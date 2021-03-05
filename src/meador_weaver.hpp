/** @file      meador_weaver.hpp
    @brief     Compute diffuse reflectance/transmittance of layers
    @copyright 2017 European Centre for Medium Range Weather Forcasts
    @license   Apache License Version 2 (see the NOTICE.md file for details)
 */
#ifndef FLOTSAM_MEADOR_WEAVER_HPP
#define FLOTSAM_MEADOR_WEAVER_HPP

#include "containers.hpp"

namespace flotsam {

  template <bool IsActive>
  void meador_weaver(const ScatteringProperties<IsActive>& prop,
		     const BeamFluxes<IsActive>& beam_fluxes,
		     DiffuseProperties<IsActive>& diffuse_props) {
    using std::max;

    typedef typename scalar<IsActive>::type areal;

    areal gamma1, gamma2, gamma3, gamma4, source, k2, k;
    areal exponential, exponential2, k_2_exponential;
    areal RT, RT_factor, exponential0;

    areal alpha1, alpha2, k_mu0, mu0, k_gamma3, k_gamma4;

    const int nz = prop.size();

    // Sources are formed from a summation so need to be initialized
    // to zero
    diffuse_props.source_dn = 0.0;
    diffuse_props.source_up = 0.0;

    for (int i = 0; i < nz; ++i) {

      if (prop.od(i) > 0.0) {


#ifndef FLOTSAM_NO_PIFM
	// Zdunkowski PIFM
       	gamma1 = 2.0 - prop.ssa(i) * (1.25 + 0.75*prop.g(i));
	gamma2 = 0.75*(prop.ssa(i) * (1.0 - prop.g(i)));
#else
	// Generalized two stream

	// In the following, we can get PIFM with
	// FLOTSAM_DIFFUSE_DIFFUSIVITY=2 and
	// FLOTSAM_BETA_PREFACTOR=3/8

	// gamma2 is gain due to scattering from the other stream
	gamma2 = (FLOTSAM_DIFFUSE_DIFFUSIVITY * FLOTSAM_BETA_PREFACTOR)
	  * prop.ssa(i)*(1.0-prop.g(i));
	// gamma1 is rate of loss due to absorption (first line) and
	// scattering to the other stream (=gamma2)
	gamma1 = FLOTSAM_DIFFUSE_DIFFUSIVITY * (1.0 - prop.ssa(i))
	  + gamma2;
#endif
	//gamma1 = (2.0 - prop.ssa(i)*(1.0 + prop.g(i)));
	//gamma2 = prop.ssa(i)*(1.0 - prop.g(i));
	//	gamma1 = gamma1 * (sqrt(3) / 2.0);
	//	gamma2 = gamma2 * (sqrt(3) / 2.0);
	//gamma1 = gamma1 * 1.33333;
	//gamma2 = gamma2 * 1.33333;

	k2 = (gamma1-gamma2) * (gamma1+gamma2);

	if (k2 > 1.0e-18) {
	  //	if (true) { k2 = max(1.0e-18,k2);
	  // Some absorption: use more complicated formulas
	  k = sqrt(k2);

	  exponential = exp(-k*prop.od(i));
	  exponential2 = exponential*exponential;
	  k_2_exponential = 2.0 * k * exponential;
	  RT = 1.0 / (k + gamma1 + (k - gamma1)*exponential2);
      
	  // First compute diffuse reflectance and transmittance

	  // Meador & Weaver (1980) Eq. 25
	  diffuse_props.reflectance(i) = gamma2 * (1.0 - exponential2) * RT;
	  // Meador & Weaver (1980) Eq. 26
	  diffuse_props.transmittance(i) = k_2_exponential * RT;
      
	  // Now compute sources as the summation of the contribution
	  // from scattering from the direct beam (one exponential
	  // term) and the downward and upward lobe fluxes (two
	  // exponential terms each)
	  for (int iexp = 0; iexp < beam_fluxes.n_components(); ++iexp) {
	    source = beam_fluxes.prefactor_source_dn(iexp,i)
	           + beam_fluxes.prefactor_source_up(iexp,i);
	    //	    std::cerr << "!!!complex " << i << " " << iexp << " " << source << "\n";
	    if (fabs(source) < 1.0e-10) {
	      continue;
	    }
	    gamma3 = beam_fluxes.prefactor_source_up(iexp,i) / source;
	    gamma4 = 1.0 - gamma3;

	    alpha1 = gamma1*gamma4 + gamma2*gamma3;
	    alpha2 = gamma1*gamma3 + gamma2*gamma4;
	    // The terms of the direct and diffuse beam propagate
	    // though the atmosphere at a particular zenith angle not
	    // always equal to the solar zenith angle
	    mu0 = prop.od(i) / beam_fluxes.exponent(iexp,i);
	    k_mu0 = k * mu0;
	    k_gamma3 = k*gamma3;
	    k_gamma4 = k*gamma4;
	    
	    // Note that we account for the fact that source fluxes
	    // are into a plane perpendicular to their direction
	    // (requiring a multiplication by mu0), and the source
	    // energy (source)
	    RT_factor = mu0 * source * RT * prop.ssa(i) / (1.0 - k_mu0*k_mu0);
	    exponential0 = beam_fluxes.expm_exponent(iexp,i);
	    // exponential0 = exp(-beam_fluxes.exponent(iexp,i));
	    
	    // Meador & Weaver (1980) Eq. 14, multiplying top & bottom by
	    // exp(-k*od) in case of very high optical depths
	    diffuse_props.source_up(i) += RT_factor 
	      * ( (1.0 - k_mu0) * (alpha2 + k_gamma3)
		  -(1.0 + k_mu0) * (alpha2 - k_gamma3)*exponential2
		  -k_2_exponential*(gamma3 - alpha2*mu0)*exponential0);
	    
	    // Meador & Weaver (1980) Eq. 15, multiplying top & bottom by
	    // exp(-k*od), minus the 1*exp(-od/mu0) term representing direct
	    // unscattered transmittance.
	    diffuse_props.source_dn(i) += RT_factor
	      * ( k_2_exponential*(gamma4 + alpha1*mu0)
		  - exponential0
		  * ( (1.0 + k_mu0) * (alpha1 + k_gamma4)
		      -(1.0 - k_mu0) * (alpha1 - k_gamma4) * exponential2) );
	    /*
	    if (iexp == 0) {
	      std::cerr << "!!! " << i << " " << diffuse_props.source_up(i)
			<< " " << diffuse_props.source_dn(i)
			<< " " << value(source*gamma3*prop.od(i)) << "\n";
	    }
	    */

	    /*
	    std::cerr << i << " " << iexp << " " << mu0 << " " << source << " "
		      << areal(RT_factor 
	      * ( (1.0 - k_mu0) * (alpha2 + k_gamma3)
		  -(1.0 + k_mu0) * (alpha2 - k_gamma3)*exponential2
		  -k_2_exponential*(gamma3 - alpha2*mu0)*exponential0))
		      << " "
		      << areal(RT_factor
	      * ( k_2_exponential*(gamma4 + alpha1*mu0)
		  - exponential0
		  * ( (1.0 + k_mu0) * (alpha1 + k_gamma4)
		      -(1.0 - k_mu0) * (alpha1 - k_gamma4) * exponential2) ))
		      << " " << diffuse_props.source_up(i) << " " << diffuse_props.source_dn(i) << "\n";
	    */
	  }
	}
	else {
	  // Conservative scattering: simpler formulas
	  
	  RT = gamma1*prop.od(i);
	  // First the diffuse reflectance and transmittance: Meador &
	  // Weaver (1980) Eq. 29
	  diffuse_props.reflectance(i)   = RT / (1.0 + RT);
	  diffuse_props.transmittance(i) = 1.0 - diffuse_props.reflectance(i);
      
	  // Now compute sources as the summation of the contribution
	  // from scattering from the direct beam (one exponential
	  // term) and the downward and upward lobe fluxes (two
	  // exponential terms each)
	  for (int iexp = 0; iexp < beam_fluxes.n_components(); ++iexp) {
	    mu0 = prop.od(i) / beam_fluxes.exponent(iexp,i);
	    exponential0 = beam_fluxes.expm_exponent(iexp,i);
	    source = beam_fluxes.prefactor_source_dn(iexp,i)
	           + beam_fluxes.prefactor_source_up(iexp,i);
	    //	    std::cerr << "!!!conservative " << i << " " << iexp << " " << source << "\n";
	    if (fabs(source) < 1.0e-10) {
	      continue;
	    }
	    gamma3 = beam_fluxes.prefactor_source_up(iexp,i) / source;
	    // Meador & Weaver (1980) Eq. 24
	    RT_factor = (RT + (gamma3-gamma1*mu0)*(1.0 - exponential0)) / (1.0 + RT);
	    // Note that we account for the fact that source fluxes
	    // are into a plane perpendicular to their direction
	    // (requiring a multiplication by mu0), and the source
	    // energy (source)
	    diffuse_props.source_up(i) += source * mu0 * RT_factor;
	    diffuse_props.source_dn(i) += source * mu0 * (1.0 - RT_factor - exponential0); 
	    //	    std::cerr << "!!!" << i << " " << iexp << " " << source << " " << mu0 << " " << RT_factor << " " << exponential0 << " " << RT << " " << gamma3 << " " << gamma1 << "\n";
	    //	    std::cerr << "!!! " << i << " " << iexp << " " << source << " " << value(mu0*RT_factor)
	    //		      << " " << value(mu0 * (1.0 - RT_factor - exponential0)) << " " << value(gamma3*prop.od(i)) << "\n";
	  }
	}
      }
      else {
	// In vacuum
	diffuse_props.reflectance(i) = 0.0;
	diffuse_props.transmittance(i) = 1.0;
      }
    }

    //    std::cerr << "diffuse_up = " << diffuse_props.source_up;
    //    std::cerr << "diffuse_dn = " << diffuse_props.source_dn;
    //    std::cerr << "TOTAL_BEAM_SRC_DN = " << sum(beam_fluxes.prefactor_source_dn,1) << "\n";
    //    std::cerr << "TOTAL_BEAM_SRC_UP = " << sum(beam_fluxes.prefactor_source_up,1) << "\n";
    //    std::cerr << "TOTAL_OD     = " << sum(prop.od) << "\n";
    //    std::cerr << "TOTAL_SOURCE = " << sum(diffuse_props.source_dn + diffuse_props.source_up) << "\n";

  }
}

#endif
