/** @file      calc_reflectance_from_fluxes.hpp
    @brief     Compute reflectance from direct, lobe and diffuse fluxes
    @copyright 2017 European Centre for Medium Range Weather Forcasts
    @license   Apache License Version 2 (see the NOTICE.md file for details)
 */
#ifndef FLOTSAM_CALC_REFLECTANCE_FROM_FLUXES_HPP
#define FLOTSAM_CALC_REFLECTANCE_FROM_FLUXES_HPP

#include "containers.hpp"

namespace flotsam {
  template <bool IsActive>
  typename scalar<IsActive>::type
  calc_reflectance_from_fluxes(const Array<1,Real,IsActive>& surface_albedo,
			       const ScatteringProperties<IsActive>& prop,
			       const BeamFluxes<IsActive>& beam,
			       const DiffuseFluxes<IsActive>& diffuse,
			       const BeamFluxes<IsActive>& ret,
			       const Array<3,Real,IsActive>& beam_to_beam,
			       Real* ref_components = 0) {

    typedef typename scalar<IsActive>::type areal;
    typedef Array<1,Real,IsActive> avector;

    static const int beam_map[5] = {0, 1, 1, 2, 2};

    areal ref         = 0.0;
    Real  ref_direct  = 0.0;
    Real  ref_lobe_dn = 0.0;
    Real  ref_lobe_up = 0.0;
    Real  ref_diffuse = 0.0;
    Real  ref_surface = 0.0;

    //    avector factor = prop.ssa / (2.0 * M_PI * (1.0 - prob.direct));
    //    avector factor = prop.ssa * prop.od / (2.0 * M_PI); // FIX
    avector factor = prop.ssa * prop.od / (4.0 * M_PI);

    // For isotropic and Rayleigh phase functions, the flux from
    // integrating over the radiance distribution ought to match the
    // TOA diffuse flux.
    areal total_od = sum(prop.od);
#ifdef FLOTSAM_DIFFUSE_SCALING_WITH_OD
    /*
    areal scaling_max = 1.0 + (FLOTSAM_DIFFUSE_SCALING_MAX-1.0) * diffuse.dn(end) / (diffuse.dn(end)+diffuse.up(end));
    areal diffuse_multiplier = 
        (total_od+FLOTSAM_DIFFUSE_SCALING_OD_SCALE*scaling_max)
      / (total_od+FLOTSAM_DIFFUSE_SCALING_OD_SCALE);
    */
    // A good fit to the analytic solution for effective diffusivity
    // is log(2/od + exp(2)), but with an optional additional term to
    // ensure the multiplier does not exceed 2
    static const Real diffuse_ref = 0.0; //FLOTSAM_DIFFUSE_MULTIPLIER_REFERENCE_OD*0.02;
    areal diffuse_multiplier
      = log(FLOTSAM_DIFFUSE_MULTIPLIER_REFERENCE_OD/sqrt(total_od*total_od + diffuse_ref*diffuse_ref)
	    + exp(2.0 * FLOTSAM_NONDIRECT_SCALING)) * 0.5;
    // But surface term is isotropic (diffuse_multiplier = 1) so need
    // to weight
    diffuse_multiplier = (diffuse.total_source * diffuse_multiplier + diffuse.up(end))
      / (diffuse.total_source + diffuse.up(end));
    
#else
    Real diffuse_multiplier = 1.0;
#endif

#ifdef FLOTSAM_LOBE_SCALING_WITH_OD
    areal sin_sza = sqrt(1.0 - 1.0/(beam.sec(0)*beam.sec(0)));
    static const Real lobe_ref = 0.0; //FLOTSAM_LOBE_MULTIPLIER_REFERENCE_OD*0.02;
    areal lobe_multiplier
      = log(FLOTSAM_LOBE_MULTIPLIER_REFERENCE_OD/sqrt(total_od*total_od + lobe_ref*lobe_ref)
	    + exp(2.0)) * 0.5;
    lobe_multiplier = 1.0 + (lobe_multiplier-1.0)*sin_sza;
#else
    Real lobe_multiplier = 1.0;
#endif

    areal single_ref;

    // Compute the reflectance due to interactions between the
    // components of the beam fluxes
    for (int i = 0; i < beam.n_components(); ++i) {
      for (int j = 0; j < ret.n_components(); ++j) {
	// std::cout << "  b2b = " << 0.5*beam_to_beam(beam_map[i],beam_map[j],__) << "\n";
	// std::cout << "  sec = " << beam.sec(i) << " " << ret.sec(j) << "\n";
	// std::cout << "  factor = " << factor*2.0/prop.od << "\n";
	// std::cout << "  beam.prefactor(i,__) = " << beam.prefactor(i,__) << "\n";
	// std::cout << "  ret.prefactor(j,__)  = " << ret.prefactor(j,__) << "\n";
	/*
	if (i == 0 && j == 0) {
	  lobe_multiplier = 1.0;
	}
	else {
	  lobe_multiplier = diffuse_multiplier;
	}
	*/
	// Save the direct reflectance
	if (i == 0 && j == 0) {
	  ref = beam.sec(i) * ret.sec(j)
	    * sum(beam_to_beam(beam_map[i], beam_map[j], __)
		*factor*beam.prefactor(i,__)*ret.prefactor(j,__)
		*(1.0 - beam.expm_exponent(i,__)*ret.expm_exponent(j,__))
		/(beam.exponent(i,__)+ret.exponent(j,__)));
	  ref_direct = value(ref);
	}
	else {
	  single_ref = lobe_multiplier * beam.sec(i) * ret.sec(j)
	    * sum(beam_to_beam(beam_map[i], beam_map[j], __)
		  *factor*beam.prefactor(i,__)*ret.prefactor(j,__)
		  *(1.0 - beam.expm_exponent(i,__)*ret.expm_exponent(j,__))
		  /(beam.exponent(i,__)+ret.exponent(j,__)));
	  ref += single_ref;
	  if (i >= 3 || j >= 3) {
	    ref_lobe_up += value(single_ref);
	  }
	  else {
	    ref_lobe_dn += value(single_ref);
	  }
	}

	/*
	std::cout << "REFLECTANCE_COMPONENT(" << i << "," << j << ") = " <<
	beam.sec(i) * ret.sec(j)
	  * sum(beam_to_beam(beam_map[i], beam_map[j], __)
		*factor*beam.prefactor(i,__)*ret.prefactor(j,__)
		*(1.0 - exp(-(beam.exponent(i,__)+ret.exponent(j,__))))
		/(beam.exponent(i,__)+ret.exponent(j,__))) << "\n";
	std::cout << "  multiplier = " << beam.sec(i) * ret.sec(j)*beam_to_beam(beam_map[i], beam_map[j], __)
	  *factor/prop.od << "\n";

	*/
      }
    }

    avector grad_diffuse_dn(prop.size());
    avector grad_diffuse_up(prop.size());

    grad_diffuse_dn = diffuse.dn(range(1,end)) - diffuse.dn(range(0,end-1));
    grad_diffuse_up = diffuse.up(range(1,end)) - diffuse.up(range(0,end-1));

    // Compute the reflectance due to the diffuse beam being received
    // by the returning beam

    // FIX multiply by factor for Rayleigh case?
    avector one_over_exponent;
    /*
    for (int i = 0; i < ret.n_components(); ++i) {
      ref += sum((ret.prefactor_source_dn(i,__)
		  *((diffuse.dn(range(0,end-1))*ret.exponent(i,__) + grad_diffuse_dn)
		    -(diffuse.dn(range(1,end)) *ret.exponent(i,__) + grad_diffuse_dn)
		    *exp(-ret.exponent(i,__)))
		  // FIX add factor 1.125 here?
		  +ret.prefactor_source_up(i,__)
		  *((diffuse.up(range(0,end-1))*ret.exponent(i,__) + grad_diffuse_up)
		    -(diffuse.up(range(1,end)) *ret.exponent(i,__) + grad_diffuse_up)
		    *exp(-ret.exponent(i,__))))
		 * prop.od*prop.od*prop.ssa / (ret.exponent(i,__)*ret.exponent(i,__))) / M_PI;
    }
    */

    for (int i = 0; i < ret.n_components(); ++i) {
      one_over_exponent = 1.0 / ret.exponent(i,__);
      /*
      areal diff_up = diffuse_multiplier * sum(
		 (
		  FLOTSAM_DIFFUSE_DN_FACTOR*
		  ret.prefactor_source_up(i,__)
		  *((diffuse.dn(range(0,end-1)) + grad_diffuse_dn*one_over_exponent)
		    -(diffuse.dn(range(1,end))  + grad_diffuse_dn*one_over_exponent)
		    *ret.expm_exponent(i,__)))
		 * prop.ssa * prop.od * one_over_exponent) / M_PI;
      areal diff_dn = diffuse_multiplier * sum(
		 (FLOTSAM_DIFFUSE_UP_FACTOR*
		  ret.prefactor_source_dn(i,__)
		  *((diffuse.up(range(0,end-1)) + grad_diffuse_up*one_over_exponent)
		    -(diffuse.up(range(1,end))  + grad_diffuse_up*one_over_exponent)
		    *ret.expm_exponent(i,__)))
		 * prop.ssa * prop.od * one_over_exponent) / M_PI;
      std::cout << "DIFFUSE BITS: " << i << " " << diff_up << " " << diff_dn << "\n";
      std::cout << "  prefactor_source_* " << ret.prefactor_source_up(i,__) << " " << ret.prefactor_source_dn(i,__) << "\n";
      std::cout << "  grad_diffuse_* " << grad_diffuse_dn << " " << grad_diffuse_up << "\n";
      std::cout << "  one_over_exponent " << one_over_exponent << " " << ret.expm_exponent(i,__) << "\n";
      std::cout << "  prop " << prop.ssa << " " << prop.od << "\n";
      */

      // Note that in the following, ret.prefactor_source_up is
      // associated with diffuse.dn because the sense of the former is
      // inverted by virtue of it being a returning probability field
      // that acts in an analogous but inverted sense from its
      // outgoing equivalent.

      ref += diffuse_multiplier * sum(
		 (
		  FLOTSAM_DIFFUSE_DN_FACTOR*
		  ret.prefactor_source_up(i,__)
		  *((diffuse.dn(range(0,end-1)) + grad_diffuse_dn*one_over_exponent)
		    -(diffuse.dn(range(1,end))  + grad_diffuse_dn*one_over_exponent)
		    *ret.expm_exponent(i,__))
		  +
		  FLOTSAM_DIFFUSE_UP_FACTOR*
		  ret.prefactor_source_dn(i,__)
		  *((diffuse.up(range(0,end-1)) + grad_diffuse_up*one_over_exponent)
		    -(diffuse.up(range(1,end))  + grad_diffuse_up*one_over_exponent)
		    *ret.expm_exponent(i,__)))
		 * prop.ssa * prop.od * one_over_exponent) / M_PI;
      /*
      std::cout << "DIFFUSE FROM DOWN " <<
	sum((ret.prefactor_source_up(i,__)
		  *((diffuse.dn(range(0,end-1)) + grad_diffuse_dn*one_over_exponent)
		    -(diffuse.dn(range(1,end))  + grad_diffuse_dn*one_over_exponent)
		    *exp(-ret.exponent(i,__))))*prop.ssa*prop.od*one_over_exponent)/M_PI
		<< "\n";
      std::cout << "DIFFUSE FROM UP   " <<
		  // FIX add factor 1.125 here?
	sum((1.125*ret.prefactor_source_dn(i,__)
		  *((diffuse.up(range(0,end-1)) + grad_diffuse_up*one_over_exponent)
		    -(diffuse.up(range(1,end))  + grad_diffuse_up*one_over_exponent)
		    *exp(-ret.exponent(i,__))))
		 * prop.ssa * prop.od * one_over_exponent) / M_PI << "\n";
		 */
    }
    ref_diffuse = value(ref) - ref_lobe_dn -ref_lobe_up - ref_direct;

    

    // Surface reflectance
    //    ref += diffuse.up(end) * ret.dn(ret.size()) / M_PI;
    //    ref += ret.dn(ret.size()) * ( (diffuse.up(end)-beam.dn(beam.size())*surface_albedo(0)) / M_PI + beam.dn(beam.size())*bdrf);


    if (surface_albedo.size() > 1){
      //      std::cerr << "brdf diff to view: " << surface_albedo(2) << " brdf hemispheric: " << surface_albedo(0) << "\n";
      //      std::cerr << "diffuse.dn(end): " << diffuse.dn(end) << "\n";
      //      std::cerr << "ret.dn(ret.size()): " << ret.dn(ret.size()) << "\n";
      ref += ret.dn(ret.size()) * ( (diffuse.dn(end))*surface_albedo(2) + beam.dn(beam.size())*surface_albedo(3));
      //      ref += ret.dn(ret.size()) * ( (diffuse.dn(end)/M_PI)*surface_albedo(2) + beam.dn(beam.size())*surface_albedo(3));
    }
    else {
      ref += diffuse.up(end) * ret.dn(ret.size()) / M_PI;
    }
    
    
    // following is to check that this amount to zero. Ok. 
    //    ref +=  diffuse.up(end)-(diffuse.dn(end)*surface_albedo(0) + beam.dn(beam.size())*surface_albedo(1));


    ref_surface = value(ref) - ref_diffuse - ref_lobe_dn - ref_lobe_up - ref_direct;

    if (ref_components) {
      ref_components[0] = ref_direct;
      ref_components[1] = ref_lobe_dn;
      ref_components[2] = ref_lobe_up;
      ref_components[3] = ref_diffuse;
      ref_components[4] = ref_surface;
    }
    /*
    std::cerr << "ref tot: " << value(ref) << "\n";
    std::cerr << "ref direct: " << ref_direct << "\n";
    std::cerr << "ref lobe: " << ref_lobe << "\n";
    std::cerr << "ref surface: " << ref_surface << "\n";
    std::cerr << "ref diffuse: " << ref_diffuse << "\n";
    */

    return ref;

  }
}

#endif
