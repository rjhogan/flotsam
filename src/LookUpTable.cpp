/** @file      LookUpTable.cpp
    @brief     Defines the struct LookUpTable
    @copyright 2017 European Centre for Medium Range Weather Forcasts
    @license   Apache License Version 2 (see the NOTICE.md file for details)
 */

// Includes the declaration for reconstruct_phase_function
#include <flotsam.hpp>

#include "LookUpTable.hpp"

//#include "phase_function_isotropic.h"
//#include "phase_function_cos_squared.h"
//#include "phase_function_exptan_25deg.h"
//#include "phase_function_gaussian_35deg_reducedsmooth.h"
#include "asymmetry_factor_from_phase_func.hpp"
//#include "phase_function_exptan_35deg_reducedsmooth.hpp"
#include "phase_function_component.hpp"
#include "angular_variance_from_phase_func.hpp"

namespace flotsam {

  /// Single global instance of LookUpTable accessible with
  /// flotsam::lut
  const LookUpTable lut;

  LookUpTable::LookUpTable() : lobe_spread(0.15),
			       forward_lobe_diffuse(FLOTSAM_FORWARD_LOBE_DIFFUSE),
			       forward_lobe_remain(FLOTSAM_FORWARD_LOBE_REMAIN*(1.0-forward_lobe_diffuse)),
			       i_lobe_component(-1) {
    Vector theta;
    Vector phi;
    int ntheta = 90, nphi = 180;
    nza = 101;

    theta = linspace(0.0, 90.0, ntheta) * M_PI / 180.0;
    phi   = linspace(0.0, 180.0, nphi)  * M_PI / 180.0;
    Vector mu0 = linspace(0.0, 1.0, nza);

    dmu = mu0(1);

    Vector weight = exp(-tan(theta) / (FLOTSAM_LOBE_HALF_WIDTH_DIFFUSIVITY_DEG*M_PI/180.0)) * sin(theta);
    weight(end) = 0.0;

    diffusivity_dn.resize(nza);
    diffusivity_up.resize(nza);
    fraction_dn.resize(nza);

    for (int iza = 0; iza < nza; ++iza) {
      Real za = acos(mu0(iza));
      Real weight_dn_sum = 0.0, transmissivity_dn_sum = 0.0;
      Real weight_up_sum = 0.0, transmissivity_up_sum = 0.0;
      // FIX
      //      Real od_ref = 2.0;//2.0 * mu0(iza);
      //      Real od_ref = 1.0;//2.0 * mu0(iza);
      //Real od_ref = 1.0 * mu0(iza);
      //Real od_ref = mu0(iza);

      //      std::cerr << "Warning: standard od_ref changed\n";
      //Real od_ref = 0.5*mu0(iza); // Standard up to 0.5.13
      Real od_ref = FLOTSAM_LUT_OPTICAL_DEPTH_REF_FACTOR * mu0(iza);

      //Real od_ref = 0.5*sqrt(mu0(iza));
      if (od_ref == 0.0) {
	od_ref = 1.0e-5;
      }
      for (int itheta = 0; itheta < ntheta; ++itheta) {
	for (int iphi = 0; iphi < nphi; ++iphi) {
	  // Zenith angle of ray
	  Real mu = cos(theta(itheta))*mu0(iza) 
	    + sin(theta(itheta))*cos(phi(iphi))*sin(za);
	  if (mu > 0.0) {
	    weight_dn_sum += weight(itheta);
	    transmissivity_dn_sum += weight(itheta) * exp(-od_ref/mu);
	  }
	  else if (mu < 0.0) {
	    weight_up_sum += weight(itheta);
	    transmissivity_up_sum += weight(itheta) * exp(od_ref/mu);
	  }
	}
      }
      fraction_dn(iza) = weight_dn_sum / (weight_up_sum + weight_dn_sum);
      diffusivity_dn(iza) = -log(transmissivity_dn_sum / weight_dn_sum) / od_ref;
      if (weight_up_sum > 0.0) {
	diffusivity_up(iza) = -log(transmissivity_up_sum / weight_up_sum) / od_ref;
      }
      else {
	diffusivity_up(iza) = diffusivity_up(iza-1);
      }
    }
#ifdef FLOTSAM_REVISION3A
    // FIX try alternative simpler diffusivity formulation 
    // diffusivity_dn = 1.0 / (0.15 + mu0*0.7);
    diffusivity_dn = 1.0 / (FLOTSAM_MIN_DIFFUSIVITY_MU_LOBE_DN
			    + mu0*(0.85-FLOTSAM_MIN_DIFFUSIVITY_MU_LOBE_DN));
    //diffusivity_dn = 1.0 / mu0; // FIXFIXFIX
    //    diffusivity_dn = 1.0 / (mu0*0.85);
    //diffusivity_dn = 1.0 / (0.1+mu0*0.75); // AltSec4
    //    diffusivity_up = 1.0 / 0.4;
    diffusivity_up = 1.0 / FLOTSAM_MU_LOBE_UP;
#endif 
#ifdef FLOTSAM_USE_OLD_COMPONENTS
    // Components up to version 0.5.15
    flotsam_phase_func_component_t id[] = { FLOTSAM_PHASE_FUNC_ISOTROPIC,
					    FLOTSAM_PHASE_FUNC_RAYLEIGH,
					    FLOTSAM_PHASE_FUNC_CONVEX_LOBE };
    set_phase_function_components(3, id);
#else
    // Improved components
    flotsam_phase_func_component_t id[] = { FLOTSAM_PHASE_FUNC_ISOTROPIC,
					    FLOTSAM_PHASE_FUNC_CONVEX_LOBE,
					    FLOTSAM_PHASE_FUNC_COS2_FORWARD,
					    FLOTSAM_PHASE_FUNC_COS2_BACKWARD,
					    FLOTSAM_PHASE_FUNC_BACKWARD,
					    FLOTSAM_PHASE_FUNC_SLOPE };
    set_phase_function_components(6, id);
#endif
  }


  void LookUpTable::set_phase_function_components(int n_components,
						  flotsam_phase_func_component_t* id)
  {
    i_lobe_component = -1;
    phase_function_component_id.resize(n_components);
    intVector mass_component_index_(n_components);
    int imass = 0;
    rayleigh_components.resize(n_components);
    rayleigh_components = 0.0;
    phase_function_component_index.resize(FLOTSAM_PHASE_FUNC_MAX_COMPONENTS);
    phase_function_component_index = -1;
    for (int i = 0; i < n_components; ++i) {
      phase_function_component_id[i] = id[i];
      phase_function_component_index[id[i]] = i;
      if (id[i] == FLOTSAM_PHASE_FUNC_CONVEX_LOBE) {
	i_lobe_component = i;
      }
      else if (id[i] == FLOTSAM_PHASE_FUNC_ISOTROPIC) {
	rayleigh_components[i] = 0.75;
      }
      else if (id[i] == FLOTSAM_PHASE_FUNC_COS2_FORWARD
	       ||id[i] == FLOTSAM_PHASE_FUNC_COS2_BACKWARD) {
	rayleigh_components[i] = 0.125;
      }
      else if (id[i] == FLOTSAM_PHASE_FUNC_RAYLEIGH) {
	rayleigh_components[i] = 0.25;
      }

      if (id[i] != FLOTSAM_PHASE_FUNC_SLOPE) {
	// This phase function component integrates to one
	mass_component_index_[imass++] = i;
      }
    }
    mass_component_index.resize(imass);
    mass_component_index = mass_component_index_(range(0,imass-1));

    // Initialize phase function components
    //nang = 181;
    nang = 91;
    // If dimension ordering is changed, also change n_phase_function_components()
    // member function:
    phase_function_components.resize(nang, n_components, FLOTSAM_NUM_SMOOTHING_VALUES+1);

    for (int i = 0; i < n_components; ++i) {
      phase_function_components(__, i, __)
	= phase_function_component(nang, id[i]);
    }

    dang = M_PI / (nang-1);

    g_components = asymmetry_factor_from_phase_func(phase_function_components(__,__,0).T());
    
    // Compute angular variances of the two forward patterns, if used
    int i_basis = phase_function_component_index[FLOTSAM_PHASE_FUNC_CONVEX_LOBE];
    if (i_basis >= 0) {
      ang_var_convex_lobe = angular_variance_from_phase_func(phase_function_components(__,i_basis,0));
    }
    i_basis = phase_function_component_index[FLOTSAM_PHASE_FUNC_COS2_FORWARD];
    if (i_basis >= 0) {
      ang_var_cos2_forward = angular_variance_from_phase_func(phase_function_components(__,i_basis,0));
    }

    ang_var_smoothing = ang_var_convex_lobe * range(0,FLOTSAM_NUM_SMOOTHING_VALUES);

    //std::cout << "g_components = " << g_components << "\n";

    //    g_back_hem_components 
    //      = asymmetry_factor_from_phase_func(phase_function_components(__,__,0).T(), true);
  }


  void LookUpTable::get_lobe_props(Real mu,           ///< Zenith angle of beam
				   Real& sec_lobe_dn, ///< Effective secant of downwelling part of lobe
				   Real& sec_lobe_up, ///< Effective secant of upwelling part of lobe
				   Real& frac_dn      ///< Fraction of lobe in downwelling part
				   ) const {

    Real real_index = mu / dmu;
    int i1 = real_index;
    if (i1 < 0) {
      i1 = 0;
    }
    else if (i1 >= nza-1) {
      i1 = nza-2;
    }
    Real wt1 = i1 + 1.0 - real_index;
    Real wt2 = 1.0 - wt1;
    sec_lobe_dn = diffusivity_dn(i1)*wt1 + diffusivity_dn(i1+1)*wt2;
    sec_lobe_up = diffusivity_up(i1)*wt1 + diffusivity_up(i1+1)*wt2;
    frac_dn = fraction_dn(i1)*wt1 + fraction_dn(i1+1)*wt2;

    // FIX
    //Ddn = 1.0 / mu0;
    //    Dup = Ddn;
    // frac_dn = 1.0;
    //Dup = 2.0;
  }

  /// Get the "raw", smoothed and twice-smoothed phase functions
  /// corresponding to the phase-function component amplitudes
  /// provided in pf_components
  Matrix LookUpTable::get_phase_function(const Vector& pf_components) const {
    Matrix M(phase_function_components.size(0), phase_function_components.size(2));
    M = 0.0;
    for (int i = 0; i < pf_components.size(); ++i) {
      M += pf_components(i) * phase_function_components(__,i,__);
    }
    return M;
  }

  /// Get the "raw", smoothed and twice-smoothed phase functions
  /// corresponding to the phase-function component amplitudes
  /// provided in pf_components, interpolated regularly in angle from
  /// 0 to 180 degrees with "nang" elements.
  Matrix reconstruct_phase_function(int nang, const Vector& pf_components) {
    Matrix pf = lut.get_phase_function(pf_components);
    Matrix pf_out(pf.size(1), nang);
    Vector ang_old = linspace(0.0,180.0,pf.size(0));
    Vector ang_new = linspace(0.0,180.0,nang);
    for (int i = 0; i < pf.size(1); ++i) {
      pf_out(i,__) = interp(ang_old, pf(__,i), ang_new);
    }
    // Lastly, account for the fact that the delta component is missing
    pf_out /= sum(pf_components(lut.mass_component_index));

    return pf_out;
  }


  /// Return indices to those components used to reconstruct the phase
  /// function that integrate to 1 (as opposed to one of them that
  /// integrates to zero)
  intVector mass_component_index() { return lut.mass_component_index; }

  /// Return a matrix containing the normalized basis functions used
  /// to parameterize phase functions
  Matrix basis_functions() { return lut.phase_function_components(__,__,0); }

}
