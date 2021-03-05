/** @file      calc_scattering_probabilities.hpp
    @brief     Calculate probabilities that direct & forward lobe are scattered into lobe or down/up hemispheres
    @copyright 2017 European Centre for Medium Range Weather Forcasts
    @license   Apache License Version 2 (see the NOTICE.md file for details)
 */
#ifndef FLOTSAM_CALC_SCATTERING_PROBABILITIES_HPP
#define FLOTSAM_CALC_SCATTERING_PROBABILITIES_HPP

#include "containers.hpp"
#include "LookUpTable.hpp"

namespace flotsam {

  template <bool IsActive>
  void calc_scattering_probabilities(Real mu, Real mu_lobe_dn,
				     Real mu_lobe_up, Real frac_lobe_dn,
				     const Array<2,Real,IsActive>& pfc,
				     ScatteringProbabilities<IsActive>& prob_direct_to,
				     ScatteringProbabilities<IsActive>& prob_lobe_dn_to,
				     ScatteringProbabilities<IsActive>& prob_lobe_up_to) {

    typedef typename scalar<IsActive>::type areal;
    typedef Array<1,Real,IsActive> avector;


    prob_direct_to.lobe       = 0.0;
    
    // Compute probability of diffuse scattering into forward/backward
    // hemisphere in beam frame of reference, then translate to
    // upward/downward hemisphere
    avector prob_direct_to_fwd(prob_direct_to.diffuse_dn.size());
    avector prob_direct_to_bwd(prob_direct_to.diffuse_up.size());
    prob_direct_to_fwd = 0.0;
    prob_direct_to_bwd = 0.0;

    // Loop through the phase function components in "pfc" and add to
    // the appropriate variable
    for (int i = 0; i < lut.n_phase_function_components(); ++i) {
      switch (lut.phase_function_component_id[i]) {
      case FLOTSAM_PHASE_FUNC_ISOTROPIC:
      case FLOTSAM_PHASE_FUNC_RAYLEIGH:
	prob_direct_to_fwd += 0.5 * pfc(__,i);
	prob_direct_to_bwd += 0.5 * pfc(__,i);
	break;
      case FLOTSAM_PHASE_FUNC_CONVEX_LOBE:
	prob_direct_to.lobe += pfc(__,i);
	break;
      case FLOTSAM_PHASE_FUNC_COS2_FORWARD:
	{
	  int i_cos2_bwd
	    = lut.phase_function_component_index[FLOTSAM_PHASE_FUNC_COS2_BACKWARD];
	  avector to_lobe;
	  if (i_cos2_bwd < 0) {
	    to_lobe = FLOTSAM_COS2_FORWARD_TO_LOBE * pfc(__,i);
	  }
	  else {
	    to_lobe = FLOTSAM_COS2_FORWARD_TO_LOBE
	      * (pfc(__,i)-pfc(__,i_cos2_bwd));
	  }
	  prob_direct_to.lobe += to_lobe;
	  prob_direct_to_fwd += (pfc(__,i) - to_lobe);
	}
	break;
      case FLOTSAM_PHASE_FUNC_COS2_BACKWARD:
      case FLOTSAM_PHASE_FUNC_BACKWARD:
	prob_direct_to_bwd += pfc(__,i);
	break;
      case FLOTSAM_PHASE_FUNC_SLOPE:
	prob_direct_to_fwd += 0.5 * pfc(__,i);
	prob_direct_to_bwd -= 0.5 * pfc(__,i);
	break;
      }
    }


#ifndef FLOTSAM_REVISION4
    // Surely this is wrong
    prob_direct_to.diffuse_dn = prob_direct_to_fwd;
    prob_direct_to.diffuse_up = prob_direct_to_bwd;
    areal weight = 0.5 + 0.5*mu_lobe_dn;
#else
    areal weight = 0.5 + 0.5*mu;
    prob_direct_to.diffuse_dn = weight*prob_direct_to_fwd + (1.0-weight)*prob_direct_to_bwd;
    prob_direct_to.diffuse_up = weight*prob_direct_to_bwd + (1.0-weight)*prob_direct_to_fwd;
    weight = 0.5 + 0.5*mu_lobe_dn;
#endif

    // Account for the possibility that we treat some of the lobe
    // scattering as going into the forward diffuse stream
    prob_direct_to.diffuse_dn 
      += (weight * flotsam::lut.forward_lobe_diffuse) * prob_direct_to.lobe;
    prob_direct_to.diffuse_up 
      += ((1.0-weight) * flotsam::lut.forward_lobe_diffuse) * prob_direct_to.lobe;
    prob_direct_to.lobe *= (1.0 - flotsam::lut.forward_lobe_diffuse);

    Real forward_lobe_remain = flotsam::lut.forward_lobe_remain * frac_lobe_dn;
    prob_lobe_dn_to.lobe = prob_direct_to.lobe * forward_lobe_remain;
    prob_lobe_up_to.lobe = prob_direct_to.lobe * forward_lobe_remain;

    avector prob_lobe_to_fwd_hem = prob_direct_to.diffuse_dn
      + (1.0-flotsam::lut.lobe_spread) * (1.0-forward_lobe_remain)*prob_direct_to.lobe;
    avector prob_lobe_to_backwd_hem = prob_direct_to.diffuse_up
      + flotsam::lut.lobe_spread*(1.0-forward_lobe_remain)*prob_direct_to.lobe;

    prob_lobe_dn_to.diffuse_dn = weight * prob_lobe_to_fwd_hem
      + (1.0 - weight) * prob_lobe_to_backwd_hem;
    prob_lobe_dn_to.diffuse_up = weight * prob_lobe_to_backwd_hem
      + (1.0 - weight) * prob_lobe_to_fwd_hem;

    // Note that mu_lobe_up should be negative here
    weight = 0.5 + 0.5*mu_lobe_up;
    // For an outgoing beam, prob_lobe_up_to.diffuse_dn is the
    // probability of the downward directed lobe being scattered into
    // the (primary) downwelling diffuse flux, whereas for a beam
    // returning to an instrument, it is proportional to the
    // probability of the upwelling diffuse flux being scattered into
    // the (primary) lobe containing radiation heading upwards.
    prob_lobe_up_to.diffuse_dn = weight * prob_lobe_to_fwd_hem
      + (1.0 - weight) * prob_lobe_to_backwd_hem;
    prob_lobe_up_to.diffuse_up = weight * prob_lobe_to_backwd_hem
      + (1.0 - weight) * prob_lobe_to_fwd_hem;
  }

}

#endif
