/** @file      calc_angular_variance.hpp
    @brief     Compute angular variance of wide lobe from phase-function components
    @copyright 2018 European Centre for Medium Range Weather Forcasts
    @license   Apache License Version 2 (see the NOTICE.md file for details)
 */
#ifndef FLOTSAM_CALC_ANGULAR_VARIANCE_HPP
#define FLOTSAM_CALC_ANGULAR_VARIANCE_HPP

#include "base.hpp"
#include "LookUpTable.hpp"

namespace flotsam {

  /// Compute angular variance of wide forward lobe from the phase function components
  ///
  /// The matrix pfc contains the fraction of the phase function
  /// contributed by each component, except the diffraction
  /// quasi-delta function in the forward direction, which is assumed
  /// to be the residual. 
  template <bool IsActive>
  void calc_angular_variance(const Array<2,Real,IsActive>& pfc, ///< Phase function components
			     Array<1,Real,IsActive>& ang_var)   ///< Output angular variance (sterad)
  {
    int i_lobe
      = lut.phase_function_component_index[FLOTSAM_PHASE_FUNC_CONVEX_LOBE];
    int i_cos2_fwd
      = lut.phase_function_component_index[FLOTSAM_PHASE_FUNC_COS2_FORWARD];
    int i_cos2_bwd
      = lut.phase_function_component_index[FLOTSAM_PHASE_FUNC_COS2_BACKWARD];

    ang_var.resize(pfc.dimension(0));
    // Default value
    ang_var = lut.ang_var_convex_lobe;

    if (i_cos2_fwd >= 0) {
      typedef typename scalar<IsActive>::type areal;
      areal denom;
      for (int i = 0; i < pfc.dimension(0); ++i) {
	if (i_cos2_bwd < 0) {
	  denom = pfc(i,i_lobe) + pfc(i,i_cos2_fwd)*FLOTSAM_COS2_FORWARD_TO_LOBE;
	  if (denom > 0.0) {
	    ang_var(i) = (lut.ang_var_convex_lobe*pfc(i,i_lobe)
			  +lut.ang_var_cos2_forward*FLOTSAM_COS2_FORWARD_TO_LOBE*pfc(i,i_cos2_fwd))
	      / denom;
	  }
	}
	else {
	  denom = pfc(i,i_lobe) + (pfc(i,i_cos2_fwd)-pfc(i,i_cos2_bwd))*FLOTSAM_COS2_FORWARD_TO_LOBE;
	  if (denom > 0.0) {
	    ang_var(i) = (lut.ang_var_convex_lobe*pfc(i,i_lobe)
			  +lut.ang_var_cos2_forward*FLOTSAM_COS2_FORWARD_TO_LOBE
			  *(pfc(i,i_cos2_fwd)-pfc(i,i_cos2_bwd)))
	      / denom;
	  }
	}
      }
    }
  }

}


#endif
