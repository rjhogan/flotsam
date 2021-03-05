/** @file      calc_asymmetry_factor.hpp
    @brief     Compute asymmetry factor from phase-function components
    @copyright 2017 European Centre for Medium Range Weather Forcasts
    @license   Apache License Version 2 (see the NOTICE.md file for details)
 */
#ifndef FLOTSAM_CALC_ASYMMETRY_FACTOR_HPP
#define FLOTSAM_CALC_ASYMMETRY_FACTOR_HPP

#include "base.hpp"
#include "LookUpTable.hpp"

namespace flotsam {

  /// Compute particulate asymmetry factor from the phase function components
  ///
  /// The matrix pfc contains the fraction of the phase function
  /// contributed by each component, except the diffraction
  /// quasi-delta function in the forward direction, which is assumed
  /// to be the residual. The quasi-delta function is assumed to have
  /// an asymmetry factor of 1.0, and a matrix multiplication is done
  /// to add the contributions from the other components.
  template <bool IsActive>
  void calc_asymmetry_factor(const Array<2,Real,IsActive>& pfc, ///< Phase function components
			     Array<1,Real,IsActive>& g)         ///< Output asymmetry factor
  {
    if (!pfc.empty()) {
      g = (1.0 - sum(pfc(__,lut.mass_component_index),1)) + (pfc ** lut.g_components);
    }
  }

}


#endif
