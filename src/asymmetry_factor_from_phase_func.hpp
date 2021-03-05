/** @file      asymmetry_factor_from_phase_func.hpp
    @brief     Compute asymmetry factor from phase function
    @copyright 2017 European Centre for Medium Range Weather Forcasts
    @license   Apache License Version 2 (see the NOTICE.md file for details)
 */
#ifndef FLOTSAM_ASYMMETRY_FACTOR_FROM_PHASE_FUNC_HPP
#define FLOTSAM_ASYMMETRY_FACTOR_FROM_PHASE_FUNC_HPP

#include "base.hpp"

namespace flotsam {

  /// Compute asymmetry factor from raw phase functions
  ///
  /// The matrix pf contains phase functions evenly spaced in angle
  /// from 0 to 180, with angle being the second, fastest varying
  /// dimension of the matrix
  Vector asymmetry_factor_from_phase_func(const Matrix& pf, ///< Phase functions
					  bool back_hem_only = false);

}


#endif
