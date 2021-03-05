/** @file      angular_variance_from_phase_func.hpp
    @brief     Compute angular variance from phase function
    @copyright 2018 European Centre for Medium Range Weather Forcasts
    @license   Apache License Version 2 (see the NOTICE.md file for details)
 */
#ifndef FLOTSAM_ANGULAR_VARIANCE_FROM_PHASE_FUNC_HPP
#define FLOTSAM_ANGULAR_VARIANCE_FROM_PHASE_FUNC_HPP

#include "base.hpp"

namespace flotsam {

  /// Compute angular variance from raw phase functions, in sterad
  ///
  /// The vector pf contains a phase function evenly spaced in angle
  /// from 0 to 180
  Real angular_variance_from_phase_func(const Vector& pf ///< Phase function
					);

}

#endif
