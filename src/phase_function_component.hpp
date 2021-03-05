/** @file      phase_function_component.hpp
    @brief     Header for computing components used for reconstructing smoothed phase functions
    @copyright 2017 European Centre for Medium Range Weather Forcasts
    @license   Apache License Version 2 (see the NOTICE.md file for details)
 */

#ifndef FLOTSAM_PHASE_FUNCTION_COMPONENT_HPP
#define FLOTSAM_PHASE_FUNCTION_COMPONENT_HPP

#include <flotsam.h>

namespace flotsam {
  using namespace adept;
  
  /// Return a phase function component and its once- and
  /// twice-smoothed versions for "nang" evenly spaced angles between
  /// 0 and pi, where icomponent is one of the list provided in
  /// flotsam.h
  Matrix phase_function_component(int nang, flotsam_phase_func_component_t icomponent);

  /// Return a phase function corresponding to smoothing by the
  /// forward lobe, using the exp(-tan(angle)/width)) parametric form. 
  Vector exptan_smoothing_kernel(int nang,        ///< Number of angles between 0 and pi
				 Real width_deg); ///< Width in degrees of the kernel

}

#endif
