/// @file      rayleigh.hpp
/// @brief     Properties of the Rayleigh-scattering atmosphere
/// @copyright 2017 European Centre for Medium Range Weather Forcasts
/// @license   Apache License Version 2 (see the NOTICE.md file for details)

#ifndef FLOTSAM_RAYLEIGH_HPP
#define FLOTSAM_RAYLEIGH_HPP

#include "base.hpp"

namespace flotsam {

  /// Compute the value of the Rayleigh phase function at the specified
  /// angle
  Real rayleigh_phase_function(Real ang ///< Scattering angle (radians)
			       );

  /// Compute the Rayleigh-scattering optical depth per Pascal of dry air
  Real rayleigh_optical_depth_per_Pa(Real wavelength ///< Wavelength (m)
				     );

}

#endif
