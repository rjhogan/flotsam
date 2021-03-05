/// @file      Channel.hpp
/// @brief     Declares the class Channel
/// @copyright 2017 European Centre for Medium Range Weather Forcasts
/// @license   Apache License Version 2 (see the NOTICE.md file for details)

#ifndef FLOTSAM_CHANNEL_HPP
#define FLOTSAM_CHANNEL_HPP

#include <vector>

#include "base.hpp"
#include "rayleigh.hpp"

namespace flotsam {

  /// Structure for holding information about a particular satellite
  /// channel
  struct Channel {
    Channel() : is_initialized(false), is_rayleigh_only(false) { }
    ~Channel() { };

    int read(const char* filename) { return FLOTSAM_FEATURE_NOT_AVAILABLE; };
    
    int rayleigh_only(Real wavelength_) {
      is_rayleigh_only = true;
      is_initialized = true;
      wavelength = wavelength_;
      rayleigh_od_per_Pa.resize(1);
      weight.resize(1);
      weight = 1.0;
      rayleigh_od_per_Pa = rayleigh_optical_depth_per_Pa(wavelength);
      return FLOTSAM_SUCCESS;
    }

    int vacuum() {
      is_rayleigh_only = true;
      is_initialized = true;
      wavelength = -1.0;
      rayleigh_od_per_Pa.resize(1);
      weight.resize(1);
      weight = 1.0;
      rayleigh_od_per_Pa = 0.0;
      return FLOTSAM_SUCCESS;
    }

    int n_g() const { return weight.size(); }

    Vector weight; ///< Weight of each g point
    Vector rayleigh_od_per_Pa; ///< Rayleigh-scattering optical depth per unit pressure for each g point
    Real wavelength;  ///< Wavelength (metres)
    bool is_initialized;
    bool is_rayleigh_only;
  };

}

#endif
