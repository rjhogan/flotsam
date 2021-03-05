/// @file      BackgroundProfile.hpp
/// @brief     Declares the class BackgroundProfile
/// @copyright 2017 European Centre for Medium Range Weather Forcasts
/// @license   Apache License Version 2 (see the NOTICE.md file for details)

#ifndef FLOTSAM_BACKGROUND_PROFILE_HPP
#define FLOTSAM_BACKGROUND_PROFILE_HPP

#include <iostream>

#include <flotsam.h>

#include "base.hpp"

namespace flotsam {

  using namespace adept;

  /// Holds the background properties of an atmospheric profile:
  /// pressure, temperature and the mixing ratios of absorbing gases
  struct BackgroundProfile {

    BackgroundProfile() { }

    /// Reset the size of the background profile
    int reset(int n);

    /// Set the pressure at layer edges, resetting the size of the
    /// arrays if zero
    int set_edge_pressure(int n, const Real* p);

    /// Set the mean layer temperature, only required for channels
    /// with gas absorption
    int set_temperature(int n, const Real* t);

    /// Set a profile of gas mass mixing ratios
    int set_gas_profile(flotsam_gas_t igas, int n, const Real* mmr);

    /// Set the mass mixing ratio of a particular gas to a constant
    /// value with height
    int set_gas_const(flotsam_gas_t igas, Real mmr);

    /// Return the number of layers
    int n_z() const { return edge_pressure.size()-1; }

    // Data
    Vector edge_pressure; ///< Pressure at layer edges (Pa), down from top
    Vector temperature;   ///< Layer-mean temperature (K)
    Matrix mixing_ratio;  ///< Mass mixing ratio (kg/kg) of absorbing gases
  };

}

#endif
