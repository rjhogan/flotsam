/// @file      BackgroundProfile.cpp
/// @brief     Defines the BackgroundProfile struct
/// @copyright 2017 European Centre for Medium Range Weather Forcasts
/// @license   Apache License Version 2 (see the NOTICE.md file for details)

#include "BackgroundProfile.hpp"

namespace flotsam {

  /// Reset the size of the background profile
  int BackgroundProfile::reset(int n) {
    if (n > 0) {
      edge_pressure.resize(n+1);
      temperature.resize(n);
      mixing_ratio.resize(FLOTSAM_MAX_GASES,n);
      edge_pressure = -1.0;
      temperature = -1.0;
      mixing_ratio = 0.0;
    }
    else {
      edge_pressure.clear();
      temperature.clear();
      mixing_ratio.clear();
    }
    return FLOTSAM_SUCCESS;
  }

  /// Set the pressure at layer edges, resetting the size of the
  /// arrays if zero
  int BackgroundProfile::set_edge_pressure(int n, const Real* p) {
    if (p[1] < p[0]) {
      return FLOTSAM_INPUTS_MUST_START_AT_TOA;
    }
    else if (edge_pressure.size() > 0 && n != edge_pressure.size()-1) {
      return FLOTSAM_INCORRECT_NUMBER_OF_LAYERS;
    }
    
    if (edge_pressure.empty()) {
      reset(n);
    }
    for (int i = 0; i <= n; ++i) {
      edge_pressure[i] = p[i];
    }
    
    return FLOTSAM_SUCCESS;
  }

  /// Set the mean layer temperature, only required for channels
  /// with gas absorption
  int BackgroundProfile::set_temperature(int n, const Real* t) {
    if (n_z() != n) {
      return FLOTSAM_INCORRECT_NUMBER_OF_LAYERS;
    }
    for (int i = 0; i <= n; ++i) {
      temperature[i] = t[i];
    }
    return FLOTSAM_SUCCESS;
  }

  /// Set a profile of gas mass mixing ratios
  int BackgroundProfile::set_gas_profile(flotsam_gas_t igas, int n,
					 const Real* mmr) {
    if (n != n_z()) {
      return FLOTSAM_INCORRECT_NUMBER_OF_LAYERS;
    }
    else if (igas < 0 || igas >= FLOTSAM_MAX_GASES) {
      return FLOTSAM_INVALID_GAS;
    }
    else {
      for (int i = 0; i < n; ++i) {
	mixing_ratio(static_cast<int>(igas),i) = mmr[i];
      }
      return FLOTSAM_SUCCESS;
    }
  }

  /// Set the mass mixing ratio of a particular gas to a constant
  /// value with height
  int BackgroundProfile::set_gas_const(flotsam_gas_t igas, Real mmr) {
    if (igas < 0 || igas >= FLOTSAM_MAX_GASES) {
      return FLOTSAM_INVALID_GAS;
    }
    else {
      mixing_ratio(static_cast<int>(igas),__) = mmr;
      return FLOTSAM_SUCCESS;
    }
  }


}
