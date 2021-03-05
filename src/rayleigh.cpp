/// @file      rayleigh.cpp
/// @brief     Properties of the Rayleigh-scattering atmosphere
/// @copyright 2017 European Centre for Medium Range Weather Forcasts
/// @license   Apache License Version 2 (see the NOTICE.md file for details)

#include "rayleigh.hpp"

using adept::Real;

/// Compute the value of the Rayleigh phase function at the specified
/// angle
///
/// Phase function normalized so that its integral over the surface of
/// the sphere is 4pi
Real
flotsam::rayleigh_phase_function(Real ang ///< Scattering angle (radians)
				 ) {
  Real cos_ang = cos(ang);
  return (1.0 + cos_ang*cos_ang) * 0.75;
}

/// Compute the Rayleigh-scattering optical depth per Pascal of dry air
///
/// Uses the method of Bucholtz, A., 1995: Rayleigh-scattering
/// calculations for the terrestrial atmosphere. Appl. Opt., 34,
/// 2765-2773.
Real 
flotsam::rayleigh_optical_depth_per_Pa(Real wavelength ///< Wavelength (m)
				       ) {

  // Physical constants
  static const Real avogadros_number = 6.02214085774e23; // mol-1
  static const Real molar_mass_dry_air = 0.028964; // kg mol-1
  static const Real accel_due_to_grav = 9.80665; // m s-2
  // Coefficients of Bucholtz method for wavlength < 0.5 microns
  static const Real A1 = 3.01577e-28;
  static const Real B1 = 3.5512;
  static const Real C1 = 1.35579;
  static const Real D1 = 0.11563;
  // Coefficients of Bucholtz method for wavlength >= 0.5 microns
  static const Real A2 = 4.01061e-28;
  static const Real B2 = 3.99668;
  static const Real C2 = 1.10298e-3;
  static const Real D2 = 2.71393e-2;

  // Convert to microns
  Real wavelength_um = wavelength * 1.0e6;
  
  // Compute Rayleigh cross-section per molecule (m2)
  Real rayleigh_xs;
  if (wavelength_um < 0.5) {
    rayleigh_xs = 1.0e-4 * A1 * pow(wavelength_um, 
				    -(B1 + C1*wavelength_um + D1/wavelength_um));
  }
  else {
    rayleigh_xs = 1.0e-4 * A2 * pow(wavelength_um, 
				    -(B2 + C2*wavelength_um + D2/wavelength_um));
  }
  return rayleigh_xs*avogadros_number / (molar_mass_dry_air*accel_due_to_grav);
}
