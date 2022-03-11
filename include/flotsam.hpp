/** @file      flotsam.hpp
    @brief     Header for C++ features of FLOTSAM solar radiance model
    @copyright 2017 European Centre for Medium Range Weather Forcasts
    @license   Apache License Version 2 (see the NOTICE.md file for details)

 */

#ifndef FLOTSAM_HPP
#define FLOTSAM_HPP

#include <flotsam.h>

#include <adept_arrays.h>

namespace flotsam {

  /// Compute the integral of the phase function over a sphere, where
  /// the input phase functions elements are assumed to be equally
  /// spaced in angle between 0 and 180 degrees.
  adept::Real
  integrate_phase_function(const adept::Vector& pf); ///< Phase function in

  /// Process raw phase functions [phase_func,angle] to compute the
  /// effective single-scattering phase function after the delta peak
  /// has been removed, and the amplitudes of the components of the
  /// phase-function decomposition. The phase functions elements are
  /// assumed to be equally spaced in angle between 0 and 180 degrees.
  void
  analyse_phase_functions(const adept::Matrix& pf,       ///< Phase functions in
			  adept::Real normalization,     ///< Integral over surface of sphere
			  adept::Matrix pf_out,          ///< Phase functions out
			  adept::Matrix pf_components);  ///< Phase function components
  /// Process raw phase functions [phase_func,angle] to compute the
  /// effective single-scattering phase function after the delta peak
  /// has been removed, and the amplitudes of the components of the
  /// phase-function decomposition. The phase functions elements need
  /// not be evenly spaced in angle, but must include 0 and pi radians.
  void
  analyse_phase_functions(const adept::Vector& angle,    ///< Scattering angle (radians)
			  const adept::Matrix& pf,       ///< Phase functions in
			  adept::Real normalization,     ///< Integral over surface of sphere
			  adept::Matrix pf_out,          ///< Phase functions out
			  adept::Matrix pf_components);  ///< Phase function components

  /// Process a single phase function to compute the effective
  /// single-scattering phase function after the delta peak has been
  /// removed, and the amplitudes of the components of the
  /// phase-function decomposition, assuming the phase function points
  /// are evenly spaced in angle
  void
  analyse_phase_function(const adept::Vector& pf,       ///< Phase function in
			 adept::Real normalization,     ///< Integral over surface of sphere
			 adept::Vector pf_out,          ///< Phase function out
			 adept::Vector pf_components);  ///< Component amplitudes out

  /// Process a single phase function to compute the effective
  /// single-scattering phase function after the delta peak has been
  /// removed, and the amplitudes of the components of the
  /// phase-function decomposition, where the scattering angle is
  /// provided explicitly
  void
  analyse_phase_function(const adept::Vector& angle,    ///< Scattering angle (radians)
			 const adept::Vector& pf,       ///< Phase function in
			 adept::Real normalization,     ///< Integral over surface of sphere
			 adept::Vector pf_out,         ///< Phase function out
			 adept::Vector pf_components); ///< Component amplitudes out
    
  /// Compress an input grid to reduce contiguous clear-sky layers
  /// into single layers
  void
  compress_grid(const adept::Vector& edge_pressure, ///< Half-level pressures
		const adept::intVector& location,   ///< Locations of cloud/aerosol
		adept::Vector& edge_pressure_out,   ///< Compressed-grid pressures
		adept::intVector& location_out);    ///< Compressed-grid locations

  /// Get the "raw", smoothed and twice-smoothed phase functions
  /// corresponding to the phase-function component amplitudes
  /// provided in pf_components, interpolated regularly in angle from
  /// 0 to 180 degrees with "nang" elements.
  adept::Matrix
  reconstruct_phase_function(int nang, const adept::Vector& pf_components);

  /// Spherical convolution of phase function "pf" (linearly spaced in
  /// angle between 0 and 180 degrees) with the phase function
  /// component indicated with "icomponent".  A negative number (the
  /// default if this argument is omitted) uses the default smoothing
  /// kernel, otherwise you can use one of the basis functions used to
  /// parameterize phase functions.
  adept::Vector
  convolve_phase_function(const adept::Vector& pf,      ///< Phase function in
			  int icomponent = -1);  ///< Index of component

  /// Return a matrix containing the normalized basis functions used
  /// to parameterize phase functions
  adept::Matrix basis_functions();

  //void generate_component_phase_functions(int nang);

  /// Return indices to those components used to reconstruct the phase
  /// function that integrate to 1 (as opposed to one of them that
  /// integrates to zero)
  adept::intVector mass_component_index();
}

#endif
