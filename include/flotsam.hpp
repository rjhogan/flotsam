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

  using namespace adept;

  /// Compute the integral of the phase function over a sphere, where
  /// the input phase functions elements are assumed to be equally
  /// spaced in angle between 0 and 180 degrees.
  Real
  integrate_phase_function(const Vector& pf); ///< Phase function in

  /// Process raw phase functions [phase_func,angle] to compute the
  /// effective single-scattering phase function after the delta peak
  /// has been removed, and the amplitudes of the components of the
  /// phase-function decomposition. The phase functions elements are
  /// assumed to be equally spaced in angle between 0 and 180 degrees.
  void
  analyse_phase_functions(const Matrix& pf,       ///< Phase functions in
			  Real normalization,     ///< Integral over surface of sphere
			  Matrix& pf_out,         ///< Phase functions out
			  Matrix& pf_components); ///< Phase function components

  /// Process a single phase function to compute the effective
  /// single-scattering phase function after the delta peak has been
  /// removed, and the amplitudes of the components of the
  /// phase-function decomposition
  void
  analyse_phase_function(const Vector& pf,       ///< Phase function in
			 Real normalization,     ///< Integral over surface of sphere
			 Vector& pf_out,         ///< Phase function out
			 Vector& pf_components); ///< Component amplitudes out


  /// Compress an input grid to reduce contiguous clear-sky layers
  /// into single layers
  void
  compress_grid(const Vector& edge_pressure, ///< Half-level pressures
		const intVector& location,   ///< Locations of cloud/aerosol
		Vector& edge_pressure_out,   ///< Compressed-grid pressures
		intVector& location_out);    ///< Compressed-grid locations

  /// Get the "raw", smoothed and twice-smoothed phase functions
  /// corresponding to the phase-function component amplitudes
  /// provided in pf_components, interpolated regularly in angle from
  /// 0 to 180 degrees with "nang" elements.
  Matrix
  reconstruct_phase_function(int nang, const Vector& pf_components);

  /// Spherical convolution of phase function "pf" (linearly spaced in
  /// angle between 0 and 180 degrees) with the phase function
  /// component indicated with "icomponent".  A negative number (the
  /// default if this argument is omitted) uses the default smoothing
  /// kernel, otherwise you can use one of the basis functions used to
  /// parameterize phase functions.
  Vector
  convolve_phase_function(const Vector& pf,      ///< Phase function in
			  int icomponent = -1);  ///< Index of component

  /// Return a matrix containing the normalized basis functions used
  /// to parameterize phase functions
  Matrix basis_functions();

  //void generate_component_phase_functions(int nang);

  /// Return indices to those components used to reconstruct the phase
  /// function that integrate to 1 (as opposed to one of them that
  /// integrates to zero)
  intVector mass_component_index();
}

#endif
