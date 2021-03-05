/// @file      bi_directional_surface_reflectivity.hpp
/// @brief     directional reflectivity properties of sea water surface
/// @copyright 2017 European Centre for Medium Range Weather Forcasts
/// @license   Apache License Version 2 (see the NOTICE.md file for details)

#ifndef FLOTSAM_BI_DIRECTIONAL_REFLECTIVITY_HPP
#define FLOTSAM_BI_DIRECTIONAL_REFLECTIVITY_HPP

#include <complex>

#include "base.hpp"

namespace flotsam {

  /// Compute the surface reflectivity at the specified
  /// angle
  void white_caps_reflectance(Real wind_speed, Real wavelength, Real& wc_reflectance, Real& wc_cover);

  Real sea_reflection_coefficient_cox_and_munk(bool apply_shadowing, int slope_dist_shape,
					       std::complex<double> cn1, 
					       std::complex<double> cn2, 
					       Real dmui, Real phii, Real dmur, Real phir,
					       Real wdir, Real wspd);

  Real sea_reflection_matrix(bool apply_shadowing, std::complex<double> cn1, 
			     std::complex<double> cn2, Real wspd,
			     Real dmui, Real phii, Real dmur, Real phir);

  void ocean_brdf_lut(bool apply_shadowing, int slope_dist_shape, Real wind_spd, 
		      Real wind_dir, Real wavelength, Real pigment_conc, 
		      Real salinity, Array3D brdf, Vector brdf_zenith, Real& brdf_hemispheric);

  Real ocean_brdf(bool apply_shadowing, int slope_dist_shape, Real wind_spd, 
		  Real wind_dir, Real wavelength,Real pigment_conc, 
		  Real salinity, Real mu_sun, Real mu_inst, Real azim);

  Real ocean_water_leaving_radiance(Real wavelength, Real pigment_conc, 
				      Real dmui, Real dmur, Real wind_spd,
				      std::complex<double> cn2);

  Real interpolate_1D(Vector xi,Vector yi,Real xo);
  Real interpolate_2D(Vector xi,Vector yi,Matrix zi,Real xo,Real yo);

  Real integrate_brdf(Vector xi,Vector yi);

  void extract_from_brdf_lut(Array3D brdf, Vector brdf_zenith, Real brdf_hemispheric, 
			     Real mu_sun, Real mu_inst, 
			     Real azim, Vector glint_components);

  Vector ocean_water_refractive_index(Real wavelength, Real salinity);


  
  
}

#endif
