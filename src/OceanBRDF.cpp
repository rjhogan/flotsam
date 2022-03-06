/// @file      OceanBRDF.cpp
/// @brief     Defines the class OceanBRDF
/// @copyright 2017 European Centre for Medium Range Weather Forcasts
/// @license   Apache License Version 2 (see the NOTICE.md file for details)

#include "base.hpp"
#include "OceanBRDF.hpp"
#include "bi_directional_surface_reflectivity.hpp"

namespace flotsam {

  /// Create look-up table 
  void OceanBRDF::create(Real wavelength_,      ///< Wavelength (m)
			 const Vector& winds_,   ///< Wind speed (m/s)
			 Real pigment_conc_mg_m3_,
			 Real salinity_ppt_,
			 bool apply_shadowing_,
			 int slope_dist_shape_) {
    clear();

    // Copy over inputs
    wavelength = wavelength_;
    winds = winds_;
    pigment_conc_mg_m3 = pigment_conc_mg_m3_;
    salinity_ppt = salinity_ppt_;
    apply_shadowing = apply_shadowing_;
    slope_dist_shape = slope_dist_shape_;
    
    // Set coordinate variables
    zenith_angles_deg = linspace(0.0, 90.0,     91); 
    dazims            = linspace(0.0, 2.0*M_PI, 240);

    // Resize look-up tables
    brdf4D.resize(winds.size(), zenith_angles_deg.size(),
		  zenith_angles_deg.size(), dazims.size());
    brdf_zenith2D.resize(winds.size(), zenith_angles_deg.size());
    brdf_hemispheric1D.resize(winds.size());
    
    Real wind_dir = 0.0;

    // Compute look-up tables for each wind speed
    for (int i = 0; i < winds.size(); i++) {
      ocean_brdf_lut(apply_shadowing,slope_dist_shape_,
		     winds[i],wind_dir,wavelength,
		     pigment_conc_mg_m3,salinity_ppt,
		     brdf4D[i],brdf_zenith2D[i],brdf_hemispheric1D(i));
    }

    is_initialized = true;
  }

  /// Extract surface reflection vector at specified angles and
  /// surface conditions
  ///
  /// @return a vector of length 4 containing: white-sky albedo,
  /// black-sky albedo (diffuse reflection from direct source),
  /// reflected radiance for a diffuse source reflected radiance for a
  /// direct source.
  Vector OceanBRDF::extract(Real mu_sun,  ///< Cosine of solar zenith angle
			    Real mu_inst, ///< Cosine of instrument zenith angle
			    Real dazim,   ///< Azimuth between sun & instrument (radians)
			    Real wind,    ///< Wind speed (m/s)
			    Real wind_dir ///< Wind direction (expressed how?)
			    ) {

    Vector albedo(4);

    // just nearest-neighbour interpolation for wind
    int w_index = 0;
    while (w_index < winds.size()-2 && wind > winds(w_index+1)) {
      w_index++;
    }
    if (w_index+1 < winds.size()
	&& fabs(wind - winds(w_index))
	> fabs(wind-winds(w_index+1))) {
      w_index++;
    }

    extract_from_brdf_lut(brdf4D[w_index], brdf_zenith2D[w_index],
			  brdf_hemispheric1D[w_index],
			  mu_sun, mu_inst, dazim, albedo);
    albedo(3) = ocean_brdf(apply_shadowing,slope_dist_shape,
			   wind, wind_dir, wavelength,
			   pigment_conc_mg_m3, salinity_ppt,
			   mu_sun, mu_inst, dazim);
    return albedo;
  }

  /// Clear all look-up tables
  void OceanBRDF::clear() {
    brdf4D.clear();
    brdf_zenith2D.clear();
    brdf_hemispheric1D.clear();
    winds.clear();
    dazims.clear();
    zenith_angles_deg.clear();
    is_initialized = false;
  }
 
}
