/// @file      OceanBRDF.hpp
/// @brief     Declares the class OceanBRDF
/// @copyright 2017 European Centre for Medium Range Weather Forcasts
/// @license   Apache License Version 2 (see the NOTICE.md file for details)

#ifndef FLOTSAM_OCEANBRDF_HPP
#define FLOTSAM_OCEANBRDF_HPP

#include <adept_arrays.h>

namespace flotsam {

  using namespace adept;

  /// Structure for holding the ocean bi-directional reflectance
  /// function look-up table at a particular wavelength
  class OceanBRDF {
  public:
    OceanBRDF() : is_initialized(false) { }
    ~OceanBRDF() { };

    /// Create look-up table 
    void create(Real wavelength_,      ///< Wavelength (m)
		const Vector& winds_,  ///< Wind speeds (m/s)
		Real pigment_conc_mg_m3 = 0.02,
		Real salinity_ppt = 35.0,
		bool apply_shadowing = true,
		int slope_dist_shape = 2);
		   
    /// Extract surface reflection vector at specified angles and
    /// surface conditions
    ///
    /// @return a vector of length 4 containing: white-sky albedo,
    /// black-sky albedo (diffuse reflection from direct source),
    /// reflected radiance for a diffuse source reflected radiance for
    /// a direct source.
    Vector extract(Real mu_sun,  ///< Cosine of solar zenith angle
		   Real mu_inst, ///< Cosine of instrument zenith angle
		   Real dazim,   ///< Azimuth between sun & instrument (radians)
		   Real wind,    ///< Wind speed (m/s)
		   Real wind_dir = 0.0 ///< Wind direction (expressed how?)
		   );

    /// Clear all look-up tables
    void clear();
    
  protected:
    Array4D brdf4D;
    Matrix  brdf_zenith2D;
    Vector  brdf_hemispheric1D;

    Vector winds;  ///< Wind speed coordinate (m/s)
    Vector dazims; ///< Sun-instrument azimuth angle coordinate (radians)
    Vector zenith_angles_deg; ///< Zenith angle coordinate (degrees)

    Real wavelength;     ///< Wavelength (m)
    Real pigment_conc_mg_m3;
    Real salinity_ppt;
    bool apply_shadowing;
    int slope_dist_shape; // shape of the glint distribution
    bool is_initialized; ///< Is the object initialized?
  };
}

#endif
