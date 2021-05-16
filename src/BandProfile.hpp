/// @file      BandProfile.hpp
/// @brief     Declares the class BandProfile
/// @copyright 2017 European Centre for Medium Range Weather Forcasts
/// @license   Apache License Version 2 (see the NOTICE.md file for details)

#ifndef FLOTSAM_BAND_PROFILE_HPP
#define FLOTSAM_BAND_PROFILE_HPP

#include <vector>

#include "base.hpp"
#include "Channel.hpp"
#include "BackgroundProfile.hpp"

namespace flotsam {

  /// Pi minus the great circle angle (radians) between two points on a sphere
  ///
  /// In terms of latitude, the formula is
  /// acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(dlon), but note
  /// that mu is the cosine of the zenith angle so it is equal to the
  /// sine of the elevation angle.
  inline
  Real great_circle_angle(Real mu_sun,  ///< Cosine of solar zenith angle
			  Real mu_inst, ///< Cosine of instrument zenith angle
			  Real azim     ///< Azimuthal separation (radians)
			  ) {
    Real cos_sun_elev  = sqrt(1.0 - mu_sun*mu_sun);
    Real cos_inst_elev = sqrt(1.0 - mu_inst*mu_inst);
    Real arg = mu_sun*mu_inst + cos_sun_elev*cos_inst_elev*cos(azim);
    if (arg < -1.0) {
      arg = -1.0;
    }
    else if (arg > 1.0) {
      arg = 1.0;
    }
    return M_PI - acos(arg);
  }

  /// Holds the optical properties of atmospheric gases for a
  /// particular atmospheric profile
  struct BandProfile {

    BandProfile() 
      : n_z(0), n_g(0) {
      //      component_g.resize(FLOTSAM_PHASE_FUNC_MAX_COMPONENTS);
      //      component_g << 0.0, 0.0, 0.834749537379482, 0.834749537379482;
    }
  
    /// Calculate gas optical properties for specified channel and
    /// background profile
    int init(const Channel& chan, const BackgroundProfile& prof);

    /// Direct initialization avoiding if the Rayleigh and absorption
    /// profiles are to be computed outside FLOTSAM, in which case the
    /// Channel and BackgroundProfile types are not needed
    int init_direct(const Vector& weight_,
		    const Matrix& od_rayleigh_,
		    const Matrix& od_gas_abs_ = Matrix()); // Default zero absorption

    /// Calculate angles and look-up tables that depend on them
    int set_geometry(Real mu_sun_, Real mu_inst_, Real azim_);

    /// Override Rayleigh optical depths in each layer
    int set_rayleigh_optical_depth(const Vector& od);

    /// Override gas absorption optical depth in each layer
    int set_gas_absorption_optical_depth(const Vector& od);

    /// Return the Rayleigh optical depth for the whole atmosphere,
    /// averaged over spectral subdivisions if there are any
    Real get_total_rayleigh_optical_depth();

    /// Return the gas-absorption optical depth for the whole
    /// atmosphere, averaged over spectral subdivisions if there are
    /// any
    Real get_total_gas_absorption_optical_depth();

    // Data

    // Atmospheric optical depth and weights
    Vector weight; ///< Weight of each g-point
    
    //    Vector component_g; ///< Asymmetry factor of each particulate phase function component

    Matrix od_gas_abs;  ///< Gas absorption optical depth
    Matrix od_rayleigh; ///< Optical depth from Rayleigh scattering only

    /// Phase function components
    ///
    /// First dimension is sun beam: 0=direct,1=lobe-dn,2=lobe-up;
    /// second dimension is instrument beam:
    /// 0=direct,1=lobe-dn,2=lobe-up; third dimension is phase
    /// function component
    Array3D pf_components;

    /// Phase function coefficients
    ///
    /// Polynomial coefficients to reconstruct the smoothing
    /// dependence of phase components. First dimension is sun beam:
    /// 0=direct,1=lobe-dn,2=lobe-up; second dimension is instrument
    /// beam: 0=direct,1=lobe-dn,2=lobe-up; third dimension is phase
    /// function component; fourth component is polynomial
    /// coefficient.
    Array4D pf_coeffs;

    Real mu_sun, mu_sun_lobe_dn, mu_sun_lobe_up, frac_sun_lobe_dn;
    Real sec_sun, sec_sun_lobe_dn, sec_sun_lobe_up;
    Real mu_inst, mu_inst_lobe_dn, mu_inst_lobe_up, frac_inst_lobe_dn;
    Real sec_inst, sec_inst_lobe_dn, sec_inst_lobe_up;
    Real azim;    ///< Azimuth angle between sun and instrument
    Matrix angle; ///< Great circle angle between sun (direct,lobe-dn,lobe-up) and instrument (same)

    Real pf_rayleigh; ///< Value of Rayleigh phase function at scattering angle
    Real max_angular_variance; ///< Maximum angular variance in look-up table
    int n_z; ///< Number of layers (whether or not they contain particulates)
    int n_g; ///< Number of g points
    bool do_sun_lobe_up; ///< Do we simulate upwelling solar lobe?
    bool do_inst_lobe_up;///< Do we simulate equivalent for return lobe?
  };
}

#endif
