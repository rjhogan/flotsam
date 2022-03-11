/// @file      c_interface.cpp
/// @brief     Implementation of C interface for FLOTSAM solar radiance model
/// @copyright 2016- European Centre for Medium Range Weather Forcasts
/// @license   Apache License Version 2 (see the NOTICE.md file for details)

#include <vector>
#include <map>

#include <flotsam.h>
#include <flotsam.hpp>

#include <adept_arrays.h>

#include "Channel.hpp"
#include "BackgroundProfile.hpp"
#include "BandProfile.hpp"
#include "OceanBRDF.hpp"
#include "LookUpTable.hpp"
#include "reflectance.hpp"

#define CATCH_CPP_ERRORS

using namespace flotsam;

static std::map<int,Channel>           channel_list;
static std::map<int,BackgroundProfile> background_profile_list;
static std::map<int,BandProfile>       band_profile_list;
static std::map<int,OceanBRDF>         ocean_brdf_list;

static const char* error_messages[] 
= {"Error thrown by C++ Adept library",
   "Incorrect number of spectral intervals",
   "Invalid number of albedo codes (must be 1 or 4)",
   "Error merging gases",
   "Automatic differentiation error",
   "Input profiles must start at top-of-atmosphere",
   "A gas required for the channel was not provided",
   "Attempted to create a look-up table with more than 3 gases",
   "Feature not available",
   "Number of layers not consistent between contexts",
   "Invalid channel, background profile or band profile context",
   "Invalid gas code (must be in range 0--4)",
   "Input variable outside physical bounds",
   "Error reading channel configuration file",
   "Channel configuration file not found",
   "Memory error"};

/// @name Support functions
/// @{ 

/// Return the next free integer ID from a map
template <typename T>
int next_free_id(const std::map<int,T>& m) {
  typename std::map<int,T>::const_iterator it;
  int i = 0;
  while (m.count(i) > 0) {
    ++i;
  }
  return i;
}

// All flotsam_* functions must be callable from C (no C++ name
// mangling)
extern "C" {

/// Return a pointer to a NULL-terminated string explaining an error
/// code
const char* flotsam_error_message(int error_code) {
  if (error_code >= 0) {
    return "Success";
  }
  else if (error_code >= FLOTSAM_FIRST_ERROR_CODE) {
    return error_messages[error_code-FLOTSAM_FIRST_ERROR_CODE];
  }
  else {
    return "Unknown error";
  }
}

/// Return the name of the gas ("H2O" etc) or NULL if not recognized
const char* flotsam_gas_name(flotsam_gas_t igas) {
  static const char* gas_names[] = {"H2O", "CO2", "O3", "CH4", "N2O"};
  if (igas >= 0 && igas < FLOTSAM_MAX_GASES) {
    return gas_names[igas];
  }
  else {
    return 0;
  }
}

/// Return the value of the phase function at the specified scattering
/// angle, assuming phase function elements to be equally spaced in
/// angle including the values at 0 and 180 degrees
flotsam_real flotsam_interp_phase_func(int n, const flotsam_real* pf,
				       flotsam_real ang) {
  flotsam_real real_index = ang * (n-1) / M_PI;
  int index = std::max(0,std::min(static_cast<int>(real_index),n-2));
  flotsam_real weight = real_index - index;
  return (1.0 - weight) * pf[index] + weight * pf[index+1];
}

/// Pi minus the great circle angle (radians) between two points on a sphere
///
/// In terms of latitude, the formula is
/// acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(dlon), but note
/// that mu is the cosine of the zenith angle so it is equal to the
/// sine of the elevation angle.
flotsam_real flotsam_scattering_angle(flotsam_real mu_sun,  ///< Cosine of solar zenith angle
				      flotsam_real mu_inst, ///< Cosine of instrument zenith angle
				      flotsam_real azim     ///< Azimuthal separation (radians)
				      ) {
  flotsam_real cos_sun_elev  = sqrt(1.0 - mu_sun*mu_sun);
  flotsam_real cos_inst_elev = sqrt(1.0 - mu_inst*mu_inst);
  flotsam_real arg = mu_sun*mu_inst + cos_sun_elev*cos_inst_elev*cos(azim);
  if (arg < -1.0) {
    arg = -1.0;
  }
  else if (arg > 1.0) {
    arg = 1.0;
  }
  return M_PI - acos(arg);
}

/// Return the number of components used to parameterize phase
/// functions (excluding the delta component)
int flotsam_n_phase_function_components() {
  return lut.n_phase_function_components();
}

/// @}

/// @name Channel functions
/// @{

int flotsam_new_channel(const char* filename) {
  int id = next_free_id(channel_list);
  channel_list[id] = Channel();
  int status = channel_list[id].read(filename);
  if (status < 0) {
    return status;
  }
  else {
    return id;
  }
}

int flotsam_new_channel_vacuum() {
  int id = next_free_id(channel_list);
  channel_list[id] = Channel();
  int status = channel_list[id].vacuum();
  if (status < 0) {
    return status;
  }
  else {
    return id;
  }
}

int flotsam_new_channel_rayleigh_only(flotsam_real wavelength) {
  int id = next_free_id(channel_list);
  channel_list[id] = Channel();
  int status = channel_list[id].rayleigh_only(wavelength);
  if (status < 0) {
    return status;
  }
  else {
    return id;
  }
}

/// Equivalent of flotsam_new_channel but for Fortran when the string
/// length is passed separately
int flotsam_new_channel_(const char* filename, long int len) {
  int id = next_free_id(channel_list);
  channel_list[id] = Channel();
  char filename_terminated[len+1]; // VLA extension
  for (int i = 0; i < len; ++i) {
    filename_terminated[i] = filename[i];
  }
  filename_terminated[len] = '\0';
  int status = channel_list[id].read(filename_terminated);
  if (status < 0) {
    return status;
  }
  else {
    return id;
  }
}

int flotsam_free_channel(int ichan) {
  if (channel_list.count(ichan) > 0) {
    channel_list.erase(ichan);
    return FLOTSAM_SUCCESS;
  }
  else {
    return FLOTSAM_INVALID_CONTEXT;
  }
}

int flotsam_get_wavelength(int ichan, flotsam_real* wavelength) {
  if (channel_list.count(ichan) > 0) {
    *wavelength = channel_list[ichan].wavelength;
    //    *wavelength_range = channel_list[ichan].wavelength_range;
    return FLOTSAM_SUCCESS;
  }
  else {
    return FLOTSAM_INVALID_CONTEXT;
  }  
}

/*
int flotsam_set_wavelength_range(int ichan, flotsam_real wavelength_range) {
  if (channel_list.count(ichan) > 0) {
    channel_list[ichan].wavelength_range = wavelength_range;
    return FLOTSAM_SUCCESS;
  }
  else {
    return FLOTSAM_INVALID_CONTEXT;
  }  
}
*/

/// @}

/// @name Background profile functions
/// @{

int flotsam_new_background_profile() {
  int id = next_free_id(background_profile_list);
  background_profile_list[id] = BackgroundProfile();
  return id;
}

int flotsam_free_background_profile(int iprof) {
  if (background_profile_list.count(iprof) > 0) {
    background_profile_list.erase(iprof);
    return FLOTSAM_SUCCESS;
  }
  else {
    return FLOTSAM_INVALID_CONTEXT;
  }  
}

int flotsam_reset_background_profile(int iprof, int n) {
  if (background_profile_list.count(iprof) > 0) {
    return background_profile_list[iprof].reset(n);
  }
  else {
    return FLOTSAM_INVALID_CONTEXT;
  }
}

int flotsam_set_edge_pressure(int iprof,
			      int n,
			      const flotsam_real* edge_pressure) {
  if (background_profile_list.count(iprof) > 0) {
    return background_profile_list[iprof].set_edge_pressure(n, edge_pressure);
  }
  else {
    return FLOTSAM_INVALID_CONTEXT;
  }
}

int flotsam_set_temperature(int iprof,
			    int n,
			    const flotsam_real* temperature) {
  if (background_profile_list.count(iprof) > 0) {
    return background_profile_list[iprof].set_temperature(n, temperature);
  }
  else {
    return FLOTSAM_INVALID_CONTEXT;
  }
}

int flotsam_set_gas_const(int iprof, flotsam_gas_t igas, flotsam_real mmr) {
  if (background_profile_list.count(iprof) > 0) {
    return background_profile_list[iprof].set_gas_const(igas, mmr);
  }
  else {
    return FLOTSAM_INVALID_CONTEXT;
  }
}

int flotsam_set_gas_profile(int iprof, flotsam_gas_t igas, int n,
			  const flotsam_real* mmr) {
  if (background_profile_list.count(iprof) > 0) {
    return background_profile_list[iprof].set_gas_profile(igas, n, mmr);
  }
  else {
    return FLOTSAM_INVALID_CONTEXT;
  }
}

/// @}

/// @name Band profile functions
/// @{

int flotsam_new_band_profile() {
  int id;
#pragma omp critical
  {
    id = next_free_id(band_profile_list);
    band_profile_list[id] = BandProfile();
  }
  return id;
}

int flotsam_free_band_profile(int iband) {
  if (band_profile_list.count(iband) > 0) {
#pragma omp critical
    band_profile_list.erase(iband);
    return FLOTSAM_SUCCESS;
  }
  else {
    return FLOTSAM_INVALID_CONTEXT;
  }  
}

int flotsam_init_band_profile(int iband, int ichan, int iprof) {
  if (band_profile_list.count(iband) == 0
      || channel_list.count(ichan) == 0
      || background_profile_list.count(iprof) == 0) {
    return FLOTSAM_INVALID_CONTEXT;
  }
  else {
    return band_profile_list[iband].init(channel_list[ichan],
					 background_profile_list[iprof]);
  }
}

int flotsam_init_band_profile_direct(int iband, int n_g, int n_z,
				     const flotsam_real* weight,
				     const flotsam_real* od_rayleigh,
				     const flotsam_real* od_gas_abs) {
  if (band_profile_list.count(iband) == 0) {
    return FLOTSAM_INVALID_CONTEXT;
  }
  else {
    const Vector weight_(const_cast<flotsam_real*>(weight), dimensions(n_g));
    const Matrix od_rayleigh_(const_cast<flotsam_real*>(od_rayleigh), dimensions(n_g,n_z));
    if (od_gas_abs) {
      const Matrix od_gas_abs_(const_cast<flotsam_real*>(od_gas_abs), dimensions(n_g,n_z));
      return band_profile_list[iband].init_direct(weight_, od_rayleigh_, od_gas_abs_);
    }
    else {
      return band_profile_list[iband].init_direct(weight_, od_rayleigh_);
    }
  } 
}

/*
int flotsam_set_phase_func_components(int iband,
				      int n_phase_func_components,
				      const flotsam_phase_func_component_t* pfc_list) {
  if (band_profile_list.count(iband) > 0) {
    return band_profile_list[iband].set_phase_func_components(n_phase_func_components,
							      pfc_list);
  }
  else {
    return FLOTSAM_INVALID_CONTEXT;
  }
}
*/

int flotsam_set_geometry(int iband,            /**< Band profile ID */
			 flotsam_real mu_sun,  /**< Zenith angle of sun */
			 flotsam_real mu_inst, /**< Zenith angle of instrument */
			 flotsam_real azim     /**< Azimuth angle between sun & instrument */
			 ) {
  if (band_profile_list.count(iband) > 0) {
    return band_profile_list[iband].set_geometry(mu_sun, mu_inst, azim);
  }
  else {
    return FLOTSAM_INVALID_CONTEXT;
  }
}

int flotsam_set_rayleigh_optical_depth(int iband,                /**< Band profile ID */
				       int n,                    /**< Number of layers */
				       const flotsam_real* od) { /**< Layer optical depth */
  if (band_profile_list.count(iband) > 0) {
    return band_profile_list[iband].set_rayleigh_optical_depth(
		       Vector(const_cast<flotsam_real*>(od),ExpressionSize<1>(n)));
  }
  else {
    return FLOTSAM_INVALID_CONTEXT;
  }
}

int flotsam_set_gas_absorption_optical_depth(int iband,                /**< Band profile ID */
					     int n,                    /**< Number of layers */
					     const flotsam_real* od) { /**< Layer optical depth */
  if (band_profile_list.count(iband) > 0) {
    return band_profile_list[iband].set_gas_absorption_optical_depth(
			     Vector(const_cast<flotsam_real*>(od),ExpressionSize<1>(n)));
  }
  else {
    return FLOTSAM_INVALID_CONTEXT;
  }
}

int flotsam_get_total_rayleigh_optical_depth(int iband,          /**< Band profile ID */
					     flotsam_real* od) { /**< Returned optical depth */
  if (band_profile_list.count(iband) > 0) {
    *od = band_profile_list[iband].get_total_rayleigh_optical_depth();
    return FLOTSAM_SUCCESS;	
  }
  else {
    return FLOTSAM_INVALID_CONTEXT;
  }
}

int flotsam_get_total_gas_absorption_optical_depth(int iband,          /**< Band profile ID */
						   flotsam_real* od) { /**< Returned optical depth */
  if (band_profile_list.count(iband) > 0) {
    *od = band_profile_list[iband].get_total_gas_absorption_optical_depth();
    return FLOTSAM_SUCCESS;	
  }
  else {
    return FLOTSAM_INVALID_CONTEXT;
  }
}


/// @}


/// @name Ocean BRDF functions

int flotsam_new_ocean_brdf(flotsam_real wavelength,    /**< Wavelength (m) */
			   int nwinds,                 /**< Number of wind speeds */
			   const flotsam_real* winds   /**< Wind speeds (m/2) */
			   ) {
  int id;
  #pragma omp critical
  {
    id = next_free_id(ocean_brdf_list);
    ocean_brdf_list[id] = OceanBRDF();
    ocean_brdf_list[id].create(wavelength, 
			       Vector(const_cast<flotsam_real*>(winds),ExpressionSize<1>(nwinds)));
  }
  return id;
}

int flotsam_new_ocean_brdf_detailed(flotsam_real wavelength,    /**< Wavelength (m) */
				    int nwinds,                 /**< Number of wind speeds */
				    const flotsam_real* winds,  /**< Wind speeds (m/2) */
				    flotsam_real pigment_conc_mg_m3,
				    flotsam_real salinity_ppt,
				    int apply_shadowing,
				    int slope_dist_shape) {
  int id = next_free_id(ocean_brdf_list);
  ocean_brdf_list[id] = OceanBRDF();
  ocean_brdf_list[id].create(wavelength, 
			     Vector(const_cast<flotsam_real*>(winds),ExpressionSize<1>(nwinds)),
			     pigment_conc_mg_m3, salinity_ppt, apply_shadowing, slope_dist_shape);

  return id;
}

int flotsam_get_ocean_albedo_components(
	  int ibrdf,              /**< Ocean BRDF ID */ 
	  flotsam_real mu_sun,    /**< Zenith angle of sun */
	  flotsam_real mu_inst,   /**< Zenith angle of instrument */
	  flotsam_real azim,      /**< Azimuth angle between sun & instrument (radians) */
	  flotsam_real wind,      /**< Wind speed (m/s) */
	  flotsam_real albedos[FLOTSAM_NUM_ALBEDO_COMPONENTS] /**< Returned albedo components */
					) {
  if (ocean_brdf_list.count(ibrdf) > 0) {
    Vector a(albedos, FLOTSAM_NUM_ALBEDO_COMPONENTS);
    a = ocean_brdf_list[ibrdf].extract(mu_sun, mu_inst, azim, wind);
    return FLOTSAM_SUCCESS;
  }
  else {
    return FLOTSAM_INVALID_CONTEXT;
  }
}

int flotsam_free_ocean_brdf(int ibrdf) {
  if (ocean_brdf_list.count(ibrdf) > 0) {
#pragma omp critical
    ocean_brdf_list.erase(ibrdf);
    return FLOTSAM_SUCCESS;
  }
  else {
    return FLOTSAM_INVALID_CONTEXT;
  }  
}


  /** @} */



/// @name Reflectance functions

int flotsam_reflectance(int iband,               /**< Band profile ID */
			int nalbedo,             /**< Number of albedo components */
			const flotsam_real* albedo, /**< Surface albedo components */
			int n,                   /**< Number of particulate layers */
			const int* loc,          /**< Location of particulate layers */
			const flotsam_real* od,  /**< Optical depth */
			const flotsam_real* ssa, /**< Single-scattering albedo */
			const flotsam_real* pf,  /**< Phase func for single scattering */
			const flotsam_real* pfc, /**< Phase function components */
			flotsam_real* reflectance/**< Returned reflectance */
			) {
  if (band_profile_list.count(iband) > 0) {
    const Vector my_albedo(const_cast<flotsam_real*>(albedo), ExpressionSize<1>(nalbedo));
    const intVector my_loc(const_cast<int*>(loc), ExpressionSize<1>(n));
    const Vector my_od(const_cast<flotsam_real*>(od), ExpressionSize<1>(n));
    const Vector my_ssa(const_cast<flotsam_real*>(ssa), ExpressionSize<1>(n));
    const Vector my_pf(const_cast<flotsam_real*>(pf), ExpressionSize<1>(n));
    const Matrix my_pfc(const_cast<flotsam_real*>(pfc),
			ExpressionSize<2>(n,lut.n_phase_function_components()));
    return flotsam::reflectance<false>(band_profile_list[iband],
				       my_albedo, my_loc, my_od, my_ssa, my_pf, my_pfc,
				       *reflectance);
  }
  else {
    return FLOTSAM_INVALID_CONTEXT;
  }
}

int flotsam_reflectance_components(int iband,    /**< Band profile ID */
			int nalbedo,             /**< Number of albedo components */
			const flotsam_real* albedo, /**< Surface albedo components */
			int n,                   /**< Number of particulate layers */
			const int* loc,          /**< Location of particulate layers */
			const flotsam_real* od,  /**< Optical depth */
			const flotsam_real* ssa, /**< Single-scattering albedo */
			const flotsam_real* pf,  /**< Phase func for single scattering */
			const flotsam_real* pfc, /**< Phase function components */
			flotsam_real data[FLOTSAM_NUM_COMPONENTS] /**< Returned data */
				   ) {
  if (band_profile_list.count(iband) > 0) {
    const Vector my_albedo(const_cast<flotsam_real*>(albedo), ExpressionSize<1>(nalbedo));
    const intVector my_loc(const_cast<int*>(loc), ExpressionSize<1>(n));
    const Vector my_od(const_cast<flotsam_real*>(od), ExpressionSize<1>(n));
    const Vector my_ssa(const_cast<flotsam_real*>(ssa), ExpressionSize<1>(n));
    const Vector my_pf(const_cast<flotsam_real*>(pf), ExpressionSize<1>(n));
    const Matrix my_pfc(const_cast<flotsam_real*>(pfc),
			ExpressionSize<2>(n,lut.n_phase_function_components()));
    WorkingData working_data;
    int status = flotsam::reflectance<false>(band_profile_list[iband],
				       my_albedo, my_loc, my_od, my_ssa, my_pf, my_pfc,
				       data[0], &working_data);
    data[1] = working_data.ref_direct;
    data[2] = working_data.ref_lobe;
    data[3] = working_data.ref_lobe_up;
    data[4] = working_data.ref_diffuse;
    data[5] = working_data.ref_surface;
    data[6] = working_data.flux_toa;
    return status;
  }
  else {
    return FLOTSAM_INVALID_CONTEXT;
  }
}

// Fortran interface with 1-based indexing for "loc"
int flotsam_reflectance_(int iband,               /**< Band profile ID */
			 int nalbedo,             /**< Number of albedo components */
			 const flotsam_real* albedo, /**< Surface albedo components */
			 int n,                   /**< Number of particulate layers */
			 const int* loc_fortran,  /**< Location of particulate layers (1-based) */
			 const flotsam_real* od,  /**< Optical depth */
			 const flotsam_real* ssa, /**< Single-scattering albedo */
			 const flotsam_real* pf,  /**< Phase func for single scattering */
			 const flotsam_real* pfc, /**< Phase function components */
			 flotsam_real* reflectance/**< Returned reflectance */
			 ) {
  std::vector<int> loc(n);
  for (int i = 0; i < n; ++i) {
    loc[i] = loc_fortran[i]-1;
  }
  return flotsam_reflectance(iband, nalbedo, albedo, n, &loc[0], od, ssa, pf, pfc,
			     reflectance);
}

int flotsam_reflectance_unindexed(int iband,     /**< Band profile ID */
			int nalbedo,             /**< Number of albedo components */
			const flotsam_real* albedo, /**< Surface albedo components */
			int n,                   /**< Number of particulate layers */
			const flotsam_real* od,  /**< Optical depth */
			const flotsam_real* ssa, /**< Single-scattering albedo */
			const flotsam_real* pf,  /**< Phase func for single scattering */
			const flotsam_real* pfc, /**< Phase function components */
			flotsam_real* reflectance/**< Returned reflectance */
			) {
  std::vector<int> loc(n);
  for (int i = 0; i < n; ++i) {
    loc[i] = i;
  }
  return flotsam_reflectance(iband, nalbedo, albedo, n, &loc[0], od, ssa, pf, pfc,
			     reflectance);
}

int flotsam_reflectance_jacobian(int iband,               /**< Band profile ID */
				 int nalbedo,             /**< Number of albedo components */
				 const flotsam_real* albedo, /**< Surface albedo components */
				 int n,                   /**< Number of particulate layers */
				 const int* loc,          /**< Location of particulate layers */
				 const flotsam_real* od,  /**< Optical depth */
				 const flotsam_real* ssa, /**< Single-scattering albedo */
				 const flotsam_real* pf,  /**< Phase func for single scattering */
				 const flotsam_real* pfc, /**< Phase function components */
				 flotsam_real* reflectance,/**< Returned reflectance */
				 flotsam_real* d_albedo,  /**< d(reflectance)/d(albedo) */
				 flotsam_real* d_od,      /**< d(reflectance)/d(od) */
				 flotsam_real* d_ssa,     /**< d(reflectance)/d(ssa) */
				 flotsam_real* d_pf,      /**< d(reflectance)/d(pf) */
				 flotsam_real* d_pfc      /**< d(reflectance)/d(pfc) */
				 ) {

  if (band_profile_list.count(iband) > 0) {
    using namespace adept;

    adept::Stack* old_stack = adept::active_stack();
    if (old_stack) {
      old_stack->deactivate();
    }
    const intVector my_loc(const_cast<int*>(loc), ExpressionSize<1>(n));

    int status = FLOTSAM_SUCCESS;

#ifdef CATCH_CPP_ERRORS
    try 
#endif
      { // Need opening bracket so that local stack is deactivated
	// when it goes out of scope
      adept::Stack stack;
      // For speed we preallocate the likely space needed
      stack.preallocate_statements(250*n);
      stack.preallocate_operations(800*n);

      BandProfile& band = band_profile_list[iband];
      int n_pfc = lut.n_phase_function_components();

      aVector a_albedo(nalbedo);
      for (int i = 0; i < nalbedo; ++i) {
	a_albedo[i] = albedo[i];
      }
      aVector a_od(n), a_ssa(n), a_pf(n);
      aMatrix a_pfc(n, n_pfc);
      for (int i = 0; i < n; ++i) {
	a_od[i]  = od[i];
	a_ssa[i] = ssa[i];
	a_pf[i]  = pf[i];
	for (int j = 0; j < n_pfc; ++j) {
	  int index = i*n_pfc + j;
	  a_pfc(i,j) = pfc[index];
	}
      }
      aReal a_reflectance;
      
      stack.new_recording();

      // The "aReal" template argument ensures that the active version
      // of the function (with automatic differentiation) is called
      status = flotsam::reflectance<true>(band,
					  a_albedo, my_loc, a_od, a_ssa,
					  a_pf, a_pfc,
					  a_reflectance);
      if (status != FLOTSAM_SUCCESS) {
	if (old_stack) {
	  stack.deactivate();
	  old_stack->activate();
	}
	return status;
      }

      *reflectance = value(a_reflectance);

      a_reflectance.set_gradient(1.0);
      stack.reverse();
      if (d_albedo) {
	// Vector pointing to external data
	Vector jac(d_albedo, ExpressionSize<1>(nalbedo));
	a_albedo.get_gradient(jac);
      }
      if (d_od) {
	Vector jac(d_od, ExpressionSize<1>(n));
	a_od.get_gradient(jac);
      }
      if (d_ssa) {
	Vector jac(d_ssa, ExpressionSize<1>(n));
	a_ssa.get_gradient(jac);
      }
      if (d_pf) {
	Vector jac(d_pf, ExpressionSize<1>(n));
	a_pf.get_gradient(jac);
      }
      if (d_pfc) {
	Matrix jac(d_pfc, ExpressionSize<2>(n,n_pfc));
	a_pfc.get_gradient(jac);
      }
    }
#ifdef CATCH_CPP_ERRORS
    catch (const std::exception& e) {
      std::cerr << "An error occurred: " << e.what() << "\n";
      status = FLOTSAM_AUTOMATIC_DIFFERENTIATION_ERROR;
    }
    catch (...) {
      std::cerr << "An error occurred\n";
      status = FLOTSAM_AUTOMATIC_DIFFERENTIATION_ERROR;
    }
#endif

    if (old_stack) {
      old_stack->activate();
    } 
    return status;

  }
  else {
    return FLOTSAM_INVALID_CONTEXT;
  }

}


/*
int flotsam_reflectance_jacobian_thread_unsafe(int iband,
			   flotsam_real mu_instrument,
			   int n,
			   const int* loc,
			   const flotsam_real* od,
			   const flotsam_real* ssa,
			   const flotsam_real* g,
			   flotsam_real* reflectance,
			   flotsam_real* d_reflectance_d_od,
			   flotsam_real* d_reflectance_d_ssa,
			   flotsam_real* d_reflectance_d_g) {
#ifndef HAVE_ADEPT
  return FLOTSAM_FEATURE_NOT_AVAILABLE;
#else
  static adept::Stack stack(false);

  if (band_profile_list.count(iband) > 0) {
    using adept::aReal;

    adept::Stack* old_stack = adept::active_stack();
    if (old_stack) {
      old_stack->deactivate();
    }
    stack.activate();

    int status = FLOTSAM_SUCCESS;
    try {
      std::vector<aReal> a_od(n), a_ssa(n), a_g(n);
      for (int i = 0; i < n; ++i) {
	a_od[i] = od[i];
	a_ssa[i] = ssa[i];
	a_g[i] = g[i];
      }
      aReal a_rad;
      
      stack.new_recording();

      status = flotsam::reflectance<aReal>(band_profile_list[iband],
					mu_instrument, n, loc, 
					&a_od[0], &a_ssa[0], &a_g[0],
					&a_rad);
      if (status != FLOTSAM_SUCCESS) {
	if (old_stack) {
	  stack.deactivate();
	  old_stack->activate();
	}
	return status;
      }

      *reflectance = value(a_rad);

      a_rad.set_gradient(1.0);
      stack.reverse();
      adept::get_gradients(&a_od[0], n, d_reflectance_d_od);
      adept::get_gradients(&a_ssa[0], n, d_reflectance_d_ssa);
      adept::get_gradients(&a_g[0], n, d_reflectance_d_g);
    }
    catch (const std::exception& e) {
      std::cerr << "An error occurred: " << e.what() << "\n";
      status = FLOTSAM_AUTOMATIC_DIFFERENTIATION_ERROR;
    }
    catch (...) {
      std::cerr << "An error occurred\n";
      status = FLOTSAM_AUTOMATIC_DIFFERENTIATION_ERROR;
    }

    stack.deactivate();
    if (old_stack) {
      old_stack->activate();
    } 
    return status;

  }
  else {
    return FLOTSAM_INVALID_CONTEXT;
  }
#endif

}
*/



int flotsam_reflectance_jacobian_thread_unsafe(int iband, /**< Band profile ID */
				 int nalbedo,             /**< Number of albedo components */
				 const flotsam_real* albedo, /**< Surface albedo components */
				 int n,                   /**< Number of particulate layers */
				 const int* loc,          /**< Location of particulate layers */
				 const flotsam_real* od,  /**< Optical depth */
				 const flotsam_real* ssa, /**< Single-scattering albedo */
				 const flotsam_real* pf,  /**< Phase func for single scattering */
				 const flotsam_real* pfc, /**< Phase function components */
				 flotsam_real* reflectance,/**< Returned reflectance */
				 flotsam_real* d_albedo,  /**< d(reflectance)/d(albedo) */
				 flotsam_real* d_od,      /**< d(reflectance)/d(od) */
				 flotsam_real* d_ssa,     /**< d(reflectance)/d(ssa) */
				 flotsam_real* d_pf,      /**< d(reflectance)/d(pf) */
				 flotsam_real* d_pfc      /**< d(reflectance)/d(pfc) */
				 ) {

  // By making this "static", the memory allocated is reusable the
  // next time the function is called... but the function is then not
  // thread-safe
  static adept::Stack stack(false);

  if (band_profile_list.count(iband) > 0) {
    using namespace adept;

    adept::Stack* old_stack = adept::active_stack();
    if (old_stack) {
      old_stack->deactivate();
    }
    const intVector my_loc(const_cast<int*>(loc), ExpressionSize<1>(n));

    int status = FLOTSAM_SUCCESS;

#ifdef CATCH_CPP_ERRORS
    try {
#endif
      // Curly bracket as we need all active variables to go out of
      // scope before deactivating the stack
      stack.activate(); { 

      BandProfile& band = band_profile_list[iband];
      int n_pfc = lut.n_phase_function_components();

      aVector a_albedo(nalbedo);
      for (int i = 0; i < nalbedo; ++i) {
	a_albedo[i] = albedo[i];
      }
      aVector a_od(n), a_ssa(n), a_pf(n);
      aMatrix a_pfc(n, n_pfc);
      for (int i = 0; i < n; ++i) {
	a_od[i]  = od[i];
	a_ssa[i] = ssa[i];
	a_pf[i]  = pf[i];
	for (int j = 0; j < n_pfc; ++j) {
	  int index = i*n_pfc + j;
	  a_pfc(i,j) = pfc[index];
	}
      }
      aReal a_reflectance;
      
      stack.new_recording();

      // The "aReal" template argument ensures that the active version
      // of the function (with automatic differentiation) is called
      status = flotsam::reflectance<true>(band,
					  a_albedo, my_loc, a_od, a_ssa,
					  a_pf, a_pfc,
					  a_reflectance);
      if (status != FLOTSAM_SUCCESS) {
	if (old_stack) {
	  stack.deactivate();
	  old_stack->activate();
	}
	return status;
      }

      *reflectance = value(a_reflectance);

      a_reflectance.set_gradient(1.0);
      stack.reverse();
      if (d_albedo) {
	// Vector pointing to external data
	Vector jac(d_albedo, ExpressionSize<1>(nalbedo));
	a_albedo.get_gradient(jac);
      }
      if (d_od) {
	Vector jac(d_od, ExpressionSize<1>(n));
	a_od.get_gradient(jac);
      }
      if (d_ssa) {
	Vector jac(d_ssa, ExpressionSize<1>(n));
	a_ssa.get_gradient(jac);
      }
      if (d_pf) {
	Vector jac(d_pf, ExpressionSize<1>(n));
	a_pf.get_gradient(jac);
      }
      if (d_pfc) {
	Matrix jac(d_pfc, ExpressionSize<2>(n,n_pfc));
	// FIX the following does not work correctly in Adept 1.9.8:
	//a_pfc.get_gradient(jac);
	// ... so we do it by row
	for (int i = 0; i < n; ++i) {
	  jac(i,__) = a_pfc(i,__).get_gradient();
	}
      }

      } stack.deactivate();
      //      std::cerr << "STACK INFO " << stack << "\n";
#ifdef CATCH_CPP_ERRORS
    }
    catch (const std::exception& e) {
      std::cerr << "An error occurred: " << e.what() << "\n";
      status = FLOTSAM_AUTOMATIC_DIFFERENTIATION_ERROR;
    }
    catch (...) {
      std::cerr << "An error occurred\n";
      status = FLOTSAM_AUTOMATIC_DIFFERENTIATION_ERROR;
    }
#endif

    if (old_stack) {
      old_stack->activate();
    } 
    return status;

  }
  else {
    return FLOTSAM_INVALID_CONTEXT;
  }

}


int flotsam_reflectance_write_internals(int iband,  /**< Band profile ID */
		     int nalbedo,                /**< Number of albedo components */
		     const flotsam_real* albedo, /**< Surface albedo components */
		     int n,                      /**< Number of particulate layers */
		     const int* loc,             /**< Location of particulate layers */
		     const flotsam_real* od,     /**< Optical depth */
		     const flotsam_real* ssa,    /**< Single-scattering albedo */
		     const flotsam_real* pf,     /**< Phase func for single scattering */
		     const flotsam_real* pfc,    /**< Phase function components */
		     flotsam_real data[FLOTSAM_NUM_COMPONENTS], /**< Returned data */
		     const char* file_name     /**< Name of file to write */
					) {
  if (band_profile_list.count(iband) > 0) {
    const Vector my_albedo(const_cast<flotsam_real*>(albedo), ExpressionSize<1>(nalbedo));
    const intVector my_loc(const_cast<int*>(loc), ExpressionSize<1>(n));
    const Vector my_od(const_cast<flotsam_real*>(od), ExpressionSize<1>(n));
    const Vector my_ssa(const_cast<flotsam_real*>(ssa), ExpressionSize<1>(n));
    const Vector my_pf(const_cast<flotsam_real*>(pf), ExpressionSize<1>(n));
    const Matrix my_pfc(const_cast<flotsam_real*>(pfc),
			ExpressionSize<2>(n,lut.n_phase_function_components()));
    WorkingData working_data;
    working_data.output_file_name = file_name;
    return flotsam::reflectance<false>(band_profile_list[iband],
				       my_albedo, my_loc, my_od, my_ssa, my_pf, my_pfc,
				       data[0], &working_data);
  }
  else {
    return FLOTSAM_INVALID_CONTEXT;
  }
}
 
int flotsam_analyse_phase_functions(int npf,  /**< Number of phase functions */
	     int nang, /**< Number of angles */
	     const flotsam_real* ang, /**< Scattering angles (radians) */
	     const flotsam_real* pf_in, /**< Phase functions in [npf,nang] */
	     const flotsam_real normalization, /**< Integral over surface of sphere */
             flotsam_real* pf_out, /**< Phase functions out [npf,nang] */
 	     flotsam_real* pf_components) { /**< Component amplitudes out [npf,FLOTSAM_NUM_COMPONENTS]*/
#ifdef CATCH_CPP_ERRORS
  try {
#endif
    const Vector my_ang(const_cast<flotsam_real*>(ang), dimensions(nang));
    const Matrix my_pf_in(const_cast<flotsam_real*>(pf_in), dimensions(npf,nang));
    Matrix my_pf_out(pf_out, dimensions(npf,nang));
    Matrix my_pf_components(pf_components, dimensions(npf,FLOTSAM_NUM_COMPONENTS));
    flotsam::analyse_phase_functions(my_ang, my_pf_in, normalization, my_pf_out, my_pf_components);
    return FLOTSAM_SUCCESS;
  }
  catch (...) {
    return FLOTSAM_ADEPT_ERROR;
  }
}


} // extern "C"
