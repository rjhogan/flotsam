!> @file      flotsam.inc -*- f90 -*-
!> @brief     Fortran header for FLOTSAM solar radiance model
!> @copyright 2021- European Centre for Medium Range Weather Forcasts
!> @license   Apache License Version 2 (see the NOTICE.md file for details)
!
!   This file provides a Fortran interface to the FLOTSAM library
!
  
! GENERAL DEFINITIONS

! This is the type used for all floating-point variables in both
! the interface and internally in the library.
#define FLOTSAM_REAL C_DOUBLE
  
! /** Gas identifiers */
! typedef enum {
!   FLOTSAM_H2O = 0,
!   FLOTSAM_CO2,
!   FLOTSAM_O3,
!   FLOTSAM_O2,
!   FLOTSAM_MAX_GASES
! } flotsam_gas_t;

!> Phase function components
integer, parameter :: FLOTSAM_PHASE_FUNC_ISOTROPIC = 0
integer, parameter :: FLOTSAM_PHASE_FUNC_RAYLEIGH = 1
integer, parameter :: FLOTSAM_PHASE_FUNC_CONVEX_LOBE = 2
integer, parameter :: FLOTSAM_PHASE_FUNC_COS2_FORWARD = 3
integer, parameter :: FLOTSAM_PHASE_FUNC_COS2_BACKWARD = 4
integer, parameter :: FLOTSAM_PHASE_FUNC_BACKWARD = 5
integer, parameter :: FLOTSAM_PHASE_FUNC_SLOPE = 6
integer, parameter :: FLOTSAM_PHASE_FUNC_MAX_COMPONENTS = 7

!> Number of ocean albedo components
integer, parameter :: FLOTSAM_NUM_ALBEDO_COMPONENTS = 4

!> Number of variables returned by flotsam_reflectance_components;
!> should not be used; rather use return value from
!> flotsam_n_phase_function_components()
!integer, parameter :: FLOTSAM_NUM_COMPONENTS = 6

!> Error codes
integer, parameter :: FLOTSAM_ADEPT_ERROR = -16
integer, parameter :: FLOTSAM_FIRST_ERROR_CODE = -15
integer, parameter :: FLOTSAM_INCORRECT_NUMBER_OF_SPECTRAL_INTERVALS = -15
integer, parameter :: FLOTSAM_INCORRECT_NUMBER_OF_ALBEDO_COMPONENTS = -14
integer, parameter :: FLOTSAM_ERROR_MERGING_GASES = -13
integer, parameter :: FLOTSAM_AUTOMATIC_DIFFERENTIATION_ERROR = -12
integer, parameter :: FLOTSAM_INPUTS_MUST_START_AT_TOA = -11
integer, parameter :: FLOTSAM_REQUIRED_GAS_NOT_PRESENT = -10
integer, parameter :: FLOTSAM_TOO_MANY_GASES_IN_LUT = -9
integer, parameter :: FLOTSAM_FEATURE_NOT_AVAILABLE = -8
integer, parameter :: FLOTSAM_INCORRECT_NUMBER_OF_LAYERS = -7
integer, parameter :: FLOTSAM_INVALID_CONTEXT = -6
integer, parameter :: FLOTSAM_INVALID_GAS = -5
integer, parameter :: FLOTSAM_INPUT_OUT_OF_PHYSICAL_BOUNDS = -4
integer, parameter :: FLOTSAM_ERROR_READING_FILE = -3
integer, parameter :: FLOTSAM_FILE_NOT_FOUND = -2
integer, parameter :: FLOTSAM_MEMORY_ERROR = -1
integer, parameter :: FLOTSAM_SUCCESS = 0


interface

  !> @name Support functions
  !> @{ 

  !> Return a pointer to a NULL-terminated string explaining an error
  !> code
  !const char* flotsam_error_message(int error_code);

  ! /** Return the name of the gas ("H2O" etc) or NULL if not
  !     recognized */
  ! const char* flotsam_gas_name(flotsam_gas_t igas);

  !> Return the value of the phase function at the specified
  !> scattering angle, assuming phase function elements to be equally
  !> spaced in angle including the values at 0 and 180 degrees */
  real(FLOTSAM_REAL) function flotsam_interp_phase_func(n, pf, ang) bind(C)
    use, intrinsic :: ISO_C_BINDING
    integer(C_INT),     intent(in), value :: n     !> Num of phase func elements
    real(FLOTSAM_REAL), intent(in)        :: pf(n) !> Phase function
    real(FLOTSAM_REAL), intent(in), value :: ang   !> Requested angle (radians)
  end function flotsam_interp_phase_func

  !> Pi minus the great circle angle (radians) between two points on a sphere.
  !> In terms of latitude, the formula is
  !> acos(sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(dlon), but note
  !> that mu is the cosine of the zenith angle so it is equal to the
  !> sine of the elevation angle.
  real(FLOTSAM_REAL) function flotsam_scattering_angle(mu_sun, mu_inst, azim) bind(C)
    use, intrinsic :: ISO_C_BINDING
    real(FLOTSAM_REAL), intent(in), value :: mu_sun   !> Cosine of solar zenith angle
    real(FLOTSAM_REAL), intent(in), value :: mu_inst  !> Cosine of instrument zenith angle 
    real(FLOTSAM_REAL), intent(in), value :: azim     !> Azimuthal separation (radians)
  end function flotsam_scattering_angle

  !> Return the number of components used to parameterize phase
  !> functions (excluding the delta component)
  integer(C_INT) function flotsam_n_phase_function_components() bind(C)
    use, intrinsic :: ISO_C_BINDING
  end function flotsam_n_phase_function_components

  !> @}

  !> @name Channel functions
  !>
  !> A satellite channel context stores look-up tables and other
  !> information about a particular channel that is independent of the
  !> properties of the atmospheric profile that it will be applied to.

  !> @{

  !> Create a new satellite channel context, configured using the
  !> information in the file provided, and return the ID for that
  !> context. The ID is positive or zero; a negative number indicates
  !> that an error occurred with the error codes as indicated above.
  ! integer function flotsam_new_channel(const char* filename) bind(C)

  !> Create new satellite channel context in a spectral region where
  !> gas absorption can be neglected and only molecular (Rayleigh)
  !> scattering is included. The ID is positive or zero; a negative
  !> number indicates that an error occurred with the error codes as
  !> indicated above.
  integer(C_INT) function flotsam_new_channel_rayleigh_only(wavelength) bind(C)
    use, intrinsic :: ISO_C_BINDING
    real(FLOTSAM_REAL), intent(in) :: wavelength
  end function flotsam_new_channel_rayleigh_only

  !> Create new satellite channel context in vacuum, or equivalently
  !> in a spectral region with no Rayleigh scattering or
  !> absorption. The ID is positive or zero; a negative number
  !> indicates that an error occurred with the error codes as
  !> indicated above.
  integer(C_INT) function flotsam_new_channel_vacuum() bind(C)
    use, intrinsic :: ISO_C_BINDING
  end function flotsam_new_channel_vacuum

  !> Free the memory associated with a satellite channel context with
  !> ID "ichan", or return FLOTSAM_INVALID_CONTEXT if not recognized
  integer(C_INT) function flotsam_free_channel(ichan) bind(C)
    use, intrinsic :: ISO_C_BINDING
    integer(C_INT), intent(in), value :: ichan
  end function flotsam_free_channel

  !> Return the wavelength of the channel, in metres
  integer(C_INT) function flotsam_get_wavelength(ichan, wavelength) bind(C)
    use, intrinsic :: ISO_C_BINDING
    integer(C_INT), intent(in), value :: ichan
    real(FLOTSAM_REAL), intent(out) :: wavelength
  end function flotsam_get_wavelength

!   /** @} */


!   /** @name Background profile functions

!       A background profile context stores the physical properties of a
!       particular atmospheric profile, such as temperature, pressure
!       and gas concentrations. It does not store cloud or aerosol
!       properties, and is independent of a particular satellite
!       channel.
!       @{
!   */

!   /** Create a new background profile context, and return an ID for
!       that context. The ID is positive or zero; a negative number
!       indicates that an error occurred with the error codes as
!       indicated above. */
!   int flotsam_new_background_profile();

!   /** Free the memory associated with background profile context with
!       ID "iprof", or return FLOTSAM_INVALID_CONTEXT if not recognized */
!   int flotsam_free_background_profile(int iprof);

!   /** Reset the size of a background profile context with ID
!       "iprof" */
!   int flotsam_reset_background_profile(
!      int iprof,                           /**< ID of the profile */
!      int n                                /**< Number of layers in profile */
! 				       );

!   /** Set pressure at layer edges in a background profile context with
!       ID "iprof" */
!   int flotsam_set_edge_pressure(
!      int iprof,                           /**< ID of the profile */
!      int n,                               /**< Number of layers in profile */
!      const flotsam_real* edge_pressure    /**< Edge pressure in Pa (n+1 elements) */
! 				     );

!   /** Set layer-mean temperature in a background profile context with
!       ID "iprof" */
!   int flotsam_set_temperature(
!      int iprof,                           /**< ID of the profile */
!      int n,                               /**< Number of layers in profile */
!      const flotsam_real* temperature      /**< Temperature in K (n elements) */
! 			      );

!   /** Set the mass mixing ratio of a well-mixed gas with gas type
!       "igas" */
!   int flotsam_set_gas_const(int iprof,
! 			    flotsam_gas_t igas,
! 			    flotsam_real mass_mixing_ratio);

!   /** Set the mass mixing ratio of a gas */
!   int flotsam_set_gas_profile(int iprof,          /**< ID of the profile */
! 			      flotsam_gas_t igas, /**< Gas ID */
! 			      int n,              /**< Number of layers (must match stored value) */
! 			      const flotsam_real* mass_mixing_ratio
! 			                          /**< Mass mixing ratio */
! 			      );

!   /** @} */

  !> @name Band profile functions

  !> A band profile context stores the optical properties of the
  !> atmosphere and surface at the wavelength of a particular
  !> satellite channel.
  
  !> Create a new band profile context, and return an ID for that
  !> context. The ID is positive or zero; a negative number indicates
  !> that an error occurred with the error codes as indicated above.
  integer(C_INT) function flotsam_new_band_profile() bind(C)
    use, intrinsic :: ISO_C_BINDING
  end function flotsam_new_band_profile

  !> Free memory associated with a band profile context, or return
  !> FLOTSAM_INVALID_CONTEXT if not recognized.
  integer(C_INT) function flotsam_free_band_profile(iband) bind(C)
    use, intrinsic :: ISO_C_BINDING
    integer(C_INT), intent(in), value :: iband
  end function flotsam_free_band_profile

!   /** Initialize a band profile context indicated by "iband" by
!       calculating the optical properties of the atmospheric profile and
!       the Planck function from the atmospheric profile indicated by
!       "iprof" at the satellite channel indicated by "ichan" */
!   int flotsam_init_band_profile(int iband, int ichan, int iprof);

  !> Initialize a band profile context indicated by "iband" by
  !> calculating the optical properties of the atmospheric profile and
  !> the Planck function from the atmospheric profile indicated by
  !> "iprof" at the satellite channel indicated by "ichan" */
  integer(C_INT) function flotsam_init_band_profile_direct(iband, n_g, n_z, &
       &  weight, od_rayleigh, od_gas_abs) bind(C)
    use, intrinsic :: ISO_C_BINDING
    integer(C_INT), intent(in), value :: iband
    integer(C_INT), intent(in), value :: n_g
    integer(C_INT), intent(in), value :: n_z
    real(FLOTSAM_REAL), intent(in) :: weight(n_g)
    real(FLOTSAM_REAL), intent(in) :: od_rayleigh(n_z,n_g)
    real(FLOTSAM_REAL), intent(in), optional :: od_gas_abs(n_z,n_g)
  end function flotsam_init_band_profile_direct

  !> Set the number and meaning of the components that will be used
  !> to describe the phase function. Note that the maximum number
  !> that can be provided is FLOTSAM_PHASE_FUNC_MAX_COMPONENTS, and
  !> FLOTSAM_PHASE_FUNC_DELTA cannot be provided in the list, as this
  !> is determined as the residual from the others.
  integer(C_INT) function flotsam_set_phase_func_components(iband, &
       & n_phase_func_components, pfc_list) bind(C)
    use, intrinsic :: ISO_C_BINDING
    integer(C_INT), intent(in) :: iband
    integer(C_INT), intent(in) :: n_phase_func_components
    integer(C_INT), intent(in) :: pfc_list(n_phase_func_components)
  end function flotsam_set_phase_func_components

  !> Set the observing geometry
  integer(C_INT) function flotsam_set_geometry(iband, mu_sun, mu_inst, azim) &
       &  bind(C)
    use, intrinsic :: ISO_C_BINDING
    integer(C_INT),     value, intent(in) :: iband   !< Band profile ID
    real(FLOTSAM_REAL), value, intent(in) :: mu_sun  !< Cosine of solar zenith angle
    real(FLOTSAM_REAL), value, intent(in) :: mu_inst !< Cosine of sensor zenith angle
    real(FLOTSAM_REAL), value, intent(in) :: azim    !< Azimuth angle between sun & instrument (radians)
  end function flotsam_set_geometry

  !> Override the Rayleigh optical depth of each layer, where the
  !> number of layers must match the layers already in the band 
  !   int flotsam_set_rayleigh_optical_depth(int iband,               !< Band profile ID 
  ! 					 int n,                   !< Number of layers 
  ! 					 const flotsam_real* od); !< Layer optical depth 

  !    Override the gas absorption optical depth of each layer, where
  !       the number of layers must match the layers already in the
  !       band 
  !   int flotsam_set_gas_absorption_optical_depth(int iband,         !< Band profile ID 
  ! 				       int n,                     !< Number of layers 
  ! 				       const flotsam_real* od);   !< Layer optical depth 
  
  !    Get the Rayleigh optical depth of the entire atmosphere 
  !   int flotsam_get_total_rayleigh_optical_depth(int iband,         !< Band profile ID 
  ! 					       flotsam_real* od); !< Returned optical depth 

  !    Get the gas absorption optical depth of the entire atmosphere 
  !   int flotsam_get_total_gas_absorption_optical_depth(int iband,         !< Band profile ID 
  ! 						     flotsam_real* od); !< Returned optical depth 

  !> @}

  !> @name Ocean BRDF functions

  !> An Ocean-BRDF context stores look-up tables needed to compute the
  !> ocean bi-directional reflectance function at a particular
  !> wavelength 

  !> @{

  !> Create a new Ocean-BRDF context for the specified wavelength and
  !> wind speeds, and return the ID for that context. The ID is
  !> positive or zero; a negative number indicates that an error
  !> occurred with the error codes as indicated above. Default values
  !> of pigment concentration and salinity are assumed.
  integer(C_INT) function flotsam_new_ocean_brdf(wavelength, nwinds, winds) bind(C)
    use, intrinsic :: ISO_C_BINDING
    real(FLOTSAM_REAL), intent(in), value :: wavelength    !< Wavelength (m)
    integer(C_INT),     intent(in), value :: nwinds        !< Number of wind speeds
    real(FLOTSAM_REAL), intent(in)        :: winds(nwinds) !< Wind speeds (m/s)
  end function flotsam_new_ocean_brdf


  !> As flotsam_new_ocean_brdf but also specifying the ocean pigment
  !> concentration, salinity, whether or not to apply shadowing 
  !> and the shape of the function describing the probability 
  !> distribution of ocean waves
  integer(C_INT) function flotsam_new_ocean_brdf_detailed(wavelength, nwinds, &
       &  winds, pigment_conc_mg_m3, salinity_ppt, apply_shadowing, slope_dist_shape) bind(C)
    use, intrinsic :: ISO_C_BINDING
    real(FLOTSAM_REAL), intent(in) :: wavelength    !< Wavelength (m)
    integer(C_INT),     intent(in) :: nwinds        !< Number of wind speeds
    real(FLOTSAM_REAL), intent(in) :: winds(nwinds) !< Wind speeds (m/s)
    real(FLOTSAM_REAL), intent(in) :: pigment_conc_mg_m3
    real(FLOTSAM_REAL), intent(in) :: salinity_ppt
    integer(C_INT),     intent(in) :: apply_shadowing  !< 1=yes, 0=no
    integer(C_INT),     intent(in) :: slope_dist_shape !< 1=anisotropic (wind dependent) 
                               !< 2=anisotropic with Gram-Charlier coefficients (default)
                               !< 3=isotropic (wind independent)
  end function flotsam_new_ocean_brdf_detailed

  !> Compute the 4 ocean albedo components (white-sky albedo,
  !> black-sky albedo = diffuse reflection from direct source,
  !> reflected radiance for a diffuse source reflected radiance for a
  !> direct source) for the specified inputs.
  integer(C_INT) function flotsam_get_ocean_albedo_components(ibrdf, mu_sun, &
       &  mu_inst, azim, wind, albedos) bind(C)
    use, intrinsic :: ISO_C_BINDING
    integer(C_INT),     intent(in), value :: ibrdf      !< Ocean BRDF ID
    real(FLOTSAM_REAL), intent(in), value :: mu_sun     !< Zenith angle of sun
    real(FLOTSAM_REAL), intent(in), value :: mu_inst    !< Zenith angle of instrument
    real(FLOTSAM_REAL), intent(in), value :: azim       !< Azimuth angle between sun & instrument (radians)
    real(FLOTSAM_REAL), intent(in), value :: wind       !< Wind speed (m/s)
    real(FLOTSAM_REAL), intent(out) :: albedos(*) !< Returned albedo components
  end function flotsam_get_ocean_albedo_components

  !> Free memory associated with an Ocean BRDF context, or return
  !> FLOTSAM_INVALID_CONTEXT if not recognized.
  integer(C_INT) function flotsam_free_ocean_brdf(ibrdf) bind(C)
    use, intrinsic :: ISO_C_BINDING
    integer(C_INT),     intent(in), value  :: ibrdf      !< Ocean BRDF ID
  end function flotsam_free_ocean_brdf

  !> @}

  !> @name Reflectance functions

  !> Compute the reflectance for the band profile indicated by
  !> "iband", for the particulate optical properties provided. These
  !> values are only stored for the layers containing particulates,
  !> the number of which is given by "n", and the location of which in
  !> flotsam_set_background_profile is given by "loc" as a list of
  !> zero-based indices. The phase function components form a matrix
  !> whose fastest varying dimension is of length
  !> n_phase_func_components. If nalbedo=1 then a Lambertian surface
  !> is assumed with albedo given in albedo[0], while if nalbedo=4
  !> then albedo[0] is the white-sky albedo, albedo[1] is the
  !> black-sky albedo (diffuse reflection from direct source),
  !> albedo[2] is the reflected radiance for a diffuse source and
  !> albedo[3] is the reflected radiance for a direct source.
  integer(C_INT) function flotsam_reflectance(iband, nalbedo, albedo, &
       &  n, loc, od, ssa, pf, pfc, reflectance) bind(C)
    use, intrinsic :: ISO_C_BINDING
    integer(C_INT), intent(in), value :: iband         !< Band profile ID 
    integer(C_INT), intent(in), value :: nalbedo       !< Number of albedo components 
    real(FLOTSAM_REAL), intent(in) :: albedo(nalbedo)  !< Surface albedo components 
    integer(C_INT), intent(in), value :: n             !< Number of particulate layers 
    integer(C_INT), intent(in) :: loc(n)               !< Location of particulate layers 
    real(FLOTSAM_REAL), intent(in) :: od(n)            !< Optical depth 
    real(FLOTSAM_REAL), intent(in) :: ssa(n)           !< Single-scattering albedo 
    real(FLOTSAM_REAL), intent(in) :: pf(n)            !< Phase func for single scattering 
    real(FLOTSAM_REAL), intent(in) :: pfc(*)         !< Phase function components 
    real(FLOTSAM_REAL), intent(out) :: reflectance !< Returned reflectance 
  end function flotsam_reflectance
  
  !    As flotsam_reflectance but return an array whose elements are
  !       [0] reflectance, [1] reflectance due to quasi-single scattering,
  !       [2] reflectance involving the lobe flux, [3] reflectance
  !       involving the diffuse flux, [4] surface contribution to
  !       reflectance, and [5] top-of-atmosphere normalized upwelling
  !       flux.  Note that elements [1] to [4] sum to element [0]. 
  !   int flotsam_reflectance_components(int iband,  !< Band profile ID 
  ! 		     int nalbedo,                !< Number of albedo components 
  ! 		     const flotsam_real* albedo, !< Surface albedo components 
  ! 		     int n,                      !< Number of particulate layers 
  ! 		     const int* loc,             !< Location of particulate layers 
  ! 		     const flotsam_real* od,     !< Optical depth 
  ! 		     const flotsam_real* ssa,    !< Single-scattering albedo 
  ! 		     const flotsam_real* pf,     !< Phase func for single scattering 
  ! 		     const flotsam_real* pfc,    !< Phase function components 
  ! 		     flotsam_real data[FLOTSAM_NUM_COMPONENTS] !< Returned data 
  ! 				     );
  
  !    As flotsam_reflectance, but with particulate variables specified at all
!       heights.  "n" must be equal to the number of layers in the band
!       profile, or a FLOTSAM_INCORRECT_NUMBER_OF_LAYERS error will be
!       returned. 
!   int flotsam_reflectance_unindexed(int iband,               !< Band profile ID 
! 				    int nalbedo,             !< Number of albedo components 
! 				    const flotsam_real* albedo, !< Surface albedo components 
! 				    int n,                   !< Number of particulate layers 
! 				    const flotsam_real* od,  !< Optical depth 
! 				    const flotsam_real* ssa, !< Single-scattering albedo 
! 				    const flotsam_real* pf,  !< Phase func for single scattering 
! 				    const flotsam_real* pfc, !< Phase function components 
! 				    flotsam_real* reflectance!< Returned reflectance 
! 				    );

!    As flotsam_reflectance, but also computes the derivatives of the
!       reflectance with respect to the input particulate optical
!       properties. This function only works if the FLOTSAM library was
!       compiled with the Adept automatic differentiation library,
!       otherwise it returns FLOTSAM_FEATURE_NOT_AVAILABLE. 
!   int flotsam_reflectance_jacobian(int iband,               !< Band profile ID 
! 				   int nalbedo,             !< Number of albedo components 
! 				   const flotsam_real* albedo, !< Surface albedo components 
! 				   int n,                   !< Number of particulate layers 
! 				   const int* loc,          !< Location of particulate layers 
! 				   const flotsam_real* od,  !< Optical depth 
! 				   const flotsam_real* ssa, !< Single-scattering albedo 
! 				   const flotsam_real* pf,  !< Phase func for single scattering 
! 				   const flotsam_real* pfc, !< Phase function components 
! 				   flotsam_real* reflectance,!< Returned reflectance 
! 				   flotsam_real* d_albedo,  !< d(reflectance)/d(albedo) 
! 				   flotsam_real* d_od,      !< d(reflectance)/d(od) 
! 				   flotsam_real* d_ssa,     !< d(reflectance)/d(ssa) 
! 				   flotsam_real* d_pf,      !< d(reflectance)/d(pf) 
! 				   flotsam_real* d_pfc      !< d(reflectance)/d(pfc) 
! 				   );
  
!    As flotsam_reflectance_jacobian, but avoiding the Adept start-up cost by
!       keeping the Adept stack in global memory - this is faster but is
!       thread unsafe 
!   int flotsam_reflectance_jacobian_thread_unsafe(int iband, !< Band profile ID 
! 				   int nalbedo,             !< Number of albedo components 
! 				   const flotsam_real* albedo, !< Surface albedo components 
! 				   int n,                   !< Number of particulate layers 
! 				   const int* loc,          !< Location of particulate layers 
! 				   const flotsam_real* od,  !< Optical depth 
! 				   const flotsam_real* ssa, !< Single-scattering albedo 
! 				   const flotsam_real* pf,  !< Phase func for single scattering 
! 				   const flotsam_real* pfc, !< Phase function components 
! 				   flotsam_real* reflectance,!< Returned reflectance 
! 				   flotsam_real* d_albedo,  !< d(reflectance)/d(albedo) 
! 				   flotsam_real* d_od,      !< d(reflectance)/d(od) 
! 				   flotsam_real* d_ssa,     !< d(reflectance)/d(ssa) 
! 				   flotsam_real* d_pf,      !< d(reflectance)/d(pf) 
! 				   flotsam_real* d_pfc      !< d(reflectance)/d(pfc) 
! 				   );

!    As flotsam_reflectance but writes an ASCII file containing the
!       results of internal calculations, such as profiles of different
!       components of the radiation distribution. 
!   int flotsam_reflectance_write_internals(int iband,  !< Band profile ID 
! 		     int nalbedo,                !< Number of albedo components 
! 		     const flotsam_real* albedo, !< Surface albedo components 
! 		     int n,                      !< Number of particulate layers 
! 		     const int* loc,             !< Location of particulate layers 
! 		     const flotsam_real* od,     !< Optical depth 
! 		     const flotsam_real* ssa,    !< Single-scattering albedo 
! 		     const flotsam_real* pf,     !< Phase func for single scattering 
! 		     const flotsam_real* pfc,    !< Phase function components 
! 		     flotsam_real data[FLOTSAM_NUM_COMPONENTS], !< Returned data 
! 		     const char* file_name     !< Name of file to write 
! 					  );

  !> @}

  !> @name Phase-function utilities

  !> @{

  !> Create a new band profile context, and return an ID for that
  !> context. The ID is positive or zero; a negative number indicates
  !> that an error occurred with the error codes as indicated above.
  integer(C_INT) function flotsam_analyse_phase_functions(npf, &
       &  nang, ang, pf_in, normalization, pf_out, pf_components) bind(C)
    use, intrinsic :: ISO_C_BINDING
    integer(C_INT), intent(in), value :: npf  !< Number of phase functions
    integer(C_INT), intent(in), value :: nang !< Number of angles
    real(FLOTSAM_REAL), intent(in) :: ang(nang) !< Scattering angle (radians)
    real(FLOTSAM_REAL), intent(in) :: pf_in(nang,npf) !< Phase function in
    real(FLOTSAM_REAL), intent(in), value :: normalization !< Integral over surf of sphere
    real(FLOTSAM_REAL), intent(out) :: pf_out(nang,npf) !< Phase function out
    real(FLOTSAM_REAL), intent(out) :: pf_components(*)
  end function flotsam_analyse_phase_functions

  !> @}


end interface
