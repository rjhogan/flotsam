/// @file      base.hpp
/// @brief     Defines basic stuff 
/// @copyright 2017-2018 European Centre for Medium Range Weather Forcasts
/// @license   Apache License Version 2 (see the NOTICE.md file for details)

//#include <exception>

#ifndef M_PI
#define M_PI 3.14159265358979323846  /* pi */
#endif

/*
  Current optimal settings for CLOUD:
  DO NOT define FLOTSAM_REVISION*
  DO NOT define FLOTSAM_NO_PIFM
  DO NOT define FLOTSAM_DIFFUSE_SCALING_WITH_OD
  FLOTSAM_FORWARD_LOBE_DIFFUSE = 0.0
  FLOTSAM_DIFFUSE_UP_FACTOR = 1.1;
  FLOTSAM_DIFFUSE_DN_FACTOR = 1.0;

  Older optimal settings for AEROSOL:
  DO define FLOTSAM_REVISION but not FLOTSAM_REVISION2
  DO NOT define FLOTSAM_NO_PIFM
  DO NOT define FLOTSAM_DIFFUSE_SCALING_WITH_OD
  FLOTSAM_FORWARD_LOBE_DIFFUSE = 0.2
  FLOTSAM_DIFFUSE_UP_FACTOR = 1.1;
  FLOTSAM_DIFFUSE_DN_FACTOR = 1.0;

  Newer settings for AEROSOL:flotsam0513g
  DO define FLOTSAM_REVISION but not FLOTSAM_REVISION2
  DO NOT define FLOTSAM_NO_PIFM
  DO define FLOTSAM_DIFFUSE_SCALING_WITH_OD
  FLOTSAM_DIFFUSE_UP_FACTOR = 1.2;
  FLOTSAM_DIFFUSE_DN_FACTOR = 1.2;
  FLOTSAM_FORWARD_LOBE_DIFFUSE = 0.2;

  Settings for AEROSOL:flotsam0514h
  DO define FLOTSAM_REVISION but not FLOTSAM_REVISION2
  DO NOT define FLOTSAM_NO_PIFM
  DO define FLOTSAM_DIFFUSE_SCALING_WITH_OD (New formula)
  FLOTSAM_DIFFUSE_UP_FACTOR = 1.0;
  FLOTSAM_DIFFUSE_DN_FACTOR = 1.0;
  FLOTSAM_NONDIRECT_SCALING = 1.15; // Replaces those above
  FLOTSAM_FORWARD_LOBE_DIFFUSE = 0.2;
  FLOTSAM_LUT_OPTICAL_DEPTH_REF_FACTOR = 0.5; // Previously hard-coded
  OMG!!! There was a bug in the plotting that removed from the total
  half of the lobe reflectance... the scheme was then tuned to this!

 */



#ifndef FLOTSAM_BASE_HPP
#define FLOTSAM_BASE_HPP

#define FLOTSAM_REVISION 1
//#define FLOTSAM_REVISION2 1 // OFF fixes non-conservation in flux profile
#define FLOTSAM_REVISION3 1
#define FLOTSAM_REVISION3A 1
#define FLOTSAM_REVISION4 1
//#define FLOTSAM_REVISION5 1 // OFF fixes reflectance at low viewing angle
#define FLOTSAM_REVISION6 1

#define FLOTSAM_NO_PIFM 1
#define FLOTSAM_DIFFUSE_DIFFUSIVITY 2
#define FLOTSAM_BETA_PREFACTOR 0.375 // PIFM value
//#define FLOTSAM_BETA_PREFACTOR 0.5

// At low optical depth there is enhanced scattering from diffuse to
// other streams
#define FLOTSAM_DIFFUSE_SCALING_WITH_OD 1
// Distribution of diffuse angles in optically thin cloud depends on
// the optical depth considered; if we set this equal to the optical
// depth of the cloud we get this "reference" optical depth
//#define FLOTSAM_DIFFUSE_MULTIPLIER_REFERENCE_OD 2.0
// But since diffuse radiation typically originates from scattering in
// the middle of the cloud, the reference optical depth perhaps should
// be doubled
#define FLOTSAM_DIFFUSE_MULTIPLIER_REFERENCE_OD 5.0

// Do something similar with lobe reflectances at low optical depth
#define FLOTSAM_LOBE_SCALING_WITH_OD 1
#define FLOTSAM_LOBE_MULTIPLIER_REFERENCE_OD 4.0

//#define FLOTSAM_USE_OLD_COMPONENTS 1

#ifdef FLOTSAM_USE_OLD_COMPONENTS
// Original values
#define FLOTSAM_LOBE_HALF_WIDTH_PHASE_FUNC_DEG 35.0
#define FLOTSAM_LOBE_HALF_WIDTH_DIFFUSIVITY_DEG 25.0
#define FLOTSAM_LOBE_HALF_WIDTH_KERNEL_DEG 20.0
#else
#define FLOTSAM_LOBE_HALF_WIDTH_PHASE_FUNC_DEG 18.0
#define FLOTSAM_LOBE_HALF_WIDTH_DIFFUSIVITY_DEG 25.0
#define FLOTSAM_LOBE_HALF_WIDTH_KERNEL_DEG 20.0
//#define FLOTSAM_LOBE_HALF_WIDTH_DIFFUSIVITY_DEG 30.0
//#define FLOTSAM_LOBE_HALF_WIDTH_KERNEL_DEG 30.0
#endif

// Not needed any more
//#define FLOTSAM_DIFFUSE_SCALING_MAX 2.0
//#define FLOTSAM_DIFFUSE_SCALING_OD_SCALE 0.4

#include <flotsam.h>

#include <adept_arrays.h>

namespace flotsam {

  using namespace adept;

  /// Select an active versus static scalar depending on bool template
  /// argument
  template <bool IsActive>
  struct scalar {
    typedef adept::Real type;
  };

  template <>
  struct scalar<true> {
    typedef adept::aReal type;
  };


  /*
  class exception : public std::exception {
  public:
    exception(const std::string& message = "Error in FLOTSAM library")
      { message_ = message; }
    virtual const char* what() const throw() { return message_.c_str(); }
    virtual ~exception() throw() { }
  protected:
    std::string message_;
  };
  */

  static const int FLOTSAM_BEAM_FLUX_TERMS = 5;
  static const int FLOTSAM_DIRECT = 0;
  static const int FLOTSAM_LOBE_DN_1 = 1;
  static const int FLOTSAM_LOBE_DN_2 = 2;
  static const int FLOTSAM_LOBE_UP_1 = 3;
  static const int FLOTSAM_LOBE_UP_2 = 4;

  static const adept::Real FLOTSAM_MAX_LAYER_OPTICAL_DEPTH = 20.0;

  // Fudge-factor accounting for enhanced scattering from diffuse
  // upward stream to lobe/direct returning beams
  static const adept::Real FLOTSAM_DIFFUSE_UP_FACTOR = 1.0;
  static const adept::Real FLOTSAM_DIFFUSE_DN_FACTOR = 1.0;

  // Preferred fudge factor
  //static const adept::Real FLOTSAM_NONDIRECT_SCALING = 1.05;
  static const adept::Real FLOTSAM_NONDIRECT_SCALING = 1.0;

  // Standard value for this was 0.5 until 0.5.13
  static const adept::Real FLOTSAM_LUT_OPTICAL_DEPTH_REF_FACTOR = 0.5;

  // Fraction of the forward lobe that is scattered into the
  // forward-hemisphere diffuse
  static const adept::Real FLOTSAM_FORWARD_LOBE_DIFFUSE = 0.0;
  //  static const adept::Real FLOTSAM_FORWARD_LOBE_DIFFUSE = 0.1;
  //  static const adept::Real FLOTSAM_FORWARD_LOBE_DIFFUSE = 0.2;
  //  static const adept::Real FLOTSAM_COS2_FORWARD_TO_LOBE = 0.5;
  static const adept::Real FLOTSAM_COS2_FORWARD_TO_LOBE = 0.5;
  //static const adept::Real FLOTSAM_COS2_FORWARD_TO_LOBE = 1.0;

  // Original value appropriate for pure exponential with angle 14deg
  //  static const adept::Real FLOTSAM_FORWARD_LOBE_REMAIN = 0.5846;
  //static const adept::Real FLOTSAM_FORWARD_LOBE_REMAIN = 0.2;
  static const adept::Real FLOTSAM_FORWARD_LOBE_REMAIN = 0.55;
  //  static const adept::Real FLOTSAM_FORWARD_LOBE_REMAIN = 0.4;
  //  static const adept::Real FLOTSAM_FORWARD_LOBE_REMAIN = 0.8;
  // Value for exp-tan function with angle 18deg
  //  static const adept::Real FLOTSAM_FORWARD_LOBE_REMAIN = 0.734424657048408;
  //static const adept::Real FLOTSAM_FORWARD_LOBE_REMAIN = 0.0;

  //static const adept::Real FLOTSAM_MU_LOBE_UP  = 0.4;
  //  static const adept::Real FLOTSAM_MIN_DIFFUSIVITY_MU_LOBE_DN = 0.15;
  static const adept::Real FLOTSAM_MU_LOBE_UP  = 0.35;
  static const adept::Real FLOTSAM_MIN_DIFFUSIVITY_MU_LOBE_DN = 0.1;

  // Lobe angular width is predicted and used to select the level of
  // smoothing to use in the version of the phase function for the
  // lobe contribution
  static const int FLOTSAM_NUM_SMOOTHING_VALUES = 5;

  // Number of polynomial coefficients to use to describe the
  // smoothing dependence
  static const int FLOTSAM_NUM_SMOOTHING_COEFFS = 4;


}

#endif
