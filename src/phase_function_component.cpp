/** @file      phase_function_component.cpp
    @brief     Compute components used for reconstructing smoothed phase functions
    @copyright 2017 European Centre for Medium Range Weather Forcasts
    @license   Apache License Version 2 (see the NOTICE.md file for details)
 */

#include "flotsam.hpp"
#include "base.hpp"

namespace flotsam {
  using namespace adept;

  /// Normalize a phase function so that its integral on the sphere
  /// (4pi steradians) is 4pi
  static Vector normalize_phase_function(const Vector& pf) {
    int nang = pf.size();
    Vector ang = linspace(0,M_PI,nang);
    Vector weight = sin(ang);
    weight(0)   = ang(1) / 8.0;
    weight(end) = ang(1) / 8.0;
    weight /= sum(weight);
    // The use of fabs here ensures sensible normalization even for
    // phase function components that integrate to zero
    return pf / dot_product(fabs(pf), weight);
  }

  /// Return a phase function corresponding to smoothing by the
  /// forward lobe, using the exp(-tan(angle)/width)) parametric form. 
  Vector exptan_smoothing_kernel(int nang,       ///< Number of angles between 0 and pi
				 Real width_deg) ///< Width in degrees of the kernel
  {
    Vector ang = linspace(0,M_PI,nang);
    Vector kernel(nang);
    kernel = 0.0;
    kernel.where(ang < M_PI*0.5) = exp(-tan(ang)/(width_deg*M_PI/180.0));
    return normalize_phase_function(kernel);
  }

  /// Return a phase function component and its smoothed versions for
  /// "nang" evenly spaced angles between 0 and pi, where icomponent
  /// is one of the list provided in flotsam.h
  Matrix
  phase_function_component(int nang, flotsam_phase_func_component_t icomponent) {
    Matrix pf(FLOTSAM_NUM_SMOOTHING_VALUES+1, nang);
    // Isotropic phase function is simply 1 everywhere
    if (icomponent == FLOTSAM_PHASE_FUNC_ISOTROPIC) {
      pf = 1.0;
    }
    else {
      Vector ang = linspace(0,M_PI,nang);
      Vector pf_tmp(nang);
      if (icomponent == FLOTSAM_PHASE_FUNC_RAYLEIGH) {
	// This is actually the cos-squared part of the Rayleigh phase
	// function
	pf_tmp = cos(ang);
	pf_tmp *= pf_tmp;
	pf(0,__) = normalize_phase_function(pf_tmp);
      }
      else if (icomponent == FLOTSAM_PHASE_FUNC_COS2_FORWARD) {
	// Cos-squared in foward hemisphere only
	pf_tmp = cos(ang);
	pf_tmp(range((nang-1)/2,end)) = 0.0;
	pf_tmp *= pf_tmp;
	pf(0,__) = normalize_phase_function(pf_tmp);
      }
      else if (icomponent == FLOTSAM_PHASE_FUNC_COS2_BACKWARD) {
	// Cos-squared in backward hemisphere only
	pf_tmp = cos(ang);
	pf_tmp(range(0,(nang-1)/2)) = 0.0;
	pf_tmp *= pf_tmp;
	pf(0,__) = normalize_phase_function(pf_tmp);
      }
      else if (icomponent == FLOTSAM_PHASE_FUNC_CONVEX_LOBE) {
	// Several functional forms are possible; with the exp-tan
	// function, widths between 18 and 40 have been tried
	pf(0,__) = exptan_smoothing_kernel(nang, FLOTSAM_LOBE_HALF_WIDTH_PHASE_FUNC_DEG);
      }
      else if (icomponent == FLOTSAM_PHASE_FUNC_BACKWARD) {
	// A component that peaks more strongly in the backward
	// direction
	Real half_pi = 0.5 * M_PI;
	pf_tmp = (ang-half_pi)*(ang-half_pi)/half_pi;
	pf_tmp(range(0,(nang-1)/2)) = 0.0;
	pf_tmp = sin(pf_tmp);
	pf(0,__) = normalize_phase_function(pf_tmp*pf_tmp);
      }
      else if (icomponent == FLOTSAM_PHASE_FUNC_SLOPE) {
	// Note that the following phase function integrates to zero
	pf(0,__) = normalize_phase_function(sin(ang)*sin(2.0*ang));
      }
      else {
	// Not recognised
	return Matrix();
      }
      // Compute the smoothed versions of the phase function
      for (int i = 1; i < FLOTSAM_NUM_SMOOTHING_VALUES+1; ++i) {
	pf(i,__) = convolve_phase_function(pf(i-1,__));
      }
    }
    return pf.T();
  }

  /*
  void generate_component_phase_functions(int nang)
  {
    Vector ang = linspace(0,M_PI,nang);
    
    Matrix pf(3,nang);
    // cos-squared
    pf(0,__) = normalize_phase_function(cos(ang)*cos(ang));
    pf(1,__) = convolve_phase_function(pf(0,__));
    pf(2,__) = convolve_phase_function(pf(1,__));

    std::cout << "phase_function_cos_squared:\n";
    std::cout << pf.T();

    // wide lobe
    Vector pf_tmp(nang);
    pf_tmp = 0.0;
    pf_tmp.where(ang < M_PI*0.5) = exp(-tan(ang)/(35.0*M_PI/180.0));
    pf(0,__) = normalize_phase_function(pf_tmp);
    pf(1,__) = convolve_phase_function(pf(0,__));
    pf(2,__) = convolve_phase_function(pf(1,__));

    std::cout << "phase_function_lobe:\n";
    std::cout << pf.T();
  }
  */


}
