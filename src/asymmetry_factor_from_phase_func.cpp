/** @file      asymmetry_factor_from_phase_func.cpp
    @brief     Compute asymmetry factor from phase function
    @copyright 2017 European Centre for Medium Range Weather Forcasts
    @license   Apache License Version 2 (see the NOTICE.md file for details)
 */

#include "asymmetry_factor_from_phase_func.hpp"

using namespace adept;

/// Compute asymmetry factor from raw phase functions
///
/// The matrix pf contains phase functions evenly spaced in angle from
/// 0 to 180, with angle being the second, fastest varying dimension
/// of the matrix
Vector
flotsam::asymmetry_factor_from_phase_func(const Matrix& pf, ///< Phase functions
					  bool back_hem_only)
{
  int nang  = pf.size(1);
  Real dang = M_PI / (nang - 1.0);
  Vector cos_mid_ang = cos(linspace(-0.5*dang, M_PI+0.5*dang, nang+1));
  // Set first value to cos(0);
  cos_mid_ang(0) = 1.0;
  // Set last value to cos(pi);
  cos_mid_ang(nang) = -1.0;
  // Area of the surface of a unit sphere corresponding to each
  // range of scattering angle
  Vector weight = -2.0*M_PI*(cos_mid_ang(range(1,end))
			   - cos_mid_ang(range(0,end-1)));
  // If we only want the backward hemisphere then set the forward
  // hemisphere weight to zero
  if (back_hem_only) {
    if (nang % 2) {
      // Odd number of angles
      int imid = (nang-1)/2;
      weight(range(0,imid-1)) = 0.0;
      weight(imid) *= 0.5;
    }
    else {
      // Even number of angles
      weight(range(0,nang/2-1)) = 0.0;
    }
  }
  // Weight multiplied by cosine of scattering angle
  Vector weight_x_cosang = weight * cos(linspace(0, M_PI, nang));

  // Unfortunately Adept-1.9.8 has a bug in matrix-vector
  // multiplication with both dimensions of matrix strided
#if ADEPT_VERSION <= 10908
  // Asymmetry factor is mean of cos(scattering angle) weighted by
  // phase function
  Matrix pf_copy;
  pf_copy = pf;
  Vector g = (pf_copy ** weight_x_cosang) / (fabs(pf_copy) ** weight);
#else
  Vector g = (pf ** weight_x_cosang) / (fabs(pf) ** weight);
#endif


  // Remove very small numbers due to rounding error
  g.where(fabs(g) < 1.0e-8) = 0.0;
  return g;
}


