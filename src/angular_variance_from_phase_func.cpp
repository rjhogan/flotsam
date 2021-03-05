/** @file      angular_variance_from_phase_func.cpp
    @brief     Compute angular variance from phase function
    @copyright 2018 European Centre for Medium Range Weather Forcasts
    @license   Apache License Version 2 (see the NOTICE.md file for details)
 */

#include "angular_variance_from_phase_func.hpp"

using namespace adept;

/// Compute angular variance from raw phase functions, in sterad
///
/// The vector pf contains a phase function evenly spaced in angle
/// from 0 to 180
Real
flotsam::angular_variance_from_phase_func(const Vector& pf ///< Phase function
					  )
{
  int nang  = pf.size();
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
  Vector ang = linspace(0,M_PI,nang);
  return sum(pf*weight*ang*ang)/sum(pf*weight);
}


