/** @file      reflectance.hpp
    @brief     Compute reflectance of an atmospheric scene
    @copyright 2017 European Centre for Medium Range Weather Forcasts
    @license   Apache License Version 2 (see the NOTICE.md file for details)
 */
#ifndef FLOTSAM_REFLECTANCE_HPP
#define FLOTSAM_REFLECTANCE_HPP

#include "BandProfile.hpp"
#include "containers.hpp"

namespace flotsam {
  using namespace adept;

  /// Compute reflectance of an atmospheric scene
  template <bool IsActive>
  int reflectance(const BandProfile& band,
		  const Array<1,Real,IsActive>& albedo,
		  const intVector& loc_particulate,
		  const Array<1,Real,IsActive>& od_particulate, 
		  const Array<1,Real,IsActive>& ssa_particulate,
		  const Array<1,Real,IsActive>& pf_particulate,
		  const Array<2,Real,IsActive>& pfc_particulate,
		  typename scalar<IsActive>::type& reflectance,
		  WorkingData* working_data = 0 ///< Store intermediate data if not NULL
		  );

}


#endif
