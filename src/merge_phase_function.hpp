/** @file      merge_phase_function.hpp
    @brief     Merge phase functions of particulate and gas scattering properties
    @copyright 2017 European Centre for Medium Range Weather Forcasts
    @license   Apache License Version 2 (see the NOTICE.md file for details)
 */
#ifndef FLOTSAM_MERGE_PHASE_FUNCTION_HPP
#define FLOTSAM_MERGE_PHASE_FUNCTION_HPP

#include "containers.hpp"
#include "LookUpTable.hpp"

namespace flotsam {

  /// Merge gas and particulate optical properties
  ///
  /// Optical depth is summed, single-scattering albedo is an average
  /// weighted by optical depth, and asymmetry factor is weighted by
  /// scattering optical depth but treating gas asymmetry factor as
  /// zero.
  template <bool IsActive>
  void merge_phase_function(const Vector& od_rayleigh, Real pf_gas,
			    const intVector& loc_particulate,
			    const ScatteringProperties<IsActive>& prop_particulate,
			    const Array<1,Real,IsActive>& pf_particulate,
			    const Array<2,Real,IsActive>& pfc_particulate,
			    Array<1,Real,IsActive>& pf,
			    Array<2,Real,IsActive>& pfc) {

    typedef Array<1,Real,IsActive> avector;

    Vector scat_od_gas = od_rayleigh(loc_particulate);
    avector scat_od_particulate = prop_particulate.od * prop_particulate.ssa;

    pf = pf_gas;

    pf(loc_particulate) = (pf_gas*scat_od_gas
			   + pf_particulate*scat_od_particulate)
      / (scat_od_gas + scat_od_particulate);

    avector weight = scat_od_particulate / (scat_od_gas + scat_od_particulate);

    // First set to Rayleigh scattering
    pfc = spread<0>(lut.rayleigh_components,pfc.dimension(0));
    
    for (int i = 0; i < loc_particulate.size(); ++i) {
      int loc = loc_particulate(i);
      pfc(loc,__) = (1.0-weight(i))*pfc(loc,__) + weight(i)*pfc_particulate(i,__);
      // Remove delta function, which is the implicit residual from one
      pfc(loc,__) /= sum(pfc(loc,lut.mass_component_index));
    }


  }

}

#endif
