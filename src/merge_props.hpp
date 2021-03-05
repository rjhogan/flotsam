/** @file      merge_props.hpp
    @brief     Merge particulate and gas scattering properties
    @copyright 2017 European Centre for Medium Range Weather Forcasts
    @license   Apache License Version 2 (see the NOTICE.md file for details)
 */
#ifndef FLOTSAM_MERGE_PROPS_HPP
#define FLOTSAM_MERGE_PROPS_HPP

#include "containers.hpp"

namespace flotsam {

  /// Merge gas and particulate optical properties
  ///
  /// Optical depth is summed, single-scattering albedo is an average
  /// weighted by optical depth, and asymmetry factor is weighted by
  /// scattering optical depth but treating gas asymmetry factor as
  /// zero.
  template <bool IsActive>
  void merge_props(const Vector& od_gas_abs,
		   const Vector& od_rayleigh,
		   const intVector& loc_particulate,
		   const ScatteringProperties<IsActive>& prop_particulate,
		   ScatteringProperties<IsActive>& prop)
  {
    // Optical depth and single-scattering albedo of gases
    prop.od = od_gas_abs + od_rayleigh;
    prop.ssa = 0.0;
    prop.ssa.where(prop.od > 0.0) = od_rayleigh / prop.od;

    // Add particulate contribution
    prop.od(loc_particulate) += prop_particulate.od;
    prop.ssa(loc_particulate) = (od_rayleigh(loc_particulate)
				 + prop_particulate.od*prop_particulate.ssa)
      / prop.od(loc_particulate);

    if (!prop.g.empty()) {
      prop.g = 0.0;
      prop.g(loc_particulate) = prop_particulate.od*prop_particulate.ssa*prop_particulate.g
	/ (prop.od(loc_particulate)*prop.ssa(loc_particulate));
    }
  }

}

#endif
