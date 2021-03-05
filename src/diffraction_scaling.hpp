/** @file      diffraction_scaling.hpp
    @brief     Perform diffraction scaling, removing the narrow forward lobe
    @copyright 2017 European Centre for Medium Range Weather Forcasts
    @license   Apache License Version 2 (see the NOTICE.md file for details)
 */
#ifndef FLOTSAM_DIFFRACTION_SCALING_HPP
#define FLOTSAM_DIFFRACTION_SCALING_HPP

#include "LookUpTable.hpp"

namespace flotsam {

  template <bool IsActive>
  void diffraction_scaling(const ScatteringProperties<IsActive>& prop,
			   const Array<2,Real,IsActive>& pfc_particulate,
			   ScatteringProperties<IsActive>& prop_prime) {

    Array<1,Real,IsActive> prob_not_delta
      = sum(pfc_particulate(__,lut.mass_component_index),1);
#ifndef FLOTSAM_REVISION
    prop_prime.od = prop.od * prob_not_delta;
    prop_prime.ssa = 1.0 + (prop.ssa - 1.0) / prob_not_delta;  // Can go negative!
    // Don't need asymmetry factor
#else
    Array<1,Real,IsActive> one_minus_ssa_delta = 1.0 - (1.0-prob_not_delta)*prop.ssa;
    prop_prime.od  = prop.od  * one_minus_ssa_delta;
    prop_prime.ssa = prop.ssa * prob_not_delta / one_minus_ssa_delta;
#endif
  }
}

#endif
