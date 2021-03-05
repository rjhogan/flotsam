/** @file      delta_eddington_scaling.hpp
    @brief     Perform delta-eddington scaling
    @copyright 2017 European Centre for Medium Range Weather Forcasts
    @license   Apache License Version 2 (see the NOTICE.md file for details)
 */
#ifndef FLOTSAM_DELTA_EDDINGTON_SCALING_HPP
#define FLOTSAM_DELTA_EDDINGTON_SCALING_HPP

namespace flotsam {

  template <bool IsActive>
  void delta_eddington_scaling(const ScatteringProperties<IsActive>& prop,
			       ScatteringProperties<IsActive>& prop_de) {

    Array<1,Real,IsActive> g2 = prop.g * prop.g;
    prop_de.od = prop.od * (1.0 - prop.ssa*g2);
    prop_de.ssa = prop.ssa * (1.0 - g2) / (1.0 - prop.ssa*g2);
    prop_de.g = prop.g / (1.0 + prop.g);

  }

}

#endif
