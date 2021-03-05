/** @file      adding.hpp
    @brief     Adding method to compute diffuse fluxes
    @copyright 2017 European Centre for Medium Range Weather Forcasts
    @license   Apache License Version 2 (see the NOTICE.md file for details)
 */
#ifndef FLOTSAM_ADDING_HPP
#define FLOTSAM_ADDING_HPP

#include "containers.hpp"

namespace flotsam {

  template <bool IsActive>
  void adding(const DiffuseProperties<IsActive>& props,
	      const Array<1,Real,IsActive>& surface_albedo,
	      const typename scalar<IsActive>::type& surface_source,
	      DiffuseFluxes<IsActive>& fluxes) {

    typedef Array<1,Real,IsActive> avector;

    int nz = props.size();

    avector albedo(nz+1), source(nz+1);
    albedo(end) = surface_albedo(0);
    source(end) = surface_source;

    avector inv_denominator(nz);

    for (int i = nz-1; i >= 0; --i) {
      inv_denominator(i) = 1.0 / (1.0 - albedo(i+1)*props.reflectance(i));
      albedo(i) = props.reflectance(i)
	+ props.transmittance(i)*props.transmittance(i)
	*albedo(i+1)*inv_denominator(i);
      source(i) = props.source_up(i) + props.transmittance(i)
	*(source(i+1) + albedo(i+1)*props.source_dn(i))*inv_denominator(i);
    }

    fluxes.dn(0) = 0.0;
    fluxes.up(0) = source(0);

    for (int i = 0; i < nz; ++i) {
      fluxes.dn(i+1) = inv_denominator(i)
	* (props.transmittance(i)*fluxes.dn(i)
	   + props.reflectance(i)*source(i+1) + props.source_dn(i));
      fluxes.up(i+1) = albedo(i+1)*fluxes.dn(i+1) + source(i+1);
    }

    // Needed in calc_reflectance_from_fluxes to compute effective
    // diffusivity
    fluxes.total_source = sum(props.source_up+props.source_dn);

  }

}

#endif
