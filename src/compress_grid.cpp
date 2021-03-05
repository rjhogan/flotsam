/** @file      compress_grid.cpp
    @brief     Compress grid to reduce clear-sky layers
    @copyright 2017 European Centre for Medium Range Weather Forcasts
    @license   Apache License Version 2 (see the NOTICE.md file for details)
 */

#include <flotsam.hpp>

using namespace adept;

/// Compress an input grid to reduce contiguous clear-sky layers
/// into single layers
void
flotsam::compress_grid(const Vector& edge_pressure, ///< Half-level pressures
		       const intVector& location,   ///< Locations of cloud/aerosol
		       Vector& edge_pressure_out,   ///< Compressed-grid pressures
		       intVector& location_out)     ///< Compressed-grid locations
{
  Vector edge_pressure_new(edge_pressure.size());
  location_out.resize(location.size());

  int iz = 1; // Index to grid used in FLOTSAM that includes
              // consolidated clear-sky layers
  
  // First constituent layer - is there a clear-sky layer above?
  edge_pressure_new(0) = edge_pressure(0);
  if (location(0) > 0) {
    // ...yes
    edge_pressure_new(1) = edge_pressure(location(0));
    ++iz;
    // Constituents start at the second layer
    location_out(0) = 1;
  }
  else {
    // Constituents start at the first layer
    location_out(0) = 0;
  }
  for (int iloc = 1; iloc < location.size(); ++iloc) {
    if (location(iloc) > location(iloc-1)+1) {
      // Insert a clear-sky layer by adding an edge pressure at the
      // base of the previous constituent layer
      edge_pressure_new(iz) = edge_pressure(location(iloc-1)+1);
      ++iz;
    }
    // Add edge pressure at top of current constituent layer
    edge_pressure_new(iz) = edge_pressure(location(iloc));
    location_out(iloc) = iz;
    ++iz;
  }
  // Add edge pressure at base of last constituent layer
  edge_pressure_new(iz) = edge_pressure(location(end)+1);
  ++iz;

  if (location(end) < edge_pressure.size()-1) {
    // Add a clear-sky layer down to the surface
    edge_pressure_new(iz) = edge_pressure(end);
    ++iz;
  }

  edge_pressure_out.clear();
  edge_pressure_out.link(edge_pressure_new(range(0,iz-1)));
}
