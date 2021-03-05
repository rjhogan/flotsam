/** @file      LookUpTable.hpp
    @brief     Declares the struct LookUpTable
    @copyright 2017 European Centre for Medium Range Weather Forcasts
    @license   Apache License Version 2 (see the NOTICE.md file for details)
 */
#ifndef FLOTSAM_LOOK_UP_TABLE_HPP
#define FLOTSAM_LOOK_UP_TABLE_HPP

#include <vector>

#include "base.hpp"

namespace flotsam {

  /// Look-up table for coefficients that depend on
  /// scattering/elevation angle
  struct LookUpTable {
    LookUpTable();
    
    /// Set phase function components to use
    void set_phase_function_components(int n_components, /// Number of components
				       flotsam_phase_func_component_t* id /// IDs of components
				       );

    /// Get the properties of the forward lobe
    void get_lobe_props(Real mu,           ///< Zenith angle of beam
			Real& sec_lobe_dn, ///< Effective secant of downwelling part of lobe
			Real& sec_lobe_up, ///< Effective secant of upwelling part of lobe
			Real& frac_dn      ///< Fraction of lobe in downwelling part
			) const;

    /// Get the "raw", smoothed and twice-smoothed phase functions
    /// corresponding to the phase-function component amplitudes
    /// provided in pf_components
    Matrix get_phase_function(const Vector& pf_components) const;

    // Return the one of the 4 or 5 phase function "basis functions"
    Vector phase_function_basis(int icomponent) const {
      if (icomponent < 0 || icomponent >= n_phase_function_components()) {
	//	throw flotsam::exception("Phase function basis index out of range");
	return Vector();
      }
      return phase_function_components(__,icomponent,0);
    }

    // Return the number of phase function components
    int n_phase_function_components() const { return phase_function_component_id.size(); }

    // Data
    Real lobe_spread;         ///< Fraction of wide scattering... FIX 
    Real forward_lobe_diffuse;///< Fraction of lobe that is scattered into forward diffuse hemisphere
    Real forward_lobe_remain; ///< Fraction of lobe remaining in lobe after a lobe-scattering event

    Vector diffusivity_dn, diffusivity_up;
    Vector fraction_dn;

    /// For each index of a phase function decomposition, store the ID
    /// of the basis function shape
    std::vector<flotsam_phase_func_component_t> phase_function_component_id;

    /// For each basis function shape, store the index to it in a
    /// phase function decomposition, or -1 if that basis function is
    /// not used
    intVector phase_function_component_index;

    /// Phase function components or basis functions, dimensioned
    /// (angle,component,smoothing), where smoothing is of length
    /// FLOTSAM_NUM_SMOOTHING_VALUES+1, smoothing==0 corresponds to no
    /// smoothing, and smoothing==n to n convolutions with the wide
    /// phase function
    Array3D phase_function_components;
    Vector ang_var_smoothing; ///< Angular variance of smoothing levels
    Vector g_components; ///< Asymmetry factor of phase function components
    //    Vector g_back_hem_components; ///< Asymmetry factor of backward hemisphere only
    ///< Index of the components with "mass" (those that integrate to
    ///1), i.e. excluding FLOTSAM_PHASE_FUNC_SLOPE which integrates to
    ///zero.
    intVector mass_component_index;
    Vector rayleigh_components; ///< Weight of each component to yield
				///Rayleigh scattering
    Real dang; ///< Angular resolution of phase_function_components
    int nang; ///< Number of angles in phase_function_components
    int i_lobe_component; ///< Index to lobe component of phase function
    Real ang_var_convex_lobe; ///< Angular variance of convex lobe basis function
    Real ang_var_cos2_forward; ///< Angular variance of forward cos-squared basis function

    Real dmu; ///< Resolution of look-up table (radians)
    int nza; ///< Number of zenith-angles in look-up table

  };

  /// Single global instance of LookUpTable accessible with
  /// flotsam::lut
  extern const LookUpTable lut;

}

#endif
