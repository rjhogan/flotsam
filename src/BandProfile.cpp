/// @file      BandProfile.cpp
/// @brief     Defines the class BandProfile
/// @copyright 2017 European Centre for Medium Range Weather Forcasts
/// @license   Apache License Version 2 (see the NOTICE.md file for details)

#include "BandProfile.hpp"
#include "LookUpTable.hpp"

namespace flotsam {

  int BandProfile::set_geometry(Real mu_sun_, Real mu_inst_, Real azim_) {
    mu_sun = mu_sun_;
    mu_inst = mu_inst_;
    azim = azim_;

    sec_sun = 1.0 / mu_sun_;
    sec_inst = 1.0 / mu_inst_;

    flotsam::lut.get_lobe_props(mu_sun, sec_sun_lobe_dn, sec_sun_lobe_up, frac_sun_lobe_dn);
    flotsam::lut.get_lobe_props(mu_inst, sec_inst_lobe_dn, sec_inst_lobe_up, frac_inst_lobe_dn);

    // FIX surely mu_sun_lobe_dn and mu_inst_lobe_dn should match the
    // mu of the sun and inst since only the sec treats diffusivity
    // effects?
#ifndef FLOTSAM_REVISION3
    mu_sun_lobe_dn  = 1.0 / sec_sun_lobe_dn;
    mu_sun_lobe_up  = 1.0 / sec_sun_lobe_up;
    mu_inst_lobe_dn = 1.0 / sec_inst_lobe_dn;
    mu_inst_lobe_up = 1.0 / sec_inst_lobe_up;
#else
    mu_sun_lobe_dn  = mu_sun;
    mu_sun_lobe_up  = FLOTSAM_MU_LOBE_UP;
    mu_inst_lobe_dn = mu_inst;
    mu_inst_lobe_up = FLOTSAM_MU_LOBE_UP;
#ifdef FLOTSAM_REVISION6
    // Account for the tendency of lobe radiation to bend towards the
    // vertical as the more off-vertical radiation has a greater
    // preference for being scattered away
    if (mu_sun_lobe_dn < 0.5) {
      mu_sun_lobe_dn = 0.15 + mu_sun*(0.35/0.5);
    }
    if (mu_inst_lobe_dn < 0.5) {
      mu_inst_lobe_dn = 0.15 + mu_inst*(0.35/0.5);
    }
#endif
#endif

    Real mu_sun_lobe_dn_eff = std::max(mu_sun, mu_sun_lobe_dn);
    Real mu_inst_lobe_dn_eff = std::max(mu_inst, mu_inst_lobe_dn);

    do_sun_lobe_up  = (frac_sun_lobe_dn  == 1.0) ? false : true;
    do_inst_lobe_up = (frac_inst_lobe_dn == 1.0) ? false : true;

    // Compute great-circle angles betwen the sun and instrument
    // considering also the lobes that may have an effective zenith
    // angle different from the direct beam angle

    //    std::cout << "mus = " << mu_sun << " " << mu_sun_lobe_dn_eff << " " << -mu_sun_lobe_up
    //	      << " " << mu_inst << " " << mu_inst_lobe_dn_eff << " " << -mu_inst_lobe_up << "\n";

    
    /*
    mu_inst_lobe_up *= -1.0;
    mu_sun_lobe_up *= -1.0;
    */

    angle.resize(3,3);
    angle(0,0) = great_circle_angle(mu_sun, mu_inst, azim);
    angle(0,1) = great_circle_angle(mu_sun, mu_inst_lobe_dn_eff, azim);
    angle(0,2) = great_circle_angle(mu_sun, -mu_inst_lobe_up, azim);
    angle(1,0) = great_circle_angle(mu_sun_lobe_dn_eff, mu_inst, azim);
    angle(1,1) = great_circle_angle(mu_sun_lobe_dn_eff, mu_inst_lobe_dn_eff, azim);
    angle(1,2) = great_circle_angle(mu_sun_lobe_dn_eff, -mu_inst_lobe_up, azim);
    angle(2,0) = great_circle_angle(-mu_sun_lobe_up, mu_inst, azim);
    angle(2,1) = great_circle_angle(-mu_sun_lobe_up, mu_inst_lobe_dn_eff, azim);
    angle(2,2) = great_circle_angle(-mu_sun_lobe_up, -mu_inst_lobe_up, azim);

    /*
    mu_inst_lobe_up *= -1.0;
    mu_sun_lobe_up *= -1.0;
    */

    //    std::cout << "angle = " << angle << "\n";

    pf_components.resize(3,3,lut.n_phase_function_components());
    // Extract the phase function components at these angles
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
	Real real_index = angle(i,j) / lut.dang;
	int index = static_cast<int>(real_index);
	if (index < 0) {
	  index = 0;
	}
	else if (index >= lut.nang-1) {
	  index = lut.nang-2;
	}
	Real weight = real_index-index;

	// The phase function components exist in three smoothnesses:
	// (0) unsmoothed, used in direct-to-direct scattering, (1)
	// smoothed once, used in direct-to-lobe and lobe-to-direct
	// scattering, and (2) smoothed twice, used in lobe-to-lobe
	// scattering
	int ismoothness = (i>0) + (j>0);
	pf_components(i,j,__)
	  = (1.0-weight) * lut.phase_function_components(index,   __, ismoothness)
	       + weight  * lut.phase_function_components(index+1, __, ismoothness);
      }
    }
    //    pf_components /= (2.0 * acos(-1.0));
    //    pf_components *= 0.5;
    pf_coeffs.resize(3,3,lut.n_phase_function_components(),
		     FLOTSAM_NUM_SMOOTHING_COEFFS);

    Matrix pf_smooth(lut.n_phase_function_components(),
		     FLOTSAM_NUM_SMOOTHING_COEFFS);
    intVector index_smooth(FLOTSAM_NUM_SMOOTHING_COEFFS);
    index_smooth << 0, 1, 3, 5;
    Vector ang_var_lut = lut.ang_var_smoothing(index_smooth);
    max_angular_variance = ang_var_lut(end);
    Matrix ang_var_powers(FLOTSAM_NUM_SMOOTHING_COEFFS,
			  FLOTSAM_NUM_SMOOTHING_COEFFS);
    ang_var_powers(0,__) = 1.0;
    ang_var_powers(1,__) = ang_var_lut;
    for (int i = 2; i < FLOTSAM_NUM_SMOOTHING_COEFFS; ++i) {
      ang_var_powers(i,__) = ang_var_lut*ang_var_powers(i-1,__);
    }
    ang_var_powers.in_place_transpose();

    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
	Real real_index = angle(i,j) / lut.dang;
	int index = static_cast<int>(real_index);
	if (index < 0) {
	  index = 0;
	}
	else if (index >= lut.nang-1) {
	  index = lut.nang-2;
	}
	Real weight = real_index-index;
	pf_smooth = (1.0-weight) * lut.phase_function_components(index,   __, index_smooth)
	               + weight  * lut.phase_function_components(index+1, __, index_smooth);
	// Solve for the polynomial coefficients
	pf_coeffs(i,j,__,__) = solve(ang_var_powers, pf_smooth.T()).T();
      }
    }

    pf_rayleigh = rayleigh_phase_function(angle(0,0));

    return FLOTSAM_SUCCESS;
  }

  int BandProfile::init(const Channel& chan, const BackgroundProfile& prof) {

    n_g = chan.n_g();
    weight.resize(n_g);
    weight = chan.weight;

    n_z = prof.n_z();

    od_gas_abs.resize(n_g, n_z);
    od_rayleigh.resize(n_g, n_z);
    
    Vector dp = prof.edge_pressure(range(1,end)) - prof.edge_pressure(range(0,end-1));

    // Could be done more efficiently with an outer product function
    for (int i = 0; i < n_g; ++i) {
      od_rayleigh(i,__) = chan.rayleigh_od_per_Pa(i)*dp;
    }

    // For the moment there is no absorption
    od_gas_abs = 0.0;

    return FLOTSAM_SUCCESS;
  }

  int BandProfile::init_direct(const Vector& weight_,
			       const Matrix& od_rayleigh_,
			       const Matrix& od_gas_abs_) {
    if (weight_.size() != od_rayleigh_.size(0)) {
      return FLOTSAM_INCORRECT_NUMBER_OF_SPECTRAL_INTERVALS;
    }
    n_g = weight_.size();
    n_z = od_rayleigh_.size(1);
    if (weight.size() != n_g) {
      weight.clear();
    }
    weight = weight_;
    if (od_rayleigh.size(0) != n_g || od_rayleigh.size(1) != n_g) {
      od_rayleigh.clear();
    }
    od_rayleigh = od_rayleigh_;
    if (od_gas_abs_.empty()) {
      od_gas_abs.resize(n_g,n_z);
      od_gas_abs = 0.0;
    }
    else {
      if (od_gas_abs.size(0) != n_g || od_gas_abs.size(1) != n_g) {
	od_gas_abs.clear();
      }
      od_gas_abs = od_gas_abs_;
    }
    return FLOTSAM_SUCCESS;
  }

  int BandProfile::set_rayleigh_optical_depth(const Vector& od) {
    if (od.size() != n_z) {
      return FLOTSAM_INCORRECT_NUMBER_OF_LAYERS;
    }
    else {
      for (int i = 0; i < n_g; ++i) {
	od_rayleigh(i,__) = od;
      }
      return FLOTSAM_SUCCESS;
    }
  }

  int BandProfile::set_gas_absorption_optical_depth(const Vector& od) {
    if (od.size() != n_z) {
      return FLOTSAM_INCORRECT_NUMBER_OF_LAYERS;
    }
    else {
      for (int i = 0; i < n_g; ++i) {
	od_gas_abs(i,__) = od;
      }
      return FLOTSAM_SUCCESS;
    }
  }

  Real BandProfile::get_total_rayleigh_optical_depth() {
    return dot_product(weight, sum(od_rayleigh,1));
  }

  Real BandProfile::get_total_gas_absorption_optical_depth() {
    return dot_product(weight, sum(od_gas_abs,1));
  }


}
