/** @file      containers.hpp
    @brief     Classes for storing intermediate radiation variables
    @copyright 2017 European Centre for Medium Range Weather Forcasts
    @license   Apache License Version 2 (see the NOTICE.md file for details)
 */

#ifndef FLOTSAM_CONTAINERS_HPP
#define FLOTSAM_CONTAINERS_HPP

#include <iostream>
#include <fstream>
#include <string>

#include "base.hpp"

namespace flotsam {

  using namespace adept;

  /// Basic optical properties for two-stream calculations
  template <bool IsActive>
  struct ScatteringProperties {
    ScatteringProperties(int n) : od(n), ssa(n), g(n) {}
    ScatteringProperties(int n, bool use_g) : od(n), ssa(n) {
      if (use_g) {
	g.resize(n);
      }
    }
    /// Constructor that links to input variables
    ScatteringProperties(Array<1,Real,IsActive> od_,
			 Array<1,Real,IsActive> ssa_,
			 Array<1,Real,IsActive> g_) 
      : ssa(ssa_), g(g_) {
      od.resize(od_.size());
      for (int i = 0; i < od_.size(); ++i) {
	if (od_(i) < FLOTSAM_MAX_LAYER_OPTICAL_DEPTH) {
	  od(i) = od_(i);
	}
	else {
	  od(i) = FLOTSAM_MAX_LAYER_OPTICAL_DEPTH
	    *(2.0 - exp(-(od_(i)-FLOTSAM_MAX_LAYER_OPTICAL_DEPTH)
			/FLOTSAM_MAX_LAYER_OPTICAL_DEPTH));
	}
      }
    }
    
    /// Return number of layers
    int size() const { return od.size(); }

    Array<1,Real,IsActive> od;  // Optical depth
    Array<1,Real,IsActive> ssa; // Single scattering albedo
    Array<1,Real,IsActive> g;   // Asymmetry factor
  };
  
  /// Probability of a scattering event leading to light ray entering
  /// different streams
  template <bool IsActive>
  struct ScatteringProbabilities {
  ScatteringProbabilities(int n) : lobe(n), diffuse_dn(n), diffuse_up(n) {}
    Array<1,Real,IsActive> lobe;      
    Array<1,Real,IsActive> diffuse_dn;
    Array<1,Real,IsActive> diffuse_up;
  };
  
  
  /// The "beam" components of the flux, specifically the quasi-direct
  /// beam and the forward lobe
  template <bool IsActive>
  struct BeamFluxes {
    typedef typename scalar<IsActive>::type areal;
    typedef Array<1,Real,IsActive> avector;
    typedef Array<2,Real,IsActive> amatrix;

    /// Number of components needed depending on whether the upwelling
    /// lobe will be simulated
    static int n_comp(bool do_lobe_up) {
      return do_lobe_up ? FLOTSAM_BEAM_FLUX_TERMS : 3;
    }

    BeamFluxes(int n, bool do_lobe_up = true)
      : prefactor(n_comp(do_lobe_up),n),
	exponent(n_comp(do_lobe_up),n),
	expm_exponent(n_comp(do_lobe_up),n),
	prefactor_source_dn(n_comp(do_lobe_up),n),
	prefactor_source_up(n_comp(do_lobe_up),n),
	// The following constructors lead to a link to a subset of
	// the prefactor and exponent matrices
	prefactor_direct(prefactor(FLOTSAM_DIRECT,__)),
	prefactor_lobe_dn_1(prefactor(FLOTSAM_LOBE_DN_1,__)),
	prefactor_lobe_dn_2(prefactor(FLOTSAM_LOBE_DN_2,__)),
	exponent_direct(exponent(FLOTSAM_DIRECT,__)),
        exponent_lobe_dn_1(exponent(FLOTSAM_LOBE_DN_1,__)),
        exponent_lobe_dn_2(exponent(FLOTSAM_LOBE_DN_2,__)),
        expm_exponent_direct(expm_exponent(FLOTSAM_DIRECT,__)),
        expm_exponent_lobe_dn_1(expm_exponent(FLOTSAM_LOBE_DN_1,__)),
        expm_exponent_lobe_dn_2(expm_exponent(FLOTSAM_LOBE_DN_2,__)),
        lobe_ang_var_dn(n),
        sec(n_comp(do_lobe_up)) {
      
      if (do_lobe_up) {
	prefactor_lobe_up_1.link(prefactor(FLOTSAM_LOBE_UP_1,__));
	prefactor_lobe_up_2.link(prefactor(FLOTSAM_LOBE_UP_2,__));
        exponent_lobe_up_1.link (exponent(FLOTSAM_LOBE_UP_1,__));
        exponent_lobe_up_2.link (exponent(FLOTSAM_LOBE_UP_2,__));
        expm_exponent_lobe_up_1.link(expm_exponent(FLOTSAM_LOBE_UP_1,__));
        expm_exponent_lobe_up_2.link(expm_exponent(FLOTSAM_LOBE_UP_2,__));
      }
    }

    Array<1,Real,IsActive> direct_dn() {
      return prefactor(FLOTSAM_DIRECT,__);
    }
    Array<1,Real,IsActive> lobe_dn() {
      return prefactor(FLOTSAM_LOBE_DN_1,__)
	+    prefactor(FLOTSAM_LOBE_DN_2,__);
    }
    Array<1,Real,IsActive> lobe_up() {
      if (do_lobe_up()) {
	return prefactor(FLOTSAM_LOBE_UP_1,__)
	  +    prefactor(FLOTSAM_LOBE_UP_2,__);
      }
      else {
	Array<1,Real,IsActive> null(prefactor.dimension(1));
	null = 0.0;
	return null;
      }
    }
    
    /// Return number of atmospheric layers
    int size() const {
      return prefactor.size(1);
    }
    /// Return number of terms needed to describe the beam fluxes:
    /// either 3 or 5
    int n_components() const {
      return prefactor.size(0);
    }

    bool do_lobe_up() const {
      return n_components() == FLOTSAM_BEAM_FLUX_TERMS;
    }

    // Downwelling flux at specified half-level
    areal dn(int i) const {
      if (i < size()) {
	// Simply sum the appropriate prefactors
	return prefactor_direct(i)+prefactor_lobe_dn_1(i)+prefactor_lobe_dn_2(i);
      }
      else if (i == size()) {
	// Surface value
	int j = i-1;
	return prefactor_direct(j)*exp(-exponent_direct(j))
	  + prefactor_lobe_dn_1(j)*exp(-exponent_lobe_dn_1(j))
	  + prefactor_lobe_dn_2(j)*exp(-exponent_lobe_dn_2(j));
      }
      else {
	abort(); // Programming error
      }
    }

    areal lobe_up_toa() const {
      if (do_lobe_up()) {
	return prefactor(FLOTSAM_LOBE_UP_1,0) + prefactor(FLOTSAM_LOBE_UP_2,0);
      }
      else {
	return 0.0;
      }
    }
    
    amatrix prefactor;
    amatrix exponent;
    amatrix expm_exponent; // Cache of exp(-exponent)
    amatrix prefactor_source_dn;
    amatrix prefactor_source_up;
    // The following are links to rows of the matrices above
    avector prefactor_direct;
    avector prefactor_lobe_dn_1;
    avector prefactor_lobe_dn_2;
    avector prefactor_lobe_up_1;
    avector prefactor_lobe_up_2;
    avector exponent_direct;
    avector exponent_lobe_dn_1;
    avector exponent_lobe_dn_2;
    avector exponent_lobe_up_1;
    avector exponent_lobe_up_2;
    avector expm_exponent_direct;
    avector expm_exponent_lobe_dn_1;
    avector expm_exponent_lobe_dn_2;
    avector expm_exponent_lobe_up_1;
    avector expm_exponent_lobe_up_2;
    avector lobe_ang_var_dn; ///< Angular variance of lobe (radians^2)
    Vector sec;
  };
  
  template <bool IsActive>
  struct DiffuseProperties {
    DiffuseProperties(int n) 
      : reflectance(n), transmittance(n), source_dn(n), source_up(n) { }
    int size() const { return reflectance.size(); }
    Array<1,Real,IsActive> reflectance, transmittance, source_dn, source_up;
  };

  /// Diffuse fluxes (should these be stored as prefactor/exponent?
  template <bool IsActive>
  struct DiffuseFluxes {
    typedef typename scalar<IsActive>::type areal;
    DiffuseFluxes(int n) : dn(n+1), up(n+1) { }
    Array<1,Real,IsActive> dn, up;
    areal total_source;
  };

  // Structure containing the intermediate data used to compute
  // reflectance
  struct WorkingData {
    WorkingData() : output_stream(&std::cout) { clear(); }
    void clear() {
      ref_direct = ref_lobe = ref_lobe_up = ref_diffuse = ref_surface = flux_toa = 0.0;
    }
    void open() {
      // "-" denotes standard output; otherwise associate
      // output_stream with a file
      if (output_file_name != "-") {
	output_file.open(output_file_name.c_str());
	output_stream = &output_file;
      }
#if ADEPT_VERSION >= 20005
      save_print_style = adept::get_array_print_style();
      adept::set_array_print_style(PRINT_STYLE_MATLAB);
#endif
    }
    template <typename T>
    void write(const char* name, T x) {
      *output_stream << name << " = " << x << ";\n";
    }
    void write_comment_line(const char* txt) {
      *output_stream << "% " << txt << "\n";
    }
    void close() {
      if (output_file_name != "-") {
	output_file.close();
      }
      output_stream = &std::cout;
#if ADEPT_VERSION >= 20005
      adept::set_array_print_style(save_print_style);
#endif
    }
    /*
    ScatteringProperties<false> prop_prime;
    ScatteringProperties<false> prop_delta_ed;
    ScatteringProbabilities<false> prob_direct_to, prob_lobe_dn_to, prob_lobe_up_to;
    BeamFluxes<false> beam_fluxes, return_fluxes;
    DiffuseProperties<false> diffuse_props;
    DiffuseFluxes<false> diffuse_fluxes;
    */
    Real ref_direct, ref_lobe, ref_lobe_up, ref_diffuse, ref_surface, flux_toa;
    std::string output_file_name;
    std::ofstream output_file;
    std::ostream* output_stream;
    adept::ArrayPrintStyle save_print_style;
  };


  template <bool IsActive>
  std::ostream& operator<<(std::ostream& os, 
			  const ScatteringProperties<IsActive>& prop) {
    os << "prop_od = " << prop.od << "\n";
    os << "prop_ssa = " << prop.ssa << "\n";
    os << "prop_g = " << prop.g << "\n";
    return os;
  }

  template <bool IsActive>
  std::ostream& operator<<(std::ostream& os, 
			  const ScatteringProbabilities<IsActive>& prob) {
    os << "prob_lobe = " << prob.lobe << "\n";
    os << "prob_diffuse_dn = " << prob.diffuse_dn << "\n";
    os << "prob_diffuse_up = " << prob.diffuse_up << "\n";
    return os;
  }

  template <bool IsActive>
  std::ostream& operator<<(std::ostream& os, 
			  const BeamFluxes<IsActive>& bf) {
    /*
    os << "beam_prefactor_direct = " << bf.prefactor_direct << "\n";
    os << "beam_prefactor_lobe_dn_1 = " << bf.prefactor_lobe_dn_1 << "\n";
    os << "beam_prefactor_lobe_dn_2 = " << bf.prefactor_lobe_dn_2 << "\n";
    os << "beam_prefactor_lobe_up_1 = " << bf.prefactor_lobe_up_1 << "\n";
    os << "beam_prefactor_lobe_up_2 = " << bf.prefactor_lobe_up_2 << "\n";
    os << "beam_exponent_direct = " << bf.exponent_direct << "\n";
    os << "beam_exponent_lobe_dn_1 = " << bf.exponent_lobe_dn_1 << "\n";
    os << "beam_exponent_lobe_dn_2 = " << bf.exponent_lobe_dn_2 << "\n";
    os << "beam_exponent_lobe_up_1 = " << bf.exponent_lobe_up_1 << "\n";
    os << "beam_exponent_lobe_up_2 = " << bf.exponent_lobe_up_2 << "\n";
    */
    os << "beam_flux_direct_dn = " << bf.prefactor(FLOTSAM_DIRECT,__) << "\n";
    os << "beam_flux_lobe_dn = "   << bf.prefactor(FLOTSAM_LOBE_DN_1,__) 
      + bf.prefactor(FLOTSAM_LOBE_DN_2,__) << "\n";
    if (bf.do_lobe_up()) {
      os << "beam_flux_lobe_up = "   << bf.prefactor(FLOTSAM_LOBE_UP_1,__)
	+ bf.prefactor(FLOTSAM_LOBE_UP_2,__) << "\n";
    }

    //    os << "beam_flux_direct_surf = " << bf.prefactor_direct(end)*exp(-bf.exponent_direct(end)) << "\n";
    //    os << "beam_flux_lobe_dn_surf = " << bf.prefactor_lobe_dn_1(end)*exp(-bf.exponent_lobe_dn_1(end))
    //      + bf.prefactor_lobe_dn_2(end)*exp(-bf.exponent_lobe_dn_2(end)) << "\n";

    os << "beam_prefactor = " << bf.prefactor << "\n";
    os << "beam_exponent = " << bf.exponent << "\n";
    os << "beam_prefactor_source_dn = " << bf.prefactor_source_dn << "\n";
    os << "beam_prefactor_source_up = " << bf.prefactor_source_up << "\n";
    os << "beam_sec = " << bf.sec << "\n";
    os << "beam_lobe_ang_var_dn = " << bf.lobe_ang_var_dn << "\n";
    return os;
  }

  template <bool IsActive>
  std::ostream& operator<<(std::ostream& os, 
			  const DiffuseProperties<IsActive>& dp) {
    os << "diffuse_props_reflectance   = " << dp.reflectance << "\n";
    os << "diffuse_props_transmittance = " << dp.transmittance << "\n";
    os << "diffuse_props_source_dn     = " << dp.source_dn << "\n";
    os << "diffuse_props_source_up     = " << dp.source_up << "\n";
    return os;
  }
 
  template <bool IsActive>
  std::ostream& operator<<(std::ostream& os, 
			  const DiffuseFluxes<IsActive>& df) {
    os << "diffuse_flux_dn = " << df.dn << "\n";
    os << "diffuse_flux_up = " << df.up << "\n";
    os << "total_source    = " << df.total_source << "\n";
    return os;
  }
}

#endif
  
