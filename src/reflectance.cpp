/** @file      reflectance.cpp
    @brief     Compute reflectance of an atmospheric scene
    @copyright 2017 European Centre for Medium Range Weather Forcasts
    @license   Apache License Version 2 (see the NOTICE.md file for details)
 */

#include "BandProfile.hpp"

#include "adding.hpp"
#include "base.hpp"
#include "calc_asymmetry_factor.hpp"
#include "calc_angular_variance.hpp"
#include "calc_beam_fluxes.hpp"
#include "calc_reflectance_from_fluxes.hpp"
#include "calc_scattering_probabilities.hpp"
#include "calc_smoothed_phase_function.hpp"
#include "containers.hpp"
#include "delta_eddington_scaling.hpp"
#include "diffraction_scaling.hpp"
#include "meador_weaver.hpp"
#include "merge_phase_function.hpp"
#include "merge_props.hpp"
#include "reflectance.hpp"
#include "bi_directional_surface_reflectivity.hpp"

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
		  WorkingData* working_data ///< Store/print intermediate data if not NULL
		  ) {

    typedef typename scalar<IsActive>::type areal;
    typedef Array<1,Real,IsActive> avector;
    typedef Array<2,Real,IsActive> amatrix;
    typedef Array<3,Real,IsActive> aarray3;

    // Check for sun below the horizon
    if (band.mu_sun <= 0.0) {
      if (working_data) {
	working_data->clear();
      }
      reflectance = 0.0;
      return FLOTSAM_SUCCESS;
    }

    // Check for albedo components
    if (albedo.size() != 1 && albedo.size() != 4) {
      std::cerr << "albedo size = " << albedo.size() << "\n";
      return FLOTSAM_INCORRECT_NUMBER_OF_ALBEDO_COMPONENTS;
    }

    // Do we write the internal data to file?
    bool write_internals = false;
    if (working_data && !working_data->output_file_name.empty()) {
      write_internals = true;
    }

    // Number of levels including particulate-free
    int n_z = band.n_z;
    int n_particulate = loc_particulate.size();

    // Merged scattering properties
    //ScatteringProperties<IsActive> prop(n_z);
    // Scaled scattering properties
    ScatteringProperties<IsActive> prop_prime(n_z, false);
    ScatteringProperties<IsActive> prop_delta_ed(n_z);

    // Probability that the direct/quasi-direct beam, and the down and
    // upwardly propagating parts of the forward lobe, are scattered
    // to the lobe or into the downward or upward hemispheres
    ScatteringProbabilities<IsActive> prob_direct_to(n_z), 
      prob_lobe_dn_to(n_z), prob_lobe_up_to(n_z);

    // Direct and lobe flux components
    BeamFluxes<IsActive> beam_fluxes(n_z,  band.do_sun_lobe_up);
    BeamFluxes<IsActive> return_fluxes(n_z,band.do_inst_lobe_up);

    // Layer-wise diffuse reflectance, transmittance and
    // upwelling/downwelling sources
    DiffuseProperties<IsActive> diffuse_props(n_z);

    // Upwelling and downwelling diffuse fluxes
    DiffuseFluxes<IsActive> diffuse_fluxes(n_z);

    // Phase function and components
    amatrix pfc(n_z, lut.n_phase_function_components());
    aarray3 beam_to_beam(3,3,n_z);
    avector pf;
    pf.link(beam_to_beam(0,0,__));

    // Scattering properties after delta-Eddington and diffraction
    // scaling
    ScatteringProperties<IsActive> prop_particulate_delta_ed(n_particulate);
    ScatteringProperties<IsActive> prop_particulate_prime(n_particulate, false);

    // Angular variance of radiation scattered into the lobe (sterad)
    avector angular_variance_lobe;

    // Calculate particulate asymmetry factor from phase function
    // components
    avector g_particulate;
    //    calc_asymmetry_factor(band.component_g, pfc_particulate, g_particulate);
    calc_asymmetry_factor(pfc_particulate, g_particulate);

    // Store unscaled particulate scattering properties
    const ScatteringProperties<IsActive> prop_particulate(od_particulate,
							  ssa_particulate,
							  g_particulate);
    // Compute scaled particulate scattering properties
    delta_eddington_scaling(prop_particulate, prop_particulate_delta_ed);
    diffraction_scaling(prop_particulate, pfc_particulate, prop_particulate_prime);

    reflectance = 0.0;

    if (working_data) {
      working_data->clear();
    }

    if (write_internals) {
      working_data->open();
      working_data->write_comment_line("This file contains inputs and outputs to the FLOTSAM radiance model for a single profile");
      working_data->write_comment_line("0. Look-up table data");
      working_data->write("lut_mu", linspace(0.0, 1.0, lut.nza));
      working_data->write("lut_diffusivity_dn", lut.diffusivity_dn);
      working_data->write("lut_diffusivity_up", lut.diffusivity_up);
      working_data->write("lut_fraction_dn", lut.fraction_dn);
      working_data->write("lut_phase_function_components", lut.phase_function_components(__,__,0));
      working_data->write("lut_phase_function_components_smooth1", lut.phase_function_components(__,__,1));
      working_data->write("lut_half_width_diffusivity", FLOTSAM_LOBE_HALF_WIDTH_DIFFUSIVITY_DEG);
      working_data->write("lut_cos2_forward_to_lobe", FLOTSAM_COS2_FORWARD_TO_LOBE);
      working_data->write("lut_angular_variance_smoothing", lut.ang_var_smoothing);

      working_data->write_comment_line("1. Geometry data");
      working_data->write("mu_sun", band.mu_sun);
      working_data->write("mu_sun_lobe_dn", band.mu_sun_lobe_dn);
      working_data->write("mu_sun_lobe_up", band.mu_sun_lobe_up);
      working_data->write("frac_sun_lobe_dn", band.frac_sun_lobe_dn);
      working_data->write("diffusivity_sun", band.sec_sun);
      working_data->write("diffusivity_sun_lobe_dn", band.sec_sun_lobe_dn);
      working_data->write("diffusivity_sun_lobe_up", band.sec_sun_lobe_up);
      working_data->write("mu_inst", band.mu_inst);
      working_data->write("mu_inst_lobe_dn", band.mu_inst_lobe_dn);
      working_data->write("mu_inst_lobe_up", band.mu_inst_lobe_up);
      working_data->write("frac_inst_lobe_dn", band.frac_inst_lobe_dn);
      working_data->write("diffusivity_inst", band.sec_inst);
      working_data->write("diffusivity_inst_lobe_dn", band.sec_inst_lobe_dn);
      working_data->write("diffusivity_inst_lobe_up", band.sec_inst_lobe_up);
      working_data->write("azimuthal_separation", band.azim);
      working_data->write("great_circle_separations", band.angle);
      working_data->write("n_layers", band.n_z);
      working_data->write("n_g_points", band.n_g);
      working_data->write("smooth_phase_function_coeffs", band.pf_coeffs(0,0,__,__));

      working_data->write_comment_line("2. Input data");
      working_data->write("albedo", albedo);
      working_data->write("location_particulate", loc_particulate);
      working_data->write("optical_depth_particulate", od_particulate);
      working_data->write("single_scattering_albedo_particulate", ssa_particulate);
      working_data->write("phase_function_particulate", pf_particulate);
      working_data->write("phase_function_components_particulate", pfc_particulate);

      working_data->write_comment_line("3. Computed particulate data");
      working_data->write("asymmetry_factor_particulate", g_particulate);
    }


    areal surface_source;

    // Loop over g point
    for (int ig = 0; ig < band.n_g; ++ig) {
      // Merge gas optical properties (assuming Rayleigh scattering
      // and absorption) with particulate properties
      merge_props(band.od_gas_abs(ig,__), band.od_rayleigh(ig,__),
		  loc_particulate, prop_particulate_delta_ed,
		  prop_delta_ed);
      merge_props(band.od_gas_abs(ig,__), band.od_rayleigh(ig,__),
		  loc_particulate, prop_particulate_prime,
		  prop_prime);

      // Do the same for phase function components
      merge_phase_function(band.od_rayleigh(ig,__), band.pf_rayleigh,
			   loc_particulate, prop_particulate,
			   pf_particulate, pfc_particulate, pf, pfc);

      // Compute the angular variance of part of the phase function
      // scattered into the lobe
      calc_angular_variance(pfc, angular_variance_lobe);

      // Compute the probability that the direct/quasi-direct beam,
      // and the forward lobe, are scattered to the lobe or into the
      // downward or upward hemispheres
      calc_scattering_probabilities(band.mu_sun, band.mu_sun_lobe_dn, band.mu_sun_lobe_up,
				    band.frac_sun_lobe_dn, pfc,
				    prob_direct_to, prob_lobe_dn_to, prob_lobe_up_to);

      // Compute the beam fluxes: the direct beam and the forward lobe
      calc_beam_fluxes(band.mu_sun, 
		       band.mu_sun, band.mu_sun_lobe_dn, band.mu_sun_lobe_up,
		       band.sec_sun, band.sec_sun_lobe_dn, band.sec_sun_lobe_up,
		       band.frac_sun_lobe_dn, prop_prime, prop_delta_ed,
		       prob_direct_to, prob_lobe_dn_to, prob_lobe_up_to,
		       angular_variance_lobe,
		       beam_fluxes);

      // FIX remove
      //      prop_delta_ed.ssa = 0.999;

      // Calculate diffuse fluxes
      meador_weaver(prop_delta_ed, beam_fluxes, diffuse_props);

      // Surface source is simply the direct and lobe downwelling
      // fluxes that are reflected up into the diffuse upwelling
      // stream at the surface. 

      //      std::cerr << "albedo size: " << albedo.size() << "\n";
      
      if (albedo.size() > 1) {
	surface_source = albedo(1) * beam_fluxes.dn(n_z);
      }
      else {
	surface_source = albedo(0) * beam_fluxes.dn(n_z);
      }

      // Adding method to get diffuse fluxes
      adding(diffuse_props, albedo, surface_source, diffuse_fluxes);

      if (write_internals) {
	working_data->write_comment_line("4. Merged properties");
	working_data->write("merged_phase_function", pf);
	working_data->write("merged_phase_function_components", pfc);
	working_data->write("primed_optical_depth", prop_prime.od);
	working_data->write("primed_single_scattering_albedo", prop_prime.ssa);
	working_data->write("primed_asymmetry_factor", prop_prime.g);
	working_data->write("delta_eddington_optical_depth", prop_delta_ed.od);
	working_data->write("delta_eddington_single_scattering_albedo", prop_delta_ed.ssa);
	working_data->write("delta_eddington_asymmetry_factor", prop_delta_ed.g);

	working_data->write_comment_line("5. Outgoing beam properties");
	working_data->write("beam_prob_direct_to_lobe",        prob_direct_to.lobe);
	working_data->write("beam_prob_direct_to_diffuse_up",  prob_direct_to.diffuse_up);
	working_data->write("beam_prob_direct_to_diffuse_dn",  prob_direct_to.diffuse_dn);
	working_data->write("beam_prob_lobe_dn_to_lobe",       prob_lobe_dn_to.lobe);
	working_data->write("beam_prob_lobe_dn_to_diffuse_up", prob_lobe_dn_to.diffuse_up);
	working_data->write("beam_prob_lobe_dn_to_diffuse_dn", prob_lobe_dn_to.diffuse_dn);
	working_data->write("beam_prob_lobe_up_to_lobe",       prob_lobe_up_to.lobe);
	working_data->write("beam_prob_lobe_up_to_diffuse_up", prob_lobe_up_to.diffuse_up);
	working_data->write("beam_prob_lobe_up_to_diffuse_dn", prob_lobe_up_to.diffuse_dn);
	working_data->write("beam_lobe_angular_variance", angular_variance_lobe);
	working_data->write("beam_flux_direct_dn", beam_fluxes.direct_dn());
	working_data->write("beam_flux_lobe_dn", beam_fluxes.lobe_dn());
	working_data->write("beam_flux_lobe_up", beam_fluxes.lobe_up());
	working_data->write("beam_lobe_angular_variance_dn", beam_fluxes.lobe_ang_var_dn);
      }

      // Compute the probability that the direct/quasi-direct return
      // beam, and the forward lobe, are scattered from the lobe or
      // the downward or upward hemispheres. Note that we overwrite
      // the prob_* data
      calc_scattering_probabilities(band.mu_inst,
				    band.mu_inst_lobe_dn, band.mu_inst_lobe_up,
				    band.frac_inst_lobe_dn, pfc,
				    prob_direct_to, prob_lobe_dn_to, prob_lobe_up_to);

      

      // Compute the return pseudo-fluxes
      calc_beam_fluxes(1.0,
		       band.mu_inst, band.mu_inst_lobe_dn, band.mu_inst_lobe_up,
		       band.sec_inst, band.sec_inst_lobe_dn, band.sec_inst_lobe_up,
		       band.frac_inst_lobe_dn,
		       //		       prop_prime, prop_delta_ed,
		       prop_prime, prop_prime, // No scaling required here
		       prob_direct_to, prob_lobe_dn_to, prob_lobe_up_to,
		       angular_variance_lobe,
		       return_fluxes);

      /*
      calc_smoothed_phase_function<IsActive>(band.pf_components,
					     pfc, beam_to_beam);
      */
      
      calc_effective_phase_functions<IsActive>(band.pf_coeffs, pfc,
					       beam_fluxes.lobe_ang_var_dn,
					       return_fluxes.lobe_ang_var_dn,
					       band.max_angular_variance,
					       beam_to_beam);
      
      Real ref[5];

      // Add contribution to reflectance
      reflectance += band.weight[ig]
	* calc_reflectance_from_fluxes(albedo,prop_prime, beam_fluxes,
				       diffuse_fluxes, return_fluxes,
				       beam_to_beam, ref);

      // Save working data (sum over spectrum)
      if (working_data) {
	working_data->ref_direct  += ref[0];
	working_data->ref_lobe    += ref[1];
	working_data->ref_lobe_up += ref[2];
	working_data->ref_diffuse += ref[3];
	working_data->ref_surface += ref[4];
	working_data->flux_toa    += value(diffuse_fluxes.up[0]);
	if (beam_fluxes.do_lobe_up()) {
	  working_data->flux_toa  += value(beam_fluxes.lobe_up_toa());
	}
      }

      if (write_internals) {
	working_data->write_comment_line("6. Diffuse properties");
	working_data->write("surface_source", surface_source);
	working_data->write("diffuse_reflectance", diffuse_props.reflectance);
	working_data->write("diffuse_transmittance", diffuse_props.transmittance);
	working_data->write("diffuse_source_up", diffuse_props.source_up);
	working_data->write("diffuse_source_dn", diffuse_props.source_dn);

	working_data->write_comment_line("7. Diffuse fluxes");
	working_data->write("diffuse_flux_dn", diffuse_fluxes.dn);
	working_data->write("diffuse_flux_up", diffuse_fluxes.up);
	working_data->write("diffuse_total_source", diffuse_fluxes.total_source);

	working_data->write_comment_line("8. Returning beam properties");
	working_data->write("return_prob_direct_to_lobe",        prob_direct_to.lobe);
	working_data->write("return_prob_direct_to_diffuse_up",  prob_direct_to.diffuse_up);
	working_data->write("return_prob_direct_to_diffuse_dn",  prob_direct_to.diffuse_dn);
	working_data->write("return_prob_lobe_dn_to_lobe",       prob_lobe_dn_to.lobe);
	working_data->write("return_prob_lobe_dn_to_diffuse_up", prob_lobe_dn_to.diffuse_up);
	working_data->write("return_prob_lobe_dn_to_diffuse_dn", prob_lobe_dn_to.diffuse_dn);
	working_data->write("return_prob_lobe_up_to_lobe",       prob_lobe_up_to.lobe);
	working_data->write("return_prob_lobe_up_to_diffuse_up", prob_lobe_up_to.diffuse_up);
	working_data->write("return_prob_lobe_up_to_diffuse_dn", prob_lobe_up_to.diffuse_dn);
	working_data->write("return_lobe_angular_variance", angular_variance_lobe);
	working_data->write("return_flux_direct_dn", return_fluxes.direct_dn());
	working_data->write("return_flux_lobe_dn", return_fluxes.lobe_dn());
	working_data->write("return_flux_lobe_up", return_fluxes.lobe_up());
	working_data->write("return_lobe_angular_variance_dn", return_fluxes.lobe_ang_var_dn);
	working_data->write("beam_to_beam", beam_to_beam);
      }

    }

    if (write_internals) {
      working_data->write_comment_line("9. Reflectances");
      working_data->write("reflectance", reflectance);
      working_data->write("reflectance_direct", working_data->ref_direct);
      working_data->write("reflectance_lobe", working_data->ref_lobe);
      working_data->write("reflectance_lobe_up", working_data->ref_lobe_up);
      working_data->write("reflectance_diffuse", working_data->ref_diffuse);
      working_data->write("reflectance_surface", working_data->ref_surface);
      working_data->write("flux_toa", working_data->flux_toa);
    }

    return FLOTSAM_SUCCESS;
  }


  // Explicit instantiations
  template 
  int reflectance<false>(const BandProfile& band,
		  const Array<1,Real,false>& albedo,
		  const intVector& loc_particulate,
		  const Array<1,Real,false>& od_particulate, 
		  const Array<1,Real,false>& ssa_particulate,
		  const Array<1,Real,false>& pf_particulate,
		  const Array<2,Real,false>& pfc_particulate,
		  scalar<false>::type& reflectance,
		  WorkingData* working_data ///< Store intermediate data if not NULL
		  );

  template
  int reflectance<true>(const BandProfile& band,
		  const Array<1,Real,true>& albedo,
		  const intVector& loc_particulate,
		  const Array<1,Real,true>& od_particulate, 
		  const Array<1,Real,true>& ssa_particulate,
		  const Array<1,Real,true>& pf_particulate,
		  const Array<2,Real,true>& pfc_particulate,
		  scalar<true>::type& reflectance,
		  WorkingData* working_data ///< Store intermediate data if not NULL
		  );

};


