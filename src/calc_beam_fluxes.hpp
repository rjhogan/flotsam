/** @file      calc_beam_fluxes.hpp
    @brief     Compute profile of beam fluxes (direct and lobe)
    @copyright 2017 European Centre for Medium Range Weather Forcasts
    @license   Apache License Version 2 (see the NOTICE.md file for details)
 */
#ifndef FLOTSAM_CALC_BEAM_FLUXES_HPP
#define FLOTSAM_CALC_BEAM_FLUXES_HPP

#include "containers.hpp"

namespace flotsam {

  /// Compute the direct and lobe flux profiles, describing them in
  /// terms of exponential functions
  template <bool IsActive>
  void calc_beam_fluxes(Real direct_toa,
			Real mu_direct, Real mu_lobe_dn, Real mu_lobe_up,
			Real sec_direct, Real sec_lobe_dn, Real sec_lobe_up,
			Real frac_lobe_dn,
			const ScatteringProperties<IsActive>& prop_prime, 
			const ScatteringProperties<IsActive>& prop_delta_ed, 
			const ScatteringProbabilities<IsActive>& prob_direct_to,
			const ScatteringProbabilities<IsActive>& prob_lobe_dn_to,
			const ScatteringProbabilities<IsActive>& prob_lobe_up_to,
			const Array<1,Real,IsActive>& angular_var_lobe,
			BeamFluxes<IsActive>& beam) {

    typedef typename scalar<IsActive>::type areal;
    typedef Array<1,Real,IsActive> avector;

    const Real frac_lobe_up = 1.0 - frac_lobe_dn;
    const int nz = prop_prime.size();

    // Copy secants over
#ifndef FLOTSAM_REVISION5
    beam.sec(FLOTSAM_DIRECT)  = sec_direct;
    beam.sec(FLOTSAM_LOBE_DN_1) = sec_lobe_dn;
    beam.sec(FLOTSAM_LOBE_DN_2) = sec_lobe_dn;
    if (beam.do_lobe_up()) {
      beam.sec(FLOTSAM_LOBE_UP_1) = sec_lobe_up;
      beam.sec(FLOTSAM_LOBE_UP_2) = sec_lobe_up;
    }
#else
    beam.sec(FLOTSAM_DIRECT)  = 1.0 / mu_direct;
    beam.sec(FLOTSAM_LOBE_DN_1) = 1.0 / mu_lobe_dn;
    beam.sec(FLOTSAM_LOBE_DN_2) = 1.0 / mu_lobe_dn;
    if (beam.do_lobe_up()) {
      beam.sec(FLOTSAM_LOBE_UP_1) = 1.0 / mu_lobe_up;
      beam.sec(FLOTSAM_LOBE_UP_2) = 1.0 / mu_lobe_up;
    }
#endif

    // TOA downwelling solar incoming radiation
    beam.prefactor_direct(0) = direct_toa; //1.0/sec_direct;

    beam.exponent_direct    = sec_direct*prop_prime.od;
    // Cache the results of exp calls for use later on
    beam.expm_exponent_direct = exp(-beam.exponent_direct);
    beam.exponent_lobe_dn_1 = beam.exponent_direct;
    beam.expm_exponent_lobe_dn_1 = beam.expm_exponent_direct;
    if (beam.do_lobe_up()) {
      beam.exponent_lobe_up_1 = beam.exponent_direct;
      beam.expm_exponent_lobe_up_1 = beam.expm_exponent_direct;
    }

    areal lobe = 0.0; // Strenth of lobe_dn at layer interface,
		      // starting at zero at TOA
    areal lobe_new;
    areal lobe_var = 0.0; // Lobe strength x angular variance
    areal lobe_var_new;

    areal lobe_var_prefactor1, lobe_var_prefactor2;
    areal lobe_var_prefactor3, lobe_var_preprefactor;

    areal lobe_exp;

    // Proceed down from TOA calculating the direct beam flux and the
    // downwelling lobe flux
    for (int i = 0; i < nz; ++i) {
      /*
      if (i == 0) {
	std::cout << "CHECK " << sec_direct << " "
		  << prob_direct_to.lobe(i) << " "
		  << frac_lobe_dn << " "
		  << prop_prime.ssa(i)
		  << beam.prefactor_direct(i) << " "
		  << sec_lobe_dn << " "
		  << prob_lobe_dn_to.lobe(i) << "\n";
      }
      */
      lobe_exp = (1.0 - prob_lobe_dn_to.lobe(i)*prop_prime.ssa(i))*sec_lobe_dn;

      using std::abs;
      if (abs(lobe_exp-sec_direct) > 1.0e-8) {
#ifndef FLOTSAM_REVISION2
      beam.prefactor_lobe_dn_1(i)
	= sec_direct*prob_direct_to.lobe(i)*frac_lobe_dn
	* prop_prime.ssa(i)*beam.prefactor_direct(i)
	/ (lobe_exp - sec_direct);
#else
	beam.prefactor_lobe_dn_1(i)
	  = sec_lobe_dn*prob_direct_to.lobe(i)*frac_lobe_dn
	  * prop_prime.ssa(i)*beam.prefactor_direct(i)
	  / (lobe_exp - sec_direct);
#endif
      }
      else {
	beam.prefactor_lobe_dn_1(i) = 0.0;
      }

      beam.prefactor_lobe_dn_2(i) 
	= lobe - beam.prefactor_lobe_dn_1(i);

      beam.exponent_lobe_dn_2(i) = lobe_exp*prop_prime.od(i);

      // Cache the results of exp calls for use later on
      beam.expm_exponent_lobe_dn_2(i) = exp(-beam.exponent_lobe_dn_2(i));
      
      // Model the variance
      //lobe_var_source_prefactor
      //= sec_direct*prop_prime.ssa(i)*(prob_direct_to.lobe(i)*frac_lobe_dn*beam.prefactor_direct(i)
      //+prob_lobe_dn_to.lobe(i)*lobe);
#ifndef FLOTSAM_REVISION2
      lobe_var_preprefactor = sec_lobe_dn*prop_prime.ssa(i)*angular_var_lobe(i);
#else
      lobe_var_preprefactor = sec_direct*prop_prime.ssa(i)*angular_var_lobe(i);
#endif
      if (abs(lobe_exp-sec_direct) > 1.0e-8) {
	lobe_var_prefactor1
	  = lobe_var_preprefactor*(frac_lobe_dn*prob_direct_to.lobe(i)*beam.prefactor_direct(i)
				   +prob_lobe_dn_to.lobe(i)*beam.prefactor_lobe_dn_1(i))
	  / (lobe_exp - sec_direct);
      }
      else {
	lobe_var_prefactor1 = 0.0;
      }
      lobe_var_prefactor2
	= lobe_var_preprefactor*(prob_lobe_dn_to.lobe(i)*beam.prefactor_lobe_dn_2(i));
      lobe_var_prefactor3 = lobe_var - lobe_var_prefactor1;
      lobe_var_new = lobe_var_prefactor1 * beam.expm_exponent_direct(i)
	+ (lobe_var_prefactor2*prop_prime.od(i) + lobe_var_prefactor3)
	* beam.expm_exponent_lobe_dn_2(i);

      // lobe_var_sum = lobe_var_prefactor;
      // std::cout << "PREL1= " << beam.prefactor_lobe_dn_1(i) << "\n";
      // std::cout << "PRE1 = " << lobe_var_prefactor << "\n";
      // lobe_var_new = lobe_var_prefactor * beam.expm_exponent_direct(i);
      
      // lobe_var_prefactor
      // 	= lobe_var_preprefactor*prob_lobe_dn_to.lobe(i)*beam.prefactor_lobe_dn_1(i)
      // 	/ (lobe_exp - beam.exponent_lobe_dn_1(i));
      // lobe_var_sum += lobe_var_prefactor;
      // std::cout << "PRE1 = " << lobe_var_prefactor << "\n";
      // lobe_var_new += lobe_var_prefactor * beam.expm_exponent_lobe_dn_1(i);
      
      // lobe_var_prefactor
      // 	= lobe_var_preprefactor*prob_lobe_dn_to.lobe(i)*beam.prefactor_lobe_dn_2(i)
      // 	/ (lobe_exp - beam.exponent_lobe_dn_2(i));
      // lobe_var_sum += lobe_var_prefactor;
      // std::cout << "PRE1 = " << lobe_var_prefactor << "\n";
      // lobe_var_new += lobe_var_prefactor * beam.expm_exponent_lobe_dn_2(i);
      
      // lobe_var_new += (lobe_var - lobe_var_sum) * beam.expm_exponent_lobe_dn_2(i);

      lobe_new
	= beam.prefactor_lobe_dn_1(i)*beam.expm_exponent_lobe_dn_1(i)
	+ beam.prefactor_lobe_dn_2(i)*beam.expm_exponent_lobe_dn_2(i);

      // Layer-mean angular variance is from the weighted-average of
      // the variances at the boundaries
      using std::max;
      beam.lobe_ang_var_dn(i) = (lobe_var + lobe_var_new)/max(1.0e-12,lobe + lobe_new);
      if (i < nz-1) {
	beam.prefactor_direct(i+1)
	  = beam.prefactor_direct(i) * beam.expm_exponent_direct(i);
	lobe = lobe_new;
	lobe_var = lobe_var_new;
      }
    }

    if (beam.do_lobe_up()) {
      lobe = 0.0; // Strength of lobe_up at the surface
      // Proceed up from surface calculating the upwelling lobe flux
      for (int i = nz-1; i >= 0; --i) {
#ifndef FLOTSAM_REVISION2
	beam.prefactor_lobe_up_1(i) = sec_direct*(prob_direct_to.lobe(i)*frac_lobe_up
						  *prop_prime.ssa(i)*beam.prefactor_direct(i))
	  / (sec_lobe_up*(1.0 - prob_lobe_up_to.lobe(i)*prop_prime.ssa(i)) + sec_direct);
#else
	beam.prefactor_lobe_up_1(i) = sec_lobe_up*(prob_direct_to.lobe(i)*frac_lobe_up
						  *prop_prime.ssa(i)*beam.prefactor_direct(i))
	  / (sec_lobe_up*(1.0 - prob_lobe_up_to.lobe(i)*prop_prime.ssa(i)) + sec_direct);
#endif
	beam.exponent_lobe_up_2(i)
	  = (prob_lobe_up_to.lobe(i)*prop_prime.ssa(i) - 1.0)*sec_lobe_up*prop_prime.od(i);
	beam.expm_exponent_lobe_up_2(i) = exp(-beam.exponent_lobe_up_2(i));
	beam.prefactor_lobe_up_2(i)
	  = (lobe - beam.prefactor_lobe_up_1(i)*beam.expm_exponent_lobe_up_1(i))
	  / beam.expm_exponent_lobe_up_2(i);
	lobe = beam.prefactor_lobe_up_1(i) + beam.prefactor_lobe_up_2(i);
      }
    }

    // Calculate source terms for diffuse fluxes

    // Firstly, Meador-Weaver uses delta-Eddington optical depth as
    // vertical coordinate so we need to change the rate of scattering
    avector factor(nz);
    factor = 1.0;
    factor.where(prop_delta_ed.od > 0.0) = prop_prime.od / prop_delta_ed.od;

    // FIX
    /*
    sec_direct  = 1.0 / mu_direct;
    sec_lobe_dn = 1.0 / mu_direct;
    sec_lobe_up = 1.0 / mu_direct;
    */

    beam.prefactor_source_dn(FLOTSAM_DIRECT,__)
      = sec_direct*factor*beam.prefactor_direct*prob_direct_to.diffuse_dn;
    beam.prefactor_source_dn(FLOTSAM_LOBE_DN_1,__)
      = sec_lobe_dn*factor*beam.prefactor_lobe_dn_1*prob_lobe_dn_to.diffuse_dn;
    beam.prefactor_source_dn(FLOTSAM_LOBE_DN_2,__)
      = sec_lobe_dn*factor*beam.prefactor_lobe_dn_2*prob_lobe_dn_to.diffuse_dn;

    beam.prefactor_source_up(FLOTSAM_DIRECT,__)
      = sec_direct*factor*beam.prefactor_direct*prob_direct_to.diffuse_up;
    beam.prefactor_source_up(FLOTSAM_LOBE_DN_1,__)
      = sec_lobe_dn*factor*beam.prefactor_lobe_dn_1*prob_lobe_dn_to.diffuse_up;
    beam.prefactor_source_up(FLOTSAM_LOBE_DN_2,__)
      = sec_lobe_dn*factor*beam.prefactor_lobe_dn_2*prob_lobe_dn_to.diffuse_up;

    if (beam.do_lobe_up()) {
      beam.prefactor_source_dn(FLOTSAM_LOBE_UP_1,__)
	= sec_lobe_up*factor*beam.prefactor_lobe_up_1*prob_lobe_up_to.diffuse_dn;
      beam.prefactor_source_dn(FLOTSAM_LOBE_UP_2,__)
	= sec_lobe_up*factor*beam.prefactor_lobe_up_2*prob_lobe_up_to.diffuse_dn;
      
      beam.prefactor_source_up(FLOTSAM_LOBE_UP_1,__)
	= sec_lobe_up*factor*beam.prefactor_lobe_up_1*prob_lobe_up_to.diffuse_up;
      beam.prefactor_source_up(FLOTSAM_LOBE_UP_2,__)
	= sec_lobe_up*factor*beam.prefactor_lobe_up_2*prob_lobe_up_to.diffuse_up;
    }

    // Check for zero optical depth, which leads to zero exponent,
    // which leads to divide by zero
    intVector v = find(prop_prime.od <= 0.0);
    if (!v.empty()) {
      beam.exponent(__,v) = 1.0e-16;
    }

  }
}

#endif
