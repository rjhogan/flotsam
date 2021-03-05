/** @file      calc_smoothed_phase_function.hpp
    @brief     Compute smoothed phase function
    @copyright 2017 European Centre for Medium Range Weather Forcasts
    @license   Apache License Version 2 (see the NOTICE.md file for details)
 */
#ifndef FLOTSAM_CALC_SMOOTHED_PHASE_FUNCTION_HPP
#define FLOTSAM_CALC_SMOOTHED_PHASE_FUNCTION_HPP

namespace flotsam {

  // template <bool IsActive>
  // static
  // void calc_powers(const Array<1,Real,IsActive>& in,
  // 		   Array<2,Real,IsActive>& out) {
  //   //    out(0,__) = 1.0;
  //   //    out(1,__) = in;
  //   //    out(2,__) = in*in;
  //   //    out(3,__) = out(2,__)*in;
  //   out(__,0) = 1.0;
  //   out(__,1) = in;
  //   out(__,2) = in*in;
  //   out(__,3) = out(__,2)*in;
  // }
#define calc_powers(in, out)			\
  out(__,0) = 1.0;				\
  out(__,1) = (in);				\
  out(__,2) = (in)*(in);			\
  out(__,3) = out(__,2)*(in) 


  // template <bool IsActive>
  // static
  // void
  // calc_effective_phase_function(const Array<2,Real,IsActive>& pfc,
  // 				const Matrix& pf_coeffs,
  // 				const Array<2,Real,IsActive>& ang_var_powers,
  // 				Array<1,Real,IsActive> out) {
  //   out = sum(pfc, pf_coeffs, ang_var_powers);
  // }
#define effective_phase_function(pfc, pf_coeffs, ang_var_powers) \
  sum((pfc ** pf_coeffs) * ang_var_powers, 1)

  template <bool IsActive>
  void calc_effective_phase_functions(const Array4D& pf_coeffs,
				      const Array<2,Real,IsActive>& pfc,
				      const Array<1,Real,IsActive>& beam_ang_var_dn,
				      const Array<1,Real,IsActive>& ret_ang_var_dn,
				      Real max_angular_variance,
				      Array<3,Real,IsActive>& b2b) {
    typedef Array<1,Real,IsActive> avector;
    typedef Array<2,Real,IsActive> amatrix;
    typedef Array<3,Real,IsActive> aarray3D;
    //    amatrix ang_var_powers(FLOTSAM_NUM_SMOOTHING_COEFFS, beam_ang_var_dn.size());
    avector ang_var;
    amatrix ang_var_powers(beam_ang_var_dn.size(), FLOTSAM_NUM_SMOOTHING_COEFFS);
    aarray3D pf_components;
    ang_var = fmin(ret_ang_var_dn, max_angular_variance);
    calc_powers(ang_var, ang_var_powers);
    b2b(0,1,__) = effective_phase_function(pfc, pf_coeffs(0,1,__,__), ang_var_powers);
    b2b(0,2,__) = effective_phase_function(pfc, pf_coeffs(0,2,__,__), ang_var_powers);
    ang_var = fmin(beam_ang_var_dn, max_angular_variance);
    calc_powers(ang_var, ang_var_powers);
    b2b(1,0,__) = effective_phase_function(pfc, pf_coeffs(1,0,__,__), ang_var_powers);
    b2b(2,0,__) = effective_phase_function(pfc, pf_coeffs(2,0,__,__), ang_var_powers);
    ang_var = fmin(beam_ang_var_dn+ret_ang_var_dn, max_angular_variance);
    calc_powers(ang_var, ang_var_powers);
    b2b(1,1,__) = effective_phase_function(pfc, pf_coeffs(1,1,__,__), ang_var_powers);
    b2b(1,2,__) = effective_phase_function(pfc, pf_coeffs(1,2,__,__), ang_var_powers);
    b2b(2,1,__) = effective_phase_function(pfc, pf_coeffs(2,1,__,__), ang_var_powers);
    b2b(2,2,__) = effective_phase_function(pfc, pf_coeffs(2,2,__,__), ang_var_powers);
  }

  template <bool IsActive>
  void calc_smoothed_phase_function(const Array3D& pf_components,
				    const Array<2,Real,IsActive>& pfc,
				    Array<3,Real,IsActive>& b2b) {
    b2b(0,1,__) = pfc ** pf_components(0,1,__); // Sun direct to inst lobe-dn
    b2b(0,2,__) = pfc ** pf_components(0,2,__); // Sun direct to inst lobe-up
    b2b(1,0,__) = pfc ** pf_components(1,0,__); // Sun lobe-dn to inst direct
    b2b(1,1,__) = pfc ** pf_components(1,1,__); // Sun lobe-dn to inst lobe-dn
    b2b(1,2,__) = pfc ** pf_components(1,2,__); // Sun lobe-dn to inst lobe-up
    b2b(2,0,__) = pfc ** pf_components(2,0,__); // Sun lobe-up to inst direct
    b2b(2,1,__) = pfc ** pf_components(2,1,__); // Sun lobe-up to inst lobe-dn
    b2b(2,2,__) = pfc ** pf_components(2,2,__); // Sun lobe-up to inst lobe-up
  }

}

#endif
