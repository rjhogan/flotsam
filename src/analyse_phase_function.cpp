/** @file      analyse_phase_functions.cpp
    @brief     Analyse phase functions to use in FLOTSAM
    @copyright 2017 European Centre for Medium Range Weather Forcasts
    @license   Apache License Version 2 (see the NOTICE.md file for details)
 */

#include <flotsam.hpp>

#include "LookUpTable.hpp"
#include "phase_function_component.hpp"
#include "asymmetry_factor_from_phase_func.hpp"
#include "calc_asymmetry_factor.hpp"

/// Compute the integral of the phase function over a sphere, where
/// the input phase functions elements are assumed to be equally
/// spaced in angle between 0 and 180 degrees.
adept::Real
flotsam::integrate_phase_function(const adept::Vector& pf) ///< Phase function in
{
  using namespace adept;

  int nang  = pf.size();
  Real dang = M_PI / (nang - 1.0);
  Vector cos_mid_ang = cos(linspace(-0.5*dang, M_PI+0.5*dang, nang+1));
  // Set first value to cos(0);
  cos_mid_ang(0) = 1.0;
  // Set last value to cos(pi);
  cos_mid_ang(nang) = -1.0;
  // Area of the surface of a unit sphere corresponding to each
  // range of scattering angle
  Vector weight = -2.0*M_PI*(cos_mid_ang(range(1,end))
			   - cos_mid_ang(range(0,end-1)));
  return dot_product(pf, weight);
}

/// Process raw phase functions [phase_func,angle] to compute the
/// effective single-scattering phase function after the delta peak
/// has been removed, and the amplitudes of the components of the
/// phase-function decomposition. The phase functions elements are
/// assumed to be equally spaced in angle between 0 and 180 degrees.
void
flotsam::analyse_phase_functions(const adept::Matrix& pf,      ///< Phase functions in
				 adept::Real normalization,    ///< Integral over surface of sphere
				 adept::Matrix& pf_out,        ///< Phase functions out
				 adept::Matrix& pf_components) ///< Component amplitudes out
{
  using namespace adept;

  static const Real four_pi = 4.0*M_PI;
  // When fitting the non-delta part of the phase function, ignore
  // 0-15 degrees, and linearly ramp down the weight given between 30
  // and 15 degrees
  static const int max_delta_angle_deg = 30;
  static const int start_rampdown_deg  = 15;
  int npf   = pf.size(0);
  int nang  = pf.size(1);
  Real dang = M_PI / (nang - 1.0);
  Vector cos_mid_ang = cos(linspace(-0.5*dang, M_PI+0.5*dang, nang+1));
  // Set first value to cos(0);
  cos_mid_ang(0) = 1.0;
  // Set last value to cos(pi);
  cos_mid_ang(nang) = -1.0;
  // Area of the surface of a unit sphere corresponding to each
  // range of scattering angle
  Vector weight = -2.0*M_PI*(cos_mid_ang(range(1,end))
			   - cos_mid_ang(range(0,end-1)));

  // We only fit to scattering angles greater than 30 degrees
  Vector weight_truncated;
  weight_truncated = weight;

  int i_max_delta = (nang*max_delta_angle_deg) / 180;
  int i_rampdown  = (nang*start_rampdown_deg)  / 180;

  weight_truncated(range(0,i_rampdown)) = 0.0;
  weight_truncated(range(i_rampdown+1,i_max_delta))
    *= linspace(0.0, 1.0, i_max_delta-i_rampdown);
  DiagMatrix weight_matrix = weight_truncated.diag_matrix();

  Matrix jacobian(nang,lut.n_phase_function_components());
  for (int i = 0; i < lut.n_phase_function_components(); ++i) {
    jacobian(__,i) = interp(linspace(0.0,M_PI,lut.nang),
			    lut.phase_function_components(__,i,0),
			    linspace(0.0,M_PI,nang));
  }

  Matrix weighted_jacobian = jacobian.T() ** weight_matrix;

  // Compute the amplitudes of the phase function components used to
  // approximate the original phase function
  pf_components.resize(npf, lut.n_phase_function_components());
  pf_components = 0.0;

  // We may need to modify first element of phase function to ensure
  // the integral equals the normalization value
  pf_out.resize(npf, nang);
  pf_out = pf;

  Vector integral = pf_out ** weight;

  for (int i = 0; i < npf; ++i) {
    if (integral(i) < normalization) {
      pf_out(i,0) += (normalization-integral(i)) / weight(0);
      integral(i) = normalization;
    }
    // Scale the phase function
    pf_out(i,__) *= (four_pi / integral(i));

    // Least squares fit to obtain phase function components
    pf_components(i,__) = solve(weighted_jacobian ** jacobian,
				weighted_jacobian ** pf_out(i,__));

    Real delta_component = 1.0 - sum(pf_components(i,lut.mass_component_index));
    if (delta_component < 1.0e-3) {
      // Delta component slightly negative: correct this by adjusting
      // the lobe component
      if (lut.i_lobe_component >= 0) {
	pf_components(i,lut.i_lobe_component) += delta_component;
      }
    }
    else {
      // Characterize delta component
      Vector pf_fit = jacobian(range(0,i_max_delta),__) ** pf_components(i,__);
      pf_fit -= pf_fit(i_max_delta);
      Real sum_fit
	= sum(weight(range(0,i_max_delta))*pf_fit);
      Real sum_delta
	= sum(weight(range(0,i_max_delta))
	      *(pf_out(i,range(0,i_max_delta))-pf_out(i,i_max_delta)));
      Vector pf_new = pf_out(i,i_max_delta)
	+ pf_fit * (sum_delta - delta_component*four_pi) / sum_fit;

      Real sum_pf = 0.0;
      Real sum_pf_ang2 = 0.0;
      for (int j = 0; j < (nang*max_delta_angle_deg)/180; ++j) {
	Real pf_diff = pf_out(i,j)-pf_new(j);
	if (pf_diff > 0.0) {
	  Real ang2 = j*M_PI/(nang-1.0);
	  ang2 *= ang2;
	  sum_pf += weight(j)*pf_diff;
	  sum_pf_ang2 += weight(j)*pf_diff*ang2;
	}
      }

      // Subtract delta component
      pf_out(i,range(0,i_max_delta)) = pf_new;

      // Scale the phase function associated with removal of the delta
      // component
      Real scaling = 1.0 / (1.0-delta_component);
      pf_out(i,__) *= scaling;

      // Angular variance of the delta component
      Real var_delta = sum_pf_ang2 / sum_pf;

      if (var_delta > 0.0) {

	Vector pf_smooth(nang);
	Vector ang = linspace(0.0, M_PI, nang);
	Vector kernel(nang);
      
	for (int k = 0; k < nang; ++k) {
	  kernel = weight * exp(-(ang-ang(k))*(ang-ang(k))/var_delta);
	  pf_smooth(k) = sum(kernel*pf_out(i,__)) / sum(kernel);
	}
#ifndef FLOTSAM_REVISION
	pf_out(i,__) = 0.5*(pf_out(i,__) + pf_smooth);
#else
	// A larger delta component means that the effective phase
	// function is smoothed more, hence the following weighting
	pf_out(i,__) = (1.0-delta_component)*pf_out(i,__)
	  + delta_component*pf_smooth;
#endif
	
      }
    }
  }

  /*
  // Check asymmetry factor of phase functions
  std::cerr << "g_true = " << asymmetry_factor_from_phase_func(pf) << "\n";
  Vector g_fit;
  calc_asymmetry_factor(pf_components, g_fit);
  std::cerr << "g_fit  = " << g_fit << "\n";
  */
}

/*
/// Original method to process raw phase functions [phase_func,angle]
/// to compute the effective single-scattering phase function after
/// the delta peak has been removed, and the amplitudes of the
/// components of the phase-function decomposition. The phase
/// functions elements are assumed to be equally spaced in angle
/// between 0 and 180 degrees.
void
flotsam::analyse_phase_functions_orig(const adept::Matrix& pf, ///< Phase functions in
				 adept::Real normalization,    ///< Integral over surface of sphere
				 adept::Matrix& pf_out,        ///< Phase functions out
				 adept::Matrix& pf_components) ///< Component amplitudes out
{
  using namespace adept;

  static const Real four_pi = 4.0*M_PI;
  int npf   = pf.size(0);
  int nang  = pf.size(1);
  Real dang = M_PI / (nang - 1.0);
  Vector cos_mid_ang = cos(linspace(-0.5*dang, M_PI+0.5*dang, nang+1));
  // Set first value to cos(0);
  cos_mid_ang(0) = 1.0;
  // Set last value to cos(pi);
  cos_mid_ang(nang) = -1.0;
  // Area of the surface of a unit sphere corresponding to each
  // range of scattering angle
  Vector weight = -2.0*M_PI*(cos_mid_ang(range(1,end))
			   - cos_mid_ang(range(0,end-1)));
  // Weight of backward hemisphere
  Vector weight_back_hem;
  weight_back_hem = weight;
  if (nang % 2) {
    // Odd number of angles
    int imid = (nang-1)/2;
    weight_back_hem(range(0,imid-1)) = 0.0;
    weight_back_hem(imid) *= 0.5;
  }
  else {
    // Even number of angles
    weight_back_hem(range(0,nang/2-1)) = 0.0;
  }

  // We may need to modify first element of phase function to ensure
  // the integral equals the normalization value
  pf_out.resize(npf, nang);
  pf_out = pf;

  Vector integral = pf_out ** weight;

  for (int i = 0; i < npf; ++i) {
    pf_out(i,0) += (normalization-integral(i)) / weight(0);
    integral(i) = normalization;
  }

  Vector frac_back_hem = (pf_out ** weight_back_hem) / integral;
  Vector g = asymmetry_factor_from_phase_func(pf_out);

  Vector g_back_hem = asymmetry_factor_from_phase_func(pf_out, true);

  // Compute the amplitudes of the phase function components used to
  // approximate the original phase function
  pf_components.resize(npf, lut.n_phase_function_components());
  pf_components = 0.0;

  // The isotropic and "Rayleigh" components are chosen such that in
  // the backward hemisphere only, the two of them have the same total
  // scattering and the same asymmetry factor as the original phase
  // function (in the backward hemisphere)
  pf_components(__, FLOTSAM_PHASE_FUNC_ISOTROPIC)
    = max(0.0,
	  2.0*frac_back_hem
	  * (g_back_hem-lut.g_back_hem_components(FLOTSAM_PHASE_FUNC_RAYLEIGH))
	  / (lut.g_back_hem_components(FLOTSAM_PHASE_FUNC_ISOTROPIC)
	     -lut.g_back_hem_components(FLOTSAM_PHASE_FUNC_RAYLEIGH)));
  pf_components(__, FLOTSAM_PHASE_FUNC_RAYLEIGH) 
    = 2.0*frac_back_hem - pf_components(__, FLOTSAM_PHASE_FUNC_ISOTROPIC);
  // The convex lobe component is chosen such that the four components
  // (the implicit one being the delta function that is a remainder of
  // the others) combine to give the same asymmetry factor as the
  // original phase function
  pf_components(__, FLOTSAM_PHASE_FUNC_CONVEX_LOBE)
    = max(0.0,
	  (g - 1.0
	   + pf_components(__, FLOTSAM_PHASE_FUNC_ISOTROPIC)
	   + pf_components(__, FLOTSAM_PHASE_FUNC_RAYLEIGH))
	  / (lut.g_components(FLOTSAM_PHASE_FUNC_CONVEX_LOBE) - 1.0));
  pf_components(__, FLOTSAM_PHASE_FUNC_CONVEX_LOBE)
    = min(pf_components(__, FLOTSAM_PHASE_FUNC_CONVEX_LOBE),
	  1.0 - pf_components(__, FLOTSAM_PHASE_FUNC_ISOTROPIC)
	  - pf_components(__, FLOTSAM_PHASE_FUNC_RAYLEIGH));

  // FIX FIX FIX check if no delta scaling...
  //pf_components(__, FLOTSAM_PHASE_FUNC_CONVEX_LOBE) = 1.0 - pf_components(__, FLOTSAM_PHASE_FUNC_ISOTROPIC) - pf_components(__, FLOTSAM_PHASE_FUNC_RAYLEIGH);

  Vector delta_component = 1.0 - sum(pf_components,1);
  
  // Create the phase functions with lobes removed

  for (int i = 0; i < npf; ++i) {
    if (delta_component(i) < 1.0e-3) {
      // No delta component: simply copy the phase function over with
      // the required normalization (integrates over surface of a
      // sphere to 4pi)
      pf_out(i,__) *= (four_pi / integral(i));
    }
    else {
      // Characterize delta component
      Real sum_pf = 0.0;
      Real sum_flat = 0.0;
      Real sum_pf_ang2 = 0.0;
      Real sum_flat_ang2 = 0.0;
      int j = 0;
      for ( ; j < nang; ++j) {
	Real ang2 = j*M_PI/(nang-1.0);
	ang2 *= ang2;
	sum_pf += weight(j)*pf_out(i,j);
	sum_flat += weight(j);
	sum_pf_ang2 += weight(j)*pf_out(i,j)*ang2;
	sum_flat_ang2 += weight(j)*ang2;

	if (sum_pf - sum_flat*pf_out(i,j+1)
	    > delta_component(i) * integral(i)) {
	  // We have reached the limit of the delta component
	  break;
	}
      }
      // Scale the phase function
      Real scaling = four_pi / (integral(i)*(1.0-delta_component(i)));
      pf_out(i,__) *= scaling;

      // Value of phase function plateau
      Real pf_flat = (sum_pf - delta_component(i) * integral(i)) / sum_flat;

      // Flatten the delta peak
      pf_out(i,range(0,j)) = pf_flat*scaling;

      // Angular variance of the delta component
      Real var_delta = (sum_pf_ang2 - sum_flat_ang2*pf_flat) / (sum_pf - sum_flat*pf_flat);
      //Real var_delta = sum_pf_ang2 / sum_pf;

      Vector pf_smooth(nang);
      Vector ang = linspace(0.0, M_PI, nang);
      Vector kernel(nang);
      
      for (int k = 0; k < nang; ++k) {
	kernel = weight * exp(-(ang-ang(k))*(ang-ang(k))/var_delta);
	pf_smooth(k) = sum(kernel*pf_out(i,__)) / sum(kernel);
      }
      //      std::cout << pf_out(i,__);
#ifndef FLOTSAM_REVISION
      pf_out(i,__) = 0.5*(pf_out(i,__) + pf_smooth);
#else
      // A larger delta component means that the effective phase
      // function is smoothed more, hence the following weighting
      pf_out(i,__) = (1.0-delta_component(i))*pf_out(i,__)
	+ delta_component(i)*pf_smooth;
#endif
    }
  }
}
*/


/// Process a single phase function to compute the effective
/// single-scattering phase function after the delta peak has been
/// removed, and the amplitudes of the components of the
/// phase-function decomposition
void
flotsam::analyse_phase_function(const adept::Vector& pf,      ///< Phase function in
				adept::Real normalization,    ///< Integral over surface of sphere
				adept::Vector& pf_out,        ///< Phase function out
				adept::Vector& pf_components) ///< Component amplitudes out
{
  using namespace adept;
  int nang = pf.size();
  Matrix pf_mat(1,nang);
  Matrix pf_mat_out(1,nang);
  Matrix pf_mat_components(1,nang);
  pf_mat(0,__) = pf;
  flotsam::analyse_phase_functions(pf_mat, normalization, pf_mat_out, pf_mat_components);
  pf_out = pf_mat_out(0,__);
  pf_components = pf_mat_components(0,__);
}

// Need Adept's "spread" function
#if ADEPT_VERSION >= 10910

/// Spherical convolution of phase function "pf" (linearly spaced in
/// angle between 0 and 180 degrees) with the phase function component
/// indicated with "icomponent".  A negative number (the default if
/// this argument is omitted) uses the default smoothing kernel,
/// otherwise you can use one of the basis functions used to
/// parameterize phase functions.
adept::Vector
flotsam::convolve_phase_function(const adept::Vector& pf, ///< Phase function in
				 int icomponent)          ///< Index of component
{
  using namespace adept;

  if (icomponent >= lut.n_phase_function_components()) {
    // Impossible value for icomponent
    return false;
  }

  int nang = pf.size();
  // Scattering angle (radians)
  Vector ang = linspace(0,M_PI,nang);
  Real dang = ang(1);
  // Longitude (radians)
  Vector lon = ang-M_PI;
  int nlon = lon.size();

  Matrix CosAng = spread<0>(cos(ang),nlon);
  Vector weighted_pf(CosAng.size());
  // 2D array linking to 1D data
  Matrix WeightedPF = weighted_pf.reshape(CosAng.size(0), CosAng.size(1));

  bool has_zero_integral = false;
  if (sum(pf) < 0.01) {
    // For phase function components that integrate to zero, we ensure
    // the correct normalization of the convolved function by adding
    // one, convolving, normalizing, then subtracting one.
    has_zero_integral = true;
    WeightedPF = spread<0>((pf+1.0)*sin(ang)*(dang*dang/M_PI),nlon);
  }
  else {
    WeightedPF = spread<0>(pf*sin(ang)*(dang*dang/M_PI),nlon);
  }

  Matrix SinAngCosLon = spread<0>(sin(ang),nlon) * spread<1>(cos(lon),nang);
  //  Matrix SinAngCosLon = outer_product(cos(lon),sin(ang));

  // Smoothing kernel phase function
  Vector pf_kernel;
  if (icomponent < 0) {
    pf_kernel = exptan_smoothing_kernel(91, FLOTSAM_LOBE_HALF_WIDTH_KERNEL_DEG);
  }
  else {
    pf_kernel = lut.phase_function_basis(icomponent);
  }
  Vector ang_kernel = linspace(0,M_PI,pf_kernel.size());
  
  // For the moment, the output resolution will match the input
  int nconv = nang;
  Vector ang_conv = ang; // Link
  Vector pf_conv(nang);

  // Long vector for interpolation
  Vector mu(CosAng.size());

  // 2D array linking to 1D data
  Matrix Mu = mu.reshape(CosAng.size(0), CosAng.size(1));
  Vector ang_great_circle;

  for (int i = 0; i < nconv; ++i) {
    Mu = CosAng * cos(ang_conv(i)) + SinAngCosLon * sin(ang_conv(i));
    mu.where(mu < -1.0) = -1.0;
    mu.where(mu >  1.0) =  1.0;
    ang_great_circle = acos(mu);
    pf_conv(i) = dot_product(weighted_pf, interp(ang_kernel, pf_kernel, ang_great_circle));
  }

  // Normalize to 4pi
  Real dang_conv = ang_conv(1);
  Vector weight_conv = sin(ang_conv);
  weight_conv(0) = dang_conv / 8.0;
  weight_conv(end) = dang_conv / 8.0;
  weight_conv /= sum(weight_conv);
  pf_conv /= dot_product(fabs(pf_conv), weight_conv);
  if (has_zero_integral) {
    pf_conv -= 1.0;
  }
  return pf_conv;
}

#endif // Adept version >= 1.9.10
