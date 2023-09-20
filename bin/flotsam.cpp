/** @file      flotsam.cpp
    @brief     Offline driver for FLOTSAM solar radiance model
    @copyright 2017- European Centre for Medium Range Weather Forcasts
    @license   Apache License Version 2 (see the NOTICE.md file for details)

*/

#include <iostream>
#include <fstream>
#include <fenv.h>

#include <flotsam.hpp>

#include "readconfig.h"
#include "Timer.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

using namespace adept;

// Check for errors when reading from config file
#define RC_CHECK_READ(status,varname) if (!(status)) \
    { std::cerr << "Error reading " << varname << "\n"; return 1; }

// Check for errors when calling FLOTSAM functions
int istatus;
#define CHECK(func) { istatus = func; if (istatus < FLOTSAM_SUCCESS) { \
			std::cerr << "Error: " << flotsam_error_message(istatus); return(istatus); } }


// Read a vector of real numbers from the config file, either with the
// name given by varname, or using <varname>_count, <varname>_begin
// and <varname>_end to construct the vector
static bool get_vector(rc_data* config, const char* varname, Vector& vec,
		       const char* units = 0) {
  int n;
  Real* var = rc_get_real_vector(config, varname, &n);
  std::string vname;
  if (var) {
    vec.resize(n);
    for (int i = 0; i < n; ++i) {
      vec[i] = var[i];
    }
    rc_free(var);
  }
  else {
    vname = std::string(varname) + "_count";
    if (!rc_assign_int(config, vname.c_str(), &n)) {
      return false;
    }
    Real val_begin, val_end;
    vname = std::string(varname) + "_begin";
    if (!rc_assign_real(config, vname.c_str(), &val_begin)) {
      return false;
    }
    vname = std::string(varname) + "_end";
    if (!rc_assign_real(config, vname.c_str(), &val_end)) {
      return false;
    }
    vec.clear();
    vec = linspace(val_begin, val_end, n);
  }
  Real scaling;
  vname = std::string(varname) + "_scaling";
  if (rc_assign_real(config, vname.c_str(), &scaling)) {
    vec *= scaling;
  }
  if (units == 0) {
    units = "";
  }
  if (vec.size() <= 8) {
    std::cerr << "  " << varname << " = " << vec << " " << units << "\n";
  }
  else {
    std::cerr << "  " << varname << " = " << vec(range(0,7)) << "... " << units << "\n";
  }

  return true;
}

// Read a phase function from a file into arrays needed by FLOTSAM
static bool read_phase_function(rc_data* config, 
				Vector& pf,
				Vector& pf_smooth,
				Vector& pf_components) {
  int len = 0;
  rc_real* pfc_raw = rc_get_real_vector(config, "phase_function_components", &len);
  char* file_name;
  static const int n_components = flotsam_n_phase_function_components();

  if (pfc_raw) {
    // Read phase function components
    if (len != n_components) {
      std::cerr << "Error: " << len << " phase function components provided, " 
		<< n_components << " needed\n";
      return false;
    }
    else {
      Array<1,rc_real> pfc(pfc_raw, dimensions(n_components)); // Create Vector sharing raw data
      pf_components.clear();
      pf_components = pfc;
      std::cerr << "Reading phase function components:\n"
		<< "  decomposition = " << pf_components << "\n"
		<< "  delta component = " << 1.0-sum(pf_components(flotsam::mass_component_index())) << "\n";
      Matrix pf_smooth_reconstructed = flotsam::reconstruct_phase_function(91, pf_components);
      pf_smooth.clear();
      pf_smooth = pf_smooth_reconstructed(0,__);
    }
    rc_free(pfc_raw);
  }
  else if ((file_name = rc_get_string(config, "phase_function_file"))) {
    // Read phase function from a file

    std::ifstream pf_file(file_name);
    if (pf_file.fail()) {
      rc_free(file_name);
      return false;
    }
    
    std::vector<Real> data;
    
    Real pf_element;
    pf_file >> pf_element;
    while (!pf_file.eof()) {
      data.push_back(pf_element);
      pf_file >> pf_element;
    };
    
    pf.resize(data.size());
    for (int i = 0; i < data.size(); ++i) {
      pf(i) = data[i];
    }
    
    Real integral = flotsam::integrate_phase_function(pf);
    
    pf_smooth.resize(pf.dimensions());
    pf_components.resize(flotsam_n_phase_function_components());
    flotsam::analyse_phase_function(pf, 4.0*M_PI, pf_smooth, pf_components);
    
    std::cerr << "Reading phase function from " << file_name << ":\n"
	      << "  " << pf.size() << " elements\n"
	      << "  integral = " << integral << " (" << integral/M_PI << " pi)\n"
	      << "  decomposition = " << pf_components << "\n"
	      << "  delta component = " << 1.0-sum(pf_components(flotsam::mass_component_index())) << "\n";
    if (rc_get_boolean(config, "use_fitted_phase_function")) {
      Matrix pf_smooth_reconstructed = flotsam::reconstruct_phase_function(91, pf_components);
      pf_smooth.clear();
      pf_smooth = pf_smooth_reconstructed(0,__);
      std::cerr << "  overriding single-scattering phase function with fit\n";
    }

    if (abs(integral - 4.0*M_PI) > 0.05) {
      std::cerr << "  WARNING: integral assumed to be 4 pi\n";
    }
    
    rc_free(file_name);
  }
  else {
    std::cerr << "Error: neither phase_function_file nor phase_function_components were specified\n";
    return false;
  }

  return true;
}
  
// Print various smoothed phase functions to standard output
static int analyse_pf(rc_data* config) {
  // Read phase function
  Vector pf_orig, pf_smooth, pf_components;
  RC_CHECK_READ(read_phase_function(config, pf_orig, pf_smooth, pf_components), "phase function");
  int nang = pf_orig.size();
  Matrix pfs(14, nang);

  pfs = 0.0;

  pfs(0,__) = linspace(0,180.0,nang);
  pfs(1,__) = pf_orig;
  pfs(2,__) = pf_smooth;
  Matrix reconstructed_pf = flotsam::reconstruct_phase_function(nang, pf_components);
  pfs(range(5,7),__) = reconstructed_pf(range(0,2),__);

  std::cerr << "First phase function convolution...\n";
  Vector pf_conv1 = flotsam::convolve_phase_function(pf_smooth);
  std::cerr << "  integral = " << flotsam::integrate_phase_function(pf_conv1) << "\n";
  pfs(3,__) = pf_conv1;

   std::cerr << "Second phase function convolution...\n";
  Vector pf_conv2 = flotsam::convolve_phase_function(pf_conv1);
  std::cerr << "  integral = " << flotsam::integrate_phase_function(pf_conv2) << "\n";
  pfs(4,__) = pf_conv2;

  pfs(8,__) = flotsam::convolve_phase_function(pfs(7,__));
  pfs(9,__) = flotsam::convolve_phase_function(pfs(8,__));
  pfs(10,__) = flotsam::convolve_phase_function(pfs(9,__));
  pfs(11,__) = flotsam::convolve_phase_function(pfs(10,__));
  pfs(12,__) = flotsam::convolve_phase_function(pfs(11,__));
  pfs(13,__) = flotsam::convolve_phase_function(pfs(12,__));

  adept::set_array_print_curly_brackets(false);
  std::cerr << "Output columns:\n"
	    << "  1. Scattering angle (degrees)\n"
	    << "  2. Original phase function\n"
	    << "  3. Phase function after narrow forward scattering removed and smoothing due to narrow forward scattering\n"
	    << "  4. Smoothing of (3) due to wide-angle forward scattering\n"
	    << "  5. Smoothing of (3) due to two instances of wide-angle forward scattering\n"
	    << "  6. Reconstruction of (3)\n"
	    << "  7. Smoothing of (6) due to wide-angle forward scattering\n"
	    << "  8. Smoothing of (6) due to two instances of wide-angle forward scattering\n";
  std::cout << pfs.T();
  return 0;
}


// Coarsen a profile of optical-depth like information in "var"
// according to the contents of "coarsening".  "var" is returned with
// a new size.  "coarsening" contains a list of numbers of layers to
// sum starting from the top of the profile.  If the total number of
// layers accounted for in "coarsening" is less than the number in
// "var" then the lowest layers are retained without being summed.
int
coarsen_profile(const intVector& coarsening, Vector& var)
{
  if (var.empty()) {
    // No action on empty vectors
    return true;
  }

  // Number of layers to remove
  int n_remove = sum(coarsening - 1);
  int n = var.size() - n_remove;
  if (n < 1) {
    return false;
  }
  Vector newvar(n);
  int i = 0;
  for (int j = 0; j < coarsening.size(); ++j) {
    newvar(j) = sum(var(range(i,i+coarsening(j)-1)));
    i += coarsening(j);
  }
  if (n > coarsening.size()) {
    newvar(range(coarsening.size(),end)) = var(range(i,end));
  }
  var.clear();
  var >>= newvar; // Steal data
  return true;
}

int
coarsen_edge_pressure(const intVector& coarsening, Vector& edge_pressure)
{
  // Number of layers to remove
  int n_remove = sum(coarsening - 1);
  int n = edge_pressure.size() - n_remove - 1;
  if (n < 1) {
    return false;
  }
  Vector new_ep(n+1);
  new_ep(0) = edge_pressure(0);
  int i = 0;
  for (int j = 0; j < coarsening.size(); ++j) {
    i += coarsening(j);
    new_ep(j+1) = edge_pressure(i);
  }
  if (n > coarsening.size()) {
    new_ep(range(coarsening.size()+1,end)) = edge_pressure(range(i+1,end));
  }
  edge_pressure.clear();
  edge_pressure >>= new_ep; // Steal data
  return true;
}



int
main(int argc, const char** argv)
{
  const Real deg_to_rad = M_PI / 180.0;

  std::cerr << "##########################\n";
#ifdef PACKAGE_STRING
  std::cerr << "  " PACKAGE_STRING "\n";
#else
  std::cerr << "  FLOTSAM\n";
#endif
  std::cerr << "##########################\n";

  // Stop on floating-point exception
  feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);

  if (argc < 2) {
    // No arguments provided
    std::cerr << "Usage: " << argv[0] << " [param1=value1] [param2=value2...] file1.cfg [file2.cfg...]\n";
    return 1;
  }

  rc_data* config;


  // Find the first config file on the command line
  int ifile = rc_get_file(argc, argv);
  if (ifile > 0) {
    config = rc_read(argv[ifile], stderr);
    if (!config) {
      std::cerr << "Error reading configuration information from " << argv[ifile] << "\n";
      return 2;
    }
  }
  else {
    config = rc_read(NULL, stderr);
  }
  
  // Read any subsequent files
  for (int i = ifile+1; i < argc; ++i) {
    rc_append(config, argv[i], stderr);
    /*
    if (!rc_append(config, argv[i], stderr)) {
      std::cerr << "Error reading configuration information from " << argv[i] << "\n";
      return 2;
    }
    */
  }


  // Supplement configuration information with command-line arguments
  rc_register_args(config, argc, argv);

#ifdef DEBUG
  std::cerr << rc_sprint(config);
#endif

  // Configure calculation mode
  enum {
    MODE_DIRECT, MODE_COMPONENTS, MODE_INTERNALS,
    MODE_JACOBIAN, MODE_JACOBIAN_THREAD_UNSAFE,
    MODE_PHASE_FUNCTION, MODE_BASIS_FUNCTIONS
  } calculation_mode = MODE_DIRECT;

  const char* mode_str = rc_get_string(config, "mode");
  if (mode_str) {
    if (mode_str == std::string("direct")) {
      calculation_mode = MODE_DIRECT;
    }
    else if (mode_str == std::string("components")) {
      calculation_mode = MODE_COMPONENTS;
    }
    else if (mode_str == std::string("internals")) {
      calculation_mode = MODE_INTERNALS;
    }
    else if (mode_str == std::string("jacobian")) {
      calculation_mode = MODE_JACOBIAN;
    }
    else if (mode_str == std::string("jacobian_thread_unsafe")) {
      calculation_mode = MODE_JACOBIAN_THREAD_UNSAFE;
    }
    else if (mode_str == std::string("phase_function")) {
      calculation_mode = MODE_PHASE_FUNCTION;
      // We analyse the phase function and print out the results,
      // rather than running FLOTSAM algorithm
      return analyse_pf(config);
    }
    else if (mode_str == std::string("basis_functions")) {
      calculation_mode = MODE_BASIS_FUNCTIONS;
#if ADEPT_VERSION >= 20005
      adept::set_array_print_style(PRINT_STYLE_PLAIN);
#endif
      std::cout << flotsam::basis_functions() << "\n";
      return 0;
    }
    else {
      std::cerr << "Calculation mode \"" << mode_str
		<< "\" not recognised (use \"direct\", \"components\", \"internals\", \"jacobian\", \"jacobian_thread_unsafe\", \"phase_function\", \"basis_functions\")\n";
      return 1;
    }
  }
  else {
    mode_str = "direct";
  }

  std::cerr << "Configuration:\n";

  std::cerr << "  mode = " << mode_str << "\n";

  // Read configuration data from config file
  Real wavelength;
  RC_CHECK_READ(rc_assign_real(config, "wavelength", &wavelength), "wavelength");

  std::cerr << "  wavelength = " << wavelength << " m\n";

  // Configure surface reflection
  Real albedo = 0.0;
  Real wind_speed = 2.0; // m/s
  enum { SURFACE_LAMBERTIAN, SURFACE_COX_MUNK } surface_type = SURFACE_LAMBERTIAN;
  char* surface_str = rc_get_string(config, "surface");
  if (surface_str) {
    if (surface_str == std::string("lambertian")) {
      surface_type = SURFACE_LAMBERTIAN;
    }
    else if (surface_str == std::string("cox_munk")) {
      surface_type = SURFACE_COX_MUNK;
    }
    else {
      std::cerr << "Surface type \"" << surface_str << "\" not recognised (use \"lambertian\" or \"cox_munk\")\n";
      return 1;
    }
  }

  switch (surface_type) {
  case SURFACE_LAMBERTIAN:
    std::cerr << "  surface = lambertian\n";
    RC_CHECK_READ(rc_assign_real(config, "albedo", &albedo), "albedo");
    std::cerr << "  albedo = " << albedo << "\n";
    break;
  case SURFACE_COX_MUNK:
    std::cerr << "  surface = cox_munk, default wind_speed = " << wind_speed << " m/s\n";
  }

  Vector sza, sensor_zenith_angle, azim, optical_depth;
  RC_CHECK_READ(get_vector(config, "sza", sza, "deg"), "sza");
  RC_CHECK_READ(get_vector(config, "sensor_zenith_angle",
			   sensor_zenith_angle, "deg"), "sensor_zenith_angle");
  RC_CHECK_READ(get_vector(config, "azim", azim, "deg"), "azim");
  RC_CHECK_READ(get_vector(config, "optical_depth", optical_depth),
		"optical_depth");

  Real single_scattering_albedo;
  RC_CHECK_READ(rc_assign_real(config, "single_scattering_albedo",
			       &single_scattering_albedo),
		"single_scattering_albedo");
  std::cerr << "  single_scattering_albedo = " << single_scattering_albedo << "\n";


  // Construct pressure profile
  int layer_count;
  Vector edge_pressure;
  Vector normalized_od;
  Vector molecular_absorption_od, molecular_scattering_od;

  if (get_vector(config, "normalized_od", normalized_od)) {
    // We have a normalized optical depth profile; now read in the
    // pressure at the base of each layer
    layer_count = normalized_od.size();
    Vector pressure_base;
    RC_CHECK_READ(get_vector(config, "pressure_base", pressure_base, "Pa"), "pressure_base");
    if (pressure_base.size() != layer_count) {
      std::cerr << "Error: normalized_od and pressure_base must have the same length\n";
      return 1;
    }
    edge_pressure.resize(layer_count+1);
    edge_pressure(0) = 0.0;
    edge_pressure(range(1,end)) = pressure_base;
    get_vector(config, "molecular_absorption_od", molecular_absorption_od);
    get_vector(config, "molecular_scattering_od", molecular_scattering_od);

    // Optionally coarsen layers
    Vector coarsening_double;
    if (get_vector(config, "coarsening", coarsening_double)) {
      intVector coarsening = coarsening_double;
      coarsen_profile(coarsening, normalized_od);
      coarsen_profile(coarsening, molecular_absorption_od);
      coarsen_profile(coarsening, molecular_scattering_od);
      coarsen_edge_pressure(coarsening, edge_pressure);
      layer_count = normalized_od.size();
    }

    std::cerr << "Sum of normalized optical depth profile is " << sum(normalized_od) << "\n";
  }
  else {
    // Read in the number of layers and the surface pressure, and
    // construct a normalized optical depth profile in which optical
    // depth is equal in each layer
    RC_CHECK_READ(rc_assign_int(config, "layer_count", &layer_count),
		  "layer_count");
    Real surface_pressure;
    RC_CHECK_READ(rc_assign_real(config, "surface_pressure",
				 &surface_pressure), "surface_pressure");
    edge_pressure = linspace(0, surface_pressure, layer_count+1);
    normalized_od.resize(layer_count);
    normalized_od = 1.0 / layer_count;
  }

  int no_rayleigh = rc_get_boolean(config, "no_rayleigh");
  if (no_rayleigh) {
    std::cerr << "Rayleigh scattering is OFF\n";
  }
  else {
    std::cerr << "Rayleigh scattering is ON\n";
  }

  int nrepeat = 1;
  rc_assign_int(config, "repeat", &nrepeat);

  if (nrepeat > 1) {
    std::cerr << "Repeating calculations " << nrepeat << " times for benchmarking\n";
  }


  // Particulates in all layers
  intVector loc(layer_count);
  loc = range(0,layer_count-1);

  // Timing
  Timer timer;
  int t_setup = timer.new_activity("channel and background-profile setup");
  //  int t_phase_function = timer.new_activity("initial phase-function processing");
  int t_geometry = timer.new_activity("geometry setup");
  int t_brdf_setup = timer.new_activity("BRDF setup");
  int t_brdf_calc = timer.new_activity("BRDF calculation");
  int t_reflectance = timer.new_activity("reflectance calculation");

  const int n_components = flotsam_n_phase_function_components();

  timer.start(t_setup);
  // Create channel (wavelength-dependent information)
  int ichan;
  if (no_rayleigh) {
    ichan = flotsam_new_channel_vacuum();
  }
  else {
    ichan = flotsam_new_channel_rayleigh_only(wavelength);
  }
  CHECK(ichan);

  // Create background profile (meteorological information)
  int iprof = flotsam_new_background_profile();
  CHECK(iprof);

  CHECK(flotsam_set_edge_pressure(iprof, layer_count, edge_pressure.data()));

  // Create band profile (channel and meteorology)
  int iband = flotsam_new_band_profile();
  CHECK(iband);
  CHECK(flotsam_init_band_profile(iband, ichan, iprof));

  if (!molecular_absorption_od.empty()) {
    CHECK(flotsam_set_gas_absorption_optical_depth(iband,
						   molecular_absorption_od.size(),
						   molecular_absorption_od.data()));
  }
  if (!molecular_scattering_od.empty()) {
    CHECK(flotsam_set_rayleigh_optical_depth(iband,
					     molecular_scattering_od.size(),
					     molecular_scattering_od.data()));
  }
  timer.stop();

  // Report atmospheric background optical depths
  Real gas_od;
  CHECK(flotsam_get_total_rayleigh_optical_depth(iband, &gas_od));
  std::cerr << "  total_rayleigh_optical_depth = " << gas_od << "\n";
  CHECK(flotsam_get_total_gas_absorption_optical_depth(iband, &gas_od));
  std::cerr << "  total_gas_absorption_optical_depth = " << gas_od << "\n";


  // Read phase function
  Vector pf_orig, pf_smooth, pf_components;
  RC_CHECK_READ(read_phase_function(config, pf_orig, pf_smooth, pf_components), "phase function");

  int ibrdf = -1;
  Vector albedo_vector;
  if (surface_type == SURFACE_LAMBERTIAN) {
    albedo_vector.resize(1);
    albedo_vector = albedo;
  }
  else { // Cox-Munk
    Real pigment_conc_mg_m3 = 0.02;
    Real salinity_ppt = 35.0;
    int no_wave_shadowing = 0;
    //shape of the glint distribution:
    //1=anisotropic (wind dependent)
    //2=anisotropic with Gram-Charlier coefficients
    //3=isotropic (wind independent)
    int slope_dist_shape = 2;

    albedo_vector.resize(FLOTSAM_NUM_ALBEDO_COMPONENTS);
    std::cerr << "Computing Cox-Munk ocean BRDF\n";

    RC_CHECK_READ(rc_assign_real(config, "wind_speed", &wind_speed),
		  "wind_speed (needed for Cox-Munk model)");
    Vector winds_lut(1);
    winds_lut << wind_speed;

    rc_assign_real(config, "pigment_conc_mg_m3", &pigment_conc_mg_m3);
    rc_assign_real(config, "salinity_ppt", &salinity_ppt);
    rc_assign_int(config, "glint_shape", &slope_dist_shape);
    no_wave_shadowing = rc_get_boolean(config, "no_wave_shadowing");

    std::cerr << "  wind_speed = " << wind_speed << " m/s\n";
    std::cerr << "  pigment_conc_mg_m3 = " << pigment_conc_mg_m3 << "\n";
    std::cerr << "  salinity_ppt = " << salinity_ppt << "\n";
    std::cerr << "  no_wave_shadowing = " << (no_wave_shadowing != 0) << "\n";
    std::cerr << "  glint_shape = " << slope_dist_shape << " ";
    if (slope_dist_shape == 1) {
      std::cerr << "(anisotropic, wind dependent)\n";
    }
    else if (slope_dist_shape == 2) {
      std::cerr << "(anisotropic with Gram-Charlier coefficients)\n";
    }
    else if (slope_dist_shape == 3) {
      std::cerr << "(isotropic, wind independent)\n";
    }
    else {
      std::cerr << "(unknown)\n";
    }

    timer.start(t_brdf_setup);
    ibrdf = flotsam_new_ocean_brdf_detailed(wavelength, winds_lut.size(), winds_lut.data(),
					    pigment_conc_mg_m3, salinity_ppt, !no_wave_shadowing,
					    slope_dist_shape);
    timer.stop();
    CHECK(ibrdf);
  }

  Vector od(layer_count);
  Vector ssa(layer_count);
  ssa = single_scattering_albedo;
  Vector pf(layer_count);
  Matrix pfc(layer_count,n_components);

  if (calculation_mode != MODE_INTERNALS) {
    // Print first line of output file describing number of
    // calculations done in each dimension
    std::cout << sza.size() << " " << sensor_zenith_angle.size() << " "
	      << azim.size() << " " << optical_depth.size();
    if (calculation_mode == MODE_COMPONENTS) {
      std::cout << " " << FLOTSAM_NUM_COMPONENTS << "\n";
    }
    else {
      std::cout << " 1\n";
    }
  }

  int ncount = 0;
  int ncount_geometry = 0;

  // Arrays to store the Jacobian
  Vector d_albedo(albedo_vector.size());
  Vector d_od(layer_count);
  Vector d_ssa(layer_count);
  Vector d_pf(layer_count);
  Matrix d_pfc(layer_count,n_components);

  d_pfc = 0.0;

  // Loop over solar zenith angle
  for (int isza = 0; isza < sza.size(); ++isza) {
    Real mu_sun = cos(sza(isza) * deg_to_rad);

    // Loop over sensor zenith angle
    for (int isens = 0; isens < sensor_zenith_angle.size(); ++isens) {
      Real mu_sens = cos(sensor_zenith_angle(isens) * deg_to_rad);

      // Loop over azimuthal separation of sun and sensor
      for (int iazim = 0; iazim < azim.size(); ++iazim) {
	Real az = azim(iazim) * deg_to_rad;

	timer.start(t_geometry);
	Real scat_ang = flotsam_scattering_angle(mu_sun, mu_sens, az);
	CHECK(flotsam_set_geometry(iband, mu_sun, mu_sens, az));
	pf = flotsam_interp_phase_func(pf_smooth.size(), pf_smooth.data(),
				       scat_ang);
	timer.stop();
	++ncount_geometry;

	if (surface_type == SURFACE_COX_MUNK) {
	  timer.start(t_brdf_calc);
	  CHECK(flotsam_get_ocean_albedo_components(ibrdf, 
		    mu_sun, mu_sens, az, wind_speed, albedo_vector.data()));
	  timer.stop();
	}

	for (int i = 0; i < n_components; ++i) {
	  pfc(__,i) = pf_components(i);
	}

	// Loop over particulate optical depth
	for (int iod = 0; iod < optical_depth.size(); ++iod) {

	  od = optical_depth(iod) * normalized_od;

	  Real reflectance[FLOTSAM_NUM_COMPONENTS];

	  if (calculation_mode != MODE_INTERNALS) {
	    // Print coordinate information
	    std::cout << sza(isza) << " "
		      << sensor_zenith_angle(isens) << " "
		      << azim(iazim) << " "
		      << optical_depth(iod) << " ";
	  }

	  if (calculation_mode == MODE_DIRECT) {
	    timer.start(t_reflectance);
	    for (int irepeat = 0; irepeat < nrepeat; ++irepeat) {
	      CHECK(flotsam_reflectance(iband, albedo_vector.size(),
					albedo_vector.data(),
					layer_count, loc.data(),
					od.data(), ssa.data(),
					pf.data(), pfc.data(),
					reflectance));
	      ++ncount;
	    }
	    timer.stop();
	    std::cout << reflectance[0] << "\n"; // Total reflectance
	  }
	  else if (calculation_mode == MODE_COMPONENTS) {
	    timer.start(t_reflectance);
	    for (int irepeat = 0; irepeat < nrepeat; ++irepeat) {
	      CHECK(flotsam_reflectance_components(iband, albedo_vector.size(),
						   albedo_vector.data(),
						   layer_count, loc.data(),
						   od.data(), ssa.data(),
						   pf.data(), pfc.data(),
						   reflectance));
	      ++ncount;
	    }
	    timer.stop();
	    std::cout << reflectance[0] << " " // Total reflectance
		      << reflectance[1] << " " // Direct component
		      << reflectance[2] << " " // Lobe dn component
		      << reflectance[3] << " " // Lobe up component
		      << reflectance[4] << " " // Diffuse component
		      << reflectance[5] << " " // Surface component
		      << reflectance[6] << "\n"; // TOA flux
	  }
	  else if (calculation_mode == MODE_INTERNALS) {
	    timer.start(t_reflectance);
	    for (int irepeat = 0; irepeat < nrepeat; ++irepeat) {
	      // Note that "-" means print internals to standard
	      // output
	      CHECK(flotsam_reflectance_write_internals(iband,
					albedo_vector.size(),
					albedo_vector.data(),
					layer_count, loc.data(),
					od.data(), ssa.data(),
					pf.data(), pfc.data(),
					reflectance, "-"));
	      ++ncount;
	    }
	    timer.stop();
	  }
	  else if (calculation_mode == MODE_JACOBIAN) {
	    timer.start(t_reflectance);
	    for (int irepeat = 0; irepeat < nrepeat; ++irepeat) {
	      CHECK(flotsam_reflectance_jacobian(iband, albedo_vector.size(),
						 albedo_vector.data(),
						 layer_count, loc.data(),
						 od.data(), ssa.data(),
						 pf.data(), pfc.data(),
						 reflectance, d_albedo.data(),
						 d_od.data(), d_ssa.data(),
						 d_pf.data(), d_pfc.data()));
	      ++ncount;
	    }
	    timer.stop();
	    std::cout << reflectance[0] << "\n";
	  }
	  else if (calculation_mode == MODE_JACOBIAN_THREAD_UNSAFE) {
	    timer.start(t_reflectance);
	    for (int irepeat = 0; irepeat < nrepeat; ++irepeat) {
	      CHECK(flotsam_reflectance_jacobian_thread_unsafe(iband, albedo_vector.size(),
						 albedo_vector.data(),
						 layer_count, loc.data(),
						 od.data(), ssa.data(),
						 pf.data(), pfc.data(),
						 reflectance, d_albedo.data(),
						 d_od.data(), d_ssa.data(),
						 d_pf.data(), d_pfc.data()));
	      ++ncount;
	    }
	    timer.stop();
	    std::cout << reflectance[0] << "\n";
	  }
	}
      }
    }
  }

  std::cerr << "Summary:\n";
  std::cerr << "  number of layers                      = " << layer_count << "\n";
  std::cerr << "  number of background/BRDF setup calls = 1\n";
  std::cerr << "  number of geometry setup calls        = " << ncount_geometry << "\n";
  if (surface_type == SURFACE_COX_MUNK) {
    std::cerr << "  number of BRDF calculations           = " << ncount_geometry << "\n";
  }
  std::cerr << "  number of reflectance calls           = " << ncount << "\n";
  std::cerr << "Timing:\n";
  std::cerr << "  time per background setup call        = "
	    << timer.timing(t_setup)*1000000.0 << " us\n";
  if (surface_type == SURFACE_COX_MUNK) {
    std::cerr << "  time per BRDF setup call              = "
	      << timer.timing(t_brdf_setup)*1000000.0 << " us\n";
  }
  std::cerr << "  time per geometry setup call          = "
	    << timer.timing(t_geometry)*1000000.0/ncount_geometry << " us\n";
  if (surface_type == SURFACE_COX_MUNK) {
    std::cerr << "  time per BRDF calculation             = "
	      << timer.timing(t_brdf_calc)*1000000.0/ncount_geometry << " us\n";
  }

  std::cerr << "  time per reflectance calculation      = "
	    << timer.timing(t_reflectance)*1000000.0/ncount << " us\n";

  std::cerr << "  total time = " << timer.timing(t_setup)+timer.timing(t_brdf_setup)+timer.timing(t_geometry)+timer.timing(t_brdf_calc)+timer.timing(t_reflectance) << "\n";
  
  timer.print_on_exit(false);
  return 0;
}
