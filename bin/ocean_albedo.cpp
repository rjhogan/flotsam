/** @file      ocean_albedo.cpp
    @brief     Offline driver for Cox-Munk ocean albedo model
    @copyright 2023- European Centre for Medium Range Weather Forcasts
    @license   Apache License Version 2 (see the NOTICE.md file for details)

*/
#include <math.h>
#include <adept_arrays.h>
#include <flotsam.h>

#define CHECK(command) istatus = (command); \
  if (istatus != FLOTSAM_SUCCESS) { \
    std::cerr << "Error: " << flotsam_error_message(istatus) << "\n"; \
    return istatus; \
  }

using namespace adept;

int
main()
{
  int istatus;
  // MODIS land channel wavelength bounds in nm
  Vector wav1_nm={620, 841, 459, 545, 1230, 1628, 2105};
  Vector wav2_nm={670, 876, 479, 565, 1250, 1652, 2155};
  // Mid-point wavelengths in metres
  //  Vector wavelengths = 1.0e-9 * 0.5 * (wav1_nm + wav2_nm);

  Vector wavelengths = 1e-9*Vector{350.0, 450.0, 550.0, 650.0, 850.0, 1000.0, 1250.0, 1400.0, 1600.0, 2100.0, 3900.0};
  
  //  Real wind = 10.0;
  // Wind speeds in m/s
  Vector winds = {0.1, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0}; 
  Vector albedo_components(FLOTSAM_NUM_ALBEDO_COMPONENTS);
  int  nmu0 = 101;
  Vector mu0 = linspace(0.0,1.0,nmu0);
  Vector white_sky_albedo(nmu0);
  Vector black_sky_albedo(nmu0);
  Vector backscatter(nmu0);

  std::cout << "% n_wavelength n_wind n_cos_solar_zenith_angle cos_solar_zenith_angles...\n";
  std::cout << wavelengths.size() << " " << winds.size() << " " << nmu0;
  for (int imu0 = 0; imu0 < nmu0; ++imu0) {
    std::cout << " " << mu0(imu0);
  }
  std::cout << "\n";
  std::cout << "% wavelength[nm] wind[m/s] white_sky_albedo black_sky_albedos...\n";
  for (int iwav = 0; iwav < wavelengths.size(); ++iwav) {
    int id = flotsam_new_ocean_brdf(wavelengths(iwav), winds.size(), winds.data());
    for (int iwind = 0; iwind < winds.size(); ++iwind) {
      for (int imu0 = 0; imu0 < nmu0; ++imu0) {
	CHECK(flotsam_get_ocean_albedo_components(id, mu0(imu0), mu0(imu0), 0, winds(iwind),
						  albedo_components.data()));
	//CHECK(flotsam_get_ocean_albedo_components(id, 0.5, mu0(imu0), 3.1415, winds(iwind),
	//					  albedo_components.data()));
	white_sky_albedo(imu0) = M_PI * albedo_components(0);
	black_sky_albedo(imu0) = albedo_components(1);
	backscatter(imu0) = albedo_components(3);
      }
      std::cout << wavelengths(iwav) << " " << winds(iwind) << " " << white_sky_albedo(0);
      for (int imu0 = 0; imu0 < nmu0; ++imu0) {
	std::cout << " " << black_sky_albedo(imu0);
	//std::cout << " " << backscatter(imu0);
      }
      std::cout << "\n";
    }
    CHECK(flotsam_free_ocean_brdf(id));
  }
}
