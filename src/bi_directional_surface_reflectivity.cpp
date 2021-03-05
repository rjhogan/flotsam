/// @file      bi_directional_surface_reflectivity.cpp
/// @brief     angular reflection function for ocean waters
/// @copyright 2017 European Centre for Medium Range Weather Forcasts
/// @license   Apache License Version 2 (see the NOTICE.md file for details)

//#include <math.h>

#include "base.hpp"
#include "bi_directional_surface_reflectivity.hpp"
//#include "ocean_reflection_parameters.h"


using adept::Real;
using adept::Vector;


void
flotsam::white_caps_reflectance(Real wind_speed, Real wavelength, Real& wc_reflectance, Real& wc_cover)
{

  //fractional cover of foam over the water surface (Monahan and Muircheartaigh, 1980)
  wc_cover  = 2.951e-6*pow(wind_speed,3.52);
  if (wc_cover > 1.0){
    wc_cover = 1.0;
  }
  //effective reflectance of the whitecaps (Koepke, 1984)
  Vector wc_reflectance_tab(39);
  Vector wc_wavelengths_tab(39);
  wc_wavelengths_tab = linspace(0.2,4.0,39);
  wc_reflectance_tab << 0.220,0.220,0.220,0.220,0.220,0.220,0.215,0.210,0.200,0.190, 
    0.175,0.155,0.130,0.080,0.100,0.105,0.100,0.080,0.045,0.055,      
    0.065,0.060,0.055,0.040,0.000,0.000,0.000,0.000,0.000,0.000,      
    0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000,0.000;
  
  wc_reflectance = interpolate_1D(wc_wavelengths_tab,wc_reflectance_tab,wavelength*1.0e6); 
  
  //  return wc_reflectance;
}


Vector 
flotsam::ocean_water_refractive_index(Real wavelength, Real salinity)
{
  //Pure water refractive index from Hale and Querry, 
  //Applied Optique, March 1973, Vol. 12,  No. 3, pp. 555-563
  
  //correction of the real part of the refr. index to take into
  //account the change in salinity. 
  //For a typical sea water with salinity 34.3 ppt and chlorinity 19.0 ppt
  //(Sverdrup) a linear correction as function of salinity is applied (McLellan).
  // REFERENCES:
  // Friedman D., Applied Optics, 1969, Vol.8, No.10, pp.2073-2078.
  // McLellan H.J., Elements of physical Oceanography, Pergamon Press, Inc.,
  //        New-York, 1965, p 129.
  // Sverdrup H.V. et al., The Oceans (Prentice-Hall, Inc., Englewood Cliffs,
  //        N.J., 1942, p 173.
  
  Vector wl(62);
  Vector tnr(62);
  Vector tni(62);
  Real ni;
  Real nr;
  
  wl << 0.250,0.275,0.300,0.325,0.345,0.375,0.400,0.425,0.445,0.475,
    0.500,0.525,0.550,0.575,0.600,0.625,0.650,0.675,0.700,0.725,
    0.750,0.775,0.800,0.825,0.850,0.875,0.900,0.925,0.950,0.975,
    1.000,1.200,1.400,1.600,1.800,2.000,2.200,2.400,2.600,2.650,
    2.700,2.750,2.800,2.850,2.900,2.950,3.000,3.050,3.100,3.150,
    3.200,3.250,3.300,3.350,3.400,3.450,3.500,3.600,3.700,3.800,
    3.900,4.000;
  tnr << 1.362,1.354,1.349,1.346,1.343,1.341,1.339,1.338,1.337,1.336,
    1.335,1.334,1.333,1.333,1.332,1.332,1.331,1.331,1.331,1.330,
    1.330,1.330,1.329,1.329,1.329,1.328,1.328,1.328,1.327,1.327,
    1.327,1.324,1.321,1.317,1.312,1.306,1.296,1.279,1.242,1.219,
    1.188,1.157,1.142,1.149,1.201,1.292,1.371,1.426,1.467,1.483,
    1.478,1.467,1.450,1.432,1.420,1.410,1.400,1.385,1.374,1.364,
    1.357,1.351;
  tni << 3.35E-08,2.35E-08,1.60E-08,1.08E-08,6.50E-09,
    3.50E-09,1.86E-09,1.30E-09,1.02E-09,9.35E-10,
    1.00E-09,1.32E-09,1.96E-09,3.60E-09,1.09E-08,
    1.39E-08,1.64E-08,2.23E-08,3.35E-08,9.15E-08,
    1.56E-07,1.48E-07,1.25E-07,1.82E-07,2.93E-07,
    3.91E-07,4.86E-07,1.06E-06,2.93E-06,3.48E-06,
    2.89E-06,9.89E-06,1.38E-04,8.55E-05,1.15E-04,
    1.10E-03,2.89E-04,9.56E-04,3.17E-03,6.70E-03,
    1.90E-02,5.90E-02,1.15E-01,1.85E-01,2.68E-01,
    2.98E-01,2.72E-01,2.40E-01,1.92E-01,1.35E-01,
    9.24E-02,6.10E-02,3.68E-02,2.61E-02,1.95E-02,
    1.32E-02,9.40E-03,5.15E-03,3.60E-03,3.40E-03,
    3.80E-03,4.60E-03;
  
  nr = interpolate_1D(wl,tnr,wavelength*1.0e6);
  ni = interpolate_1D(wl,tni,wavelength*1.0e6);
  
  if (salinity >= 0.0) {
    nr += 0.006*(salinity/34.3);
  }else{
    //default to 34ppt
    nr += 0.006;
  }
  
  Vector n(2);
  n(0) = nr;
  n(1) = ni;

  return n;
  
}


Real
flotsam::ocean_water_leaving_radiance(Real wavelength, Real pigment_conc, 
				      Real dmui, Real dmur, Real wind_spd,
				      std::complex<double> cn2){
  //Adapted from the routine morcasiwat.f90 from 6SV. See 6S documentation
  //    ! Spectral diffuse attenuation coefficient of Case I Waters as Predicted 
  //    ! by MOREL within the spectral range 400-700nm (1988, Journal of Geophysical 
  //    ! Research, Vol.93, No C9, pp 10749-10768)
  //    !
  //    ! input parameters:	wavelength (IN MICROMETERS)
  //    !			pigment concentration
  //    ! output parameter:	reflectance of water
  //    !
  //    ! According Morel,1988, we use:
  //    !
  //    ! Kd	spectral value of the attenuation coefficient for 
  //    !	 downwelling irradiance
  //    !	 with: Kd=Kw+Xc*pigment_conc**e
  //    ! Kw	spectral value of the diffuse attenuation coefficient 
  //    !	 for pure oceanic water
  //    ! Xc, e	spectral coefficients to compute the diffuse attenuation 
  //    !	 coefficient for pigment
  //    ! bb	total backscattering coefficient
  //    !	 with: bb=0.5*bw+bbt*b
  //    ! bw	spectral value of the molecular scattering coefficient of water
  //    ! bbt,b	parameters to compute the scattering coefficients of pigments
  //    !
  //    ! R2	reflectance of water below the surface
  //    !	 with: R2=(0.33/u)*(bb/Kd)	where u is depending of R2

  Real Kw;
  Real Kd;
  Vector tKw(61);
  Vector tXc(61);
  Vector te(61);
  Vector tbw(61);
  Real Xc;
  Real e;
  Real bw;
  Real bb;
  Real b;
  Real bbt;
  Real u1;
  Real R1;
  Real u2;
  Real err;
  int iwl;
  Real R2=0.0;

  tKw << 0.0209,0.0200,0.0196,0.0189,0.0183,
         0.0182,0.0171,0.0170,0.0168,0.0166,
         0.0168,0.0170,0.0173,0.0174,0.0175,
         0.0184,0.0194,0.0203,0.0217,0.0240,
         0.0271,0.0320,0.0384,0.0445,0.0490,
         0.0505,0.0518,0.0543,0.0568,0.0615,
         0.0640,0.0640,0.0717,0.0762,0.0807,
         0.0940,0.1070,0.1280,0.1570,0.2000,
         0.2530,0.2790,0.2960,0.3030,0.3100,
         0.3150,0.3200,0.3250,0.3300,0.3400,
         0.3500,0.3700,0.4050,0.4180,0.4300,
         0.4400,0.4500,0.4700,0.5000,0.5500,
         0.6500;
  tXc << 0.1100,0.1110,0.1125,0.1135,0.1126,
         0.1104,0.1078,0.1065,0.1041,0.0996,
         0.0971,0.0939,0.0896,0.0859,0.0823,
         0.0788,0.0746,0.0726,0.0690,0.0660,
         0.0636,0.0600,0.0578,0.0540,0.0498,
         0.0475,0.0467,0.0450,0.0440,0.0426,
         0.0410,0.0400,0.0390,0.0375,0.0360,
         0.0340,0.0330,0.0328,0.0325,0.0330,
         0.0340,0.0350,0.0360,0.0375,0.0385,
         0.0400,0.0420,0.0430,0.0440,0.0445,
         0.0450,0.0460,0.0475,0.0490,0.0515,
         0.0520,0.0505,0.0440,0.0390,0.0340,
	 0.0300;
  te << 0.668,0.672,0.680,0.687,0.693,
	0.701,0.707,0.708,0.707,0.704,
	0.701,0.699,0.700,0.703,0.703,
	0.703,0.703,0.704,0.702,0.700,
	0.700,0.695,0.690,0.685,0.680,
	0.675,0.670,0.665,0.660,0.655,
	0.650,0.645,0.640,0.630,0.623,
	0.615,0.610,0.614,0.618,0.622,
	0.626,0.630,0.634,0.638,0.642,
	0.647,0.653,0.658,0.663,0.667,
	0.672,0.677,0.682,0.687,0.695,
	0.697,0.693,0.665,0.640,0.620,
	0.600;
  tbw << 0.0076,0.0072,0.0068,0.0064,0.0061,
         0.0058,0.0055,0.0052,0.0049,0.0047,
         0.0045,0.0043,0.0041,0.0039,0.0037,
         0.0036,0.0034,0.0033,0.0031,0.0030,
         0.0029,0.0027,0.0026,0.0025,0.0024,
         0.0023,0.0022,0.0022,0.0021,0.0020,
         0.0019,0.0018,0.0018,0.0017,0.0017,
         0.0016,0.0016,0.0015,0.0015,0.0014,
         0.0014,0.0013,0.0013,0.0012,0.0012,
         0.0011,0.0011,0.0010,0.0010,0.0010,
         0.0010,0.0009,0.0008,0.0008,0.0008,
         0.0007,0.0007,0.0007,0.0007,0.0007,
         0.0007;

  //Contribution from water leaving radiance
  Matrix tdsbnd(6,5); //lut of transmissivity for incident radiation
  Matrix tdvbnd(6,5); //lut of transmissivity for backscattered radiation
  Vector angbnd(5); //zenith angles vaules used to compute tdsbnd/tdvbnd lut
  Vector wsbnd(6);  //surface wind speed vaules used to compute tdsbnd/tdvbnd lut

  angbnd << 0.0,45.0,60.0,75.0,85.0;
  wsbnd << 1.0,3.0,5.0,7.0,9.0,20.0;
  tdsbnd << 0.9787803,0.9706900,0.9479931,0.9690591,0.9980542,
    0.9787738,0.9698871,0.9404608,0.9275920,0.9602273,
    0.9787626,0.9691746,0.9385692,0.9058769,0.9114283,
    0.9787467,0.9685547,0.9381815,0.8951812,0.8713799,
    0.9787264,0.9680276,0.9384519,0.8899654,0.8417820,
    0.9785573,0.9666586,0.9430056,0.8892645,0.7800314;

  tdvbnd << 0.9787764,0.9692680,0.9225163,0.8048478,0.7294627,
    0.9787535,0.9637051,0.9069787,0.8479503,0.8137348,
    0.9787106,0.9564344,0.9044844,0.8678726,0.8453338,
    0.9786453,0.9495727,0.9052351,0.8797889,0.8629867,
    0.9785548,0.9438773,0.9068328,0.8878716,0.8745421,
    0.9775019,0.9288712,0.9153687,0.9091171,0.9036854;

  wavelength *= 1.0e6;

  if (wavelength >= 0.4 && wavelength <= 0.7){ 

    iwl = 1 + static_cast<int>((wavelength-0.4)/0.005);
    Kw  = tKw(iwl);
    Xc  = tXc(iwl);
    e   = te(iwl);
    bw  = tbw(iwl);
  

    if (fabs(pigment_conc) < 0.0001){
      bb = 0.5*bw;
      Kd = Kw;
    }else{
      b   = 0.3*pow(pigment_conc,0.62);
      bbt = 0.002+0.02*(0.5-0.25*log10(pigment_conc))*0.55/wavelength;
      bb  = 0.5*bw+bbt*b;
      Kd  = Kw+Xc*pow(pigment_conc,e);
    } 
      // FIX: these variables may be used uninitialized
      //std::cerr << "Kd=" << Kd << " Kw=" << Kw << " Xc=" << Xc << " pig=" << pigment_conc << " e=" << e << "\n";
      
    u1 = 0.75;
    R1 = 0.33*bb/u1/Kd;
    
    err = 1.0;
    while (err > 0.0001) {
      u2  = 0.9*(1.0-R1)/(1.0+2.25*R1);
      R2  = 0.33*bb/u2/Kd;
      err = fabs((R2-R1)/R2);
      R1 = R2;
    }
    
  }
  //Isotropic back-scattered radiation from within water body
  Real Rwater=R2;

  Real tds; //transmissivity for incident radiation
  Real tdv; //transmissivity for backscattered radiation

  //Factor to take into account the change in solid angle 
  //from under to above the surface.
  //Taking into accunt the total internal reflection,
  //the effective transmissivity for the up-welling radiation
  //can be approximated as tdv/(abs(n2)*abs(n2)) 
  Real one_over_refr_idx2 = (1.0/(abs(cn2)*abs(cn2)));

  //downward water reflectance coefficient for upwelling radiance
  //at the water-air boundary. 
  //Practically this is (1-tdv/(abs(n2)*abs(n2))) but
  //a constant broad-band value can be assumed with negligible
  //loss of accuracy (Austin, 1974)
  Real a = 0.485;

  //wind speed index
  int w_index=0;
  while (w_index < wsbnd.size()-2
	 && wind_spd >= wsbnd(w_index+1)) {
    w_index++;
  }

  //incident angle index
  int i_index=0;
  while (i_index < angbnd.size()-2
	 && acos(dmui)/M_PI*180. >= angbnd(i_index+1)) {
    i_index++;
  }

  //viewing angle index
  int r_index=0;
  while (r_index < angbnd.size()-2
	 && acos(dmur)/M_PI*180. >= angbnd(r_index+1)) {
    r_index++;
  }

  //interpolation...perhaps more accurate than nearest neighbour
  //but id does not compare to DISORT which uses nearest neighbour
  //tds = interpolate_2D(wsbnd,angbnd,tdsbnd,wind_spd,acos(dmui)/M_PI*180.);
  //tdv = interpolate_2D(wsbnd,angbnd,tdvbnd,wind_spd,acos(dmur)/M_PI*180.);

  //std::cout << "w_index: " << w_index << " r_index: " << r_index << " tetar " << acos(dmur)/M_PI*180. << " tdvbnd: " << tdvbnd(w_index,r_index) << "\n";
  //std::cerr << "tdsbnd: " << tdsbnd(w_index,i_index) << " tdvbnd: " << tdvbnd(w_index,r_index) << "\n";
  tds = tdsbnd(w_index,i_index);
  tdv = tdvbnd(w_index,r_index);
  //water-leaving radiance. It depends on incident/scattering angles
  Real Rsw = one_over_refr_idx2 * tds * tdv * Rwater / (1.0-a*Rwater);

  
  return Rsw;


}


Real
flotsam::sea_reflection_coefficient_cox_and_munk(bool apply_shadowing, int slope_dist_shape,
						 std::complex<double> n1, std::complex<double> n2, 
						 Real dmui, Real phii, Real dmur, Real phir,
						 Real wdir, Real wspd){
  //Compute the angular reflectance factor following Cox and Munk. included the contribution from 
  //anisotropy due to wind direction and shadowing

  //Basic equations from Cox and Munk 1954 (and translated in modern notation by e.g. 
  //Feng et al. IEEE TRANSACTIONS ON GEOSCIENCE AND REMOTE SENSING, VOL. 54, NO. 2, FEBRUARY 2016)
  //Wave shadowing effect from Tsang et al. (1985)
  //
  //Inputs: complex refractive indices of air (n1) and water (n2)
  //        incident and scattering zenith and azimuth angles
  //        surface wind speed (m/s) and direction (0 deg is North)
  //
  //Output: dimensionsless surface reflectance factor (defined as pi*I/(E*mu) )
  //        Divide by pi to obtain radiance reflectance in sr-1
  //Author A. Bozzo 2017
  
  if (fabs(dmui-1.0) < 1.0e-9) {
    dmui=0.999999999999;
  }
  if (fabs(dmur-1.0) < 1.0e-9) {
    dmur=0.999999999999;
  }

  
  Real phi_relative = -phii + phir;
  Real phi_wind_relative = -phii + wdir;
  Real thetai = acos(dmui);
  Real thetar = acos(dmur);
  Real sini = sin(thetai);
  Real sinr = sin(thetar);
  //surface slope
  Real Zx = ( -sinr*sin(phi_relative) ) / ( dmui + dmur);
  Real Zy = ( sini + sinr*cos(phi_relative) ) / ( dmui + dmur );
  //adjustment relative to wind direction
  Real Zxp = cos(phi_wind_relative)*Zx + sin(phi_wind_relative)*Zy; //cross wind
  Real Zyp = -sin(phi_wind_relative)*Zx + cos(phi_wind_relative)*Zy; //up wind

  //for anisotropic slope distribution
  Real sigmax2 = 0.003+0.00192*wspd; //cross wind
  Real sigmay2 = 0.00316*wspd; //upwind
  //for isotropic slope distribution
  Real sigma2=0.5*(0.003+0.00512*wspd); //wind direction independent

  Real zeta2 = Zxp*Zxp/sigmax2;
  Real eta2 = Zyp*Zyp/sigmay2; 
  //Real zeta = Zxp/sqrt(sigmax2);
  Real eta = Zyp/sqrt(sigmay2); 
  Real exponent;
  Real slope_dist;
  if (slope_dist_shape == 1 || slope_dist_shape == 2) {
    //anisotropic slope distribution
    exponent = exp( -(zeta2 + eta2)/2. );
    slope_dist = exponent / ( 2.*M_PI*sqrt(sigmax2)*sqrt(sigmay2) );
  }else if (slope_dist_shape == 3){
    //isotropic gaussian distribution
    exponent = exp( -0.5*(Zx*Zx + Zy*Zy)/sigma2 );
    slope_dist = exponent / ( 2.*M_PI*sigma2 );
  }else{
    std::cerr << "Error: type of slope distribution \"" << slope_dist_shape << "\" for sun glint not recognized\n";
    abort();
    //    std::cerr << "set to anisotropic with Gran-Charlier coefficients \n";
    //    exponent = exp( -(zeta2 + eta2)/2. );
    //    slope_dist = exponent / ( 2.*M_PI*sqrt(sigmax2)*sqrt(sigmay2) );
  }

  //Gram-Charlier series for modified peakedness and skewness of
  //the distribution
  Real c21 = 0.01-0.0086*wspd; 
  Real c03 = 0.04-0.0330*wspd; 
  Real c40 = 0.4;
  Real c04 = 0.23;
  Real c22 = 0.12;

  Real t1=-0.5*c21*(zeta2-1.)*eta;
  Real t2=-1./6.*c03*eta*(eta2-3.);
  Real t3=1./24.*c40*(zeta2*zeta2-6.*zeta2+3.);
  Real t4=1./24.*c04*(eta2*eta2-6.*eta2+3.);
  Real t5=1./4.*c22*(zeta2-1)*(eta2-1);

  if (slope_dist_shape == 2){
    //anisotropic slope distribution with 
    //Gram-Charlier coefficients
    slope_dist *= (1.+t1+t2+t3+t4+t5);
  }


  //geometric factors for local tangent plane and local incidence angle
  //2*theta=angle between incident and reflected beam, theta=angle of 
  //incidence relative to the local normal
  Real cos2theta = dmui*dmur + sini*sinr*cos(phi_relative);
  if (cos2theta >= 1.) {
    cos2theta = 0.9999999999;
  }
  if (cos2theta <= -1.) {
    cos2theta = -0.999999999;
  }
  Real mui_i = sqrt(0.5*(1+cos2theta));
  Real sin_i = sqrt(0.5*(1-cos2theta));
  //beta=facet tilt
  Real cosbeta = ( dmui + dmur )/sqrt(2.+2.*cos2theta);
  //alternative(equivalent) computation
  //Real tilt = atan(sqrt(Zx*Zx+Zy*Zy));
  //Real cosbeta = cos(tilt);
  Real cosbeta4 = cosbeta*cosbeta*cosbeta*cosbeta;

  //****************************************************************************
  //Fresnel reflection coefficients (relative to local normal)

  std::complex<double> refr_idx_factor=(n1/n2*sin_i)*(n1/n2*sin_i);
  std::complex<double> Rpar = ( (n1*mui_i - n2*sqrt(1.-refr_idx_factor)) /
  				(n1*mui_i + n2*sqrt(1.-refr_idx_factor)) );
  std::complex<double> Rper = ( (-n2*mui_i + n1*sqrt(1.-refr_idx_factor)) /
  				(n2*mui_i + n1*sqrt(1.-refr_idx_factor)) );
  Real Rf = 0.5*(abs(Rpar)*abs(Rpar) + abs(Rper)*abs(Rper));


  //*************************************************************************


  //sun glint
  Real rhogl = M_PI*slope_dist*Rf/(4.*dmui*dmur*cosbeta4);


  if (apply_shadowing) {
  //Shadowing effect (significant only for large incidence/viewing angles)
  //this is taken from Mishchenko implementation and it does not 
  //depend on wind direction

  //mean square slope
    Real s1=sqrt(2.0*sigma2/M_PI);
    Real s3=1.0/(sqrt(2.0*sigma2));
    Real s2=s3*s3;
    
    Real dcot = dmui/sqrt(1.0-dmui*dmui);
    Real tt1 = exp(-dcot*dcot*s2);
    Real tt2 = 1.0-erf(dcot*s3);
    Real shadow_i = 0.5*(s1*tt1/dcot-tt2); 
    dcot = dmur/sqrt(1.0-dmur*dmur);
    tt1 = exp(-dcot*dcot*s2);
    tt2 = 1.0-erf(dcot*s3);
    Real shadow_r = 0.5*(s1*tt1/dcot-tt2); 
    
    Real shad_i_r=1.0/(1.0+shadow_i+shadow_r);
    rhogl *= shad_i_r; 
  }
   
  return rhogl;

}


Real
flotsam::sea_reflection_matrix(bool apply_shadowing, std::complex<double> cn1, 
			       std::complex<double> cn2, Real wind_spd,
			       Real dmui, Real phii, Real dmur, Real phir){
  // Compute the angular reflectivity from a generic rough surface.
  //#adapted from RMATR subrouine from Mishchenko and Travis code
  //# M. I. Mishchenko and L. D. Travis, Satellite retrieval
  //#   of aerosol properties over the ocean using polarization as well as
  //#   intensity of reflected sunlight.  J. Geophys. Res. 102, 16989-
  //#   17013 (1997).

  if (fabs(dmui-1.0) < 1.0e-9) {
    dmui=0.999999999999;
  }
  if (fabs(dmur-1.0) < 1.0e-9) {
    dmur=0.999999999999;
  }
  //#MEAN SQUARE SURFACE SLOPE (EQ. (18) IN  JGR PAPER)
  Real sigma2=0.5*(0.003+0.00512*wind_spd);
  

  Real dcosi = cos(phii+M_PI);
  Real dsini=sin(phii+M_PI);
  Real dcosr=cos(phir);
  Real dsinr=sin(phir);
  Real dsi=sqrt(1.0-dmui*dmui);
  Real dsr=sqrt(1.0-dmur*dmur);
  Real vi1=dsi*dcosi;
  Real vi2=dsi*dsini;
  Real vi3=-dmui;
  Real vr1=dsr*dcosr;
  Real vr2=dsr*dsinr;
  Real vr3=dmur;


//#!    LOCAL SURFACE NORMAL FOR SPECULAR REFLECTION

  Real unit1=vi1-vr1;
  Real unit2=vi2-vr2;
  Real unit3=vi3-vr3;
  Real fact1=unit1*unit1 + unit2*unit2 + unit3*unit3;
  Real factor=sqrt(1.0/fact1);

//#!    FRESNEL REFLECTION COEFFICIENTS

  Real xi1=factor*(unit1*vi1+unit2*vi2+unit3*vi3);
  std::complex<double> cxi2=1.0 - (1.0-xi1*xi1)*cn1*cn1/(cn2*cn2);
  cxi2=sqrt(cxi2);
  std::complex<double> c1=cn1*xi1;
  std::complex<double> c2=cn2*cxi2;
  std::complex<double> crper=(c1-c2)/(c1+c2);
 //#      crper=1.0
  c1=cn2*xi1;
  c2=cn1*cxi2;
  std::complex<double> crpar=(c1-c2)/(c1+c2);
 //#      crpar=1

 //#!    CALCULATION OF THE AMPLITUDE SCATTERING MATRIX

  Real ti1=-dmui*dcosi;
  Real ti2=-dmui*dsini;
  Real ti3=-dsi;

  Real tr1=dmur*dcosr;
  Real tr2=dmur*dsinr;
  Real tr3=-dsr;

  Real pi1=-dsini;
  Real pi2=dcosi;
  Real pi3=0.0;

  Real pr1=-dsinr;
  Real pr2=dcosr;
  Real pr3=0.0;

  Real pikr=pi1*vr1+pi2*vr2+pi3*vr3;
  Real prki=pr1*vi1+pr2*vi2+pr3*vi3;
  Real tikr=ti1*vr1+ti2*vr2+ti3*vr3;
  Real trki=tr1*vi1+tr2*vi2+tr3*vi3;

  Real e1=pikr*prki;
  Real e2=tikr*trki;
  Real e3=tikr*prki;
  Real e4=pikr*trki;

  std::complex<double> cf11=e1*crper+e2*crpar;
  std::complex<double> cf12=-e3*crper+e4*crpar;
  std::complex<double> cf21=-e4*crper+e3*crpar;
  std::complex<double> cf22=e2*crper+e1*crpar;

 //#!   CALCULATION OF THE STOKES REFLECTION MATRIX WITH
 //#!   SLOPE DISTRIBUTION

  Real vp1=vi2*vr3-vi3*vr2;
  Real vp2=vi3*vr1-vi1*vr3;
  Real vp3=vi1*vr2-vi2*vr1;
  Real dmod=vp1*vp1+vp2*vp2+vp3*vp3;

  dmod=dmod*dmod;
    
  Real rdz2=unit3*unit3;
  Real rdz4=rdz2*rdz2;

 //#!   CALCULATION OF THE SLOPE DISTRIBUTION FUNCTION
  Real dcoeff=1.0/(4.0*dmui*dmur*dmod*rdz4*2.0*sigma2);
  Real dex= -(unit1*unit1 + unit2*unit2)/(2.0*sigma2*rdz2);
  dex=exp(dex);
  dcoeff=dcoeff*fact1*fact1*dex;




  Real af=0.5*dcoeff;
  Real af11=abs(cf11);
  Real af12=abs(cf12);
  Real af21=abs(cf21);
  Real af22=abs(cf22);
  af11=af11*af11;
  af12=af12*af12;
  af21=af21*af21;
  af22=af22*af22;

 //#first of stokes parameter in a 4x4 reflection matrix
  Real r=(af11+af12+af21+af22)*af;

  if (apply_shadowing){
  //Shadowing effect (significant only for large incidence/viewing angles)

  Real s1=sqrt(2.0*sigma2/M_PI);
  Real s3=1.0/(sqrt(2.0*sigma2));
  Real s2=s3*s3;

  Real dcot = dmui/sqrt(1.0-dmui*dmui);
  Real t1 = exp(-dcot*dcot*s2);
  Real t2 = 1.0-erf(dcot*s3);
  Real shadow_i = 0.5*(s1*t1/dcot-t2); 
  dcot = dmur/sqrt(1.0-dmur*dmur);
  t1 = exp(-dcot*dcot*s2);
  t2 = 1.0-erf(dcot*s3);
  Real shadow_r = 0.5*(s1*t1/dcot-t2); 


  Real shad_i_r=1.0/(1.0+shadow_i+shadow_r);
  r *= shad_i_r; 
  }


  return r;

}

Real
flotsam::ocean_brdf(bool apply_shadowing, int slope_dist_shape, Real wind_spd, 
		    Real wind_dir,Real wavelength, Real pigment_conc, 
		    Real salinity, Real mu_sun, Real mu_inst, Real azim){

  //#*********************************      
  //#input parameters      
  //#*********************************           

  Real brdf;
  //#incidence azimuth left at 0.0 because we use the relative azimuth
  Real phii=0.0/180.*M_PI; 

  Vector n=ocean_water_refractive_index(wavelength,salinity);

  std::complex<double> cn1(1.0, 0.0); // #air complex refr index
  std::complex<double> cn2(n(0),n(1));

  //***************************************************
  //Lambertian contribution from white caps

  Real wc_reflectance = 0.0;
  Real wc_cover = 0.0;
  white_caps_reflectance(wind_spd,wavelength,wc_reflectance,wc_cover);

  //***************************************************
  //***************************************************
  //computation of BRDF for input geometry

  Real dmui=mu_sun;
  Real dmur=mu_inst;

  //Contribution from water leaving radiance  
  Real Rsw = ocean_water_leaving_radiance(wavelength, pigment_conc, dmui,dmur,wind_spd,cn2);


  //Mishchenko and Travis reflection matrix. 
  //The first element of the reflection matrix matches Cox and Munk 
  //with isotropic slope distribution. Shadowing already applied
  //Real r=sea_reflection_matrix(apply_shadowing,cn1,cn2,wind_spd,dmui,phii,dmur,azim);
  
  //Cox and Munk anisotropic with shadowing already applied
  Real r=sea_reflection_coefficient_cox_and_munk(apply_shadowing,slope_dist_shape,
						 cn1,cn2,dmui,phii,dmur,azim,wind_dir,wind_spd);

  //Total brdf with contribution from white caps and water-leaving radiance
  //the water leaving radiance is weighted by (1-wc_cover*wc_reflectance). This measn that 
  //we consider the white cap area as not opaque with some radiation crossing into
  //the water.
  //Here brdf is expressed as dimensionsless reflectance factor. The scaling by pi 
  //is to transform into reflectance with dimension sr-1
  brdf=((1.-wc_cover)*r + wc_cover*wc_reflectance + (1.-wc_cover*wc_reflectance)*Rsw)/M_PI;
  
  return brdf;
  
  //  std::cerr << "brdf: " << brdf << " glint: " << r << " wc_cover: " << wc_cover << " wc_reflectance: " << wc_reflectance << " Rsw: " << Rsw << "\n";
}



void
flotsam::ocean_brdf_lut(bool apply_shadowing, int slope_dist_shape, Real wind_spd, 
			Real wind_dir,Real wavelength, Real pigment_conc, Real salinity, 
			Array3D brdf, Vector brdf_zenith, Real& brdf_hemispheric){

  //#*********************************      
  //#input parameters      
  //#*********************************           

  int nazim=brdf.dimension(2); //#number of azimuth angles
  int nzen=brdf.dimension(0); //#number of zenith angles

  Vector zenith_angles(nzen);
  Vector phir(nazim);
  //LUT computed up to 85 deg to avoid problems at large incidence/viewing angles
  //the shadowing correction solves most problems but some
  //approximations for back-scattered radiation from beneath the 
  //surface are not defined for zenith>85
  //The cap at 85deg marginally affects the integral for the 
  //computation of the hemispheirc reflectance
  zenith_angles = linspace(0.0,90.0,nzen);
  phir = linspace(0.0,2.*M_PI,nazim);
  //#incidence azimuth left at 0.0 because we use the relative azimuth
  Real phii=0.0/180.*M_PI; 

  Vector mu_values = cos(zenith_angles/180.*M_PI);

  Vector n=ocean_water_refractive_index(wavelength,salinity);

  std::complex<double> cn1(1.0, 0.0); // #air complex refr index
  std::complex<double> cn2(n(0),n(1));

  //***************************************************
  //Lambertian contribution from white caps

  Real wc_reflectance = 0.0;
  Real wc_cover = 0.0;
  white_caps_reflectance(wind_spd,wavelength,wc_reflectance,wc_cover);

  //***************************************************
  //***************************************************
  //computation of the full BRDF matrix LUT

  for (int i = 0; i < zenith_angles.size(); i++)
    {
      Real dmui=mu_values(i);
      for (int j = 0; j < zenith_angles.size(); j++)
	{
	  Real dmur=mu_values(j);
	  //Contribution from water leaving radiance  
	  Real Rsw = ocean_water_leaving_radiance(wavelength, pigment_conc, 
						  dmui,dmur,
						  wind_spd,cn2);
	  for (int k = 0; k < nazim; k++)
	    {
	      //Mishchenko and Travis reflection matrix. 
	      //The first element of the reflection matrix matches Cox and Munk 
	      //with isotropic slope distribution. Shadowing already applied
	      //Real r=sea_reflection_matrix(apply_shadowing,cn1,cn2,wind_spd,dmui,phii,dmur,phir(k));

	      //Cox and Munk anisotropic with shadowing already applied
	      Real r=sea_reflection_coefficient_cox_and_munk(apply_shadowing,slope_dist_shape,
							     cn1,cn2,dmui,phii,dmur,phir(k),wind_dir,wind_spd);

	      //Total brdf with contribution from white caps and water-leaving radiance
	      //the water leaving radiance is weighted by (1-wc_cover*wc_reflectance). This measn that 
	      //we consider the white cap area as not-opaque with some radiation crossing into
	      //the water.
	      //Here brdf is expressed as dimensionsless reflectance factor. The scaling by pi 
	      //is to transform into reflectance with dimension sr-1
	      brdf(i,j,k)=((1.-wc_cover)*r + wc_cover*wc_reflectance + (1.-wc_cover*wc_reflectance)*Rsw)/M_PI;

	    }
	}
    }



  //Computation of the integral of brdf over phi and mu to obtain the 
  //hemispheric value and the diffuse-to-instrument or
  //sun-to-diffuse terms

  Matrix r_intg1(mu_values.size(),mu_values.size());
  Vector r_intg2(mu_values.size());
  for (int i = 0; i < mu_values.size(); i++){
    for (int j = 0; j < mu_values.size(); j++){
       r_intg1(i,j)=integrate_brdf(phir,brdf(i,j,__));
    }
    r_intg2(i)=integrate_brdf(mu_values,r_intg1(i,__)*mu_values);
  }

  //factor 2 because the integration goes only from 0 to pi/2 but we need to
  //integrate in the full hemisphere (the azimuth part was done already in 
  //the first integration over the viewing hemisphere)
  brdf_hemispheric=2./M_PI*integrate_brdf(mu_values,r_intg2*mu_values);

  brdf_zenith = r_intg2;

}

void
flotsam::extract_from_brdf_lut(Array3D brdf, Vector brdf_zenith, Real brdf_hemispheric,
			       Real mu_sun, Real mu_inst, Real azim, Vector glint_components){


  int nazim=brdf.dimension(2); //#number of azimuth angles in lut
  int nzen=brdf.dimension(0); //#number of zenith angles in lut
  Vector lut_zens(nzen);
  Vector lut_azis(nazim);
  lut_zens = linspace(0.0,90.0,nzen);
  lut_azis = linspace(0.0,2.*M_PI,nazim);
  

  //Real brdf_dir_to_inst=0.0;
  Real brdf_diffuse_to_inst=0.0;
  Real brdf_dir_to_diffuse=0.0;
  //  Vector r_int2(lut_zens.size());

  Vector xi=lut_zens*M_PI/180.;
  Real xsun=acos(mu_sun);
  Real xinst=acos(mu_inst);

  //  for (int i = 0; i < lut_zens.size()-1; i++){
  //    Matrix brdf_lut=brdf(i,__,__); 
  //    r_int2(i)=interpolate_2D(xi,lut_azis,brdf_lut,xinst,azim);
  //  }


  brdf_diffuse_to_inst=interpolate_1D(xi,brdf_zenith,xinst);
  brdf_dir_to_diffuse=interpolate_1D(xi,brdf_zenith,xsun);
  //  brdf_dir_to_inst=interpolate_1D(xi,r_int2,xsun);

  //incident diffuse hemispheric flux into diffuse hemispheric reflected flux
  //(white sky albedo)
  glint_components(0)=brdf_hemispheric;

  //incident direct beam flux into diffuse hemispheric reflected flux
  //(black sky albedo)
  glint_components(1)=brdf_dir_to_diffuse;

  //incident diffuse hemispheric flux into reflected radiance
  glint_components(2)=brdf_diffuse_to_inst/M_PI;

  //incident direct beam flux into reflected radiance
  //  glint_components(3)=brdf_dir_to_inst;

  /*
  //  std::cout << "azi: " << azim << " mu_inst: " << mu_inst << " mu_sun: " << 
  mu_sun << " brdf direct to view: " << glint_components(3) << "\n";
  std::cerr << "brdf diff to view: " << brdf_diffuse_to_inst << "\n";
  std::cerr << "brdf direct to diff: " << brdf_dir_to_diffuse << "\n";
  std::cerr << "brdf hemispheric: " << brdf_hemispheric << "\n";
  std::cerr << "brdf direct to view: " << brdf_dir_to_inst << "\n";
  */


}

Real
flotsam::interpolate_1D(Vector xi,Vector yi,Real xo)
{

  //Linear interpolation of y on x
  
  Real yo=NAN;

  //extrapolation not allowed  
  if (xi(0) < xi(end)) {
    if (xo < xi(0) || xo > xi(end)) {
      std::cerr << "error reading values from BRDF lut";
      std::cerr << "xi(0): " << xi(0) << " xi(end): " << xi(end) << " xo: " << xo << "\n";
    }
  }else{
    if (xo < xi(end) || xo > xi(0)) {
      std::cerr << "error reading values from BRDF lut";
      std::cerr << "xi(0): " << xi(0) << " xi(end): " << xi(end) << " xo: " << xo << "\n";
    }
  }
    
  Real x_end=xi(end);
    
    if (xo==xi(0)){
      yo=yi(0);
      //      std::cerr << "x(0): " << xi(0) << " x0: " << xo << "\n"; 
    }else if (xo==xi(end)){
      //      std::cerr << "x(end): " << xi(end) << " x0: " << xo << "\n"; 
      yo=yi(end);
    }
    else{
      for (int i = 0; i < xi.size()-1; i++){
	if (xi(0) < x_end) {
	  if (xi(i) <= xo && xo < xi(i+1)) {
	    yo=yi(i)+(xo-xi(i))*(yi(i+1)-yi(i))/(xi(i+1)-xi(i));
	    //	    std::cerr << "x(i): " << xi(i) << " x(i+1):  " << xi(i+1) << " x0: " << xo << "\n"; 
	    break;
	  }
	}else if (xi(0) > x_end) {
	  if (xi(i) >= xo && xo > xi(i+1)) {
	    yo=yi(i)+(xo-xi(i))*(yi(i+1)-yi(i))/(xi(i+1)-xi(i));
	    //	    std::cerr << "x(i): " << xi(i) << " x(i+1):  " << xi(i+1) << " x0: " << xo  << "\n"; 
	    break;
	  }
	}else{
	  std::cerr << "error in interpolation of BRDF LUT!! \n";
	}
      }
    }
    
    return yo;
    
    
}

Real
flotsam::interpolate_2D(Vector xi,Vector yi,Matrix zi,Real xo,Real yo)
{

  //Bilinear interpolation of (z) on (x) and (y)
  // input grid from lut must be monotonically increasing

  int x_index = 0;
  int y_index = 0;

  
  
  Real zo;

  //extrapolation not allowed
  //check x range
    if (xi(0) < xi(end)) {
      if (xo < xi(0) || xo > xi(end)) {
	/*	std::cerr << "warning while reading values from BRDF lut \n";
	std::cerr << "xi(0): " << xi(0) << " xi(end): " << xi(end) << " xo: " << xo << "\n";
	std::cerr << "not doing extrapolation, using the extremes of the input grid \n";*/
	if (xo < xi(0)){
	  xo = xi(0);
	}
	if (xo > xi(end)){
	  xo = xi(end);
	}
    }
    }else{
      std::cerr << "error in bi-linear interpolation of brdf: input grid not monotonically increasing \n";
    }
  //check y range
    if (yi(0) < yi(end)) {
      if (yo < yi(0) || yo > yi(end)) {
	/*       	std::cerr << "warning while reading values from BRDF lut \n";
     	std::cerr << "yi(0): " << yi(0) << " yi(end): " << yi(end) << " yo: " << yo << "\n";
	std::cerr << "not doing extrapolation, using the extremes of the input grid \n";*/
	if (yo < yi(0)){
	  yo = yi(0);
	}
	if (yo > yi(end)){
	  yo = yi(end);
	}
    }
    }else{
      std::cerr << "error in bi-linear interpolation of brdf: input grid not monotonically increasing \n";
    }
    
    
    
    if (xo==xi(0) && yo==yi(0)){
      zo=zi(0,0);
    }else if (xo==xi(end) && yo==yi(end)){
      zo=zi(end,end);
    }
    else{
      // Bilinear interpolation 
       while (x_index < xi.size()-2 && xo > xi(x_index+1)) {
	x_index++;
      }
      while (y_index < yi.size()-2 && yo > yi(y_index+1)) {
	y_index++;
      }
    }
    Real inv_denom = 1.0 / (xi(x_index+1)-xi(x_index));
    Real weight0_ = (xi(x_index+1)-xo)*inv_denom;
    Real weight1_ = (xo-xi(x_index))*inv_denom;
    inv_denom = 1.0 / (yi(y_index+1)-yi(y_index));
    Real weight_0 = (yi(y_index+1)-yo)*inv_denom;
    Real weight_1 = (yo-yi(y_index))*inv_denom;
    Real weight00 = weight0_*weight_0;
    Real weight01 = weight0_*weight_1;
    Real weight10 = weight1_*weight_0;
    Real weight11 = weight1_*weight_1;

    zo = weight00*zi(x_index, y_index)
      + weight01*zi(x_index, y_index+1)
      + weight10*zi(x_index+1, y_index)
      + weight11*zi(x_index+1, y_index+1);

    return zo;    
 
    
}
  
  
Real
flotsam::integrate_brdf(Vector xi,Vector yi)
{
  //Simple trapezoidal rule

  Real integral;
  integral=0.0;
  for (int i = 0; i < xi.size()-1; i++){
    integral += (xi(i+1)-xi(i))*(yi(i+1)+yi(i));
  }
  integral = fabs(integral)*0.5;

  //Simpson 1/3 rule 
  //(assuming x monotonically increasing)
  // Real h = xi(1)-xi(0);
  // Real s = yi(0) + yi(end);

  // for (int i = 1; i < xi.size()-1; i=i+2){
  //   s += 4 * yi(i);
  // }
  // for (int i = 2; i < xi.size()-2; i=i+2){
  //   s += 2 * yi(i);
  // }


  // s = s * h / 3.;
  // std::cerr << "trap int: " << integral << " simps int: " << s << "\n";

  return integral;
  // return s;


}
