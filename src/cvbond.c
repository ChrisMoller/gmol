#include <stdio.h>
#include <math.h>

#define AN_H	1
#define AN_C	6
#define AN_N	7
#define AN_O	8
#define AN_F	9
#define AN_Si	14
#define AN_P	15
#define AN_S	16
#define AN_Cl	17
#define AN_Br	35
#define AN_I	53

#define combo(l,r)	(((l) << 8) | (r))

// from http://chem.libretexts.org/Textbook_Maps/General_Chemistry_Textbook_Maps/Map%3A_Chemistry%3A_The_Central_Science_(Brown_et_al.)/08._Basic_Concepts_of_Chemical_Bonding/8.8%3A_Strength_of_Covalent_Bonds

/***
    max =  432
    min =  149
 ***/

// from http://www.fourmilab.ch/documents/specrend/specrend.c

/* A colour system is defined by the CIE x and y coordinates of
   its three primary illuminants and the x and y coordinates of
   the white point. */

struct colourSystem {
    char *name;     	    	    /* Colour system name */
    double xRed, yRed,	    	    /* Red x, y */
           xGreen, yGreen,  	    /* Green x, y */
           xBlue, yBlue,    	    /* Blue x, y */
           xWhite, yWhite,  	    /* White point x, y */
	   gamma;   	    	    /* Gamma correction for system */
};

/* White point chromaticities. */

#define IlluminantC     0.3101, 0.3162	    	/* For NTSC television */
#define IlluminantD65   0.3127, 0.3291	    	/* For EBU and SMPTE */
#define IlluminantE 	0.33333333, 0.33333333  /* CIE equal-energy illuminant */

/*  Gamma of nonlinear correction.

    See Charles Poynton's ColorFAQ Item 45 and GammaFAQ Item 6 at:
    
       http://www.poynton.com/ColorFAQ.html
       http://www.poynton.com/GammaFAQ.html
 
*/

#define GAMMA_REC709	0		/* Rec. 709 */

static struct colourSystem
                  /* Name                  xRed    yRed    xGreen  yGreen  xBlue  yBlue    White point        Gamma   */
SMPTEsystem =  { "SMPTE",              0.630,  0.340,  0.310,  0.595,  0.155,  0.070,  IlluminantD65,  GAMMA_REC709 };


static void
xyz_to_rgb(struct colourSystem *cs,
	   double xc, double yc, double zc,
	   double *r, double *g, double *b)
{
    double xr, yr, zr, xg, yg, zg, xb, yb, zb;
    double xw, yw, zw;
    double rx, ry, rz, gx, gy, gz, bx, by, bz;
    double rw, gw, bw;

    xr = cs->xRed;    yr = cs->yRed;    zr = 1 - (xr + yr);
    xg = cs->xGreen;  yg = cs->yGreen;  zg = 1 - (xg + yg);
    xb = cs->xBlue;   yb = cs->yBlue;   zb = 1 - (xb + yb);

    xw = cs->xWhite;  yw = cs->yWhite;  zw = 1 - (xw + yw);

    /* xyz -> rgb matrix, before scaling to white. */
    
    rx = (yg * zb) - (yb * zg);  ry = (xb * zg) - (xg * zb);  rz = (xg * yb) - (xb * yg);
    gx = (yb * zr) - (yr * zb);  gy = (xr * zb) - (xb * zr);  gz = (xb * yr) - (xr * yb);
    bx = (yr * zg) - (yg * zr);  by = (xg * zr) - (xr * zg);  bz = (xr * yg) - (xg * yr);

    /* White scaling factors.
       Dividing by yw scales the white luminance to unity, as conventional. */
       
    rw = ((rx * xw) + (ry * yw) + (rz * zw)) / yw;
    gw = ((gx * xw) + (gy * yw) + (gz * zw)) / yw;
    bw = ((bx * xw) + (by * yw) + (bz * zw)) / yw;

    /* xyz -> rgb matrix, correctly scaled to white. */
    
    rx = rx / rw;  ry = ry / rw;  rz = rz / rw;
    gx = gx / gw;  gy = gy / gw;  gz = gz / gw;
    bx = bx / bw;  by = by / bw;  bz = bz / bw;

    /* rgb of the desired point */
    
    *r = (rx * xc) + (ry * yc) + (rz * zc);
    *g = (gx * xc) + (gy * yc) + (gz * zc);
    *b = (bx * xc) + (by * yc) + (bz * zc);
}


/*                          SPECTRUM_TO_XYZ

    Calculate the CIE X, Y, and Z coordinates corresponding to
    a light source with spectral distribution given by  the
    function SPEC_INTENS, which is called with a series of
    wavelengths between 380 and 780 nm (the argument is 
    expressed in meters), which returns emittance at  that
    wavelength in arbitrary units.  The chromaticity
    coordinates of the spectrum are returned in the x, y, and z
    arguments which respect the identity:

            x + y + z = 1.
*/

static void
spectrum_to_xyz(double (*spec_intens)(double wavelength),
		double *x, double *y, double *z)
{
    int i;
    double lambda, X = 0, Y = 0, Z = 0, XYZ;

    /* CIE colour matching functions xBar, yBar, and zBar for
       wavelengths from 380 through 780 nanometers, every 5
       nanometers.  For a wavelength lambda in this range:

            cie_colour_match[(lambda - 380) / 5][0] = xBar
            cie_colour_match[(lambda - 380) / 5][1] = yBar
            cie_colour_match[(lambda - 380) / 5][2] = zBar

	To save memory, this table can be declared as floats
	rather than doubles; (IEEE) float has enough 
	significant bits to represent the values. It's declared
	as a double here to avoid warnings about "conversion
	between floating-point types" from certain persnickety
	compilers. */

    static double cie_colour_match[81][3] = {
        {0.0014,0.0000,0.0065}, {0.0022,0.0001,0.0105}, {0.0042,0.0001,0.0201},
        {0.0076,0.0002,0.0362}, {0.0143,0.0004,0.0679}, {0.0232,0.0006,0.1102},
        {0.0435,0.0012,0.2074}, {0.0776,0.0022,0.3713}, {0.1344,0.0040,0.6456},
        {0.2148,0.0073,1.0391}, {0.2839,0.0116,1.3856}, {0.3285,0.0168,1.6230},
        {0.3483,0.0230,1.7471}, {0.3481,0.0298,1.7826}, {0.3362,0.0380,1.7721},
        {0.3187,0.0480,1.7441}, {0.2908,0.0600,1.6692}, {0.2511,0.0739,1.5281},
        {0.1954,0.0910,1.2876}, {0.1421,0.1126,1.0419}, {0.0956,0.1390,0.8130},
        {0.0580,0.1693,0.6162}, {0.0320,0.2080,0.4652}, {0.0147,0.2586,0.3533},
        {0.0049,0.3230,0.2720}, {0.0024,0.4073,0.2123}, {0.0093,0.5030,0.1582},
        {0.0291,0.6082,0.1117}, {0.0633,0.7100,0.0782}, {0.1096,0.7932,0.0573},
        {0.1655,0.8620,0.0422}, {0.2257,0.9149,0.0298}, {0.2904,0.9540,0.0203},
        {0.3597,0.9803,0.0134}, {0.4334,0.9950,0.0087}, {0.5121,1.0000,0.0057},
        {0.5945,0.9950,0.0039}, {0.6784,0.9786,0.0027}, {0.7621,0.9520,0.0021},
        {0.8425,0.9154,0.0018}, {0.9163,0.8700,0.0017}, {0.9786,0.8163,0.0014},
        {1.0263,0.7570,0.0011}, {1.0567,0.6949,0.0010}, {1.0622,0.6310,0.0008},
        {1.0456,0.5668,0.0006}, {1.0026,0.5030,0.0003}, {0.9384,0.4412,0.0002},
        {0.8544,0.3810,0.0002}, {0.7514,0.3210,0.0001}, {0.6424,0.2650,0.0000},
        {0.5419,0.2170,0.0000}, {0.4479,0.1750,0.0000}, {0.3608,0.1382,0.0000},
        {0.2835,0.1070,0.0000}, {0.2187,0.0816,0.0000}, {0.1649,0.0610,0.0000},
        {0.1212,0.0446,0.0000}, {0.0874,0.0320,0.0000}, {0.0636,0.0232,0.0000},
        {0.0468,0.0170,0.0000}, {0.0329,0.0119,0.0000}, {0.0227,0.0082,0.0000},
        {0.0158,0.0057,0.0000}, {0.0114,0.0041,0.0000}, {0.0081,0.0029,0.0000},
        {0.0058,0.0021,0.0000}, {0.0041,0.0015,0.0000}, {0.0029,0.0010,0.0000},
        {0.0020,0.0007,0.0000}, {0.0014,0.0005,0.0000}, {0.0010,0.0004,0.0000},
        {0.0007,0.0002,0.0000}, {0.0005,0.0002,0.0000}, {0.0003,0.0001,0.0000},
        {0.0002,0.0001,0.0000}, {0.0002,0.0001,0.0000}, {0.0001,0.0000,0.0000},
        {0.0001,0.0000,0.0000}, {0.0001,0.0000,0.0000}, {0.0000,0.0000,0.0000}
    };

    for (i = 0, lambda = 380; lambda < 780.1; i++, lambda += 5) {
        double Me;

        Me = (*spec_intens)(lambda);
        X += Me * cie_colour_match[i][0];
        Y += Me * cie_colour_match[i][1];
        Z += Me * cie_colour_match[i][2];
    }
    XYZ = (X + Y + Z);
    *x = X / XYZ;
    *y = Y / XYZ;
    *z = Z / XYZ;
}

/*                            BB_SPECTRUM

    Calculate, by Planck's radiation law, the emittance of a black body
    of temperature bbTemp at the given wavelength (in metres).  */

static double bbTemp = 5000;                 /* Hidden temperature argument
                                         to BB_SPECTRUM. */
static double
bb_spectrum(double wavelength)
{
    double wlm = wavelength * 1e-9;   /* Wavelength in meters */

    return (3.74183e-16 * pow(wlm, -5.0)) /
           (exp(1.4388e-2 / (wlm * bbTemp)) - 1.0);
}

/*  Built-in test program which displays the x, y, and Z and RGB
    values for black body spectra from 1000 to 10000 degrees kelvin.
    When run, this program should produce the following output:

    Temperature       x      y      z       R     G     B
    -----------    ------ ------ ------   ----- ----- -----
       1000 K      0.6528 0.3444 0.0028   1.000 0.007 0.000 (Approximation)
       1500 K      0.5857 0.3931 0.0212   1.000 0.126 0.000 (Approximation)
       2000 K      0.5267 0.4133 0.0600   1.000 0.234 0.010
       2500 K      0.4770 0.4137 0.1093   1.000 0.349 0.067
       3000 K      0.4369 0.4041 0.1590   1.000 0.454 0.151
       3500 K      0.4053 0.3907 0.2040   1.000 0.549 0.254
       4000 K      0.3805 0.3768 0.2428   1.000 0.635 0.370
       4500 K      0.3608 0.3636 0.2756   1.000 0.710 0.493
       5000 K      0.3451 0.3516 0.3032   1.000 0.778 0.620
       5500 K      0.3325 0.3411 0.3265   1.000 0.837 0.746
       6000 K      0.3221 0.3318 0.3461   1.000 0.890 0.869
       6500 K      0.3135 0.3237 0.3628   1.000 0.937 0.988
       7000 K      0.3064 0.3166 0.3770   0.907 0.888 1.000
       7500 K      0.3004 0.3103 0.3893   0.827 0.839 1.000
       8000 K      0.2952 0.3048 0.4000   0.762 0.800 1.000
       8500 K      0.2908 0.3000 0.4093   0.711 0.766 1.000
       9000 K      0.2869 0.2956 0.4174   0.668 0.738 1.000
       9500 K      0.2836 0.2918 0.4246   0.632 0.714 1.000
      10000 K      0.2807 0.2884 0.4310   0.602 0.693 1.000
*/

#if 0
int main()
{ 
    double t, x, y, z, r, g, b;
    struct colourSystem *cs = &SMPTEsystem;

    printf("Temperature       x      y      z       R     G     B\n");
    printf("-----------    ------ ------ ------   ----- ----- -----\n");
    
    for (t = 1000; t <= 10000; t+= 500) {
        bbTemp = t;
        spectrum_to_xyz(bb_spectrum, &x, &y, &z);
        xyz_to_rgb(cs, x, y, z, &r, &g, &b);
        printf("  %5.0f K      %.4f %.4f %.4f   ", t, x, y, z);
        if (constrain_rgb(&r, &g, &b)) {
	    norm_rgb(&r, &g, &b);
            printf("%.3f %.3f %.3f (Approximation)\n", r, g, b);
        } else {
	    norm_rgb(&r, &g, &b);
            printf("%.3f %.3f %.3f\n", r, g, b);
        }
    }
    return 0;
}
#endif

void
get_rgb (double t, double *red, double *green, double *blue)
{
  double x, y, z, r, g, b;
  
  struct colourSystem *cs = &SMPTEsystem;
  bbTemp = 5000.0 + (t - 0.8) * 20000.0;

  spectrum_to_xyz(bb_spectrum, &x, &y, &z);
  xyz_to_rgb(cs, x, y, z, &r, &g, &b);
  if (red)   *red = r;
  if (green) *green = g;
  if (blue)  *blue = b;
}



#define BASE_STRENGTH 200.0
double
get_strength (int a, int b)
{
  double strength = BASE_STRENGTH;
  
  switch(combo (a, b)) {
  case combo (AN_H, AN_H):
    strength = 432.0;
    break;
  case combo (AN_H, AN_C):
  case combo (AN_C, AN_H):
    strength = 411.0;
    break;
  case combo (AN_H, AN_Si):
  case combo (AN_Si, AN_H):
    strength = 318.0;
    break;
  case combo (AN_H, AN_N):
  case combo (AN_N, AN_H):
    strength = 386.0;
    break;
  case combo (AN_H, AN_P):
  case combo (AN_P, AN_H):
    strength = 322.0;
    break;
  case combo (AN_H, AN_O):
  case combo (AN_O, AN_H):
    strength = 459.0;
    break;
  case combo (AN_H, AN_S):
  case combo (AN_S, AN_H):
    strength = 363.0;
    break;
  case combo (AN_H, AN_F):
  case combo (AN_F, AN_H):
    strength = 565.0;
    break;
  case combo (AN_H, AN_Cl):
  case combo (AN_Cl, AN_H):
    strength = 428.0;
    break;
  case combo (AN_H, AN_Br):
  case combo (AN_Br, AN_H):
    strength = 362.0;
    break;
  case combo (AN_H, AN_I):
  case combo (AN_I, AN_H):
    strength = 295.0;
    break;
  case combo (AN_C, AN_C):
    strength = 346.0;
    break;
  case combo (AN_C, AN_Si):
  case combo (AN_Si, AN_C):
    strength = 318.0;
    break;
  case combo (AN_C, AN_N):
  case combo (AN_Si, AN_N):
    strength = 305.0;
    break;
  case combo (AN_C, AN_O):
  case combo (AN_O, AN_C):
    strength = 358.0;
    break;
  case combo (AN_C, AN_S):
  case combo (AN_S, AN_C):
    strength = 272.0;
    break;
  case combo (AN_C, AN_F):
  case combo (AN_F, AN_C):
    strength = 485.0;
    break;
  case combo (AN_C, AN_Cl):
  case combo (AN_Cl, AN_C):
    strength = 327.0;
    break;
  case combo (AN_C, AN_Br):
  case combo (AN_Br, AN_C):
    strength = 285.0;
    break;
  case combo (AN_C, AN_I):
  case combo (AN_I, AN_C):
    strength = 318.0;
    break;
  case combo (AN_Si, AN_Si):
    strength = 362.0;
    break;
  case combo (AN_Si, AN_O):
  case combo (AN_O, AN_Si):
    strength = 452.0;
    break;
  case combo (AN_N, AN_N):
    strength = 167.0;
    break;
  case combo (AN_N, AN_O):
  case combo (AN_O, AN_N):
    strength = 201.0;
    break;
  case combo (AN_N, AN_F):
  case combo (AN_F, AN_N):
    strength = 283.0;
    break;
  case combo (AN_N, AN_Cl):
  case combo (AN_Cl, AN_N):
    strength = 313.0;
    break;
  case combo (AN_N, AN_Br):
  case combo (AN_Br, AN_N):
    strength = 243.0;
    break;
  case combo (AN_P, AN_P):
    strength = 201.0;
    break;
  case combo (AN_O, AN_O):
    strength = 142.0;
    break;
  case combo (AN_O, AN_F):
  case combo (AN_F, AN_O):
    strength = 190.0;
    break;
  case combo (AN_O, AN_Cl):
  case combo (AN_Cl, AN_O):
    strength = 218.0;
    break;
  case combo (AN_O, AN_Br):
  case combo (AN_Br, AN_O):
    strength = 201.0;
    break;
  case combo (AN_O, AN_I):
  case combo (AN_I, AN_O):
    strength = 201.0;
    break;
  case combo (AN_S, AN_S):
    strength = 226.0;
    break;
  case combo (AN_S, AN_F):
  case combo (AN_F, AN_S):
    strength = 284.0;
    break;
  case combo (AN_S, AN_Cl):
  case combo (AN_Cl, AN_S):
    strength = 255.0;
    break;
  case combo (AN_S, AN_Br):
  case combo (AN_Br, AN_S):
    strength = 218.0;
    break;
  case combo (AN_F, AN_F):
    strength = 155.0;
    break;
  case combo (AN_F, AN_Cl):
  case combo (AN_Cl, AN_F):
    strength = 249.0;
    break;
  case combo (AN_F, AN_Br):
  case combo (AN_Br, AN_F):
    strength = 249.0;
    break;
  case combo (AN_F, AN_I):
  case combo (AN_I, AN_F):
    strength = 278.0;
    break;
  case combo (AN_Cl, AN_Cl):
    strength = 240.0;
    break;
  case combo (AN_Cl, AN_Br):
  case combo (AN_Br, AN_Cl):
    strength = 216.0;
    break;
  case combo (AN_Cl, AN_I):
  case combo (AN_I, AN_Cl):
    strength = 208.0;
    break;
  case combo (AN_Br, AN_Br):
    strength = 190.0;
    break;
  case combo (AN_Br, AN_I):
  case combo (AN_I, AN_Br):
    strength = 175.0;
    break;
  case combo (AN_I, AN_I):
    strength = 149.0;
    break;
  }
  return strength / BASE_STRENGTH;;
}
