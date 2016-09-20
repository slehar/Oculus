//This code is copyright Picturesolve.com
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
//#include <fftw.h>
#include <time.h>
//#include <fftw-int.h>

#ifdef ARRAYSTRUCTBLOCK
#else
#define ARRAYSTRUCTBLOCK 1
/*this is where array_type is defined with two integers and a pointer to the array*/
typedef struct{
	int xsize, ysize;
	fftw_complex *ptr;
}array_type;
#endif


//
//prototyping of alex's deblurring routines
void load_array_complex(array_type *);
void print_array_complex(array_type *);
void scale_array_complex(array_type *);
void shift_array_complex(array_type *);
void hanning_window_complex(array_type *);
int psf_oof_complex(array_type *, float);
void psf_lin_complex(array_type *, float, float);
int minimum_twonums (int, int);
int maximum_twonums (int, int);
fftw_complex conj_c (fftw_complex);
fftw_complex add_c (fftw_complex, fftw_complex);
fftw_complex prod_c (fftw_complex, fftw_complex);
fftw_complex div_c (fftw_complex, fftw_complex);
int mapgen_complex(array_type *, array_type *, float);
void scale_array_complex_convert(array_type *, unsigned char *);
void mult_array_complex(array_type *, array_type *,array_type *);
double set_dc_complex(array_type *, double, int);
void powerspectrum(array_type *, unsigned char *, int, float);
int convolve(array_type *, array_type *, array_type *, int);
void normalize_array_complex(array_type *);
int deblur(array_type *, array_type *, array_type *, float);
void add_noise_complex(array_type *, float);
void make_grayscale_complex(array_type *);
void accumulate_array_complex(array_type *, array_type *);
void zero_array_complex(array_type *);
void magnitude_phase_complex(array_type *);
void complex_to_byte_array(array_type *, unsigned char *, int );
void average_power_spectrum(array_type *, int , int , array_type *, int, int *, int *);
void forward_fft(array_type *);
void inverse_fft(array_type *);
void value_array_complex(array_type *, double , double);
void stretch_complex_array(array_type *, double , double , int );
void convert_to_logthreshold_complex(array_type *, double );
void cepstrum(array_type *, int , int , array_type *, int );
void invert_realpart_complex(array_type *);
void divide_element(array_type *, float );
void clip_realpart_complex(array_type *, double, double);
void maxmin_array(array_type *, double *, double *, int );
void cepest(array_type *, int , int , array_type *, int );
void move_array(array_type *, array_type *);
void edge(array_type *, array_type *, array_type *);
void patgen(array_type *, array_type *, int , int , int );
int cnvolv_spatial_to_spatial(array_type *, array_type *, array_type *, int );
void convert_to_exponential_complex(array_type *);
void add_array_complex(array_type *, double );
void mult_const_complex(array_type *, double );
int deblur1(array_type *, array_type *, array_type *, float );
double getsum_array_complex(array_type *, int );
void normalize_justreals(array_type *);
void maxmin_array_special(array_type *, double *, double *, int , int *, int *);
double get_dc_complex(array_type *, int );
void pattern_removal(array_type *, array_type *, array_type *, int , int , int ,float);
void mult_array_complex_special(array_type *, array_type *, int , int , array_type *);
void remove_nyquist(array_type *, int, int);
void mult_array_otherreal(array_type *, array_type *, array_type *);
void hanning_window_box(array_type *, int , int );
void border_padding(array_type *, int , int , double );
void quadrant_swap(array_type *);
void punch_box(array_type *, int , int , double );
void punch_hole(array_type *, double , double );
void psf_gen(array_type *,float ,float, int);
void psf_linthick_complex(array_type *array, float, float, float);
void psf_genthick(array_type *array, float, float, float, int, float, float);
int psf_ellipse_complex(array_type *array, float, float, float);
int oneIfWithinEllipse( int , int , int, int , int , int , float  );
int radius_squared( int , int );
float theta_1( int , int );
int psf_gaussellipse(array_type *array, float, float, float, float);
int psf_invert_gaussellipse(array_type *array, float, float, float, float, float);
float gaussianEllipseValue( float, int , int , int , int , int , int , float  );
void rescale_array_complex(array_type *array, float);
void make_fuzzy(array_type *array, array_type *brray, float);
