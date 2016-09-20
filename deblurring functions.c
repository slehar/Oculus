#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fftw.h>
//#include <fftw.h>
#include <time.h>
//#include <fftw-int.h>
#include "deblurring functions.h"



//--------------This code is copyright Picturesolve.com-----------------
//-----This is the collection of Alex's Restoration modules-------------
//-----	   												   -------------
//								||
//							   \  /
//								\/
//
//    "Everything should be as simple as possible, but not Simpler!"  
//						Albert Einstein
//
//
//	General notes:  
/*The approach taken here is that fourier domain filtering requires the extensive 
use of in-place Complex FFT's.  Thus complex arithmetic lies at the heart of 
a lot of these operations.  We need to pass images of complex type around, of arbitrary size
and we define a structure which carries ix, iy, and a pointer to the complex image data.  
This data structure is very extensively used in the code below, sometimes carrying redundancy
in the imaginary part for simplicity.  */

void load_array_complex(array_type *array)
/***********************************************************
*	load_array(array)	Alex Lehar								
*	This routine takes a struct as an argument, which contains 
*	Two integers and a pointer to a 2D complex (double) array
*.  It loads the Reals 
*	with input data, and zeros the Imaginaries.
*************************************************************/
{
	int i,j,xsize,ysize;
	fftw_complex putin;
	double a_real, a_imag;
	a_real=0.0;/*initial value*/
	a_imag=0.0;
	c_re(putin)=a_real;
	c_im(putin)=a_imag;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;
	
	/*load the array*/
	for (j=0;j<xsize;j++){/*for every line*/
		for (i=0;i<ysize; i++){
			*(array->ptr+j*ysize+i) = putin;				/*assign complex*/	
		}
	}
	a_real=1.0;
	c_re(putin)=a_real;
	/*load the array again*/
	for (j=0;j<xsize;j++){				/*for every line*/
		for (i=0;i<ysize; i++){
			//r=rand();
			c_re(putin)=(double)a_real;
			*(array->ptr+j*ysize+i) = putin;				/*assign complex*/	
		}
	}
		

}

void shift_array_complex(array_type *array)
/*************************************************************
*	take the input dynamically resizable input complex pair array
*	and 'shift' its quadrants for display of dc at center
*	WARNING...this only works for even dimensioned arrays
**************************************************************/
{
	int i,j,xsize, ysize;
	double scale;
	fftw_complex c;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;

	
	/*calculate the values*/
	for (j=0;j<xsize;j++){
		for(i=0;i<ysize; i++){
			scale=pow((-1.0),(double)(i+j));
			c=*(array->ptr+j*ysize+i);
			c_re(c)=c_re(c)*scale;
			c_im(c)=c_im(c)*scale;
			*(array->ptr+j*ysize+i)=c;
		}
	}
}	

void print_array_complex(array_type *array)
/*************************************************************
*	interpret a dynamically resizable 2-D complex'array' and print
*	its contents
**************************************************************/
/*this function will print an array whose structure is defined*/
{
	int i,j,xsize, ysize;
	double a,b,magn;
	fftw_complex c;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;
	
	/*print out the values*/
	for (j=0;j<xsize;j++){
		/*printf("%2d %2d", xsize, ysize);*/
		for(i=0;i<ysize; i++){
			c=*(array->ptr+j*ysize+i);
			a=c_re(c);
			b=c_im(c);
			magn=sqrt(a*a+b*b);
			/*printf("   %2.2f %2.2f", a, b);*/
			printf("  %2.2f", magn); 
		}
		printf("\n");
	}
}	

void scale_array_complex(array_type *array)
/*************************************************************
*	take the input dynamically resizable input complex pair array
*	and scale it by 1/(xsize*ysize)
**************************************************************/
{
	int i,j,xsize, ysize;
	double scale;
	fftw_complex c;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;
	scale=1.0/(((double)xsize)*((double)ysize));
	
	/*calculate the values*/
	for (j=0;j<xsize;j++){
		for(i=0;i<ysize; i++){
			c=*(array->ptr+j*ysize+i);
			c_re(c)=c_re(c)*scale;
			c_im(c)=c_im(c)*scale;
			*(array->ptr+j*ysize+i)=c;
		}
	}
}	



void hanning_window_complex(array_type *array)
/************************************************************
*	alex lehar	15 oct 01
*	apply a hanning window to a rectangular complex
*	array in the structure array_type which contains dimensions
*	and data pointer.  
*	note that in its present form, it does not quite taper to zero
*	at the center edges.  
***************************************************************/

{
	int i,j,x,y;
	float r,d,a,b,weight;
	fftw_complex c;
	#define PI 3.1415927
	
	/*unpack struct into local variables*/
	x= array->xsize;
	y= array->ysize;
	
	/*go through the array*/
	/*radius of the circumscribed circle*/
	r=sqrt(x*x+y*y)/2.0;
	for (j=0; j<x; j++) {
		for (i=0;i<y; i++){
			/*distance between center and examined point*/
			a=pow((x/2-j),2);
			b=pow((y/2-i),2);
			d=sqrt(a+b);
			/*cosine weighting factor*/
			weight=1.0-cos((PI*(r-d)/2.0)/r);
			c=*(array->ptr+j*y+i);
			c_re(c)=c_re(c)*weight;
			c_im(c)=c_im(c)*weight;	/*weight both real and imaginary parts*/
			*(array->ptr+j*y+i)=c;	/*final reassignment*/
		}
	}
}

int psf_oof_complex(array_type *array,float radius)
/************************************************************
*	alex lehar	17 oct 01
*	Insert into the center of a rectangular complex array structure
*	a flat-top disc function of specified radius, into the real part.
***************************************************************/
{
	int i,j,x,y,z;
	float d,a,b;
	fftw_complex c;
	#define PI 3.1415927
	
	/*unpack struct into local variables*/
	x= array->xsize;
	y= array->ysize;
	z=minimum_twonums(x,y);	/*compute the minimum*/
	/*make sure the radius is not too big for the array*/
	if(2.0*radius < (float)z)
		{
		/*go through the array*/
		/*and set the necessary elements in the circle.*/
		for (j=0; j<x; j++) {
			for (i=0;i<y; i++){
				/*distance between center and examined point*/
				a=pow((x/2-j),2);
				b=pow((y/2-i),2);
				d=sqrt(a+b);
				c=*(array->ptr+j*y+i);
				if	(d < radius)
					{
					c_re(c)=1.0;
					c_im(c)=0.0;
					*(array->ptr+j*y+i)=c;	/*reassignment*/
					}
				else
					{
					c_re(c)=0.0;
					c_im(c)=0.0;
					*(array->ptr+j*y+i)=c;	/*reassignment*/
					}

			}
		}
		}
	else
		return (0);
	}
   
	

int minimum_twonums (int a, int b)
/***************************
*alex lehar 18 oct 01
*return the minumum of two integers
*****************************/
{
	if (a < b)
		return (a);
	else
		return (b);
}

int maximum_twonums (int a, int b)
/***************************
*alex lehar 27 dec 01
*return the maximum of two integers
*****************************/
{
	if (a < b)
		return (b);
	else
		return (a);
}

fftw_complex conj_c (fftw_complex e)
/*********************************************
*	alex lehar 19 oct 01
*	return the complex conjugate of a double 
*	precision complex number
**********************************************/
{
	fftw_complex f;
	double a,b;
	a=c_re(e);
	b=c_im(e);
	c_re(f)=a;
	c_im(f)=-b;
	return (f);
}

fftw_complex add_c (fftw_complex e, fftw_complex f)
/****************************************************
*	alex lehar	19 oct 01
*	return the sum of two double precision complex
*	numbers, using relationship:
*	(a+bj)+(c+dj)=(a+c)+ j(b+d)
*****************************************************/
{
	double a,b,c,d;
	a=c_re(e);
	b=c_im(e);
	c=c_re(f);
	d=c_im(f);
	/*put result in e*/
	c_re(e)=(a+c);
	c_im(e)=(b+d);
	return (e);
}

fftw_complex prod_c (fftw_complex e, fftw_complex f)
/****************************************************
*	alex lehar	19 oct 01
*	return the product of two double precision complex
*	numbers, using relationship:
*	(a+bj)(c+dj)=(ac-bd)+ j(ad+bc)
*****************************************************/
{
	double a,b,c,d;
	a=c_re(e);
	b=c_im(e);
	c=c_re(f);
	d=c_im(f);
	/*put result in e*/
	c_re(e)=(a*c-b*d);
	c_im(e)=(a*d+b*c);
	return (e);
}

fftw_complex div_c (fftw_complex e, fftw_complex f)
/******************************************************
*	alex lehar 19 oct 01
*	return the result of dividing two double precision complex numbers
*	using the relationship:
*					(ac+bd)			(bc-ad) 
*	(a+bj)/(c+dj)=	-------		+	-------		j
*					(c^2 +d^2)		(c^2+d^2)
*******************************************************/
{
	double a,b,c,d;
	a=c_re(e);
	b=c_im(e);
	c=c_re(f);
	d=c_im(f);
	/*put result in e*/
	c_re(e)=(a*c+b*d)/(c*c+d*d);
	c_im(e)=(b*c-a*d)/(c*c+d*d);
	return (e);
}

int mapgen_complex(array_type *array, array_type *brray, float snr)
/************************************************************
*	alex lehar 19 oct 01
*	take the input double precision complex input 'array'
*	assumed to contain the point spread function, using the signal 
*	to noise ratio 'snr', compute the inverse linear MAP restoration
*	filter.  Inverse Transform this to create a Spatial Domain
*	Inverse Restoration Kernel and put the result into 'brray'.  
*	Assume array and brray are identical in size.  The form of this 
*	filter in the Fourier domain is:
*
*				(H*(u,v)+C)
*	K(u,v)	=	------------------    where C=1/snr
*				(H*(u,v).H(u,v)+C)
*
*	note: that array now contains Transform K(u,v) and brray shifted inverse result
*	modified 11 nov, to do the whole job, spatial domain to spatial domain
*************************************************************/
{
	array_type *crray;				/*required internal buffer*/
	int i,j,x,y,x1,y1;
	fftw_complex c,k,s,t,u,w;
	double invsnr;
	invsnr=1.0/snr;
	/*set complex invsnr*/
	c_re(w)=invsnr;
	c_im(w)=0.;
	
	/*unpack struct into local variables*/
	x= array->xsize;
	y= array->ysize;
	x1=brray->xsize;
	y1=brray->ysize;
	
	/*allocate memory for the internal buffer crray*/
	crray = (array_type*)malloc(sizeof(array_type));
	crray->xsize=x;
	crray->ysize=y;
	crray->ptr = (fftw_complex *)malloc(x*y*sizeof(fftw_complex));
	
	if ((x+y) != (x1+y1))
		return (0);	/*escape if images incompatible*/
	else
		{
		move_array(array, crray);			/*to preserve psf in array*/
		shift_array_complex(brray);			/*for quadrant shifting*/
		/*transform the input psf*/
		forward_fft(crray);					/*shifted form*/
		
		/*create from this the Fourier MAP filter*/
			 for (j=0; j<x; j++) {
			 	for (i=0; i<y; i++) {
			 		c=*(crray->ptr+j*y+i);	/*read the element*/
			 		s=add_c(conj_c(c),w);	/*first term in MAP{see journal p 10} */
			 		t=prod_c(conj_c(c),c);	/*intermediate term*/
			 		u=add_c(w,t);			/*	"				*/
			 		k=div_c(s,u);			/*	"				*/
			 		*(brray->ptr+j*y+i)=k;	/*final assignment*/
			 	}
			 }
		move_array(brray, array);			/*put Transform into array for inspection*/		
		inverse_fft(brray);
		///normalize_array_complex(brray);		/*to give correct scaling*/
		}
		
	free(crray->ptr);
	free(crray);							/*release the local memory*/
	}
	
void scale_array_complex_convert(array_type *array, unsigned char *bytAr)
/*************************************************************
*	take the input dynamically resizable input complex pair array
*	and scale the real part 0-255, and then convert data to an
*	identically dimensioned byte array located at *bytAr
**************************************************************/

{
	int i,j,xsize, ysize;
	double r,s,t,min,max;
	fftw_complex c;
	min=10000000.;
	max=-1000000.;	/*starting values*/
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;
	/*find the extrema*/
	for (j=0; j<xsize; j++){
		for(i=0; i<ysize; i++){
			c=*(array->ptr+j*ysize+i);
				r=c_re(c);
				max = (r > max) ? r : max;
				min = (r < min) ? r : min;
		}
	}
	
	s=(max-min);			
				
	/*convert the pixels*/
	for (j=0;j<xsize;j++){
		for(i=0;i<ysize; i++){
			c=*(array->ptr+j*ysize+i);
			r=c_re(c);
			t=255.*(r-min)/s;	/*make sure it is scaled*/
			*(bytAr+j*ysize+i)=(unsigned char)t;
		}
	}
}	

void mult_array_complex(array_type *array, array_type *brray, array_type *crray)
/*************************************************************
*	Alex Lehar	25 oct 01
*	take the input dynamically resizable input complex pair array
*	and (complex) multiply each element with that of brray, to produce
*	crray.  For the time being assume all arrays are the same size.
**************************************************************/

{
	int i,j,xsize, ysize;
	fftw_complex a,b,c;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;

	
	/*calculate the values*/
	for (j=0;j<xsize;j++){
		for(i=0;i<ysize; i++){
			a=*(array->ptr+j*ysize+i);
			b=*(brray->ptr+j*ysize+i);
			c=prod_c(a,b);
			*(crray->ptr+j*ysize+i)=c;
		}
	}
}	


double set_dc_complex(array_type *array, double dc, int icen)
/*************************************************************
*	Alex Lehar	27 Oct 01
*	take the input dynamically resizable input complex pair array
*	and set the d.c. term.  To deal with the possibility that 
*	the array might be shifted or unshifted, look for the point
*	of maximum magnitude, identify that as dc, and set that.  
*	return the value of the dc term before reset
*	Revision 20 Nov 01:  That is prone to error: Introduce extra
*	parameter icen.
*	arguments
*	*array		complex data structure
*	dc			value to set dc at
*	icen		=0 if at array center, =1 if at topleft corner
**************************************************************/
{
	int i,j,xsize, ysize;
	int k,l;						/*indices for 'found' dc term*/
	fftw_complex c;
	double max,r,s,t;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;

	
	/*set dc's indices and value*/
	if (icen == 0)
		{
		k=ysize/2;
		l=xsize/2;
		}
	else
		{
		k=0;
		l=0;
		}
	c=*(array->ptr+l*ysize+k);
	max=c_re(c);	
	c_re(c)=dc;
	c_im(c)=0.;
	*(array->ptr+l*ysize+k)=c;
	return max;
}	

void powerspectrum(array_type *array, unsigned char *bytAr, int logornot, float level)
/*************************************************************
*	alex lehar	28 oct 01
*	take the input dynamically resizable input complex pair array
*	and compute the power spectrum, scale scale it 0-255, and then 
*	transfer data to an identically dimensioned byte array located 
*	at *bytAr.  If logornot=1 magnitude, if logornot=2, log magnitude
*	result. Clip the data at low threshold level of LEVEL.
*	arguments
*	*array		pointer to complex array structure containing fourier transform
*	*bytAr		pointer to byte array containing 0-255 output values
*	logornot	flag 1 or 2, if Magnitude or Log Magnitude is needed
*	level		lower magnitude threshold to set on data
**************************************************************/

{
	int i,j,xsize, ysize;
	float r,s,t,u,v,min,max;
	fftw_complex c;
	min=1e12;
	max=-1e12;	/*starting values*/
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;
	
	/*examine one of two possibilities, magnitude or logmagnitude*/
	switch (logornot) 
	{
		case 1:		/*just the magnitude*/
			{
				/*find the extrema*/
				for (j=0; j<xsize; j++){
					for(i=0; i<ysize; i++){
						c=*(array->ptr+j*ysize+i);
							r=c_re(c);
							s=c_im(c);
							t=sqrt(r*r+s*s);			/*magnitude*/
							max=(t > max) ? t : max;	/*concise if form*/
							min=(t < min) ? t : min;
					}
				}
				
				u=(max-min);			
							
				/*compute output pixels*/
				for (j=0;j<xsize;j++){
					for(i=0;i<ysize; i++){
						c=*(array->ptr+j*ysize+i);
						r=c_re(c);
						s=c_im(c);
						t=(sqrt(r*r+s*s) > level) ? sqrt(r*r+s*s) : level; /*cut off low thresh*/
						v=255.*(t-min)/u;	/*make sure it is scaled*/
						*(bytAr+j*ysize+i)=(unsigned char)v;
					}
				}
			}
			break;

		case 2:		/*log magnitude*/
			{
				/*find the extrema*/
				for (j=0; j<xsize; j++){
					for(i=0; i<ysize; i++){
						c=*(array->ptr+j*ysize+i);
							r=c_re(c);
							s=c_im(c);
							t=log(sqrt(r*r+s*s));			/*magnitude floor*/
							max=(t > max) ? t : max;		/*concise if form*/
							min=(t < min) ? t : min;
					}
				}
				
				u=(max-min);			
							
				/*compute output pixels*/
				for (j=0;j<xsize;j++){
					for(i=0;i<ysize; i++){
						c=*(array->ptr+j*ysize+i);
						r=c_re(c);
						s=c_im(c);
						t=(log(sqrt(r*r+s*s)) > level) ? log(sqrt(r*r+s*s)) : log(level);/*magnitude floor*/
						v=255.*(t-min)/u;									 /*make sure it is scaled*/
						*(bytAr+j*ysize+i)=(unsigned char)v;
					}
				}
			}
			break;
		}
	}	

	int convolve(array_type *array, array_type *brray, array_type *crray, int window)
	/************************************************************
	*	alex lehar 29 oct 01
	*	take the input double precision complex input 'array'
	*	and similar 'brray', each containing a spatial domain data
	*	in the Real part, and zero in Imaginary; and convolve them
	*	to produce an output 'crray', containing the spatial domain
	*	result in its real part.  'brray' is normalized first.
	*	All this is accomplished by
	*	multiplication in the intermediate Fourier domain according
	*	to the Convolution Theorem. 
	*	If Window=0, use no windowing, if 1, use Hanning Window. 
	*	note that 'array' is corrupted in the process!
	*************************************************************/
	{
		int xa,ya,xb,yb,xc,yc,flags;
		fftw_direction isign;
		fftwnd_plan p;
		flags  = FFTW_IN_PLACE;
		
		/*unpack struct into local variables*/
		xa= array->xsize;
		ya= array->ysize;
		xb=brray->xsize;
		yb=brray->ysize;
		xc=crray->xsize;
		yc=crray->ysize;
		
		if (((xa+ya) + (xb+yb) + (xc+yc) ) != 3*(xa+ya))	/*incompatible sizes*/
			return (0);	/*escape if images incompatible*/
		else
			{
			if(window == 0)
				{
				isign=FFTW_FORWARD;  				/*Forward Fourier Transform*/
				p= fftw2d_create_plan(xa,ya,isign,flags);
				fftwnd_one(p, array->ptr, NULL); 	/*Transform*/
				fftwnd_one(p, brray->ptr, NULL); 	/*Transform*/
				normalize_justreals(brray);		/*normalization*/
				mult_array_complex(array, brray, crray);	/*multiply=convolution*/
				isign=FFTW_BACKWARD;				/*Inverse Fourier Transform*/
				p= fftw2d_create_plan(xa,ya,isign,flags);
				shift_array_complex(crray);
				fftwnd_one(p, crray->ptr, NULL);
				scale_array_complex(crray);			/*scale back*/
				}
			else if (window ==1)
				{
				hanning_window_complex(array);		/*apply window*/
				isign=FFTW_FORWARD;  				/*Forward Fourier Transform*/
				p= fftw2d_create_plan(xa,ya,isign,flags);
				fftwnd_one(p, array->ptr, NULL); 	/*Transform*/
				hanning_window_complex(brray);		/*apply window*/
				fftwnd_one(p, brray->ptr, NULL); 	/*Transform*/
				normalize_justreals(brray);			/*normalization*/
				mult_array_complex(array, brray, crray);	/*multiply=convolution*/
				isign=FFTW_BACKWARD;				/*Inverse Fourier Transform*/
				p= fftw2d_create_plan(xa,ya,isign,flags);
				shift_array_complex(crray);
				fftwnd_one(p, crray->ptr, NULL);
				scale_array_complex(crray);			/*scale back*/
				}
			}
	}

void normalize_array_complex(array_type *array)
/*************************************************************
*	alex lehar 30 oct 01
*	take the input dynamically resizable input complex pair array
*	and normalize it by 1/(sum of elements)
*	modified 9 nov, to involve sums of magnitudes.
**************************************************************/
{
	int i,j,xsize,ysize;
	fftw_complex c;
	double sum,r,s;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;
	
	/*initialize the sum*/
	sum=0.;
		
	/*find the sum for scale*/
	for (j=0;j<xsize;j++){
		for(i=0;i<ysize; i++){
			c=*(array->ptr+j*ysize+i);
			r=c_re(c);
			s=c_im(c);
			sum=sum+sqrt(r*r+s*s);		/*accumulate magnitudes*/
		}
	}
	
	/*do the normalization*/
	for (j=0;j<xsize;j++){
		for(i=0;i<ysize; i++){
			c=*(array->ptr+j*ysize+i);
			r=c_re(c);
			s=c_im(c);
			c_re(c)=r/sum;
			c_im(c)=s/sum;
			*(array->ptr+j*ysize+i)=c;
		}
	}
	
}	

	int deblur(array_type *array, array_type *brray, array_type *crray, float snr)
	/************************************************************
	*	alex lehar 30 oct 01
	*	array	Input Picture Complex Array Structure
	*	brray	Input PSF Complex Array Structure
	*			Output contains Restored Picture
	*	crray	Output Linear MAP Kernel.
	*	snr		Estimated Signal to Noise Ratio (100 typical)
	*	implements Fourier Domain Restoration, using the form of the 
	*	Convolution operation Spatial-Fourier-BacktoSpatial
	*	WARNING.  You might have to use windowing option (=1) in 
	*	final convolution to control ringing in the picture.  
	*************************************************************/
	{
		int xa,ya,xb,yb,xc,yc;
		
		/*unpack struct into local variables*/
		xa= array->xsize;
		ya= array->ysize;
		xb=brray->xsize;
		yb=brray->ysize;
		xc=crray->xsize;
		yc=crray->ysize;

		
		if (((xa+ya) + (xb+yb) + (xc+yc) ) != 3*(xa+ya))	/*incompatible sizes*/
			return (0);										/*escape if images incompatible*/
		else
			{
			mapgen_complex(brray, crray, snr);				/*MAP kernel in crray*/
			//convolve(array, crray, brray, 0);				/*spatial domain convolution*/
															/*result sits in brray overwriting psf*/
			cnvolv_spatial_to_spatial(array, crray, brray, 0);/*more up to date version*/
			}
	}

void add_noise_complex(array_type *array, float snr)
/************************************************************
*	alex lehar	15 oct 01
*	assume that complex array structure pointed to by 'array'
*	contains real data.  Go through the array and add 'snr' 
*	random noise to each pixel. 'snr' is really noise-to-signal ratio here
***************************************************************/
{
	int i,j,x,y,seed;
	double r,s,t, hi, lo, v;
	long int M;
	fftw_complex c;
	seed=10000;		//choose a seed value
	srand(seed);	//initialize random number generator
	
	maxmin_array(array, &hi, &lo, 0);	//get the range of real values
	v=hi-lo;
	
	/*unpack struct into local variables*/
	x= array->xsize;
	y= array->ysize;
	
	/*go through the array*/
	for (j=0; j<x; j++) {
		for (i=0;i<y; i++){
			c=*(array->ptr+j*y+i);
			r=c_re(c);
			s=((double)rand()/((double)(RAND_MAX)+(double)(1)));
			s=r+v*snr*(s-0.5);
			c_re(c)=s;
			*(array->ptr+j*y+i)=c;	/*final reassignment*/
		}
	}
}

void make_grayscale_complex(array_type *array)
/************************************************************
*	alex lehar	15 oct 01
*	overwrites a grayscale ramp set in the middle of the complex
*	array, for testing purposes.
***************************************************************/
{
	int i,j,x,y,yh,k;
	fftw_complex c;
	
	/*unpack struct into local variables*/
	x= array->xsize;
	y= array->ysize;
	yh=y/2;

	
	/*go through the array, just central region*/
	for (j=x-10; j<x; j++) {
		for (i= (yh-128),k=0;i < (yh+127); i++, k++){	/*k indexes grey ramp*/
			c_re(c)=(double)k;
			c_im(c)=0.;
			*(array->ptr+j*y+i)=c;	/*final reassignment*/
		}
	}
}

void accumulate_array_complex(array_type *array, array_type *brray)
/*************************************************************
*	Alex Lehar	25 oct 01
*	take the input dynamically resizable input complex pair array
*	and add the contents of 'array' into those of 'brray'
**************************************************************/
{
	int i,j,xsize,ysize;
	fftw_complex a,b,c;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;

	
	/*calculate the values*/
	for (j=0; j<xsize; j++){
		for(i=0; i<ysize; i++){
			a=*(array->ptr+j*ysize+i);
			b=*(brray->ptr+j*ysize+i);
			c=add_c(a,b);
			*(brray->ptr+j*ysize+i)=c;
		}
	}
}	

void zero_array_complex(array_type *array)
/*************************************************************
*	take the input dynamically resizable input complex pair array
*	and zero it.
**************************************************************/
{
	int i,j,xsize, ysize;
	fftw_complex c;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;
	
	/*set the values*/
	for (j=0;j<xsize;j++){
		for(i=0;i<ysize; i++){
			c_re(c)=0.0;
			c_im(c)=0.0;
			*(array->ptr+j*ysize+i)=c;
		}
	}
}	

void magnitude_phase_complex(array_type *array)
/*************************************************************
*	alex lehar	1 nov 01
*	take the input dynamically resizable input complex pair array
*	compute the magnitude {sqrt(x^2+y^2)} of each element, and store the result
*	in the real part.  Phase {atan(y/x)} is stored in the imaginary part.
**************************************************************/
{
	int i,j,xsize, ysize;
	double x,y,mag,pha;
	fftw_complex c;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;
	
	/*set the values*/
	for (j=0;j<xsize;j++){
		for(i=0;i<ysize; i++){
			c=*(array->ptr+j*ysize+i);
			x=c_re(c);
			y=c_im(c);
			mag=sqrt(x*x+y*y);
			pha=atan(y/x);
			c_re(c)=mag;
			c_im(c)=pha;
			*(array->ptr+j*ysize+i)=c;
		}
	}
}	

void complex_to_byte_array(array_type *array, unsigned char *bytAr, int realorimag)
/**********************************************************************************
*	alex lehar	1 nov 01
*	Input an input dynamically resizable input complex pair array, pointed to by
*	*array,  and choose whether to deal with the real or imaginary part, for 
*	realorimag= 0 or 1 respectively.  Then scale the values found, and 
*	transfer them to a similarly dimensioned byte array located at *bytAr.
***********************************************************************************/
{
	int i,j,xsize,ysize;
	double r,u,v,min,max;
	fftw_complex c;
	min=1e12;
	max=-1e12;	/*starting values*/
		
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;
	
	
	/*examine one of two possibilities, realorimag = 0 or 1  */
	switch (realorimag) 
	{
		case 0:		/*MAGNITUDE TRANSFER*/
			{
				/*find the extrema*/
				for (j=0; j<xsize; j++){
					for(i=0; i<ysize; i++){
						c=*(array->ptr+j*ysize+i);
							r=c_re(c);
							max=(r > max) ? r : max;	/*concise if form*/
							min=(r < min) ? r : min;
					}
				}
				
				u=(max-min);			
							
				/*compute output pixels*/
				for (j=0;j<xsize;j++){
					for(i=0;i<ysize; i++){
						c=*(array->ptr+j*ysize+i);
						r=c_re(c);
						v=255.*(r-min)/u;	/*make sure it is scaled*/
						*(bytAr+j*ysize+i)=(unsigned char)v;
					}
				}
			}
			break;

		case 2:		/*IMAGINARY TRANSFER */
			{
				/*find the extrema*/
				for (j=0; j<xsize; j++){
					for(i=0; i<ysize; i++){
						c=*(array->ptr+j*ysize+i);
							r=c_im(c);
							max=(r > max) ? r : max;	/*concise if form*/
							min=(r < min) ? r : min;
					}
				}
				
				u=(max-min);			
							
				/*compute output pixels*/
				for (j=0;j<xsize;j++){
					for(i=0;i<ysize; i++){
						c=*(array->ptr+j*ysize+i);
						r=c_im(c);
						v=255.*(r-min)/u;	/*make sure it is scaled*/
						*(bytAr+j*ysize+i)=(unsigned char)v;
					}
				}
			}
			break;

		}
	}
void forward_fft(array_type *array)
/**********************************************
*	alex lehar 2 nov 01
*	set up and perform forward in-place fft
*	on complex array structure 'array'
*	confirmed works 2 nov 01
***********************************************/

{
	int xa, ya, flags;
	fftw_direction isign;
	fftwnd_plan p;
	flags = FFTW_IN_PLACE;
	xa= array->xsize;
	ya= array->ysize;
	
	isign=FFTW_FORWARD;
	p=fftw2d_create_plan(xa,ya,isign,flags);
	fftwnd_one(p, array->ptr, NULL);
	//scale_array_complex(array);
}


void inverse_fft(array_type *array)
/**********************************************
*	alex lehar 2 nov 01
*	set up and perform inverse in-place fft
*	on complex array structure 'array'.
*	scale the result to remove nxm bias factor
***********************************************/

{
	int xa, ya,flags;
	fftw_direction isign;
	fftwnd_plan p;
	flags = FFTW_IN_PLACE;
	xa= array->xsize;
	ya= array->ysize;
	//scale_array_complex(array);
	isign=FFTW_BACKWARD;
	p=fftw2d_create_plan(xa,ya,isign,flags);
	fftwnd_one(p, array->ptr, NULL);
	scale_array_complex(array);
}

void average_power_spectrum(array_type *array, int xboxn, int yboxn, array_type *drray, int ovlp,
	int *jx, int *iy)
	/************************************************************
	*	alex lehar 1 nov 01
	*	array	Input Picture Complex Array Structure
	*	xboxn	number of non-overlapping boxes in vertical   x dir'n
	*	yxoxn	number of non-overlapping boxes in horizontal y dir'n
	*	drray	Contains the output, in same format as 'array', but since
	*			the output is smaller, it is embedded into the center
	*			of the larger format.  
	*	ovlp	=0 if non-overlapping, =1 if 50% overlap
	*	*jx		returned address of box row numbers
	*	*iy		returned address of box column numbers
	*
	*	worked first time 5 nov 01, non-overlapping case
	*			also overlapping case!
	*	modified 23 nov to return box dimensions by reference.
	*************************************************************/
	{
		array_type *brray, *crray;							/*necessary internal buffers*/
		/*indices are;i,j the big array &drray; e,f small arrays brray & crray;, m,n the box counters*/
		int i,j,xa,ya,xb,yb,xc,yc,m,n,e,f;
		fftw_complex c;
		
		/*unpack struct into local variables*/
		xa= array->xsize;
		ya= array->ysize;
		
		/*create the box sub-array 'brray' and accumulator buffer 'crray'*/
		/*allocate the array struct*/
		/*this allocates only the space for two ints and a pointer*/
		brray = (array_type*)malloc(sizeof(array_type));/*for the output*/
		crray = (array_type*)malloc(sizeof(array_type));/*for the accumulator*/
		
		brray->xsize = xa/xboxn;							/*determine box size from image*/
		brray->ysize = ya/yboxn;							/*               "             */
		
		*jx=xa/xboxn;										/*box sizes for passing out*/
		*iy=ya/yboxn;
		
		crray->xsize = xa/xboxn;
		crray->ysize = ya/yboxn;

		/*confirm size of the box*/
		xb=brray->xsize;
		yb=brray->ysize;
		
		/*confirm size of the box*/
		xc=crray->xsize;
		yc=crray->ysize;
		
		/*allocate the memory for the box*/		
		brray->ptr = (fftw_complex *)malloc(xb*yb*sizeof(fftw_complex));
		crray->ptr = (fftw_complex *)malloc(xc*yc*sizeof(fftw_complex));

		/*zero the buffers*/
		zero_array_complex(crray);
		zero_array_complex(brray);
		/*set up the loop which will compute power spectrum in each block, load data
		into brray, and then accumulate them in crray, thus averaging them*/
		if (ovlp == 0)											/*non-overlapping case*/
			{
				/*loop through each box (tiles butting against eachother)*/
				for (m=0; m<yboxn; m++){
					for (n=0; n<xboxn; n++){
						/*transfer data from array to brray*/
						for (j= m*xb,f=0; j<(m+1)*xb; j++,f++){
							for (i= n*yb,e=0; i<(n+1)*yb; i++,e++){
								c=*(array->ptr+j*ya+i);				/*read pixel in big array*/
								*(brray->ptr+yb*f+e)=c;				/*write pixel to small box*/
							}
						}
						hanning_window_complex(brray);				/*do the edging*/
						shift_array_complex(brray);
						forward_fft(brray); 						/*FFT the Box*/
						magnitude_phase_complex(brray);				/*polar form*/
						accumulate_array_complex(brray,crray);		/*save the sum in buffer*/
					}
				}
			}
		else
			{
				/*loop through more overlapping boxes*/
				for (m=0; m<(2*yboxn-1); m++){						/*because there are 2*m-1 of them*/
					for (n=0; n<(2*xboxn-1); n++){					/*		"			2*n-1	"	 */
						/*transfer data from array to brray*/
						for (j= m*xb/2,f=0; f<xb; j++,f++){
							for (i= n*yb/2,e=0; e<yb; i++,e++){
								c=*(array->ptr+j*ya+i);				/*read pixel in big array*/
								*(brray->ptr+yb*f+e)=c;				/*write pixel to small box*/
							}
						}
						hanning_window_complex(brray);				/*do the edging*/
						shift_array_complex(brray);
						forward_fft(brray); 						/*FFT the Box*/
						magnitude_phase_complex(brray);				/*polar form*/
						accumulate_array_complex(brray,crray);		/*save the sum in buffer*/
					}
				}
			}
		
		
		/*having completed the accumulation process, now embed the result into the input image
		format of the output array*/
		value_array_complex(drray, 0.,0.);					/*first zero the target array*/
		convert_to_logthreshold_complex(crray,.001);
		//stretch_complex_array(crray, 255., 0., 0);
		/*it is very important to get the positioning of this right because 
		of subsequent processes involving matching of different sized image 
		operations!*/

		for (j=(xa/2)-(xb/2),f=0; f<xb;j++,f++){
			for(i=(ya/2)-(yb/2),e=0; e<yb; i++,e++){
				c=*(drray->ptr+j*ya+i);
				*(drray->ptr+j*ya+i)=*(crray->ptr+f*yc+e);
			}
		}
	free(brray->ptr);
	free(crray->ptr);											/*release temporary memory allocation*/		
	free(brray);
	free(crray);											/*release temporary memory allocation*/		
		
	}


void value_array_complex(array_type *array, double real, double imaginary)
/*************************************************************
*	take the input dynamically resizable input complex pair array
*	and set the real and imaginary values of each element to given
*	values.
**************************************************************/
{
	int i,j,xsize, ysize;
	fftw_complex c;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;
	
	/*set the values*/
	for (j=0;j<xsize;j++){
		for(i=0;i<ysize; i++){
			c_re(c)=real;
			c_im(c)=imaginary;
			*(array->ptr+j*ysize+i)=c;
		}
	}
}	

void stretch_complex_array(array_type *array, double high, double low, int r_or_i)
/**********************************************************************************
*	alex lehar	5 nov 01
*	Input an input dynamically resizable input complex pair array, pointed to by
*	*array,  and choose whether to deal with the real or imaginary part, for 
*	r_or_i= 0 or 1 respectively.  Then scale the values found, so that either the
*	real or imaginary parts respectively scale from 'low' to 'high'
***********************************************************************************/
{
	int i,j,xsize,ysize;
	double r,u,v,min,max;
	fftw_complex c,d,e;
	min=1e12;
	max=-1e12;	/*starting values*/
	v=high-low;	/*desired range*/
		
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;
	
	
	/*examine one of two possibilities, realorimag = 0 or 1  */
	switch (r_or_i) 
	{
		case 0:		/*REAL SCALED*/
			{
				/*find the extrema*/
				for (j=0; j<xsize; j++){
					for(i=0; i<ysize; i++){
						c=*(array->ptr+j*ysize+i);
							r=c_re(c);
							max=(r > max) ? r : max;	/*concise if form*/
							min=(r < min) ? r : min;
					}
				}
				
				u=(max-min);			
							
				/*compute output pixels*/
				for (j=0;j<xsize;j++){
					for(i=0;i<ysize; i++){
						c=*(array->ptr+j*ysize+i);
						c_re(d)=v*(c_re(c)-min)/u;
						c_im(d)=0.0;
						*(array->ptr+j*ysize+i)=d;
					}
				}
			}
			break;

		case 2:		/*IMAGINARY SCALED */
			{
				/*find the extrema*/
				for (j=0; j<xsize; j++){
					for(i=0; i<ysize; i++){
						c=*(array->ptr+j*ysize+i);
							r=c_im(c);
							max=(r > max) ? r : max;	/*concise if form*/
							min=(r < min) ? r : min;
					}
				}
				
				u=(max-min);			
							
				/*compute output pixels*/
				for (j=0;j<xsize;j++){
					for(i=0;i<ysize; i++){
						c=*(array->ptr+j*ysize+i);
						c_re(d)=0.0;
						c_im(d)=v*(c_im(c)-min)/u;
						*(array->ptr+j*ysize+i)=d;
					}
				}
			}
			break;

		}
	}

void convert_to_logthreshold_complex(array_type *array, double thresh)
/*************************************************************
*	alex lehar	5 nov 01
*	take the input dynamically resizable input complex pair array
*	convert all reals and imaginaries to their (natural) logarithm, applying
*	a lower threshold 'thresh', to prevent -infinity problem
*
*	modifications: (1)
*	try setting all imaginaries to zero instead
**************************************************************/
{
	int i,j,xsize, ysize;
	double x,y,p,q;
	fftw_complex c,d;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;
	
	/*set the values*/
	for (j=0;j<xsize;j++){
		for(i=0;i<ysize; i++){
			c=*(array->ptr+j*ysize+i);
			x=c_re(c);
			y=c_im(c);
				if (x > thresh)
					p=log(x);
				else
					p=thresh;
				if (y > thresh)
					q=log(y);
				else
					q=thresh;
					
			c_re(d)=p;
			c_im(d)=q;
			*(array->ptr+j*ysize+i)=d;
		}
	}
}	

void cepstrum(array_type *array, int xboxn, int yboxn, array_type *drray, int ovlp)
	/************************************************************
	*	alex lehar 1 nov 01
	*	this variant accumulates the average power spectrum first, and 
	*	then computes the 'power cepstrum' from that at the end. A lot of
	*	careful thresholding and clipping goes on, to try to distill the 
	*	contribution of 'zeros' structures in the power spectrum.
	*	array	Input Picture Complex Array Structure
	*	xboxn	number of non-overlapping boxes in vertical   x dir'n
	*	yxoxn	number of non-overlapping boxes in horizontal y dir'n
	*	drray	Contains the output, in same format as 'array', but since
	*			the output is smaller, it is embedded into the center
	*			of the larger format.  
	*	ovlp	=0 if non-overlapping, =1 if 50% overlap
	*	based on average_power_spectrum.
	*	worked first 9 nov 01.
	*************************************************************/
	{
		array_type *brray, *crray;							/*necessary internal buffers*/
		/*indices are;i,j the big array &drray; e,f small arrays brray & crray;, m,n the box counters*/
		int i,j,xa,ya,xb,yb,xc,yc,m,n,e,f,total_boxes;
		fftw_complex c;
		double dc, hilev, lolev,range;
		
		/*unpack struct into local variables*/
		xa= array->xsize;
		ya= array->ysize;
		
		/*create the box sub-array 'brray' and accumulator buffer 'crray'*/
		/*allocate the array struct*/
		/*this allocates only the space for two ints and a pointer*/
		brray = (array_type*)malloc(sizeof(array_type));/*for the output*/
		crray = (array_type*)malloc(sizeof(array_type));/*for the accumulator*/
		
		brray->xsize = xa/xboxn;							/*determine box size from image*/
		brray->ysize = ya/yboxn;							/*               "             */
		
		crray->xsize = xa/xboxn;
		crray->ysize = ya/yboxn;

		/*confirm size of the box*/
		xb=brray->xsize;
		yb=brray->ysize;
		
		/*confirm size of the box*/
		xc=crray->xsize;
		yc=crray->ysize;
		
		/*allocate the memory for the box*/		
		brray->ptr = (fftw_complex *)malloc(xb*yb*sizeof(fftw_complex));
		crray->ptr = (fftw_complex *)malloc(xc*yc*sizeof(fftw_complex));

		/*zero the buffers*/
		zero_array_complex(crray);
		zero_array_complex(brray);
		/*set up the loop which will compute power spectrum in each block, load data
		into brray, and then accumulate them in crray, thus averaging them*/
		if (ovlp == 0)											/*non-overlapping case*/
			{
				/*loop through each box (tiles butting against eachother)*/
				for (m=0; m<yboxn; m++){
					for (n=0; n<xboxn; n++){
						/*transfer data from array to brray*/
						for (j= m*xb,f=0; j<(m+1)*xb; j++,f++){
							for (i= n*yb,e=0; i<(n+1)*yb; i++,e++){
								c=*(array->ptr+j*ya+i);				/*read pixel in big array*/
								*(brray->ptr+yb*f+e)=c;				/*write pixel to small box*/
							}
						}
						hanning_window_complex(brray);				/*do the edging*/
						shift_array_complex(brray);
						forward_fft(brray); 						/*FFT the Box*/
						magnitude_phase_complex(brray);				/*polar form*/
						accumulate_array_complex(brray,crray);		/*save the sum in buffer*/
					}
				}
			}
		else
			{
				/*loop through more overlapping boxes*/
				for (m=0; m<(2*yboxn-1); m++){						/*because there are 2*m-1 of them*/
					for (n=0; n<(2*xboxn-1); n++){					/*		"			2*n-1	"	 */
						/*transfer data from array to brray*/
						for (j= m*xb/2,f=0; f<xb; j++,f++){
							for (i= n*yb/2,e=0; e<yb; i++,e++){
								c=*(array->ptr+j*ya+i);				/*read pixel in big array*/
								*(brray->ptr+yb*f+e)=c;				/*write pixel to small box*/
							}
						}
						hanning_window_complex(brray);				/*do the edging*/
						shift_array_complex(brray);
						forward_fft(brray); 						/*FFT the Box*/
						magnitude_phase_complex(brray);				/*polar form*/
						accumulate_array_complex(brray,crray);		/*save the sum in buffer*/
					}
				}
			}
		
		/*record the total number of boxes processed*/
		total_boxes=m*n;
		divide_element(crray, (float)(total_boxes));			/*effectuate uggh! the average*/
		
		/*having completed the accumulation process, now embed the result into the input image
		format of the output array*/
		value_array_complex(drray, 0.,0.);					/*first zero the target array*/

		convert_to_logthreshold_complex(crray,1.0e-34);		/*log & threshold critical for zeros detection*/

		invert_realpart_complex(crray);						/*inversion*/
		maxmin_array(crray, &hilev, &lolev, 0);				/*find extrema*/
		range=hilev-lolev;									/*range of values*/
		clip_realpart_complex(crray, lolev+(.5*range),
				hilev-(.00001*range));							/*this clip level is very CRITICAL */

		stretch_complex_array(crray, 255., 0.,0);			/*scale to normal range*/

		shift_array_complex(crray);
		inverse_fft(crray);
		//dc=set_dc_complex(crray, 0.,0);
		magnitude_phase_complex(crray);
		convert_to_logthreshold_complex(crray, 0.001);
		dc=set_dc_complex(crray, 0.,0);
		stretch_complex_array(crray, 255., 0., 0);

		for (j=(xa/2)-(xb/2),f=0; f<xb;j++,f++){
			for(i=(ya/2)-(yb/2),e=0; e<yb; i++,e++){
				c=*(drray->ptr+j*ya+i);
				*(drray->ptr+j*ya+i)=*(crray->ptr+f*yc+e);
			}
		}
	free(brray->ptr);
	free(crray->ptr);											/*release temporary memory allocation*/		
	free(brray);
	free(crray);											/*release temporary memory allocation*/		
		
	}

void invert_realpart_complex(array_type *array)
/*************************************************************
*	alex lehar	5 nov 01
*	take the input dynamically resizable input complex pair array
*	and just invert all the real values. (specialized application 
*	in cepstrum calculation)
*
*	modification (1)
**************************************************************/
{
	int i,j,xsize, ysize;
	double x,y,p;
	fftw_complex c,d;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;
	
	/*set the values*/
	for (j=0;j<xsize;j++){
		for(i=0;i<ysize; i++){
			c=*(array->ptr+j*ysize+i);
			x=c_re(c);
			y=c_im(c);
			p=-x;				/*do the inversion*/
			c_re(d)=p;
			c_im(d)=0.;			/*(1) changed from setting to y*/
			*(array->ptr+j*ysize+i)=d;
		}
	}
}	


void divide_element(array_type *array, float e)
/************************************************************
*	alex lehar	15 oct 01
*	divide each element of complex array structure pointed to by 'array'
*	by 'e', treating as complex divide for real and imaginary
*	parts.
***************************************************************/
{
	int i,j,x,y;
	double r,s,t,inv;
	fftw_complex b,c,d;
	t=32768.;
	
	/*unpack struct into local variables*/
	x= array->xsize;
	y= array->ysize;
	
	/*go through the array*/
	for (j=0; j<x; j++) {
		for (i=0;i<y; i++){
			c=*(array->ptr+j*y+i);
			c_re(d)=e;
			c_im(d)=0.0;
			b=div_c(c,d);
			*(array->ptr+j*y+i)=b;	/*final reassignment*/
		}
	}
}
void clip_realpart_complex(array_type *array, double lolevel, double hilevel)
/*************************************************************
*	alex lehar	8 nov 01
*	take the input dynamically resizable input complex pair array
*	and operate just on the real part, to clip all values between 
*	lolevel and hilevel.  Leave imaginary 
*	part alone.  Useful in cepstral preprocessing. 
**************************************************************/
{
	int i,j,xsize, ysize;
	double x;
	fftw_complex c;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;
	
	/*set the values*/
	for (j=0;j<xsize;j++){
		for(i=0;i<ysize; i++){
			c=*(array->ptr+j*ysize+i);
			x=c_re(c);
			if (x < lolevel)		/*do the low clipping*/
				x=lolevel;
			else if (x > hilevel)
				x=hilevel;			/*do the high side*/

			c_re(c)=x;
			*(array->ptr+j*ysize+i)=c;
		}
	}
}	
#if 0
void auto_clipper(array_type *array, double lowpercentile, double hipercentile)
/*************************************************************
*	this is work in progress....
*	alex lehar	8 nov 01
*	take the input dynamically resizable input complex pair array
*	and operate just on the real part, to clip all values below a calculated
*	level (low) and clip all values above a calculated level (high).
*	'Low' and 'High' are determined to by 'lowpercentile' and 'hipercentile' 
*	respectively, based on the histogram of 'n' bins.  Leave imaginary 
*	part alone.  Useful in cepstral preprocessing. 
*	arguments
*	*array			pointer to complex array structure containing data
*	lowpercentile	value ranging from 0. to 100.
*	hipercentile	value ranging from 0. to 100.
**************************************************************/
{
	int i,j,xsize, ysize, n;
	double low, high;
	double x;
	fftw_complex c;
	
	n=1000;						/*number of bins*/
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;
	
	
	
	/*set the values*/
	for (j=0;j<xsize;j++){
		for(i=0;i<ysize; i++){
			c=*(array->ptr+j*ysize+i);
			x=c_re(c);
			if (x < level)		/*do the clipping*/
				x=level;

			c_re(c)=x;
			*(array->ptr+j*ysize+i)=c;
		}
	}
}	
#endif 


void maxmin_array(array_type *array, double *high, double *low, int r_or_i)
/**********************************************************************************
*	alex lehar	8 nov 01
*	Input an input dynamically resizable input complex pair array, pointed to by
*	*array,  and choose whether to deal with the real or imaginary part, for 
*	r_or_i= 0 or 1 respectively.  Then determine the max (high) and min (low)
*	values, and return these values to the addresses of 'high' and 'low'
*	method of calling is: 
*			maxmin_array(array, &high, &low);
***********************************************************************************/
{
	int i,j,xsize,ysize;
	double r,min,max;
	fftw_complex c;
	min=1e34;
	max=-1e34;	/*starting values*/
		
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;
	
	
	/*examine one of two possibilities, realorimag = 0 or 1  */
	switch (r_or_i) 
	{
		case 0:		/*REAL SCALED*/
			{
				/*find the extrema*/
				for (j=0; j<xsize; j++){
					for(i=0; i<ysize; i++){
						c=*(array->ptr+j*ysize+i);
							r=c_re(c);
							max=(r > max) ? r : max;	/*concise if form*/
							min=(r < min) ? r : min;
					}
				}
				
				*high=max;								/*report results*/
				*low=min;			
			}
			break;

		case 2:		/*IMAGINARY SCALED */
			{
				/*find the extrema*/
				for (j=0; j<xsize; j++){
					for(i=0; i<ysize; i++){
						c=*(array->ptr+j*ysize+i);
							r=c_im(c);
							max=(r > max) ? r : max;	/*concise if form*/
							min=(r < min) ? r : min;
					}
				}
				
				*high=max;
				*low=min;		
			}
			break;

		}
	}
	
void edge(array_type *pic, array_type *psf, array_type *edg)
/***********************************************************
*	alex lehar 12 nov 01
*	take an input picture 'pic' and a blur point spread function
*	'psf', and create an edged picture, which combines the unaltered
*	'pic' at its center and a doubly blurred picture 'blr' formed 
*	by convolving 'pic' with 'psf'.  The transition from edge to 
*	center is controlled by some measure of the autocorrelation of the 
*	point spread function.  This is actually, the autocorrelation of y's
*	projection b1(x), and the autocorrelation of x's projection  b2(y). For this
*	routine it is necessary to generate temporary workspace internally
*	for 2-D 'blr' and the (2 x 3) 1-D autocorrelation projections. see journal p 24. 
*	note that we have to keep a spare copy of 'pic' on hand in an additional
*	buffer, 'pik'
*
*	edg		=		alpha D + (1-alpha) {pic*psf}   where
*
*		alpha(x,y)= (1-b1(x))*(1-b2(y))
*
************************************************************/
{
	array_type *blr,*buf,*pik,*b1,*b2,*b3,*b4,*b5,*b6;		/*internal buffers*/
	double sum, pix,q,r;
	int xpic,ypic,i,j;
	fftw_complex alpha,beta,g,h,s;
	
	/*find out size of the input array*/
	xpic=pic->xsize;
	ypic=pic->ysize;
	
	/*arrange internal buffers*/
	b1 = (array_type*)malloc(sizeof(array_type));		/*for y-projection*/
	b2 = (array_type*)malloc(sizeof(array_type));		/*for x-projection*/
	b3 = (array_type*)malloc(sizeof(array_type));		/*spare for y-projection*/
	b4 = (array_type*)malloc(sizeof(array_type));		/*spare for x-projection*/
	b5 = (array_type*)malloc(sizeof(array_type));		/*spare for y-projection*/
	b6 = (array_type*)malloc(sizeof(array_type));		/*spare for x-projection*/

	blr = (array_type*)malloc(sizeof(array_type));		/*for doubly blurred image*/
	buf = (array_type*)malloc(sizeof(array_type));		/*for internal image buffer*/
	pik = (array_type*)malloc(sizeof(array_type));		/*for internal image buffer*/

	
	/*main memory allocation*/
	b1->ptr = (fftw_complex *)malloc(xpic*sizeof(fftw_complex));	/*1-D array*/
	b2->ptr = (fftw_complex *)malloc(ypic*sizeof(fftw_complex));	/*		"  */
	b3->ptr = (fftw_complex *)malloc(xpic*sizeof(fftw_complex));	/*1-D array*/
	b4->ptr = (fftw_complex *)malloc(ypic*sizeof(fftw_complex));	/*		"  */
	b5->ptr = (fftw_complex *)malloc(xpic*sizeof(fftw_complex));	/*1-D array*/
	b6->ptr = (fftw_complex *)malloc(ypic*sizeof(fftw_complex));	/*		"  */


	blr->ptr = (fftw_complex *)malloc(xpic*ypic*sizeof(fftw_complex));
	buf->ptr = (fftw_complex *)malloc(xpic*ypic*sizeof(fftw_complex));
	pik->ptr = (fftw_complex *)malloc(xpic*ypic*sizeof(fftw_complex));
	
	/*set buffer sizes*/
	b1->xsize = xpic;
	b1->ysize = 1;
	b2->xsize = 1;
	b2->ysize = ypic;
	b3->xsize = xpic;
	b3->ysize = 1;
	b4->xsize = 1;
	b4->ysize = ypic;
	b5->xsize = xpic;
	b5->ysize = 1;
	b6->xsize = 1;
	b6->ysize = ypic;

	blr->xsize = xpic;
	blr->ysize = ypic;
	pik->xsize = xpic;
	pik->ysize = ypic;
	
	/*zero some buffers*/
	zero_array_complex(b1);				/*makes sure imaginaries are zero too*/
	zero_array_complex(b2);
	zero_array_complex(b3);
	zero_array_complex(b4);
	zero_array_complex(b5);
	zero_array_complex(b6);

	/*save the psf in buf*/
	move_array(psf, buf);				/*necessary because convolve overwrites psf*/
	move_array(pic, pik);				/*necessary because pic gets corrupted on convolve*/
		
	/*create doubly blurred picture*/
	cnvolv_spatial_to_spatial(pic, psf, blr, 0);							/*result is in blr*/
	stretch_complex_array(blr, 255.0, 0., 0);
	stretch_complex_array(pik, 255.0, 0., 0);			/*make sure both pictures properly scaled*/
	
	/*make the y-projection*/
	for (j=0; j<xpic; j++){
		sum=0.;											
		for (i=0; i<ypic; i++){						/*sum for projection element*/
			s=*(buf->ptr+j*ypic+i);
			pix=c_re(s);
			sum=sum+pix;
		}
		c_re(*(b1->ptr+j))=sum;
	}
		
	/*make the x-projection*/	
	for (i=0; i<ypic; i++){
		sum=0.;											
		for (j=0; j<xpic; j++){						/*sum for projection element*/
			s=*(buf->ptr+j*ypic+i);
			pix=c_re(s);
			sum=sum+pix;
		}
		c_re(*(b2->ptr+i))=sum;
	}
	
	
	/*make projection copies*/
	move_array(b1, b3);
	move_array(b2, b4);
	
	
	/*create autocorrelation of the projections*/
	cnvolv_spatial_to_spatial(b1, b3, b5, 0);
	cnvolv_spatial_to_spatial(b2, b4, b6, 0);

	/*swap the quadrants using trick of transform and shifting*/
	forward_fft(b5);
	forward_fft(b6);
	shift_array_complex(b5);
	shift_array_complex(b6);
	inverse_fft(b5);
	inverse_fft(b6);		


	/*scale the projections*/
	stretch_complex_array(b5, 1.0, 0., 0);
	stretch_complex_array(b6, 1.0, 0., 0);

	
	/*merge pic and blr arrays governed by autocorrelation projections*/
	for (j=0; j<xpic; j++){
		for (i=0; i<ypic; i++){
			s=*(b5->ptr+j);
			q=1.0-c_re(s);
			s=*(b6->ptr+i);
			r=1.0-c_re(s);
			c_re(alpha)=q*r;
			c_im(alpha)=0.;
			c_re(beta)=1.0-(q*r);
			c_im(beta)=0.;
			g=prod_c(alpha, *(pik->ptr+j*ypic+i));
			h=prod_c(beta, *(blr->ptr+j*ypic+i));
			*(edg->ptr+j*ypic+i)=add_c(g,h);
		}
	}
	
	free(b1->ptr);
	free(b2->ptr);											/*release the memory*/
	free(b3->ptr);
	free(b4->ptr);											
	free(b5->ptr);
	free(b6->ptr);											
	free(buf->ptr);
	free(blr->ptr);
	free(pik->ptr);
	free(b1);
	free(b2);											/*release the memory*/
	free(b3);
	free(b4);											
	free(b5);
	free(b6);											
	free(buf);
	free(blr);
	free(pik);
}	



void cepest(array_type *array, int xboxn, int yboxn, array_type *drray, int ovlp)
	/************************************************************
	*	alex lehar 9 nov 01
	*	This computes the Cepstrum (Transform of the log of the Transform of windowed
	*	chunks of the picture) for each of many sub-images, and then
	*	averages these together to get the result.   A lot of
	*	careful thresholding and clipping goes on, to try to distill the 
	*	contribution of 'zeros' structures in the power spectrum.
	*	array	Input Picture Complex Array Structure
	*	xboxn	number of non-overlapping boxes in vertical   x dir'n
	*	yxoxn	number of non-overlapping boxes in horizontal y dir'n
	*	drray	Contains the output, in same format as 'array', but since
	*			the output is smaller, it is embedded into the center
	*			of the larger format.  
	*	ovlp	=0 if non-overlapping, =1 if 50% overlap.  You pretty much 
	*			always want to overlap to get maximum information, particularly
	*			since the windowing destroys stuff at the edges.  
	*	based on average_power_spectrum.
	*	worked first 9 nov 01.
	*************************************************************/
	{
		array_type *brray, *crray;							/*necessary internal buffers*/
		/*indices are;i,j the big array &drray; e,f small arrays brray & crray;, m,n the box counters*/
		int i,j,xa,ya,xb,yb,xc,yc,m,n,e,f,total_boxes;
		fftw_complex c;
		double dc, hilev,lolev,range;
		
		/*unpack struct into local variables*/
		xa= array->xsize;
		ya= array->ysize;
		
		/*create the box sub-array 'brray' and accumulator buffer 'crray'*/
		/*allocate the array struct*/
		/*this allocates only the space for two ints and a pointer*/
		brray = (array_type*)malloc(sizeof(array_type));/*for the output*/
		crray = (array_type*)malloc(sizeof(array_type));/*for the accumulator*/
		
		brray->xsize = xa/xboxn;							/*determine box size from image*/
		brray->ysize = ya/yboxn;							/*               "             */
		
		crray->xsize = xa/xboxn;
		crray->ysize = ya/yboxn;

		/*confirm size of the box*/
		xb=brray->xsize;
		yb=brray->ysize;
		
		/*confirm size of the box*/
		xc=crray->xsize;
		yc=crray->ysize;
		
		/*allocate the memory for the box*/		
		brray->ptr = (fftw_complex *)malloc(xb*yb*sizeof(fftw_complex));
		crray->ptr = (fftw_complex *)malloc(xc*yc*sizeof(fftw_complex));

		/*zero the buffers*/
		zero_array_complex(crray);
		zero_array_complex(brray);
		/*set up the loop which will compute power spectrum in each block, load data
		into brray, and then accumulate them in crray, thus averaging them*/
		if (ovlp == 0)											/*non-overlapping case*/
			{
				/*loop through each box (tiles butting against eachother)*/
				for (m=0; m<yboxn; m++){
					for (n=0; n<xboxn; n++){
						/*transfer data from array to brray*/
						for (j= m*xb,f=0; j<(m+1)*xb; j++,f++){
							for (i= n*yb,e=0; i<(n+1)*yb; i++,e++){
								c=*(array->ptr+j*ya+i);				/*read pixel in big array*/
								*(brray->ptr+yb*f+e)=c;				/*write pixel to small box*/
							}
						}
						hanning_window_complex(brray);				/*do the edging*/
						shift_array_complex(brray);
						forward_fft(brray); 						/*FFT the Box*/
						magnitude_phase_complex(brray);				/*polar form*/
						convert_to_logthreshold_complex(brray, 1.0e-34);/*log and clip*/
						invert_realpart_complex(brray);				/*inversion*/
						maxmin_array(brray, &hilev, &lolev, 0);		/*find extrema*/
						range=hilev-lolev;							/*range of values*/
						clip_realpart_complex(brray, lolev+(.45*range),
							hilev-(.2*range));						/*CRITICAL clipping*/
						stretch_complex_array(brray, 255., 0., 0);	/*scale to normal range*/
						shift_array_complex(brray);
						inverse_fft(brray);
						magnitude_phase_complex(brray);
						convert_to_logthreshold_complex(brray, .001);
						accumulate_array_complex(brray,crray);		/*save the sum in buffer*/
					}
				}
			}
		else
			{
				/*loop through more overlapping boxes*/
				for (m=0; m<(2*yboxn-1); m++){						/*because there are 2*m-1 of them*/
					for (n=0; n<(2*xboxn-1); n++){					/*		"			2*n-1	"	 */
						/*transfer data from array to brray*/
						for (j= m*xb/2,f=0; f<xb; j++,f++){
							for (i= n*yb/2,e=0; e<yb; i++,e++){
								c=*(array->ptr+j*ya+i);				/*read pixel in big array*/
								*(brray->ptr+yb*f+e)=c;				/*write pixel to small box*/
							}
						}
						hanning_window_complex(brray);				/*do the edging*/
						shift_array_complex(brray);
						forward_fft(brray); 						/*FFT the Box*/
						magnitude_phase_complex(brray);				/*polar form*/
						convert_to_logthreshold_complex(brray, 1.0e-34);/*log and clip*/
						invert_realpart_complex(brray);				/*inversion*/
						maxmin_array(brray, &hilev, &lolev, 0);		/*find extrema*/
						range=hilev-lolev;							/*range of values*/
						clip_realpart_complex(brray, lolev+(.45*range),
							hilev-(.2*range));						/*CRITICAL clipping*/
						stretch_complex_array(brray, 255., 0., 0);	/*scale to normal range*/
						shift_array_complex(brray);
						inverse_fft(brray);
						magnitude_phase_complex(brray);
						convert_to_logthreshold_complex(brray, .001);
							
						accumulate_array_complex(brray,crray);		/*save the sum in buffer*/
					}
				}
			}
		
		/*record the total number of boxes processed*/
		total_boxes=m*n;
		divide_element(crray, (float)(total_boxes));			/*effectuate uggh! the average*/
		
		/*having completed the accumulation process, now embed the result into the input image
		format of the output array*/
		value_array_complex(drray, 0.,0.);					/*first zero the target array*/

		dc=set_dc_complex(crray, 0., 0);					/*dc is at center*/
		stretch_complex_array(crray, 255., 0., 0);

		for (j=(xa/2)-(xb/2),f=0; f<xb;j++,f++){
			for(i=(ya/2)-(yb/2),e=0; e<yb; i++,e++){
				c=*(drray->ptr+j*ya+i);
				*(drray->ptr+j*ya+i)=*(crray->ptr+f*yc+e);
			}
		}
	free(brray->ptr);
	free(crray->ptr);											/*release temporary memory allocation*/		
	free(brray);
	free(crray);											/*release temporary memory allocation*/		
		
	}
	
	void move_array(array_type *array, array_type *brray)
/*************************************************************
*	alex lehar	9 nov 01
*	copy the data in 'array' to 'brray'
**************************************************************/
{
	int i,j,xsize, ysize;
	fftw_complex c;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;
	
	/*set the values*/
	for (j=0;j<xsize;j++){
		for(i=0;i<ysize; i++){
			c=*(array->ptr+j*ysize+i);
			*(brray->ptr+j*ysize+i)=c;
		}
	}
}	


void patgen(array_type *array, array_type *crray, int xboxn, int yboxn, int ovlp)
/************************************************************
*	alex lehar 15 nov 01
*	create a spatial domain filter to remove any regular
*	pattern from input image array, and put the filter into crray.
*	This routine will examine the Fourier Domain magnitude
*	and then find the energy 'peaks' by differencing with the 
*	prevailing annular average value. (Accumulated in array
*	aada) This signal will then
*	be appropriately inverted and clipped, and then inverse
*	transformed to produce filter kernel.
*	input array is not corrupted
*	arguments:
*	array			input array (remains uncorrupted)
*	crray			output spatial domain kernel, same size, 
*					zero-filled as necessary
*	xboxn			number of horizontal tiles
*	yboxn			number of vertical tiles
*	ovlp			0 if nonoverlap, 1 if overlapped for power 
*					spectrum computation.
*************************************************************/
{
	array_type  *aada;						/*required internal buffer*/
	int i,j,x,y,hi,hj,n,m,jx,iy, bigimage_area, smallimage_area,hjx,hiy,ii,jj;
	float fi,fj;
	double p,q,r,sum,hilev,lolev,range,dc,xmax,chi, area_ratio;
	fftw_complex c;
	
	/*unpack struct into local variables*/
	x= array->xsize;
	y= array->ysize;
	hi=y/2;
	hj=x/2;
	fi=((float)y)/2.;
	fj=((float)x)/2.;								/*center point*/
	
	/*allocate memory for aada*/
	aada = (array_type*)malloc(sizeof(array_type));
	/*annular average distribution array*/
	n=(int)(sqrt(hi*hi+hj*hj));		/*dimension of half diagonal to image, 
										add arbitrary 4 to be safe*/	
	aada->xsize=n;
	aada->ysize=1;
	aada->ptr = (fftw_complex *)malloc((n+4)*sizeof(fftw_complex));	/*+4 is for safety elbow room*/
	
	/*compute the average power spectrum*/
	average_power_spectrum(array, xboxn, yboxn, crray, ovlp, &jx, &iy);
	hjx=(jx)/2;
	hiy=(iy)/2;						/*corresponding half-sizes*/
	
	/*note that size of the box is smaller than main image*/
	bigimage_area=x*y;
	smallimage_area=jx*iy;
	area_ratio=((double)bigimage_area)/((double)smallimage_area);
		
	/*Go through the Power Spectrum and build up the Annular Average Distribution
	Array 'aada'  Note devilish clever use of real and imaginary parts of the aada. 
	Modification 23 nov 01: but don't examine beyond the limits of the central box!*/
		zero_array_complex(aada);		/*initialize*/
		 for (j=hj-hjx,jj=0; jj<jx; jj++,j++) {				/*within the box*/
		 	for (i=hi-hiy,ii=0; ii<iy; ii++,i++) {				/*within the box*/
		 		c=*(crray->ptr+j*y+i);	/*read the (magnitude) element*/
		 		p=c_re(c);				/*it's in the real part*/
		 		r=(double)sqrt((float)(i-fi)*(i-fi)+(j-fj)*(j-fj));/*radius from center*/
		 		m=rint(r);				/*decide which bin to put into*/
		 		q=c_re(*(aada->ptr+m));	/*what's already in the bin*/
		 		q=q+p;					/*accumulate in the right bin*/
		 		c_re(*(aada->ptr+m))=q;
		 		c_im(*(aada->ptr+m))=c_im(*(aada->ptr+m))+1.0;/*craftily use complex part to 
		 								record how many increments are in the bin for use in
		 								later averaging calculation*/
		 	}
		 }
	/*now go through aada to compute the annular averages*/
	for (m=0; m<n; m++) {
		p=c_re(*(aada->ptr+m));
		q=c_im(*(aada->ptr+m));			/*there better not be any zero q's!!*/
		c_re(*(aada->ptr+m))=p/q;		/*compute the averages*/
	}
		
	/*now go through the Power Spectrum and subtract the radial average to expose
	any peaks.  Modify to just look in the box area*/
		 for (j=hj-hjx,jj=0; jj<jx; jj++,j++) {				/*within the box*/
		 	for (i=hi-hiy,ii=0; ii<iy; ii++,i++) {				/*within the box*/
		 		c=*(crray->ptr+j*y+i);	/*read the (magnitude) element*/
		 		p=c_re(c);				/*it's in the real part*/
		 		r=(double)sqrt((float)((i-fi)*(i-fi)+(j-fj)*(j-fj)));/*radius from center*/
		 		m=rint(r);				/*decide which bin to refer to*/
		 		q=c_re(*(aada->ptr+m));	/*what's already in the bin*/
		 		p=p-q;					/*subtract the annular average*/
		 		c_re(*(crray->ptr+j*y+i))=p;/*assign the new value*/
		 	}
		 }
		 

	/*Now take the annular-average-difference signal and condition in various ways
	to create a fourier domain pattern attenuation filter*/
	mult_const_complex(crray, area_ratio);/*scaling compensation*/
	dc=get_dc_complex(crray,0);			/*dc is at center*/
	maxmin_array_special(crray, &hilev, &lolev, 0, &jx, &iy);	/*find extrema in central box*/
	range=hilev-lolev;					/*calculate range of values*/
	add_array_complex(crray, -dc);		/*subtract dc to 'normalize' kernel*/
	clip_realpart_complex(crray, 0.2*hilev, hilev);/*clip*/
	add_array_complex(crray, -.2*hilev);/*bias dc to zero*/
	xmax=.8*hilev;						/*new max*/
xi:	chi=35./xmax;						/*scale factor*/
	mult_const_complex(crray, chi);		/*scaling for 35db gain*/
	invert_realpart_complex(crray);		/*perform inversion*/
	convert_to_exponential_complex(crray);/*exponentiate.  Values here range from
										max of 1.0 to min of near zero.  The DC is at
										the center.  The filter waits in crray*/

	shift_array_complex(crray);			/*quadrant shifting*/
	set_dc_complex(crray, 0., 0);		/*set dc to zero*/
	inverse_fft(crray);					/*back to spatial domain*/

		
	free(aada->ptr);							/*release the local memory*/
	free(aada);							/*release the local memory*/
	}
	
	int cnvolv_spatial_to_spatial(array_type *array, array_type *brray, array_type *crray, int window)
	/************************************************************
	*	alex lehar 17 nov 01
	*	take the input double precision complex input 'array'
	*	and similar 'brray', each containing a spatial domain data
	*	in the Real part, and zero in Imaginary; and convolve them
	*	to produce an output 'crray', containing the spatial domain
	*	result in its real part.  'brray' is normalized first.
	*	All this is accomplished by
	*	multiplication in the intermediate Fourier domain according
	*	to the Convolution Theorem. 
	*	If Window=0, use no windowing, if 1, use Hanning Window. 
	*	note that 'array' is corrupted in the process!
	*	revised from 'convolve', to clean up use of fft.
	*	WARNING  seems not to work for odd sized images. investigate!
	*************************************************************/
	{
		int xa,ya,xb,yb,xc,yc;
		
		/*unpack struct into local variables*/
		xa= array->xsize;
		ya= array->ysize;
		xb=brray->xsize;
		yb=brray->ysize;
		xc=crray->xsize;
		yc=crray->ysize;
		
		if (((xa+ya) + (xb+yb) + (xc+yc) ) != 3*(xa+ya))	/*incompatible sizes*/
			return (0);	/*escape if images incompatible*/
		else
			{
			if(window == 0)
				{
				forward_fft(array); 				/*Transform*/
				forward_fft(brray); 				/*Transform*/
				//normalize_array_complex(brray);
				//normalize_justreals(brray);			/*normalization*/
				mult_array_complex(array, brray, crray);	/*multiply=convolution*/
				shift_array_complex(crray);
				inverse_fft(crray);
				//scale_array_complex(crray);			/*scale back*/
				}
			else if (window ==1)
				{
				hanning_window_complex(array);		/*apply window*/
				forward_fft(array); 				/*Transform*/
				hanning_window_complex(brray);		/*apply window*/
				forward_fft(brray); 				/*Transform*/
				//normalize_array_complex(brray);
				//normalize_justreals(brray);			/*normalization*/
				mult_array_complex(array, brray, crray);	/*multiply=convolution*/
				shift_array_complex(crray);
				inverse_fft(crray);
				//scale_array_complex(crray);			/*scale back*/
				}
			}
	}

void convert_to_exponential_complex(array_type *array)
/*************************************************************
*	alex lehar	18 nov 01
*	take the input dynamically resizable input complex pair array
*	convert all reals to their exponentials
**************************************************************/
{
	int i,j,xsize, ysize;
	double x,y,p,q;
	fftw_complex c,d;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;
	
	/*set the values*/
	for (j=0;j<xsize;j++){
		for(i=0;i<ysize; i++){
			c=*(array->ptr+j*ysize+i);
			x=c_re(c);
			y=0.;
			p=exp(x);
			c_re(d)=p;
			c_im(d)=y;
			*(array->ptr+j*ysize+i)=d;
		}
	}
}	

void add_array_complex(array_type *array, double r)
/*************************************************************
*	Alex Lehar	19 nov 01
*	take the input dynamically resizable input complex pair array
*	and (complex) add constant r to each realpart element
**************************************************************/

{
	int i,j,xsize, ysize;
	fftw_complex a;
	double c,d;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;

	
	/*calculate the values*/
	for (j=0;j<xsize;j++){
		for(i=0;i<ysize; i++){
			a=*(array->ptr+j*ysize+i);
			c_re(a)=c_re(a)+r;
			*(array->ptr+j*ysize+i)=a;
		}
	}
}	

void mult_const_complex(array_type *array, double r)
/*************************************************************
*	Alex Lehar	19 nov 01
*	take the input dynamically resizable input complex pair array
*	and multiply each realpart element by r.
**************************************************************/

{
	int i,j,xsize, ysize;
	fftw_complex a;
	double c,d;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;

	
	/*calculate the values*/
	for (j=0;j<xsize;j++){
		for(i=0;i<ysize; i++){
			a=*(array->ptr+j*ysize+i);
			c_re(a)=c_re(a)*r;
			*(array->ptr+j*ysize+i)=a;
		}
	}
}	
	int deblur1(array_type *array, array_type *brray, array_type *crray, float snr)
	/************************************************************
	*	alex lehar 30 oct 01
	*	variant of deblur, to use windowing on the map kernel.
	*	array	Input Picture Complex Array Structure
	*	brray	Input PSF Complex Array Structure
	*			Output contains Restored Picture
	*	crray	Output Linear MAP Kernel.
	*	snr		Estimated Signal to Noise Ratio (100 typical)
	*	implements Fourier Domain Restoration, using the form of the 
	*	Convolution operation Spatial-Fourier-BacktoSpatial
	*	WARNING.  You might have to use windowing option (=1) in 
	*	final convolution to control ringing in the picture.  
	*************************************************************/
	{
		int xa,ya,xb,yb,xc,yc;
		
		/*unpack struct into local variables*/
		xa= array->xsize;
		ya= array->ysize;
		xb=brray->xsize;
		yb=brray->ysize;
		xc=crray->xsize;
		yc=crray->ysize;

		
		if (((xa+ya) + (xb+yb) + (xc+yc) ) != 3*(xa+ya))	/*incompatible sizes*/
			return (0);										/*escape if images incompatible*/
		else
			{
			mapgen_complex(brray, crray, snr);				/*MAP kernel in crray*/
			convolve(array, crray, brray, 1);				/*spatial domain convolution with windowing*/
															/*result sits in brray overwriting psf*/
			}
	}


double getsum_array_complex(array_type *array, int iwhich)
/*************************************************************
*	alex lehar 21 nov 01
*	take the input dynamically resizable input complex pair array
*	and find out what the sum of the components are in either
*	the real iwhich=0, or imaginary iwhich >= 1 parts.
**************************************************************/
{
	int i,j,xsize,ysize;
	fftw_complex c;
	double sumr,sumi,r,s;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;
	
	/*initialize the sums*/
	sumr=0.;
	sumi=0.;
		
	/*go get'em*/
	for (j=0;j<xsize;j++){
		for(i=0;i<ysize; i++){
			c=*(array->ptr+j*ysize+i);
			r=c_re(c);
			s=c_im(c);
			if (iwhich == 0)
				{
				sumr=sumr+r;
				}
			else
				{
				sumi=sumi+s;
				}
		}
	}
	if (iwhich == 0) return sumr;
	else return sumi;	
}	

void normalize_justreals(array_type *array)
/*************************************************************
*	alex lehar 21 nov 01
*	take the input dynamically resizable input complex pair array
*	and normalize just the real part, by 1/(sum of elements)
*	the complex value is also normalized to the same realsum.
**************************************************************/
{
	int i,j,xsize,ysize;
	fftw_complex c;
	double sum,r,s;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;
	
	/*initialize the sum*/
	sum=0.;
		
	/*find the sum for scale*/
	for (j=0;j<xsize;j++){
		for(i=0;i<ysize; i++){
			c=*(array->ptr+j*ysize+i);
			r=c_re(c);
			sum=sum+r;		/*accumulate reals*/
		}
	}
	
	/*do the normalization*/
	for (j=0;j<xsize;j++){
		for(i=0;i<ysize; i++){
			c=*(array->ptr+j*ysize+i);
			r=c_re(c);
			s=c_im(c);
			c_re(c)=r/sum;
			c_im(c)=s/sum;	/*just for the heck of it*/
			*(array->ptr+j*ysize+i)=c;
		}
	}
	
}	
void maxmin_array_special(array_type *array, double *high, double *low, int r_or_i, 
	int *jx, int *iy)
/**********************************************************************************
*	alex lehar	23 nov 01
*	special modification on maxmin_array, to introduce a couple of extra arguments
*	to enable the selective examination of central box region defined by jx, iy.
*
*	Input an input dynamically resizable input complex pair array, pointed to by
*	*array,  and choose whether to deal with the real or imaginary part, for 
*	r_or_i= 0 or 1 respectively.  Then determine the max (high) and min (low)
*	values, and return these values to the addresses of 'high' and 'low'
*	method of calling is: 
*			maxmin_array(array, &high, &low, &jx, &iy);
***********************************************************************************/
{
	int i,j,xsize,ysize,hx,hy,ii,jj,hjx,hiy;
	double r,min,max;
	fftw_complex c;
	min=1e34;
	max=-1e34;	/*starting values*/
		
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;
	hx=xsize/2;					/*half point parameters*/
	hy=ysize/2;
	hjx=*jx/2;
	hiy=*iy/2;
	
	/*examine one of two possibilities, realorimag = 0 or 1  */
	switch (r_or_i) 
	{
		case 0:		/*REAL SCALED*/
			{
				/*find the extrema*/
				for (j=hx-hjx, jj=0; jj<*jx; jj++,j++){		/*look just within box*/
					for(i=hy-hiy, ii=0; ii<*iy; ii++,i++){	/*look just within box*/
						c=*(array->ptr+j*ysize+i);
							r=c_re(c);
							max=(r > max) ? r : max;	/*concise if form*/
							min=(r < min) ? r : min;
					}
				}
				
				*high=max;								/*report results*/
				*low=min;			
			}
			break;

		case 2:		/*IMAGINARY SCALED */
			{
				/*find the extrema*/
				for (j=hx-hjx, jj=0; jj<*jx; jj++,j++){		/*look just within box*/
					for(i=hy-hiy, ii=0; ii<*iy; ii++,i++){	/*look just within box*/
						c=*(array->ptr+j*ysize+i);
							r=c_im(c);
							max=(r > max) ? r : max;	/*concise if form*/
							min=(r < min) ? r : min;
					}
				}
				
				*high=max;
				*low=min;		
			}
			break;

		}
	}
	
double get_dc_complex(array_type *array, int icen)
/*************************************************************
*	Alex Lehar	23 nov 01
*	just return the dc value
*	*array		complex data structure
*	dc			value to set dc at
*	icen		=0 if at array center, =1 if at topleft corner
**************************************************************/
{
	int i,j,xsize, ysize;
	int k,l;						/*indices for 'found' dc term*/
	fftw_complex c;
	double max,r,s,t;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;

	
	/*set dc's indices and value*/
	if (icen == 0)
		{
		k=ysize/2;
		l=xsize/2;
		}
	else
		{
		k=0;
		l=0;
		}
	c=*(array->ptr+l*ysize+k);
	max=c_re(c);	
	return max;
}	

void pattern_removal(array_type *array, array_type *crray, array_type *brray, 
	int xboxn, int yboxn, int ovlp, float spike)
/************************************************************
*	alex lehar 25 nov 01
*	create a fourier domain filter to remove any regular
*	pattern from input image array, and put the transform-domain 
*	filter into crray, and resulting picture into brray.
*	This routine will examine the Fourier Domain magnitude
*	and then find the energy 'peaks' by differencing with the 
*	prevailing annular average value. (Accumulated in array
*	aada) This signal will then
*	be appropriately inverted and clipped, and then inverse
*	transformed to produce filter kernel. 'spike' is the fraction
*	of the spike which is effectively cut out.
*	input array is transformed.
*	arguments:
*	array			input array (gets transformed)
*	crray			output fourier domain kernel, same size, 
*					zero-filled as necessary
*	brray			output filtered image
*	xboxn			number of horizontal tiles
*	yboxn			number of vertical tiles
*	ovlp			0 if nonoverlap, 1 if overlapped for power 
*					spectrum computation.
*	spike			spike attenuation usually set at 0.8 added 17 dec.
*************************************************************/
{
	array_type  *aada;						/*required internal buffer*/
	int i,j,x,y,hi,hj,n,m,jx,iy, bigimage_area, smallimage_area,hjx,hiy,ii,jj;
	float fi,fj;
	double p,q,r,sum,hilev,lolev,range,dc,xmax,chi, area_ratio;
	fftw_complex c;
	
	/*unpack struct into local variables*/
	x= array->xsize;
	y= array->ysize;
	hi=y/2;
	hj=x/2;
	fi=((float)y)/2.;
	fj=((float)x)/2.;								/*center point*/
	
	/*allocate memory for aada*/
	aada = (array_type*)malloc(sizeof(array_type));
	/*annular average distribution array*/
	n=(int)(sqrt(hi*hi+hj*hj));		/*dimension of half diagonal to image, 
										add arbitrary 4 to be safe*/	
	aada->xsize=n;
	aada->ysize=1;
	aada->ptr = (fftw_complex *)malloc((n+4)*sizeof(fftw_complex));	/*+4 is for safety elbow room*/
	
	/*compute the average power spectrum*/
	average_power_spectrum(array, xboxn, yboxn, crray, ovlp, &jx, &iy);
	/*apply hanning window*/
	//hanning_window_box(crray, jx, iy);
	hjx=(jx)/2;
	hiy=(iy)/2;						/*corresponding half-sizes*/
	
	/*note that size of the box is smaller than main image*/
	bigimage_area=x*y;
	smallimage_area=jx*iy;
	area_ratio=((double)bigimage_area)/((double)smallimage_area);
		
	/*Go through the Power Spectrum and build up the Annular Average Distribution
	Array 'aada'  Note devilish clever use of real and imaginary parts of the aada. 
	Modification 23 nov 01: but don't examine beyond the limits of the central box!*/
		zero_array_complex(aada);		/*initialize*/
		 for (j=hj-hjx,jj=0; jj<jx; jj++,j++) {				/*within the box*/
		 	for (i=hi-hiy,ii=0; ii<iy; ii++,i++) {				/*within the box*/
		 		c=*(crray->ptr+j*y+i);	/*read the (magnitude) element*/
		 		p=c_re(c);				/*it's in the real part*/
		 		r=(double)sqrt((float)(i-fi)*(i-fi)+(j-fj)*(j-fj));/*radius from center*/
		 		m=rint(r);				/*decide which bin to put into*/
		 		q=c_re(*(aada->ptr+m));	/*what's already in the bin*/
		 		q=q+p;					/*accumulate in the right bin*/
		 		c_re(*(aada->ptr+m))=q;
		 		c_im(*(aada->ptr+m))=c_im(*(aada->ptr+m))+1.0;/*craftily use complex part to 
		 								record how many increments are in the bin for use in
		 								later averaging calculation*/
		 	}
		 }
	/*now go through aada to compute the annular averages*/
	for (m=0; m<n; m++) {
		p=c_re(*(aada->ptr+m));
		q=c_im(*(aada->ptr+m));			/*there better not be any zero q's!!*/
		c_re(*(aada->ptr+m))=p/q;		/*compute the averages*/
	}
		
	/*now go through the Power Spectrum and subtract the radial average to expose
	any peaks.  Modify to just look in the box area*/
		 for (j=hj-hjx,jj=0; jj<jx; jj++,j++) {				/*within the box*/
		 	for (i=hi-hiy,ii=0; ii<iy; ii++,i++) {				/*within the box*/
		 		c=*(crray->ptr+j*y+i);	/*read the (magnitude) element*/
		 		p=c_re(c);				/*it's in the real part*/
		 		r=(double)sqrt((float)((i-fi)*(i-fi)+(j-fj)*(j-fj)));/*radius from center*/
		 		m=rint(r);				/*decide which bin to refer to*/
		 		q=c_re(*(aada->ptr+m));	/*what's already in the bin*/
		 		p=p-q;					/*subtract the annular average*/
		 		c_re(*(crray->ptr+j*y+i))=p;/*assign the new value*/
		 	}
		 }
		 

	/*Now take the annular-average-difference signal and condition in various ways
	to create a fourier domain pattern attenuation filter*/
	mult_const_complex(crray, area_ratio);/*scaling compensation*/
	dc=get_dc_complex(crray,0);			/*dc is at center*/
	maxmin_array_special(crray, &hilev, &lolev, 0, &jx, &iy);	/*find extrema in central box*/
	range=hilev-lolev;					/*calculate range of values*/
	add_array_complex(crray, -dc);		/*subtract dc to 'normalize' kernel*/
										/*spike proportion used to create filter, .8 usually*/
	clip_realpart_complex(crray, (1.0-spike)*hilev, hilev);/*clip*/
	add_array_complex(crray, -(1.0-spike)*hilev);/*bias dc to zero*/
	xmax=spike*hilev;					/*new max*/
xx:	chi=35./xmax;						/*scale factor*/
	mult_const_complex(crray, chi);		/*scaling for 35db gain*/
	invert_realpart_complex(crray);		/*perform inversion*/
	convert_to_exponential_complex(crray);/*exponentiate.  Values here range from
										max of 1.0 to min of near zero.  The DC is at
										the center.  The filter waits in crray*/
	maxmin_array_special(crray, &hilev, &lolev, 0, &jx, &iy);	/*find extrema in central box*/
	border_padding(crray,jx,iy,hilev);	/*important because edge frequencies don't exist  
										now in the matching crray, set to 1. */
	shift_array_complex(array);			/*prepare for fft*/
	forward_fft(array);					/*transform of array is now center-dc'd*/
	border_padding(array, jx, iy, 0.);	/*set borders to 0.*/

	mult_array_complex(array, crray, brray);/*do the filtering*/
	quadrant_swap(brray);
	//hanning_window_complex(brray);	/*window the whole image*/
	inverse_fft(brray);	
	free(aada->ptr);					/*release the local memory*/
	free(aada);							/*release the local memory*/
	}
	
void remove_nyquist(array_type *array, int jx, int iy)
/*************************************************************
*	alex lehar 25 nov 01
*	take the input dynamically resizable input complex pair array
*	and remove nyquist terms in the fourier domain. Assume that 
*	the transform is dc centered. Nyquist all those terms beyond 
*	the box defined in the center, dimensions jx, iy.
*
*	modified 27 nov 01... to pad the area outside the box with a
*	gray value, corresponding to the average of the values within
*	the box.  
**************************************************************/
{
	int i,j,xsize,ysize,ii,jj, hj, hi, hjx, hiy, gapx, gapy;
	fftw_complex c;
	double sum,n,average;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;
	
	hj=xsize/2;		/*image half size*/
	hi=ysize/2;
	hjx=jx/2;		/*box half size*/
	hiy=iy/2;
	gapx=(double)((xsize-jx)/2); /*strip width*/
	gapy=(double)((ysize-iy)/2);
	
#if 0
	/*find the average value in the box*/
	sum=0.;
	n=0.;
	for (j=gapx,jj=0; jj<jx; j++, jj++){
		for (i=gapy, ii=0; ii<iy; i++, ii++){
			n=n+1.;
			c=*(array->ptr+j*ysize+i);
			sum=sum+c_re(c);
		}
	}
	average=sum/n;
#endif	
	average=1.0;		/*override for pattern removal application*/	
	/*set constants*/
	c_re(c)=average;
	c_im(c)=0.;
		
	
	/*set borders to zero*/
		/*do the top swath*/
		for (j=0, jj=0; jj<gapx; j++,jj++)
			{
				for(i=0;i<ysize; i++){
					*(array->ptr+j*ysize+i)=c;
				}
			}
		/*do the bottom swath*/
		for (j=hj+hjx,jj=0; jj<(gapx); j++,jj++)
			{
				for(i=0;i<ysize; i++){
					*(array->ptr+j*ysize+i)=c;
				}
			}

		/*do the left swath*/
		for (i=0, ii=0; ii<gapy; i++,ii++)
			{
				for(j=0;j<xsize; j++){
					*(array->ptr+j*ysize+i)=c;
				}
			}
		/*do the right swath*/
		for (i=hi+hiy,ii=0; ii<(gapy); i++, ii++)
			{
				for(j=0;j<xsize; j++){
					*(array->ptr+j*ysize+i)=c;
				}
			}

}	

void mult_array_complex_special(array_type *array, array_type *brray, 
	int jx, int iy, array_type *crray)
/*************************************************************
*	Alex Lehar	27 nov 01
*	take the input dynamically resizable input complex pair array
*	and (complex) multiply each element with that of brray, to produce
*	crray.  In this routine we allow the combination of unequally
*	sized regions of array and brray.  brray has a center which 
*	has box dimensions jx, iy, and we only modify array pixels
*	within this box span.
*	arguments:
*	*array			input array (big)
*	*brray			input array	(small)
*	jx				internal x box size of brray
*	iy				internal y box size of brray
**************************************************************/

{
	int i,j,x,y,hx,hy,hjx,hiy,xgap,ygap;
	fftw_complex a,b,c;
	
	/*unpack struct into local variables*/
	x= array->xsize;
	y= array->ysize;
	hx=x/2;
	hy=y/2;
	hjx=jx/2;
	hiy=iy/2;
	xgap=(x-jx)/2;
	ygap=(y-iy)/2;

	/*calculate the values differently in and out of box*/
	for (j=0;j<x;j++){
		for(i=0;i<y; i++){
			if (((j > xgap) || j > (xgap+jx)) && ((i > ygap) || i > (ygap+iy))) /*in box*/
				{
				a=*(array->ptr+j*y+i);
				b=*(brray->ptr+j*y+i);
				c=prod_c(a,b);
				*(crray->ptr+j*y+i)=c;
				}
			else															/*out of box*/
				{
				a=*(array->ptr+j*y+i);
				*(crray->ptr+j*y+i)=a;
				}
		}
	}
}	

void mult_array_otherreal(array_type *array, array_type *brray, array_type *crray)
/*************************************************************
*	Alex Lehar	25 oct 01
*	take the input dynamically resizable input complex pair array
*	and multiply real and imaginary parts of 'array' with the 
*	real part of 'brray'
**************************************************************/

{
	int i,j,xsize, ysize;
	fftw_complex a,b,c;
	double bb;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;

	
	/*calculate the values*/
	for (j=0;j<xsize;j++){
		for(i=0;i<ysize; i++){
			a=*(array->ptr+j*ysize+i);
			b=*(brray->ptr+j*ysize+i);
			bb=c_re(b);
			c_re(c)=c_re(a)*bb;
			c_im(c)=c_im(a)*bb;
			*(crray->ptr+j*ysize+i)=c;
		}
	}
}	

void hanning_window_box(array_type *array, int jx, int iy)
/*************************************************************
*	alex lehar 28 nov 01
*	take the input dynamically resizable input complex pair array
*	and apply a hanning window attenuation to the central box
*	defined by jx & iy
**************************************************************/
{
	int i,j,x,y,ii,jj, hj, hi, hjx, hiy, gapx, gapy;
	float r,d,a,b,weight;
	fftw_complex c;
	#define PI 3.1415927
	
	/*unpack struct into local variables*/
	x= array->xsize;
	y= array->ysize;
	
	hj=x/2;		/*image half size*/
	hi=y/2;
	hjx=jx/2;		/*box half size*/
	hiy=iy/2;
	gapx=(double)((x-jx)/2); /*strip width*/
	gapy=(double)((y-iy)/2);
	

	/*go through the array*/
	/*radius of the circumscribed circle*/
	r=sqrt(x*x+y*y)/2.0;
	for (j=gapx,jj=0; jj<jx; j++, jj++){
		for (i=gapy, ii=0; ii<iy; i++, ii++){
			/*distance between center and examined point*/
			a=pow((x/2-j),2);
			b=pow((y/2-i),2);
			d=sqrt(a+b);
			/*cosine weighting factor*/
			weight=1.0-cos((PI*(r-d)/2.0)/r);
			c=*(array->ptr+j*y+i);
			c_re(c)=c_re(c)*weight;
			c_im(c)=c_im(c)*weight;	/*weight both real and imaginary parts*/
			*(array->ptr+j*y+i)=c;	/*final reassignment*/
		}
	}

}	

void border_padding(array_type *array, int jx, int iy, double value)
/*************************************************************
*	alex lehar 25 nov 01
*	take the input dynamically resizable input complex pair array
*	and  to pad the area outside the box with a
*	gray 'value'. If the value is negative, then 'value' is
*	determined by the average value of the pixels inside the box
*	defined by jx, iy.  If positive, then determined by 
*	value itself. Border is only set in both real and imaginary 
*	parts.
**************************************************************/
{
	int i,j,xsize,ysize,ii,jj, hj, hi, hjx, hiy, gapx, gapy;
	fftw_complex c;
	double sum,n,average;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;
	
	hj=xsize/2;		/*image half size*/
	hi=ysize/2;
	hjx=jx/2;		/*box half size*/
	hiy=iy/2;
	gapx=(double)((xsize-jx)/2); /*strip width*/
	gapy=(double)((ysize-iy)/2);
	
if (value < 0.)	/*automatic determination*/
	{
		/*find the average value in the box*/
		sum=0.;
		n=0.;
		for (j=gapx,jj=0; jj<jx; j++, jj++){
			for (i=gapy, ii=0; ii<iy; i++, ii++){
				n=n+1.;
				c=*(array->ptr+j*ysize+i);
				sum=sum+c_re(c);
			}
		}
		average=sum/n;
	}
else
	{	
		average=value;		/*specified value*/	
	}
	
	/*set constants*/
	c_re(c)=average;
	c_im(c)=average;
		
	
	/*set borders to zero*/
		/*do the top swath*/
		for (j=0, jj=0; jj<gapx; j++,jj++)
			{
				for(i=0;i<ysize; i++){
					*(array->ptr+j*ysize+i)=c;
				}
			}
		/*do the bottom swath*/
		for (j=hj+hjx,jj=0; jj<(gapx); j++,jj++)
			{
				for(i=0;i<ysize; i++){
					*(array->ptr+j*ysize+i)=c;
				}
			}

		/*do the left swath*/
		for (i=0, ii=0; ii<gapy; i++,ii++)
			{
				for(j=0;j<xsize; j++){
					*(array->ptr+j*ysize+i)=c;
				}
			}
		/*do the right swath*/
		for (i=hi+hiy,ii=0; ii<(gapy); i++, ii++)
			{
				for(j=0;j<xsize; j++){
					*(array->ptr+j*ysize+i)=c;
				}
			}

}	

void quadrant_swap(array_type *array)
/*************************************************************
*	Alex Lehar	30 nov 01
*	take the input dynamically resizable input complex pair array
*	and swap the quadrants.  Do so using the shift property of 
*	the Fourier Transform.  Only works for even sized array.
**************************************************************/

{
	forward_fft(array);
	shift_array_complex(array);
	inverse_fft(array);
}	

void punch_box(array_type *array, int jx, int iy, double value)
/*************************************************************
*	alex lehar 2 dec 01
*	take the input dynamically resizable input complex pair array
*	and  overwrite a box with a
*	gray 'value'. If the value is negative, then 'value' is
*	determined by the average value of the pixels inside the centrally
*	located box defined by jx, iy.  If positive, then determined by 
*	value itself. Box area is set in both real and imaginary 
*	parts.
**************************************************************/
{
	int i,j,xsize,ysize,ii,jj, hj, hi, hjx, hiy, gapx, gapy;
	fftw_complex c;
	double sum,n,average;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;
	
	hj=xsize/2;		/*image half size*/
	hi=ysize/2;
	hjx=jx/2;		/*box half size*/
	hiy=iy/2;
	gapx=(double)((xsize-jx)/2); /*strip width*/
	gapy=(double)((ysize-iy)/2);
	
if (value < 0.)	/*automatic determination*/
	{
		/*find the average value in the box*/
		sum=0.;
		n=0.;
		for (j=gapx,jj=0; jj<jx; j++, jj++){
			for (i=gapy, ii=0; ii<iy; i++, ii++){
				n=n+1.;
				c=*(array->ptr+j*ysize+i);
				sum=sum+c_re(c);
			}
		}
		average=sum/n;
	}
else
	{	
		average=value;		/*specified value*/	
	}
	
	/*set constants*/
	c_re(c)=average;
	c_im(c)=average;
		
	
	/*set central box pixel values*/
	for (j=gapx,jj=0; jj<jx; j++, jj++){
		for (i=gapy, ii=0; ii<iy; i++, ii++){
			*(array->ptr+j*ysize+i)=c;
		}
	}

}	

void punch_hole(array_type *array, double radius, double value)
/*************************************************************
*	alex lehar 2 dec 01
*	take the input dynamically resizable input complex pair array
*	and  overwrite a disc with a
*	gray 'value'. If the value is negative, then 'value' is
*	determined by the average value of the pixels inside the centrally
*	located circle defined by radius.  If positive, then determined by 
*	value itself. Disc area is set in both real and imaginary 
*	parts.
**************************************************************/
{
	int i,j,xsize,ysize,ii,jj, hj, hi, gapx, gapy, id;
	fftw_complex c;
	double sum,n,average,d;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;
	
	hj=xsize/2;		/*image half size*/
	hi=ysize/2;
	id=(int)radius*2;
	
	gapx=(double)((xsize-(int)radius)/2); /*strip width*/
	gapy=(double)((ysize-(int)radius)/2);
	
if (value < 0.)	/*automatic determination*/
	{
		/*find the average value in the disc*/
		sum=0.;
		n=0.;
		for (j=gapx,jj=0; jj<id; j++, jj++){
			for (i=gapy, ii=0; ii<id; i++, ii++){
				d=sqrt((j-hj)*(j-hj)+(i-hi)*(i-hi));
				if(d <= radius)
					{
					n=n+1.;
					c=*(array->ptr+j*ysize+i);
					sum=sum+c_re(c);
					}
			}
		}
		average=sum/n;
	}
else
	{	
		average=value;		/*specified value*/	
	}
	
	/*set constants*/
	c_re(c)=average;
	c_im(c)=average;
		
	
	/*set disc values*/
		for (j=gapx,jj=0; jj<id; j++, jj++){
			for (i=gapy, ii=0; ii<id; i++, ii++){
				d=sqrt((j-hj)*(j-hj)+(i-hi)*(i-hi));
				if(d <= radius)
					{
					*(array->ptr+j*ysize+i)=c;
					}
			}
		}

}	

void psf_lin_complex(array_type *array, float length, float angle)
/************************************************************
*	alex lehar	22 oct 01
*	Insert into the center of a rectangular complex array structure
*	a linear wall function of specified length and angle (radians), into the real part.
*	with special provision made for vertical case of infinite gradient
*	modified extensively by Matt Lehar to fix special case of vertical and horiz.
***************************************************************/
{
	int i,j,x,c_x,y,c_y,ia,z1,z2,z3,z4,z5,z6,   /* */ didItTrip;
	float x1,y1,x2,y2,x0,y0,l2,k,a,b,c,theta;
	fftw_complex cc;
	#define PI 3.1415927
	#define THICK .51	/*this is the thickness of the line*/
	#define OCULUSMAXFLOAT 1.0E36
	
	/*unpack struct into local variables*/
	x= array->xsize;
	y= array->ysize;
	l2=length/2.0; /*halflength*/
	c_x=x/2;	/*center point*/
	c_y=y/2;
	
	didItTrip = 0;
	ia = 0;
	
	/*here we have to worry about relative sizes of x1,y1 and x2,y2 for the 
	different gradient scenarios.*/

	
	if 	( angle == 0. || fabsf(angle) >= PI )
		{
		ia=1;
		if ( fabsf(angle) >= PI )
			{
			angle = 0.;
			}
		x1 = c_x-l2;
		x2 = c_x+l2;
		y1 = c_y;
		y2 = c_y;
		
		k=(y1-y2)/(x1-x2);		/*gradient*/
		c=(y2*x1-y1*x2)/(x1-x2);/*intercept see p 11*/
		}
		
	if	( (0. < fabsf(angle)) && (PI/2. > fabsf(angle)))
		{
		ia=1;					/*flags small angle case*/
		x1=c_x-(l2*cos(angle));
		y1=c_y-(l2*sin(angle));/*one end of line*/
		x2=c_x+(l2*cos(angle));
		y2=c_y+(l2*sin(angle));/*other end of line*/
		
		k=(y1-y2)/(x1-x2);		/*gradient*/
		c=(y2*x1-y1*x2)/(x1-x2);/*intercept see p 11*/
		}
			
	if ((PI/2. < fabsf(angle)) && (PI > fabsf(angle)))
		{
		ia=2;					/*flags large angle case*/
		theta=PI-angle;
		x1=c_x-(l2*cos(theta));
		y1=c_y+(l2*sin(theta));/*other end of line*/
		x2=c_x+(l2*cos(theta));
		y2=c_y-(l2*sin(theta));/*one end of line*/
		
		k=(y1-y2)/(x1-x2);		/*gradient*/
		c=(y2*x1-y1*x2)/(x1-x2);/*intercept see p 11*/
		}
		
	if ( angle == PI/2.0 || fabsf(x2-x1) < 1.0 )  /*deal with special case*/
		{
		ia=1;
		//didItTrip = 2;
		x1 = c_x;
		x2 = c_x;
		y1 = c_y-l2;
		y2 = c_y+l2;
		
		//ia=3;					/*infinite gradient case*/
		k=OCULUSMAXFLOAT;
		c=OCULUSMAXFLOAT;
		}
		
	
	/*go through the array*/
	/*and set the necessary elements in the spanning box.*/
	for (j=0; j<x; j++) {	/*-------------------------------------------outer for*/
		for (i=0;i<y; i++){	 /*--------------------------------------------inner for*/
			/*step through each point and if within line terminus
			range, and if within line thickness set to 1.0, other-
			wise set to zero*/
				{		/*for block*/
					
					/* In the following 4 lines, the bounding box is set to be at least as
							thick as the desired line */
					if( ia == 1 )
						{
						z1=maximum_twonums(((float)j > (((x1+x2)/2)-THICK)), ((float)j > x1)); /*x past line start*/
						z2=maximum_twonums(((float)j < (((x1+x2)/2)+THICK)), ((float)j < x2)); /*x before line end*/
						z3=maximum_twonums(((float)i > (((y1+y2)/2)-THICK)), ((float)i > y1)); /*y past line start*/
						z4=maximum_twonums(((float)i < (((y1+y2)/2)+THICK)), ((float)i < y2)); /*y before line end*/
						z5=fabsf(x2-x1) < 1.;	/*steep gradient extreme*/
						}
						
					if( ia == 2 )
						{
						z1=maximum_twonums(((float)j > (((x1+x2)/2)-THICK)), ((float)j > x1)); /*x past line start*/
						z2=maximum_twonums(((float)j < (((x1+x2)/2)+THICK)), ((float)j < x2)); /*x before line end*/
						z3=maximum_twonums(((float)i > (((y1+y2)/2)-THICK)), ((float)i < y1)); /*y past line start*/
						z4=maximum_twonums(((float)i < (((y1+y2)/2)+THICK)), ((float)i > y2)); /*y before line end*/
						z5=fabsf(x2-x1) < 1.;	/*steep gradient extreme*/
						}
					
					if(((z1 && z2) && (z3 && z4))) /*point is in box range*/
						{
							if ( fabsf(x2-x1) >= 1.0)
								{
								a=(float)(k*((float)i-c)+(float)j)/(1.0+k*k);
								x0=a;
								y0=a*k+c;
								b=(float)sqrt((x0-(float)j)*
								(x0-(float)j)+(y0-(float)i)*(y0-(float)i));	/*distance*/
								}
							else
								{
								b=0.;		/*special case of vertical*/
								}
							if (b < THICK)
								{
								//didItTrip = 1;
								c_re(cc)=1.0;
								c_im(cc)=0.0;
								*(array->ptr+j*y+i)=cc;	/*reassignment*/
								}
							else
								{
								c_re(cc)=0.0;
								c_im(cc)=0.0;
								*(array->ptr+j*y+i)=cc;	/*reassignment*/
								}
						}	
						
				else
					{
					c_re(cc)=0.0;
					c_im(cc)=0.0;
					*(array->ptr+j*y+i)=cc;		/*reassignment*/
					}
				}/*for block*/
			}/*---------------------------------------------------------inner for*/
		}/*--------------------------------------------------------outer for*/
	}	

void psf_gen(array_type *array, float length, float angle, int whichone)
/************************************************************
*	alex lehar 3 jan 2002
*	generate either a linear or circular focus blur psf in 
*	a complex array, depending on whether WHICHONE is 1 or 0
*************************************************************/
	{
	if (whichone ==0)
		{
		psf_oof_complex(array, length);
		}
	else
		{
		psf_lin_complex(array, length, angle);
		}
	}
#if 0	
void leave_hole(array_type *array, double radius, double value)
/*************************************************************
*	alex lehar 5 jan 02
*	take the input dynamically resizable input complex pair array
*	and  overwrite outside a disc with a
*	gray 'value'. If the value is negative, then 'value' is
*	determined by the average value of the pixels inside the centrally
*	located circle defined by radius.  If positive, then determined by 
*	value itself. Disc area is set in both real and imaginary 
*	parts.
**************************************************************/
{
	int i,j,xsize,ysize,ii,jj, hj, hi, gapx, gapy, id;
	fftw_complex c;
	double sum,n,average,d;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;
	
	hj=xsize/2;		/*image half size*/
	hi=ysize/2;
	id=(int)radius*2;
	
	gapx=(double)((xsize-(int)radius)/2); /*strip width*/
	gapy=(double)((ysize-(int)radius)/2);
	
if (value < 0.)	/*automatic determination*/
	{
		/*find the average value in the disc*/
		sum=0.;
		n=0.;
		for (j=gapx,jj=0; jj<id; j++, jj++){
			for (i=gapy, ii=0; ii<id; i++, ii++){
				d=sqrt((j-hj)*(j-hj)+(i-hi)*(i-hi));
				if(d <= radius)
					{
					n=n+1.;
					c=*(array->ptr+j*ysize+i);
					sum=sum+c_re(c);
					}
			}
		}
		average=sum/n;
	}
else
	{	
		average=value;		/*specified value*/	
	}
	
	/*set constants*/
	c_re(c)=average;
	c_im(c)=average;
		
	
	/*set disc values*/
		for (j=gapx,jj=0; jj<id; j++, jj++){
			for (i=gapy, ii=0; ii<id; i++, ii++){
				d=sqrt((j-hj)*(j-hj)+(i-hi)*(i-hi));
				if(d <= radius)
					{
					*(array->ptr+j*ysize+i)=c;
					}
			}
		}

}	
#endif


// debug stuff
void delta_array_complex(array_type *array)
/*************************************************************
*   9 feb 05	alex lehar  Picturesolve.com
*	take the input dynamically resizable input complex pair array
*	and zero it. Then load a single value of 1.0 in the real part
*   at the center of the array, to form a "delta" function.  
**************************************************************/
{
	int i,j,xsize, ysize, xh, yh;
	fftw_complex c;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;
	
	//calculate the half points
	xh=xsize/2;
	yh=ysize/2;
	c_re(c)=0.0; //zero the complex value
	c_im(c)=0.0;
	
	/*set the values*/
	for (j=0;j<xsize;j++){
		for(i=0;i<ysize; i++){
			*(array->ptr+j*ysize+i)=c;
		}
	}
	c_re(c)=1.0;	//make the delta
	*(array->ptr+xh*ysize+yh)=c;
}

void rampup_array_complex(array_type *array)
/*************************************************************
*   9 feb 05	alex lehar  Picturesolve.com
*	take the input dynamically resizable input complex pair array
*	and enter data starting 0 increment up in raster fashion
*   through whole picture for test purposes.
**************************************************************/
{
	int i,j,xsize, ysize, k;
	fftw_complex c;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;
	
	c_re(c)=0.0; //zero the complex value
	c_im(c)=0.0;
	k=0;
	
	/*set the values*/
	for (j=0;j<xsize;j++){
		for(i=0;i<ysize; i++){
			k=k+1;
			c_re(c)=(float)k;
			*(array->ptr+j*ysize+i)=c;
		}
	}
}

// end debug stuff



void psf_linthick_complex(array_type *array, float length, float angle, float thick)
/************************************************************
*	alex lehar	22 oct 01
*   modified from psf_lin_complex, to add another argument for thickness, "thick"
*   14 feb 05
*
*	Insert into the center of a rectangular complex array structure
*	a linear wall function of specified length and angle (radians), into the real part.
*	with special provision made for vertical case of infinite gradient
*	modified extensively by Matt Lehar to fix special case of vertical and horiz.
***************************************************************/
{
	int i,j,x,c_x,y,c_y,ia,z1,z2,z3,z4,z5,z6,   /* */ didItTrip;
	float x1,y1,x2,y2,x0,y0,l2,k,a,b,c,theta;
	fftw_complex cc;
	#define PI 3.1415927
	#define OCULUSMAXFLOAT 1.0E36
	
	/*unpack struct into local variables*/
	x= array->xsize;
	y= array->ysize;
	l2=length/2.0; /*halflength*/
	c_x=x/2;	/*center point*/
	c_y=y/2;
	
	didItTrip = 0;
	ia = 0;
	
	/*here we have to worry about relative sizes of x1,y1 and x2,y2 for the 
	different gradient scenarios.*/

	
	if 	( angle == 0. || fabsf(angle) >= PI )
		{
		ia=1;
		if ( fabsf(angle) >= PI )
			{
			angle = 0.;
			}
		x1 = c_x-l2;
		x2 = c_x+l2;
		y1 = c_y;
		y2 = c_y;
		
		k=(y1-y2)/(x1-x2);		/*gradient*/
		c=(y2*x1-y1*x2)/(x1-x2);/*intercept see p 11*/
		}
		
	if	( (0. < fabsf(angle)) && (PI/2. > fabsf(angle)))
		{
		ia=1;					/*flags small angle case*/
		x1=c_x-(l2*cos(angle));
		y1=c_y-(l2*sin(angle));/*one end of line*/
		x2=c_x+(l2*cos(angle));
		y2=c_y+(l2*sin(angle));/*other end of line*/
		
		k=(y1-y2)/(x1-x2);		/*gradient*/
		c=(y2*x1-y1*x2)/(x1-x2);/*intercept see p 11*/
		}
			
	if ((PI/2. < fabsf(angle)) && (PI > fabsf(angle)))
		{
		ia=2;					/*flags large angle case*/
		theta=PI-angle;
		x1=c_x-(l2*cos(theta));
		y1=c_y+(l2*sin(theta));/*other end of line*/
		x2=c_x+(l2*cos(theta));
		y2=c_y-(l2*sin(theta));/*one end of line*/
		
		k=(y1-y2)/(x1-x2);		/*gradient*/
		c=(y2*x1-y1*x2)/(x1-x2);/*intercept see p 11*/
		}
		
	if ( angle == PI/2.0 || fabsf(x2-x1) < 1.0 )  /*deal with special case*/
		{
		ia=1;
		//didItTrip = 2;
		x1 = c_x;
		x2 = c_x;
		y1 = c_y-l2;
		y2 = c_y+l2;
		
		//ia=3;					/*infinite gradient case*/
		k=OCULUSMAXFLOAT;
		c=OCULUSMAXFLOAT;
		}
		
	
	/*go through the array*/
	/*and set the necessary elements in the spanning box.*/
	for (j=0; j<x; j++) {	/*-------------------------------------------outer for*/
		for (i=0;i<y; i++){	 /*--------------------------------------------inner for*/
			/*step through each point and if within line terminus
			range, and if within line thickness set to 1.0, other-
			wise set to zero*/
				{		/*for block*/
					
					/* In the following 4 lines, the bounding box is set to be at least as
							thick as the desired line */
					if( ia == 1 )
						{
						z1=maximum_twonums(((float)j > (((x1+x2)/2)-thick)), ((float)j > x1)); /*x past line start*/
						z2=maximum_twonums(((float)j < (((x1+x2)/2)+thick)), ((float)j < x2)); /*x before line end*/
						z3=maximum_twonums(((float)i > (((y1+y2)/2)-thick)), ((float)i > y1)); /*y past line start*/
						z4=maximum_twonums(((float)i < (((y1+y2)/2)+thick)), ((float)i < y2)); /*y before line end*/
						z5=fabsf(x2-x1) < 1.;	/*steep gradient extreme*/
						}
						
					if( ia == 2 )
						{
						z1=maximum_twonums(((float)j > (((x1+x2)/2)-thick)), ((float)j > x1)); /*x past line start*/
						z2=maximum_twonums(((float)j < (((x1+x2)/2)+thick)), ((float)j < x2)); /*x before line end*/
						z3=maximum_twonums(((float)i > (((y1+y2)/2)-thick)), ((float)i < y1)); /*y past line start*/
						z4=maximum_twonums(((float)i < (((y1+y2)/2)+thick)), ((float)i > y2)); /*y before line end*/
						z5=fabsf(x2-x1) < 1.;	/*steep gradient extreme*/
						}
					
					if(((z1 && z2) && (z3 && z4))) /*point is in box range*/
						{
							if ( fabsf(x2-x1) >= 1.0)
								{
								a=(float)(k*((float)i-c)+(float)j)/(1.0+k*k);
								x0=a;
								y0=a*k+c;
								b=(float)sqrt((x0-(float)j)*
								(x0-(float)j)+(y0-(float)i)*(y0-(float)i));	/*distance*/
								}
							else
								{
								b=0.;		/*special case of vertical*/
								}
							if (b < thick)
								{
								//didItTrip = 1;
								c_re(cc)=1.0;
								c_im(cc)=0.0;
								*(array->ptr+j*y+i)=cc;	/*reassignment*/
								}
							else
								{
								c_re(cc)=0.0;
								c_im(cc)=0.0;
								*(array->ptr+j*y+i)=cc;	/*reassignment*/
								}
						}	
						
				else
					{
					c_re(cc)=0.0;
					c_im(cc)=0.0;
					*(array->ptr+j*y+i)=cc;		/*reassignment*/
					}
				}/*for block*/
			}/*---------------------------------------------------------inner for*/
		}/*--------------------------------------------------------outer for*/
	}	

void psf_genthick(array_type *array, float length, float angle, float thick, int whichone, float maxmod, float nstds)
/************************************************************
*   modified from psf_gen, to accommodate extra argument "thick"
*   14 feb 05
*   added ellipse option 15 feb 05.  Written by Matt.
*   add further value for whichone =3 to generate grayscale
*   ellipse, with gaussian distribution
*
*	alex lehar 3 jan 2002
*	generate either a linear or circular focus blur psf in 
*	a complex array, depending on whether WHICHONE is 1 or 0
*	maxmod is a floating value 0.to 1. to modulate psf center
*	21 apr 05, add nstds, number of standard deviations for intensity modulation
*************************************************************/
	{
	if (whichone ==0)
		{
		psf_oof_complex(array, length);
		}
	else if (whichone ==1)
		{
		psf_linthick_complex(array, length, angle, thick);
		}
	else if (whichone ==2)
		{
		psf_ellipse_complex(array, length, thick, angle);
		}
	else if (whichone ==3)
		{
		psf_gaussellipse(array, length, thick, angle, nstds);
		}
	else
		{
		psf_invert_gaussellipse(array, length, thick, angle, nstds, maxmod);
		}
	}

int psf_ellipse_complex(array_type *array, float length, float width, float theta)
/************************************************************
*	alex & matt lehar	15 feb 05
*	Insert into the center of a rectangular complex array structure
*	an ellipse of specified proportion and orientation, into the real part.
***************************************************************/
{
	int i,j,x,y,z,test;
	float d,a,b;
	fftw_complex c;
	#define PI 3.1415927
	
	//redefine width as an aspect ratio
	width=length*width;
	
	/*unpack struct into local variables*/
	x= array->xsize;
	y= array->ysize;
	z=minimum_twonums(x,y);	/*compute the minimum*/
	/*make sure the radius is not too big for the array*/
	if(2.0*length < (float)z)
		{
		/*go through the array*/
		/*and set the necessary elements in the circle.*/
		for (j=0; j<x; j++) {
			for (i=0;i<y; i++){
				/*distance between center and examined point*/
				test=oneIfWithinEllipse(j,i,x,y,length,width,theta);
				c=*(array->ptr+j*y+i);
				if	(test > 0)
					{
					c_re(c)=1.0;
					c_im(c)=0.0;
					*(array->ptr+j*y+i)=c;	/*reassignment*/
					}
				else
					{
					c_re(c)=0.0;
					c_im(c)=0.0;
					*(array->ptr+j*y+i)=c;	/*reassignment*/
					}

			}
		}
		}
	else
		return (0);
	}
int psf_gaussellipse(array_type *array, float length, float width, float theta, float nstd)
/************************************************************
*	alex & matt lehar	20 feb 05
*	Insert into the center of a rectangular complex array structure
*	an ellipse of specified proportion and orientation, into the real part.
*   the ellipse has a gaussian distribution intensity variation
***************************************************************/
{
	int i,j,x,y,z;
	float d,a,b,test;
	fftw_complex c;
	#define PI 3.1415927
	
	//redefine width as an aspect ratio
	width=length*width;

	
	/*unpack struct into local variables*/
	x= array->xsize;
	y= array->ysize;
	z=minimum_twonums(x,y);	/*compute the minimum*/
	
	/*make sure the radius is not too big for the array*/
	if(2.0*length < (float)z)
		{
		/*go through the array*/
		/*and set the necessary elements in the circle.*/
		for (j=0; j<x; j++) {
			for (i=0;i<y; i++){
				/*distance between center and examined point*/
				test=gaussianEllipseValue(nstd,j,i,x,y,length,width,theta);
				c=*(array->ptr+j*y+i);
					c_re(c)=test;
					c_im(c)=0.0;
					*(array->ptr+j*y+i)=c;	/*reassignment*/
			}
		}
		}
	else
		return (0);
	}
	
int psf_invert_gaussellipse(array_type *array, float length, float width, float theta, float nstd, float max_mod)
/************************************************************
*	alex & matt lehar	20 feb 05
*	Insert into the center of a rectangular complex array structure
*	an ellipse of specified proportion and orientation, into the real part.
*   the ellipse has a gaussian distribution intensity variation
*	max_mod is floating point value 0. to 1. giving totality of blackness at center of ellipse
*	nstd is std deviations between ctr and rim.
***************************************************************/
{
	int i,j,x,y,z;
	float d,a,b,test;
	fftw_complex c;
	#define PI 3.1415927
	
	//redefine width as an aspect ratio
	width=length*width;

	
	/*unpack struct into local variables*/
	x= array->xsize;
	y= array->ysize;
	z=minimum_twonums(x,y);	/*compute the minimum*/
	
	/*make sure the radius is not too big for the array*/
	if(2.0*length < (float)z)
		{
		/*go through the array*/
		/*and set the necessary elements in the circle.*/
		for (j=0; j<x; j++) {
			for (i=0;i<y; i++){
				/*distance between center and examined point*/
				test=oneIfWithinEllipse(j,i,x,y,length,width,theta)*(255.0-max_mod*gaussianEllipseValue(nstd,j,i,x,y,length,width,theta));
				c=*(array->ptr+j*y+i);
					c_re(c)=test;
					c_im(c)=0.0;
					*(array->ptr+j*y+i)=c;	/*reassignment*/
			}
		}
		}
	else
		return (0);
	}	
	
	
// written by Matt Lehar, 2.20.05
// returns a number between 0 and 255 corresponding to the value of the Gaussian at the given point
// the argument "numStdDevs" provides the number of standard deviations between the center and the rim of the ellipse
//
float gaussianEllipseValue( float numStdDevs, int x_prime, int y_prime, int n_x, int n_y, int length, int width, float rotation )
  {
  int x = x_prime - n_x/2;
  int y = -y_prime + n_y/2;
  
  if (x==0 && y==0) { return 255.0; } //special case of center point

  float R_o = sqrt( 1 / (1/pow(length,2) * pow(cos( theta_1(x,y)+rotation ),2) + 1/pow(width,2) * pow(sin( theta_1(x,y)+rotation ),2)));

  return floor( 256.0 * exp(-radius_squared(x,y) / (2*pow(R_o/numStdDevs,2))));
  }


//this code is written by matt lehar, for establishing the test criterion for the ellipse
//	
int oneIfWithinEllipse( int x_prime, int y_prime, int n_x, int n_y, int length, int width, float rotation )
  {
  //length=20;
  //width=15;
  //rotation=0;
  
  int x = x_prime - n_x/2;
  int y = -y_prime + n_y/2;
  
  if (x==0 && y==0) { return 1; } //special case of center point
//printstatements
		//printf("hello from matt \n");
		//printf("   %2.2f %2.2f \n", extent, angle);
		//printf("   %d %d \n", theta_1(3,3), width);
  if ( radius_squared(x,y) <= 1 / (1/pow(length,2) * pow(cos( theta_1(x,y)+rotation ),2) + 1/pow(width,2) * pow(sin( theta_1(x,y)+rotation ),2)))
    {
    return 1;
    }
  else { return 0; }
  }


int radius_squared( int x, int y )
  {
  return pow(x,2) + pow(y,2);
  }

float theta_1( int x, int y )
  {
  float cos_result = acos( x / (sqrt(radius_squared( x, y ))));
  float sin_result = asin( y / (sqrt(radius_squared( x, y ))));
  
  if (y>0) {return cos_result; }
  else if (x>0) { return sin_result; }
  else { return acos( -x / (sqrt(radius_squared( x, y )))) + PI; }
  }

void make_fuzzy(array_type *array, array_type *brray, float radius)
/****************************************************************************
*   make_fuzzy	Alex Lehar  25 Feb 05
*   take an array in input "array", convolve it with a disc of "radius" pixels
*   and put the result in "brray"
*****************************************************************************/
{
	//declare internal buffer
	array_type *crray;
	int x,y;
	x= array->xsize;
	y= array->ysize;
	
	//allocate memory for the internal buffer crray
	crray = (array_type*)malloc(sizeof(array_type));
	crray->xsize=x;
	crray->ysize=y;
	crray->ptr = (fftw_complex *)malloc(x*y*sizeof(fftw_complex));

	//generate the convolving disc
	psf_genthick(brray, radius, 0., 1.0, 0,1.,0.0);
	
	//do the convolution
	cnvolv_spatial_to_spatial(array, brray, crray, 0);
	move_array(crray, brray);
	
	//release local memory
	free (crray->ptr);
	free (crray);
}	

void rescale_array_complex(array_type *array, float rescale)
/*************************************************************
*	alex lehar 29 may 05
*	take the input dynamically resizable input complex pair array
*	and re-scale it by "rescale"
*	modified 9 nov, to involve sums of magnitudes.
**************************************************************/
{
	int i,j,xsize,ysize;
	fftw_complex c;
	double r,s;
	
	/*unpack struct into local variables*/
	xsize= array->xsize;
	ysize= array->ysize;
	
	
	/*do the rescaling*/
	for (j=0;j<xsize;j++){
		for(i=0;i<ysize; i++){
			c=*(array->ptr+j*ysize+i);
			r=c_re(c);
			s=c_im(c);
			c_re(c)=r*rescale;
			c_im(c)=s*rescale;
			*(array->ptr+j*ysize+i)=c;
		}
	}
	
}	
