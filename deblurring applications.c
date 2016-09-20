//this code is copyright Picturesolve.com
//#include	"ZJPEGFile.h"
//#include	"ZGWorld.h"
//#include	"MacZoop.h"

//#include	<ImageCompression.h>
//#include	<FixMath.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fftw.h>
#include <time.h>
//#include <fftw-int.h>
#include <deblurring functions.h>
#include <deblurring applications.h>
//#include	"ZExpParser.h"
//copyToWindowArray 

/*
void copyToWindowArray(unsigned char *inputArray1, 
						unsigned char *inputArray2, 
						unsigned char *inputArray3, 
						int colorImage,
						char * baseAddressImage,
						int heightInputArray, int widthInputArray,
						int lengthPixelMapLine)
	{
	int i,j,ii,offset,offsetInputArray;
	offsetInputArray=0;
	offset=0;
	//sets the pixel values in the image to be displayed to the scaled values as determined
	//by the fft array -- baseAddressImage is the array used by the window
	for (j=0;j<heightInputArray;j++)
		{
		for (i=offsetInputArray,ii=offset;i<offsetInputArray+widthInputArray;ii+=3,i+=1)
			{
			unsigned char pixelvalue;	
			if(colorImage==1)
				{
				pixelvalue=(unsigned char)inputArray1[i];
				baseAddressImage[ii+1]=( char)pixelvalue;
				pixelvalue=(unsigned char)inputArray2[i];
				baseAddressImage[ii+2]=( char)pixelvalue;
				pixelvalue=(unsigned char)inputArray3[i];
				baseAddressImage[ii+3]=( char)pixelvalue;
				}
			else
				{
				pixelvalue=(unsigned char)inputArray1[i];
				baseAddressImage[ii+1]=( char)pixelvalue;
				baseAddressImage[ii+2]=( char)pixelvalue;
				baseAddressImage[ii+3]=( char)pixelvalue;
				}
			}
		offset+=lengthPixelMapLine;
		offsetInputArray+=widthInputArray;
		}
	}
*/
void copyToWindowArray(unsigned char *inputArray1, 
					   unsigned char *inputArray2, 
					   unsigned char *inputArray3, 
					   int processColors,
					   char *baseAddressImage,
					   int samplesPerPixel, int height, int width)
{
	int i;
	for(i=0; i<width*height; i++)
	{
		char *pixel=baseAddressImage+samplesPerPixel*i;
		if(processColors)
		{
			pixel[0]=(unsigned char)inputArray1[i];
			pixel[1]=(unsigned char)inputArray2[i];
			pixel[2]=(unsigned char)inputArray3[i];
		}
		else
		{
			if(samplesPerPixel==3)
			{
				pixel[0]=(unsigned char)inputArray1[i];
				pixel[1]=(unsigned char)inputArray1[i];
				pixel[2]=(unsigned char)inputArray1[i];
			}
			else if(samplesPerPixel==1)
				pixel[0]=(unsigned char)inputArray1[i];
		}
	}
}

void processImage(short samplesPerPixel, short width, short height,
				  char *baseAddressImage, char *baseAddressArray, DEBLURINPUT *deblurInput)
{
	//ZExpParser	ep;
	short temp;
	long i,ii,j;
	int flags;
	fftw_complex putin;

	int		Int_1, Int_2, Int_3, Int_4;
	float	Float_1,Float_2,Float_3,Float_4,Float_5,Float_6,Float_7,
			Float_8,Float_9,Float_10,Float_11,Float_12,Float_13,Float_14;
	int		processColors,numberColors,l;
	
	//declarations for stretching
	double hilev, lolev, v;
	
	
	// take account of any mouse-driven change in the value "a"
	
	//ep.SetSymbolValue( "a", aa );
	/*arrange for clock*/
	//double t;
	//fftw_time begin, end;
	
	//fftwnd_plan p;
	array_type *array;	/*this will contain the image*/
	array_type *brray;	/*this will contain the psf*/
	array_type *crray;	/*this will contain the output*/
	array_type *drray;	/*this will contain extra working array*/
	
	/*declarations for transform*/
	fftw_direction isign;
#if 1	
	//read the input from the modeless dialog for use if wanted
	Int_1=deblurInput->inputint1;
	Int_2=deblurInput->inputint2;
	Int_3=deblurInput->inputint3;
	Int_4=deblurInput->inputint4;
	Float_1=deblurInput->inputfloat1;
	Float_2=deblurInput->inputfloat2;
	Float_3=deblurInput->inputfloat3;
	Float_4=deblurInput->inputfloat4;
	Float_5=deblurInput->inputfloat5;
	Float_6=deblurInput->inputfloat6;
	Float_7=deblurInput->inputfloat7;
	Float_8=deblurInput->inputfloat8;
	Float_9=deblurInput->inputfloat9;
	Float_10=deblurInput->inputfloat10;
	Float_11=deblurInput->inputfloat11;
	Float_12=deblurInput->inputfloat12;
	Float_13=deblurInput->inputfloat13;
	Float_14=deblurInput->inputfloat14;
#endif
	
	//printf("%f, %f, %f, %f, %f, %f, %f\n", Float_1, Float_2, Float_3, Float_4, Float_5, Float_6, Float_7);
	
	if(Int_4==1)
		processColors=1;
	else
		processColors=0;
	
  /*   fftw_time begin, end;*/
	flags  = FFTW_IN_PLACE;
		/*load the variables*/
	isign=FFTW_FORWARD;  /*this gives the Forward Transform*/
	unsigned char *workArray1;
	unsigned char *workArray2;
	unsigned char *workArray3;
	unsigned char *workArray4;
	//unsigned char *workBrray;		/*for the block power spectrum data*/
	//unsigned long nn[3];
	long widthFFT,heightFFT,offset,offsetFFT;
	short lengthPixelMapLine;
	lengthPixelMapLine=width;

#if 0
	for (i=0;i<width*height*3;i++)
	{
		baseAddressImage[i]^=255;
	}
#endif
	//baseAddressImage=malloc(10000);
//-----------------------------------------------------------------deblurring
	//just set the size for fft as the image size
	widthFFT=width;
	heightFFT=height;	
	/* allocate the array struct **/
	/*this allocates only the space for two ints and a pointer*/
	array = (array_type*)malloc(sizeof(array_type));
	brray = (array_type*)malloc(sizeof(array_type));/*for the psf*/
	crray = (array_type*)malloc(sizeof(array_type));/*for the output*/
	drray = (array_type*)malloc(sizeof(array_type));/*for extra working array*/
	
	/*set the size variables*/
	array->xsize = heightFFT; //Alex swapped these around
	array->ysize = widthFFT;
	/*and for the extra arrays too*/
	brray->xsize = heightFFT;
	brray->ysize = widthFFT;
	crray->xsize = heightFFT;
	crray->ysize = widthFFT;
	/*for working array*/
	drray->xsize = heightFFT;
	drray->ysize = widthFFT;
	
	//nn[2]=heightFFT;
	array->ptr = (fftw_complex *)malloc(widthFFT*heightFFT*sizeof(fftw_complex));
	brray->ptr = (fftw_complex *)malloc(widthFFT*heightFFT*sizeof(fftw_complex)); /*for psf*/
	crray->ptr = (fftw_complex *)malloc(widthFFT*heightFFT*sizeof(fftw_complex)); /*for output*/
	drray->ptr = (fftw_complex *)malloc(widthFFT*heightFFT*sizeof(fftw_complex)); /*for extra working array*/

	workArray1 = (unsigned char *)malloc(widthFFT*heightFFT);
	workArray2 = (unsigned char *)malloc(widthFFT*heightFFT);
	workArray3 = (unsigned char *)malloc(widthFFT*heightFFT);
	workArray4 = (unsigned char *)malloc(widthFFT*heightFFT);
	//alex swapped heightFFT, and widthFFT parameters around here
	//p= fftw2d_create_plan(heightFFT,widthFFT,isign,flags);
#if 0
	temp=(*ppix)->rowBytes;
	lengthPixelMapLine=temp&0x3FFF;
#endif
#if 1
	if(processColors==1)
		numberColors=3;
	else
		numberColors=1;
	for(l=0;l<numberColors;l++)
	{
		//loops through image, averages RGB components of pixel, and set the fft array real value to
		//the average value. Sets the imaginary value to 0
		/*	
		offsetFFT=0;
		offset=0;
		for (j=0;j<heightFFT;j++)
			{
			for (i=offsetFFT,ii=offset;i<offsetFFT+widthFFT;ii+=3,i+=1)
				{
				float r;	// alex changed this from int to float
				if(processColors==1)
					{
					if(l==0)
						r=(float)((unsigned char)baseAddressImage[ii+1]);
					else if(l==1)
						r=(float)((unsigned char)baseAddressImage[ii+2]);
					else if(l==2)
						r=(float)((unsigned char)baseAddressImage[ii+3]);
					}
				else
					r=(float)((unsigned char)baseAddressImage[ii+1]+(unsigned char)baseAddressImage[ii+2]+(unsigned char)baseAddressImage[ii+3])/3.0;
				c_re(putin)=(double)r;		
				c_im(putin)=0.0;
				*(array->ptr+i) = putin;				//assign complex
				}
			offset+=lengthPixelMapLine;
			offsetFFT+=widthFFT;
			}
		 */
		int i;
		for(i=0; i<width*height; i++)
		{
			char *pixel=baseAddressImage+samplesPerPixel*i;
			float r;
			if(processColors)
				r = (float)pixel[l];
			else
			{
				if(samplesPerPixel==3)
					r = ((float)pixel[0]+(float)pixel[1]+(float)pixel[2])/3.0;
				else if(samplesPerPixel==1)
					r = (float)pixel[0];
			}
			c_re(putin)=(double)r;		
			c_im(putin)=0.0;
			*(array->ptr+i) = putin;	//assign complex
		}
#endif


	

//-----------------------------------------------------------------deblurring		
		//alex,put in code to manipulate workArray, which is organized as an ordinary byte image,
		//that is, each pixel is a single byte, which takes on a value from 0 to 255, and those 
		//bytes are organized as usual in sequential rows of pixels


//alex||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||testing

/*these routines can be switched on INDIVIDUALLY, and all make the following
assumptions:
	~	The input array pointer is *array
	~	The output image displayed, has the same size as *array
		and is designated *brray        */

//insertion point for new stuff
#define CEPSTRUM 0
#define ARTIFICIAL_CEPSTRUM 1
#define ARTIFICIAL_AVERAGE_POWER_SPECTRUM 2
#define PSF_POWERSPECTRUM 3
#define AVERAGE_POWER_SPECTRUM 4
#define ADD_NOISE 5
#define	PSF_BLUR_IMAGE 6
#define TWO_PSFS_CONVOLVE 7
#define SIMULATED_RESTORATION 8
#define EDGING_DEMO 9
#define	DEBLUR 10
#define DEBLUR_WITH_WINDOWS 11
#define SHOW_MAP_KERNEL 12
#define	DISPLAY_PSF 13
#define DISPLAY_PATGEN_MAGNITUDE 14
#define PATTERN_REMOVAL 15
#define FFT 16
#define TEST_NYQUIST 17
#define BANDPASS_FILTERING 18
#define STRETCH 19


	if (Int_3==CEPSTRUM)
		{

		//---0 demonstrate Cepstral Estimation
		cepest(array, Int_1,Int_2, brray, 1);
		}
		
	else if (Int_3==ARTIFICIAL_CEPSTRUM)
		{
		//---1 cepstrum of artificially blurred image:
		//blur an image with a point spread function
		#define PI 3.1415927
		float extent, angle;
		extent=Float_1;
		angle=PI*.7;
		psf_gen(brray, extent, angle, 0);
		cnvolv_spatial_to_spatial(array, brray, crray, 0);
		cepest(crray, Int_1,Int_2, brray, 1);
		}

	else if (Int_3==ARTIFICIAL_AVERAGE_POWER_SPECTRUM)
		{
		//---2 average power spectrum of artificially blurred image:
		//blur an image with a point spread function
		#define PI 3.1415927
		float extent, angle;
		extent=30.;
		angle=PI*.7;
		int jx, iy;
		double high, low;
		psf_oof_complex(brray, extent);
		cnvolv_spatial_to_spatial(array, brray, crray, 0);
		average_power_spectrum(crray, Int_1, Int_2, brray, 1, &jx, &iy);
		maxmin_array_special(brray, &high, &low, 0, &jx, &iy);
		border_padding(brray, jx, iy, -4);
		stretch_complex_array(brray, high, low, 0);
		}
		
	else if (Int_3==PSF_POWERSPECTRUM)
		{
		//---3 demonstrate power spectrum of point spread function
		/*adding a thickness parameter in the case of linear motion*/
		#define PI 3.1415927
		float extent, angle, thick;
		extent=Float_1;
		angle=PI*Float_2/180.;
		thick=Float_4;
		//printf("values are\n");
		//printf("   %2.2f %2.2f \n", extent, angle);
		psf_genthick(brray, extent, angle, thick, Int_2, Float_7, Float_8);

		//check to see if want to convolve psf with disc to blur its edge
			if(Int_1 == 4)
				{
				make_fuzzy(brray, array, Float_7);
				move_array(array, brray);
				}
			else
				{
				}

		shift_array_complex(brray);
		forward_fft(brray);
		magnitude_phase_complex(brray);
		convert_to_logthreshold_complex(brray, .0001);
		}


	else if (Int_3==AVERAGE_POWER_SPECTRUM)
		{
		//---4 demonstrate average power spectrum of an image
		int jx, iy;
		average_power_spectrum(array, Int_1, Int_2, brray, 1, &jx, &iy);
		}


	else if (Int_3==ADD_NOISE)
		{
		//---5 demonstration addition of noise
		float snr;
		snr=Float_7;
		//printf("snr is \n");
		//printf("   %2.2f  \n", snr);
		add_noise_complex(array, snr);
		move_array(array, brray);
		}


	else if (Int_3==PSF_BLUR_IMAGE)
		{
		//---6 blur an image with a point spread function
		#define PI 3.1415927
		float extent, angle, thick;
		extent=Float_1;
		angle=PI*Float_2/180.;
		thick=Float_4;
		//psf_gen(brray, extent, angle, Int_2);
		psf_genthick(brray, extent, angle, thick, Int_2, Float_7, Float_8);
		cnvolv_spatial_to_spatial(array, brray, crray, 0);
		move_array(crray, brray);
		}
		

	else if (Int_3==TWO_PSFS_CONVOLVE)
		{
		//---7 convolve two different psf's together.
		#define PI 3.1415927
		float extent, angle, thick;
		extent=Float_1;
		angle=PI*Float_2/180.;
		thick=Float_4;
		psf_genthick(array, extent, angle, thick, Int_2, Float_7, Float_8);//choose the starting psf(any one)
		psf_genthick(brray, Float_7, angle, 1.0, 0, Float_7, Float_8);//choose the convolving function radius
		cnvolv_spatial_to_spatial(brray, array, crray, 0);
		move_array(crray, brray);
		}


	else if (Int_3==SIMULATED_RESTORATION)
		{
		//---8 simulated restoration
		//worked 13 nov 01
		#define PI 3.1415927
		float extent, angle, snr, thick;
		extent=Float_1;
		angle=PI*Float_2/180.;
		snr=Float_3;
		thick=Float_4;
		//generate psf for artificial blurring
		psf_genthick(brray, extent, angle, thick, Int_2, Float_7, Float_8);//generate psf
		//add_noise_complex(array, snr);		//add some noise
		//artificially blur
		cnvolv_spatial_to_spatial(array, brray, crray,0);//artificially blur input
		move_array(crray, array);
		psf_genthick(brray, extent, angle, thick, Int_2, Float_7, Float_8);//generate psf
		edge(array, brray, crray);			 //edge the picture
		psf_genthick(brray, extent, angle, thick, Int_2, Float_7, Float_8);//regenerate psf
		//check to see if want to convolve psf with disc to blur its edge
			if(Int_1 == 4)
				{
				make_fuzzy(brray, array, Float_7);
				move_array(array, brray);
				}
			else
				{
				}
		//now do the restoration
		deblur(crray, brray, array, snr);	 //nitty gritty
		}


	else if (Int_3==EDGING_DEMO)
		{
		//---9 Edging Demonstration--------------------------------------------------
		//worked 13 nov 01
		#define PI 3.1415927
		float extent, angle, thick;
		extent=Float_1;
		angle=PI*Float_2/180.;
		thick=Float_4;
		//psf_gen(brray, extent, angle, Int_2);
		psf_genthick(brray, extent, angle, thick, Int_2, Float_7, Float_8);//generate psf

		edge(array, brray, crray);			//edge the picture
		move_array(crray, brray);			//replace original picture
		}

	else if (Int_3==DEBLUR)
		{
		//---10 restore a picture----------------------------------------------------

		#define PI 3.1415927
		float extent, angle, snr, thick;
		extent=Float_1;
		angle=PI*Float_2/180.;
		snr=Float_3;
		thick=Float_4;
		
		psf_genthick(brray, extent, angle, thick, Int_2, Float_7, Float_8);//generate psf
		
		//check to see if want to convolve psf with disc to blur its edge
			if(Int_1 == 4)
				{
				make_fuzzy(brray, drray, Float_7);
				move_array(drray, brray);
				}
			else
				{
				}
		
		
		edge(array, brray, crray);			 //edge the picture
		psf_genthick(brray, extent, angle, thick, Int_2, Float_7, Float_8);//regenerate psf
		//check to see if want to convolve psf with disc to blur its edge
			if(Int_1 == 4)
				{
				make_fuzzy(brray, drray, Float_7);
				move_array(drray, brray);
				}
			else
				{
				}
		//now do the restoration
		//try this to fix apparent scaling problem 29 may 05
		//this fixes scaling effects, but there is still a problem in the restoration.
		normalize_array_complex(brray);
		rescale_array_complex(brray, Float_13);//for testing
		
		deblur(crray, brray, array, snr);	 //nitty gritty
		}

	else if (Int_3==DEBLUR_WITH_WINDOWS)
		{
		//---11 restore a picture--------------------------------------------------------------
		#define PI 3.1415927
		float extent, angle, snr;
		extent=4.20;
		angle=PI*.62;
		snr=100.;
		psf_oof_complex(brray, extent);
		edge(array, brray, crray);
		psf_oof_complex(brray, extent);
		deblur1(crray, brray, array, snr);
		}

	else if (Int_3==SHOW_MAP_KERNEL)
		{
		//---12 display restoration kernel transform--------------------------------------------
		#define PI 3.1415927
		float extent, snr, angle,thick;
		extent=Float_1;
		angle=PI*Float_2/180.;
		snr=Float_3;
		thick=Float_4;
		psf_genthick(array, extent, angle,thick, Int_2, Float_7, Float_8);
		
		//check to see if want to convolve psf with disc to blur its edge
			if(Int_1 == 4)
				{
				make_fuzzy(array, brray, Float_7);
				move_array(brray, array);
				}
			else
				{
				}
		
		
		shift_array_complex(array);
		mapgen_complex(array,brray,snr);
		magnitude_phase_complex(array);
		convert_to_logthreshold_complex(array, .000001);
		move_array(array, brray);
		}

	else if (Int_3==DISPLAY_PSF)
		{
		//---13 display a psf
		#define PI 3.1415927
		float extent, angle, thick;
		extent=Float_1;
		angle=PI*Float_2/180.;
		thick=Float_4;
		psf_genthick(brray, extent, angle, thick, Int_2, Float_7, Float_8);
		if(Int_1 == 4)
			{
			make_fuzzy(brray, array, Float_7);
			move_array(array, brray);
			}
		else
			{
			}
		}
		
	else if (Int_3==DISPLAY_PATGEN_MAGNITUDE)
		{
		//---14 pattern removal filter display demonstration
		///patgen(array, brray, 4, 4, 1);		//generate filter
		pattern_removal(array, crray, brray, Int_1, Int_2, 1, Float_1);
		magnitude_phase_complex(crray);
		move_array(crray, brray);
		}


	else if (Int_3==PATTERN_REMOVAL)
		//---15 pattern removal
		{
		pattern_removal(array, crray, brray, Int_1, Int_2,1, Float_1);
		}

	else if (Int_3==FFT)
		{
		//---16 fft demo
		//forward_fft(array);
		//inverse_fft(array);
		//move_array(array, brray);
		//---16 fft demo
		//normalize_array_complex(array);
		//delta_array_complex(array);
		//rampup_array_complex(array);
		//printf("delta function in array \n");
		//print_array_complex(array);
		forward_fft(array);
		//printf("fft in array \n");
		//print_array_complex(array);
		inverse_fft(array);
		//printf("inverst fft in array \n");
		//print_array_complex(array);
		move_array(array, brray);
		}
		
	//#elif TEST_NYQUIST
	else if (Int_3==TEST_NYQUIST)
		{
		//---17 test the 'nyquist removal' routine
		double hilev, lolev;
		shift_array_complex(array);
		forward_fft(array);
		magnitude_phase_complex(array);
		convert_to_logthreshold_complex(array, .1);
		maxmin_array(array, &hilev, &lolev, 0);
		border_padding(array, 128, 128, -3);
		move_array(array, brray);
		}
		
	//#elif BANDPASS_FILTERING
	else if (Int_3==BANDPASS_FILTERING)
		{
		//---18 Try out some simple fourier filters
#if 0
		#define PI 3.1415927
		float extent, angle, thick, nstd, max_mod;
		extent=Float_1;
		angle=PI*Float_2/180.;
		thick=Float_4;
		nstd=Float_3;
#endif		
		forward_fft(array);
		quadrant_swap(array);
		//border_padding(array, Int_1, Int_2, 0.);
		punch_hole(array, Float_1, Float_2);
#if 0
		psf_invert_gaussellipse(array, extent, thick, angle, nstd, max_mod);
#endif
		quadrant_swap(array);
		inverse_fft(array);
		move_array(array, brray);
		}
	//#elif STRETCHING
	else if (Int_3==STRETCH)
		{
		//---19 Try out STRETCHING OPERATION
		stretch_complex_array(array, Float_1, Float_2, 0);
		move_array(array, brray);
		}

	//--------------------------------------------------------------------	
	//#else
	else 
		{
		move_array(array, brray);
		}
	//#endif
	
	//choose to invert or not invert output image, according to the value of Float_14
	
		if (Float_14 > 0.5)
			{
			invert_realpart_complex(brray);
			stretch_complex_array(brray, 255., 0., 0);
			}
		else
			{
			stretch_complex_array(brray, 255., 0., 0);
			}
			
	//clip to stretch interactively according to values set in Float_5 and Float_6
		maxmin_array(brray, &hilev, &lolev, 0);//find extrema
		v=hilev-lolev;
		clip_realpart_complex(brray, lolev+(Float_5*v),
		hilev-(Float_6*v));					//clip the top&bottom		

			
		if(l==0)
			complex_to_byte_array(brray, workArray1, 0);
		else if (l==1)
			complex_to_byte_array(brray, workArray2, 0);
		else 
			complex_to_byte_array(brray, workArray3, 0);
	}//closes numberColors loop

//alex||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||testing
	
		//if you wish to process color, the barrays below should be 
		//replaced in order by complex arrays based respectively on 
		//red, green and blue data

	copyToWindowArray(workArray1, workArray2, workArray3, processColors,
					  baseAddressImage, samplesPerPixel, height, width);
	
	// test
	//move_array(brray, drray);
	complex_to_byte_array(drray, workArray4, 0);
	if(baseAddressArray)
		copyToWindowArray(workArray4, workArray4, workArray4, 0,
						  baseAddressArray, samplesPerPixel, height, width);

done:
	//fftwnd_destroy_plan(p);
	free(workArray1);
	free(workArray2);
	free(workArray3);
	free(array->ptr);
	free(brray->ptr);
	free(crray->ptr);
	free(drray->ptr);
	free(array);
	free(brray);
	free(crray);
	free(drray);
}
