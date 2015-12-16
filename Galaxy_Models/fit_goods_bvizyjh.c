/*
3DCOLORFIT: Niraj Welikala (University of Oxford)
Last version: 07/12/2015
*/

/*
Model fit colour gradient measurement - works on catalog of objects in GOODS-S
*/
#include <pthread.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <unistd.h>
#include <fitsio.h>
#include <complex.h>
#include <fftw3.h>
#include <string.h>

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
//#include <gsl/gsl_complex.h>

#define pi M_PI
#define c 3e8 // ms-1
#define pi M_PI

#define Npix_side_acs 333
#define Npix_side_wfc3 168
#define pixelsize_arcsec_acs 0.03
#define pixelsize_arcsec_wfc3 0.06

#define zeropoint_filter1 25.673
#define zeropoint_filter2 26.486
#define zeropoint_filter3 25.654
#define zeropoint_filter4 24.862
#define zeropoint_filter5 26.2687
#define zeropoint_filter6 26.2303
#define zeropoint_filter7 25.9463

//Approximate exposure times in full depth stacked images
#define exposuretime_filter5 6000  // secs
#define exposuretime_filter6 6000 // secs
#define exposuretime_filter7 6000 //secs

#define sb_limit_bulge 0.00001404 //rmax=4.5re
#define sb_limit_disk 0.0111  //rmax=4.5h

/* Dust correction from Graham & Worely (2008) */
/* c_lambda */
#define l_lambda 1.5852e-09
#define m_lambda -5.24305e-05
#define n_lambda 1.45466

/* d_lambda */
#define u_lambda -2.75675e-09
#define v_lambda 3.88187e-05
#define w_lambda 0.23887

/* correlation noise estimated from 10 galaxies in file sample_corrfactor */
#define corr_correlations_filter1 0.642 //0.608274
#define corr_correlations_filter2 0.692 //0.637484
#define corr_correlations_filter3 0.643 //0.597494
#define corr_correlations_filter4 0.696 //0.647182
#define corr_correlations_filter5 0.850 // UPDATE
#define corr_correlations_filter6 0.857  //UPDATE
#define corr_correlations_filter7 0.839  //UPDATE


#define Nparams 5 // No. of free parameters in model
#define nfilters_acs 4 // No. of ACS bands
#define nfilters_wfc3 3 // No. of WFC3 bands

#define ntypes 6 // No. of morphology types of models
#define npoints_bulge 10000 // number of points in base bulge model that is read in
#define npoints_disk 10000 // number of points in base disk model that is read in
#define Nmaxiterations 2000 //2000 // No. of free parameters in model

// accuracy tolerance for GSL simplex routine
double simplex_tol = 1.0e-5;

//colour factor filter2:filter 1
int itype;
int npixels_acs, npixels_wfc3;

void getimagedim(char *imagename);
void readimage(char* filename, double *);


double *apix, *simdata_filter1, *simdata_filter2, *simdata_filter3, *simdata_filter4,*simdata_filter5, *simdata_filter6, *simdata_filter7, *psf_image_filter1, *psf_image_filter2, *psf_image_filter3, *psf_image_filter4,*psf_image_filter5, *psf_image_filter6, *psf_image_filter7, *simrms_filter1, *simrms_filter2, *simrms_filter3, *simrms_filter4, *simrms_filter5, *simrms_filter6, *simrms_filter7, *model_conv_acs, *model_conv_wfc3, *model_unconv_acs, *model_unconv_wfc3, *psf_image_acs, *psf_image_wfc3;
double *in_acs, *in_wfc3, *identity_acs, *identity_wfc3, *final_acs, *final_wfc3;
fftw_complex *inTrans_acs, *inTrans_wfc3, *identityTrans_acs, *identityTrans_wfc3, *FinalFFT_acs, *FinalFFT_wfc3;
fftw_plan plan1_acs, plan1_wfc3, plan2_acs, plan2_wfc3, plan3_acs, plan3_wfc3;

double *x_bulge,*y_bulge,*z_bulge, *I_xyz_bulge_filter1, *I_xyz_disk_filter1, *x_disk, *y_disk, *z_disk, f2f1_ratio_disk_type, f2f1_ratio_bulge_type;
double pixelsize_kpc_acs, pixelsize_kpc_wfc3, arcsecstokpc;
int ncols_var, nrows_var;
int ncols_psf_acs, ncols_psf_wfc3, nrows_psf_acs, nrows_psf_wfc3;
long ns;
int w, h;
double Nfreedom;
double amp_etype_type_filter1, amp_etype_type_filter2, amp_disk_type_filter1, amp_disk_type_filter2, amp_bulge_type_filter1, amp_bulge_type_filter2;
double amp_etype_type_filter3, amp_etype_type_filter4, amp_disk_type_filter3, amp_disk_type_filter4, amp_bulge_type_filter3, amp_bulge_type_filter4;
double amp_etype_type_filter5, amp_etype_type_filter6, amp_disk_type_filter5, amp_disk_type_filter6, amp_bulge_type_filter5, amp_bulge_type_filter6, amp_bulge_type_filter7, amp_disk_type_filter7;


/* Make unconvolved model on ACS pixel grid */
void make2Dgalaxymodelunconvolved_acs (int componenttype, int filtertype, double scaling, double e1, double e2, double xcentre, double ycentre, double arcsecstokpc, double weight_component)
{

  double c_lambda, d_lambda;

  if (filtertype==0) { c_lambda = (l_lambda * 4350.0 * 4350.0 ) + (m_lambda * 4350.0) + n_lambda  ;  d_lambda = (u_lambda * 4350.0 * 4350.0) +  (v_lambda * 4350.0) + w_lambda ;  }
  if (filtertype==1) { c_lambda = (l_lambda * 6060.0 * 6060.0 ) + (m_lambda * 6060.0) + n_lambda  ; d_lambda = (u_lambda * 6060.0 *6060.0) + (v_lambda * 6060.0) + w_lambda ; } //V
  if (filtertype==2) { c_lambda = (l_lambda * 7750.0 * 7750.0 ) + (m_lambda * 7750.0) + n_lambda  ; d_lambda = (u_lambda * 7750.0 *7750.0) + (v_lambda * 7750.0) + w_lambda ; } //I
  if (filtertype==3) { c_lambda = (l_lambda * 8500.0 * 8500.0 ) + (m_lambda * 8500.0) + n_lambda  ; d_lambda = (u_lambda * 8500.0 *8500.0) + (v_lambda * 8500.0) + w_lambda ; } //Z

  //double *model_unconv;
  double x_final_bulge, y_final_bulge, x_final_disk, y_final_disk, summodel;
  double disk_flux_tot = 0.0;
  double bulge_flux_tot = 0.0;
  double pixelsize_kpc, norm_disk, norm_bulge;

  /*
  printf("xcentre: %lf\n",xcentre);
  printf("ycentre: %lf\n",ycentre);
  printf("scaling: %lf\n", scaling);
  printf("e1: %lf\n", e1);
  printf("e2: %lf\n", e2);
  */

  int i_index=0, j_index=0;
  int i=0, j=0, k=0, ncols, nrows, pix, npixels,ii,jj,ifz=0;
  int nrows_tot, ncols_tot, ii_o, jj_o, pix_o, pix_tot,jj_o_psf, ii_o_psf, pix_o_psf;
  double xmin, xmax, ymin, ymax;
  double theta_rad, posangle_rad, filterfactor;

  filterfactor = 1.0;

  nrows_tot=Npix_side_acs; ncols_tot=Npix_side_acs;

  h = nrows_tot; w = nrows_tot;

  pixelsize_kpc = pixelsize_arcsec_acs * arcsecstokpc; // pixel size in kpc

/* Set boundaries xmin, xmax, ymin, ymax as N=3 times larger than xmin or xmax: CHANGE THIS IN NEXT ITERATION TO OUTER RADIUS CUT OFF */
  xmin = -(ncols_tot/2) * pixelsize_kpc;
  xmax = (ncols_tot/2) * pixelsize_kpc;
  ymin = -(nrows_tot/2) * pixelsize_kpc;
  ymax = (nrows_tot/2) * pixelsize_kpc;


  /*
  printf("Npix_side: %d\n",Npix_side);
  printf("pixelsize_arcsec arcsecstokpc: %lf %lf\n", pixelsize_arcsec, arcsecstokpc);
  printf("xmin: %lf\n",xmin);
  printf("xmax: %lf\n",xmax);
  printf("ymin: %lf\n",ymin);
  printf("ymax: %lf\n",ymax);
  */

  npixels = nrows_tot * ncols_tot;

  //printf("npixels: %d\n", npixels);

  theta_rad = asin(sqrt((e1*e1) + (e2*e2)));
  posangle_rad  =  0.5 * atan2(e2,e1) ;

  double theta_e1e2;

  //printf("theta_deg: %lf\n",theta_rad * (180.0/M_PI));
  //printf("posangle_deg: %lf\n", posangle_rad * (180.0/M_PI));

  double x_scaled, y_scaled, z_scaled, x_inclined, y_inclined;

  /* Initalize arrays  */

   for (pix=0; pix<npixels; pix++)
    {
      model_unconv_acs[pix] = 0.0; //unconvolved
      model_conv_acs[pix]  = 0.0; //final convolved
    }


  /* disk */
  if (componenttype==1) {

    //double *x_final_disk, *y_final_disk;
    disk_flux_tot = 0.0;
    norm_disk = 0.0;


    /* Scale, incline and rotate galaxy with these parameters */
	 for (ifz=0; ifz<npoints_disk;ifz++){

	   /* Scale x, y, z by same factors */

	   x_scaled= (c_lambda - (d_lambda*(2.5*log10(cos(theta_rad))))) * x_disk[ifz] * 5 * scaling;  //kpc
	   y_scaled = (c_lambda - (d_lambda*(2.5*log10(cos(theta_rad))))) *  y_disk[ifz] * 5 * scaling ;  //kpc
	   z_scaled = (c_lambda - (d_lambda*(2.5*log10(cos(theta_rad))))) * z_disk[ifz] * 5 * scaling;  //kpc


   /* Inclination: Rotate the object about x axis. This rotates axis of symmetry of galaxy from z to z' */
  /* look at projected image in plane of sky (perpendicular to line of sight) -> take x',y' */

	   x_inclined = x_scaled ;
	   y_inclined = (cos(theta_rad) * y_scaled) - (sin(theta_rad) * z_scaled) ;

    /* Now apply rotation about new z axis (line of sight) -- in plane of sky */
    /* Added a positional offset */
	   //printf("xcentre, ycentre: %lf %lf\n",xcentre,ycentre);
	   x_final_disk = ((cos(posangle_rad) * x_inclined) - ( sin(posangle_rad) * y_inclined  )) ;
	   y_final_disk = ((sin(posangle_rad) * x_inclined) + ( cos(posangle_rad) * y_inclined ))  ;


      /* find which pixel the particle falls into and do gridding of mass */

	   i_index = ( (int) (((x_final_disk+xcentre-xmin)/pixelsize_kpc) + 1) ) - 1; // pixel in x direction
	   j_index = ( (int) (((y_final_disk+ycentre-ymin)/pixelsize_kpc) + 1) ) - 1;  // pixel in y direction
	   pix = (ncols_tot *j_index) + i_index;

	   if (i_index>=0 && i_index<ncols_tot && j_index>=0 && j_index<nrows_tot)
	      {
	     	model_unconv_acs[pix] = model_unconv_acs[pix] + (I_xyz_disk_filter1[ifz]); // add the flux to the pixel
		disk_flux_tot = disk_flux_tot  + (I_xyz_disk_filter1[ifz]) ;
	      }


	 }

	 //printf("disk_flux_tot: %lf\n",disk_flux_tot);

  /* Normalize */

  for (pix=0; pix<npixels; pix++)
      {
	model_unconv_acs[pix] =  model_unconv_acs[pix]/disk_flux_tot;
      }

  }


 if (componenttype==0) {

    bulge_flux_tot=0.0;


/*re-initialize galaxy_image_filter1 and 2 */


	 /* Scale, incline and rotate galaxy with these parameters */
    for (ifz=0; ifz<npoints_bulge;ifz++){

	   /* Scale x, y, z by same factors */

	   x_scaled = x_bulge[ifz] * 5 * scaling ;  //from kpc
	   y_scaled = y_bulge[ifz] * 5 * scaling ;  //from kpc
	   z_scaled = z_bulge[ifz] * 5 * scaling ;  //from kpc


   /* Inclination: Rotate the object about x axis. This rotates axis of symmetry of galaxy from z to z' */
  /* look at projected image in plane of sky (perpendicular to line of sight) -> take x',y' */

	   x_inclined = x_scaled;
	   y_inclined = (cos(theta_rad) * y_scaled) - (sin(theta_rad) * z_scaled);

    /* Now apply random rotation about new z axis (line of sight) -- in plane of sky */

	   x_final_bulge = ((cos(posangle_rad) * x_inclined) - ( sin(posangle_rad) * y_inclined  ))  ;
	   y_final_bulge = ((sin(posangle_rad) * x_inclined) + ( cos(posangle_rad) * y_inclined ))  ;

	   i_index = ( (int) (((x_final_bulge+xcentre-xmin)/pixelsize_kpc) + 1) ) - 1; // pixel in x direction
	   j_index = ( (int) (((y_final_bulge+ycentre-ymin)/pixelsize_kpc) + 1) ) - 1;  // pixel in y direction
	   pix = (ncols_tot *j_index) + i_index;

	   if (i_index>=0 && i_index<ncols_tot && j_index>=0 && j_index<nrows_tot)
	      {
		model_unconv_acs[pix] = model_unconv_acs[pix] + (I_xyz_bulge_filter1[ifz]); // add the flux to the pixel
		bulge_flux_tot = bulge_flux_tot  + (I_xyz_bulge_filter1[ifz]) ;
	       }


	}


     for (pix=0;pix<npixels;pix++)
       {
	 model_unconv_acs[pix] = model_unconv_acs[pix]/bulge_flux_tot ;
       }

 }

}



/* Make unconvolved model on WFC3 pixel grid */
void make2Dgalaxymodelunconvolved_wfc3(int componenttype, int filtertype, double scaling, double e1, double e2, double xcentre, double ycentre, double arcsecstokpc, double weight_component)
{

  double c_lambda, d_lambda;

  if (filtertype==4) { c_lambda = (l_lambda * 10500.0 * 10500.0 ) + (m_lambda*10500.0) + n_lambda  ; d_lambda = (u_lambda*10500.0*10500.0) + (v_lambda*10500.0) + w_lambda ; } //Y
  if (filtertype==5) { c_lambda = (l_lambda * 12500.0 * 12500.0 ) + (m_lambda*12500.0) + n_lambda  ; d_lambda = (u_lambda*12500.0*12500.0) + (v_lambda*12500.0) + w_lambda ; } //J
  if (filtertype==6) { c_lambda = (l_lambda * 16500.0 * 16500.0 ) + (m_lambda*16500.0) + n_lambda  ; d_lambda = (u_lambda*16500.0*16500.0) + (v_lambda*16500.0) + w_lambda ; } //H

  double x_final_bulge, y_final_bulge, x_final_disk, y_final_disk, summodel;
  double disk_flux_tot = 0.0;
  double bulge_flux_tot = 0.0;
  double pixelsize_kpc, norm_disk, norm_bulge;

  /*
  printf("xcentre: %lf\n",xcentre);
  printf("ycentre: %lf\n",ycentre);
  printf("scaling: %lf\n", scaling);
  printf("e1: %lf\n", e1);
  printf("e2: %lf\n", e2);
  */

  int i_index=0, j_index=0;
  int i=0, j=0, k=0, ncols, nrows, pix, npixels,ii,jj,ifz=0;
  int nrows_tot, ncols_tot, ii_o, jj_o, pix_o, pix_tot,jj_o_psf, ii_o_psf, pix_o_psf;
  double xmin, xmax, ymin, ymax;
  double theta_rad, posangle_rad, filterfactor;

  filterfactor = 1.0;

  nrows_tot=Npix_side_wfc3; ncols_tot=Npix_side_wfc3;

  h = nrows_tot; w = nrows_tot;

  pixelsize_kpc = pixelsize_arcsec_wfc3 * arcsecstokpc; // pixel size in kpc

/* Set boundaries xmin, xmax, ymin, ymax as N=3 times larger than xmin or xmax: CHANGE THIS IN NEXT ITERATION TO OUTER RADIUS CUT OFF */
  xmin = -(ncols_tot/2) * pixelsize_kpc;
  xmax = (ncols_tot/2) * pixelsize_kpc;
  ymin = -(nrows_tot/2) * pixelsize_kpc;
  ymax = (nrows_tot/2) * pixelsize_kpc;


  /*
  printf("Npix_side: %d\n",Npix_side);
  printf("pixelsize_arcsec arcsecstokpc: %lf %lf\n", pixelsize_arcsec, arcsecstokpc);
  printf("xmin: %lf\n",xmin);
  printf("xmax: %lf\n",xmax);
  printf("ymin: %lf\n",ymin);
  printf("ymax: %lf\n",ymax);
  */

  npixels = nrows_tot * ncols_tot;

  //printf("npixels: %d\n", npixels);

  theta_rad = asin(sqrt((e1*e1) + (e2*e2)));
  posangle_rad  =  0.5 * atan2(e2,e1) ;

  double theta_e1e2;

  //printf("theta_deg: %lf\n",theta_rad * (180.0/M_PI));
  //printf("posangle_deg: %lf\n", posangle_rad * (180.0/M_PI));

  double x_scaled, y_scaled, z_scaled, x_inclined, y_inclined;

  /* Initalize arrays  */

   for (pix=0; pix<npixels; pix++)
    {
      model_unconv_wfc3[pix] = 0.0; //unconvolved
      model_conv_wfc3[pix]  = 0.0; //final convolved
    }


  /* disk */
  if (componenttype==1) {

    //double *x_final_disk, *y_final_disk;
    disk_flux_tot = 0.0;
    norm_disk = 0.0;


    /* Scale, incline and rotate galaxy with these parameters */
	 for (ifz=0; ifz<npoints_disk;ifz++){

	   /* Scale x, y, z by same factors */

	   x_scaled=  (c_lambda - (d_lambda*(2.5*log10(cos(theta_rad)))))  * x_disk[ifz] * 5 * scaling;  //kpc
	   y_scaled = (c_lambda - (d_lambda*(2.5*log10(cos(theta_rad))))) * y_disk[ifz] * 5 * scaling ;  //kpc
	   z_scaled =  (c_lambda - (d_lambda*(2.5*log10(cos(theta_rad))))) * z_disk[ifz] * 5 * scaling;  //kpc


   /* Inclination: Rotate the object about x axis. This rotates axis of symmetry of galaxy from z to z' */
  /* look at projected image in plane of sky (perpendicular to line of sight) -> take x',y' */

	   x_inclined = x_scaled ;
	   y_inclined = (cos(theta_rad) * y_scaled) - (sin(theta_rad) * z_scaled) ;

    /* Now apply rotation about new z axis (line of sight) -- in plane of sky */
    /* Added a positional offset */
	   //printf("xcentre, ycentre: %lf %lf\n",xcentre,ycentre);
	   x_final_disk = ((cos(posangle_rad) * x_inclined) - ( sin(posangle_rad) * y_inclined  )) ;
	   y_final_disk = ((sin(posangle_rad) * x_inclined) + ( cos(posangle_rad) * y_inclined ))  ;


      /* find which pixel the particle falls into and do gridding of mass */

	   i_index = ( (int) (((x_final_disk+xcentre-xmin)/pixelsize_kpc) + 1) ) - 1; // pixel in x direction
	   j_index = ( (int) (((y_final_disk+ycentre-ymin)/pixelsize_kpc) + 1) ) - 1;  // pixel in y direction
	   pix = (ncols_tot *j_index) + i_index;

	   if (i_index>=0 && i_index<ncols_tot && j_index>=0 && j_index<nrows_tot)
	      {
	     	model_unconv_wfc3[pix] = model_unconv_wfc3[pix] + (I_xyz_disk_filter1[ifz]); // add the flux to the pixel
		disk_flux_tot = disk_flux_tot  + (I_xyz_disk_filter1[ifz]) ;
	      }


     }

	 //printf("disk_flux_tot: %lf\n",disk_flux_tot);

  /* Normalize */

  for (pix=0; pix<npixels; pix++)
      {
	model_unconv_wfc3[pix] =  model_unconv_wfc3[pix]/disk_flux_tot;
      }

  }


 if (componenttype==0) {

    bulge_flux_tot=0.0;


/*re-initialize galaxy_image_filter1 and 2 */


	 /* Scale, incline and rotate galaxy with these parameters */
    for (ifz=0; ifz<npoints_bulge;ifz++){

	   /* Scale x, y, z by same factors */

	   x_scaled =  x_bulge[ifz] * 5 * scaling ;  //from kpc
	   y_scaled =  y_bulge[ifz] * 5 * scaling ;  //from kpc
	   z_scaled =  z_bulge[ifz] * 5 * scaling ;  //from kpc


   /* Inclination: Rotate the object about x axis. This rotates axis of symmetry of galaxy from z to z' */
  /* look at projected image in plane of sky (perpendicular to line of sight) -> take x',y' */

	   x_inclined = x_scaled;
	   y_inclined = (cos(theta_rad) * y_scaled) - (sin(theta_rad) * z_scaled);

    /* Now apply random rotation about new z axis (line of sight) -- in plane of sky */

	   x_final_bulge = ((cos(posangle_rad) * x_inclined) - ( sin(posangle_rad) * y_inclined  ))  ;
	   y_final_bulge = ((sin(posangle_rad) * x_inclined) + ( cos(posangle_rad) * y_inclined ))  ;

	   i_index = ( (int) (((x_final_bulge+xcentre-xmin)/pixelsize_kpc) + 1) ) - 1; // pixel in x direction
	   j_index = ( (int) (((y_final_bulge+ycentre-ymin)/pixelsize_kpc) + 1) ) - 1;  // pixel in y direction
	   pix = (ncols_tot *j_index) + i_index;

	   if (i_index>=0 && i_index<ncols_tot && j_index>=0 && j_index<nrows_tot)
	      {
		model_unconv_wfc3[pix] = model_unconv_wfc3[pix] + (I_xyz_bulge_filter1[ifz]); // add the flux to the pixel
		bulge_flux_tot = bulge_flux_tot  + (I_xyz_bulge_filter1[ifz]) ;
	       }


	}


    //printf("bulge_flux_tot: %lf\n",bulge_flux_tot);


     for (pix=0;pix<npixels;pix++)
       {
	 model_unconv_wfc3[pix] = model_unconv_wfc3[pix]/bulge_flux_tot ;
       }

 }

}


// make PSF convolved galaxy model at ACS image
void make2Dgalaxymodelconvolved_acs (int componenttype, int filtertype, double scaling, double e1, double e2, double xcentre, double ycentre, double arcsecstokpc, double weight_component)
{

   double c_lambda, d_lambda;

  if (filtertype==0) { c_lambda = (l_lambda * 4350.0 * 4350.0 ) + (m_lambda*4350.0) + n_lambda  ; d_lambda = (u_lambda*4350.0*4350.0) + (v_lambda*4350.0) + w_lambda ; } //B
  if (filtertype==1) { c_lambda = (l_lambda * 6060.0 * 6060.0 ) + (m_lambda*6060.0) + n_lambda  ; d_lambda = (u_lambda*6060.0*6060.0) + (v_lambda*6060.0) + w_lambda ; } //V
  if (filtertype==2) { c_lambda = (l_lambda * 7750.0 * 7750.0 ) + (m_lambda*7750.0) + n_lambda  ; d_lambda = (u_lambda*7750.0*7750.0) + (v_lambda*7750.0) + w_lambda ; } //I
  if (filtertype==3) { c_lambda = (l_lambda * 8500.0 * 8500.0 ) + (m_lambda*8500.0) + n_lambda  ; d_lambda = (u_lambda*8500.0*8500.0) + (v_lambda*8500.0) + w_lambda ; } //Z

  int pheight, pwidth, imheight, imwidth, psfsize,imsize, padfactor, padwidth, padheight, padsize, complexwidth, complexheight, complexsize,padcomplexsize;
  double x_final_bulge, y_final_bulge, x_final_disk, y_final_disk, summodel;
  double disk_flux_tot = 0.0;
  double bulge_flux_tot = 0.0;
  double norm_disk;
  double pixelsize_kpc;

  /*
  printf("xcentre: %lf\n",xcentre);
  printf("ycentre: %lf\n",ycentre);
  printf("scaling: %lf\n", scaling);
  printf("e1: %lf\n", e1);
  printf("e2: %lf\n", e2);
  printf("f2f1_ratio_bulge, f2f1_ratio_disk: %lf %lf\n",f2f1_ratio_bulge,f2f1_ratio_disk);
  */

  int i_index=0, j_index=0;
  int i=0, j=0, k=0, ncols, nrows, pix, npixels,ii,jj,ifz=0;
  int nrows_tot, ncols_tot, ii_o, jj_o, pix_o, pix_tot,jj_o_psf, ii_o_psf, pix_o_psf;
  double xmin, xmax, ymin, ymax;
  double theta_rad, posangle_rad, filterfactor;

  filterfactor = 1.0;

  nrows_tot=Npix_side_acs; ncols_tot=Npix_side_acs;

  h = nrows_tot; w = ncols_tot;

  pixelsize_kpc = pixelsize_arcsec_acs * arcsecstokpc; // pixel size in kpc

  xmin = -(ncols_tot/2) * pixelsize_kpc;
  xmax = (ncols_tot/2) * pixelsize_kpc;
  ymin = -(nrows_tot/2) * pixelsize_kpc;
  ymax = (nrows_tot/2) * pixelsize_kpc;

  npixels = nrows_tot * ncols_tot;
  //printf("h,w, npixels: %d %d %d\n",h,w,npixels);

  theta_rad = asin(sqrt((e1*e1) + (e2*e2)));
  posangle_rad  =  0.5 * atan2(e2,e1) ;

  double theta_e1e2;

  double x_scaled, y_scaled, z_scaled, x_inclined, y_inclined;

   /* Initalize arrays  */

   for (pix=0; pix<npixels; pix++)
    {
      model_unconv_acs[pix] = 0.0; //unconvolved
      model_conv_acs[pix]  = 0.0; //final convolved
    }



  /* disk */
  if (componenttype==1) {

    //double *x_final_disk, *y_final_disk;
    disk_flux_tot = 0.0;
    norm_disk = 0.0;


    /* Scale, incline and rotate galaxy with these parameters */
	 for (ifz=0; ifz<npoints_disk;ifz++){

	   /* Scale x, y, z by same factors */

	   x_scaled= (c_lambda - (d_lambda*(2.5*log10(cos(theta_rad))))) * x_disk[ifz] * 5 * scaling;  //kpc
	   y_scaled = (c_lambda - (d_lambda*(2.5*log10(cos(theta_rad))))) * y_disk[ifz] * 5 * scaling ;  //kpc
	   z_scaled = (c_lambda - (d_lambda*(2.5*log10(cos(theta_rad))))) *z_disk[ifz] * 5 * scaling;  //kpc


   /* Inclination: Rotate the object about x axis. This rotates axis of symmetry of galaxy from z to z' */
  /* look at projected image in plane of sky (perpendicular to line of sight) -> take x',y' */

	   x_inclined = x_scaled ;
	   y_inclined = (cos(theta_rad) * y_scaled) - (sin(theta_rad) * z_scaled) ;

    /* Now apply rotation about new z axis (line of sight) -- in plane of sky */
    /* Added a positional offset */
	   //printf("xcentre, ycentre: %lf %lf\n",xcentre,ycentre);
	   x_final_disk = ((cos(posangle_rad) * x_inclined) - ( sin(posangle_rad) * y_inclined  )) ;
	   y_final_disk = ((sin(posangle_rad) * x_inclined) + ( cos(posangle_rad) * y_inclined ))  ;


      /* find which pixel the particle falls into and do gridding of mass */

	   i_index = ( (int) (((x_final_disk+xcentre-xmin)/pixelsize_kpc) + 1) ) - 1; // pixel in x direction
	   j_index = ( (int) (((y_final_disk+ycentre-ymin)/pixelsize_kpc) + 1) ) - 1;  // pixel in y direction
	   pix = (ncols_tot *j_index) + i_index;

	   if (i_index>=0 && i_index<ncols_tot && j_index>=0 && j_index<nrows_tot)
	      {
		model_unconv_acs[pix] = model_unconv_acs[pix] + (I_xyz_disk_filter1[ifz]); // add the flux to the pixel
		disk_flux_tot = disk_flux_tot  + (I_xyz_disk_filter1[ifz]) ;
	      }


     }

	 //printf("disk_flux_tot: %lf\n",disk_flux_tot);

  /* Normalize */

    for (pix=0; pix<npixels; pix++)
      {
	model_unconv_acs[pix] =  model_unconv_acs[pix]/disk_flux_tot;
      }


  }


 if (componenttype==0) {

    bulge_flux_tot=0.0;


/*re-initialize galaxy_image_filter1 and 2 */


	 /* Scale, incline and rotate galaxy with these parameters */
    for (ifz=0; ifz<npoints_bulge;ifz++){

	   /* Scale x, y, z by same factors */

	   x_scaled =  x_bulge[ifz] * 5 * scaling ;  //from kpc
	   y_scaled =  y_bulge[ifz] * 5 * scaling ;  //from kpc
	   z_scaled =  z_bulge[ifz] * 5 * scaling ;  //from kpc


   /* Inclination: Rotate the object about x axis. This rotates axis of symmetry of galaxy from z to z' */
  /* look at projected image in plane of sky (perpendicular to line of sight) -> take x',y' */

	   x_inclined = x_scaled;
	   y_inclined = (cos(theta_rad) * y_scaled) - (sin(theta_rad) * z_scaled);

    /* Now apply random rotation about new z axis (line of sight) -- in plane of sky */

	   x_final_bulge = ((cos(posangle_rad) * x_inclined) - ( sin(posangle_rad) * y_inclined  ))  ;
	   y_final_bulge = ((sin(posangle_rad) * x_inclined) + ( cos(posangle_rad) * y_inclined ))  ;

	   i_index = ( (int) (((x_final_bulge+xcentre-xmin)/pixelsize_kpc) + 1) ) - 1; // pixel in x direction
	   j_index = ( (int) (((y_final_bulge+ycentre-ymin)/pixelsize_kpc) + 1) ) - 1;  // pixel in y direction
	   pix = (ncols_tot *j_index) + i_index;

	   if (i_index>=0 && i_index<ncols_tot && j_index>=0 && j_index<nrows_tot)
	      {
		 model_unconv_acs[pix] = model_unconv_acs[pix] + (I_xyz_bulge_filter1[ifz]); // add the flux to the pixel
		 bulge_flux_tot = bulge_flux_tot  + (I_xyz_bulge_filter1[ifz]) ;
	       }


	}


    //printf("bulge_flux_tot: %lf\n",bulge_flux_tot);

   for (pix=0;pix<npixels;pix++)
       {
	 model_unconv_acs[pix] = model_unconv_acs[pix]/bulge_flux_tot ;
       }

 }

   /* -------------------------- */
   /* Do PSF convolution */
   /* -------------------------- */

   /* get the PSF into array */
  i=0;
  for (i=0;i<h*w;i++)
     {
       if (filtertype==0) {psf_image_acs[i]=psf_image_filter1[i] ;}
       if (filtertype==1) {psf_image_acs[i]=psf_image_filter2[i] ;}
       if (filtertype==2) {psf_image_acs[i]=psf_image_filter3[i] ;}
       if (filtertype==3) {psf_image_acs[i]=psf_image_filter4[i] ;}
      }


  /* populate real arrays  */

  for (i=0; i<h*w; i++)
   {
     identity_acs[i] = psf_image_acs[i];
     in_acs[i] = model_unconv_acs[i];
     //fprintf(testreal,"i identity in %d %lf %lf\n",i,identity[i],in[i]);
     //printf("i identity in %d %lf %lf\n",i,identity[i],in[i]);
     //if (in[i]!=0.0) {printf("zero!");}
   }


  //printf("populated real arrays\n");

   //Execute them

  fftw_execute(plan1_acs);
  fftw_execute(plan2_acs);

  int iiii=0;
  for(iiii=0; iiii<(w/2+1)*h; iiii++)
    {

      //fprintf(testfft,"ii inTrans idensityTrans %d %lf %lf\n",ii,inTrans[ii],identityTrans[ii]);
      // printf("ii inTrans idensityTrans %d %lf %lf\n",ii,inTrans[ii],identityTrans[ii]);
      FinalFFT_acs[iiii] = inTrans_acs[iiii] * identityTrans_acs[iiii] ;
      //printf("ii inTrans idensityTrans FinalFFT %d %lf %lf %lf\n",iiii,inTrans[iiii],identityTrans[iiii],FinalFFT[iiii]);
    }

  fftw_execute(plan3_acs);



 i=0;j=0;ii=0;jj=0;

  /* Swap quandrants */

 for (i=0; i<h; i++)
   {
     ii = i + h/2;
     if (ii>=h) ii -= h;
     for (j=0; j<w; j++)
       {
	 jj = j + w/2;
	 if (jj>=w) jj -= w;
	 model_conv_acs[jj + w*ii] = final_acs[j + w*i];
       }
   }

 double tot_flux = 0.0;

 /*Normalize */

 for (i=0;i<h*w;i++)
   tot_flux = tot_flux + model_conv_acs[i];

 for (i=0;i<h*w;i++)
   model_conv_acs[i] = model_conv_acs[i] * filterfactor * weight_component/tot_flux ;

}





// make PSF convolved galaxy model at WFC3 image
void make2Dgalaxymodelconvolved_wfc3 (int componenttype, int filtertype, double scaling, double e1, double e2, double xcentre, double ycentre, double arcsecstokpc, double weight_component)
{
   double c_lambda, d_lambda;

  if (filtertype==4) { c_lambda = (l_lambda * 10500.0 * 10500.0 ) + (m_lambda*10500.0) + n_lambda  ; d_lambda = (u_lambda*10500.0*10500.0) + (v_lambda*10500.0) + w_lambda ; } //Y
  if (filtertype==5) { c_lambda = (l_lambda * 12500.0 * 12500.0 ) + (m_lambda*12500.0) + n_lambda  ; d_lambda = (u_lambda*12500.0*12500.0) + (v_lambda*12500.0) + w_lambda ; } //J
  if (filtertype==6) { c_lambda = (l_lambda * 16500.0 * 16500.0 ) + (m_lambda*16500.0) + n_lambda  ; d_lambda = (u_lambda*16500.0*16500.0) + (v_lambda*16500.0) + w_lambda ; } //H


  int pheight, pwidth, imheight, imwidth, psfsize,imsize, padfactor, padwidth, padheight, padsize, complexwidth, complexheight, complexsize,padcomplexsize;
  double x_final_bulge, y_final_bulge, x_final_disk, y_final_disk, summodel;
  double disk_flux_tot = 0.0;
  double bulge_flux_tot = 0.0;
  double norm_disk;
  double pixelsize_kpc;

  /*
  printf("xcentre: %lf\n",xcentre);
  printf("ycentre: %lf\n",ycentre);
  printf("scaling: %lf\n", scaling);
  printf("e1: %lf\n", e1);
  printf("e2: %lf\n", e2);
  printf("f2f1_ratio_bulge, f2f1_ratio_disk: %lf %lf\n",f2f1_ratio_bulge,f2f1_ratio_disk);
  */

  int i_index=0, j_index=0;
  int i=0, j=0, k=0, ncols, nrows, pix, npixels,ii,jj,ifz=0;
  int nrows_tot, ncols_tot, ii_o, jj_o, pix_o, pix_tot,jj_o_psf, ii_o_psf, pix_o_psf;
  double xmin, xmax, ymin, ymax;
  double theta_rad, posangle_rad, filterfactor;

  filterfactor = 1.0;

  nrows_tot=Npix_side_wfc3; ncols_tot=Npix_side_wfc3;

  h = nrows_tot; w = ncols_tot;

  pixelsize_kpc = pixelsize_arcsec_wfc3 * arcsecstokpc; // pixel size in kpc

  xmin = -(ncols_tot/2) * pixelsize_kpc;
  xmax = (ncols_tot/2) * pixelsize_kpc;
  ymin = -(nrows_tot/2) * pixelsize_kpc;
  ymax = (nrows_tot/2) * pixelsize_kpc;

  npixels = nrows_tot * ncols_tot;
  //printf("h,w, npixels: %d %d %d\n",h,w,npixels);

  theta_rad = asin(sqrt((e1*e1) + (e2*e2)));
  posangle_rad  =  0.5 * atan2(e2,e1) ;

  double theta_e1e2;

  double x_scaled, y_scaled, z_scaled, x_inclined, y_inclined;

   /* Initalize arrays  */

   for (pix=0; pix<npixels; pix++)
    {
      model_unconv_wfc3[pix] = 0.0; //unconvolved
      model_conv_wfc3[pix]  = 0.0; //final convolved
    }



  /* disk */
  if (componenttype==1) {

    //double *x_final_disk, *y_final_disk;
    disk_flux_tot = 0.0;
    norm_disk = 0.0;


    /* Scale, incline and rotate galaxy with these parameters */
	 for (ifz=0; ifz<npoints_disk;ifz++){

	   /* Scale x, y, z by same factors */

	   x_scaled= (c_lambda - (d_lambda*(2.5*log10(cos(theta_rad)))))  * x_disk[ifz] * 5 * scaling;  //kpc
	   y_scaled = (c_lambda - (d_lambda*(2.5*log10(cos(theta_rad)))))  * y_disk[ifz] * 5 * scaling ;  //kpc
	   z_scaled = (c_lambda - (d_lambda*(2.5*log10(cos(theta_rad)))))  * z_disk[ifz] * 5 * scaling;  //kpc


   /* Inclination: Rotate the object about x axis. This rotates axis of symmetry of galaxy from z to z' */
  /* look at projected image in plane of sky (perpendicular to line of sight) -> take x',y' */

	   x_inclined = x_scaled ;
	   y_inclined = (cos(theta_rad) * y_scaled) - (sin(theta_rad) * z_scaled) ;

    /* Now apply rotation about new z axis (line of sight) -- in plane of sky */
    /* Added a positional offset */
	   //printf("xcentre, ycentre: %lf %lf\n",xcentre,ycentre);
	   x_final_disk = ((cos(posangle_rad) * x_inclined) - ( sin(posangle_rad) * y_inclined  )) ;
	   y_final_disk = ((sin(posangle_rad) * x_inclined) + ( cos(posangle_rad) * y_inclined ))  ;


      /* find which pixel the particle falls into and do gridding of mass */

	   i_index = ( (int) (((x_final_disk+xcentre-xmin)/pixelsize_kpc) + 1) ) - 1; // pixel in x direction
	   j_index = ( (int) (((y_final_disk+ycentre-ymin)/pixelsize_kpc) + 1) ) - 1;  // pixel in y direction
	   pix = (ncols_tot *j_index) + i_index;

	   if (i_index>=0 && i_index<ncols_tot && j_index>=0 && j_index<nrows_tot)
	      {
		model_unconv_wfc3[pix] = model_unconv_wfc3[pix] + (I_xyz_disk_filter1[ifz]); // add the flux to the pixel
		disk_flux_tot = disk_flux_tot  + (I_xyz_disk_filter1[ifz]) ;
	      }


     }

	 //printf("disk_flux_tot: %lf\n",disk_flux_tot);

  /* Normalize */

    for (pix=0; pix<npixels; pix++)
      {
	model_unconv_wfc3[pix] =  model_unconv_wfc3[pix]/disk_flux_tot;
      }


  }


 if (componenttype==0) {

    bulge_flux_tot=0.0;


/*re-initialize galaxy_image_filter1 and 2 */


	 /* Scale, incline and rotate galaxy with these parameters */
    for (ifz=0; ifz<npoints_bulge;ifz++){

	   /* Scale x, y, z by same factors */

	   x_scaled = x_bulge[ifz] * 5 * scaling ;  //from kpc
	   y_scaled = y_bulge[ifz] * 5 * scaling ;  //from kpc
	   z_scaled = z_bulge[ifz] * 5 * scaling ;  //from kpc


   /* Inclination: Rotate the object about x axis. This rotates axis of symmetry of galaxy from z to z' */
  /* look at projected image in plane of sky (perpendicular to line of sight) -> take x',y' */

	   x_inclined = x_scaled;
	   y_inclined = (cos(theta_rad) * y_scaled) - (sin(theta_rad) * z_scaled);

    /* Now apply random rotation about new z axis (line of sight) -- in plane of sky */

	   x_final_bulge = ((cos(posangle_rad) * x_inclined) - ( sin(posangle_rad) * y_inclined  ))  ;
	   y_final_bulge = ((sin(posangle_rad) * x_inclined) + ( cos(posangle_rad) * y_inclined ))  ;

	   i_index = ( (int) (((x_final_bulge+xcentre-xmin)/pixelsize_kpc) + 1) ) - 1; // pixel in x direction
	   j_index = ( (int) (((y_final_bulge+ycentre-ymin)/pixelsize_kpc) + 1) ) - 1;  // pixel in y direction
	   pix = (ncols_tot *j_index) + i_index;

	   if (i_index>=0 && i_index<ncols_tot && j_index>=0 && j_index<nrows_tot)
	      {
		 model_unconv_wfc3[pix] = model_unconv_wfc3[pix] + (I_xyz_bulge_filter1[ifz]); // add the flux to the pixel
		 bulge_flux_tot = bulge_flux_tot  + (I_xyz_bulge_filter1[ifz]) ;
	       }


	}


    //printf("bulge_flux_tot: %lf\n",bulge_flux_tot);

   for (pix=0;pix<npixels;pix++)
       {
	 model_unconv_wfc3[pix] = model_unconv_wfc3[pix]/bulge_flux_tot ;
       }

 }

   /* -------------------------- */
   /* Do PSF convolution */
   /* -------------------------- */

   /* get the PSF into array */
  i=0;
  for (i=0;i<h*w;i++)
     {
       if (filtertype==4) {psf_image_wfc3[i]=psf_image_filter5[i] ;}
       if (filtertype==5) {psf_image_wfc3[i]=psf_image_filter6[i] ;}
       if (filtertype==6) {psf_image_wfc3[i]=psf_image_filter7[i] ;}
      }


  /* populate real arrays  */

  for (i=0; i<h*w; i++)
   {
     identity_wfc3[i] = psf_image_wfc3[i];
     in_wfc3[i] = model_unconv_wfc3[i];
     //fprintf(testreal,"i identity in %d %lf %lf\n",i,identity[i],in[i]);
     //printf("i identity in %d %lf %lf\n",i,identity[i],in[i]);
     //if (in[i]!=0.0) {printf("zero!");}
   }


  //printf("populated real arrays\n");

   //Execute them

  fftw_execute(plan1_wfc3);
  fftw_execute(plan2_wfc3);

  int iiii=0;
  for(iiii=0; iiii<(w/2+1)*h; iiii++)
    {

      //fprintf(testfft,"ii inTrans idensityTrans %d %lf %lf\n",ii,inTrans[ii],identityTrans[ii]);
      // printf("ii inTrans idensityTrans %d %lf %lf\n",ii,inTrans[ii],identityTrans[ii]);
      FinalFFT_wfc3[iiii] = inTrans_wfc3[iiii] * identityTrans_wfc3[iiii] ;
      //printf("ii inTrans idensityTrans FinalFFT %d %lf %lf %lf\n",iiii,inTrans[iiii],identityTrans[iiii],FinalFFT[iiii]);
    }

  fftw_execute(plan3_wfc3);



 i=0;j=0;ii=0;jj=0;

  /* Swap quandrants */

 for (i=0; i<h; i++)
   {
     ii = i + h/2;
     if (ii>=h) ii -= h;
     for (j=0; j<w; j++)
       {
	 jj = j + w/2;
	 if (jj>=w) jj -= w;
	 model_conv_wfc3[jj + w*ii] = final_wfc3[j + w*i];
       }
   }

 double tot_flux = 0.0;

 /*Normalize */

 for (i=0;i<h*w;i++)
   tot_flux = tot_flux + model_conv_wfc3[i];

 //printf("tot_flux: %lf\n",tot_flux);

 for (i=0;i<h*w;i++)
   {
     model_conv_wfc3[i] = model_conv_wfc3[i] * filterfactor * weight_component/tot_flux ;
   }

}





// evaluate function to measure likelihood of data given the model, assuming the
// maximum-likelihood value for the normalisation
double jointlikel(gsl_vector *v, void *params)
{


  //printf("Entered jointlikel () 1\n");
  int j, ns, comp;
  int ll=0;
  double sumx, sumy, sumxx, sumyy, sumzz, sumxy, sumxz, sumzy, denom, lval;
  double  amp_bulge_filter1, amp_bulge_filter2, amp_bulge_filter3, amp_bulge_filter4, amp_bulge_filter5, amp_bulge_filter6, amp_bulge_filter7, amp_disk_filter1, amp_disk_filter2, amp_disk_filter3, amp_disk_filter4, amp_disk_filter5, amp_disk_filter6, amp_disk_filter7;
  double amp_etype_filter1, amp_etype_filter2, amp_etype_filter3, amp_etype_filter4, amp_etype_filter5, amp_etype_filter6, amp_etype_filter7, penalty_amp_bulge_filter1=0.0, penalty_amp_bulge_filter2=0.0, penalty_amp_bulge_filter3=0.0, penalty_amp_bulge_filter4=0.0, penalty_amp_bulge_filter5=0.0, penalty_amp_bulge_filter6=0.0, penalty_amp_bulge_filter7=0.0, penalty_amp_disk_filter1=0.0, penalty_amp_disk_filter2=0.0, penalty_amp_disk_filter3=0.0, penalty_amp_disk_filter4=0.0, penalty_amp_disk_filter5=0.0, penalty_amp_disk_filter6=0.0, penalty_amp_disk_filter7=0.0, penalty_amp_etype_filter1=0.0, penalty_amp_etype_filter2=0.0, penalty_amp_etype_filter3=0.0, penalty_amp_etype_filter4=0.0, penalty_amp_etype_filter5=0.0, penalty_amp_etype_filter6=0.0, penalty_amp_etype_filter7=0.0;

  double xcentre, ycentre, scaling, e1, e2;
  double theta_rad_fit, posangle_rad_fit, theta_deg_fit, posangle_deg_fit;
  double *data_filter1, *model_disk_filter1, *model_bulge_filter1;
  double *data_filter2, *model_disk_filter2, *model_bulge_filter2;
  double *data_filter3, *model_disk_filter3, *model_bulge_filter3;
  double *data_filter4, *model_disk_filter4, *model_bulge_filter4;
  double *data_filter5, *model_disk_filter5, *model_bulge_filter5;
  double *data_filter6, *model_disk_filter6, *model_bulge_filter6;
  double *data_filter7, *model_disk_filter7, *model_bulge_filter7;

  denom=0.0;

  data_filter1 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
  model_bulge_filter1=(double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
  model_disk_filter1 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));

  data_filter2 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
  model_bulge_filter2=(double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
  model_disk_filter2 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));

  data_filter3 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
  model_bulge_filter3=(double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
  model_disk_filter3 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));

  data_filter4 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
  model_bulge_filter4=(double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
  model_disk_filter4 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));

  data_filter5 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
  model_bulge_filter5=(double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
  model_disk_filter5 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));

  data_filter6 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
  model_bulge_filter6=(double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
  model_disk_filter6 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));

  data_filter7 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
  model_bulge_filter7=(double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
  model_disk_filter7 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));


  int *p = (int*)params;
  comp = *p;

  scaling = gsl_vector_get(v, 0);
  e1 = gsl_vector_get(v, 1);
  e2 = gsl_vector_get(v, 2);
  xcentre = gsl_vector_get(v, 3);
  ycentre = gsl_vector_get(v, 4);

  /* parameter ranges */

 /* Apply penalty if parameter values are outside allowed ranges */

  double penalty_scaling = 0;
  double penalty_e1 = 0;
  double penalty_e2 =0;
  double penalty_e1e2quad =0;
  double e1_diff, e2_diff, e1e2quad_diff;
  double theta_e1e2;
  double theta_deg, posangle_deg;
  double penalty_xcentre = 0.0;
  double penalty_ycentre = 0.0;
  double penalty_bulge_total=0.0;
  double penalty_amp_bulge=0.0;
  double penalty_amp_disk=0.0;
  double penalty_amp_etype=0.0;


  if (sqrt(e2*e2 + e1*e1) > 1.0) {penalty_e1e2quad = pow(npixels_acs-Nparams-1,4) * pow(sqrt(e2*e2 + e1*e1)-1,2) ; theta_e1e2 = atan2(e2,e1); e1 = 0.99 * cos(theta_e1e2); e2 = 0.99 * sin(theta_e1e2) ;}
  if (scaling<=0.) { penalty_scaling = pow(npixels_acs-Nparams-1,4) * (scaling) * (scaling) ; }
  //if (abs(xcentre)>= Npix_side/2.0) {penalty_xcentre =pow(ns-Nparams,6) * pow(xcentre - (Npix_side/2.0),2); }
  //if (abs(ycentre)>= Npix_side/2.0) {penalty_ycentre =pow(ns-Nparams,6) * pow(ycentre - (Npix_side/2.0),2); }

  /* copy contents from global variable to local ones before they are replaced by next call to make*/

  make2Dgalaxymodelconvolved_acs(0,0, scaling, e1, e2, xcentre, ycentre, arcsecstokpc,1.0);

   for (ll=0;ll<(npixels_acs);ll++)
     model_bulge_filter1[ll]= model_conv_acs[ll];

   make2Dgalaxymodelconvolved_acs(1,0, scaling, e1, e2, xcentre, ycentre, arcsecstokpc,1.0);

   for (ll=0;ll<(npixels_acs);ll++)
     model_disk_filter1[ll]=model_conv_acs[ll];

   make2Dgalaxymodelconvolved_acs(0,1, scaling, e1, e2, xcentre, ycentre,arcsecstokpc,1.0);

   for (ll=0;ll<(npixels_acs);ll++)
     model_bulge_filter2[ll]= model_conv_acs[ll];

   make2Dgalaxymodelconvolved_acs(1,1, scaling, e1, e2, xcentre, ycentre,  arcsecstokpc,1.0);

   for (ll=0;ll<(npixels_acs);ll++)
     model_disk_filter2[ll]=model_conv_acs[ll];

    make2Dgalaxymodelconvolved_acs(0,2, scaling, e1, e2, xcentre, ycentre,arcsecstokpc,1.0);

   for (ll=0;ll<(npixels_acs);ll++)
     model_bulge_filter3[ll]= model_conv_acs[ll];

   make2Dgalaxymodelconvolved_acs(1,2, scaling, e1, e2, xcentre, ycentre,  arcsecstokpc,1.0);

   for (ll=0;ll<(npixels_acs);ll++)
     model_disk_filter3[ll]=model_conv_acs[ll];

   make2Dgalaxymodelconvolved_acs(0,3, scaling, e1, e2, xcentre, ycentre, arcsecstokpc,1.0);

   for (ll=0;ll<(npixels_acs);ll++)
     model_bulge_filter4[ll]= model_conv_acs[ll];

   make2Dgalaxymodelconvolved_acs(1,3, scaling, e1, e2, xcentre, ycentre, arcsecstokpc,1.0);

   for (ll=0;ll<(npixels_acs);ll++)
     model_disk_filter4[ll]=model_conv_acs[ll];

   make2Dgalaxymodelconvolved_wfc3(0,4, scaling, e1, e2, xcentre, ycentre, arcsecstokpc,1.0);

   for (ll=0;ll<(npixels_wfc3);ll++)
     model_bulge_filter5[ll]= model_conv_wfc3[ll];

   make2Dgalaxymodelconvolved_wfc3(1,4, scaling, e1, e2, xcentre, ycentre,  arcsecstokpc,1.0);

   for (ll=0;ll<(npixels_wfc3);ll++)
     model_disk_filter5[ll]=model_conv_wfc3[ll];


   make2Dgalaxymodelconvolved_wfc3(0,5, scaling, e1, e2, xcentre, ycentre, arcsecstokpc,1.0);

   for (ll=0;ll<(npixels_wfc3);ll++)
     model_bulge_filter6[ll]= model_conv_wfc3[ll];

   make2Dgalaxymodelconvolved_wfc3(1,5, scaling, e1, e2, xcentre, ycentre,  arcsecstokpc,1.0);

   for (ll=0;ll<(npixels_wfc3);ll++)
     model_disk_filter6[ll]=model_conv_wfc3[ll];

   make2Dgalaxymodelconvolved_wfc3(0,6, scaling, e1, e2, xcentre, ycentre, arcsecstokpc,1.0);

   for (ll=0;ll<(npixels_wfc3);ll++)
     model_bulge_filter7[ll]= model_conv_wfc3[ll];

   make2Dgalaxymodelconvolved_wfc3(1,6, scaling, e1, e2, xcentre, ycentre, arcsecstokpc,1.0);

   for (ll=0;ll<(npixels_wfc3);ll++)
     model_disk_filter7[ll]=model_conv_wfc3[ll];



  for (j=0; j<npixels_acs; j++)
   {
     model_disk_filter1[j]= model_disk_filter1[j]/simrms_filter1[j];
     model_bulge_filter1[j]= model_bulge_filter1[j]/simrms_filter1[j];
     data_filter1[j] = simdata_filter1[j]/simrms_filter1[j];
     model_disk_filter2[j]= model_disk_filter2[j]/simrms_filter2[j];
     model_bulge_filter2[j]= model_bulge_filter2[j]/simrms_filter2[j];
     data_filter2[j] = simdata_filter2[j]/simrms_filter2[j];
     model_disk_filter3[j]= model_disk_filter3[j]/simrms_filter3[j];
     model_bulge_filter3[j]= model_bulge_filter3[j]/simrms_filter3[j];
     data_filter3[j] = simdata_filter3[j]/simrms_filter3[j];
     model_disk_filter4[j]= model_disk_filter4[j]/simrms_filter4[j];
     model_bulge_filter4[j]= model_bulge_filter4[j]/simrms_filter4[j];
     data_filter4[j] = simdata_filter4[j]/simrms_filter4[j];

     //printf("model_disk_filter1[j], model_bulge_filter1[j] %lf %lf \n",model_disk_filter1[j],model_bulge_filter1[j]);
   }

   for (j=0; j<npixels_wfc3; j++)
   {
     model_disk_filter5[j]= model_disk_filter5[j]/simrms_filter5[j];
     model_bulge_filter5[j]= model_bulge_filter5[j]/simrms_filter5[j];
     data_filter5[j] = simdata_filter5[j]/simrms_filter5[j];
     model_disk_filter6[j]= model_disk_filter6[j]/simrms_filter6[j];
     model_bulge_filter6[j]= model_bulge_filter6[j]/simrms_filter6[j];
     data_filter6[j] = simdata_filter6[j]/simrms_filter6[j];
     model_disk_filter7[j]= model_disk_filter7[j]/simrms_filter7[j];
     model_bulge_filter7[j]= model_bulge_filter7[j]/simrms_filter7[j];
     data_filter7[j] = simdata_filter7[j]/simrms_filter7[j];

   }


 /* Compute likelihood for E-type separately - bulge only */

    if (itype==0)
      {

	lval = 0.0;
	sumx=sumy=sumxx=sumyy=sumxy=0.0;

	for (j=0; j<npixels_acs; j++)
	  {

	    sumx += model_bulge_filter1[j];
	    sumy += data_filter1[j];
	    sumxx += pow(model_bulge_filter1[j],2);
	    sumyy += pow(data_filter1[j],2);
	    sumxy += data_filter1[j]*model_bulge_filter1[j];
	  }

  	denom = sumxx;


	if (denom != 0.)
	  {
	    amp_etype_filter1 = sumxy/sumxx;
	    amp_bulge_filter1=amp_etype_filter1;

	    for (j=0;j<npixels_acs;j++)
	      {
		lval += pow((data_filter1[j]- (amp_etype_filter1 * model_bulge_filter1[j])),2);

	      }

	    if (amp_etype_filter1<0.0) { penalty_amp_etype_filter1 = 1e7;}

	  }
	else
	  {
	    lval += 9e99;
	  }


	/*2nd filter*/
	sumx=0.0;
	sumy=0.0;
	sumxx=0.0;
	sumyy=0.0;
	sumxy=0.0;

	for (j=0; j<npixels_acs; j++)
	  {

	    sumx += model_bulge_filter2[j];
	    sumy += data_filter2[j];
	    sumxx += pow(model_bulge_filter2[j],2);
	    sumyy += pow(data_filter2[j],2);
	    sumxy += data_filter2[j]*model_bulge_filter2[j];
	  }

  	denom = sumxx;


	if (denom != 0.)
	  {
	    amp_etype_filter2=sumxy/sumxx;
	    amp_bulge_filter2=amp_etype_filter2;

	    for (j=0;j<npixels_acs;j++)
	      {
		lval += pow((data_filter2[j]- (amp_etype_filter2 * model_bulge_filter2[j])),2);

	      }

	     if (amp_etype_filter2<0.0) { penalty_amp_etype_filter2 = 1e7;}

	  }
	else
	  {
	    lval += 9e99;
	  }


	/* 3rd filter */
	sumx=0.0;
	sumy=0.0;
	sumxx=0.0;
	sumyy=0.0;
	sumxy=0.0;

	for (j=0; j<npixels_acs; j++)
	  {

	    sumx += model_bulge_filter3[j];
	    sumy += data_filter3[j];
	    sumxx += pow(model_bulge_filter3[j],2);
	    sumyy += pow(data_filter3[j],2);
	    sumxy += data_filter3[j]*model_bulge_filter3[j];
	  }

  	denom = sumxx;


	if (denom != 0.)
	  {
	    amp_etype_filter3=sumxy/sumxx;
	    amp_bulge_filter3=amp_etype_filter3;

	    for (j=0;j<npixels_acs;j++)
	      {
		lval += pow((data_filter3[j]- (amp_etype_filter3 * model_bulge_filter3[j])),2);

	      }
	    if (amp_etype_filter3<0.0) { penalty_amp_etype_filter3 = 1e7;}

	  }
	else
	  {
	    lval += 9e99;
	  }

	/* 4th filter */

	sumx=0.0;
	sumy=0.0;
	sumxx=0.0;
	sumyy=0.0;
	sumxy=0.0;

	for (j=0; j<npixels_acs; j++)
	  {

	    sumx += model_bulge_filter4[j];
	    sumy += data_filter4[j];
	    sumxx += pow(model_bulge_filter4[j],2);
	    sumyy += pow(data_filter4[j],2);
	    sumxy += data_filter4[j]*model_bulge_filter4[j];
	  }

  	denom = sumxx;


	if (denom != 0.)
	  {
	    amp_etype_filter4=sumxy/sumxx;
	    amp_bulge_filter4=amp_etype_filter4;

	    for (j=0;j<npixels_acs;j++)
	      {
		lval += pow((data_filter4[j]- (amp_etype_filter4 * model_bulge_filter4[j])),2);

	      }

	     if (amp_etype_filter4<0.0) { penalty_amp_etype_filter4 = 1e7;}

	  }
	else
	  {
	    lval += 9e99;
	  }

	sumx=0.0;
	sumy=0.0;
	sumxx=0.0;
	sumyy=0.0;
	sumxy=0.0;


	/* 5th filter */

	sumx=0.0;
	sumy=0.0;
	sumxx=0.0;
	sumyy=0.0;
	sumxy=0.0;

	for (j=0; j<npixels_wfc3; j++)
	  {

	    sumx += model_bulge_filter5[j];
	    sumy += data_filter5[j];
	    sumxx += pow(model_bulge_filter5[j],2);
	    sumyy += pow(data_filter5[j],2);
	    sumxy += data_filter5[j]*model_bulge_filter5[j];
	  }

  	denom = sumxx;


	if (denom != 0.)
	  {
	    amp_etype_filter5=sumxy/sumxx;
	    amp_bulge_filter5=amp_etype_filter5;

	    for (j=0;j<npixels_wfc3;j++)
	      {
		lval += pow((data_filter5[j]- (amp_etype_filter5 * model_bulge_filter5[j])),2);

	      }

	     if (amp_etype_filter5<0.0) { penalty_amp_etype_filter5 = 1e7;}

	  }
	else
	  {
	    lval += 9e99;
	  }

	sumx=0.0;
	sumy=0.0;
	sumxx=0.0;
	sumyy=0.0;
	sumxy=0.0;


       /* 6th filter */

	sumx=0.0;
	sumy=0.0;
	sumxx=0.0;
	sumyy=0.0;
	sumxy=0.0;

	for (j=0; j<npixels_wfc3; j++)
	  {

	    sumx += model_bulge_filter6[j];
	    sumy += data_filter6[j];
	    sumxx += pow(model_bulge_filter6[j],2);
	    sumyy += pow(data_filter6[j],2);
	    sumxy += data_filter6[j]*model_bulge_filter6[j];
	  }

  	denom = sumxx;


	if (denom != 0.)
	  {
	    amp_etype_filter6=sumxy/sumxx;
	    amp_bulge_filter6=amp_etype_filter6;

	    for (j=0;j<npixels_wfc3;j++)
	      {
		lval += pow((data_filter6[j]- (amp_etype_filter6 * model_bulge_filter6[j])),2);

	      }

	     if (amp_etype_filter6<0.0) { penalty_amp_etype_filter6 = 1e7;}

	  }
	else
	  {
	    lval += 9e99;
	  }

	sumx=0.0;
	sumy=0.0;
	sumxx=0.0;
	sumyy=0.0;
	sumxy=0.0;


	/* 7th filter */

	sumx=0.0;
	sumy=0.0;
	sumxx=0.0;
	sumyy=0.0;
	sumxy=0.0;

	for (j=0; j<npixels_wfc3; j++)
	  {

	    sumx += model_bulge_filter7[j];
	    sumy += data_filter7[j];
	    sumxx += pow(model_bulge_filter7[j],2);
	    sumyy += pow(data_filter7[j],2);
	    sumxy += data_filter7[j]*model_bulge_filter7[j];
	  }

  	denom = sumxx;


	if (denom != 0.)
	  {
	    amp_etype_filter7=sumxy/sumxx;
	    amp_bulge_filter7=amp_etype_filter7;

	    for (j=0;j<npixels_wfc3;j++)
	      {
		lval += pow((data_filter7[j]- (amp_etype_filter7 * model_bulge_filter7[j])),2);

	      }

	     if (amp_etype_filter7<0.0) { penalty_amp_etype_filter7 = 1e7;}

	  }
	else
	  {
	    lval += 9e99;
	  }

	sumx=0.0;
	sumy=0.0;
	sumxx=0.0;
	sumyy=0.0;
	sumxy=0.0;

	//printf("E-type amp, chi2 %lf %lf\n",amp_etype,lval);
	/* apply penalties for unphysical ranges */
	lval = lval + penalty_e1e2quad + penalty_scaling + penalty_amp_etype_filter1 + penalty_amp_etype_filter2 + penalty_amp_etype_filter3 + penalty_amp_etype_filter4 + penalty_amp_etype_filter5 + penalty_amp_etype_filter6 + penalty_amp_etype_filter7;

      }



     /* Continue for other types */
   if (itype>0)
      {


	denom=0.0;
	lval=0.0;
	sumxx=0.0;
	sumzz=0.0;
	sumyy=0.0;
	sumxy=0.0;
	sumzy=0.0;
	sumxz=0.0;

	/* 1 filter */
	for (j=0; j<npixels_acs; j++)
	  {

	    sumxx += pow(model_disk_filter1[j],2);
	    sumzz += pow(model_bulge_filter1[j],2);
	    sumyy += pow(data_filter1[j],2);
	    sumxy += data_filter1[j]*model_disk_filter1[j];
	    sumzy += data_filter1[j]*model_bulge_filter1[j];
	    sumxz += model_disk_filter1[j]*model_bulge_filter1[j];
	  }

	denom = sumxx*sumzz-sumxz*sumxz;

	if (denom != 0.0)
	  {

	    amp_disk_filter1 = (sumzz*sumxy-sumxz*sumzy)/denom; //disk
	    amp_bulge_filter1 = (sumxx*sumzy-sumxz*sumxy)/denom;  // bulge

	    //printf("amp(bulge), amp(disk) amp(bulge)/amp(total): %lf %lf %lf\n", amp_bulge, amp_disk,amp_bulge/(amp_bulge+amp_disk));

	    for (j=0;j<npixels_acs;j++)
	      {
		lval = lval + pow(data_filter1[j]- ((amp_disk_filter1 * model_disk_filter1[j]) + (amp_bulge_filter1 * model_bulge_filter1[j]) ), 2);

	      }


	  }
	else
	  {
	    printf("filter 1 sumxx: %lf\n",sumxx);
	    printf("filter 1 sumzz: %lf\n",sumzz);
	    printf("filter 1 sumxz: %lf\n",sumxz);
	    lval += 9e99 ;
	  }

	sumxx=0.0;
	sumyy=0.0;
	sumxy=0.0;
	sumzy=0.0;
	sumxz=0.0;
	sumzz=0.0;
	denom=0.0;

	/* 2nd filter */
	 for (j=0; j<npixels_acs; j++)
	  {

	    sumxx += pow(model_disk_filter2[j],2);
	    sumzz += pow(model_bulge_filter2[j],2);
	    sumyy += pow(data_filter2[j],2);
	    sumxy += data_filter2[j]*model_disk_filter2[j];
	    sumzy += data_filter2[j]*model_bulge_filter2[j];
	    sumxz += model_disk_filter2[j]*model_bulge_filter2[j];
	  }

	denom = sumxx*sumzz-sumxz*sumxz;

	if (denom != 0.0)
	  {

	    amp_disk_filter2 = (sumzz*sumxy-sumxz*sumzy)/denom; //disk
	    amp_bulge_filter2 = (sumxx*sumzy-sumxz*sumxy)/denom;  // bulge

	    for (j=0;j<npixels_acs;j++)
	      {
		lval = lval + pow(data_filter2[j]- ((amp_disk_filter2 * model_disk_filter2[j]) + (amp_bulge_filter2 * model_bulge_filter2[j]) ) ,2);

	      }


	  }
	else
	  {
	    printf("filter 2 sumxx: %lf\n",sumxx);
	    printf("filter 2 sumzz: %lf\n",sumzz);
	    printf("filter 2 sumxz: %lf\n",sumxz);
	    lval += 9e99 ;
	  }

	sumxx=0.0;
	sumyy=0.0;
	sumxy=0.0;
	sumzy=0.0;
	sumxz=0.0;
	sumzz=0.0;

	/* 3rd filter */

	for (j=0; j<npixels_acs; j++)
	  {

	    sumxx += pow(model_disk_filter3[j],2);
	    sumzz += pow(model_bulge_filter3[j],2);
	    sumyy += pow(data_filter3[j],2);
	    sumxy += data_filter3[j]*model_disk_filter3[j];
	    sumzy += data_filter3[j]*model_bulge_filter3[j];
	    sumxz += model_disk_filter3[j]*model_bulge_filter3[j];
	  }

	denom = sumxx*sumzz-sumxz*sumxz;

	if (denom != 0.0)
	  {

	    amp_disk_filter3 = (sumzz*sumxy-sumxz*sumzy)/denom; //disk
	    amp_bulge_filter3 = (sumxx*sumzy-sumxz*sumxy)/denom;  // bulge

	    for (j=0;j<npixels_acs;j++)
	      {
		lval = lval + pow(data_filter3[j]- ((amp_disk_filter3 * model_disk_filter3[j]) + (amp_bulge_filter3 * model_bulge_filter3[j]) ) ,2);

	      }


	  }
	else
	  {
	    printf("filter 3 sumxx: %lf\n",sumxx);
	    printf("filter 3 sumzz: %lf\n",sumzz);
	    printf("filter 3 sumxz: %lf\n",sumxz);
	    lval += 9e99 ;
	  }

	sumxx=0.0;
	sumyy=0.0;
	sumxy=0.0;
	sumzy=0.0;
	sumxz=0.0;
	sumzz=0.0;

	/* 4th filter */

	 for (j=0; j<npixels_acs; j++)
	  {

	    sumxx += pow(model_disk_filter4[j],2);
	    sumzz += pow(model_bulge_filter4[j],2);
	    sumyy += pow(data_filter4[j],2);
	    sumxy += data_filter4[j]*model_disk_filter4[j];
	    sumzy += data_filter4[j]*model_bulge_filter4[j];
	    sumxz += model_disk_filter4[j]*model_bulge_filter4[j];
	  }

	denom = sumxx*sumzz-sumxz*sumxz;

	if (denom != 0.0)
	  {

	    amp_disk_filter4 = (sumzz*sumxy-sumxz*sumzy)/denom; //disk
	    amp_bulge_filter4 = (sumxx*sumzy-sumxz*sumxy)/denom;  // bulge

	    for (j=0;j<npixels_acs;j++)
	      {
		lval = lval + pow(data_filter4[j]- ((amp_disk_filter4 * model_disk_filter4[j]) + (amp_bulge_filter4 * model_bulge_filter4[j]) ) ,2);

	      }


	  }
	else
	  {
	    printf("filter 4 sumxx: %lf\n",sumxx);
	    printf("filter 4 sumzz: %lf\n",sumzz);
	    printf("filter 4 sumxz: %lf\n",sumxz);
	    lval += 9e99 ;
	  }

	sumxx=0.0;
	sumyy=0.0;
	sumxy=0.0;
	sumzy=0.0;
	sumxz=0.0;
	sumzz=0.0;


	 /* 5th filter */

	 for (j=0; j<npixels_wfc3; j++)
	  {

	    sumxx += pow(model_disk_filter5[j],2);
	    sumzz += pow(model_bulge_filter5[j],2);
	    sumyy += pow(data_filter5[j],2);
	    sumxy += data_filter5[j]*model_disk_filter5[j];
	    sumzy += data_filter5[j]*model_bulge_filter5[j];
	    sumxz += model_disk_filter5[j]*model_bulge_filter5[j];
	  }

	denom = sumxx*sumzz-sumxz*sumxz;

	if (denom != 0.0)
	  {

	    amp_disk_filter5 = (sumzz*sumxy-sumxz*sumzy)/denom; //disk
	    amp_bulge_filter5 = (sumxx*sumzy-sumxz*sumxy)/denom;  // bulge

	    for (j=0;j<npixels_wfc3;j++)
	      {
		lval = lval + pow(data_filter5[j]- ((amp_disk_filter5 * model_disk_filter5[j]) + (amp_bulge_filter5 * model_bulge_filter5[j]) ) ,2);

	      }


	  }
	else
	  {
	    printf("filter 5 sumxx: %lf\n",sumxx);
	    printf("filter 5 sumzz: %lf\n",sumzz);
	    printf("filter 5 sumxz: %lf\n",sumxz);
	    lval += 9e99 ;
	  }

	sumxx=0.0;
	sumyy=0.0;
	sumxy=0.0;
	sumzy=0.0;
	sumxz=0.0;
	sumzz=0.0;


	 /* 6th filter */

	 for (j=0; j<npixels_wfc3; j++)
	  {

	    sumxx += pow(model_disk_filter6[j],2);
	    sumzz += pow(model_bulge_filter6[j],2);
	    sumyy += pow(data_filter6[j],2);
	    sumxy += data_filter6[j]*model_disk_filter6[j];
	    sumzy += data_filter6[j]*model_bulge_filter6[j];
	    sumxz += model_disk_filter6[j]*model_bulge_filter6[j];
	  }

	denom = sumxx*sumzz-sumxz*sumxz;

	if (denom != 0.0)
	  {

	    amp_disk_filter6 = (sumzz*sumxy-sumxz*sumzy)/denom; //disk
	    amp_bulge_filter6 = (sumxx*sumzy-sumxz*sumxy)/denom;  // bulge

	    for (j=0;j<npixels_wfc3;j++)
	      {
		lval = lval + pow(data_filter6[j]- ((amp_disk_filter6 * model_disk_filter6[j]) + (amp_bulge_filter6 * model_bulge_filter6[j]) ) ,2);

	      }


	  }
	else
	  {
	    printf("filter 6 sumxx: %lf\n",sumxx);
	    printf("filter 6 sumzz: %lf\n",sumzz);
	    printf("filter 6 sumxz: %lf\n",sumxz);
	    lval += 9e99 ;
	  }

	sumxx=0.0;
	sumyy=0.0;
	sumxy=0.0;
	sumzy=0.0;
	sumxz=0.0;
	sumzz=0.0;



	/* 7th filter */

	 for (j=0; j<npixels_wfc3; j++)
	  {

	    sumxx += pow(model_disk_filter7[j],2);
	    sumzz += pow(model_bulge_filter7[j],2);
	    sumyy += pow(data_filter7[j],2);
	    sumxy += data_filter7[j]*model_disk_filter7[j];
	    sumzy += data_filter7[j]*model_bulge_filter7[j];
	    sumxz += model_disk_filter7[j]*model_bulge_filter7[j];
	  }

	denom = sumxx*sumzz-sumxz*sumxz;

	if (denom != 0.0)
	  {

	    amp_disk_filter7 = (sumzz*sumxy-sumxz*sumzy)/denom; //disk
	    amp_bulge_filter7 = (sumxx*sumzy-sumxz*sumxy)/denom;  // bulge

	    for (j=0;j<npixels_wfc3;j++)
	      {
		lval = lval + pow(data_filter7[j]- ((amp_disk_filter7 * model_disk_filter7[j]) + (amp_bulge_filter7 * model_bulge_filter7[j]) ) ,2);

	      }


	  }
	else
	  {
	    printf("filter 7 sumxx: %lf\n",sumxx);
	    printf("filter 7 sumzz: %lf\n",sumzz);
	    printf("filter 7 sumxz: %lf\n",sumxz);
	    lval += 9e99 ;
	  }

	sumxx=0.0;
	sumyy=0.0;
	sumxy=0.0;
	sumzy=0.0;
	sumxz=0.0;
	sumzz=0.0;


	if (amp_bulge_filter1<0.0) { penalty_amp_bulge_filter1 = 1e7; }
	if (amp_disk_filter1<0.0) { penalty_amp_disk_filter1 = 1e7; }
	if (amp_bulge_filter2<0.0) { penalty_amp_bulge_filter2 = 1e7; }
	if (amp_disk_filter2<0.0) { penalty_amp_disk_filter2 = 1e7; }
	if (amp_bulge_filter3<0.0) { penalty_amp_bulge_filter3 = 1e7; }
	if (amp_disk_filter3<0.0) { penalty_amp_disk_filter3 = 1e7; }
	if (amp_bulge_filter4<0.0) { penalty_amp_bulge_filter4 = 1e7; }
	if (amp_disk_filter4<0.0) { penalty_amp_disk_filter4 = 1e7; }
	if (amp_bulge_filter5<0.0) { penalty_amp_bulge_filter5 = 1e7; }
	if (amp_disk_filter5<0.0) { penalty_amp_disk_filter5 = 1e7; }
	if (amp_bulge_filter6<0.0) { penalty_amp_bulge_filter6 = 1e7; }
	if (amp_disk_filter6<0.0) { penalty_amp_disk_filter6 = 1e7; }
	if (amp_bulge_filter7<0.0) { penalty_amp_bulge_filter7 = 1e7; }
	if (amp_disk_filter7<0.0) { penalty_amp_disk_filter7 = 1e7; }



	amp_bulge_type_filter1 = amp_bulge_filter1;
	amp_disk_type_filter1 = amp_disk_filter1;
	amp_bulge_type_filter2 = amp_bulge_filter2;
	amp_disk_type_filter2 = amp_disk_filter2;
	amp_bulge_type_filter3 = amp_bulge_filter3;
	amp_disk_type_filter3 = amp_disk_filter3;
	amp_bulge_type_filter4 = amp_bulge_filter4;
	amp_disk_type_filter4 = amp_disk_filter4;
	amp_bulge_type_filter5 = amp_bulge_filter5;
	amp_disk_type_filter5 = amp_disk_filter5;
	amp_bulge_type_filter6 = amp_bulge_filter6;
	amp_disk_type_filter6 = amp_disk_filter6;
	amp_bulge_type_filter7 = amp_bulge_filter7;
	amp_disk_type_filter7 = amp_disk_filter7;

	lval = lval + penalty_e1e2quad + penalty_scaling + penalty_amp_bulge_filter1 + penalty_amp_disk_filter1 + penalty_amp_bulge_filter2 + penalty_amp_disk_filter2 + penalty_amp_bulge_filter3 + penalty_amp_disk_filter3 + penalty_amp_bulge_filter4 + penalty_amp_disk_filter4  + penalty_amp_bulge_filter5 + penalty_amp_disk_filter5 +  penalty_amp_bulge_filter6 + penalty_amp_disk_filter6 + penalty_amp_bulge_filter7 + penalty_amp_disk_filter7 ;// + penalty_xcentre + penalty_ycentre ;

	//printf("amp_bulge_filter1 amp_bulge_filter2 amp_disk_filter1 amp_disk_filter2 ratio_amp_bulge ratio_amp_disk %lf %lf %lf %lf %lf %lf\n",amp_bulge_filter1, amp_bulge_filter2, amp_disk_filter1, amp_disk_filter2, amp_bulge_filter2/amp_bulge_filter1, amp_disk_filter2/amp_disk_filter1);

      }

  free(data_filter1);
  free(model_disk_filter1);
  free(model_bulge_filter1);
  free(data_filter2);
  free(model_disk_filter2);
  free(model_bulge_filter2);
  free(data_filter3);
  free(model_disk_filter3);
  free(model_bulge_filter3);
  free(data_filter4);
  free(model_disk_filter4);
  free(model_bulge_filter4);
  free(data_filter5);
  free(model_disk_filter5);
  free(model_bulge_filter5);
  free(data_filter6);
  free(model_disk_filter6);
  free(model_bulge_filter6);
  free(data_filter7);
  free(model_disk_filter7);
  free(model_bulge_filter7);

  return lval;

}

double (*my_f)();


double fitjointmodel(double scaling, double e1,double e2, double xcentre, double ycentre, double scaling_step, double e1_step, double e2_step, double xcentre_step, double ycentre_step, double *rvals)
{


  double chi2;
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *ss, *x;
  gsl_multimin_function minex_func;

  size_t iter = 0;
  int status, j, comp=0
;
  double size;

  //printf("ENTERED fitjointmodel() 2");


  /* Give Starting values of parameters */

  x = gsl_vector_alloc (Nparams);

  gsl_vector_set (x, 0, scaling); //set parameters
  gsl_vector_set (x, 1, e1); //set parameters
  gsl_vector_set (x, 2, e2); //set parameters
  gsl_vector_set (x, 3, xcentre); //set parameters
  gsl_vector_set (x, 4, ycentre); //set parameters

  //printf("ENTERED fitjointmodel() 3");

  /* Set initial step sizes to 0.5 in scaling factor (for scale length in xyz) and 5 deg angles */
  ss = gsl_vector_alloc (Nparams);


  gsl_vector_set (ss, 0, scaling_step);  //
  gsl_vector_set (ss, 1, e1_step);  //
  gsl_vector_set (ss, 2, e2_step);  //
  gsl_vector_set (ss, 3, xcentre_step);  //
  gsl_vector_set (ss, 4, ycentre_step);  //

  // printf("ENTERED fitjointmodel() 4");


/*  specify function */
  my_f = &jointlikel;


/* Initialize method and iterate */
   minex_func.n = Nparams;
   minex_func.f = my_f;
   minex_func.params = (void *)&comp;

/*   //printf("ok3\n"); */


   s = gsl_multimin_fminimizer_alloc (T, Nparams);

/*   //printf("ok4\n"); */

   gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

/*   //printf("ok5\n"); */

   // printf("ENTERED fitjointmodel() 5");


   do
     {
       iter++;
       status = gsl_multimin_fminimizer_iterate(s);
       //printf("ENTERED fitjointmodel() 6");

       if (status)
 	break;

       size = gsl_multimin_fminimizer_size (s);
       status = gsl_multimin_test_size (size, simplex_tol);
       // printf("ENTERED fitjointmodel() 7");

       if (status == GSL_SUCCESS)

	{
 	  printf ("converged to minimum at\n");

	  printf ("%5d %f %f %f %f %f %g %lf\n",
 	      (int)iter,
 	      gsl_vector_get (s->x, 0),
	      gsl_vector_get (s->x, 1),
	      gsl_vector_get (s->x, 2),
	      gsl_vector_get (s->x, 3),
	      gsl_vector_get (s->x, 4),
 	      s->fval, size);

 	  fflush(stdout);
     }
     }
   while (status == GSL_CONTINUE && iter < Nmaxiterations);


   printf("Number of iterations %zu\n",iter);


   gsl_vector_free(x);
   gsl_vector_free(ss);

   for (j=0; j<Nparams; j++)
     rvals[j] = gsl_vector_get (s->x, j);

   chi2 = (double) s->fval;  // store chi2

   gsl_multimin_fminimizer_free (s);

   return chi2;
}






/* get image dimensions (so we can allocate array sizes) */

 void getimagedim(char *imagename)

 {
   fitsfile *afptr;    /* FITS file pointers */
    int status = 0;             /* CFITSIO status value must be initialized to zero */
     int anaxis, bnaxis;
     long anaxes[2] = { 1, 1 }, fpixel[2] = {1, 1}, bnaxes[2] = {1, 1};
     int size, bsize;

     /* open input image */

     fits_open_file(&afptr, imagename, READONLY, &status);

/* read dimensions */
     fits_get_img_dim(afptr, &anaxis, &status);
     fits_get_img_size(afptr, 2, anaxes, &status);

     if (status) {
         fits_report_error(stderr, status);      /* print error message */
 	exit(EXIT_FAILURE);
     }

     if (anaxis != 2) {
         printf
             ("Error: images with other than 2 dimensions are not supported\n");
         exit(EXIT_FAILURE);
     }

     size = anaxes[0] * anaxes[1];

     ncols_var = anaxes[0];
     nrows_var = anaxes[1];

/* close main image file */

     fits_close_file(afptr, &status);

     if (status) {
       fits_report_error(stderr, status);      /* print error message */
         exit(EXIT_FAILURE);
     }

 }


/* read image  */

  void readimage(char* imagename, double *apix)
 {

   //// char* imagename = malloc(500*sizeof(char));

   /////strcpy(imagename,"psf.fits");

    fitsfile *afptr;    /* FITS file pointers */
     int status = 0;             /* CFITSIO status value must be initialized to zero */
     int anaxis, bnaxis;
     long anaxes[2] = { 1, 1 }, fpixel[2] = {1, 1}, bnaxes[2] = {1, 1};
     int size, bsize;

 /* open input image*/
     fits_open_file(&afptr, imagename, READONLY, &status);



 /* read dimensions*/
     fits_get_img_dim(afptr, &anaxis, &status);
     fits_get_img_size(afptr, 2, anaxes, &status);

     if (status) {
         fits_report_error(stderr, status);      /* print error message */
 	exit(EXIT_FAILURE);
     }

     if (anaxis != 2) {
         printf
             ("Error: images with other than 2 dimensions are not supported\n");
         exit(EXIT_FAILURE);
     }

     size = anaxes[0] * anaxes[1];

     printf("anaxes[0]: %d\n",anaxes[0]);
     printf("anaxes[1]: %d\n",anaxes[1]);

     printf("image size: %d\n",size);

     /* read input data into image array */

       if (fits_read_pix(afptr, TDOUBLE, fpixel, size, NULL, apix, NULL, &status)) {
         printf(" error reading pixel data \n");
         exit(EXIT_FAILURE);
     }


 /* close main image file */

     fits_close_file(afptr, &status);

     if (status) {
       fits_report_error(stderr, status);      /* print error message */
         exit(EXIT_FAILURE);
     }


 }


unsigned long int random_seed()
{
  unsigned int seed, eavailable;
  struct timeval tv;
  FILE *devrandom, *devurandom, *entfile;

  devrandom = fopen("/dev/random","r");
  devurandom = fopen("/dev/urandom","r");

  entfile = fopen("/proc/sys/kernel/random/entropy_avail","r");
  eavailable = 0;
  if (entfile != NULL)
    {
      fread(&eavailable,sizeof(eavailable),1,entfile);
      fclose(entfile);
    }
  printf("available entropy %u\n", eavailable);

  if (devrandom != NULL && eavailable > 0)
    {
      fread(&seed,sizeof(seed),1,devrandom);
      printf("Got seed %u from /dev/random\n",seed);
      fclose(devrandom);
    }
  else if (devurandom != NULL)
    {
      fread(&seed,sizeof(seed),1,devurandom);
      printf("Got seed %u from /dev/urandom\n",seed);
      fclose(devurandom);
    }
  else
    {
      gettimeofday(&tv,0);
      seed = tv.tv_sec + tv.tv_usec;
      printf("Got seed %u from gettimeofday()\n",seed);
    }

  return(seed);

}




int main(int argc, char *argv[])
 {

   int nfreedom_mask;
   double chi2_check;
   char objidstring[4];
   int obji=0, nobjects, dummy1, h_acs, w_acs, h_wfc3, w_wfc3;
   FILE *fpcat, *cgfile_bulge_disk, *cgfile_in_out, *chi2file;
   double sum_weights_inside_filter1, sum_weights_inside_filter2, sum_weights_outside_filter1, sum_weights_outside_filter2;
   double sum_weights_inside_filter3, sum_weights_inside_filter4, sum_weights_outside_filter3, sum_weights_outside_filter4;
   double sum_weights_inside_filter5, sum_weights_inside_filter6, sum_weights_outside_filter5, sum_weights_outside_filter6;
   double sum_weights_inside_filter7, sum_weights_outside_filter7;
   double *mask_filter1, *mask_filter2, *mask_filter3, *mask_filter4, *mask_filter5, *mask_filter6, *mask_filter7;

   int objid;
   double ra, dec, dummy2, dummy3,dummy4,dummy5,dummy6,dummy7,dummy8,dummy9,dummy10,dummy11,dummy12,dummy13,dummy14,dummy15, dummy16, dummy17, zphot, R50_data_filter1, R50_data_filter2, R50_data_filter3, R50_data_filter4, R50_bestmodel_filter1, R50_bestmodel_filter2,  R50_bestmodel_filter3,  R50_bestmodel_filter4, R50_catalog, age_estimated, age_grid, R50;
   double BAB,VAB,RAB,IAB,ZAB,M_B,M_V,M_R,M_I,M_Z;
   double *simdata_orig_filter1, *simdata_orig_filter2, *simdata_orig_filter3, *simdata_orig_filter4, *simdata_orig_filter5, *simdata_orig_filter6, *simdata_orig_filter7, *simnoise_filter1, *simnoise_filter2, *simnoise_filter3,  *simnoise_filter4,  *simnoise_filter5,  *simnoise_filter6,  *simnoise_filter7, *weight_filter1, *weight_filter2, *weight_filter3, *weight_filter4, *weight_filter5, *weight_filter6, *weight_filter7;
   int iter=0;
   double sumflux_bulge_filter1 = 0.0;
   double sumflux_disk_filter1 = 0.0;
   double sumflux_bulge_filter2 = 0.0;
   double sumflux_disk_filter2 = 0.0;
   double sumflux_bulge_filter3 = 0.0;
   double sumflux_disk_filter3 = 0.0;
   double sumflux_bulge_filter4 = 0.0;
   double sumflux_disk_filter4 = 0.0;
   double sumflux_bulge_filter5 = 0.0;
   double sumflux_disk_filter5 = 0.0;
   double sumflux_bulge_filter6 = 0.0;
   double sumflux_disk_filter6 = 0.0;
   double sumflux_bulge_filter7 = 0.0;
   double sumflux_disk_filter7 = 0.0;

   double sumflux_data_galaxy_filter1=0.0;
   double sumflux_data_inside_filter1=0.0;
   double sumflux_data_outside_filter1=0.0;
   double sumflux_data_galaxy_filter2=0.0;
   double sumflux_data_inside_filter2=0.0;
   double sumflux_data_outside_filter2=0.0;
   double sumflux_data_galaxy_filter3=0.0;
   double sumflux_data_inside_filter3=0.0;
   double sumflux_data_outside_filter3=0.0;
   double sumflux_data_galaxy_filter4=0.0;
   double sumflux_data_inside_filter4=0.0;
   double sumflux_data_outside_filter4=0.0;
   double sumflux_data_galaxy_filter5=0.0;
   double sumflux_data_inside_filter5=0.0;
   double sumflux_data_outside_filter5=0.0;
   double sumflux_data_galaxy_filter6=0.0;
   double sumflux_data_inside_filter6=0.0;
   double sumflux_data_outside_filter6=0.0;
   double sumflux_data_galaxy_filter7=0.0;
   double sumflux_data_inside_filter7=0.0;
   double sumflux_data_outside_filter7=0.0;

   double sumflux_model_galaxy_filter1=0.0;
   double sumflux_model_inside_filter1=0.0;
   double sumflux_model_outside_filter1=0.0;
   double sumflux_model_galaxy_filter2=0.0;
   double sumflux_model_inside_filter2=0.0;
   double sumflux_model_outside_filter2=0.0;
   double sumflux_model_galaxy_filter3=0.0;
   double sumflux_model_inside_filter3=0.0;
   double sumflux_model_outside_filter3=0.0;
   double sumflux_model_galaxy_filter4=0.0;
   double sumflux_model_inside_filter4=0.0;
   double sumflux_model_outside_filter4=0.0;
   double sumflux_model_galaxy_filter5=0.0;
   double sumflux_model_inside_filter5=0.0;
   double sumflux_model_outside_filter5=0.0;
   double sumflux_model_galaxy_filter6=0.0;
   double sumflux_model_inside_filter6=0.0;
   double sumflux_model_outside_filter6=0.0;
   double sumflux_model_galaxy_filter7=0.0;
   double sumflux_model_inside_filter7=0.0;
   double sumflux_model_outside_filter7=0.0;

   double sum_weights_galaxy_filter1=0.0;
   double sum_weights_galaxy_filter2=0.0;
   double sum_weights_galaxy_filter3=0.0;
   double sum_weights_galaxy_filter4=0.0, sum_weights_galaxy_filter5=0.0, sum_weights_galaxy_filter6=0.0, sum_weights_galaxy_filter7=0.0;

   double sumrms_data_inside_filter1=0.0;
   double sumrms_data_outside_filter1=0.0;
   double sumrms_data_galaxy_filter1=0.0;
   double sumrms_data_inside_filter2=0.0;
   double sumrms_data_outside_filter2=0.0;
   double sumrms_data_galaxy_filter2=0.0;
   double sumrms_data_inside_filter3=0.0;
   double sumrms_data_outside_filter3=0.0;
   double sumrms_data_galaxy_filter3=0.0;
   double sumrms_data_inside_filter4=0.0;
   double sumrms_data_outside_filter4=0.0;
   double sumrms_data_galaxy_filter4=0.0;
   double sumrms_data_inside_filter5=0.0;
   double sumrms_data_outside_filter5=0.0;
   double sumrms_data_galaxy_filter5=0.0;
   double sumrms_data_inside_filter6=0.0;
   double sumrms_data_outside_filter6=0.0;
   double sumrms_data_galaxy_filter6=0.0;
    double sumrms_data_inside_filter7=0.0;
   double sumrms_data_outside_filter7=0.0;
   double sumrms_data_galaxy_filter7=0.0;

   double mag_data_inside_filter1=0.0;
   double mag_data_outside_filter1=0.0;
   double mag_data_galaxy_filter1=0.0;
   double mag_data_inside_filter2=0.0;
   double mag_data_outside_filter2=0.0;
   double mag_data_galaxy_filter2=0.0;
   double mag_data_inside_filter3=0.0;
   double mag_data_outside_filter3=0.0;
   double mag_data_galaxy_filter3=0.0;
   double mag_data_inside_filter4=0.0;
   double mag_data_outside_filter4=0.0;
   double mag_data_galaxy_filter4=0.0;
   double mag_data_inside_filter5=0.0;
   double mag_data_outside_filter5=0.0;
   double mag_data_galaxy_filter5=0.0;
   double mag_data_inside_filter6=0.0;
   double mag_data_outside_filter6=0.0;
   double mag_data_galaxy_filter6=0.0;
   double mag_data_inside_filter7=0.0;
   double mag_data_outside_filter7=0.0;
   double mag_data_galaxy_filter7=0.0;

   double mag_model_inside_filter1=0.0;
   double mag_model_outside_filter1=0.0;
   double mag_model_galaxy_filter1=0.0;
   double mag_model_inside_filter2=0.0;
   double mag_model_outside_filter2=0.0;
   double mag_model_galaxy_filter2=0.0;
   double mag_model_inside_filter3=0.0;
   double mag_model_outside_filter3=0.0;
   double mag_model_galaxy_filter3=0.0;
   double mag_model_inside_filter4=0.0;
   double mag_model_outside_filter4=0.0;
   double mag_model_galaxy_filter4=0.0;
   double mag_model_inside_filter5=0.0;
   double mag_model_outside_filter5=0.0;
   double mag_model_galaxy_filter5=0.0;
   double mag_model_inside_filter6=0.0;
   double mag_model_outside_filter6=0.0;
   double mag_model_galaxy_filter6=0.0;
   double mag_model_inside_filter7=0.0;
   double mag_model_outside_filter7=0.0;
   double mag_model_galaxy_filter7=0.0;

   double sumflux_bulge_inside_filter1=0.0;
   double sumflux_disk_inside_filter1=0.0;
   double sumflux_bulge_outside_filter1=0.0;
   double sumflux_disk_outside_filter1=0.0;
   double sumflux_bulge_inside_filter2=0.0;
   double sumflux_disk_inside_filter2=0.0;
   double sumflux_bulge_outside_filter2=0.0;
   double sumflux_disk_outside_filter2=0.0;
   double sumflux_bulge_inside_filter3=0.0;
   double sumflux_disk_inside_filter3=0.0;
   double sumflux_bulge_outside_filter3=0.0;
   double sumflux_disk_outside_filter3=0.0;
   double sumflux_bulge_inside_filter4=0.0;
   double sumflux_disk_inside_filter4=0.0;
   double sumflux_bulge_outside_filter4=0.0;
   double sumflux_disk_outside_filter4=0.0;
   double sumflux_bulge_inside_filter5=0.0;
   double sumflux_disk_inside_filter5=0.0;
   double sumflux_bulge_outside_filter5=0.0;
   double sumflux_disk_outside_filter5=0.0;
   double sumflux_bulge_inside_filter6=0.0;
   double sumflux_disk_inside_filter6=0.0;
   double sumflux_bulge_outside_filter6=0.0;
   double sumflux_disk_outside_filter6=0.0;
   double sumflux_bulge_inside_filter7=0.0;
   double sumflux_disk_inside_filter7=0.0;
   double sumflux_bulge_outside_filter7=0.0;
   double sumflux_disk_outside_filter7=0.0;

   double colour_data_galaxy=0.0;
   double colour_model_galaxy=0.0;
   double colour_data_inside, colour_data_outside,colour_model_inside, colour_model_outside;
   int tt=0; int ss=0;
   double totalflux_inside_summed_filter1 = 0.0;
   double totalflux_inside_bulge_filter1 = 0.0;
   double totalflux_inside_disk_filter1 = 0.0;
   double totalflux_outside_summed_filter1 = 0.0;
   double totalflux_outside_bulge_filter1 = 0.0;
   double totalflux_outside_disk_filter1 = 0.0;
   double totalflux_inside_summed_filter2 = 0.0;
   double totalflux_inside_bulge_filter2 = 0.0;
   double totalflux_inside_disk_filter2 = 0.0;
   double totalflux_outside_summed_filter2 =0.0;
   double totalflux_outside_bulge_filter2 = 0.0;
   double totalflux_outside_disk_filter2 = 0.0;
   double totalflux_inside_summed_filter3 = 0.0;
   double totalflux_inside_bulge_filter3 = 0.0;
   double totalflux_inside_disk_filter3 = 0.0;
   double totalflux_outside_summed_filter3 = 0.0;
   double totalflux_outside_bulge_filter3 = 0.0;
   double totalflux_outside_disk_filter3 = 0.0;
   double totalflux_inside_summed_filter4 = 0.0;
   double totalflux_inside_bulge_filter4 = 0.0;
   double totalflux_inside_disk_filter4 = 0.0;
   double totalflux_outside_summed_filter4 =0.0;
   double totalflux_outside_bulge_filter4 = 0.0;
   double totalflux_outside_disk_filter4 = 0.0;
   double totalflux_outside_summed_filter5 =0.0;
   double totalflux_outside_bulge_filter5 = 0.0;
   double totalflux_outside_disk_filter5 = 0.0;
   double totalflux_outside_summed_filter6 =0.0;
   double totalflux_outside_bulge_filter6 = 0.0;
   double totalflux_outside_disk_filter6 = 0.0;
   double totalflux_outside_summed_filter7 =0.0;
   double totalflux_outside_bulge_filter7 = 0.0;
   double totalflux_outside_disk_filter7 = 0.0;

   double penalty_amp_etype_filter1=0.0, penalty_amp_etype_filter2=0.0, penalty_amp_bulge_filter1=0.0, penalty_amp_bulge_filter2=0.0, penalty_amp_disk_filter1=0.0, penalty_amp_disk_filter2=0.0;
   double penalty_amp_etype_filter3=0.0, penalty_amp_etype_filter4=0.0, penalty_amp_bulge_filter3=0.0, penalty_amp_bulge_filter4=0.0, penalty_amp_disk_filter3=0.0, penalty_amp_disk_filter4=0.0;
   double penalty_amp_etype_filter5=0.0, penalty_amp_etype_filter6=0.0, penalty_amp_bulge_filter5=0.0, penalty_amp_bulge_filter6=0.0, penalty_amp_disk_filter5=0.0, penalty_amp_disk_filter6=0.0;
   double penalty_amp_etype_filter7=0.0, penalty_amp_bulge_filter7=0.0,  penalty_amp_disk_filter7=0.0;
   double btratio_filter1, btratio_filter2, btratio_filter3, btratio_filter4, btratio_filter5, btratio_filter6, btratio_filter7;

   double vi_central_pixel_data, vierr_central_pixel_data, vi_central_pixel_model_summed_unconvolved, vi_central_pixel_model_summed_convolved, vi_central_pixel_model_bulge_convolved, vi_central_pixel_model_bulge_unconvolved, vi_central_pixel_model_disk_convolved, vi_central_pixel_model_disk_unconvolved;
   int  ipoints;
   double *data_filter1, *model_disk_filter1, *model_bulge_filter1, *data_filter2, *model_disk_filter2, *model_bulge_filter2;
   double *data_filter3, *model_disk_filter3, *model_bulge_filter3, *data_filter4, *model_disk_filter4, *model_bulge_filter4;
   double *data_filter5, *model_disk_filter5, *model_bulge_filter5, *data_filter6, *model_disk_filter6, *model_bulge_filter6;
   double *data_filter7,  *model_disk_filter7, *model_bulge_filter7;

   double datamax_filter1, datamax_filter2, residual_filter1, residual_filter2, residualmax_filter1, residualmax_filter2;
   double datamax_filter3, datamax_filter4, residual_filter3, residual_filter4, residualmax_filter3, residualmax_filter4;
   double datamax_filter5, datamax_filter6, residual_filter5, residual_filter6, residualmax_filter5, residualmax_filter6;
   double datamax_filter7, residual_filter7, residualmax_filter7;

   double lval,sumx, sumy,sumxx,sumzz,sumyy,sumxy,sumzy,sumxz,amp_disk,amp_etype_filter1, amp_etype_filter2,denom,amp_bulge_filter1, amp_disk_filter1, amp_bulge_filter2, amp_disk_filter2, amp_bulge_global_filter1, amp_bulge_global_filter2, amp_disk_global_filter1, amp_disk_global_filter2;
   double amp_etype_filter3, amp_etype_filter4,amp_bulge_filter3, amp_disk_filter3, amp_bulge_filter4, amp_disk_filter4, amp_bulge_global_filter3, amp_bulge_global_filter4, amp_disk_global_filter3, amp_disk_global_filter4;
   double amp_etype_filter5, amp_etype_filter6,amp_bulge_filter5, amp_disk_filter5, amp_bulge_filter6, amp_disk_filter6, amp_bulge_global_filter5, amp_bulge_global_filter6, amp_disk_global_filter5, amp_disk_global_filter6;
   double amp_etype_filter7, amp_bulge_filter7, amp_disk_filter7, amp_bulge_global_filter7,  amp_disk_global_filter7;
   double penalty_e1e2quad=0.0, penalty_scaling=0.0;
   double f2f1_bulge_e, f2f1_disk_e, f2f1_bulge_s0, f2f1_disk_s0, f2f1_bulge_sa, f2f1_disk_sa, f2f1_bulge_sb, f2f1_disk_sb, f2f1_bulge_sc, f2f1_disk_sc, f2f1_bulge_sd, f2f1_disk_sd;
   double dL,age_universe;
   double *cum_model_flux_radius_filter1, *cum_model_flux_radius_filter2, *cum_model_flux_radius_filter3, *cum_model_flux_radius_filter4, *cum_model_flux_radius_filter5, *cum_model_flux_radius_filter6, *cum_model_flux_radius_filter7 ;
   double *cum_data_flux_radius_filter1, *cum_data_flux_radius_filter2, *cum_data_flux_radius_filter3, *cum_data_flux_radius_filter4, *cum_data_flux_radius_filter5, *cum_data_flux_radius_filter6, *cum_data_flux_radius_filter7;
   double *cum_data_err_radius_filter1, *cum_data_err_radius_filter2, *cum_data_err_radius_filter3, *cum_data_err_radius_filter4, *cum_data_err_radius_filter5, *cum_data_err_radius_filter6, *cum_data_err_radius_filter7 ;
   double r0, radius_step, radius_max, fluxtot_data_filter1, fluxtot_data_filter2,fluxtot_model_filter1, fluxtot_model_filter2, radius_pix, cum_model_flux_istep_filter1, cum_model_flux_istep_filter2, cum_data_flux_istep_filter1, cum_data_flux_istep_filter2, r50_data_filter1, r50_data_filter2, r50_model_filter1, r50_model_filter2, maxvalue_filter1, maxvalue_filter2, cum_data_err_istep_filter1, cum_data_err_istep_filter2;
   double fluxtot_data_filter3, fluxtot_data_filter4,fluxtot_model_filter3, fluxtot_model_filter4, cum_model_flux_istep_filter3, cum_model_flux_istep_filter4, cum_data_flux_istep_filter3, cum_data_flux_istep_filter4, r50_data_filter3, r50_data_filter4, r50_model_filter3, r50_model_filter4, maxvalue_filter3, maxvalue_filter4, cum_data_err_istep_filter3, cum_data_err_istep_filter4;
   double fluxtot_data_filter5, fluxtot_data_filter6,fluxtot_model_filter5, fluxtot_model_filter6, cum_model_flux_istep_filter5, cum_model_flux_istep_filter6, cum_data_flux_istep_filter5, cum_data_flux_istep_filter6, r50_data_filter5, r50_data_filter6, r50_model_filter5, r50_model_filter6, maxvalue_filter5, maxvalue_filter6, cum_data_err_istep_filter5, cum_data_err_istep_filter6;
   double fluxtot_data_filter7,fluxtot_model_filter7, cum_model_flux_istep_filter7, cum_data_flux_istep_filter7, r50_data_filter7, r50_model_filter7, maxvalue_filter7, cum_data_err_istep_filter7;

   int x0_filter1,y0_filter1, x0_filter2, y0_filter2, x0_filter3,y0_filter3, x0_filter4, y0_filter4, x0_filter5,y0_filter5, x0_filter6, y0_filter6, x0_filter7, y0_filter7, nsteps, istep, istep_r50_data_filter1, istep_r50_data_filter2, istep_r50_model_filter1, istep_r50_model_filter2,istep_r50_data_filter3, istep_r50_data_filter4, istep_r50_model_filter3, istep_r50_model_filter4, istep_r50_data_filter5, istep_r50_data_filter6, istep_r50_model_filter5, istep_r50_model_filter6,istep_r50_data_filter7, istep_r50_model_filter7, x0_bestmodel_acs, y0_bestmodel_acs, x0_bestmodel_wfc3, y0_bestmodel_wfc3, pix_bestmodel_acs, pix_bestmodel_wfc3;

   double datafluxinner_filter1, datafluxinner_filter2, datafluxouter_filter1, datafluxouter_filter2, dataerrinner_filter1, dataerrinner_filter2, dataerrouter_filter1, dataerrouter_filter2;
   double datafluxinner_filter3, datafluxinner_filter4, datafluxouter_filter3, datafluxouter_filter4, dataerrinner_filter3, dataerrinner_filter4, dataerrouter_filter3, dataerrouter_filter4;
   double datafluxinner_filter5, datafluxinner_filter6, datafluxouter_filter5, datafluxouter_filter6, dataerrinner_filter5, dataerrinner_filter6, dataerrouter_filter5, dataerrouter_filter6;
   double datafluxinner_filter7, datafluxouter_filter7, dataerrinner_filter7, dataerrouter_filter7;
   double maginner_data_filter1, maginner_data_filter2, magouter_data_filter1, magouter_data_filter2, magerrinner_data_filter1, magerrinner_data_filter2, magerrouter_data_filter1, magerrouter_data_filter2, maginner_data_filter3, maginner_data_filter4, magouter_data_filter3, magouter_data_filter4, magerrinner_data_filter3, magerrinner_data_filter4, magerrouter_data_filter3, magerrouter_data_filter4,CG_aperture_data, CGerr_aperture_data, CG_aperture_data_err;
   double maginner_data_filter5, maginner_data_filter6, magouter_data_filter5, magouter_data_filter6, magerrinner_data_filter5, magerrinner_data_filter6, magerrouter_data_filter5, magerrouter_data_filter6, maginner_data_filter7, magouter_data_filter7,  magerrinner_data_filter7, magerrouter_data_filter7;

   double modelfluxinner_filter1, modelfluxinner_filter2, modelfluxouter_filter1, modelfluxouter_filter2,modelfluxinner_filter3, modelfluxinner_filter4, modelfluxouter_filter3, modelfluxouter_filter4 ;
   double maginner_model_filter1, maginner_model_filter2, magouter_model_filter1, magouter_model_filter2, maginner_model_filter3, maginner_model_filter4, magouter_model_filter3, magouter_model_filter4, CG_aperture_model;
   double modelfluxinner_filter5, modelfluxinner_filter6, modelfluxouter_filter5, modelfluxouter_filter6;
   double maginner_model_filter7, magouter_model_filter7;

   double totalflux_unconv_bulge_filter1=0.0;
   double totalflux_unconv_disk_filter1=0.0;
   double totalflux_conv_bulge_filter1=0.0;
   double totalflux_conv_disk_filter1=0.0;
   double totalflux_unconv_bulge_filter2=0.0;
   double totalflux_unconv_disk_filter2=0.0;
   double totalflux_conv_bulge_filter2=0.0;
   double totalflux_conv_disk_filter2=0.0;
   double totalflux_unconv_bulge_filter3=0.0;
   double totalflux_unconv_disk_filter3=0.0;
   double totalflux_conv_bulge_filter3=0.0;
   double totalflux_conv_disk_filter3=0.0;
   double totalflux_unconv_bulge_filter4=0.0;
   double totalflux_unconv_disk_filter4=0.0;
   double totalflux_conv_bulge_filter4=0.0;
   double totalflux_conv_disk_filter4=0.0;
   double totalflux_unconv_bulge_filter5=0.0;
   double totalflux_unconv_disk_filter5=0.0;
   double totalflux_conv_bulge_filter5=0.0;
   double totalflux_conv_disk_filter5=0.0;
   double totalflux_unconv_bulge_filter6=0.0;
   double totalflux_unconv_disk_filter6=0.0;
   double totalflux_conv_bulge_filter6=0.0;
   double totalflux_conv_disk_filter6=0.0;
   double totalflux_unconv_bulge_filter7=0.0;
   double totalflux_unconv_disk_filter7=0.0;
   double totalflux_conv_bulge_filter7=0.0;
   double totalflux_conv_disk_filter7=0.0;


   FILE *fp1, *fp2;
   double f2f1_ratio_bulge_global, f2f1_ratio_disk_global;
   double *psfpix_filter1, *psfpix_filter2, *psfpix_filter3,*psfpix_filter4, *psfpix_filter5, *psfpix_filter6, *psfpix_filter7;
   double xdummy, ydummy, zdummy, Idummy_filter1, norm_disk_filter1, norm_disk_filter2, norm_bulge_filter1, norm_bulge_filter2,Idummy_filter3, norm_disk_filter3, norm_disk_filter4, norm_disk_filter5, norm_disk_filter6, norm_disk_filter7, norm_bulge_filter3, norm_bulge_filter4, norm_bulge_filter5, norm_bulge_filter6, norm_bulge_filter7;
   int i=0, j=0, comp, ftsize, ig, testnumber;
   double xcentre, ycentre,scaling, theta_deg, theta_rad, posangle_rad, posangle_deg, scaling_best,e1_best, e2_best,theta_deg_best, theta_rad_best, posangle_rad_best, posangle_deg_best,  mass_disk_chosen, mass_disk, mass_disk_best, mass_bulge_chosen, e1, e2,  mass_bulge, mass_bulge_best, psftot1, psftot2, psftot3,psftot4,psftot5, psftot6, psftot7,exposuretime_filter1, exposuretime_filter2, exposuretime_filter3, exposuretime_filter4, xcentre_best, ycentre_best;
   double xcentre_step, scaling_step, e1_step, e2_step, ycentre_step;
   double xcentre_step_best, scaling_step_best, e1_step_best, e2_step_best, ycentre_step_best;
   double chi2_best_type, scaling_best_global, e1_best_global, e2_best_global, xcentre_best_global,ycentre_best_global,  theta_best_global, posangle_best_global;
   double chi2_best_global, chi2_masked;
   int itype_best_global;
   double angle, fwhm, psfe1, psfe2, snorm, shear;
   double *wpsf, *sumpsf;
   double rvals1[Nparams],rvals2[Nparams],rvals3[Nparams],rvals4[Nparams],rvals5[Nparams];
   double Ro;
   double x_scaled, y_scaled, z_scaled, x_inclined, y_inclined;
   int k=0, ncols, nrows, pix, i_index,j_index,ii,jj;
   int nrows_tot_acs, ncols_tot_acs, nrows_tot_wfc3, ncols_tot_wfc3, ii_o, jj_o, pix_o, pix_tot,jj_o_psf, ii_o_psf, pix_o_psf;
   double xmin, xmax, ymin, ymax;
   double disk_flux_tot_filter1 = 0.0, disk_flux_tot_filter2=0.0, disk_flux_tot_filter3 = 0.0, disk_flux_tot_filter4 = 0.0, disk_flux_tot_filter5=0.0, disk_flux_tot_filter6=0.0, disk_flux_tot_filter7=0.0 ;

   double *model_conv_disk_filter1, *model_conv_bulge_filter1, *model_conv_summed_filter1;
   double *model_conv_disk_filter2, *model_conv_bulge_filter2, *model_conv_summed_filter2;
   double *model_conv_disk_filter3, *model_conv_bulge_filter3, *model_conv_summed_filter3;
   double *model_conv_disk_filter4, *model_conv_bulge_filter4, *model_conv_summed_filter4;
   double *model_conv_disk_filter5, *model_conv_bulge_filter5, *model_conv_summed_filter5;
   double *model_conv_disk_filter6, *model_conv_bulge_filter6, *model_conv_summed_filter6;
   double *model_conv_disk_filter7, *model_conv_bulge_filter7, *model_conv_summed_filter7;

   double *model_conv_residual_filter1, *model_conv_residual_filter2, *model_conv_residual_filter3, *model_conv_residual_filter4, *model_conv_residual_filter5, *model_conv_residual_filter6, *model_conv_residual_filter7;
   double *model_unconv_disk_filter1, *model_unconv_bulge_filter1, *model_unconv_summed_filter1;
   double *model_unconv_disk_filter2, *model_unconv_bulge_filter2, *model_unconv_summed_filter2;
   double *model_unconv_disk_filter3, *model_unconv_bulge_filter3, *model_unconv_summed_filter3;
   double *model_unconv_disk_filter4, *model_unconv_bulge_filter4, *model_unconv_summed_filter4;
   double *model_unconv_disk_filter5, *model_unconv_bulge_filter5, *model_unconv_summed_filter5;
   double *model_unconv_disk_filter6, *model_unconv_bulge_filter6, *model_unconv_summed_filter6;
   double *model_unconv_disk_filter7, *model_unconv_bulge_filter7, *model_unconv_summed_filter7;

   double mag_bulge_filter1, mag_bulge_filter2, mag_disk_filter1, mag_disk_filter2,mag_bulge_filter3, mag_bulge_filter4, mag_bulge_filter5, mag_bulge_filter6, mag_bulge_filter7, mag_disk_filter3, mag_disk_filter4, mag_disk_filter5, mag_disk_filter6, mag_disk_filter7;
   double colour_bulge, colour_disk, CG_bulge_disk;
   double x_disk_dummy, y_disk_dummy, z_disk_dummy, I_xyz_disk_filter1_dummy, dummy5_disk, dummy5_bulge;
   double x_bulge_dummy, y_bulge_dummy, z_bulge_dummy, I_xyz_bulge_filter1_dummy;

   char* imagename = malloc(500*sizeof(char));
   char *galaxytype = malloc(20*sizeof(char));
   char *modelfilename_disk, *modelfilename_bulge;
   modelfilename_disk = (char *) calloc(500, sizeof(char));
   modelfilename_bulge = (char *) calloc(500, sizeof(char));

   char *datafitsname_filter1 = malloc(500*sizeof(char));
   char *whtfitsname_filter1 = malloc(500*sizeof(char));
   char *psffitsname_filter1 = malloc(200*sizeof(char));

   char *datafitsname_filter2 = malloc(500*sizeof(char));
   char *whtfitsname_filter2 = malloc(500*sizeof(char));
   char *psffitsname_filter2 = malloc(200*sizeof(char));

   char *datafitsname_filter3 = malloc(500*sizeof(char));
   char *whtfitsname_filter3 = malloc(500*sizeof(char));
   char *psffitsname_filter3 = malloc(200*sizeof(char));

   char *datafitsname_filter4 = malloc(500*sizeof(char));
   char *whtfitsname_filter4 = malloc(500*sizeof(char));
   char *psffitsname_filter4 = malloc(200*sizeof(char));

   char *datafitsname_filter5 = malloc(500*sizeof(char));
   char *whtfitsname_filter5 = malloc(500*sizeof(char));
   char *psffitsname_filter5 = malloc(200*sizeof(char));

   char *datafitsname_filter6 = malloc(500*sizeof(char));
   char *whtfitsname_filter6 = malloc(500*sizeof(char));
   char *psffitsname_filter6 = malloc(200*sizeof(char));

   char *datafitsname_filter7 = malloc(500*sizeof(char));
   char *whtfitsname_filter7 = malloc(500*sizeof(char));
   char *psffitsname_filter7 = malloc(200*sizeof(char));

   char *catfilename =  malloc(200*sizeof(char));
   char *catfilename_bulge_disk = malloc(200*sizeof(char));
   char *catfilename_in_out = malloc(200*sizeof(char));

   char *zootype = malloc(200*sizeof(char));

   char *inputmapspath = malloc(500*sizeof(char));
   char *outputmapspath = malloc(500*sizeof(char));
   char *filename_data_filter1 = malloc(500*sizeof(char));
   char *filename_data_filter2 = malloc(500*sizeof(char));
   char *filename_model_colour = malloc(500*sizeof(char));
   char *filename_model_unconv_filter1 = malloc(500*sizeof(char));
   char *filename_model_unconv_filter2 = malloc(500*sizeof(char));
   char *filename_model_conv_filter1 = malloc(500*sizeof(char));
   char *filename_model_conv_filter2 = malloc(500*sizeof(char));
   char *filename_model_residual_filter1 = malloc(500*sizeof(char));
   char *filename_model_residual_filter2 = malloc(500*sizeof(char));

   char *filename_data_filter3 = malloc(500*sizeof(char));
   char *filename_data_filter4 = malloc(500*sizeof(char));
   char *filename_model_unconv_filter3 = malloc(500*sizeof(char));
   char *filename_model_unconv_filter4 = malloc(500*sizeof(char));
   char *filename_model_conv_filter3 = malloc(500*sizeof(char));
   char *filename_model_conv_filter4 = malloc(500*sizeof(char));
   char *filename_model_residual_filter3 = malloc(500*sizeof(char));
   char *filename_model_residual_filter4 = malloc(500*sizeof(char));

   char *filename_data_filter5 = malloc(500*sizeof(char));
   char *filename_data_filter6 = malloc(500*sizeof(char));
   char *filename_model_unconv_filter5 = malloc(500*sizeof(char));
   char *filename_model_unconv_filter6 = malloc(500*sizeof(char));
   char *filename_model_conv_filter5 = malloc(500*sizeof(char));
   char *filename_model_conv_filter6 = malloc(500*sizeof(char));
   char *filename_model_residual_filter5 = malloc(500*sizeof(char));
   char *filename_model_residual_filter6 = malloc(500*sizeof(char));

   char *filename_data_filter7 = malloc(500*sizeof(char));
   char *filename_model_unconv_filter7 = malloc(500*sizeof(char));
   char *filename_model_conv_filter7 = malloc(500*sizeof(char));
   char *filename_model_residual_filter7 = malloc(500*sizeof(char));

   char *filename_mask_filter1 = malloc(500*sizeof(char));
   char *filename_mask_filter2 = malloc(500*sizeof(char));
   char *filename_mask_filter3 = malloc(500*sizeof(char));
   char *filename_mask_filter4 = malloc(500*sizeof(char));
   char *filename_mask_filter5 = malloc(500*sizeof(char));
   char *filename_mask_filter6 = malloc(500*sizeof(char));
   char *filename_mask_filter7 = malloc(500*sizeof(char));


   double totflux_filter1, totflux_filter2,totflux_filter3, totflux_filter4, totflux_filter5, totflux_filter6, totflux_filter7;
   double sumfluxsquared_filter1 = 0.0;
   double sumfluxsquared_filter2 = 0.0;
   double sumfluxsquared_filter3 = 0.0;
   double sumfluxsquared_filter4 = 0.0, sumfluxsquared_filter5=0.0, sumfluxsquared_filter6=0.0, sumfluxsquared_filter7=0.0;


   if (argc!=2)
     {
       printf("USAGE program cat_filename\n");
       exit(0);
     }

   /*---------------------------- */
   /* Inputs (fixed) */
   /*---------------------------- */

   strcpy(psffitsname_filter1,"psf_sect33_b.fits");
   strcpy(psffitsname_filter2,"psf_sect33_v.fits");
   strcpy(psffitsname_filter3,"psf_sect33_i.fits");
   strcpy(psffitsname_filter4,"psf_sect33_z.fits");
   strcpy(psffitsname_filter5,"psf_wfc3_y.fits");
   strcpy(psffitsname_filter6,"psf_wfc3_j.fits");
   strcpy(psffitsname_filter7,"psf_wfc3_h.fits");

   /* Input file */

   //strcpy(catfilename, "sample_bvizyjh_");
   strcat(catfilename,argv[1]);
   fpcat=fopen(catfilename,"r");

   printf("%s\n",catfilename);


   /* Output file */
   strcpy(catfilename_bulge_disk, "cg_bulge_disk_dust_");
   strcat(catfilename_bulge_disk,argv[1]);
   cgfile_bulge_disk=fopen(catfilename_bulge_disk,"a");

   strcpy(catfilename_in_out, "cg_in_out_dust_");
   strcat(catfilename_in_out,argv[1]);
   cgfile_in_out = fopen(catfilename_in_out, "a");

   chi2file=fopen("chi2_bvizyjh_bulge_disk_dust","a");

   /* Input maps path */
   strcpy(inputmapspath,"/users/welikala/data/GOODS/");

   Nfreedom = (nfilters_acs*Npix_side_acs*Npix_side_acs) + (nfilters_wfc3*Npix_side_wfc3*Npix_side_wfc3) - Nparams;
   nrows_tot_acs = Npix_side_acs;
   ncols_tot_acs = Npix_side_acs;
   npixels_acs = nrows_tot_acs * ncols_tot_acs;
   nrows_tot_wfc3 = Npix_side_wfc3;
   ncols_tot_wfc3 = Npix_side_wfc3;
   npixels_wfc3 = nrows_tot_wfc3 * ncols_tot_wfc3;

   /* Define radius of annuli for aperture colour gradient measurement */
   radius_step = 0.06; ///arcsecs
   radius_max = 3.0; //Npix_side * pixelsize_arcsec/2;
   nsteps = radius_max/radius_step;



   cum_model_flux_radius_filter1 = (double *) malloc ((nsteps) * sizeof(double));
   cum_model_flux_radius_filter2 = (double *) malloc ((nsteps) * sizeof(double));
   cum_data_flux_radius_filter1 = (double *) malloc ((nsteps) * sizeof(double));
   cum_data_flux_radius_filter2 = (double *) malloc ((nsteps) * sizeof(double));
   cum_data_err_radius_filter1 = (double *) malloc ((nsteps) * sizeof(double));
   cum_data_err_radius_filter2 = (double *) malloc ((nsteps) * sizeof(double));

   cum_model_flux_radius_filter3 = (double *) malloc ((nsteps) * sizeof(double));
   cum_model_flux_radius_filter4 = (double *) malloc ((nsteps) * sizeof(double));
   cum_data_flux_radius_filter3 = (double *) malloc ((nsteps) * sizeof(double));
   cum_data_flux_radius_filter4 = (double *) malloc ((nsteps) * sizeof(double));
   cum_data_err_radius_filter3 = (double *) malloc ((nsteps) * sizeof(double));
   cum_data_err_radius_filter4 = (double *) malloc ((nsteps) * sizeof(double));

   model_unconv_acs = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_unconv_wfc3 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_conv_acs = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_conv_wfc3 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   psf_image_acs = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   psf_image_wfc3 =  (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));

   model_conv_disk_filter1 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_conv_bulge_filter1 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_conv_disk_filter2 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_conv_bulge_filter2 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_conv_disk_filter3 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_conv_bulge_filter3 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_conv_disk_filter4 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_conv_bulge_filter4 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_conv_disk_filter5 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   model_conv_bulge_filter5 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   model_conv_disk_filter6 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   model_conv_bulge_filter6 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   model_conv_disk_filter7 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   model_conv_bulge_filter7 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));


   model_conv_summed_filter1 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_conv_summed_filter2 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_conv_summed_filter3 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_conv_summed_filter4 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_conv_summed_filter5 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   model_conv_summed_filter6 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   model_conv_summed_filter7 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));

   model_conv_residual_filter1 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_conv_residual_filter2 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_conv_residual_filter3 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_conv_residual_filter4 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_conv_residual_filter5 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   model_conv_residual_filter6 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   model_conv_residual_filter7 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));

   model_unconv_disk_filter1 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_unconv_bulge_filter1 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_unconv_disk_filter2 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_unconv_bulge_filter2 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_unconv_disk_filter3 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_unconv_bulge_filter3 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_unconv_disk_filter4 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_unconv_bulge_filter4 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_unconv_disk_filter5 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   model_unconv_bulge_filter5 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   model_unconv_disk_filter6 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   model_unconv_bulge_filter6 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   model_unconv_disk_filter7 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   model_unconv_bulge_filter7 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));

   model_unconv_summed_filter1 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_unconv_summed_filter2 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_unconv_summed_filter3 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_unconv_summed_filter4 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_unconv_summed_filter5 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   model_unconv_summed_filter6 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   model_unconv_summed_filter7 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));

   data_filter1 =  (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   data_filter2 =  (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   data_filter3 =  (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   data_filter4 =  (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   data_filter5 =  (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   data_filter6 =  (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   data_filter7 =  (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));

   model_disk_filter1 =  (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_bulge_filter1 =  (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_disk_filter2 =  (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_bulge_filter2 =  (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_disk_filter3 =  (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_bulge_filter3 =  (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_disk_filter4 =  (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_bulge_filter4 =  (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   model_disk_filter5 =  (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   model_bulge_filter5 =  (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   model_disk_filter6 =  (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   model_bulge_filter6 =  (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   model_disk_filter7 =  (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   model_bulge_filter7 =  (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));

   psf_image_filter1 =  (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   psf_image_filter2 =  (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   psf_image_filter3 =  (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   psf_image_filter4 =  (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   psf_image_filter5 =  (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   psf_image_filter6 =  (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   psf_image_filter7 =  (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));


   x_disk = (double *) malloc ((npoints_disk) * sizeof(double));
   y_disk = (double *) malloc ((npoints_disk) * sizeof(double));
   z_disk = (double *) malloc ((npoints_disk) * sizeof(double));
   I_xyz_disk_filter1 = (double *) malloc ((npoints_disk) * sizeof(double));
   x_bulge = (double *) malloc ((npoints_bulge) * sizeof(double));
   y_bulge = (double *) malloc ((npoints_bulge) * sizeof(double));
   z_bulge = (double *) malloc ((npoints_bulge) * sizeof(double));
   I_xyz_bulge_filter1 = (double *) malloc ((npoints_bulge) * sizeof(double));

   in_acs = (double*) fftw_malloc(sizeof(double)*Npix_side_acs*Npix_side_acs);
   identity_acs = (double*) fftw_malloc(sizeof(double)*Npix_side_acs*Npix_side_acs);
   final_acs=(double*) fftw_malloc(sizeof(double)*Npix_side_acs*Npix_side_acs);
   in_wfc3 = (double*) fftw_malloc(sizeof(double)*Npix_side_wfc3*Npix_side_wfc3);
   identity_wfc3 = (double*) fftw_malloc(sizeof(double)*Npix_side_wfc3*Npix_side_wfc3);
   final_wfc3=(double*) fftw_malloc(sizeof(double)*Npix_side_wfc3*Npix_side_wfc3);

   inTrans_acs = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Npix_side_acs*(Npix_side_acs/2+1));
   identityTrans_acs = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Npix_side_acs*(Npix_side_acs/2+1));
   FinalFFT_acs = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Npix_side_acs*(Npix_side_acs/2+1));
   inTrans_wfc3 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Npix_side_wfc3*(Npix_side_wfc3/2+1));
   identityTrans_wfc3 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Npix_side_wfc3*(Npix_side_wfc3/2+1));
   FinalFFT_wfc3 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Npix_side_wfc3*(Npix_side_wfc3/2+1));

   //Initialize forward plans
   plan1_acs  = fftw_plan_dft_r2c_2d(Npix_side_acs, Npix_side_acs, in_acs, inTrans_acs, FFTW_ESTIMATE);
   plan2_acs = fftw_plan_dft_r2c_2d(Npix_side_acs, Npix_side_acs, identity_acs, identityTrans_acs, FFTW_ESTIMATE);
   plan3_acs  = fftw_plan_dft_c2r_2d(Npix_side_acs, Npix_side_acs,  FinalFFT_acs, final_acs, FFTW_ESTIMATE);
   plan1_wfc3  = fftw_plan_dft_r2c_2d(Npix_side_wfc3, Npix_side_wfc3, in_wfc3, inTrans_wfc3, FFTW_ESTIMATE);
   plan2_wfc3 = fftw_plan_dft_r2c_2d(Npix_side_wfc3, Npix_side_wfc3, identity_wfc3, identityTrans_wfc3, FFTW_ESTIMATE);
   plan3_wfc3  = fftw_plan_dft_c2r_2d(Npix_side_wfc3, Npix_side_wfc3,  FinalFFT_wfc3, final_wfc3, FFTW_ESTIMATE);

   /* Local variables x,y,z for each type */

   double ** x_bulge_type, **y_bulge_type, **z_bulge_type, **I_xyz_bulge_filter1_type;
   double ** x_disk_type, **y_disk_type, **z_disk_type, **I_xyz_disk_filter1_type;
   int iii=0;

   x_bulge_type = (double**)malloc(sizeof(double*)*npoints_bulge);
   for (iii = 0; iii < npoints_bulge; iii++)
    x_bulge_type[iii] = (double*)malloc(sizeof(double)*ntypes);

   y_bulge_type = (double**)malloc(sizeof(double*)*npoints_bulge);
   for (iii = 0; iii < npoints_bulge; iii++)
     y_bulge_type[iii] = (double*)malloc(sizeof(double)*ntypes);

   z_bulge_type = (double**)malloc(sizeof(double*)*npoints_bulge);
   for (iii = 0; iii < npoints_bulge; iii++)
     z_bulge_type[iii] = (double*)malloc(sizeof(double)*ntypes);

   I_xyz_bulge_filter1_type = (double**)malloc(sizeof(double*)*npoints_bulge);
   for (iii = 0; iii < npoints_bulge; iii++)
     I_xyz_bulge_filter1_type[iii] = (double*)malloc(sizeof(double)*ntypes);

   x_disk_type = (double**)malloc(sizeof(double*)*npoints_disk);
   for (iii = 0; iii < npoints_disk; iii++)
     x_disk_type[iii] = (double*)malloc(sizeof(double)*ntypes);

   y_disk_type = (double**)malloc(sizeof(double*)*npoints_disk);
   for (iii = 0; iii < npoints_disk; iii++)
     y_disk_type[iii] = (double*)malloc(sizeof(double)*ntypes);

   z_disk_type = (double**)malloc(sizeof(double*)*npoints_disk);
   for (iii = 0; iii < npoints_disk; iii++)
     z_disk_type[iii] = (double*)malloc(sizeof(double)*ntypes);

   I_xyz_disk_filter1_type = (double**)malloc(sizeof(double*)*npoints_disk);
   for (iii = 0; iii < npoints_disk; iii++)
     I_xyz_disk_filter1_type[iii] = (double*)malloc(sizeof(double)*ntypes);


   h_acs = Npix_side_acs;
   w_acs = Npix_side_acs;
   h_wfc3 = Npix_side_wfc3;
   w_wfc3 = Npix_side_wfc3;

   /* Initialize arrays */

   for (pix=0;pix<npixels_acs;pix++)
     {

       model_unconv_bulge_filter1[pix] = 0.0;
       model_unconv_disk_filter1[pix] = 0.0;
       model_unconv_bulge_filter2[pix] = 0.0;
       model_unconv_disk_filter2[pix] = 0.0;
       model_conv_bulge_filter1[pix] = 0.0;
       model_conv_disk_filter1[pix] = 0.0;
       model_conv_bulge_filter2[pix] = 0.0;
       model_conv_disk_filter2[pix] = 0.0;

       model_unconv_summed_filter1[pix] =  0.0;
       model_unconv_summed_filter2[pix] =  0.0;
       model_conv_summed_filter1[pix] =    0.0;
       model_conv_summed_filter2[pix] =    0.0;
       model_conv_residual_filter1[pix] =  0.0;
       model_conv_residual_filter2[pix] =  0.0;

       model_unconv_bulge_filter3[pix] = 0.0;
       model_unconv_disk_filter3[pix] = 0.0;
       model_unconv_bulge_filter4[pix] = 0.0;
       model_unconv_disk_filter4[pix] = 0.0;
       model_conv_bulge_filter3[pix] = 0.0;
       model_conv_disk_filter3[pix] = 0.0;
       model_conv_bulge_filter4[pix] = 0.0;
       model_conv_disk_filter4[pix] = 0.0;

       model_unconv_summed_filter3[pix] =  0.0;
       model_unconv_summed_filter4[pix] =  0.0;
       model_conv_summed_filter3[pix] =    0.0;
       model_conv_summed_filter4[pix] =    0.0;
       model_conv_residual_filter3[pix] =  0.0;
       model_conv_residual_filter4[pix] =  0.0;

     }

   for (pix=0;pix<npixels_wfc3;pix++)
     {

       model_unconv_bulge_filter5[pix] = 0.0;
       model_unconv_disk_filter5[pix] = 0.0;
       model_unconv_bulge_filter6[pix] = 0.0;
       model_unconv_disk_filter6[pix] = 0.0;
       model_unconv_bulge_filter7[pix] = 0.0;
       model_unconv_disk_filter7[pix] = 0.0;

       model_conv_bulge_filter5[pix] = 0.0;
       model_conv_disk_filter5[pix] = 0.0;
       model_conv_bulge_filter6[pix] = 0.0;
       model_conv_disk_filter6[pix] = 0.0;
       model_conv_bulge_filter7[pix] = 0.0;
       model_conv_disk_filter7[pix] = 0.0;

       model_unconv_summed_filter5[pix] =  0.0;
       model_unconv_summed_filter6[pix] =  0.0;
       model_unconv_summed_filter7[pix] =  0.0;

       model_conv_summed_filter5[pix] =    0.0;
       model_conv_summed_filter6[pix] =    0.0;
       model_conv_summed_filter7[pix] =    0.0;

       model_conv_residual_filter5[pix] =  0.0;
       model_conv_residual_filter6[pix] =  0.0;
       model_conv_residual_filter7[pix] =  0.0;
     }


   /* ---------------------------------------- */
   /* Read in galaxy models (morphology types) */
   /* ---------------------------------------- */
       int kk=0;
       int ll=0;

       for (itype=0; itype<ntypes; itype++)
	 {
	   strcpy(modelfilename_disk, "3dmodel_disk_");
	   strcpy(modelfilename_bulge, "3dmodel_bulge_");

	   if (itype==0) { strcpy(galaxytype,"e");}
	   if (itype==1) { strcpy(galaxytype,"s0");}
	   if (itype==2) { strcpy(galaxytype,"sa");}
	   if (itype==3) { strcpy(galaxytype,"sb");}
	   if (itype==4) { strcpy(galaxytype,"sc");}
	   if (itype==5) { strcpy(galaxytype,"sd");}

	   printf("********** GALAXY TYPE %s\n",galaxytype);

	   strcat(modelfilename_disk, galaxytype);
	   strcat(modelfilename_disk, ".dat");
	   printf("********** MODEL FILE %s\n",modelfilename_disk);

	   strcat(modelfilename_bulge, galaxytype);
	   strcat(modelfilename_bulge, ".dat");
	   printf("********** MODEL FILE %s\n",modelfilename_bulge);

	   fp1 = fopen(modelfilename_disk,"r");
	   fp2 = fopen(modelfilename_bulge,"r");

	   kk=0;
	   while(fscanf(fp1, "%lf %lf %lf %lf ",&x_disk_dummy, &y_disk_dummy, &z_disk_dummy, &I_xyz_disk_filter1_dummy ) == 4 )
	     {
	       x_disk_type[kk][itype]=x_disk_dummy;
	       y_disk_type[kk][itype]=y_disk_dummy;
               z_disk_type[kk][itype]=z_disk_dummy;
               I_xyz_disk_filter1_type[kk][itype]=I_xyz_disk_filter1_dummy;

	       kk++;
	     }

	   fclose(fp1);


	   ll=0;
	   while(fscanf(fp2, "%lf %lf %lf %lf ", &x_bulge_dummy, &y_bulge_dummy, &z_bulge_dummy, &I_xyz_bulge_filter1_dummy ) == 4 )
	     {
	       x_bulge_type[ll][itype]=x_bulge_dummy;
               y_bulge_type[ll][itype]=y_bulge_dummy;
               z_bulge_type[ll][itype]=z_bulge_dummy;
               I_xyz_bulge_filter1_type[ll][itype]=I_xyz_bulge_filter1_dummy;
	       //printf("x,y,z,I: %lf %lf %lf %lf\n",x_bulge_type[ll][itype],y_bulge_type[ll][itype],z_bulge_type[ll][itype],I_xyz_bulge_filter1_type[ll][itype]);
	       ll++;
	     }

	   fclose(fp2);



	 } // loop over and store galaxy models



   /*---------------------------- */
  /* Read in the PSF (assumed unpadded and not same size as final image which is Npix_side x Npix_side ) */
  /*---------------------------- */
 /* Pad galaxy and PSF image to size Npix_side x Npix_side */

/* first get image dimensions and read image - psf image 1 */
  getimagedim(psffitsname_filter1);
  psfpix_filter1 = (double *) calloc(ncols_var*nrows_var, sizeof(double));
  readimage(psffitsname_filter1, psfpix_filter1);

/* same for psf image 2 */

  getimagedim(psffitsname_filter2);
  psfpix_filter2 = (double *) calloc(ncols_var*nrows_var, sizeof(double));
  readimage(psffitsname_filter2, psfpix_filter2);

/* same for psf image 3 */

  getimagedim(psffitsname_filter3);
  psfpix_filter3 = (double *) calloc(ncols_var*nrows_var, sizeof(double));
  readimage(psffitsname_filter3, psfpix_filter3);

/* same for psf image 4 */

  getimagedim(psffitsname_filter4);
  psfpix_filter4 = (double *) calloc(ncols_var*nrows_var, sizeof(double));
  readimage(psffitsname_filter4, psfpix_filter4);

  ncols_psf_acs = ncols_var;
  nrows_psf_acs = nrows_var;


/* same for psf image 5 */

  getimagedim(psffitsname_filter5);
  psfpix_filter5 = (double *) calloc(ncols_var*nrows_var, sizeof(double));
  readimage(psffitsname_filter5, psfpix_filter5);

  /* same for psf image 6 */

  getimagedim(psffitsname_filter6);
  psfpix_filter6 = (double *) calloc(ncols_var*nrows_var, sizeof(double));
  readimage(psffitsname_filter6, psfpix_filter6);

  /* same for psf image 7 */

  getimagedim(psffitsname_filter7);
  psfpix_filter7 = (double *) calloc(ncols_var*nrows_var, sizeof(double));
  readimage(psffitsname_filter7, psfpix_filter7);

  ncols_psf_wfc3 = ncols_var;
  nrows_psf_wfc3 = nrows_var;


  /* do padding of PSF images */

  pix_tot = 0; pix_o_psf = 0.0; jj_o_psf = 0; ii_o_psf=0;

  for (jj=0;jj<nrows_tot_acs;jj++)
    for (ii=0;ii<ncols_tot_acs;ii++)
     {
       pix_tot = (ncols_tot_acs *jj) + ii;
       psf_image_filter1[pix_tot] = 0.0;
       psf_image_filter2[pix_tot] = 0.0;
       psf_image_filter3[pix_tot] = 0.0;
       psf_image_filter4[pix_tot] = 0.0;

       if ( (ii>=(ncols_tot_acs-ncols_psf_acs)/2) && (ii<(ncols_tot_acs + ncols_psf_acs)/2) && (jj>=(nrows_tot_acs-nrows_psf_acs)/2) && (jj<(nrows_tot_acs + nrows_psf_acs)/2)  )
	 {

	   jj_o_psf =  jj - ( (nrows_tot_acs-nrows_psf_acs)/2 )  ;
	   ii_o_psf =  ii - ((ncols_tot_acs-ncols_psf_acs)/2 )   ;
	   pix_o_psf = (ncols_psf_acs*(jj_o_psf)) + ii_o_psf;
	   psf_image_filter1[pix_tot] = psfpix_filter1[pix_o_psf];
	   psf_image_filter2[pix_tot] = psfpix_filter2[pix_o_psf];
	   psf_image_filter3[pix_tot] = psfpix_filter3[pix_o_psf];
	   psf_image_filter4[pix_tot] = psfpix_filter4[pix_o_psf];

	 }

     }


  pix_tot = 0; pix_o_psf = 0; jj_o_psf = 0; ii_o_psf=0;

  for (jj=0;jj<nrows_tot_wfc3;jj++)
    for (ii=0;ii<ncols_tot_wfc3;ii++)
     {
       pix_tot = (ncols_tot_wfc3 *jj) + ii;
       psf_image_filter5[pix_tot] = 0.0;
       psf_image_filter6[pix_tot] = 0.0;
       psf_image_filter7[pix_tot] = 0.0;

       if ( (ii>=(ncols_tot_wfc3-ncols_psf_wfc3)/2) && (ii<(ncols_tot_wfc3 + ncols_psf_wfc3)/2) && (jj>=(nrows_tot_wfc3-nrows_psf_wfc3)/2) && (jj<(nrows_tot_wfc3 + nrows_psf_wfc3)/2)  )
	 {

	   jj_o_psf =  jj - ( (nrows_tot_wfc3-nrows_psf_wfc3)/2 )  ;
	   ii_o_psf =  ii - ((ncols_tot_wfc3-ncols_psf_wfc3)/2 )   ;
	   pix_o_psf = (ncols_psf_wfc3*(jj_o_psf)) + ii_o_psf;
	   psf_image_filter5[pix_tot] = psfpix_filter5[pix_o_psf];
	   psf_image_filter6[pix_tot] = psfpix_filter6[pix_o_psf];
	   psf_image_filter7[pix_tot] = psfpix_filter7[pix_o_psf];

	 }

     }


  /* Normalize the PSFs to 1 */

  psftot1 = 0.0;
  psftot2 = 0.0;
  psftot3 = 0.0;
  psftot4 = 0.0;
  psftot5 = 0.0;
  psftot6 = 0.0;
  psftot7 = 0.0;

  /* ACS */
  for (pix=0; pix<npixels_acs;pix++)
    {
      psftot1 = psftot1 + psf_image_filter1[pix];
      psftot2 = psftot2 + psf_image_filter2[pix];
      psftot3 = psftot3 + psf_image_filter3[pix];
      psftot4 = psftot4 + psf_image_filter4[pix];

    }

 for (pix=0; pix<npixels_acs;pix++)
    {
      psf_image_filter1[pix] = psf_image_filter1[pix]/psftot1;
      psf_image_filter2[pix] = psf_image_filter2[pix]/psftot2;
      psf_image_filter3[pix] = psf_image_filter3[pix]/psftot3;
      psf_image_filter4[pix] = psf_image_filter4[pix]/psftot4;

   }

 /* WFC3 */
for (pix=0; pix<npixels_wfc3;pix++)
    {
      psftot5 = psftot5 + psf_image_filter5[pix];
      psftot6 = psftot6 + psf_image_filter6[pix];
      psftot7 = psftot7 + psf_image_filter7[pix];
    }

 for (pix=0; pix<npixels_wfc3;pix++)
    {
      psf_image_filter5[pix] = psf_image_filter5[pix]/psftot5;
      psf_image_filter6[pix] = psf_image_filter6[pix]/psftot6;
      psf_image_filter7[pix] = psf_image_filter7[pix]/psftot7;

   }



   /*---------------------------- */
   /* Read in the data (assumed this is the specified Npix_side x Npix_side) */
   /*---------------------------- */
   double *countrate_filter1, *countrate_filter2,  *countrate_filter3, *countrate_filter4, *countrate_filter5, *countrate_filter6, *countrate_filter7, *wht_filter1, *wht_filter2,*wht_filter3, *wht_filter4, *wht_filter5, *wht_filter6, *wht_filter7;
   countrate_filter1=(double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   countrate_filter2=(double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   countrate_filter3=(double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   countrate_filter4=(double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   countrate_filter5=(double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   countrate_filter6=(double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   countrate_filter7=(double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));

   wht_filter1=(double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   wht_filter2=(double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   wht_filter3=(double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   wht_filter4=(double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   wht_filter5 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   wht_filter6 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   wht_filter7 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));

   mask_filter1=(double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   mask_filter2=(double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   mask_filter3=(double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   mask_filter4=(double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   mask_filter5 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   mask_filter6 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   mask_filter7 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));

   simdata_orig_filter1 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   simdata_orig_filter2 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   simdata_orig_filter3 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   simdata_orig_filter4 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   simdata_orig_filter5 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   simdata_orig_filter6 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   simdata_orig_filter7 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));


   simdata_filter1 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   simrms_filter1 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   simnoise_filter1 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   weight_filter1 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));

   simdata_filter2 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   simrms_filter2 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   simnoise_filter2 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   weight_filter2 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));

   simdata_filter3 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   simrms_filter3 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   simnoise_filter3 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   weight_filter3 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));

   simdata_filter4 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   simrms_filter4 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   simnoise_filter4 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));
   weight_filter4 = (double *) malloc ((Npix_side_acs*Npix_side_acs) * sizeof(double));

   simdata_filter5 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   simrms_filter5 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   simnoise_filter5 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   weight_filter5 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));

   simdata_filter6 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   simrms_filter6 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   simnoise_filter6 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   weight_filter6 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));

   simdata_filter7 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   simrms_filter7 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   simnoise_filter7 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));
   weight_filter7 = (double *) malloc ((Npix_side_wfc3*Npix_side_wfc3) * sizeof(double));


   /* ----------------------- */
   /* Read catalog of objects */
   /* ----------------------- */

   /* GOODS ACS catalog matched with z-band SEXtractor catalog - applied cuts in R50>0.3" and I<26.0*/

   obji=0;



   while(fscanf(fpcat,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",&objid,&ra,&dec,&zphot,&R50,&VAB,&IAB,&M_I,&exposuretime_filter1, &exposuretime_filter2, &exposuretime_filter3, &exposuretime_filter4,&arcsecstokpc,&dL,&age_universe)==15)

       {


       R50 = R50 * pixelsize_arcsec_acs;  //ACS size in I band in arc seconds
       sprintf(objidstring,"%d",objid);

       strcpy(datafitsname_filter1,"/users/welikala/data/GOODS/ACS/stamps/");
       strcpy(datafitsname_filter2,"/users/welikala/data/GOODS/ACS/stamps/");
       strcpy(datafitsname_filter3,"/users/welikala/data/GOODS/ACS/stamps/");
       strcpy(datafitsname_filter4,"/users/welikala/data/GOODS/ACS/stamps/");
       strcpy(whtfitsname_filter1,"/users/welikala/data/GOODS/ACS/stamps/");
       strcpy(whtfitsname_filter2,"/users/welikala/data/GOODS/ACS/stamps/");
       strcpy(whtfitsname_filter3,"/users/welikala/data/GOODS/ACS/stamps/");
       strcpy(whtfitsname_filter4,"/users/welikala/data/GOODS/ACS/stamps/");

       strcat(datafitsname_filter1,objidstring);
       strcat(datafitsname_filter1,"_b.fits");
       printf("%s\n",datafitsname_filter1);
       strcat(datafitsname_filter2,objidstring);
       strcat(datafitsname_filter2,"_v.fits");
       printf("%s\n",datafitsname_filter2);
       strcat(datafitsname_filter3,objidstring);
       strcat(datafitsname_filter3,"_i.fits");
       printf("%s\n",datafitsname_filter3);
       strcat(datafitsname_filter4,objidstring);
       strcat(datafitsname_filter4,"_z.fits");
       printf("%s\n",datafitsname_filter4);

       strcat(whtfitsname_filter1,objidstring);
       strcat(whtfitsname_filter1,"_wht_b.fits");
       strcat(whtfitsname_filter2,objidstring);
       strcat(whtfitsname_filter2,"_wht_v.fits");
       strcat(whtfitsname_filter3,objidstring);
       strcat(whtfitsname_filter3,"_wht_i.fits");
       strcat(whtfitsname_filter4,objidstring);
       strcat(whtfitsname_filter4,"_wht_z.fits");

       strcpy(datafitsname_filter5,"/users/welikala/data/GOODS/WFC3/f105w/Ystamp");
       strcpy(whtfitsname_filter5,"/users/welikala/data/GOODS/WFC3/f105w/Ystampwht");
       strcpy(datafitsname_filter6,"/users/welikala/data/GOODS/WFC3/f125w/Jstamp");
       strcpy(whtfitsname_filter6,"/users/welikala/data/GOODS/WFC3/f125w/Jstampwht");
       strcpy(datafitsname_filter7,"/users/welikala/data/GOODS/WFC3/f160w/Hstamp");
       strcpy(whtfitsname_filter7,"/users/welikala/data/GOODS/WFC3/f160w/Hstampwht");

       strcat(datafitsname_filter5,objidstring);
       strcat(datafitsname_filter5,"_0.fits");
       strcat(datafitsname_filter6,objidstring);
       strcat(datafitsname_filter6,"_0.fits");
       strcat(datafitsname_filter7,objidstring);
       strcat(datafitsname_filter7,"_0.fits");
       strcat(whtfitsname_filter5,objidstring);
       strcat(whtfitsname_filter5,"_0.fits");
       strcat(whtfitsname_filter6,objidstring);
       strcat(whtfitsname_filter6,"_0.fits");
       strcat(whtfitsname_filter7,objidstring);
       strcat(whtfitsname_filter7,"_0.fits");


	/* Put check if fits files are readable, otherwise skip to next object?*/

       readimage(datafitsname_filter1,countrate_filter1);
       readimage(datafitsname_filter2,countrate_filter2);
       readimage(datafitsname_filter3,countrate_filter3);
       readimage(datafitsname_filter4,countrate_filter4);
       readimage(datafitsname_filter5,countrate_filter5);
       readimage(datafitsname_filter6,countrate_filter6);
       readimage(datafitsname_filter7,countrate_filter7);

       readimage(whtfitsname_filter1, wht_filter1);
       readimage(whtfitsname_filter2, wht_filter2);
       readimage(whtfitsname_filter3, wht_filter3);
       readimage(whtfitsname_filter4, wht_filter4);
       readimage(whtfitsname_filter5, wht_filter5);
       readimage(whtfitsname_filter6, wht_filter6);
       readimage(whtfitsname_filter7, wht_filter7);

       /* Calibrate the data to units of fluxes in nJy */

       for (pix=0;pix<npixels_acs;pix++)
         {
           simdata_filter1[pix] = countrate_filter1[pix] * pow(10,(48.6+zeropoint_filter1)/(-2.5)) * 1e23 * 1e9 ;
           simdata_filter2[pix] = countrate_filter2[pix] * pow(10,(48.6+zeropoint_filter2)/(-2.5)) * 1e23 * 1e9;
	   simdata_filter3[pix] = countrate_filter3[pix] * pow(10,(48.6+zeropoint_filter3)/(-2.5)) * 1e23 * 1e9 ;
           simdata_filter4[pix] = countrate_filter4[pix] * pow(10,(48.6+zeropoint_filter4)/(-2.5)) * 1e23 * 1e9;
	   simrms_filter1[pix] = sqrt( (pow(corr_correlations_filter1,2)/wht_filter1[pix]) + (countrate_filter1[pix]/exposuretime_filter1) ) * pow(10,(48.6+zeropoint_filter1)/(-2.5)) * 1e23 * 1e9 ;
	   simrms_filter2[pix] = sqrt( (pow(corr_correlations_filter2,2)/wht_filter2[pix]) + (countrate_filter2[pix]/exposuretime_filter2) ) * pow(10,(48.6+zeropoint_filter2)/(-2.5)) * 1e23 * 1e9 ;
	   simrms_filter3[pix] = sqrt( (pow(corr_correlations_filter3,2)/wht_filter3[pix]) + (countrate_filter3[pix]/exposuretime_filter3) ) * pow(10,(48.6+zeropoint_filter3)/(-2.5)) * 1e23 * 1e9 ;
	   simrms_filter4[pix] = sqrt( (pow(corr_correlations_filter4,2)/wht_filter4[pix]) + (countrate_filter4[pix]/exposuretime_filter4) ) * pow(10,(48.6+zeropoint_filter4)/(-2.5)) * 1e23 * 1e9 ;
	 }

       for (pix=0;pix<npixels_wfc3;pix++)
         {
	   simdata_filter5[pix] = countrate_filter5[pix] * pow(10,(48.6+zeropoint_filter5)/(-2.5)) * 1e23 * 1e9;
	   simdata_filter6[pix] = countrate_filter6[pix] * pow(10,(48.6+zeropoint_filter6)/(-2.5)) * 1e23 * 1e9 ;
           simdata_filter7[pix] = countrate_filter7[pix] * pow(10,(48.6+zeropoint_filter7)/(-2.5)) * 1e23 * 1e9;
	   simrms_filter5[pix] = sqrt( (pow(corr_correlations_filter5,2)/wht_filter5[pix]) + (countrate_filter5[pix]/exposuretime_filter5) ) * pow(10,(48.6+zeropoint_filter5)/(-2.5)) * 1e23 * 1e9 ;
	   simrms_filter6[pix] = sqrt( (pow(corr_correlations_filter6,2)/wht_filter6[pix]) + (countrate_filter6[pix]/exposuretime_filter6) ) * pow(10,(48.6+zeropoint_filter6)/(-2.5)) * 1e23 * 1e9 ;
	   simrms_filter7[pix] = sqrt( (pow(corr_correlations_filter7,2)/wht_filter7[pix]) + (countrate_filter7[pix]/exposuretime_filter7) ) * pow(10,(48.6+zeropoint_filter7)/(-2.5)) * 1e23 * 1e9 ;

	 }


       chi2_best_global = 9e99;

       sumflux_data_inside_filter1=0.0;
       sumflux_data_outside_filter1=0.0;
       sumflux_data_galaxy_filter1=0.0;
       sumrms_data_inside_filter1=0.0;
       sumrms_data_outside_filter1=0.0;

       sumflux_data_inside_filter2=0.0;
       sumflux_data_outside_filter2=0.0;
       sumflux_data_galaxy_filter2=0.0;
       sumrms_data_inside_filter2=0.0;
       sumrms_data_outside_filter2=0.0;

       sumflux_data_inside_filter3=0.0;
       sumflux_data_outside_filter3=0.0;
       sumflux_data_galaxy_filter3=0.0;
       sumrms_data_inside_filter3=0.0;
       sumrms_data_outside_filter3=0.0;

       sumflux_data_inside_filter4=0.0;
       sumflux_data_outside_filter4=0.0;
       sumflux_data_galaxy_filter4=0.0;
       sumrms_data_inside_filter4=0.0;
       sumrms_data_outside_filter4=0.0;

       sumflux_data_inside_filter5=0.0;
       sumflux_data_outside_filter5=0.0;
       sumflux_data_galaxy_filter5=0.0;
       sumrms_data_inside_filter5=0.0;
       sumrms_data_outside_filter5=0.0;

       sumflux_data_inside_filter6=0.0;
       sumflux_data_outside_filter6=0.0;
       sumflux_data_galaxy_filter6=0.0;
       sumrms_data_inside_filter6=0.0;
       sumrms_data_outside_filter6=0.0;

       sumflux_data_inside_filter7=0.0;
       sumflux_data_outside_filter7=0.0;
       sumflux_data_galaxy_filter7=0.0;
       sumrms_data_inside_filter7=0.0;
       sumrms_data_outside_filter7=0.0;

       sumflux_model_inside_filter1=0.0;
       sumflux_model_outside_filter1=0.0;
       sumflux_model_galaxy_filter1=0.0;
       sumflux_bulge_filter1 = 0.0;
       sumflux_disk_filter1 = 0.0;

       sumflux_model_inside_filter2=0.0;
       sumflux_model_outside_filter2=0.0;
       sumflux_model_galaxy_filter2=0.0;
       sumflux_bulge_filter2 = 0.0;
       sumflux_disk_filter2 = 0.0;

       sumflux_model_inside_filter3=0.0;
       sumflux_model_outside_filter3=0.0;
       sumflux_model_galaxy_filter3=0.0;
       sumflux_bulge_filter3 = 0.0;
       sumflux_disk_filter3 = 0.0;

       sumflux_model_inside_filter4=0.0;
       sumflux_model_outside_filter4=0.0;
       sumflux_model_galaxy_filter4=0.0;
       sumflux_bulge_filter4 = 0.0;
       sumflux_disk_filter4 = 0.0;

       sumflux_model_inside_filter5=0.0;
       sumflux_model_outside_filter5=0.0;
       sumflux_model_galaxy_filter5=0.0;
       sumflux_bulge_filter5 = 0.0;
       sumflux_disk_filter5 = 0.0;

       sumflux_model_inside_filter6=0.0;
       sumflux_model_outside_filter6=0.0;
       sumflux_model_galaxy_filter6=0.0;
       sumflux_bulge_filter6 = 0.0;
       sumflux_disk_filter6 = 0.0;

       sumflux_model_inside_filter7=0.0;
       sumflux_model_outside_filter7=0.0;
       sumflux_model_galaxy_filter7=0.0;
       sumflux_bulge_filter7 = 0.0;
       sumflux_disk_filter7 = 0.0;


  /* ---------------------------------------------------------------------------------------------------- */
  /* FITTING */
  /* Loop over each morphology type and do joint fit to model and disk and get best fit model parameters   */
  /*------------------------------------------------------------------------------------------------------ */

  /* initialize */

 for (ipoints=0;ipoints<npoints_disk;ipoints++)
    {
      x_disk[ipoints]=0.0;
      y_disk[ipoints]=0.0;
      z_disk[ipoints]=0.0;
      I_xyz_disk_filter1[ipoints]=0.0;
    }
  for (ipoints=0;ipoints<npoints_bulge;ipoints++)
    {
      x_bulge[ipoints]=0.0;
      y_bulge[ipoints]=0.0;
      z_bulge[ipoints]=0.0;
      I_xyz_bulge_filter1[ipoints]=0.0;
    }


  for (itype=0;itype<ntypes;itype++)
    {

      /* Define colours for each morphology type for each galaxy  */


      for (ipoints=0;ipoints<npoints_disk;ipoints++)
	{
	  x_disk[ipoints]=x_disk_type[ipoints][itype];
	  y_disk[ipoints]=y_disk_type[ipoints][itype];
	  z_disk[ipoints]=z_disk_type[ipoints][itype];
	  I_xyz_disk_filter1[ipoints]=I_xyz_disk_filter1_type[ipoints][itype];
	}

      for (ipoints=0;ipoints<npoints_bulge;ipoints++)
	{
	  x_bulge[ipoints]=x_bulge_type[ipoints][itype];
	  y_bulge[ipoints]=y_bulge_type[ipoints][itype];
	  z_bulge[ipoints]=z_bulge_type[ipoints][itype];
	  I_xyz_bulge_filter1[ipoints]=I_xyz_bulge_filter1_type[ipoints][itype];
	}



      /* Do fitting in filter 1 */

      scaling=0.5; //0.5
      e1 = -0.5;  //-0.5
      e2= 0.5;   //0.5
      xcentre = 0.0; //0
      ycentre = 0.0;   //

      scaling_step=0.5;
      e1_step = 0.5;
      e2_step=0.5;
      xcentre_step=0.5;
      ycentre_step =0.5;

      printf("ITYPE: %d\n",itype);


      chi2_best_type=fitjointmodel(scaling,e1,e2,xcentre,ycentre,scaling_step,e1_step,e2_step,xcentre_step,ycentre_step,rvals1);


      scaling_best = rvals1[0];
      e1_best = rvals1[1];
      e2_best = rvals1[2];
      xcentre_best=rvals1[3];
      ycentre_best=rvals1[4];

      scaling_step_best =0.05;
      e1_step_best=0.05;
      e2_step_best=0.05;
      xcentre_step_best=0.05;
      ycentre_step_best=0.05;


      chi2_best_type=fitjointmodel(scaling_best,e1_best, e2_best, xcentre_best, ycentre_best, scaling_step_best,e1_step_best,e2_step_best,xcentre_step_best, ycentre_step_best, rvals2);

      scaling_best=rvals2[0];
      e1_best = rvals2[1];
      e2_best = rvals2[2];
      xcentre_best = rvals2[3];
      ycentre_best = rvals2[4];

      /* Third iteration of simplex - skip here */

      /*
      scaling_step_best=0.05;
      e1_step_best=0.05;
      e2_step_best=0.05;
      xcentre_step_best=0.05;
      ycentre_step_best=0.05;

      chi2_best_type=fitjointmodel_filter1(scaling_best,e1_best,e2_best, xcentre_best, ycentre_best, scaling_step_best,e1_step_best,e2_step_best, xcentre_step_best, ycentre_step_best, rvals3);

      scaling_best=rvals3[0];
      e1_best = rvals3[1];
      e2_best = rvals3[2];
      xcentre_best = rvals3[3];
      ycentre_best = rvals3[4];
      */

      if (chi2_best_type < chi2_best_global) {
        chi2_best_global = chi2_best_type;
	amp_bulge_global_filter1 = amp_bulge_type_filter1;
	amp_disk_global_filter1 = amp_disk_type_filter1;
	amp_bulge_global_filter2 = amp_bulge_type_filter2;
	amp_disk_global_filter2 = amp_disk_type_filter2;
	amp_bulge_global_filter3 = amp_bulge_type_filter3;
	amp_disk_global_filter3 = amp_disk_type_filter3;
	amp_bulge_global_filter4 = amp_bulge_type_filter4;
	amp_disk_global_filter4 = amp_disk_type_filter4;
	amp_bulge_global_filter5 = amp_bulge_type_filter5;
	amp_disk_global_filter5 = amp_disk_type_filter5;
	amp_bulge_global_filter6 = amp_bulge_type_filter6;
	amp_disk_global_filter6 = amp_disk_type_filter6;
	amp_bulge_global_filter7 = amp_bulge_type_filter7;
	amp_disk_global_filter7 = amp_disk_type_filter7;

        itype_best_global = itype;
        scaling_best_global = scaling_best;
        e1_best_global = e1_best;
        e2_best_global = e2_best;
        xcentre_best_global = xcentre_best;
        ycentre_best_global = ycentre_best;
      }

    } // loop over galaxy models


      /* Find centre of image in ACS pixel coordinates */

      pixelsize_kpc_acs = pixelsize_arcsec_acs * arcsecstokpc;
      xmin = (-ncols_tot_acs/2.0)*pixelsize_kpc_acs;  //kpc
      ymin = (-nrows_tot_acs/2.0)*pixelsize_kpc_acs;  //kpc
      x0_bestmodel_acs = (int) ((xcentre_best_global-xmin)/pixelsize_kpc_acs) ; // ACS pixel in x direction
      y0_bestmodel_acs = (int) ((ycentre_best_global-ymin)/pixelsize_kpc_acs) ;  // ACS pixel in y direction
      pix_bestmodel_acs = (ncols_tot_acs * y0_bestmodel_acs) + x0_bestmodel_acs;

      /* Find centre of image in WFC3 pixel coordinates */

      pixelsize_kpc_wfc3 = pixelsize_arcsec_wfc3 * arcsecstokpc;
      xmin = (-ncols_tot_wfc3/2.0)*pixelsize_kpc_wfc3;  //kpc
      ymin = (-nrows_tot_wfc3/2.0)*pixelsize_kpc_wfc3;  //kpc
      x0_bestmodel_wfc3 = (int) ((xcentre_best_global-xmin)/pixelsize_kpc_wfc3) ; // WFC3 pixel in x direction
      y0_bestmodel_wfc3 = (int) ((ycentre_best_global-ymin)/pixelsize_kpc_wfc3) ;  // WFC3 pixel in y direction
      pix_bestmodel_wfc3 = (ncols_tot_wfc3 * y0_bestmodel_wfc3) + x0_bestmodel_wfc3;


      /*
      theta_rad_best = asin(sqrt(e1_best*e1_best + e2_best*e2_best));
      posangle_rad_best  =  0.5 * atan2(e2_best,e1_best) ;
      theta_deg_best = theta_rad_best * (180.0/M_PI);
      posangle_deg_best = posangle_rad_best * (180.0/M_PI);
      */

      /*----------------------------------- */
      /*----------------------------------- */
      /* Make best-fit models */
      /* ------------------------------------ */
      /*----------------------------------- */

      for (ipoints=0;ipoints<npoints_disk;ipoints++)
        {
          x_disk[ipoints]=x_disk_type[ipoints][itype_best_global];
          y_disk[ipoints]=y_disk_type[ipoints][itype_best_global];
          z_disk[ipoints]=z_disk_type[ipoints][itype_best_global];
          I_xyz_disk_filter1[ipoints]=I_xyz_disk_filter1_type[ipoints][itype_best_global];
        }

      for (ipoints=0;ipoints<npoints_bulge;ipoints++)
        {
          x_bulge[ipoints]=x_bulge_type[ipoints][itype_best_global];
          y_bulge[ipoints]=y_bulge_type[ipoints][itype_best_global];
          z_bulge[ipoints]=z_bulge_type[ipoints][itype_best_global];
          I_xyz_bulge_filter1[ipoints]=I_xyz_bulge_filter1_type[ipoints][itype_best_global];
        }

      /* ACS bands */

      /* filter 1*/
      make2Dgalaxymodelunconvolved_acs(0,0,scaling_best_global, e1_best_global, e2_best_global, xcentre_best_global,ycentre_best_global, arcsecstokpc,1.0);
      for (pix=0;pix<npixels_acs;pix++)
	model_unconv_bulge_filter1[pix] = model_unconv_acs[pix];

      make2Dgalaxymodelunconvolved_acs(1,0, scaling_best_global, e1_best_global, e2_best_global, xcentre_best_global,ycentre_best_global, arcsecstokpc, 1.0);
      for (pix=0;pix<npixels_acs;pix++)
	model_unconv_disk_filter1[pix] = model_unconv_acs[pix];

      make2Dgalaxymodelconvolved_acs(0,0,scaling_best_global, e1_best_global, e2_best_global, xcentre_best_global,ycentre_best_global, arcsecstokpc, 1.0);
      for (pix=0;pix<npixels_acs;pix++)
        model_conv_bulge_filter1[pix] = model_conv_acs[pix];

      make2Dgalaxymodelconvolved_acs(1,0, scaling_best_global, e1_best_global, e2_best_global,xcentre_best_global,ycentre_best_global,  arcsecstokpc, 1.0);
      for (pix=0;pix<npixels_acs;pix++)
        model_conv_disk_filter1[pix] = model_conv_acs[pix];

      /* filter 2 */
     make2Dgalaxymodelunconvolved_acs(0,1,scaling_best_global, e1_best_global, e2_best_global,xcentre_best_global,ycentre_best_global,arcsecstokpc, 1.0);
     for (pix=0;pix<npixels_acs;pix++)
       model_unconv_bulge_filter2[pix] = model_unconv_acs[pix];

      make2Dgalaxymodelunconvolved_acs(1,1,scaling_best_global, e1_best_global, e2_best_global,xcentre_best_global,ycentre_best_global, arcsecstokpc, 1.0);
      for (pix=0;pix<npixels_acs;pix++)
	model_unconv_disk_filter2[pix] = model_unconv_acs[pix];

      make2Dgalaxymodelconvolved_acs(0,1,scaling_best_global, e1_best_global, e2_best_global,xcentre_best_global, ycentre_best_global, arcsecstokpc, 1.0);
      for (pix=0;pix<npixels_acs;pix++)
	model_conv_bulge_filter2[pix] = model_conv_acs[pix];

      make2Dgalaxymodelconvolved_acs(1,1, scaling_best_global, e1_best_global, e2_best_global,xcentre_best_global,ycentre_best_global, arcsecstokpc, 1.0);
      for (pix=0;pix<npixels_acs;pix++)
	model_conv_disk_filter2[pix] = model_conv_acs[pix];

      /* filter 3 */

     make2Dgalaxymodelunconvolved_acs(0,2,scaling_best_global, e1_best_global, e2_best_global,xcentre_best_global,ycentre_best_global, arcsecstokpc, 1.0);
     for (pix=0;pix<npixels_acs;pix++)
       model_unconv_bulge_filter3[pix] = model_unconv_acs[pix];

      make2Dgalaxymodelunconvolved_acs(1,2,scaling_best_global, e1_best_global, e2_best_global,xcentre_best_global,ycentre_best_global, arcsecstokpc, 1.0);
      for (pix=0;pix<npixels_acs;pix++)
	model_unconv_disk_filter3[pix] = model_unconv_acs[pix];

      make2Dgalaxymodelconvolved_acs(0,2,scaling_best_global, e1_best_global, e2_best_global,xcentre_best_global, ycentre_best_global, arcsecstokpc, 1.0);
      for (pix=0;pix<npixels_acs;pix++)
	model_conv_bulge_filter3[pix] = model_conv_acs[pix];

      make2Dgalaxymodelconvolved_acs(1,2, scaling_best_global, e1_best_global, e2_best_global,xcentre_best_global,ycentre_best_global, arcsecstokpc, 1.0);
      for (pix=0;pix<npixels_acs;pix++)
	model_conv_disk_filter3[pix] = model_conv_acs[pix];

      /* Make best fit models in filter 4 */

      make2Dgalaxymodelunconvolved_acs(0,3,scaling_best_global, e1_best_global, e2_best_global,xcentre_best_global,ycentre_best_global,  arcsecstokpc, 1.0);
     for (pix=0;pix<npixels_acs;pix++)
       model_unconv_bulge_filter4[pix] = model_unconv_acs[pix];

      make2Dgalaxymodelunconvolved_acs(1,3,scaling_best_global, e1_best_global, e2_best_global,xcentre_best_global,ycentre_best_global, arcsecstokpc, 1.0);
      for (pix=0;pix<npixels_acs;pix++)
	model_unconv_disk_filter4[pix] = model_unconv_acs[pix];

      make2Dgalaxymodelconvolved_acs(0,3,scaling_best_global, e1_best_global, e2_best_global,xcentre_best_global, ycentre_best_global, arcsecstokpc, 1.0);
      for (pix=0;pix<npixels_acs;pix++)
	model_conv_bulge_filter4[pix] = model_conv_acs[pix];

      make2Dgalaxymodelconvolved_acs(1,3, scaling_best_global, e1_best_global, e2_best_global,xcentre_best_global,ycentre_best_global, arcsecstokpc, 1.0);
      for (pix=0;pix<npixels_acs;pix++)
	model_conv_disk_filter4[pix] = model_conv_acs[pix];


      /* WFC3 bands */

      /* filter 5 */

      make2Dgalaxymodelunconvolved_wfc3(0,4,scaling_best_global, e1_best_global, e2_best_global,xcentre_best_global,ycentre_best_global,  arcsecstokpc, 1.0);
     for (pix=0;pix<npixels_wfc3;pix++)
       model_unconv_bulge_filter5[pix] = model_unconv_wfc3[pix];

      make2Dgalaxymodelunconvolved_wfc3(1,4,scaling_best_global, e1_best_global, e2_best_global,xcentre_best_global,ycentre_best_global, arcsecstokpc, 1.0);
      for (pix=0;pix<npixels_wfc3;pix++)
	model_unconv_disk_filter5[pix] = model_unconv_wfc3[pix];

      make2Dgalaxymodelconvolved_wfc3(0,4,scaling_best_global, e1_best_global, e2_best_global,xcentre_best_global, ycentre_best_global,arcsecstokpc, 1.0);
      for (pix=0;pix<npixels_wfc3;pix++)
	model_conv_bulge_filter5[pix] = model_conv_wfc3[pix];

      make2Dgalaxymodelconvolved_wfc3(1,4, scaling_best_global, e1_best_global, e2_best_global,xcentre_best_global,ycentre_best_global,arcsecstokpc, 1.0);
      for (pix=0;pix<npixels_wfc3;pix++)
	model_conv_disk_filter5[pix] = model_conv_wfc3[pix];


      /* filter 6 */

        make2Dgalaxymodelunconvolved_wfc3(0,5,scaling_best_global, e1_best_global, e2_best_global,xcentre_best_global,ycentre_best_global, arcsecstokpc, 1.0);
     for (pix=0;pix<npixels_wfc3;pix++)
       model_unconv_bulge_filter6[pix] = model_unconv_wfc3[pix];

      make2Dgalaxymodelunconvolved_wfc3(1,5,scaling_best_global, e1_best_global, e2_best_global,xcentre_best_global,ycentre_best_global, arcsecstokpc, 1.0);
      for (pix=0;pix<npixels_wfc3;pix++)
	model_unconv_disk_filter6[pix] = model_unconv_wfc3[pix];

      make2Dgalaxymodelconvolved_wfc3(0,5,scaling_best_global, e1_best_global, e2_best_global,xcentre_best_global, ycentre_best_global, arcsecstokpc, 1.0);
      for (pix=0;pix<npixels_wfc3;pix++)
	model_conv_bulge_filter6[pix] = model_conv_wfc3[pix];

      make2Dgalaxymodelconvolved_wfc3(1,5, scaling_best_global, e1_best_global, e2_best_global,xcentre_best_global,ycentre_best_global, arcsecstokpc, 1.0);
      for (pix=0;pix<npixels_wfc3;pix++)
	model_conv_disk_filter6[pix] = model_conv_wfc3[pix];



      /* filter 7 */

        make2Dgalaxymodelunconvolved_wfc3(0,6,scaling_best_global, e1_best_global, e2_best_global,xcentre_best_global,ycentre_best_global, arcsecstokpc, 1.0);
     for (pix=0;pix<npixels_wfc3;pix++)
       model_unconv_bulge_filter7[pix] = model_unconv_wfc3[pix];

      make2Dgalaxymodelunconvolved_wfc3(1,6,scaling_best_global, e1_best_global, e2_best_global,xcentre_best_global,ycentre_best_global,arcsecstokpc, 1.0);
      for (pix=0;pix<npixels_wfc3;pix++)
	model_unconv_disk_filter7[pix] = model_unconv_wfc3[pix];

      make2Dgalaxymodelconvolved_wfc3(0,6,scaling_best_global, e1_best_global, e2_best_global,xcentre_best_global, ycentre_best_global, arcsecstokpc, 1.0);
      for (pix=0;pix<npixels_wfc3;pix++)
	model_conv_bulge_filter7[pix] = model_conv_wfc3[pix];

      make2Dgalaxymodelconvolved_wfc3(1,6, scaling_best_global, e1_best_global, e2_best_global,xcentre_best_global,ycentre_best_global, arcsecstokpc, 1.0);
      for (pix=0;pix<npixels_wfc3;pix++)
	model_conv_disk_filter7[pix] = model_conv_wfc3[pix];


      /* --------------------------------------------------------------------------------------- */
      /* Joint fitting for morphology and colours (amplitudes) */
      /* --------------------------------------------------------------------------------------- */

      for (j=0; j<npixels_acs; j++)
	{
	  model_disk_filter1[j]= model_conv_disk_filter1[j]/simrms_filter1[j];
	  model_bulge_filter1[j]= model_conv_bulge_filter1[j]/simrms_filter1[j];
	  model_disk_filter2[j]= model_conv_disk_filter2[j]/simrms_filter2[j];
	  model_bulge_filter2[j]= model_conv_bulge_filter2[j]/simrms_filter2[j];
	  model_disk_filter3[j]= model_conv_disk_filter3[j]/simrms_filter3[j];
	  model_bulge_filter3[j]= model_conv_bulge_filter3[j]/simrms_filter3[j];
	  model_disk_filter4[j]= model_conv_disk_filter4[j]/simrms_filter4[j];
	  model_bulge_filter4[j]= model_conv_bulge_filter4[j]/simrms_filter4[j];

	  data_filter1[j] = simdata_filter1[j]/simrms_filter1[j];
	  data_filter2[j] = simdata_filter2[j]/simrms_filter2[j];
	  data_filter3[j] = simdata_filter3[j]/simrms_filter3[j];
	  data_filter4[j] = simdata_filter4[j]/simrms_filter4[j];

	}


  for (j=0; j<npixels_wfc3; j++)
	{

	  model_disk_filter5[j]= model_conv_disk_filter5[j]/simrms_filter5[j];
	  model_bulge_filter5[j]= model_conv_bulge_filter5[j]/simrms_filter5[j];
	  model_disk_filter6[j]= model_conv_disk_filter6[j]/simrms_filter6[j];
	  model_bulge_filter6[j]= model_conv_bulge_filter6[j]/simrms_filter6[j];
	  model_disk_filter7[j]= model_conv_disk_filter7[j]/simrms_filter7[j];
	  model_bulge_filter7[j]= model_conv_bulge_filter7[j]/simrms_filter7[j];

	  data_filter5[j] = simdata_filter5[j]/simrms_filter5[j];
	  data_filter6[j] = simdata_filter6[j]/simrms_filter6[j];
	  data_filter7[j] = simdata_filter7[j]/simrms_filter7[j];

	}




      /* ----------------------------------------------------- */
      /* Get amplitudes for filter 1 and filter 2 separately  */
      /* ---------------------------------------------------- */
      amp_etype_filter1=0.0;
      amp_bulge_filter1=0.0;
      amp_disk_filter1=0.0;
      amp_etype_filter2=0.0;
      amp_bulge_filter2=0.0;
      amp_disk_filter2=0.0;

      amp_etype_filter3=0.0;
      amp_bulge_filter3=0.0;
      amp_disk_filter3=0.0;
      amp_etype_filter4=0.0;
      amp_bulge_filter4=0.0;
      amp_disk_filter4=0.0;

      amp_etype_filter5=0.0;
      amp_bulge_filter5=0.0;
      amp_disk_filter5=0.0;
      amp_etype_filter6=0.0;
      amp_bulge_filter6=0.0;
      amp_disk_filter6=0.0;
      amp_etype_filter7=0.0;
      amp_bulge_filter7=0.0;
      amp_disk_filter7=0.0;

      denom=0.0;

      /* Compute likelihood for E-type separately - bulge only */

      if (itype_best_global==0)
	{

	  sumx=sumy=sumxx=sumyy=sumxy=0.0;
	  lval = 0.0;

	  for (j=0; j<npixels_acs; j++)
	    {

	      sumx += model_bulge_filter1[j];
	      sumy += data_filter1[j];
	      sumxx += pow(model_bulge_filter1[j],2);
	      sumyy += pow(data_filter1[j],2);
	      sumxy += data_filter1[j]*model_bulge_filter1[j];
	    }

	  denom = sumxx;

	  if (denom != 0.)
	    {
	      amp_etype_filter1 = sumxy/sumxx;
	      amp_bulge_filter1=amp_etype_filter1;

	      for (j=0;j<npixels_acs;j++)
		{
		  lval += pow((data_filter1[j]- (amp_etype_filter1 * model_bulge_filter1[j])),2);

		}

	    }
	  else
	    {
	      lval += 9e99;
	    }
	  sumx=sumy=sumxx=sumyy=sumxy=0.0;

	  //printf("amp_bulge (E-TYPE) lval: %6.3e %6.3e\n",amp_etype_filter1, lval);


	  /* 2nd filter */

          for (j=0; j<npixels_acs; j++)
            {

              sumx += model_bulge_filter2[j];
              sumy += data_filter2[j];
              sumxx += pow(model_bulge_filter2[j],2);
              sumyy += pow(data_filter2[j],2);
              sumxy += data_filter2[j]*model_bulge_filter2[j];
            }

          denom = sumxx;

          if (denom != 0.)
            {
              amp_etype_filter2 = sumxy/sumxx;
              amp_bulge_filter2=amp_etype_filter2;
              for (j=0;j<npixels_acs;j++)
                {
                  lval += pow((data_filter2[j]- (amp_etype_filter2 * model_bulge_filter2[j])),2);

                }

            }
          else
            {
              lval += 9e99;
            }
          sumx=sumy=sumxx=sumyy=sumxy=0.0;

	  /* 3rd filter */
	  for (j=0; j<npixels_acs; j++)
            {

              sumx += model_bulge_filter3[j];
              sumy += data_filter3[j];
              sumxx += pow(model_bulge_filter3[j],2);
              sumyy += pow(data_filter3[j],2);
              sumxy += data_filter3[j]*model_bulge_filter3[j];
            }

          denom = sumxx;

          if (denom != 0.)
            {
              amp_etype_filter3 = sumxy/sumxx;
              amp_bulge_filter3=amp_etype_filter3;
              for (j=0;j<npixels_acs;j++)
                {
                  lval += pow((data_filter3[j]- (amp_etype_filter3 * model_bulge_filter3[j])),2);

                }

            }
          else
            {
              lval += 9e99;
            }
          sumx=sumy=sumxx=sumyy=sumxy=0.0;

	  /* 4th filter */
	  for (j=0; j<npixels_acs; j++)
            {

              sumx += model_bulge_filter4[j];
              sumy += data_filter4[j];
              sumxx += pow(model_bulge_filter4[j],2);
              sumyy += pow(data_filter4[j],2);
              sumxy += data_filter4[j]*model_bulge_filter4[j];
            }

          denom = sumxx;

          if (denom != 0.)
            {
              amp_etype_filter4 = sumxy/sumxx;
              amp_bulge_filter4=amp_etype_filter4;
              for (j=0;j<npixels_acs;j++)
                {
                  lval += pow((data_filter4[j]- (amp_etype_filter4 * model_bulge_filter4[j])),2);

                }

            }
          else
            {
              lval += 9e99;
            }
          sumx=sumy=sumxx=sumyy=sumxy=0.0;

	  /* filter 5 */

	  for (j=0; j<npixels_wfc3; j++)
            {

              sumx += model_bulge_filter5[j];
              sumy += data_filter5[j];
              sumxx += pow(model_bulge_filter5[j],2);
              sumyy += pow(data_filter5[j],2);
              sumxy += data_filter5[j]*model_bulge_filter5[j];
            }

          denom = sumxx;

          if (denom != 0.)
            {
              amp_etype_filter5 = sumxy/sumxx;
              amp_bulge_filter5=amp_etype_filter5;
              for (j=0;j<npixels_wfc3;j++)
                {
                  lval += pow((data_filter5[j]- (amp_etype_filter5 * model_bulge_filter5[j])),2);

                }

            }
          else
            {
              lval += 9e99;
            }
          sumx=sumy=sumxx=sumyy=sumxy=0.0;


	  /* filter 6 */

	  for (j=0; j<npixels_wfc3; j++)
            {

              sumx += model_bulge_filter6[j];
              sumy += data_filter6[j];
              sumxx += pow(model_bulge_filter6[j],2);
              sumyy += pow(data_filter6[j],2);
              sumxy += data_filter6[j]*model_bulge_filter6[j];
            }

          denom = sumxx;

          if (denom != 0.)
            {
              amp_etype_filter6 = sumxy/sumxx;
              amp_bulge_filter6=amp_etype_filter6;
              for (j=0;j<npixels_wfc3;j++)
                {
                  lval += pow((data_filter6[j]- (amp_etype_filter6 * model_bulge_filter6[j])),2);

                }

            }
          else
            {
              lval += 9e99;
            }
          sumx=sumy=sumxx=sumyy=sumxy=0.0;


	   /* filter 7 */

	  for (j=0; j<npixels_wfc3; j++)
            {

              sumx += model_bulge_filter7[j];
              sumy += data_filter7[j];
              sumxx += pow(model_bulge_filter7[j],2);
              sumyy += pow(data_filter7[j],2);
              sumxy += data_filter7[j]*model_bulge_filter7[j];
            }

          denom = sumxx;

          if (denom != 0.)
            {
              amp_etype_filter7 = sumxy/sumxx;
              amp_bulge_filter7=amp_etype_filter7;
              for (j=0;j<npixels_wfc3;j++)
                {
                  lval += pow((data_filter7[j]- (amp_etype_filter7 * model_bulge_filter7[j])),2);

                }

            }
          else
            {
              lval += 9e99;
            }
          sumx=sumy=sumxx=sumyy=sumxy=0.0;



          //printf("E-type amp, chi2 %lf %lf\n",amp_etype,lval);
	  /* apply penalties for unphysical ranges */
	  if (amp_etype_filter1<0.0) {penalty_amp_etype_filter1=1e7;}
	  if (amp_etype_filter2<0.0) {penalty_amp_etype_filter2=1e7;}
	  if (amp_etype_filter3<0.0) {penalty_amp_etype_filter3=1e7;}
	  if (amp_etype_filter4<0.0) {penalty_amp_etype_filter4=1e7;}
	  if (amp_etype_filter5<0.0) {penalty_amp_etype_filter5=1e7;}
	  if (amp_etype_filter6<0.0) {penalty_amp_etype_filter6=1e7;}
	  if (amp_etype_filter7<0.0) {penalty_amp_etype_filter7=1e7;}


          lval = lval + penalty_e1e2quad + penalty_scaling + penalty_amp_etype_filter1 + penalty_amp_etype_filter2 + penalty_amp_etype_filter3 + penalty_amp_etype_filter4 + penalty_amp_etype_filter5 + penalty_amp_etype_filter6 + penalty_amp_etype_filter7;

	  /* Set disk amplitudes to zero by definition */
	  amp_disk_filter1 = 0.0;
	  amp_disk_filter2 = 0.0;
	  amp_disk_filter3 = 0.0;
	  amp_disk_filter4 = 0.0;
	  amp_disk_filter5 = 0.0;
	  amp_disk_filter6 = 0.0;
	  amp_disk_filter7 = 0.0;

	  //printf("amp_etype: %lf %lf %lf %lf\n", amp_etype_filter1, amp_etype_filter2, amp_etype_filter3, amp_etype_filter4);

	}



      /* Continue for other types */
      if (itype_best_global>0)
	{

	  lval=sumxx=sumzz=sumyy=sumxy=sumzy=sumxz=0.0;

	  /* 1 filter */
	  for (j=0; j<npixels_acs; j++)
	    {

	      sumxx += pow(model_disk_filter1[j],2);
	      sumzz += pow(model_bulge_filter1[j],2);
	      sumyy += pow(data_filter1[j],2);
	      sumxy += data_filter1[j]*model_disk_filter1[j];
	      sumzy += data_filter1[j]*model_bulge_filter1[j];
	      sumxz += model_disk_filter1[j]*model_bulge_filter1[j];
	    }


	  denom = sumxx*sumzz-sumxz*sumxz;


	  if (denom != 0.0)
	    {

	      amp_disk_filter1 = (sumzz*sumxy-sumxz*sumzy)/denom; //disk
	      amp_bulge_filter1 = (sumxx*sumzy-sumxz*sumxy)/denom;  // bulge

	      for (j=0;j<npixels_acs;j++)
		{
		  lval = lval + pow(data_filter1[j]- ((amp_disk_filter1 * model_disk_filter1[j]) + (amp_bulge_filter1 * model_bulge_filter1[j]) ) ,2);

		}


	    }
	  else
	    {
	      printf("sumxx: %lf\n",sumxx);
	      printf("sumzz: %lf\n",sumzz);
	      printf("sumxz: %lf\n",sumxz);
	      lval += 9e99 ;
	    }

	  sumxx=0.0;
	  sumyy=0.0;
	  sumxy=0.0;
	  sumzy=0.0;
	  sumxz=0.0;
	  sumzz=0.0;

	  /* 2nd filter */
          for (j=0; j<npixels_acs; j++)
            {

              sumxx += pow(model_disk_filter2[j],2);
              sumzz += pow(model_bulge_filter2[j],2);
              sumyy += pow(data_filter2[j],2);
              sumxy += data_filter2[j]*model_disk_filter2[j];
              sumzy += data_filter2[j]*model_bulge_filter2[j];
              sumxz += model_disk_filter2[j]*model_bulge_filter2[j];
            }


          denom = sumxx*sumzz-sumxz*sumxz;

          if (denom != 0.0)
            {

              amp_disk_filter2 = (sumzz*sumxy-sumxz*sumzy)/denom; //disk
	      amp_bulge_filter2 = (sumxx*sumzy-sumxz*sumxy)/denom;  // bulge
              for (j=0;j<npixels_acs;j++)
                {
                  lval = lval + pow(data_filter2[j]- ((amp_disk_filter2 * model_disk_filter2[j]) + (amp_bulge_filter2 * model_bulge_filter2[j]) ) ,2);

                }


            }
          else
            {
              printf("sumxx: %lf\n",sumxx);
              printf("sumzz: %lf\n",sumzz);
              printf("sumxz: %lf\n",sumxz);
              lval += 9e99 ;
            }

          sumxx=0.0;
          sumyy=0.0;
          sumxy=0.0;
          sumzy=0.0;
          sumxz=0.0;
          sumzz=0.0;

	  /* 3rd filter */
	   for (j=0; j<npixels_acs; j++)
            {

              sumxx += pow(model_disk_filter3[j],2);
              sumzz += pow(model_bulge_filter3[j],2);
              sumyy += pow(data_filter3[j],2);
              sumxy += data_filter3[j]*model_disk_filter3[j];
              sumzy += data_filter3[j]*model_bulge_filter3[j];
              sumxz += model_disk_filter3[j]*model_bulge_filter3[j];
            }


          denom = sumxx*sumzz-sumxz*sumxz;

          if (denom != 0.0)
            {

              amp_disk_filter3 = (sumzz*sumxy-sumxz*sumzy)/denom; //disk
	      amp_bulge_filter3 = (sumxx*sumzy-sumxz*sumxy)/denom;  // bulge
              for (j=0;j<npixels_acs;j++)
                {
                  lval = lval + pow(data_filter3[j]- ((amp_disk_filter3 * model_disk_filter3[j]) + (amp_bulge_filter3 * model_bulge_filter3[j]) ) ,2);

                }


            }
          else
            {
              printf("sumxx: %lf\n",sumxx);
              printf("sumzz: %lf\n",sumzz);
              printf("sumxz: %lf\n",sumxz);
              lval += 9e99 ;
            }

          sumxx=0.0;
          sumyy=0.0;
          sumxy=0.0;
          sumzy=0.0;
          sumxz=0.0;
          sumzz=0.0;

	  /* 4th filter */

	   for (j=0; j<npixels_acs; j++)
            {

              sumxx += pow(model_disk_filter4[j],2);
              sumzz += pow(model_bulge_filter4[j],2);
              sumyy += pow(data_filter4[j],2);
              sumxy += data_filter4[j]*model_disk_filter4[j];
              sumzy += data_filter4[j]*model_bulge_filter4[j];
              sumxz += model_disk_filter4[j]*model_bulge_filter4[j];
            }


          denom = sumxx*sumzz-sumxz*sumxz;

          if (denom != 0.0)
            {

              amp_disk_filter4 = (sumzz*sumxy-sumxz*sumzy)/denom; //disk
	      amp_bulge_filter4 = (sumxx*sumzy-sumxz*sumxy)/denom;  // bulge
              for (j=0;j<npixels_acs;j++)
                {
                  lval = lval + pow(data_filter4[j]- ((amp_disk_filter4 * model_disk_filter4[j]) + (amp_bulge_filter4 * model_bulge_filter4[j]) ) ,2);

                }


            }
          else
            {
              printf("sumxx: %lf\n",sumxx);
              printf("sumzz: %lf\n",sumzz);
              printf("sumxz: %lf\n",sumxz);
              lval += 9e99 ;
            }

          sumxx=0.0;
          sumyy=0.0;
          sumxy=0.0;
          sumzy=0.0;
          sumxz=0.0;
          sumzz=0.0;

	   /* 5th filter */

	   for (j=0; j<npixels_wfc3; j++)
            {

              sumxx += pow(model_disk_filter5[j],2);
              sumzz += pow(model_bulge_filter5[j],2);
              sumyy += pow(data_filter5[j],2);
              sumxy += data_filter5[j]*model_disk_filter5[j];
              sumzy += data_filter5[j]*model_bulge_filter5[j];
              sumxz += model_disk_filter5[j]*model_bulge_filter5[j];
            }


          denom = sumxx*sumzz-sumxz*sumxz;

          if (denom != 0.0)
            {

              amp_disk_filter5 = (sumzz*sumxy-sumxz*sumzy)/denom; //disk
	      amp_bulge_filter5 = (sumxx*sumzy-sumxz*sumxy)/denom;  // bulge
              for (j=0;j<npixels_wfc3;j++)
                {
                  lval = lval + pow(data_filter5[j]- ((amp_disk_filter5 * model_disk_filter5[j]) + (amp_bulge_filter5 * model_bulge_filter5[j]) ) ,2);

                }


            }
          else
            {
              printf("sumxx: %lf\n",sumxx);
              printf("sumzz: %lf\n",sumzz);
              printf("sumxz: %lf\n",sumxz);
              lval += 9e99 ;
            }

          sumxx=0.0;
          sumyy=0.0;
          sumxy=0.0;
          sumzy=0.0;
          sumxz=0.0;
          sumzz=0.0;

	   /* 6th filter */

	   for (j=0; j<npixels_wfc3; j++)
            {

              sumxx += pow(model_disk_filter6[j],2);
              sumzz += pow(model_bulge_filter6[j],2);
              sumyy += pow(data_filter6[j],2);
              sumxy += data_filter6[j]*model_disk_filter6[j];
              sumzy += data_filter6[j]*model_bulge_filter6[j];
              sumxz += model_disk_filter6[j]*model_bulge_filter6[j];
            }


          denom = sumxx*sumzz-sumxz*sumxz;

          if (denom != 0.0)
            {

              amp_disk_filter6 = (sumzz*sumxy-sumxz*sumzy)/denom; //disk
	      amp_bulge_filter6 = (sumxx*sumzy-sumxz*sumxy)/denom;  // bulge
              for (j=0;j<npixels_wfc3;j++)
                {
                  lval = lval + pow(data_filter6[j]- ((amp_disk_filter6 * model_disk_filter6[j]) + (amp_bulge_filter6 * model_bulge_filter6[j]) ) ,2);

                }


            }
          else
            {
              printf("sumxx: %lf\n",sumxx);
              printf("sumzz: %lf\n",sumzz);
              printf("sumxz: %lf\n",sumxz);
              lval += 9e99 ;
            }

          sumxx=0.0;
          sumyy=0.0;
          sumxy=0.0;
          sumzy=0.0;
          sumxz=0.0;
          sumzz=0.0;

	   /* 7th filter */

	   for (j=0; j<npixels_wfc3; j++)
            {

              sumxx += pow(model_disk_filter7[j],2);
              sumzz += pow(model_bulge_filter7[j],2);
              sumyy += pow(data_filter7[j],2);
              sumxy += data_filter7[j]*model_disk_filter7[j];
              sumzy += data_filter7[j]*model_bulge_filter7[j];
              sumxz += model_disk_filter7[j]*model_bulge_filter7[j];
            }


          denom = sumxx*sumzz-sumxz*sumxz;

          if (denom != 0.0)
            {

              amp_disk_filter7 = (sumzz*sumxy-sumxz*sumzy)/denom; //disk
	      amp_bulge_filter7 = (sumxx*sumzy-sumxz*sumxy)/denom;  // bulge
              for (j=0;j<npixels_wfc3;j++)
                {
                  lval = lval + pow(data_filter7[j]- ((amp_disk_filter7 * model_disk_filter7[j]) + (amp_bulge_filter7 * model_bulge_filter7[j]) ) ,2);

                }


            }
          else
            {
              printf("sumxx: %lf\n",sumxx);
              printf("sumzz: %lf\n",sumzz);
              printf("sumxz: %lf\n",sumxz);
              lval += 9e99 ;
            }

          sumxx=0.0;
          sumyy=0.0;
          sumxy=0.0;
          sumzy=0.0;
          sumxz=0.0;
          sumzz=0.0;

	  if (amp_bulge_filter1<0.0) {penalty_amp_bulge_filter1=1e7;}
	  if (amp_disk_filter1<0.0) {penalty_amp_disk_filter1=1e7;}
	  if (amp_bulge_filter2<0.0) {penalty_amp_bulge_filter2=1e7;}
	  if (amp_disk_filter2<0.0) {penalty_amp_disk_filter2=1e7;}
	  if (amp_bulge_filter3<0.0) {penalty_amp_bulge_filter3=1e7;}
	  if (amp_disk_filter3<0.0) {penalty_amp_disk_filter3=1e7;}
	  if (amp_bulge_filter4<0.0) {penalty_amp_bulge_filter4=1e7;}
	  if (amp_disk_filter4<0.0) {penalty_amp_disk_filter4=1e7;}

	  if (amp_bulge_filter5<0.0) {penalty_amp_bulge_filter5=1e7;}
	  if (amp_disk_filter5<0.0) {penalty_amp_disk_filter5=1e7;}
	  if (amp_bulge_filter6<0.0) {penalty_amp_bulge_filter6=1e7;}
	  if (amp_disk_filter6<0.0) {penalty_amp_disk_filter6=1e7;}
	  if (amp_bulge_filter7<0.0) {penalty_amp_bulge_filter7=1e7;}
	  if (amp_disk_filter7<0.0) {penalty_amp_disk_filter7=1e7;}

	  lval = lval + penalty_e1e2quad + penalty_scaling + penalty_amp_bulge_filter1 + penalty_amp_disk_filter1 +  penalty_amp_bulge_filter2 + penalty_amp_disk_filter2 + penalty_amp_bulge_filter3 + penalty_amp_disk_filter3 + penalty_amp_bulge_filter4  + penalty_amp_disk_filter4 + penalty_amp_bulge_filter5 + penalty_amp_disk_filter5 + penalty_amp_bulge_filter6 + penalty_amp_disk_filter6 + penalty_amp_bulge_filter7 + penalty_amp_disk_filter7 ;// + penalty_xcentre + penalty_ycentre

	}


        /* ------------------------------------*/
      /* Recompute chi2 with mask/isophotes */
      /* ------------------------------------*/

       lval = 0.0;
       sumxx=0.0;
       sumyy=0.0;
       sumxy=0.0;
       sumzy=0.0;
       sumxz=0.0;
       sumzz=0.0;
       denom=0.0;
       nfreedom_mask=0;

       for (j=0;j<npixels_acs;j++)
	 {

	   mask_filter1[j]=0.0;
	   mask_filter2[j]=0.0;
	   mask_filter3[j]=0.0;
	   mask_filter4[j] = 0.0;

	   if (itype_best_global==0) {
	   if (model_conv_bulge_filter1[j]>(model_conv_bulge_filter1[pix_bestmodel_acs]*sb_limit_bulge)) {mask_filter1[j]=1.0;}
	   if (model_conv_bulge_filter2[j]>(model_conv_bulge_filter2[pix_bestmodel_acs]*sb_limit_bulge)) {mask_filter2[j]=1.0;}
	   if (model_conv_bulge_filter3[j]>(model_conv_bulge_filter3[pix_bestmodel_acs]*sb_limit_bulge)) {mask_filter3[j]=1.0;}
	   if (model_conv_bulge_filter4[j]>(model_conv_bulge_filter4[pix_bestmodel_acs]*sb_limit_bulge)) {mask_filter4[j]=1.0;}
	   }

	   if (itype_best_global>0) {
	   if (model_conv_disk_filter1[j]>(model_conv_disk_filter1[pix_bestmodel_acs]*sb_limit_disk)) {mask_filter1[j]=1.0;}
	   if (model_conv_disk_filter2[j]>(model_conv_disk_filter2[pix_bestmodel_acs]*sb_limit_disk)) {mask_filter2[j]=1.0;}
	   if (model_conv_disk_filter3[j]>(model_conv_disk_filter3[pix_bestmodel_acs]*sb_limit_disk)) {mask_filter3[j]=1.0;}
	   if (model_conv_disk_filter4[j]>(model_conv_disk_filter4[pix_bestmodel_acs]*sb_limit_disk)) {mask_filter4[j]=1.0;}
	   }

	 }


       for (j=0;j<npixels_wfc3;j++)
	 {

	   mask_filter5[j] = 0.0;
	   mask_filter6[j] = 0.0;
	   mask_filter7[j] = 0.0;

	   if (itype_best_global==0) {
	   if (model_conv_bulge_filter5[j]>(model_conv_bulge_filter5[pix_bestmodel_wfc3]*sb_limit_bulge)) {mask_filter5[j]=1.0;}
	   if (model_conv_bulge_filter6[j]>(model_conv_bulge_filter6[pix_bestmodel_wfc3]*sb_limit_bulge)) {mask_filter6[j]=1.0;}
	   if (model_conv_bulge_filter7[j]>(model_conv_bulge_filter7[pix_bestmodel_wfc3]*sb_limit_bulge)) {mask_filter7[j]=1.0;}


	   }

	   if (itype_best_global>0) {

	   if (model_conv_disk_filter5[j]>(model_conv_disk_filter5[pix_bestmodel_wfc3]*sb_limit_disk)) {mask_filter5[j]=1.0;}
	   if (model_conv_disk_filter6[j]>(model_conv_disk_filter6[pix_bestmodel_wfc3]*sb_limit_disk)) {mask_filter6[j]=1.0;}
	   if (model_conv_disk_filter7[j]>(model_conv_disk_filter7[pix_bestmodel_wfc3]*sb_limit_disk)) {mask_filter7[j]=1.0;}
	   }

	 }


       /* Compute likelihood for E-type separately - bulge only */

      if (itype_best_global==0)
	{

	  sumx=sumy=sumxx=sumyy=sumxy=0.0;
	  lval = 0.0;

	  for (j=0; j<npixels_acs; j++)
	    {
	       if (model_conv_bulge_filter1[j]>(model_conv_bulge_filter1[pix_bestmodel_acs]*sb_limit_bulge))
		{

	      sumx += model_bulge_filter1[j];
	      sumy += data_filter1[j];
	      sumxx += pow(model_bulge_filter1[j],2);
	      sumyy += pow(data_filter1[j],2);
	      sumxy += data_filter1[j]*model_bulge_filter1[j];
		}
	    }

	  denom = sumxx;

	  if (denom != 0.)
	    {
	      amp_etype_filter1 = sumxy/sumxx;
	      amp_bulge_filter1=amp_etype_filter1;

	      for (j=0;j<npixels_acs;j++)
		{
		   if (model_conv_bulge_filter1[j]>(model_conv_bulge_filter1[pix_bestmodel_acs]*sb_limit_bulge))
		{
		  lval += pow((data_filter1[j]- (amp_etype_filter1 * model_bulge_filter1[j])),2);
		  nfreedom_mask++;
		}

		}

	    }
	  else
	    {
	      lval += 9e99;
	    }
	  sumx=sumy=sumxx=sumyy=sumxy=0.0;

	  //printf("amp_bulge (E-TYPE) lval: %6.3e %6.3e\n",amp_etype_filter1, lval);


	  /* 2nd filter */

          for (j=0; j<npixels_acs; j++)
            {

	      if (model_conv_bulge_filter2[j]>(model_conv_bulge_filter2[pix_bestmodel_acs]*sb_limit_bulge))
		{
              sumx += model_bulge_filter2[j];
              sumy += data_filter2[j];
              sumxx += pow(model_bulge_filter2[j],2);
              sumyy += pow(data_filter2[j],2);
              sumxy += data_filter2[j]*model_bulge_filter2[j];
		}
            }

          denom = sumxx;

          if (denom != 0.)
            {
              amp_etype_filter2 = sumxy/sumxx;
              amp_bulge_filter2=amp_etype_filter2;
              for (j=0;j<npixels_acs;j++)
                {
		   if (model_conv_bulge_filter2[j]>(model_conv_bulge_filter2[pix_bestmodel_acs]*sb_limit_bulge))
		{
                  lval += pow((data_filter2[j]- (amp_etype_filter2 * model_bulge_filter2[j])),2);
		  nfreedom_mask++;

		}
                }

            }
          else
            {
              lval += 9e99;
            }
          sumx=sumy=sumxx=sumyy=sumxy=0.0;

	  /* 3rd filter */
	  for (j=0; j<npixels_acs; j++)
            {
	       if (model_conv_bulge_filter3[j]>(model_conv_bulge_filter3[pix_bestmodel_acs]*sb_limit_bulge))
		{
              sumx += model_bulge_filter3[j];
              sumy += data_filter3[j];
              sumxx += pow(model_bulge_filter3[j],2);
              sumyy += pow(data_filter3[j],2);
              sumxy += data_filter3[j]*model_bulge_filter3[j];
		}
            }

          denom = sumxx;

          if (denom != 0.)
            {
              amp_etype_filter3 = sumxy/sumxx;
              amp_bulge_filter3=amp_etype_filter3;
              for (j=0;j<npixels_acs;j++)
                {
		   if (model_conv_bulge_filter3[j]>(model_conv_bulge_filter3[pix_bestmodel_acs]*sb_limit_bulge))
		{
                  lval += pow((data_filter3[j]- (amp_etype_filter3 * model_bulge_filter3[j])),2);
		  nfreedom_mask++;
		}

                }

            }
          else
            {
              lval += 9e99;
            }
          sumx=sumy=sumxx=sumyy=sumxy=0.0;

	  /* 4th filter */
	  for (j=0; j<npixels_acs; j++)
            {

	       if (model_conv_bulge_filter4[j]>(model_conv_bulge_filter4[pix_bestmodel_acs]*sb_limit_bulge))
		{
              sumx += model_bulge_filter4[j];
              sumy += data_filter4[j];
              sumxx += pow(model_bulge_filter4[j],2);
              sumyy += pow(data_filter4[j],2);
              sumxy += data_filter4[j]*model_bulge_filter4[j];
		}
            }

          denom = sumxx;

          if (denom != 0.)
            {
              amp_etype_filter4 = sumxy/sumxx;
              amp_bulge_filter4=amp_etype_filter4;
              for (j=0;j<npixels_acs;j++)
                {
		   if (model_conv_bulge_filter4[j]>(model_conv_bulge_filter4[pix_bestmodel_acs]*sb_limit_bulge))
		{
                  lval += pow((data_filter4[j]- (amp_etype_filter4 * model_bulge_filter4[j])),2);
		  nfreedom_mask++;

		}

                }

            }
          else
            {
              lval += 9e99;
            }
          sumx=sumy=sumxx=sumyy=sumxy=0.0;


	  /* 5th filter */
	  for (j=0; j<npixels_wfc3; j++)
            {

	       if (model_conv_bulge_filter5[j]>(model_conv_bulge_filter5[pix_bestmodel_wfc3]*sb_limit_bulge))
		{
              sumx += model_bulge_filter5[j];
              sumy += data_filter5[j];
              sumxx += pow(model_bulge_filter5[j],2);
              sumyy += pow(data_filter5[j],2);
              sumxy += data_filter5[j]*model_bulge_filter5[j];
		}
            }

          denom = sumxx;

          if (denom != 0.)
            {
              amp_etype_filter5 = sumxy/sumxx;
              amp_bulge_filter5=amp_etype_filter5;
              for (j=0;j<npixels_wfc3;j++)
                {
		   if (model_conv_bulge_filter5[j]>(model_conv_bulge_filter5[pix_bestmodel_wfc3]*sb_limit_bulge))
		{
                  lval += pow((data_filter5[j]- (amp_etype_filter5 * model_bulge_filter5[j])),2);
		  nfreedom_mask++;

		}

                }

            }
          else
            {
              lval += 9e99;
            }
          sumx=sumy=sumxx=sumyy=sumxy=0.0;


	  /* 6th filter */
	  for (j=0; j<npixels_wfc3; j++)
            {

	       if (model_conv_bulge_filter6[j]>(model_conv_bulge_filter6[pix_bestmodel_wfc3]*sb_limit_bulge))
		{
              sumx += model_bulge_filter6[j];
              sumy += data_filter6[j];
              sumxx += pow(model_bulge_filter6[j],2);
              sumyy += pow(data_filter6[j],2);
              sumxy += data_filter6[j]*model_bulge_filter6[j];
		}
            }

          denom = sumxx;

          if (denom != 0.)
            {
              amp_etype_filter6 = sumxy/sumxx;
              amp_bulge_filter6=amp_etype_filter6;
              for (j=0;j<npixels_wfc3;j++)
                {
		   if (model_conv_bulge_filter6[j]>(model_conv_bulge_filter6[pix_bestmodel_wfc3]*sb_limit_bulge))
		{
                  lval += pow((data_filter6[j]- (amp_etype_filter6 * model_bulge_filter6[j])),2);
		  nfreedom_mask++;

		}

                }

            }
          else
            {
              lval += 9e99;
            }
          sumx=sumy=sumxx=sumyy=sumxy=0.0;


	  /* 7th filter */
	  for (j=0; j<npixels_wfc3; j++)
            {

	       if (model_conv_bulge_filter7[j]>(model_conv_bulge_filter7[pix_bestmodel_wfc3]*sb_limit_bulge))
		{
              sumx += model_bulge_filter7[j];
              sumy += data_filter7[j];
              sumxx += pow(model_bulge_filter7[j],2);
              sumyy += pow(data_filter7[j],2);
              sumxy += data_filter7[j]*model_bulge_filter7[j];
		}
            }

          denom = sumxx;

          if (denom != 0.)
            {
              amp_etype_filter7 = sumxy/sumxx;
              amp_bulge_filter7=amp_etype_filter7;
              for (j=0;j<npixels_wfc3;j++)
                {
		   if (model_conv_bulge_filter7[j]>(model_conv_bulge_filter7[pix_bestmodel_wfc3]*sb_limit_bulge))
		{
                  lval += pow((data_filter7[j]- (amp_etype_filter7 * model_bulge_filter7[j])),2);
		  nfreedom_mask++;

		}

                }

            }
          else
            {
              lval += 9e99;
            }
          sumx=sumy=sumxx=sumyy=sumxy=0.0;

	  if (amp_bulge_filter1<0.0) {penalty_amp_bulge_filter1=1e7;}
	  if (amp_disk_filter1<0.0) {penalty_amp_disk_filter1=1e7;}
	  if (amp_bulge_filter2<0.0) {penalty_amp_bulge_filter2=1e7;}
	  if (amp_disk_filter2<0.0) {penalty_amp_disk_filter2=1e7;}
	  if (amp_bulge_filter3<0.0) {penalty_amp_bulge_filter3=1e7;}
	  if (amp_disk_filter3<0.0) {penalty_amp_disk_filter3=1e7;}
	  if (amp_bulge_filter4<0.0) {penalty_amp_bulge_filter4=1e7;}
	  if (amp_disk_filter4<0.0) {penalty_amp_disk_filter4=1e7;}

	  if (amp_bulge_filter5<0.0) {penalty_amp_bulge_filter5=1e7;}
	  if (amp_disk_filter5<0.0) {penalty_amp_disk_filter5=1e7;}
	  if (amp_bulge_filter6<0.0) {penalty_amp_bulge_filter6=1e7;}
	  if (amp_disk_filter6<0.0) {penalty_amp_disk_filter6=1e7;}
	  if (amp_bulge_filter7<0.0) {penalty_amp_bulge_filter7=1e7;}
	  if (amp_disk_filter7<0.0) {penalty_amp_disk_filter7=1e7;}

	  lval = lval + penalty_e1e2quad + penalty_scaling + penalty_amp_bulge_filter1 + penalty_amp_disk_filter1 +  penalty_amp_bulge_filter2 + penalty_amp_disk_filter2 + penalty_amp_bulge_filter3 + penalty_amp_disk_filter3 + penalty_amp_bulge_filter4  + penalty_amp_disk_filter4 + penalty_amp_bulge_filter5 + penalty_amp_disk_filter5 + penalty_amp_bulge_filter6 + penalty_amp_disk_filter6 + penalty_amp_bulge_filter7 + penalty_amp_disk_filter7 ;



	}



      /* Continue for other types */
      if (itype_best_global>0)
	{

	  lval=sumxx=sumzz=sumyy=sumxy=sumzy=sumxz=0.0;

	  /* 1 filter */
	  for (j=0; j<npixels_acs; j++)
	    {

	      if (model_conv_disk_filter1[j]>(model_conv_disk_filter1[pix_bestmodel_acs]*sb_limit_disk))
		{
	      sumxx += pow(model_disk_filter1[j],2);
	      sumzz += pow(model_bulge_filter1[j],2);
	      sumyy += pow(data_filter1[j],2);
	      sumxy += data_filter1[j]*model_disk_filter1[j];
	      sumzy += data_filter1[j]*model_bulge_filter1[j];
	      sumxz += model_disk_filter1[j]*model_bulge_filter1[j];
		}
	    }


	  denom = sumxx*sumzz-sumxz*sumxz;


	  if (denom != 0.0)
	    {

	      amp_disk_filter1 = (sumzz*sumxy-sumxz*sumzy)/denom; //disk
	      amp_bulge_filter1 = (sumxx*sumzy-sumxz*sumxy)/denom;  // bulge

	      for (j=0;j<npixels_acs;j++)
		{
		  if (model_conv_disk_filter1[j]>(model_conv_disk_filter1[pix_bestmodel_acs]*sb_limit_disk))
		{
		  lval = lval + pow(data_filter1[j]- ((amp_disk_filter1 * model_disk_filter1[j]) + (amp_bulge_filter1 * model_bulge_filter1[j]) ) ,2);
		  nfreedom_mask++;
		}

		}


	    }
	  else
	    {
	      printf("sumxx: %lf\n",sumxx);
	      printf("sumzz: %lf\n",sumzz);
	      printf("sumxz: %lf\n",sumxz);
	      lval += 9e99 ;
	    }

	  sumxx=0.0;
	  sumyy=0.0;
	  sumxy=0.0;
	  sumzy=0.0;
	  sumxz=0.0;
	  sumzz=0.0;

	  /* 2nd filter */
          for (j=0; j<npixels_acs; j++)
            {
	      if (model_conv_disk_filter2[j]>(model_conv_disk_filter2[pix_bestmodel_acs]*sb_limit_disk))
		{
              sumxx += pow(model_disk_filter2[j],2);
              sumzz += pow(model_bulge_filter2[j],2);
              sumyy += pow(data_filter2[j],2);
              sumxy += data_filter2[j]*model_disk_filter2[j];
              sumzy += data_filter2[j]*model_bulge_filter2[j];
              sumxz += model_disk_filter2[j]*model_bulge_filter2[j];
		}
            }


          denom = sumxx*sumzz-sumxz*sumxz;

          if (denom != 0.0)
            {

              amp_disk_filter2 = (sumzz*sumxy-sumxz*sumzy)/denom; //disk
	      amp_bulge_filter2 = (sumxx*sumzy-sumxz*sumxy)/denom;  // bulge
              for (j=0;j<npixels_acs;j++)
                {
		  if (model_conv_disk_filter2[j]>(model_conv_disk_filter2[pix_bestmodel_acs]*sb_limit_disk))
		{
                  lval = lval + pow(data_filter2[j]- ((amp_disk_filter2 * model_disk_filter2[j]) + (amp_bulge_filter2 * model_bulge_filter2[j]) ) ,2);
		  nfreedom_mask++;
		}

                }


            }
          else
            {
              printf("sumxx: %lf\n",sumxx);
              printf("sumzz: %lf\n",sumzz);
              printf("sumxz: %lf\n",sumxz);
              lval += 9e99 ;
            }

          sumxx=0.0;
          sumyy=0.0;
          sumxy=0.0;
          sumzy=0.0;
          sumxz=0.0;
          sumzz=0.0;

	  /* 3rd filter */
	   for (j=0; j<npixels_acs; j++)
            {

	      if (model_conv_disk_filter3[j]>(model_conv_disk_filter3[pix_bestmodel_acs]*sb_limit_disk))
		{
              sumxx += pow(model_disk_filter3[j],2);
              sumzz += pow(model_bulge_filter3[j],2);
              sumyy += pow(data_filter3[j],2);
              sumxy += data_filter3[j]*model_disk_filter3[j];
              sumzy += data_filter3[j]*model_bulge_filter3[j];
              sumxz += model_disk_filter3[j]*model_bulge_filter3[j];
		}
            }


          denom = sumxx*sumzz-sumxz*sumxz;

          if (denom != 0.0)
            {

              amp_disk_filter3 = (sumzz*sumxy-sumxz*sumzy)/denom; //disk
	      amp_bulge_filter3 = (sumxx*sumzy-sumxz*sumxy)/denom;  // bulge
              for (j=0;j<npixels_acs;j++)
                {
		  if (model_conv_disk_filter3[j]>(model_conv_disk_filter3[pix_bestmodel_acs]*sb_limit_disk))
		{
                  lval = lval + pow(data_filter3[j]- ((amp_disk_filter3 * model_disk_filter3[j]) + (amp_bulge_filter3 * model_bulge_filter3[j]) ) ,2);
		  nfreedom_mask++;
		}

                }


            }
          else
            {
              printf("sumxx: %lf\n",sumxx);
              printf("sumzz: %lf\n",sumzz);
              printf("sumxz: %lf\n",sumxz);
              lval += 9e99 ;
            }

          sumxx=0.0;
          sumyy=0.0;
          sumxy=0.0;
          sumzy=0.0;
          sumxz=0.0;
          sumzz=0.0;

	  /* 4th filter */

	   for (j=0; j<npixels_acs; j++)
            {

	      if (model_conv_disk_filter4[j]>(model_conv_disk_filter4[pix_bestmodel_acs]*sb_limit_disk))
		{
              sumxx += pow(model_disk_filter4[j],2);
              sumzz += pow(model_bulge_filter4[j],2);
              sumyy += pow(data_filter4[j],2);
              sumxy += data_filter4[j]*model_disk_filter4[j];
              sumzy += data_filter4[j]*model_bulge_filter4[j];
              sumxz += model_disk_filter4[j]*model_bulge_filter4[j];
		}
            }


          denom = sumxx*sumzz-sumxz*sumxz;

          if (denom != 0.0)
            {

              amp_disk_filter4 = (sumzz*sumxy-sumxz*sumzy)/denom; //disk
	      amp_bulge_filter4 = (sumxx*sumzy-sumxz*sumxy)/denom;  // bulge
              for (j=0;j<npixels_acs;j++)
                {
		  if (model_conv_disk_filter4[j]>(model_conv_disk_filter4[pix_bestmodel_acs]*sb_limit_disk))
		{
                  lval = lval + pow(data_filter4[j]- ((amp_disk_filter4 * model_disk_filter4[j]) + (amp_bulge_filter4 * model_bulge_filter4[j]) ) ,2);
		  nfreedom_mask++;
		}

                }


            }
          else
            {
              printf("sumxx: %lf\n",sumxx);
              printf("sumzz: %lf\n",sumzz);
              printf("sumxz: %lf\n",sumxz);
              lval += 9e99 ;
            }

          sumxx=0.0;
          sumyy=0.0;
          sumxy=0.0;
          sumzy=0.0;
          sumxz=0.0;
          sumzz=0.0;

	  /* 5th filter */

	   for (j=0; j<npixels_wfc3; j++)
            {

	      if (model_conv_disk_filter5[j]>(model_conv_disk_filter5[pix_bestmodel_wfc3]*sb_limit_disk))
		{
              sumxx += pow(model_disk_filter5[j],2);
              sumzz += pow(model_bulge_filter5[j],2);
              sumyy += pow(data_filter5[j],2);
              sumxy += data_filter5[j]*model_disk_filter5[j];
              sumzy += data_filter5[j]*model_bulge_filter5[j];
              sumxz += model_disk_filter5[j]*model_bulge_filter5[j];
		}
            }


          denom = sumxx*sumzz-sumxz*sumxz;

          if (denom != 0.0)
            {

              amp_disk_filter5 = (sumzz*sumxy-sumxz*sumzy)/denom; //disk
	      amp_bulge_filter5 = (sumxx*sumzy-sumxz*sumxy)/denom;  // bulge
              for (j=0;j<npixels_wfc3;j++)
                {
		  if (model_conv_disk_filter5[j]>(model_conv_disk_filter5[pix_bestmodel_wfc3]*sb_limit_disk))
		{
                  lval = lval + pow(data_filter5[j]- ((amp_disk_filter5 * model_disk_filter5[j]) + (amp_bulge_filter5 * model_bulge_filter5[j]) ) ,2);
		  nfreedom_mask++;
		}

                }


            }
          else
            {
              printf("sumxx: %lf\n",sumxx);
              printf("sumzz: %lf\n",sumzz);
              printf("sumxz: %lf\n",sumxz);
              lval += 9e99 ;
            }

          sumxx=0.0;
          sumyy=0.0;
          sumxy=0.0;
          sumzy=0.0;
          sumxz=0.0;
          sumzz=0.0;


	   /* 6th filter */

	   for (j=0; j<npixels_wfc3; j++)
            {

	      if (model_conv_disk_filter6[j]>(model_conv_disk_filter6[pix_bestmodel_wfc3]*sb_limit_disk))
		{
              sumxx += pow(model_disk_filter6[j],2);
              sumzz += pow(model_bulge_filter6[j],2);
              sumyy += pow(data_filter6[j],2);
              sumxy += data_filter6[j]*model_disk_filter6[j];
              sumzy += data_filter6[j]*model_bulge_filter6[j];
              sumxz += model_disk_filter6[j]*model_bulge_filter6[j];
		}
            }


          denom = sumxx*sumzz-sumxz*sumxz;

          if (denom != 0.0)
            {

              amp_disk_filter6 = (sumzz*sumxy-sumxz*sumzy)/denom; //disk
	      amp_bulge_filter6 = (sumxx*sumzy-sumxz*sumxy)/denom;  // bulge
              for (j=0;j<npixels_wfc3;j++)
                {
		  if (model_conv_disk_filter6[j]>(model_conv_disk_filter6[pix_bestmodel_wfc3]*sb_limit_disk))
		{
                  lval = lval + pow(data_filter6[j]- ((amp_disk_filter6 * model_disk_filter6[j]) + (amp_bulge_filter6 * model_bulge_filter6[j]) ) ,2);
		  nfreedom_mask++;
		}

                }


            }
          else
            {
              printf("sumxx: %lf\n",sumxx);
              printf("sumzz: %lf\n",sumzz);
              printf("sumxz: %lf\n",sumxz);
              lval += 9e99 ;
            }

          sumxx=0.0;
          sumyy=0.0;
          sumxy=0.0;
          sumzy=0.0;
          sumxz=0.0;
          sumzz=0.0;


	   /* 7th filter */

	   for (j=0; j<npixels_wfc3; j++)
            {

	      if (model_conv_disk_filter7[j]>(model_conv_disk_filter7[pix_bestmodel_wfc3]*sb_limit_disk))
		{
              sumxx += pow(model_disk_filter7[j],2);
              sumzz += pow(model_bulge_filter7[j],2);
              sumyy += pow(data_filter7[j],2);
              sumxy += data_filter7[j]*model_disk_filter7[j];
              sumzy += data_filter7[j]*model_bulge_filter7[j];
              sumxz += model_disk_filter7[j]*model_bulge_filter7[j];
		}
            }


          denom = sumxx*sumzz-sumxz*sumxz;

          if (denom != 0.0)
            {

              amp_disk_filter7 = (sumzz*sumxy-sumxz*sumzy)/denom; //disk
	      amp_bulge_filter7 = (sumxx*sumzy-sumxz*sumxy)/denom;  // bulge
              for (j=0;j<npixels_wfc3;j++)
                {
		  if (model_conv_disk_filter7[j]>(model_conv_disk_filter7[pix_bestmodel_wfc3]*sb_limit_disk))
		{
                  lval = lval + pow(data_filter7[j]- ((amp_disk_filter7 * model_disk_filter7[j]) + (amp_bulge_filter7 * model_bulge_filter7[j]) ) ,2);
		  nfreedom_mask++;
		}

                }


            }
          else
            {
              printf("sumxx: %lf\n",sumxx);
              printf("sumzz: %lf\n",sumzz);
              printf("sumxz: %lf\n",sumxz);
              lval += 9e99 ;
            }


          sumxx=0.0;
          sumyy=0.0;
          sumxy=0.0;
          sumzy=0.0;
          sumxz=0.0;
          sumzz=0.0;


	  if (amp_bulge_filter1<0.0) {penalty_amp_bulge_filter1=1e7;}
	  if (amp_disk_filter1<0.0) {penalty_amp_disk_filter1=1e7;}
	  if (amp_bulge_filter2<0.0) {penalty_amp_bulge_filter2=1e7;}
	  if (amp_disk_filter2<0.0) {penalty_amp_disk_filter2=1e7;}
	  if (amp_bulge_filter3<0.0) {penalty_amp_bulge_filter3=1e7;}
	  if (amp_disk_filter3<0.0) {penalty_amp_disk_filter3=1e7;}
	  if (amp_bulge_filter4<0.0) {penalty_amp_bulge_filter4=1e7;}
	  if (amp_disk_filter4<0.0) {penalty_amp_disk_filter4=1e7;}

	  if (amp_bulge_filter5<0.0) {penalty_amp_bulge_filter5=1e7;}
	  if (amp_disk_filter5<0.0) {penalty_amp_disk_filter5=1e7;}
	  if (amp_bulge_filter6<0.0) {penalty_amp_bulge_filter6=1e7;}
	  if (amp_disk_filter6<0.0) {penalty_amp_disk_filter6=1e7;}
	  if (amp_bulge_filter7<0.0) {penalty_amp_bulge_filter7=1e7;}
	  if (amp_disk_filter7<0.0) {penalty_amp_disk_filter7=1e7;}

	  lval = lval + penalty_e1e2quad + penalty_scaling + penalty_amp_bulge_filter1 + penalty_amp_disk_filter1 +  penalty_amp_bulge_filter2 + penalty_amp_disk_filter2 + penalty_amp_bulge_filter3 + penalty_amp_disk_filter3 + penalty_amp_bulge_filter4  + penalty_amp_disk_filter4 + penalty_amp_bulge_filter5 + penalty_amp_disk_filter5 + penalty_amp_bulge_filter6 + penalty_amp_disk_filter6 + penalty_amp_bulge_filter7 + penalty_amp_disk_filter7 ;

	}


      chi2_masked = lval;

      nfreedom_mask = nfreedom_mask - Nparams;


       /* -------------------------------------------------------------------- */
       /* Measure data and model colour gradients from these single band fits */
       /* -------------------------------------------------------------------- */

      totalflux_unconv_bulge_filter1=0.0;
      totalflux_unconv_disk_filter1=0.0;
      totalflux_unconv_bulge_filter2=0.0;
      totalflux_unconv_disk_filter2=0.0;
      totalflux_unconv_bulge_filter3=0.0;
      totalflux_unconv_disk_filter3=0.0;
      totalflux_unconv_bulge_filter4=0.0;
      totalflux_unconv_disk_filter4=0.0;
      totalflux_unconv_bulge_filter5=0.0;
      totalflux_unconv_disk_filter5=0.0;
      totalflux_unconv_bulge_filter6=0.0;
      totalflux_unconv_disk_filter6=0.0;
      totalflux_unconv_bulge_filter7=0.0;
      totalflux_unconv_disk_filter7=0.0;


      totalflux_conv_bulge_filter1=0.0;
      totalflux_conv_disk_filter1=0.0;
      totalflux_conv_bulge_filter2=0.0;
      totalflux_conv_disk_filter2=0.0;
      totalflux_conv_bulge_filter3=0.0;
      totalflux_conv_disk_filter3=0.0;
      totalflux_conv_bulge_filter4=0.0;
      totalflux_conv_disk_filter4=0.0;
      totalflux_conv_bulge_filter5=0.0;
      totalflux_conv_disk_filter5=0.0;
      totalflux_conv_bulge_filter6=0.0;
      totalflux_conv_disk_filter6=0.0;
      totalflux_conv_bulge_filter7=0.0;
      totalflux_conv_disk_filter7=0.0;


      if (itype_best_global>0)
	{

	  for (pix=0;pix<npixels_acs;pix++)
	    {

	      model_unconv_disk_filter1[pix] = amp_disk_filter1 * model_unconv_disk_filter1[pix];
	      model_unconv_bulge_filter1[pix] = amp_bulge_filter1 * model_unconv_bulge_filter1[pix];
	      model_conv_disk_filter1[pix] = amp_disk_filter1 * model_conv_disk_filter1[pix];
	      model_conv_bulge_filter1[pix] = amp_bulge_filter1 * model_conv_bulge_filter1[pix];
	      model_unconv_summed_filter1[pix] =  (model_unconv_disk_filter1[pix]) + (model_unconv_bulge_filter1[pix]);
	      model_conv_summed_filter1[pix] =    (model_conv_disk_filter1[pix]) + (model_conv_bulge_filter1[pix]);

	      totalflux_unconv_bulge_filter1 = totalflux_unconv_bulge_filter1 + (model_unconv_bulge_filter1[pix]);
	      totalflux_unconv_disk_filter1 = totalflux_unconv_disk_filter1 + (model_unconv_disk_filter1[pix]);
	      totalflux_conv_bulge_filter1 = totalflux_conv_bulge_filter1 + (model_conv_bulge_filter1[pix]);
              totalflux_conv_disk_filter1 = totalflux_conv_disk_filter1 + (model_conv_disk_filter1[pix]);

	      model_unconv_disk_filter2[pix] = amp_disk_filter2 * model_unconv_disk_filter2[pix];
	      model_unconv_bulge_filter2[pix] = amp_bulge_filter2 * model_unconv_bulge_filter2[pix];
	      model_conv_disk_filter2[pix] = amp_disk_filter2 * model_conv_disk_filter2[pix];
	      model_conv_bulge_filter2[pix] = amp_bulge_filter2 * model_conv_bulge_filter2[pix];
	      model_unconv_summed_filter2[pix] =  (model_unconv_disk_filter2[pix]) + (model_unconv_bulge_filter2[pix]);
	      model_conv_summed_filter2[pix] =    (model_conv_disk_filter2[pix]) + (model_conv_bulge_filter2[pix]);

	      totalflux_unconv_bulge_filter2 = totalflux_unconv_bulge_filter2 + (model_unconv_bulge_filter2[pix]);
	      totalflux_unconv_disk_filter2 = totalflux_unconv_disk_filter2 + (model_unconv_disk_filter2[pix]);
	      totalflux_conv_bulge_filter2 = totalflux_conv_bulge_filter2 + (model_conv_bulge_filter2[pix]);
              totalflux_conv_disk_filter2 = totalflux_conv_disk_filter2 + (model_conv_disk_filter2[pix]);

	      model_unconv_disk_filter3[pix] = amp_disk_filter3 * model_unconv_disk_filter3[pix];
	      model_unconv_bulge_filter3[pix] = amp_bulge_filter3 * model_unconv_bulge_filter3[pix];
	      model_conv_disk_filter3[pix] = amp_disk_filter3 * model_conv_disk_filter3[pix];
	      model_conv_bulge_filter3[pix] = amp_bulge_filter3 * model_conv_bulge_filter3[pix];
	      model_unconv_summed_filter3[pix] =  (model_unconv_disk_filter3[pix]) + (model_unconv_bulge_filter3[pix]);
	      model_conv_summed_filter3[pix] =    (model_conv_disk_filter3[pix]) + (model_conv_bulge_filter3[pix]);

	      totalflux_unconv_bulge_filter3 = totalflux_unconv_bulge_filter3 + (model_unconv_bulge_filter3[pix]);
	      totalflux_unconv_disk_filter3 = totalflux_unconv_disk_filter3 + (model_unconv_disk_filter3[pix]);
	      totalflux_conv_bulge_filter3 = totalflux_conv_bulge_filter3 + (model_conv_bulge_filter3[pix]);
              totalflux_conv_disk_filter3 = totalflux_conv_disk_filter3 + (model_conv_disk_filter3[pix]);

	      model_unconv_disk_filter4[pix] = amp_disk_filter4 * model_unconv_disk_filter4[pix];
	      model_unconv_bulge_filter4[pix] = amp_bulge_filter4 * model_unconv_bulge_filter4[pix];
	      model_conv_disk_filter4[pix] = amp_disk_filter4 * model_conv_disk_filter4[pix];
	      model_conv_bulge_filter4[pix] = amp_bulge_filter4 * model_conv_bulge_filter4[pix];
	      model_unconv_summed_filter4[pix] =  (model_unconv_disk_filter4[pix]) + (model_unconv_bulge_filter4[pix]);
	      model_conv_summed_filter4[pix] =    (model_conv_disk_filter4[pix]) + (model_conv_bulge_filter4[pix]);

	      totalflux_unconv_bulge_filter4 = totalflux_unconv_bulge_filter4 + (model_unconv_bulge_filter4[pix]);
	      totalflux_unconv_disk_filter4 = totalflux_unconv_disk_filter4 + (model_unconv_disk_filter4[pix]);
	      totalflux_conv_bulge_filter4 = totalflux_conv_bulge_filter4 + (model_conv_bulge_filter4[pix]);
              totalflux_conv_disk_filter4 = totalflux_conv_disk_filter4 + (model_conv_disk_filter4[pix]);

	      model_conv_residual_filter1[pix] = (simdata_filter1[pix] - model_conv_summed_filter1[pix])/simrms_filter1[pix] ;
      	      model_conv_residual_filter2[pix] = (simdata_filter2[pix] - model_conv_summed_filter2[pix])/simrms_filter2[pix] ;
	      model_conv_residual_filter3[pix] = (simdata_filter3[pix] - model_conv_summed_filter3[pix])/simrms_filter3[pix] ;
      	      model_conv_residual_filter4[pix] = (simdata_filter4[pix] - model_conv_summed_filter4[pix])/simrms_filter4[pix] ;


	    }


	   for (pix=0;pix<npixels_wfc3;pix++)
	    {

	      model_unconv_disk_filter5[pix] = amp_disk_filter5 * model_unconv_disk_filter5[pix];
	      model_unconv_bulge_filter5[pix] = amp_bulge_filter5 * model_unconv_bulge_filter5[pix];
	      model_conv_disk_filter5[pix] = amp_disk_filter5 * model_conv_disk_filter5[pix];
	      model_conv_bulge_filter5[pix] = amp_bulge_filter5 * model_conv_bulge_filter5[pix];
	      model_unconv_summed_filter5[pix] =  (model_unconv_disk_filter5[pix]) + (model_unconv_bulge_filter5[pix]);
	      model_conv_summed_filter5[pix] =    (model_conv_disk_filter5[pix]) + (model_conv_bulge_filter5[pix]);

	      totalflux_unconv_bulge_filter5 = totalflux_unconv_bulge_filter5 + (model_unconv_bulge_filter5[pix]);
	      totalflux_unconv_disk_filter5 = totalflux_unconv_disk_filter5 + (model_unconv_disk_filter5[pix]);
	      totalflux_conv_bulge_filter5 = totalflux_conv_bulge_filter5 + (model_conv_bulge_filter5[pix]);
              totalflux_conv_disk_filter5 = totalflux_conv_disk_filter5 + (model_conv_disk_filter5[pix]);

	      model_unconv_disk_filter6[pix] = amp_disk_filter6 * model_unconv_disk_filter6[pix];
	      model_unconv_bulge_filter6[pix] = amp_bulge_filter6 * model_unconv_bulge_filter6[pix];
	      model_conv_disk_filter6[pix] = amp_disk_filter6 * model_conv_disk_filter6[pix];
	      model_conv_bulge_filter6[pix] = amp_bulge_filter6 * model_conv_bulge_filter6[pix];
	      model_unconv_summed_filter6[pix] =  (model_unconv_disk_filter6[pix]) + (model_unconv_bulge_filter6[pix]);
	      model_conv_summed_filter6[pix] =    (model_conv_disk_filter6[pix]) + (model_conv_bulge_filter6[pix]);

	      totalflux_unconv_bulge_filter6 = totalflux_unconv_bulge_filter6 + (model_unconv_bulge_filter6[pix]);
	      totalflux_unconv_disk_filter6 = totalflux_unconv_disk_filter6 + (model_unconv_disk_filter6[pix]);
	      totalflux_conv_bulge_filter6 = totalflux_conv_bulge_filter6 + (model_conv_bulge_filter6[pix]);
              totalflux_conv_disk_filter6 = totalflux_conv_disk_filter6 + (model_conv_disk_filter6[pix]);

	      model_unconv_disk_filter7[pix] = amp_disk_filter7 * model_unconv_disk_filter7[pix];
	      model_unconv_bulge_filter7[pix] = amp_bulge_filter7 * model_unconv_bulge_filter7[pix];
	      model_conv_disk_filter7[pix] = amp_disk_filter7 * model_conv_disk_filter7[pix];
	      model_conv_bulge_filter7[pix] = amp_bulge_filter7 * model_conv_bulge_filter7[pix];
	      model_unconv_summed_filter7[pix] =  (model_unconv_disk_filter7[pix]) + (model_unconv_bulge_filter7[pix]);
	      model_conv_summed_filter7[pix] =    (model_conv_disk_filter7[pix]) + (model_conv_bulge_filter7[pix]);

	      totalflux_unconv_bulge_filter7 = totalflux_unconv_bulge_filter7 + (model_unconv_bulge_filter7[pix]);
	      totalflux_unconv_disk_filter7 = totalflux_unconv_disk_filter7 + (model_unconv_disk_filter7[pix]);
	      totalflux_conv_bulge_filter7 = totalflux_conv_bulge_filter7 + (model_conv_bulge_filter7[pix]);
              totalflux_conv_disk_filter7 = totalflux_conv_disk_filter7 + (model_conv_disk_filter7[pix]);

	      model_conv_residual_filter5[pix] = (simdata_filter5[pix] - model_conv_summed_filter5[pix])/simrms_filter5[pix] ;
      	      model_conv_residual_filter6[pix] = (simdata_filter6[pix] - model_conv_summed_filter6[pix])/simrms_filter6[pix] ;
	      model_conv_residual_filter7[pix] = (simdata_filter7[pix] - model_conv_summed_filter7[pix])/simrms_filter7[pix] ;

	    }

	}


      if (itype_best_global==0) /* E type - only bulge so set this to zero */
	{

	  for (pix=0;pix<npixels_acs;pix++)
	    {

	      model_unconv_bulge_filter1[pix] = amp_etype_filter1 * model_unconv_bulge_filter1[pix];
	      model_conv_bulge_filter1[pix] = amp_etype_filter1 * model_conv_bulge_filter1[pix];
	      model_unconv_summed_filter1[pix] = model_unconv_bulge_filter1[pix];
	      model_conv_summed_filter1[pix] = model_conv_bulge_filter1[pix];

	      totalflux_conv_bulge_filter1 = totalflux_conv_bulge_filter1 + (model_conv_bulge_filter1[pix]);
	      totalflux_conv_disk_filter1 = 0.0;
	      model_conv_residual_filter1[pix] = (simdata_filter1[pix] - model_conv_summed_filter1[pix])/simrms_filter1[pix] ;

	      model_unconv_bulge_filter2[pix] = amp_etype_filter2 * model_unconv_bulge_filter2[pix];
	      model_conv_bulge_filter2[pix] = amp_etype_filter2 * model_conv_bulge_filter2[pix];
	      model_unconv_summed_filter2[pix] = model_unconv_bulge_filter2[pix];
	      model_conv_summed_filter2[pix] = model_conv_bulge_filter2[pix];

	      totalflux_conv_bulge_filter2 = totalflux_conv_bulge_filter2 + (model_conv_bulge_filter2[pix]);
	      totalflux_conv_disk_filter2 = 0.0;
	      model_conv_residual_filter2[pix] = (simdata_filter2[pix] - model_conv_summed_filter2[pix])/simrms_filter2[pix] ;

	      model_unconv_bulge_filter3[pix] = amp_etype_filter3 * model_unconv_bulge_filter3[pix];
	      model_conv_bulge_filter3[pix] = amp_etype_filter3 * model_conv_bulge_filter3[pix];
	      model_unconv_summed_filter3[pix] = model_unconv_bulge_filter3[pix];
	      model_conv_summed_filter3[pix] = model_conv_bulge_filter3[pix];

	      totalflux_conv_bulge_filter3 = totalflux_conv_bulge_filter3 + (model_conv_bulge_filter3[pix]);
	      totalflux_conv_disk_filter3 = 0.0;
	      model_conv_residual_filter3[pix] = (simdata_filter3[pix] - model_conv_summed_filter3[pix])/simrms_filter3[pix] ;

	      model_unconv_bulge_filter4[pix] = amp_etype_filter4 * model_unconv_bulge_filter4[pix];
	      model_conv_bulge_filter4[pix] = amp_etype_filter4 * model_conv_bulge_filter4[pix];
	      model_unconv_summed_filter4[pix] = model_unconv_bulge_filter4[pix];
	      model_conv_summed_filter4[pix] = model_conv_bulge_filter4[pix];

	      totalflux_conv_bulge_filter4 = totalflux_conv_bulge_filter4 + (model_conv_bulge_filter4[pix]);
	      totalflux_conv_disk_filter4 = 0.0;
	      model_conv_residual_filter4[pix] = (simdata_filter4[pix] - model_conv_summed_filter4[pix])/simrms_filter4[pix] ;

	    }


	   for (pix=0;pix<npixels_wfc3;pix++)
	    {

	      model_unconv_bulge_filter5[pix] = amp_etype_filter5 * model_unconv_bulge_filter5[pix];
	      model_conv_bulge_filter5[pix] = amp_etype_filter5 * model_conv_bulge_filter5[pix];
	      model_unconv_summed_filter5[pix] = model_unconv_bulge_filter5[pix];
	      model_conv_summed_filter5[pix] = model_conv_bulge_filter5[pix];

	      totalflux_conv_bulge_filter5 = totalflux_conv_bulge_filter5 + (model_conv_bulge_filter5[pix]);
	      totalflux_conv_disk_filter5 = 0.0;
	      model_conv_residual_filter5[pix] = (simdata_filter5[pix] - model_conv_summed_filter5[pix])/simrms_filter5[pix] ;

	      model_unconv_bulge_filter6[pix] = amp_etype_filter6 * model_unconv_bulge_filter6[pix];
	      model_conv_bulge_filter6[pix] = amp_etype_filter6 * model_conv_bulge_filter6[pix];
	      model_unconv_summed_filter6[pix] = model_unconv_bulge_filter6[pix];
	      model_conv_summed_filter6[pix] = model_conv_bulge_filter6[pix];

	      totalflux_conv_bulge_filter6 = totalflux_conv_bulge_filter6 + (model_conv_bulge_filter6[pix]);
	      totalflux_conv_disk_filter6 = 0.0;
	      model_conv_residual_filter6[pix] = (simdata_filter6[pix] - model_conv_summed_filter6[pix])/simrms_filter6[pix] ;

	      model_unconv_bulge_filter7[pix] = amp_etype_filter7 * model_unconv_bulge_filter7[pix];
	      model_conv_bulge_filter7[pix] = amp_etype_filter7 * model_conv_bulge_filter7[pix];
	      model_unconv_summed_filter7[pix] = model_unconv_bulge_filter7[pix];
	      model_conv_summed_filter7[pix] = model_conv_bulge_filter7[pix];

	      totalflux_conv_bulge_filter7 = totalflux_conv_bulge_filter7 + (model_conv_bulge_filter7[pix]);
	      totalflux_conv_disk_filter7 = 0.0;
	      model_conv_residual_filter7[pix] = (simdata_filter7[pix] - model_conv_summed_filter7[pix])/simrms_filter7[pix] ;

	    }


	}





	  /* ----------------------------------- */
	  /* Measure data and model CG           */
	  /* ---------------------------------- */

	  /* Sum over weights for DATA = S_i(no noise)/sigma_i */

	  tt=0; ss=0;
	  sum_weights_inside_filter1 = 0.0;
          sum_weights_outside_filter1 = 0.0;
	  sum_weights_galaxy_filter1 = 0.0;
	  sum_weights_inside_filter2 = 0.0;
          sum_weights_outside_filter2 = 0.0;
	  sum_weights_galaxy_filter2 = 0.0;
	  sum_weights_inside_filter3 = 0.0;
          sum_weights_outside_filter3 = 0.0;
	  sum_weights_galaxy_filter3 = 0.0;
	  sum_weights_inside_filter4 = 0.0;
          sum_weights_outside_filter4 = 0.0;
	  sum_weights_galaxy_filter4 = 0.0;
	  sum_weights_inside_filter5 = 0.0;
          sum_weights_outside_filter5 = 0.0;
	  sum_weights_galaxy_filter5 = 0.0;
	  sum_weights_inside_filter6 = 0.0;
          sum_weights_outside_filter6 = 0.0;
	  sum_weights_galaxy_filter6 = 0.0;
	  sum_weights_inside_filter7 = 0.0;
          sum_weights_outside_filter7 = 0.0;
	  sum_weights_galaxy_filter7 = 0.0;

	  sumflux_data_inside_filter1 = 0.0;
	  sumflux_data_outside_filter1 = 0.0;
	  sumflux_model_inside_filter1 = 0.0;
	  sumflux_model_outside_filter1 = 0.0;
	  sumflux_model_galaxy_filter1 = 0.0;
	  sumflux_data_galaxy_filter1=0.0;

	  sumflux_data_inside_filter2 = 0.0;
	  sumflux_data_outside_filter2 = 0.0;
	  sumflux_model_inside_filter2 = 0.0;
	  sumflux_model_outside_filter2 = 0.0;
	  sumflux_model_galaxy_filter2 = 0.0;
	  sumflux_data_galaxy_filter2=0.0;

	  sumflux_data_inside_filter3 = 0.0;
	  sumflux_data_outside_filter3 = 0.0;
	  sumflux_model_inside_filter3 = 0.0;
	  sumflux_model_outside_filter3 = 0.0;
	  sumflux_model_galaxy_filter3 = 0.0;
	  sumflux_data_galaxy_filter3=0.0;

	  sumflux_data_inside_filter4 = 0.0;
	  sumflux_data_outside_filter4 = 0.0;
	  sumflux_model_inside_filter4 = 0.0;
	  sumflux_model_outside_filter4 = 0.0;
	  sumflux_model_galaxy_filter4 = 0.0;
	  sumflux_data_galaxy_filter4=0.0;

	  sumflux_data_inside_filter5 = 0.0;
	  sumflux_data_outside_filter5 = 0.0;
	  sumflux_model_inside_filter5 = 0.0;
	  sumflux_model_outside_filter5 = 0.0;
	  sumflux_model_galaxy_filter5 = 0.0;
	  sumflux_data_galaxy_filter5=0.0;

	  sumflux_data_inside_filter6 = 0.0;
	  sumflux_data_outside_filter6 = 0.0;
	  sumflux_model_inside_filter6 = 0.0;
	  sumflux_model_outside_filter6 = 0.0;
	  sumflux_model_galaxy_filter6 = 0.0;
	  sumflux_data_galaxy_filter6=0.0;

	  sumflux_data_inside_filter7 = 0.0;
	  sumflux_data_outside_filter7 = 0.0;
	  sumflux_model_inside_filter7 = 0.0;
	  sumflux_model_outside_filter7 = 0.0;
	  sumflux_model_galaxy_filter7 = 0.0;
	  sumflux_data_galaxy_filter7=0.0;

	  sumflux_bulge_filter1 = 0.0;
	  sumflux_disk_filter1 = 0.0;
	  sumflux_bulge_filter2 = 0.0;
	  sumflux_disk_filter2 = 0.0;
	  sumflux_bulge_filter3 = 0.0;
	  sumflux_disk_filter3 = 0.0;
	  sumflux_bulge_filter4 = 0.0;
	  sumflux_disk_filter4 = 0.0;
	  sumflux_bulge_filter5 = 0.0;
	  sumflux_disk_filter5 = 0.0;
	  sumflux_bulge_filter6 = 0.0;
	  sumflux_disk_filter6 = 0.0;
	  sumflux_bulge_filter7 = 0.0;
	  sumflux_disk_filter7 = 0.0;

	  sumrms_data_inside_filter1 = 0.0;
	  sumrms_data_outside_filter1 = 0.0;
	  sumrms_data_galaxy_filter1 = 0.0;
	  sumrms_data_inside_filter2 = 0.0;
	  sumrms_data_outside_filter2 = 0.0;
	  sumrms_data_galaxy_filter2 = 0.0;
	  sumrms_data_inside_filter3 = 0.0;
	  sumrms_data_outside_filter3 = 0.0;
	  sumrms_data_galaxy_filter3 = 0.0;
	  sumrms_data_inside_filter4 = 0.0;
	  sumrms_data_outside_filter4 = 0.0;
	  sumrms_data_galaxy_filter4 = 0.0;
	  sumrms_data_inside_filter5 = 0.0;
	  sumrms_data_outside_filter5 = 0.0;
	  sumrms_data_galaxy_filter5 = 0.0;
	  sumrms_data_inside_filter6 = 0.0;
	  sumrms_data_outside_filter6 = 0.0;
	  sumrms_data_galaxy_filter6 = 0.0;
	  sumrms_data_inside_filter7 = 0.0;
	  sumrms_data_outside_filter7 = 0.0;
	  sumrms_data_galaxy_filter7 = 0.0;

	  /* ACS */
	  for (tt=0;tt<ncols_tot_acs;tt++) {
            for (ss=0;ss<nrows_tot_acs;ss++){
	      pix = (ncols_tot_acs *ss) + tt;
	      radius_pix = sqrt( pow(tt-x0_bestmodel_acs,2) + pow(ss-y0_bestmodel_acs,2) ) * pixelsize_arcsec_acs ;  ///arcsec

	      weight_filter1[pix] = model_conv_summed_filter1[pix]/pow(simrms_filter1[pix],2);
	      weight_filter2[pix] = model_conv_summed_filter2[pix]/pow(simrms_filter2[pix],2);
	      weight_filter3[pix] = model_conv_summed_filter3[pix]/pow(simrms_filter3[pix],2);
	      weight_filter4[pix] = model_conv_summed_filter4[pix]/pow(simrms_filter4[pix],2);


	      if (radius_pix <= (0.3*R50))
                {
		  sum_weights_inside_filter1 += weight_filter1[pix] ;
		  sum_weights_inside_filter2 += weight_filter2[pix] ;
		  sum_weights_inside_filter3 += weight_filter3[pix] ;
		  sum_weights_inside_filter4 += weight_filter4[pix] ;

		}

              if ( (radius_pix > (0.3*R50)) && (radius_pix <= (R50)))
		{
                  sum_weights_outside_filter1 += weight_filter1[pix] ;
                  sum_weights_outside_filter2 += weight_filter2[pix] ;
		  sum_weights_outside_filter3 += weight_filter3[pix] ;
                  sum_weights_outside_filter4 += weight_filter4[pix] ;

		}


		  sum_weights_galaxy_filter1 += weight_filter1[pix] ;
		  sum_weights_galaxy_filter2 += weight_filter2[pix] ;
		  sum_weights_galaxy_filter3 += weight_filter3[pix] ;
		  sum_weights_galaxy_filter4 += weight_filter4[pix] ;


	    }}

	  /* WFC3 */
	  for (tt=0;tt<ncols_tot_wfc3;tt++) {
            for (ss=0;ss<nrows_tot_wfc3;ss++){
	      pix = (ncols_tot_wfc3 *ss) + tt;
	      radius_pix = sqrt( pow(tt-x0_bestmodel_wfc3,2) + pow(ss-y0_bestmodel_wfc3,2) ) * pixelsize_arcsec_wfc3 ;  ///arcsec

	      weight_filter5[pix] = model_conv_summed_filter5[pix]/pow(simrms_filter5[pix],2);
	      weight_filter6[pix] = model_conv_summed_filter6[pix]/pow(simrms_filter6[pix],2);
	      weight_filter7[pix] = model_conv_summed_filter7[pix]/pow(simrms_filter7[pix],2);


	      if (radius_pix <= (0.3*R50))
                {
		  sum_weights_inside_filter5 += weight_filter5[pix] ;
		  sum_weights_inside_filter6 += weight_filter6[pix] ;
		  sum_weights_inside_filter7 += weight_filter7[pix] ;

		}

              if ( (radius_pix > (0.3*R50)) && (radius_pix <= (R50)))
		{
                  sum_weights_outside_filter5 += weight_filter5[pix] ;
                  sum_weights_outside_filter6 += weight_filter6[pix] ;
		  sum_weights_outside_filter7 += weight_filter7[pix] ;
		}


		  sum_weights_galaxy_filter5 += weight_filter5[pix] ;
		  sum_weights_galaxy_filter6 += weight_filter6[pix] ;
		  sum_weights_galaxy_filter7 += weight_filter7[pix] ;

	    }}

	  /* ACS */
	   for (tt=0;tt<ncols_tot_acs;tt++) {
            for (ss=0;ss<nrows_tot_acs;ss++)
            {

              pix = (ncols_tot_acs *ss) + tt;
              radius_pix = sqrt( pow(tt-x0_bestmodel_acs,2) + pow(ss-y0_bestmodel_acs,2) ) * pixelsize_arcsec_acs ;  ///arcsec

              if (radius_pix <= (0.3*R50))
                {

                  sumflux_data_inside_filter1  += (weight_filter1[pix] * simdata_filter1[pix]);
                  sumflux_data_inside_filter2 += (weight_filter2[pix] * simdata_filter2[pix]);
		  sumflux_data_inside_filter3  += (weight_filter3[pix] * simdata_filter3[pix]);
                  sumflux_data_inside_filter4 += (weight_filter4[pix] * simdata_filter4[pix]);

		  sumflux_model_inside_filter1  += (model_conv_summed_filter1[pix]);
                  sumflux_model_inside_filter2 += (model_conv_summed_filter2[pix]);
		  sumflux_model_inside_filter3  += (model_conv_summed_filter3[pix]);
                  sumflux_model_inside_filter4 += (model_conv_summed_filter4[pix]);

                  sumrms_data_inside_filter1 += pow(simrms_filter1[pix],2)  ;
                  sumrms_data_inside_filter2 += pow(simrms_filter2[pix],2)  ;
		  sumrms_data_inside_filter3 += pow(simrms_filter3[pix],2)  ;
                  sumrms_data_inside_filter4 += pow(simrms_filter4[pix],2)  ;
                }

              if ( (radius_pix > (0.3*R50)) && (radius_pix <= R50))
                {

                  sumflux_data_outside_filter1 += (weight_filter1[pix] * simdata_filter1[pix]);
                  sumflux_data_outside_filter2 += (weight_filter2[pix] * simdata_filter2[pix]);
		  sumflux_data_outside_filter3 += (weight_filter3[pix] * simdata_filter3[pix]);
                  sumflux_data_outside_filter4 += (weight_filter4[pix] * simdata_filter4[pix]);

		  sumflux_model_outside_filter1  += (model_conv_summed_filter1[pix]);
                  sumflux_model_outside_filter2 += (model_conv_summed_filter2[pix]);
		  sumflux_model_outside_filter3  += (model_conv_summed_filter3[pix]);
                  sumflux_model_outside_filter4 += (model_conv_summed_filter4[pix]);

                  sumrms_data_outside_filter1 += pow(simrms_filter1[pix],2)  ;
                  sumrms_data_outside_filter2 += pow(simrms_filter2[pix],2)  ;
		  sumrms_data_outside_filter3 += pow(simrms_filter3[pix],2)  ;
                  sumrms_data_outside_filter4 += pow(simrms_filter4[pix],2)  ;


                }

	          sumflux_data_galaxy_filter1 += (weight_filter1[pix] * simdata_filter1[pix]);
                  sumflux_data_galaxy_filter2 += (weight_filter2[pix] * simdata_filter2[pix]);
		  sumflux_data_galaxy_filter3 += (weight_filter3[pix] * simdata_filter3[pix]);
                  sumflux_data_galaxy_filter4 += (weight_filter4[pix] * simdata_filter4[pix]);

		  sumflux_bulge_filter1 += (model_conv_bulge_filter1[pix]) ;
		  sumflux_bulge_filter2 += (model_conv_bulge_filter2[pix]) ;
		  sumflux_bulge_filter3 += (model_conv_bulge_filter3[pix]) ;
		  sumflux_bulge_filter4 += (model_conv_bulge_filter4[pix]) ;

		  sumflux_disk_filter1 += (model_conv_disk_filter1[pix]) ;
                  sumflux_disk_filter2 += (model_conv_disk_filter2[pix]) ;
		  sumflux_disk_filter3 += (model_conv_disk_filter3[pix]) ;
                  sumflux_disk_filter4 += (model_conv_disk_filter4[pix]) ;

		  sumflux_model_galaxy_filter1  += (model_unconv_summed_filter1[pix]);
                  sumflux_model_galaxy_filter2 += (model_unconv_summed_filter2[pix]);
		  sumflux_model_galaxy_filter3  += (model_unconv_summed_filter3[pix]);
                  sumflux_model_galaxy_filter4 += (model_unconv_summed_filter4[pix]);

		  sumrms_data_galaxy_filter1 += pow(simrms_filter1[pix],2)  ;
                  sumrms_data_galaxy_filter2 += pow(simrms_filter2[pix],2)  ;
		  sumrms_data_galaxy_filter3 += pow(simrms_filter3[pix],2)  ;
                  sumrms_data_galaxy_filter4 += pow(simrms_filter4[pix],2)  ;

	    }}

	      /* WFC3 */

	    for (tt=0;tt<ncols_tot_wfc3;tt++) {
            for (ss=0;ss<nrows_tot_wfc3;ss++)
            {

              pix = (ncols_tot_wfc3 *ss) + tt;
              radius_pix = sqrt( pow(tt-x0_bestmodel_wfc3,2) + pow(ss-y0_bestmodel_wfc3,2) ) * pixelsize_arcsec_wfc3 ;  ///arcsec

              if (radius_pix <= (0.3*R50))
                {

                  sumflux_data_inside_filter5  += (weight_filter5[pix] * simdata_filter5[pix]);
                  sumflux_data_inside_filter6 += (weight_filter6[pix] * simdata_filter6[pix]);
		  sumflux_data_inside_filter7  += (weight_filter7[pix] * simdata_filter7[pix]);

		  sumflux_model_inside_filter5  += (model_conv_summed_filter5[pix]);
                  sumflux_model_inside_filter6 += (model_conv_summed_filter6[pix]);
		  sumflux_model_inside_filter7  += (model_conv_summed_filter7[pix]);

                  sumrms_data_inside_filter5 += pow(simrms_filter5[pix],2)  ;
                  sumrms_data_inside_filter6 += pow(simrms_filter6[pix],2)  ;
		  sumrms_data_inside_filter7 += pow(simrms_filter7[pix],2)  ;
                }

              if ( (radius_pix > (0.3*R50)) && (radius_pix <= R50))
                {

                  sumflux_data_outside_filter5 += (weight_filter5[pix] * simdata_filter5[pix]);
                  sumflux_data_outside_filter6 += (weight_filter6[pix] * simdata_filter6[pix]);
		  sumflux_data_outside_filter7 += (weight_filter7[pix] * simdata_filter7[pix]);

		  sumflux_model_outside_filter5  += (model_conv_summed_filter5[pix]);
                  sumflux_model_outside_filter6 += (model_conv_summed_filter6[pix]);
		  sumflux_model_outside_filter7  += (model_conv_summed_filter7[pix]);

                  sumrms_data_outside_filter5 += pow(simrms_filter5[pix],2)  ;
                  sumrms_data_outside_filter6 += pow(simrms_filter6[pix],2)  ;
		  sumrms_data_outside_filter7 += pow(simrms_filter7[pix],2)  ;

                }



		  sumflux_data_galaxy_filter5 += (weight_filter5[pix] * simdata_filter5[pix]);
                  sumflux_data_galaxy_filter6 += (weight_filter6[pix] * simdata_filter6[pix]);
		  sumflux_data_galaxy_filter7 += (weight_filter7[pix] * simdata_filter7[pix]);

		  sumflux_bulge_filter5 += (model_conv_bulge_filter5[pix]) ;
		  sumflux_bulge_filter6 += (model_conv_bulge_filter6[pix]) ;
		  sumflux_bulge_filter7 += (model_conv_bulge_filter7[pix]) ;

		  sumflux_disk_filter5 += (model_conv_disk_filter5[pix]) ;
                  sumflux_disk_filter6 += (model_conv_disk_filter6[pix]) ;
		  sumflux_disk_filter7 += (model_conv_disk_filter7[pix]) ;

		  sumflux_model_galaxy_filter5  += (model_unconv_summed_filter5[pix]);
                  sumflux_model_galaxy_filter6 += (model_unconv_summed_filter6[pix]);
		  sumflux_model_galaxy_filter7  += (model_unconv_summed_filter7[pix]);

		  sumrms_data_galaxy_filter5 += pow(simrms_filter5[pix],2)  ;
                  sumrms_data_galaxy_filter6 += pow(simrms_filter6[pix],2)  ;
		  sumrms_data_galaxy_filter7 += pow(simrms_filter7[pix],2)  ;

	    }}



        sumflux_data_inside_filter1 = sumflux_data_inside_filter1/sum_weights_inside_filter1;
        sumflux_data_inside_filter2 = sumflux_data_inside_filter2/sum_weights_inside_filter2;
	sumflux_data_inside_filter3 = sumflux_data_inside_filter3/sum_weights_inside_filter3;
        sumflux_data_inside_filter4 = sumflux_data_inside_filter4/sum_weights_inside_filter4;
	sumflux_data_inside_filter5 = sumflux_data_inside_filter5/sum_weights_inside_filter5;
	sumflux_data_inside_filter6 = sumflux_data_inside_filter6/sum_weights_inside_filter6;
	sumflux_data_inside_filter7 = sumflux_data_inside_filter7/sum_weights_inside_filter7;

	sumflux_data_outside_filter1 = sumflux_data_outside_filter1/sum_weights_outside_filter1;
        sumflux_data_outside_filter2 = sumflux_data_outside_filter2/sum_weights_outside_filter2;
	sumflux_data_outside_filter3 = sumflux_data_outside_filter3/sum_weights_outside_filter3;
        sumflux_data_outside_filter4 = sumflux_data_outside_filter4/sum_weights_outside_filter4;
	sumflux_data_outside_filter5 = sumflux_data_outside_filter5/sum_weights_outside_filter5;
	sumflux_data_outside_filter6= sumflux_data_outside_filter6/sum_weights_outside_filter6;
	sumflux_data_outside_filter7 = sumflux_data_outside_filter7/sum_weights_outside_filter7;

        sumflux_data_galaxy_filter1 = sumflux_data_galaxy_filter1/sum_weights_galaxy_filter1;
        sumflux_data_galaxy_filter2 = sumflux_data_galaxy_filter2/sum_weights_galaxy_filter2;
	sumflux_data_galaxy_filter3 = sumflux_data_galaxy_filter1/sum_weights_galaxy_filter1;
        sumflux_data_galaxy_filter4 = sumflux_data_galaxy_filter2/sum_weights_galaxy_filter2;
	sumflux_data_galaxy_filter5 = sumflux_data_galaxy_filter5/sum_weights_galaxy_filter5;
	sumflux_data_galaxy_filter6 = sumflux_data_galaxy_filter6/sum_weights_galaxy_filter6;
	sumflux_data_galaxy_filter7 = sumflux_data_galaxy_filter7/sum_weights_galaxy_filter7;

	sumflux_model_inside_filter1 = sumflux_model_inside_filter1;
        sumflux_model_inside_filter2 = sumflux_model_inside_filter2;
	sumflux_model_inside_filter3 = sumflux_model_inside_filter3;
        sumflux_model_inside_filter4 = sumflux_model_inside_filter4;
	sumflux_model_inside_filter5 = sumflux_model_inside_filter5;
	sumflux_model_inside_filter6 = sumflux_model_inside_filter6;
	sumflux_model_inside_filter7 = sumflux_model_inside_filter7;

	sumflux_model_outside_filter1 = sumflux_model_outside_filter1;
        sumflux_model_outside_filter2 = sumflux_model_outside_filter2;
	sumflux_model_outside_filter3 = sumflux_model_outside_filter3;
        sumflux_model_outside_filter4 = sumflux_model_outside_filter4;
	sumflux_model_outside_filter5 = sumflux_model_outside_filter5;
	sumflux_model_outside_filter6 = sumflux_model_outside_filter6;
	sumflux_model_outside_filter7 = sumflux_model_outside_filter7;


        sumflux_model_galaxy_filter1 = sumflux_model_galaxy_filter1;
        sumflux_model_galaxy_filter2 = sumflux_model_galaxy_filter2;
	sumflux_model_galaxy_filter3 = sumflux_model_galaxy_filter3;
        sumflux_model_galaxy_filter4 = sumflux_model_galaxy_filter4;
	sumflux_model_galaxy_filter5 = sumflux_model_galaxy_filter5;
	sumflux_model_galaxy_filter6 = sumflux_model_galaxy_filter6;
	sumflux_model_galaxy_filter7 = sumflux_model_galaxy_filter7;


        sumrms_data_inside_filter1 = sqrt(sumrms_data_inside_filter1);
        sumrms_data_inside_filter2 = sqrt(sumrms_data_inside_filter2);
	sumrms_data_inside_filter3 = sqrt(sumrms_data_inside_filter3);
        sumrms_data_inside_filter4 = sqrt(sumrms_data_inside_filter4);


        sumrms_data_outside_filter1 = sqrt(sumrms_data_outside_filter1);
        sumrms_data_outside_filter2 = sqrt(sumrms_data_outside_filter2);
	sumrms_data_outside_filter3 = sqrt(sumrms_data_outside_filter3);
        sumrms_data_outside_filter4 = sqrt(sumrms_data_outside_filter4);
	sumrms_data_outside_filter5 = sqrt(sumrms_data_outside_filter5);
	sumrms_data_outside_filter6 = sqrt(sumrms_data_outside_filter6);
	sumrms_data_outside_filter7 = sqrt(sumrms_data_outside_filter7);

	sumrms_data_galaxy_filter1 = sqrt(sumrms_data_galaxy_filter1);
        sumrms_data_galaxy_filter2 = sqrt(sumrms_data_galaxy_filter2);
	sumrms_data_galaxy_filter3 = sqrt(sumrms_data_galaxy_filter3);
        sumrms_data_galaxy_filter4 = sqrt(sumrms_data_galaxy_filter4);
	sumrms_data_galaxy_filter5 = sqrt(sumrms_data_galaxy_filter5);
	sumrms_data_galaxy_filter6 = sqrt(sumrms_data_galaxy_filter6);
	sumrms_data_galaxy_filter7 = sqrt(sumrms_data_galaxy_filter7);


        mag_data_inside_filter1 = -2.5 * log10(sumflux_data_inside_filter1) + 31.426;
        mag_data_inside_filter2= -2.5 * log10(sumflux_data_inside_filter2) + 31.426;
	mag_data_inside_filter3 = -2.5 * log10(sumflux_data_inside_filter3) + 31.426;
        mag_data_inside_filter4= -2.5 * log10(sumflux_data_inside_filter4) + 31.426;
	mag_data_inside_filter5 = -2.5 * log10(sumflux_data_inside_filter5) + 31.426;
	mag_data_inside_filter6 = -2.5 * log10(sumflux_data_inside_filter6) + 31.426;
	mag_data_inside_filter7 = -2.5 * log10(sumflux_data_inside_filter7) + 31.426;

	mag_data_outside_filter1= -2.5 * log10(sumflux_data_outside_filter1) + 31.426;
        mag_data_outside_filter2= -2.5 * log10(sumflux_data_outside_filter2) + 31.426;
	mag_data_outside_filter3= -2.5 * log10(sumflux_data_outside_filter3) + 31.426;
        mag_data_outside_filter4= -2.5 * log10(sumflux_data_outside_filter4) + 31.426;
	mag_data_outside_filter5= -2.5 * log10(sumflux_data_outside_filter5) + 31.426;
	mag_data_outside_filter6= -2.5 * log10(sumflux_data_outside_filter6) + 31.426;
	mag_data_outside_filter7= -2.5 * log10(sumflux_data_outside_filter7) + 31.426;

        mag_data_galaxy_filter1 = -2.5 * log10(sumflux_data_galaxy_filter1) + 31.426;
        mag_data_galaxy_filter2= -2.5 * log10(sumflux_data_galaxy_filter2) + 31.426;
        mag_data_galaxy_filter3 = -2.5 * log10(sumflux_data_galaxy_filter3) + 31.426;
        mag_data_galaxy_filter4= -2.5 * log10(sumflux_data_galaxy_filter4) + 31.426;
        mag_data_galaxy_filter5 = -2.5 * log10(sumflux_data_galaxy_filter5) + 31.426;
        mag_data_galaxy_filter6 = -2.5 * log10(sumflux_data_galaxy_filter6) + 31.426;
        mag_data_galaxy_filter7 = -2.5 * log10(sumflux_data_galaxy_filter7) + 31.426;

	mag_model_inside_filter1 = -2.5 * log10(sumflux_model_inside_filter1) + 31.426;
        mag_model_inside_filter2= -2.5 * log10(sumflux_model_inside_filter2) + 31.426;
	mag_model_inside_filter3 = -2.5 * log10(sumflux_model_inside_filter3) + 31.426;
        mag_model_inside_filter4= -2.5 * log10(sumflux_model_inside_filter4) + 31.426;
	mag_model_inside_filter5 = -2.5 * log10(sumflux_model_inside_filter5) + 31.426;
	mag_model_inside_filter6 = -2.5 * log10(sumflux_model_inside_filter6) + 31.426;
	mag_model_inside_filter7 = -2.5 * log10(sumflux_model_inside_filter7) + 31.426;


        mag_model_outside_filter1= -2.5 * log10(sumflux_model_outside_filter1) + 31.426;
        mag_model_outside_filter2= -2.5 * log10(sumflux_model_outside_filter2) + 31.426;
	mag_model_outside_filter3= -2.5 * log10(sumflux_model_outside_filter3) + 31.426;
        mag_model_outside_filter4= -2.5 * log10(sumflux_model_outside_filter4) + 31.426;
	mag_model_outside_filter5= -2.5 * log10(sumflux_model_outside_filter5) + 31.426;
	mag_model_outside_filter6= -2.5 * log10(sumflux_model_outside_filter6) + 31.426;
	mag_model_outside_filter7= -2.5 * log10(sumflux_model_outside_filter7) + 31.426;

        mag_model_galaxy_filter1 = -2.5 * log10(sumflux_model_galaxy_filter1) + 31.426;
        mag_model_galaxy_filter2= -2.5 * log10(sumflux_model_galaxy_filter2) + 31.426;
	mag_model_galaxy_filter3 = -2.5 * log10(sumflux_model_galaxy_filter3) + 31.426;
        mag_model_galaxy_filter4= -2.5 * log10(sumflux_model_galaxy_filter4) + 31.426;
	mag_model_galaxy_filter5 = -2.5 * log10(sumflux_model_galaxy_filter5) + 31.426;
	mag_model_galaxy_filter6 = -2.5 * log10(sumflux_model_galaxy_filter6) + 31.426;
	mag_model_galaxy_filter7 = -2.5 * log10(sumflux_model_galaxy_filter7) + 31.426;


	mag_bulge_filter1 = -2.5 * log10(sumflux_bulge_filter1) + 31.426;
        mag_bulge_filter2 = -2.5 * log10(sumflux_bulge_filter2) + 31.426;
	mag_bulge_filter3 = -2.5 * log10(sumflux_bulge_filter3) + 31.426;
        mag_bulge_filter4 = -2.5 * log10(sumflux_bulge_filter4) + 31.426;
	mag_bulge_filter5 = -2.5 * log10(sumflux_bulge_filter5) + 31.426;
	mag_bulge_filter6 = -2.5 * log10(sumflux_bulge_filter6) + 31.426;
	mag_bulge_filter7 = -2.5 * log10(sumflux_bulge_filter7) + 31.426;

        mag_disk_filter1 = -2.5 * log10(sumflux_disk_filter1) + 31.426;
        mag_disk_filter2 = -2.5 * log10(sumflux_disk_filter2) + 31.426;
	mag_disk_filter3 = -2.5 * log10(sumflux_disk_filter3) + 31.426;
        mag_disk_filter4 = -2.5 * log10(sumflux_disk_filter4) + 31.426;
	mag_disk_filter5 = -2.5 * log10(sumflux_disk_filter5) + 31.426;
	mag_disk_filter6 = -2.5 * log10(sumflux_disk_filter6) + 31.426;
	mag_disk_filter7 = -2.5 * log10(sumflux_disk_filter7) + 31.426;

	btratio_filter1 = amp_bulge_filter1/(amp_bulge_filter1+amp_disk_filter1);
	btratio_filter2 = amp_bulge_filter2/(amp_bulge_filter2+amp_disk_filter2);
	btratio_filter3 = amp_bulge_filter3/(amp_bulge_filter3+amp_disk_filter3);
	btratio_filter4 = amp_bulge_filter4/(amp_bulge_filter4+amp_disk_filter4);
	btratio_filter5 = amp_bulge_filter5/(amp_bulge_filter5+amp_disk_filter5);
	btratio_filter6 = amp_bulge_filter6/(amp_bulge_filter6+amp_disk_filter6);
	btratio_filter7 = amp_bulge_filter7/(amp_bulge_filter7+amp_disk_filter7);


	//if (sumflux_data_inside_filter1 > 0.0 && sumflux_data_inside_filter2>0.0 && sumflux_data_outside_filter1>0.0 && sumflux_data_outside_filter2>0.0 && sumflux_model_inside_filter1>0.0 && sumflux_model_inside_filter2>0.0 && sumflux_model_outside_filter1>0.0 && sumflux_model_outside_filter2>0.0 && sumflux_data_inside_filter3 > 0.0 && sumflux_data_inside_filter4>0.0 && sumflux_data_outside_filter3>0.0 && sumflux_data_outside_filter4>0.0 && sumflux_model_inside_filter3>0.0 && sumflux_model_inside_filter4>0.0 && sumflux_model_outside_filter3>0.0 && sumflux_model_outside_filter4>0.0)
	//{

	   fprintf(cgfile_bulge_disk,"%d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",objid,ra,dec,zphot,dL,age_universe,arcsecstokpc,R50,VAB, IAB, M_I, itype_best_global, scaling_best_global, e1_best_global, e2_best_global, xcentre_best_global, ycentre_best_global, sumflux_bulge_filter1, sumflux_bulge_filter2, sumflux_bulge_filter3, sumflux_bulge_filter4, sumflux_bulge_filter5, sumflux_bulge_filter6, sumflux_bulge_filter7,sumflux_disk_filter1, sumflux_disk_filter2,  sumflux_disk_filter3, sumflux_disk_filter4, sumflux_disk_filter5, sumflux_disk_filter6, sumflux_disk_filter7, sumrms_data_galaxy_filter1, sumrms_data_galaxy_filter2, sumrms_data_galaxy_filter3, sumrms_data_galaxy_filter4, sumrms_data_galaxy_filter5, sumrms_data_galaxy_filter6, sumrms_data_galaxy_filter7, btratio_filter1, btratio_filter2, btratio_filter3, btratio_filter4, btratio_filter5, btratio_filter5, btratio_filter6, btratio_filter7, chi2_best_global/Nfreedom);  // 47 arguments


	   fprintf(cgfile_in_out,"%d %lf %lf %lf %lf %lf %lf %lf %lf %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", objid,ra,dec,zphot,arcsecstokpc,R50,VAB, IAB, M_I, sumflux_data_inside_filter1, sumflux_model_inside_filter1, sumflux_data_outside_filter1, sumflux_model_outside_filter1, sumflux_data_inside_filter2, sumflux_model_inside_filter2, sumflux_data_outside_filter2, sumflux_model_outside_filter2, sumflux_data_inside_filter3, sumflux_model_inside_filter3, sumflux_data_outside_filter3, sumflux_model_outside_filter3, sumflux_data_inside_filter4, sumflux_model_inside_filter4, sumflux_data_outside_filter4, sumflux_model_outside_filter4, sumflux_data_inside_filter5, sumflux_model_inside_filter5, sumflux_data_outside_filter5, sumflux_model_outside_filter5, sumflux_data_inside_filter6, sumflux_model_inside_filter6, sumflux_data_outside_filter6, sumflux_model_outside_filter6, sumflux_data_inside_filter7, sumflux_model_inside_filter7, sumflux_data_outside_filter7, sumflux_model_outside_filter7, sumrms_data_galaxy_filter1, sumrms_data_galaxy_filter2, sumrms_data_galaxy_filter3, sumrms_data_galaxy_filter4, sumrms_data_galaxy_filter5, sumrms_data_galaxy_filter6, sumrms_data_galaxy_filter7);  //45 arguments

	    fprintf(chi2file,"%d %10.7lf %10.7lf %lf %d %lf %lf\n",objid,ra,dec,zphot,itype_best_global,chi2_best_global/Nfreedom, chi2_masked/nfreedom_mask);


	   // }


	/* -------------------------------- */
	/* Write output maps as FITS images */
	/* -------------------------------- */
	 /* FITS file names */

	 char *rootname;
	 rootname = (char *) calloc(300, sizeof(char));
	 bzero(rootname, 300);

	 char *stampname;
	 stampname = (char *) calloc(300, sizeof(char));
	 double *otestout;

	 strcpy(stampname, rootname);
	 strcat(stampname, ".fits");

	 int anaxis = 2;
	 long fpixel[3] = { 1, 1, 1 }, lrow;
	 long anaxes[3];
	 anaxes[0] = Npix_side_acs;
	 anaxes[1] = Npix_side_acs;

	 int bitpix = -64;
	 fitsfile *tfptr;
	 char *pstr;
	 int status = 0;


	strcpy(outputmapspath,"results/maps_bvizyjh_dust/");

	strcpy(filename_model_conv_filter1,outputmapspath);
	strcat(filename_model_conv_filter1,objidstring);
	strcat(filename_model_conv_filter1,"model_conv_summed_B.fits");

	strcpy(filename_model_conv_filter2,outputmapspath);
        strcat(filename_model_conv_filter2,objidstring);
	strcat(filename_model_conv_filter2,"model_conv_summed_V.fits");

	strcpy(filename_model_conv_filter3,outputmapspath);
        strcat(filename_model_conv_filter3,objidstring);
	strcat(filename_model_conv_filter3,"model_conv_summed_I.fits");

	strcpy(filename_model_conv_filter4,outputmapspath);
        strcat(filename_model_conv_filter4,objidstring);
	strcat(filename_model_conv_filter4,"model_conv_summed_Z.fits");

	strcpy(filename_model_conv_filter5,outputmapspath);
        strcat(filename_model_conv_filter5,objidstring);
	strcat(filename_model_conv_filter5,"model_conv_summed_Y.fits");

	strcpy(filename_model_conv_filter6,outputmapspath);
        strcat(filename_model_conv_filter6,objidstring);
	strcat(filename_model_conv_filter6,"model_conv_summed_J.fits");

	strcpy(filename_model_conv_filter7,outputmapspath);
        strcat(filename_model_conv_filter7,objidstring);
	strcat(filename_model_conv_filter7,"model_conv_summed_H.fits");


	strcpy(filename_model_residual_filter1,outputmapspath);
	strcat(filename_model_residual_filter1,objidstring);
	strcat(filename_model_residual_filter1,"model_conv_residual_B.fits");

	strcpy(filename_model_residual_filter2,outputmapspath);
        strcat(filename_model_residual_filter2,objidstring);
	strcat(filename_model_residual_filter2,"model_conv_residual_V.fits");

	strcpy(filename_model_residual_filter3,outputmapspath);
	strcat(filename_model_residual_filter3,objidstring);
	strcat(filename_model_residual_filter3,"model_conv_residual_I.fits");

	strcpy(filename_model_residual_filter4,outputmapspath);
        strcat(filename_model_residual_filter4,objidstring);
	strcat(filename_model_residual_filter4,"model_conv_residual_Z.fits");

	strcpy(filename_model_residual_filter5,outputmapspath);
        strcat(filename_model_residual_filter5,objidstring);
	strcat(filename_model_residual_filter5,"model_conv_residual_Y.fits");

	strcpy(filename_model_residual_filter6,outputmapspath);
	strcat(filename_model_residual_filter6,objidstring);
	strcat(filename_model_residual_filter6,"model_conv_residual_J.fits");

	strcpy(filename_model_residual_filter7,outputmapspath);
        strcat(filename_model_residual_filter7,objidstring);
	strcat(filename_model_residual_filter7,"model_conv_residual_H.fits");



	strcpy(filename_mask_filter1,outputmapspath);
	strcat(filename_mask_filter1,objidstring);
	strcat(filename_mask_filter1,"mask_B.fits");

	strcpy(filename_mask_filter2,outputmapspath);
        strcat(filename_mask_filter2,objidstring);
	strcat(filename_mask_filter2,"mask_V.fits");

	strcpy(filename_mask_filter3,outputmapspath);
	strcat(filename_mask_filter3,objidstring);
	strcat(filename_mask_filter3,"mask_I.fits");

	strcpy(filename_mask_filter4,outputmapspath);
        strcat(filename_mask_filter4,objidstring);
	strcat(filename_mask_filter4,"mask_Z.fits");

	strcpy(filename_mask_filter5,outputmapspath);
        strcat(filename_mask_filter5,objidstring);
	strcat(filename_mask_filter5,"mask_Y.fits");

	strcpy(filename_mask_filter6,outputmapspath);
	strcat(filename_mask_filter6,objidstring);
	strcat(filename_mask_filter6,"mask_J.fits");

	strcpy(filename_mask_filter7,outputmapspath);
        strcat(filename_mask_filter7,objidstring);
	strcat(filename_mask_filter7,"mask_H.fits");


	    fits_create_file(&tfptr, filename_model_conv_filter1, &status);
            fits_create_img(tfptr, bitpix, anaxis, anaxes, &status);
            if (status) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }
            if (fits_write_pix
                (tfptr, TDOUBLE, fpixel, npixels_acs, model_conv_summed_filter1, &status)) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
		       stampname);
                    fits_report_error(stderr, status);
            }

            fits_close_file(tfptr, &status);

            if (status) {
              fits_report_error(stderr, status);
              exit(2);
            }


	    fits_create_file(&tfptr,filename_model_conv_filter2, &status);
            fits_create_img(tfptr, bitpix, anaxis, anaxes, &status);
            if (status) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }
            if (fits_write_pix
                (tfptr, TDOUBLE, fpixel, npixels_acs, model_conv_summed_filter2, &status)) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }

            fits_close_file(tfptr, &status);

            if (status) {
              fits_report_error(stderr, status);
              exit(2);
            }

	       fits_create_file(&tfptr,filename_model_conv_filter3, &status);
            fits_create_img(tfptr, bitpix, anaxis, anaxes, &status);
            if (status) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }
            if (fits_write_pix
                (tfptr, TDOUBLE, fpixel, npixels_acs, model_conv_summed_filter3, &status)) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }

            fits_close_file(tfptr, &status);

            if (status) {
              fits_report_error(stderr, status);
              exit(2);
            }


	    fits_create_file(&tfptr,filename_model_conv_filter4, &status);
            fits_create_img(tfptr, bitpix, anaxis, anaxes, &status);
            if (status) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }
            if (fits_write_pix
                (tfptr, TDOUBLE, fpixel, npixels_acs, model_conv_summed_filter4, &status)) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }

            fits_close_file(tfptr, &status);

            if (status) {
              fits_report_error(stderr, status);
              exit(2);
            }



	    fits_create_file(&tfptr, filename_model_residual_filter1, &status);
	    fits_create_img(tfptr, bitpix, anaxis, anaxes, &status);
	    if (status) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
	    }
	    if (fits_write_pix
		(tfptr, TDOUBLE, fpixel, npixels_acs, model_conv_residual_filter1, &status)) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
	    }

	    fits_close_file(tfptr, &status);

	    if (status) {
	      fits_report_error(stderr, status);
	      exit(2);
	    }


	    fits_create_file(&tfptr, filename_model_residual_filter2, &status);
            fits_create_img(tfptr, bitpix, anaxis, anaxes, &status);
            if (status) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }
            if (fits_write_pix
                (tfptr, TDOUBLE, fpixel, npixels_acs, model_conv_residual_filter2, &status)) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }

            fits_close_file(tfptr, &status);

            if (status) {
              fits_report_error(stderr, status);
              exit(2);
            }


        fits_create_file(&tfptr, filename_model_residual_filter3, &status);
	    fits_create_img(tfptr, bitpix, anaxis, anaxes, &status);
	    if (status) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
	    }
	    if (fits_write_pix
		(tfptr, TDOUBLE, fpixel, npixels_acs, model_conv_residual_filter3, &status)) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
	    }

	    fits_close_file(tfptr, &status);

	    if (status) {
	      fits_report_error(stderr, status);
	      exit(2);
	    }


	    fits_create_file(&tfptr, filename_model_residual_filter4, &status);
            fits_create_img(tfptr, bitpix, anaxis, anaxes, &status);
            if (status) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }
            if (fits_write_pix
                (tfptr, TDOUBLE, fpixel, npixels_acs, model_conv_residual_filter4, &status)) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }

            fits_close_file(tfptr, &status);

            if (status) {
              fits_report_error(stderr, status);
              exit(2);
            }



	    fits_create_file(&tfptr, filename_mask_filter1, &status);
            fits_create_img(tfptr, bitpix, anaxis, anaxes, &status);
            if (status) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }
            if (fits_write_pix
                (tfptr, TDOUBLE, fpixel, npixels_acs, mask_filter1, &status)) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }

            fits_close_file(tfptr, &status);

            if (status) {
              fits_report_error(stderr, status);
              exit(2);
            }


	    fits_create_file(&tfptr, filename_mask_filter2, &status);
            fits_create_img(tfptr, bitpix, anaxis, anaxes, &status);
            if (status) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }
            if (fits_write_pix
                (tfptr, TDOUBLE, fpixel, npixels_acs, mask_filter2, &status)) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }

            fits_close_file(tfptr, &status);

            if (status) {
              fits_report_error(stderr, status);
              exit(2);
            }


	    fits_create_file(&tfptr, filename_mask_filter3, &status);
            fits_create_img(tfptr, bitpix, anaxis, anaxes, &status);
            if (status) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }
            if (fits_write_pix
                (tfptr, TDOUBLE, fpixel, npixels_acs, mask_filter3, &status)) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }

            fits_close_file(tfptr, &status);

            if (status) {
              fits_report_error(stderr, status);
              exit(2);
            }



	    fits_create_file(&tfptr, filename_mask_filter4, &status);
            fits_create_img(tfptr, bitpix, anaxis, anaxes, &status);
            if (status) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }
            if (fits_write_pix
                (tfptr, TDOUBLE, fpixel, npixels_acs, mask_filter4, &status)) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }

            fits_close_file(tfptr, &status);

            if (status) {
              fits_report_error(stderr, status);
              exit(2);
            }


	    /* WFC3 */

	    anaxes[0] = Npix_side_wfc3;
	    anaxes[1] = Npix_side_wfc3;


	    fits_create_file(&tfptr,filename_model_conv_filter5, &status);
            fits_create_img(tfptr, bitpix, anaxis, anaxes, &status);
            if (status) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }
            if (fits_write_pix
                (tfptr, TDOUBLE, fpixel, npixels_wfc3, model_conv_summed_filter5, &status)) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }

            fits_close_file(tfptr, &status);

            if (status) {
              fits_report_error(stderr, status);
              exit(2);
            }


	     fits_create_file(&tfptr,filename_model_conv_filter6, &status);
            fits_create_img(tfptr, bitpix, anaxis, anaxes, &status);
            if (status) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }
            if (fits_write_pix
                (tfptr, TDOUBLE, fpixel, npixels_wfc3, model_conv_summed_filter6, &status)) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }

            fits_close_file(tfptr, &status);

            if (status) {
              fits_report_error(stderr, status);
              exit(2);
            }


	     fits_create_file(&tfptr,filename_model_conv_filter7, &status);
            fits_create_img(tfptr, bitpix, anaxis, anaxes, &status);
            if (status) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }
            if (fits_write_pix
                (tfptr, TDOUBLE, fpixel, npixels_wfc3, model_conv_summed_filter7, &status)) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }

            fits_close_file(tfptr, &status);

            if (status) {
              fits_report_error(stderr, status);
              exit(2);
            }



	    fits_create_file(&tfptr,filename_mask_filter5, &status);
            fits_create_img(tfptr, bitpix, anaxis, anaxes, &status);
            if (status) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }
            if (fits_write_pix
                (tfptr, TDOUBLE, fpixel, npixels_wfc3, mask_filter5, &status)) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }

            fits_close_file(tfptr, &status);

            if (status) {
              fits_report_error(stderr, status);
              exit(2);
            }

	    	    fits_create_file(&tfptr,filename_mask_filter6, &status);
            fits_create_img(tfptr, bitpix, anaxis, anaxes, &status);
            if (status) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }
            if (fits_write_pix
                (tfptr, TDOUBLE, fpixel, npixels_wfc3, mask_filter6, &status)) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }

            fits_close_file(tfptr, &status);

            if (status) {
              fits_report_error(stderr, status);
              exit(2);
            }

	    	    fits_create_file(&tfptr,filename_mask_filter7, &status);
            fits_create_img(tfptr, bitpix, anaxis, anaxes, &status);
            if (status) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }
            if (fits_write_pix
                (tfptr, TDOUBLE, fpixel, npixels_wfc3, mask_filter7, &status)) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
            }

            fits_close_file(tfptr, &status);

            if (status) {
              fits_report_error(stderr, status);
              exit(2);
            }




	     fits_create_file(&tfptr, filename_model_residual_filter5, &status);
	    fits_create_img(tfptr, bitpix, anaxis, anaxes, &status);
	    if (status) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
	    }
	    if (fits_write_pix
		(tfptr, TDOUBLE, fpixel, npixels_wfc3, model_conv_residual_filter5, &status)) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
	    }

	    fits_close_file(tfptr, &status);

	    if (status) {
	      fits_report_error(stderr, status);
	      exit(2);
	    }


	     fits_create_file(&tfptr, filename_model_residual_filter6, &status);
	    fits_create_img(tfptr, bitpix, anaxis, anaxes, &status);
	    if (status) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
	    }
	    if (fits_write_pix
		(tfptr, TDOUBLE, fpixel, npixels_wfc3, model_conv_residual_filter6, &status)) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
	    }

	    fits_close_file(tfptr, &status);

	    if (status) {
	      fits_report_error(stderr, status);
	      exit(2);
	    }


	     fits_create_file(&tfptr, filename_model_residual_filter7, &status);
	    fits_create_img(tfptr, bitpix, anaxis, anaxes, &status);
	    if (status) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
	    }
	    if (fits_write_pix
		(tfptr, TDOUBLE, fpixel, npixels_wfc3, model_conv_residual_filter7, &status)) {
                    printf
                      (" error writing postagestamp test pixel data %s \n",
                       stampname);
                    fits_report_error(stderr, status);
	    }

	    fits_close_file(tfptr, &status);

	    if (status) {
	      fits_report_error(stderr, status);
	      exit(2);
	    }



	obji++;

       } // loop over object catalog


   fclose(fpcat);
   fclose(cgfile_bulge_disk);
   fclose(cgfile_in_out);
   fclose(chi2file);

  /* ------------------------------------------- */
  /* free memory*/
  /* ------------------------------------------- */


  free(x_disk);
  free(y_disk);
  free(z_disk);
  free(x_bulge);
  free(y_bulge);
  free(z_bulge);
  free(I_xyz_disk_filter1);
  free(I_xyz_bulge_filter1);


  free(model_bulge_filter1);
  free(model_disk_filter1);
  free(data_filter1);
  free(model_bulge_filter2);
  free(model_disk_filter2);
  free(data_filter2);
  free(model_bulge_filter3);
  free(model_disk_filter3);
  free(data_filter3);
  free(model_bulge_filter4);
  free(model_disk_filter4);
  free(data_filter4);
  free(model_bulge_filter5);
  free(model_disk_filter5);
  free(data_filter5);
  free(model_bulge_filter6);
  free(model_disk_filter6);
  free(data_filter6);
  free(model_bulge_filter7);
  free(model_disk_filter7);
  free(data_filter7);

  free(model_conv_disk_filter1);
  free(model_conv_bulge_filter1);
  free(model_conv_summed_filter1);
  free(model_conv_disk_filter2);
  free(model_conv_bulge_filter2);
  free(model_conv_summed_filter2);
  free(model_conv_disk_filter3);
  free(model_conv_bulge_filter3);
  free(model_conv_summed_filter3);
  free(model_conv_disk_filter4);
  free(model_conv_bulge_filter4);
  free(model_conv_summed_filter4);
  free(model_conv_disk_filter5);
  free(model_conv_bulge_filter5);
  free(model_conv_summed_filter5);
  free(model_conv_disk_filter6);
  free(model_conv_bulge_filter6);
  free(model_conv_summed_filter6);
  free(model_conv_disk_filter7);
  free(model_conv_bulge_filter7);
  free(model_conv_summed_filter7);

  free(model_unconv_disk_filter1);
  free(model_unconv_bulge_filter1);
  free(model_unconv_summed_filter1);
  free(model_conv_residual_filter1);
  free(model_unconv_disk_filter2);
  free(model_unconv_bulge_filter2);
  free(model_unconv_summed_filter2);
  free(model_conv_residual_filter2);
  free(model_unconv_disk_filter3);
  free(model_unconv_bulge_filter3);
  free(model_unconv_summed_filter3);
  free(model_conv_residual_filter3);
  free(model_unconv_disk_filter4);
  free(model_unconv_bulge_filter4);
  free(model_unconv_summed_filter4);
  free(model_conv_residual_filter4);
  free(model_unconv_disk_filter5);
  free(model_unconv_bulge_filter5);
  free(model_unconv_summed_filter5);
  free(model_conv_residual_filter5);
  free(model_unconv_disk_filter6);
  free(model_unconv_bulge_filter6);
  free(model_unconv_summed_filter6);
  free(model_conv_residual_filter6);
  free(model_unconv_disk_filter7);
  free(model_unconv_bulge_filter7);
  free(model_unconv_summed_filter7);
  free(model_conv_residual_filter7);

  free(cum_model_flux_radius_filter1);
  free(cum_data_flux_radius_filter1);
  free(cum_data_err_radius_filter1);
  free(cum_model_flux_radius_filter2);
  free(cum_data_flux_radius_filter2);
  free(cum_data_err_radius_filter2);
  free(cum_model_flux_radius_filter3);
  free(cum_data_flux_radius_filter3);
  free(cum_data_err_radius_filter3);
  free(cum_model_flux_radius_filter4);
  free(cum_data_flux_radius_filter4);
  free(cum_data_err_radius_filter4);
  free(cum_model_flux_radius_filter5);
  free(cum_data_flux_radius_filter5);
  free(cum_data_err_radius_filter5);
  free(cum_model_flux_radius_filter6);
  free(cum_data_flux_radius_filter6);
  free(cum_data_err_radius_filter6);
  free(cum_model_flux_radius_filter7);
  free(cum_data_flux_radius_filter7);
  free(cum_data_err_radius_filter7);

  free(simdata_filter1);
  free(simdata_filter2);
  free(simdata_filter3);
  free(simdata_filter4);
  free(simdata_filter5);
  free(simdata_filter6);
  free(simdata_filter7);

  free(simrms_filter1);
  free(simrms_filter2);
  free(simrms_filter3);
  free(simrms_filter4);
  free(simrms_filter5);
  free(simrms_filter6);
  free(simrms_filter7);

  free(simdata_orig_filter1);
  free(simdata_orig_filter2);
  free(simdata_orig_filter3);
  free(simdata_orig_filter4);
  free(simdata_orig_filter5);
  free(simdata_orig_filter6);
  free(simdata_orig_filter7);

  free(simnoise_filter1);
  free(simnoise_filter2);
  free(simnoise_filter3);
  free(simnoise_filter4);
  free(simnoise_filter5);
  free(simnoise_filter6);
  free(simnoise_filter7);

  free(weight_filter1);
  free(weight_filter2);
  free(weight_filter3);
  free(weight_filter4);
  free(weight_filter5);
  free(weight_filter6);
  free(weight_filter7);

  free(wht_filter1);
  free(wht_filter2);
  free(wht_filter3);
  free(wht_filter4);
  free(wht_filter5);
  free(wht_filter6);
  free(wht_filter7);

  free(mask_filter1);
  free(mask_filter2);
  free(mask_filter3);
  free(mask_filter4);
  free(mask_filter5);
  free(mask_filter6);
  free(mask_filter7);


  free(psf_image_filter1);
  free(psf_image_filter2);
  free(psf_image_filter3);
  free(psf_image_filter4);
  free(psf_image_filter5);
  free(psf_image_filter6);
  free(psf_image_filter7);

  free(filename_model_unconv_filter1);
  free(filename_model_conv_filter1);
  free(filename_model_residual_filter1);
  free(filename_model_unconv_filter2);
  free(filename_model_conv_filter2);
  free(filename_model_residual_filter2);
  free(filename_model_unconv_filter3);
  free(filename_model_conv_filter3);
  free(filename_model_residual_filter3);
  free(filename_model_unconv_filter4);
  free(filename_model_conv_filter4);
  free(filename_model_residual_filter4);
  free(filename_model_unconv_filter5);
  free(filename_model_conv_filter5);
  free(filename_model_residual_filter5);
  free(filename_model_unconv_filter6);
  free(filename_model_conv_filter6);
  free(filename_model_residual_filter6);
  free(filename_model_unconv_filter7);
  free(filename_model_conv_filter7);
  free(filename_model_residual_filter7);

  free(filename_data_filter1);
  free(filename_data_filter2);
  free(filename_data_filter3);
  free(filename_data_filter4);
  free(filename_data_filter5);
  free(filename_data_filter6);
  free(filename_data_filter7);

  free(datafitsname_filter1);
  free(datafitsname_filter2);
  free(datafitsname_filter3);
  free(datafitsname_filter4);
  free(datafitsname_filter5);
  free(datafitsname_filter6);
  free(datafitsname_filter7);

  free(whtfitsname_filter1);
  free(whtfitsname_filter2);
  free(whtfitsname_filter3);
  free(whtfitsname_filter4);
  free(whtfitsname_filter5);
  free(whtfitsname_filter6);
  free(whtfitsname_filter7);

  free(psffitsname_filter1);
  free(psffitsname_filter2);
  free(psffitsname_filter3);
  free(psffitsname_filter4);
  free(psffitsname_filter5);
  free(psffitsname_filter6);
  free(psffitsname_filter7);


  free(catfilename);
  free(catfilename_bulge_disk);
  free(catfilename_in_out);


  iii=0;
  for (iii = 0; iii < npoints_bulge; iii++)
    free(x_bulge_type[iii]);
  free(x_bulge_type);
  for (iii = 0; iii < npoints_bulge; iii++)
    free(y_bulge_type[iii]);
  free(y_bulge_type);
  for  (iii = 0; iii < npoints_bulge; iii++)
    free(z_bulge_type[iii]);
  free(z_bulge_type);
for (iii = 0; iii < npoints_bulge; iii++)
    free(I_xyz_bulge_filter1_type[iii]);
  free(I_xyz_bulge_filter1_type);

  for (iii = 0; iii < npoints_disk; iii++)
    free(x_disk_type[iii]);
  free(x_disk_type);
  for (iii = 0; iii < npoints_disk; iii++)
    free(y_disk_type[iii]);
  free(y_disk_type);
  for (iii = 0; iii < npoints_disk; iii++)
    free(z_disk_type[iii]);
  free(z_disk_type);
  for (iii = 0; iii < npoints_disk; iii++)
    free(I_xyz_disk_filter1_type[iii]);
  free(I_xyz_disk_filter1_type);

  free(model_unconv_acs);
  free(psf_image_acs);
  free(model_conv_acs);
  free(model_unconv_wfc3);
  free(psf_image_wfc3);
  free(model_conv_wfc3);

  free(apix);
  free(psfpix_filter1);
  free(psfpix_filter2);
  free(psfpix_filter3);
  free(psfpix_filter4);
  free(psfpix_filter5);
  free(psfpix_filter6);
  free(psfpix_filter7);

  fftw_destroy_plan(plan1_acs);
  fftw_destroy_plan(plan2_acs);
  fftw_destroy_plan(plan3_acs);
  fftw_free(in_acs);
  fftw_free(identity_acs);
  fftw_free(final_acs);
  fftw_free(inTrans_acs);
  fftw_free(identityTrans_acs);
  fftw_free(FinalFFT_acs);

  fftw_destroy_plan(plan1_wfc3);
  fftw_destroy_plan(plan2_wfc3);
  fftw_destroy_plan(plan3_wfc3);
  fftw_free(in_wfc3);
  fftw_free(identity_wfc3);
  fftw_free(final_wfc3);
  fftw_free(inTrans_wfc3);
  fftw_free(identityTrans_wfc3);
  fftw_free(FinalFFT_wfc3);

  free(imagename);
  free(galaxytype);
  free(modelfilename_disk);
  free(modelfilename_bulge);


  return 0;

 }



