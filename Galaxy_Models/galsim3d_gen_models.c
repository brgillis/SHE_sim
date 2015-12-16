#include <pthread.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <string.h>
#include <unistd.h>
#include <fitsio.h>
#include <complex.h>
#include <fftw3.h>

/*gsl header files */
#include </opt/local/include/gsl/gsl_integration.h>
#include </opt/local/include/gsl/gsl_spline.h>
#include </opt/local/include/gsl/gsl_linalg.h>
#include </opt/local/include/gsl/gsl_errno.h>
#include </opt/local/include/gsl/gsl_matrix.h>
#include </opt/local/include/gsl/gsl_rng.h>


//CCfits libraries
#ifdef _MSC_VER
#include "MSconfig.h" // for truncation warning
#endif

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

/* general constants  */
#define hubble 0.7
#define PI 3.1415926
#define c 3e8 // ms-1
#define h 6.64e-27  //ergs s

/* sampling points from distributions */
/*radial */

#define zcut_z0_ratio 4.5  // z_cut/z0 (?) 
#define rcut_r0_ratio 4.5 // r_cut/r0  (van de Kruit & Searl 1981)

//#define rpoints 5001 /* to generate the profile to start with */
#define rpoints_inf 100001  /*for integration */
#define kpc_unit 0.01//0.0005  // kpc in the simulation ***HAS TO BE LESS THAN pixelscale_kpc of instrument that's input ****
/*z*/
//#define zpoints 5001   /* to generate the profile to start with */
#define zpoints_inf 100001 /* for integration */


/*routines */

double cum_integral(double r[], double rstart, double rmax, double y[], int npoints);   // find cumulative integral at each r point
double convolvesed(char sedfile[],char filterfile[],double);


int main(int argc, char *argv[])
{

  FILE *filterlist;
  FILE *filtername;
  FILE *newfilter;
  FILE *filterfluxes;


  /* for random number */
  /*--------------------- */
  const gsl_rng_type * T;
  gsl_rng * random;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  random = gsl_rng_alloc (T);
  /*--------------------- */
 
  if (argc!=8) {
    printf("USAGE type r0_bulge r0_disk z0_disk bulge_compression Npoints n_Sersic");
    exit(0);
  }
  

  FILE *f_r_file; 
  f_r_file = fopen("func-r.dat", "w");
  FILE *cum_f_r_file; 
  cum_f_r_file = fopen("cumfunc-r.dat", "w");
  FILE *file2d;
  file2d=fopen("xypoints.dat","w");

  FILE *f_z_file; 
  f_z_file = fopen("func-z.dat", "w");
  FILE *cum_f_z_file; 
  cum_f_z_file = fopen("cumfunc-z.dat", "w");
 
  /* first x-y-z points and flux */
  FILE *file3d;
  
  /* ---------------------------*/
  /* Declare variables */
  /* ---------------------------*/
 
  int ir, npixels, F_nintervals, zpoints, rpoints;
  double nsersic,bn;
  double *nparticles,*radius, *I_r, *F_r, *radius_inf, *I_r_inf,F_r_inf;
  int *npoints_radshell, *npoints_zshell;
  double *z, *z_new, *I_z, *F_z, *F_z_new, *z_inf, *I_z_inf, F_z_inf,*intensity; 
  double I_z0 = 1.0; 
  double redshift, d_L, r0bulge_kpc_i, r0disk_kpc_i, z0_kpc_i, squeeze_bulge_zfactor;
  double I_disk_r0_i, I_bulge_r0_i;
  char* galaxytype= malloc(5*sizeof(char)); 
  char* sedfilename = malloc(12*sizeof(char));
  
  double *z_store, *I_z_new, *I_r_new, *I_xyz; 
  double F_interval, F_new;
  double *radius_store, *theta_store, *phi_store, *x, *y; 
  int *zindex;
  
   
   I_disk_r0_i = 1.0; 
   I_bulge_r0_i = 1.0; 
   
 
  /* ---------------------------*/
  /* Read parameters for generating galaxy  */
  /* ---------------------------*/
   
   strcpy(galaxytype,argv[1]);
   r0bulge_kpc_i = atof(argv[2]); /* bulge Re in kpc */
   r0disk_kpc_i = atof(argv[3]); /* disk scale length kpc  */
   z0_kpc_i = atof(argv[4]);
   squeeze_bulge_zfactor = atof(argv[5]); /* for bulges, squeeze factor */
   F_nintervals = atoi(argv[6]);
   nsersic = atof(argv[7]);

   bn = (1.9992 * nsersic) - 0.3271 ; // Capaccioli 1989
  
   /*  Cut off radii in r and z */
   if (strncmp(galaxytype,"disk",4)==0) {
     rpoints = round(rcut_r0_ratio * r0disk_kpc_i/kpc_unit) +1   ;  }
  if (strncmp(galaxytype,"bulge",5)==0) {
     rpoints = round(rcut_r0_ratio * r0bulge_kpc_i/kpc_unit) + 1  ;  }

  //rpoints = 10000;  // hack NW 
    
    zpoints = rpoints;
   //zpoints = round(zcut_z0_ratio * z0_kpc_i/kpc_unit)    ; 

   printf("rpoints: %d\n",rpoints);
   printf("zpoints: %d\n",zpoints);

   radius = (double *) malloc ((rpoints) * sizeof(double));
   radius_inf = (double *) malloc ((rpoints_inf) * sizeof(double));
   I_r = (double *) malloc ((rpoints) * sizeof(double));
   I_r_inf = (double *) malloc ((rpoints_inf) * sizeof(double));
   F_r = (double *) malloc ((rpoints) * sizeof(double));

   z = (double *) malloc ((zpoints) * sizeof(double));
   z_inf = (double *) malloc ((zpoints_inf) * sizeof(double));
   I_z = (double *) malloc ((zpoints) * sizeof(double));
   I_z_inf = (double *) malloc ((zpoints_inf) * sizeof(double));
   F_z = (double *) malloc ((zpoints) * sizeof(double));
 
   npoints_radshell = (int *) malloc ((rpoints) * sizeof(int));
   npoints_zshell =  (int *) malloc ((rpoints) * sizeof(int));
   
   z_store = (double *) malloc ((F_nintervals) * sizeof(double));
   I_z_new = (double *) malloc ((F_nintervals) * sizeof(double));
   zindex = (int *) malloc ((F_nintervals) * sizeof(int));
   radius_store = (double *) malloc ((F_nintervals) * sizeof(double));
   theta_store = (double *) malloc ((F_nintervals) * sizeof(double));
   phi_store = (double *) malloc ((F_nintervals) * sizeof(double));
   x = (double *) malloc ((F_nintervals) * sizeof(double));
   y = (double *) malloc ((F_nintervals) * sizeof(double));
   I_r_new = (double *) malloc ((F_nintervals) * sizeof(double));
   I_xyz = (double *) malloc ((F_nintervals) * sizeof(double));


  char outfilename[80];
  strcat(outfilename, "3dmodel_");
  strcat(outfilename,galaxytype);
  strcat(outfilename,".dat");
  
  file3d=fopen(outfilename,"w");


/* ---------------------------*/
/* Start the simulation  */
/* ---------------------------*/

/* ---------------------------*/
  /* Set starting grid and input mass profiles*/
/* ---------------------------*/

  for (ir=0;ir<rpoints;ir++) {
    radius[ir] = ir * kpc_unit; 

    if (strncmp(galaxytype,"disk",4)==0) {
      I_r[ir] = I_disk_r0_i * exp(-radius[ir]/r0disk_kpc_i); } // r0_disk_pic_i is the SCALE LENGTH OF THE DISK (corresponds to h in Graham & Worley 2008)
    if (strncmp(galaxytype,"bulge",5)==0) {
      //I_r[ir] = I_bulge_r0_i * exp( -log(10) * 0.4 * 1.086 * bn * ( pow( (radius[ir]/r0bulge_kpc_i),1/nsersic) -1 ) ) ; }  // r0_bulge_kpc_i is the Effective radius/half light radius of the bulge (r_e in Graham & Worely 2008)
      I_r[ir] = I_bulge_r0_i * exp( -bn * ( pow( (radius[ir]/r0bulge_kpc_i),1/nsersic) -1 ) ) ; }  // r0_bulge_kpc_i is the Effective radius/half light radius of the bulge (r_e in Graham & Worely 2008)
      fprintf(f_r_file, "%lf %lf\n",radius[ir], I_r[ir]);
  }


/* ---------------------------*/
  /* evaluate cumulative pdf F(r) */
  /* Integrate the pdfs p(r): evaluate indefinite integral F_r_inf */
  /* For disk, p(r) dr = I(r) 2 pi r dr*/
  /* For bulge, p(r) dr = I(r) 4 pi r^2 dr*/  
/* ---------------------------*/
    
  
   for (ir=0;ir<rpoints_inf;ir++) {
    radius_inf[ir] = ir * kpc_unit; 
    if (strncmp(galaxytype,"disk",4)==0) {
      I_r_inf[ir] = I_disk_r0_i * exp(-radius_inf[ir]/r0disk_kpc_i) * 2 * PI * radius_inf[ir]; } 
    if (strncmp(galaxytype,"bulge",5)==0) {
      //I_r_inf[ir] = I_bulge_r0_i * exp( -log(10) * 0.4 * 1.086 * bn * ( pow( (radius_inf[ir]/r0bulge_kpc_i),1/nsersic) -1 ) ) * 2 * PI * radius_inf[ir] ;
      I_r_inf[ir] = I_bulge_r0_i * exp( -bn * ( pow( (radius_inf[ir]/r0bulge_kpc_i),1/nsersic) -1 ) ) * 2 * PI * radius_inf[ir] ;
    }

     }

   F_r_inf = cum_integral(radius_inf, 0.0, radius_inf[rpoints_inf-1], I_r_inf, rpoints_inf); 
  

   for (ir=0;ir<rpoints;ir++) {
     if (strncmp(galaxytype,"disk",4)==0) {
      I_r[ir] = I_disk_r0_i * exp(-radius[ir]/r0disk_kpc_i) * 2 * PI * radius[ir]; }
     if (strncmp(galaxytype,"bulge",5)==0) {
       //I_r[ir] = I_bulge_r0_i * exp( -log(10) * 0.4 * 1.086 * bn * ( pow( (radius[ir]/r0bulge_kpc_i),1/nsersic) -1 ) ) * 2 * PI * radius[ir] ;
       I_r[ir] = I_bulge_r0_i * exp( - bn * ( pow( (radius[ir]/r0bulge_kpc_i),1/nsersic) -1 ) ) * 2 * PI * radius[ir];
     }
     F_r[ir] = cum_integral(radius, 0.0, radius[ir], I_r,rpoints)/F_r_inf; 
     fprintf(cum_f_r_file, "%lf %lf\n", radius[ir], F_r[ir]);
 
  }

/* -----------------------------------*/
/* Now draw radius points from equally spaced intervals in cumulative pdf F(r) */
/* new grid equally spaced in F(r) from 0 to 1 */
/* -----------------------------------*/

  F_interval = (F_r[rpoints-1] - F_r[0])/F_nintervals;

  int ifr;
  double sumflux_bulge = 0, sumflux_disk ;
   
  for (ifr=0;ifr<F_nintervals;ifr++) {
      /* choose mid-point of the integral */
    
      F_new = (F_interval * (0.5 + ifr)) ;
      //F_new =  ((double)rand() / (double)RAND_MAX) ;
   
      /* Use high resolution grid of F_r_new and get the r_new at this F_new value */
   
      for (ir = 0; ir < rpoints; ir++)
	{
 	 
	  if ( (F_new >= F_r[ir]) && (F_new < F_r[ir+1]) ) 
	    {
	      radius_store[ifr] = radius[ir] ;// store that r value
	      theta_store[ifr] =  gsl_rng_uniform (random) * 2 * PI;  // randomly chose a theta for that r value from 0 to 2pi radians
	      phi_store[ifr] =  acos( (2* gsl_rng_uniform (random)) -1);  // cos(phi) from -1 to +1 -> phi from 0 to PI radians 
	      	      
	      /* get corresponding x and y */
	      if (strncmp(galaxytype,"disk",4)==0) {

	      x[ifr] = radius_store[ifr] * cos(theta_store[ifr]);
	      y[ifr] = radius_store[ifr] * sin(theta_store[ifr]);
	      I_r_new[ifr] = (I_disk_r0_i * exp(-radius_store[ifr])/r0disk_kpc_i); 
		
	      }

	      if (strncmp(galaxytype,"bulge",5)==0) 
		  {
		 
		    x[ifr] = radius_store[ifr] * cos(theta_store[ifr]) ;//; * sin(phi_store[ifr]) ;
		    y[ifr] = radius_store[ifr] * sin(theta_store[ifr]); //;* sin(phi_store[ifr]) ;
		    //z_store[ifr] = squeeze_bulge_zfactor * radius_store[ifr] * cos(phi_store[ifr]) ;			     
		    //I_r_new[ifr] = I_bulge_r0_i * exp( -7.6692 * ( pow( (radius_store[ifr]/r0bulge_kpc_i),0.25) -1 ) ); 
		    //I_r_new[ifr] = ( I_bulge_r0_i * exp( -log(10) * 0.4 * 1.086 * bn * ( pow( (radius_store[ifr]/r0bulge_kpc_i),1/nsersic) -1 ) ) ) ; 
		    I_r_new[ifr] = I_bulge_r0_i * exp( -bn * ( pow( (radius_store[ifr]/r0bulge_kpc_i),1/nsersic) -1 ) )  ; 
		    I_xyz[ifr] = I_r_new[ifr]; // * 1e23 * (1+redshift) * L_solar/(4 * PI * 3.08568025e24 * 3.08568025e24 * d_L * d_L) ;

		    sumflux_bulge += I_xyz[ifr];
		    

		  }
	
	    }

	}


  }


if (strncmp(galaxytype,"bulge",5)==0) 
{
  for (ifr=0;ifr<F_nintervals;ifr++)
    {

      I_xyz[ifr] = I_xyz[ifr]/sumflux_bulge;
      fprintf(file3d, "%lf %lf %lf\n", x[ifr], y[ifr], I_xyz[ifr]);
    
    }
 }
  

  

 /* ------------------------------ */
  /* For disk, now put in z direction  */
  /* ------------------------------ */

 /* evaluate indefinite integral F_z_inf */


  if (strncmp(galaxytype,"disk",4)==0) {

  int iz; 

  double zmin_inf = (-kpc_unit) * (zpoints_inf-1)/2;
  double zmax_inf = kpc_unit * (zpoints_inf-1)/2;

  for (iz=0; iz< zpoints_inf;iz++) {
      z_inf[iz] = zmin_inf + (iz * kpc_unit); 

      //I_z_inf[iz] = I_z0 * pow(2.0/( exp(-fabs(z_inf[iz])/(2*z0_kpc_i)) + exp(fabs(z_inf[iz])/(2*z0_kpc_i)) )  ,2); // I(z) = sech^2(z) 

      I_z_inf[iz] = I_z0/pow( cosh(z_inf[iz]/z0_kpc_i) ,2);

     }


  F_z_inf = cum_integral(z_inf, zmin_inf , zmax_inf, I_z_inf, zpoints_inf); 
 

  /* original grid and input profile (sech(z)^2)*/

   double zmin =  -(kpc_unit) * (zpoints-1)/2;
   double zmax =  kpc_unit * (zpoints+1)/2 ;

  for (iz=0;iz<zpoints;iz++) {
        z[iz] = zmin + (iz * kpc_unit); 

	//I_z[iz] =  I_z0 * pow(2.0/( exp(-fabs(z[iz])/(2*z0_kpc_i)) + exp(fabs(z[iz])/(2*z0_kpc_i)) )  ,2);  // sech^2 profile
	I_z[iz] = I_z0/pow( cosh(z[iz]/z0_kpc_i) ,2)  ;
	fprintf(f_z_file, "%lf %lf\n", z[iz], I_z[iz]);
          
  }


  /* evaluate cumulative pdf F(z) */

  for (iz=0;iz<zpoints;iz++) {

    F_z[iz] = cum_integral(z, zmin, z[iz], I_z, zpoints)/F_z_inf; 
    fprintf(cum_f_z_file, "%lf %lf\n", z[iz], F_z[iz]);
 
  }


  
  /* new grid equally spaced in F(r) from 0 to 1 */
   
  F_interval = 0.0;
  F_interval = (F_z[zpoints-1] - F_z[0])/F_nintervals;
 
  int ifz; 
  F_new = 0.0;

  
  for (ifz=0;ifz<F_nintervals;ifz++) {
    /* choose mid-point of the integral */
    
    F_new = (F_interval * (0.5 + ifz)) ;
     
    for (iz = 0; iz < zpoints; iz++)
	{
 	 
	  if ( (F_new >= F_z[iz]) && (F_new < F_z[iz+1]) ) 
	    {
	        z_store[ifz] = z[iz] ;// store that r value

		//I_z_new[ifz] = I_z0 * pow(2.0/( exp(-fabs(z_store[ifz])/(2*z0_kpc_i)) + exp(fabs(z_store[ifz])/(2*z0_kpc_i)) )  ,2); // I(z) = sech^2(z)  ; 
		I_z_new[ifz] = (I_z0/pow( cosh(z_store[ifz]/z0_kpc_i),2)); 

	        zindex[ifz] = ifz; 

	    } 

	}


  }
  


  /* Randomize the order of z */
  if (F_nintervals > 1) {
        int i;
        for (i = F_nintervals - 1; i > -1; i--) {
            int j = (unsigned int) (drand48()*(i+1));
            int t = zindex[j];
            zindex[j] = zindex[i];
            zindex[i] = t;
	 }
    }


 /* Now map from r(x,y) to z &output 3d galaxy coordinates and intensity */
  
  for (ifz=0; ifz<F_nintervals;ifz++){
    
    I_xyz[ifz] = I_r_new[ifz] * I_z_new[zindex[ifz]] ; 
    

     }



   /*------------------------------------------------- */
   /* Adjust fluxes in two bands and print out results */
   /*------------------------------------------------- */

 
 for (ifz=0; ifz<F_nintervals;ifz++){
    
   I_xyz[ifz] = I_r_new[ifz] * I_z_new[zindex[ifz]] ; // * (1+redshift) * 1e23 * L_solar/(4 * PI * 3.08568025e24 * 3.08568025e24 * d_L * d_L); 
   sumflux_disk +=I_xyz[ifz];
   
     }


  for (ifz=0; ifz<F_nintervals;ifz++)
    {       
      I_xyz[ifz] =  I_xyz[ifz]/sumflux_disk;
      fprintf(file3d, "%lf %lf %lf %lf \n", x[ifz], y[ifz], z_store[zindex[ifz]],  I_xyz[ifz]); 
     
    }


  }  /* end if block for DISK GALAXY only */

 
    gsl_rng_free (random);


    free(npoints_radshell);
    free(npoints_zshell);
    
  return 0;
 
}










double cum_integral(double radius[], double rstart, double rmax, double func[], int npoints)
    {

      double result;


 gsl_integration_workspace * w 
    = gsl_integration_workspace_alloc (npoints);

 {
      gsl_interp_accel *acc 
	= gsl_interp_accel_alloc ();
      gsl_spline *spline 
	= gsl_spline_alloc (gsl_interp_cspline, npoints);

      gsl_spline_init (spline, radius,func, npoints);

      result = gsl_spline_eval_integ (spline,rstart,rmax, acc);

      gsl_spline_free (spline);
      gsl_interp_accel_free (acc);
 
   }
 
gsl_integration_workspace_free (w);


 return result;

    }




