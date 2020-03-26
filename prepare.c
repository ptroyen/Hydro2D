/*
 ============================================================================
 Name        : prepare.c
 Author      : Sagar Pokharel
 Version     :
 Copyright   :
 Description : 
 ============================================================================
 */

 // Standard Libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fenv.h>

// Code Headers
#include "2Dsolver.h"
#include "gridgen.h"
#include "thermo.h"


// UNIVERSAL CONSTANTS
#define ATM 101325.0
#define mili 1.0e-3
#define micro 1.0e-6
#define nano 1.0e-9
#define pico 1.0e-12
#define femto 1.0e-15

// SPECIFIC CONSTANTS
#define GAMMA 1.4
#define c_v 0.718*1000
// #define RHO 1.225 // Air
// #define RHO 1.2506 // Nitrogen
#define RHO 1.1740


// BASIC INITIALIZATION
#define ncx 20
#define ncy 30
#define nx  ncx+1
#define ny  ncy+1
#define nc_row ncy
#define nc_col ncx

#define PI 3.14159265358979323846
#define R0 8.314*1000
#define MW_O2 32.0
#define MW_N2 28.0134
#define MW_M (0.8*MW_N2+0.2*MW_O2)
#define R (R0/MW_M)


//
void initialize(int nc_row , int nc_col, double r[][nc_col] , double u[][nc_col], double v[][nc_col], double p[][nc_col], double tempr[][nc_col]){

    double p0[2] , u0[2] , v0[2] , r0[2] ;

    // Initial Step Conditions

    // SOD Problem
//// Right Going
//    p0[0] = 1.0 ; 		p0[1] = 0.1;
//    u0[0] = 0.0 ;		    u0[1] = 0.0 ;
//    v0[0] = 0.0 ; 		v0[1] = 0.0;
//    r0[0] = 1.0 ; 		r0[1] = 0.125;

//// Left Going
//    p0[0] = 0.1 ; 		p0[1] = 1.0;
//    u0[0] = 0.0 ;		u0[1] = 0.0 ;
//    v0[0] = 0.0 ; 		v0[1] = 0.0;
//    r0[0] = 0.125 ; 	r0[1] = 1.0;


    // Initial Step Conditions For SEDOV - 2D

//    p0[0] = 1*1e-6 ; 	    p0[1] = 1*1e-6;
//    u0[0] = 0.0 ;		    u0[1] = 0.0 ;
//    v0[0] = 0.0 ; 		v0[1] = 0.0;
//    r0[0] = 1.0 ; 		r0[1] = 1.0;

    p0[0] = 1.0 * ATM ; 	p0[1] = 1.0 * ATM;
    u0[0] = 0.0 ;		    u0[1] = 0.0 ;
    v0[0] = 0.0 ; 		    v0[1] = 0.0;
    r0[0] = 1.0 * RHO ;     r0[1] = 1.0 * RHO;


    for (i=0; i < nc_row ; i++){
    		for (j=0; j < nc_col ; j++){

    		// Left State
    			if (j < nc_col/2){
    				p[i][j] = p0[0];
    				u[i][j] = u0[0];
    				v[i][j] = v0[0];
    				r[i][j] = r0[0];
    			}
    		// Right State
    			else{
    				p[i][j] = p0[1];
					u[i][j] = u0[1];
					v[i][j] = v0[1];
					r[i][j] = r0[1];
    			}

                tempr[i][j] = p[i][j] / (r[i][j] * R ) ;

            }
    }
/*
   	// For transforming from left -- right to up -- down
   	    for (j=0; j < nc_col ; j++){
   		    for (i=0; i < nc_row ; i++){
             // Down
       			if (i < nc_col/2){

   				    p[i][j] = p0[0];
   				    u[i][j] = u0[0];
   				    v[i][j] = v0[0];
   				    r[i][j] = r0[0];
   			    }
             // UP
   			    else{
   				    p[i][j] = p0[1];
					    u[i][j] = u0[1];
					    v[i][j] = v0[1];
					    r[i][j] = r0[1];
   			    }
   		    }
   	    }

*/
}
//

//
void get_input(char case[], int *nc_row , int *nc_col, double *XL , double *XR , double *YL, double *YR, double *t0 , double *tend, double *CFL){

    char title[16];
    FILE *fpini;    

    sprintf(title, "%s.inp", case);
    printf("\nReading Input From Case \n");
    fpini = fopen(title, "r");
        fscanf(fpini, "%d", nc_col); fscanf(fpini, "%d", nc_row);        
        fscanf(fpini, "\n%lf", XL); fscanf(fpini, "%lf", XR);
        fscanf(fpini, "\n%lf", YL); fscanf(fpini, "%lf", YR);
        fscanf(fpini, "\n%lf", t0); fscanf(fpini, "%lf", tend);
        fscanf(fpini, "\n%lf", CFL);
     fclose(fpini);
}
//

//
void write_case(char case[], int nc_row , int nc_col, double XL , double XR , double YL, double YR, double t0, double tend,
                double x[] , double y[], double r[][nc_col] , double u[][nc_col], double v[][nc_col], double p[][nc_col], double tempr[][nc_col] ){

    char title[16];
    FILE *fpini;    

    int i, j ;
    double cx[nc_col], cy[nc_row] ;

    // Parameters written
    sprintf(title, "%s.inp", case);
    printf("\nWriting Case \n");
    fpini = fopen(title, "w+");
        fprintf(fpini, "%d", nc_col); fprintf(fpini, "%d", nc_row);        
        fprintf(fpini, "\n%lf", XL); fprintf(fpini, "%lf", XR);
        fprintf(fpini, "\n%lf", YL); fprintf(fpini, "%lf", YR);
        fprintf(fpini, "\n%lf", t0); fprintf(fpini, "%lf", tend);
     fclose(fpini);


     //Write mesh [x,y] not centroids
    sprintf(title, "%s_mesh.inp", case);
    fpini = fopen(title, "w+");
    for (i=0; i < nc_col+1 ; i++){ fprintf(fpini,"%lf \t",x[i]);}
    fprintf(fpini,"\n");
    for (j=0; j < nc_row+1 ; j++){ fprintf(fpini,"%lf \t",y[j]);}
    fclose(fpini);

     // Centroids
    for(i=0; i<nc_col; i++){cx[i] = 0.5*(x[i] + x[i+1]);}
    for(i=0; i<nc_row; i++){cy[i] = 0.5*(y[i] + y[i+1]);}

     // Write the solutions
    sprintf(title, "%s_data.inp", case);
    fpini = fopen(title, "w+");
    for (i=0; i < nc_row ; i++){
        for (j=0; j < nc_col ; j++){
                fprintf(fpini,"%lf		%lf		%lf		%lf		%lf		%lf		%lf \n"
                        ,cx[j],cy[i],r[i][j],u[i][j],v[i][j],p[i][j],tempr[i][j]);
        }
    }
    fclose(fpini);
}
//

//
void save_data(char case[], int nc_row , int nc_col, double time,
                double cx[] , double cy[], double r[][nc_col] , double u[][nc_col], double v[][nc_col], double p[][nc_col], double tempr[][nc_col] ){

    char title[16];
    FILE *fpini;    

     // Write the solutions
    sprintf(title, "%s_data.inp", case);
    fpini = fopen(title, "w+");
    fprintf(fpini,"%lf",time);
    for (i=0; i < nc_row ; i++){
        for (j=0; j < nc_col ; j++){
                fprintf(fpini,"%lf		%lf		%lf		%lf		%lf		%lf		%lf \n"
                        ,cx[j],cy[i],r[i][j],u[i][j],v[i][j],p[i][j],tempr[i][j]);
        }
    }
    fclose(fpini);
}
//


//
void load_case(char case[], int nc_row , int nc_col, double cx[] , double cy[], double r[][nc_col] , double u[][nc_col], double v[][nc_col], double p[][nc_col], double tempr[][nc_col]){

    char title[16];
    FILE *fpini;    

    sprintf(title, "%s_data.inp", case);
    printf("\nReading Case \n");
    fpini = fopen(title, "r");
        for (i=0; i < nc_row ; i++){
            for (j=0; j < nc_col ; j++){
                    fscanf(fpini,"%lf		%lf		%lf		%lf		%lf		%lf		%lf \n"
                            ,cx[j],cy[i],r[i][j],u[i][j],v[i][j],p[i][j],tempr[i][j]);
            }
    }
     fclose(fpini);
}
//


