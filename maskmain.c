/*
 ============================================================================
 Name        : main.c
 Author      : Sagar Pokharel
 Version     :
 Copyright   :
 Description : 2D Euler Solver with Energy Addition
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

#define PI 3.14159265358979323846
#define R0 8.314*1000
#define MW_O2 32.0
#define MW_N2 28.0134
#define MW_M (0.8*MW_N2+0.2*MW_O2)
#define R (R0/MW_M)

int main(char case_name[]){

    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

    int i, j , k ;
    int ncx , ncy, nx, ny, nc_row , nc_col ;
    double  XL , XR , YL , YR, CFL;
    double  t0, time, tend , dt;

    // Get Input from file
    get_input(case_name , &ncy , &ncx , &XL , &XR, &YL , &YR, &t0, &tend, &CFL);
    nx = ncx+1; ny= ncy+1 ; nc_row = ncy ; nc_col = ncx;

    // MESH PARAMETERS
    double  dx,dy,x[nx],y[ny],cx[ncx],cy[ncy] ;

    // SOLUTION STORAGE VARIABLES // Free at the end
    double (*r)[nc_col] = malloc(sizeof(double[nc_row][nc_col])),
            (*u)[nc_col] = malloc(sizeof(double[nc_row][nc_col])),
            (*v)[nc_col] = malloc(sizeof(double[nc_row][nc_col])),
            (*p)[nc_col] = malloc(sizeof(double[nc_row][nc_col])),
            (*E)[nc_col] = malloc(sizeof(double[nc_row][nc_col])),
            (*a)[nc_col] = malloc(sizeof(double[nc_row][nc_col])),
            (*speed)[nc_col] = malloc(sizeof(double[nc_row][nc_col])),
            (*tempr)[nc_col] = malloc(sizeof(double[nc_row][nc_col])),
            (*source_accu)[nc_col] = malloc(sizeof(double[nc_row][nc_col]));

    double (*q)[nc_row][nc_col] = malloc(sizeof(double[4][nc_row][nc_col])),
            (*qfinal)[nc_row][nc_col] = malloc(sizeof(double[4][nc_row][nc_col])),
            (*dudt)[nc_row][nc_col] = malloc(sizeof(double[4][nc_row][nc_col])),
            (*accu_a)[nc_row][nc_col] = malloc(sizeof(double[4][nc_row][nc_col])),
            (*source_geom)[nc_row][nc_col] = malloc(sizeof(double[4][nc_row][nc_col]));


    // Load Case Data if Stored
    if (t0 > 0.0){
        load_case(case_name , nc_row, nc_col, x, y, r,u,v,p,tempr);
    }else{
    // Generate Grid and Initialize Data
    grid(nx,ny,XL,XR,YL,YR,x,y);
    initialize( nc_row, nc_col,r,u,v,p,tempr);
    }

    // Centroids
    for(i=0; i<ncx; i++){cx[i] = 0.5*(x[i] + x[i+1]);}
    for(i=0; i<ncy; i++){cy[i] = 0.5*(y[i] + y[i+1]);}

    // NDIM = 0 Cartesian, 1 Cylindrical, 2 Spherical
    int NDIM = 0 ;

    double stored_energy ;
    char title[16];
    static int count = 0 ;
    FILE *fp , *fpini;

    // INITIALIZATION OF VARIABLES  - DONE
    printf("\nVARIABLES INITIALIZED \n");

///********************************************************************
    // [rho , rho*u , rho*v , E]
    concalculate(nc_row, nc_col,r,u,v,p,GAMMA,q);
    concalculate(nc_row, nc_col,r,u,v,p,GAMMA,qfinal);
    printf("CONSERVATIVE VARIABLES INITIALIZED\n");

    // Changes Gamma as well with temperature
    // thermo_concalculate(nc_row, nc_col,r,u,v,p,q);
    // thermo_concalculate(nc_row, nc_col,r,u,v,p,qfinal);

    time = t0;
    stored_energy = 0.0 ;
///*******************************************************************
    /// INTEGRATION LOOP STARTS

    while (time < tend) {

        // calculate the timestep
        dt = CFLmaintain(nc_row, nc_col, r,u,v,p,GAMMA,CFL,x,y,count,time);
        // dt = thermo_CFLmaintain(nc_row, nc_col, r,u,v,p,CFL,x,y,count);

        if (time+dt > tend ){ dt = tend - time ; }

        // Set DUDT to zero before integration
        for (i=0; i < nc_row ; i++){
            for (j=0; j < nc_col ; j++){
                for (k=0; k < 4 ; k++){
                        dudt[k][i][j] = 0.0 ;
                        accu_a[k][i][j] = 0.0 ;
                        //source_accu[k][i][j] = 0.0 ;
                }
            }
        }

    // Convective Flux Calculation del(F)/del(V) = DUDT
    // Store on accu_a
        Flux_M(nc_row, nc_col, q,x,y,dt,GAMMA, accu_a);
        // thermo_Flux_M(nc_row, nc_col, q,p,x,y,dt,GAMMA, accu_a);

    // Source Flux Calculate
    // store on source_accu : only for energy equation , see source_accu initialization
        Source(nc_row, nc_col, q, x, y, time,dt, c_v, source_accu );
        // Geom_F(nc_row, nc_col,q,GAMMA,NDIM,x,y,source_geom);

    // Now Integrate the equation
        for (i=0; i < nc_row ; i++){
            for (j=0; j < nc_col ; j++){
                for (k=0; k < 4 ; k++){
                    dudt[k][i][j] =  - accu_a[k][i][j] - source_geom[k][i][j];
                        if (k == 3){
                            dudt[3][i][j] = dudt[3][i][j] + source_accu[i][j] ; // + add diffusive flux
                        }
                    qfinal[k][i][j] = q[k][i][j] + dt * dudt[k][i][j] ;
                }
            }
        }


// FOR TVD RK, USE THIS PART ELSE FOR EULER COMMENT IT
// ----------------------------------- start
                    // thermo_prmcalculate(nc_row, nc_col, qfinal,GAMMA,r,u,v,p,tempr,speed);
                    // thermo_Flux_M(nc_row, nc_col, qfinal,p,x,y,dt,GAMMA, accu_a);
                    Flux_M(nc_row, nc_col, qfinal,x,y,dt,GAMMA, accu_a);

                   Source(nc_row, nc_col, qfinal, x, y, time,dt, c_v, source_accu );

                   for (i=0; i < nc_row ; i++){
				            for (j=0; j < nc_col ; j++){
					            for (k=0; k < 4 ; k++){
                                 dudt[k][i][j] =  - accu_a[k][i][j] ;
                                     if (k == 3){
                                          dudt[3][i][j] = dudt[3][i][j] + source_accu[i][j] ; // + add difussive flux
                                         }
                           // Integration in time here
                           qfinal[k][i][j] = 0.5 * ( q[k][i][j] +  dt * dudt[k][i][j]) + 0.5 * qfinal[k][i][j] ;

						        }
				            }
                   }
// -------------- END-TVD RK ---------------


// Sum of energy added
if ( time < 40.0*nano){
 for (i=0; i < nc_row ; i++){
                for (j=0; j < nc_col ; j++){
                    for (k=0; k < 4 ; k++){
                            stored_energy = stored_energy+ qfinal[3][i][j] - q[3][i][j] ;
                    }
                }
        }

   //added energy total
    fp = fopen("output/energy.txt", "w+");
    fprintf(fp,"Energy added = %0.7f \n",stored_energy);
   fclose(fp);
}

        // CALCULATE UPDATED VALUE AFTER THE TIME STEP
		prmcalculate(nc_row, nc_col, qfinal,GAMMA,r,u,v,p,speed);
        // thermo_prmcalculate(nc_row, nc_col, qfinal,GAMMA,r,u,v,p,tempr,speed);

        // UPDATE CONSERVATIVE VARIABLES FOR NEXT INTEGRATION STEP
         for (i=0; i < nc_row ; i++){
                for (j=0; j < nc_col ; j++){
                    for (k=0; k < 4 ; k++){
                            q[k][i][j] = qfinal[k][i][j] ;
                    }
                       tempr[i][j] =  (p[i][j]/(r[i][j] * R)) ;
                }
        }

		//  UPDAtE TIME FOR NEXT STEP
		count = count + 1;
		time = time + dt ;

        // DDISPLAY MID CELL'S VALUE IN CONSOLE
        printf("\n TIME : %0.15f",time);
        printf("\t DT : %0.18f \n",dt);
		printf("%lf		%lf		%0.10f		%0.10f		%0.10f		%0.10f		%0.10f \n",
                cx[(int)(nc_col/2)],cy[(int)(nc_row/2)],r[(int)(nc_row/2)][(int)(nc_col/2)],
                u[(int)(nc_row/2)][(int)(nc_col/2)],v[(int)(nc_row/2)][(int)(nc_col/2)],
                p[(int)(nc_row/2)][(int)(nc_col/2)],tempr[(int)(nc_row/2)][(int)(nc_col/2)]);


        if ( count%5 == 0 && time < 45.0e-9 ){
            sprintf(title, "output/add%d.txt", count);
            save_data(title , nc_row, nc_col,time, cx, cy, r,u,v,p,tempr);
        }

        if ( count % 500 == 0){
            sprintf(title, "output/%d.txt", count);
            save_data(title , nc_row, nc_col,time, cx, cy, r,u,v,p,tempr);
        }

        if ( time - tend == 1.0e-37){
            sprintf(title, "%s_out.txt", case_name);
            save_data(title , nc_row, nc_col,time, cx, cy, r,u,v,p,tempr);
        }

    }
 /// INTEGRATION LOOP ENDS
///*******************************************************************************
// Save the simulation as a case so that you can run from this later
    write_case(case_name, nc_row, nc_col, XL,XR,YL,YR,t0,tend,x,y,r,u,v,p,tempr);

    
    //free

    return 0;
}




