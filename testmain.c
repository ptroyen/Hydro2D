/*
 ============================================================================
 Name        : main.c
 Author      : Sagar Pokharel
 Version     :
 Copyright   :
 Description : testmain.c
 ============================================================================
 */

 // Standard Libraries
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fenv.h>

// Code Headers
#include "2Dsolver.h"
#include "thermo.h"
#include "gridgen.h"


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
#define ncx 100
#define ncy 100
#define nx  ncx+1
#define ny  ncy+1
#define nc_row ncy
#define nc_col ncx

#define PI 3.14159265358979323846
#define R0 8.314*1000
#define MW_O2 32.0
#define MW_N2 28.0134
#define MW_M (0.80*MW_N2+0.20*MW_O2)
#define R (R0/MW_M)

int main(){

    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

    // printf("R0 = %f , R = %f\n" , R0 , R);


    /*
    List of main variables used
    */

    // BASIC INITIALIZATION
    // const int ncx = 400 , ncy = 400 ;
    // const int nx = ncx+1,ny= ncy+1,nc_row = ncy,nc_col = ncx;
    // NDIM = 0 Cartesian, 1 Cylindrical, 2 Spherical
    int NDIM = 0 ;
    double  l=8 * mili ,lx = 30.0 * mili,ly = 10.0 * mili,CFL = 0.15;
    static double  t0, time, tend , dt;
    int i, j , k ;

    // Initial Conditions Storage Variables
    // LEFT and RIGHT States
    static double  p0[2], u0[2], v0[2], r0[2] ;

    double stored_energy ;
    char title[16];

    static int count = 0 ;
    FILE *fp , *fpini;

    // INITIALIZATION OF VARIABLES  - DONE
        printf("\nVARIABLES INITIALIZED \n");


///***********************************************************

    /// DEFINE BASIC PARAMETERS
    // domain length // 20mm for laser energy deposition
    //l = 20 * mili ;
    //lx = l;        ly = l;

    // Number of Cells
    //ncx = 400;     ncy = 400 ;

    // Time
    t0 = 0.0;
    // tend = 30.0 * micro ; // Laser tear shape
    tend = 600 * micro ; // Laser tear shape
    //CFL NO.
    //CFL = 0.1;

///*************************************************************

        printf("BASIC PARAMETERS SET\n");


    // MESH DEFINITION STARTS

    // ROW and COLUMN Lengths
    //nc_row = ncy;    nc_col = ncx;
    // Number of interfaces
    //nx = ncx+1 ;      ny = ncy+1 ;

    // STORAGE SIZE FIXED
    // MESH PARAMETERS
    static double  dx,dy,x[nx],y[ny],cx[ncx],cy[ncy] ;

    // SOLUTION STORAGE VARIABLES
    static double  p[nc_row][nc_col] , u[nc_row][nc_col] , v[nc_row][nc_col] , r[nc_row][nc_col] ;
    static double  speed[nc_row][nc_col],tempr[nc_row][nc_col], E[nc_row][nc_col], a[nc_row][nc_col] ;
    static double  q[4][nc_row][nc_col],qfinal[4][nc_row][nc_col],
            dudt[4][nc_row][nc_col],accu_a[4][nc_row][nc_col], source_accu[nc_row][nc_col], source_geom[4][nc_row][nc_col];
             // add new flux accumulation definition and pass later


	double exp_ratio = 1.4;


    // DEFINE the MESH
    //-------------------------------------------------------
    // dx = l/(nx-1);    dy = l/(ny-1);
    // for(i=0; i<nx; i++){x[i] = i*dx;}
    // for(i=0; i<ny; i++){y[i] = i*dy;}


    // USING GRIDGENERATOR ------------

   // void grid(int IMAX , int JMAX ,double XL, double XR, double YL, double YU, double x[], double y[]);

   grid(nx,ny,0.0,lx,0.0,ly,x,y);
    //  Centroids
    for(i=0; i<ncx; i++){cx[i] = 0.5*(x[i] + x[i+1]);}
    for(i=0; i<ncy; i++){cy[i] = 0.5*(y[i] + y[i+1]);}


    //---------------------------------------------------------



        printf("MESH DEFINED\n");




     //// FILL IN INITIAL VALUES IN MESH

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


// Laser Energy Deposition Test -- Standard/Normal

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


//    	// For transforming from left -- right to up -- down
//    	    for (j=0; j < nc_col ; j++){
//    		    for (i=0; i < nc_row ; i++){
//              // Down
//        			if (i < nc_col/2){
//
//    				    p[i][j] = p0[0];
//    				    u[i][j] = u0[0];
//    				    v[i][j] = v0[0];
//    				    r[i][j] = r0[0];
//    			    }
//              // UP
//    			    else{
//    				    p[i][j] = p0[1];
//					    u[i][j] = u0[1];
//					    v[i][j] = v0[1];
//					    r[i][j] = r0[1];
//    			    }
//    		    }
//    	    }


    // CALCULATE CONSERVATIVE VARIABLES FOR THE EQUATIONS TO SOLVE :
    // [rho , rho*u , rho*v , E]
    // concalculate(nc_row, nc_col,r,u,v,p,GAMMA,q);
    // concalculate(nc_row, nc_col,r,u,v,p,GAMMA,qfinal);

    // Changes Gamma as well with temperature
    thermo_concalculate(nc_row, nc_col,r,u,v,p,q);
    thermo_concalculate(nc_row, nc_col,r,u,v,p,qfinal);

        printf("CONSERVATIVE VARIABLES INITIALIZED\n");


    	// INITIAL CONDITIONS WRITE TO FILE
    fpini = fopen("initial.txt", "w+");
    for (i=0; i < nc_row ; i++){
        for (j=0; j < nc_col ; j++){
                fprintf(fpini,"%lf		%lf		%lf		%lf		%lf		%lf		%lf \n"
                        ,cx[j],cy[i],r[i][j],u[i][j],v[i][j],p[i][j],q[3][i][j]);
        }
    }
        printf("INITIAL VARIABLES WRITTEN ON initial.txt\n");



    time = t0;

stored_energy = 0.0 ;
///*******************************************************************
    /// INTEGRATION LOOP STARTS

    while (time < tend) {

        // calculate the timestep
        // dt = CFLmaintain(nc_row, nc_col, r,u,v,p,GAMMA,CFL,x,y,count);
        dt = thermo_CFLmaintain(nc_row, nc_col, r,u,v,p,CFL,x,y,count);

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

       // void Geom_F(int nc_row, int nc_col,double q[4][nc_row][nc_col],double gamma,int NDIM, double x[], double y[], double source_geom[4][nc_row][nc_col]){
 
        // Geom_F(nc_row, nc_col,q,GAMMA,NDIM,x,y,source_geom);


    // Now Integrate the equation
        for (i=0; i < nc_row ; i++){
            for (j=0; j < nc_col ; j++){
                for (k=0; k < 4 ; k++){
                    dudt[k][i][j] =  - accu_a[k][i][j] ;//- source_geom[k][i][j];
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

                   Flux_M(nc_row, nc_col, qfinal,x,y,dt,GAMMA, accu_a);
                //    thermo_Flux_M(nc_row, nc_col, qfinal,p,x,y,dt,GAMMA, accu_a);
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
// -------------- END-TVD RK


// Sum of energy added
if ( time < 40*nano){
 for (i=0; i < nc_row ; i++){
                for (j=0; j < nc_col ; j++){
                    for (k=0; k < 4 ; k++){
                            stored_energy = stored_energy+ qfinal[3][i][j] - q[3][i][j] ;
                    }
                }
        }
}

        // CALCULATE UPDATED VALUE AFTER THE TIME STEP
		prmcalculate(nc_row, nc_col, qfinal,GAMMA,r,u,v,p,speed);

        // thermo_prmcalculate(nc_row, nc_col, qfinal,GAMMA,r,u,v,p,tempr,speed);

		// CHECK IF Values are negative HERE


        // UPDATE CONSERVATIVE VARIABLES FOR NEXT INTEGRATION STEP
         for (i=0; i < nc_row ; i++){
                for (j=0; j < nc_col ; j++){
                    for (k=0; k < 4 ; k++){
                            q[k][i][j] = qfinal[k][i][j] ;
                    }
                }
        }


		//temperature(nc_row, nc_col, q,c_v, tempr);


		//  UPDAtE TIME FOR NEXT STEP
		count = count + 1;
		time = time + dt ;

        // DDISPLAY MID CELL'S VALUE IN CONSOLE
        printf("\n TIME : %0.15f",time);
        printf("\t DT : %0.18f \n",dt);
		printf("%lf		%lf		%0.10f		%0.10f		%0.10f		%0.10f		%0.10f \n",
                cx[(int)(nc_col/2)],cy[(int)(nc_row/2)],r[(int)(nc_row/2)][(int)(nc_col/2)],
                u[(int)(nc_row/2)][(int)(nc_col/2)],v[(int)(nc_row/2)][(int)(nc_col/2)],
                p[(int)(nc_row/2)][(int)(nc_col/2)],q[3][(int)(nc_row/2)][(int)(nc_col/2)]);


        // WRITE DATA TO FILE
        if ( count == 1){
           fp = fopen("output/testout2.txt", "w+");
           //fprintf(fp,"Time =  %0.16f\n", time) ;
            // Centroid_X		Centroid_Y		Density		Velocity_U		Velocity_V		Pressure_P		Energy_E
           for (i=0; i < nc_row ; i++){
                for (j=0; j < nc_col ; j++){
                    fprintf(fp,"%0.12f   %0.12f   %0.12f  %0.12f  %0.12f  %0.12f  %0.12f\n",
                            cx[j],cy[i],r[i][j],u[i][j],v[i][j],p[i][j],q[3][i][j]);
                }
           }
           fclose(fp);
        }
        if ( count%10 == 0 && count < 200 ){
            sprintf(title, "output/add%d.txt", count);
            fp = fopen(title, "w+");
            fprintf(fp,"Time =  %0.16f\n", time) ;
            for (i=0; i < nc_row ; i++){
   				for (j=0; j < nc_col ; j++){
                    fprintf(fp,"%0.12f   %0.12f   %0.12f  %0.12f  %0.12f  %0.12f  %0.12f\n",
                        cx[j],cy[i],r[i][j],u[i][j],v[i][j],p[i][j],q[3][i][j]);
                }
           }
           fclose(fp);
        }

        if ( count % 2000 == 0){
            
            sprintf(title, "output/%d.txt", count);
            fp = fopen(title, "w+");
            fprintf(fp,"Time =  %0.16f\n", time);
            for (i=0; i < nc_row ; i++){
   				for (j=0; j < nc_col ; j++){
                    fprintf(fp,"%0.12f   %0.12f   %0.12f  %0.12f  %0.12f  %0.12f  %0.12f\n",
                        cx[j],cy[i],r[i][j],u[i][j],v[i][j],p[i][j],q[3][i][j]);
                }
           }
           fclose(fp);
        }

 
    }

    /// INTEGRATION LOOP ENDS
///*******************************************************************************

    // FINAL TIMESTEP RESULTS WRITE ON FILE
    fp = fopen("output/testout.txt", "w+");
    //fprintf(fp,"Time =  %0.16f\n", time) ;
    for (i=0; i < nc_row ; i++){
        for (j=0; j < nc_col ; j++){
            fprintf(fp,"%0.12f   %0.12f   %0.12f  %0.12f  %0.12f  %0.12f  %0.12f\n",
                    cx[j],cy[i],r[i][j],u[i][j],v[i][j],p[i][j],q[3][i][j]);
        }
   }
   fclose(fp);



   //added energy total
    fp = fopen("output/energy.txt", "w+");
    //fprintf(fp,"Time =  %0.16f\n", time) ;
    fprintf(fp,"Energy added = %0.7f \n",stored_energy);
  
   fclose(fp);

    return 0;
}


