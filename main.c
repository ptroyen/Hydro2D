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

// Code Headers
#include "2Dsolver.h"


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
#define RHO 1.225


// BASIC INITIALIZATION
#define ncx 300
#define ncy 300
#define nx  ncx+1
#define ny  ncy+1
#define nc_row ncy
#define nc_col ncx

int main(){

    /*
    List of main variables used
    */

    // BASIC INITIALIZATION
    // const int ncx = 400 , ncy = 400 ;
    // const int nx = ncx+1,ny= ncy+1,nc_row = ncy,nc_col = ncx;
    double  l=20 * mili ,lx = l,ly = l,CFL = 0.35;
    static double  t0, time, tend , dt;
    int i, j , k ;

    // Initial Conditions Storage Variables
    // LEFT and RIGHT States
    static double  p0[2], u0[2], v0[2], r0[2] ;

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
    tend = 500.0 * micro ;
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
            dudt[4][nc_row][nc_col],accu_a[4][nc_row][nc_col], source_accu[nc_row][nc_col];
             // add new flux accumulation definition and pass later



    // DEFINE the MESH
    dx = l/(nx-1);    dy = l/(ny-1);
    for(i=0; i<nx; i++){x[i] = i*dx;}
    for(i=0; i<ny; i++){y[i] = i*dy;}
    //  Centroids
    for(i=0; i<ncx; i++){cx[i] = 0.5*(x[i] + x[i+1]);}
    for(i=0; i<ncy; i++){cy[i] = 0.5*(y[i] + y[i+1]);}

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
    concalculate(nc_row, nc_col,r,u,v,p,GAMMA,q);
    concalculate(nc_row, nc_col,r,u,v,p,GAMMA,qfinal);

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


///*******************************************************************
    /// INTEGRATION LOOP STARTS

    while (time < tend) {

        // calculate the timestep
        dt = CFLmaintain(nc_row, nc_col, r,u,v,p,GAMMA,CFL,x,y,count);
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

    // Source Flux Calculate
    // store on source_accu : only for energy equation , see source_accu initialization
        Source(nc_row, nc_col, q, x, y, time,dt, c_v, source_accu );


    // Now Integrate the equation
        for (i=0; i < nc_row ; i++){
            for (j=0; j < nc_col ; j++){
                for (k=0; k < 4 ; k++){
                    dudt[k][i][j] =  - accu_a[k][i][j] ;
                        if (k == 3){
                            dudt[3][i][j] = dudt[3][i][j] + source_accu[i][j] ; // + add diffusive flux
                        }
                    qfinal[k][i][j] = q[k][i][j] + dt * dudt[k][i][j] ;
                }
            }
        }


//// FOR TVD RK, USE THIS PART ELSE FOR EULER COMMENT IT
//// ----------------------------------- start
//
//                    Flux_M(nc_row, nc_col, qfinal,x,y,dt,GAMMA, accu_a);
//                    Source(nc_row, nc_col, qfinal, x, y, time,dt, c_v, source_accu );
//                    for (i=0; i < nc_row ; i++){
//				            for (j=0; j < nc_col ; j++){
//					            for (k=0; k < 4 ; k++){
//                                  dudt[k][i][j] =  - accu_a[k][i][j] ;
//                                      if (k == 3){
//                                           dudt[3][i][j] = dudt[3][i][j] + source_accu[i][j] ; // + add difussive flux
//                                          }
//
//                            // Integration in time here
//                            qfinal[k][i][j] = 0.5 * ( q[k][i][j] +  dt * dudt[k][i][j]) + 0.5 * qfinal[k][i][j] ;
//
//						        }
//				            }
//                    }
//// -------------- END-TVD RK

        // CALCULATE UPDATED VALUE AFTER THE TIME STEP
		prmcalculate(nc_row, nc_col, qfinal,GAMMA,r,u,v,p,speed);


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
           fp = fopen("testout2.txt", "w+");
           //fprintf(fp,"Time =  %0.16f\n", time) ;
            // Centroid_X		Centroid_Y		Density		Velocity_U		Velocity_V		Pressure_P		Energy_E
           for (i=0; i < nc_row ; i++){
                for (j=0; j < nc_col ; j++){
                    fprintf(fp,"%0.7f   %0.7f   %0.10f  %0.10f  %0.10f  %0.10f  %0.10f\n",
                            cx[j],cy[i],r[i][j],u[i][j],v[i][j],p[i][j],q[3][i][j]);
                }
           }
           fclose(fp);
        }
        if ( count == 500){
            fp = fopen("testout40.txt", "w+");
            fprintf(fp,"Time =  %0.16f\n", time) ;
            for (i=0; i < nc_row ; i++){
   				for (j=0; j < nc_col ; j++){
                    fprintf(fp,"%0.7f   %0.7f   %0.10f  %0.10f  %0.10f  %0.10f  %0.10f\n",
                        cx[j],cy[i],r[i][j],u[i][j],v[i][j],p[i][j],q[3][i][j]);
                }
           }
           fclose(fp);
        }




        if ( count % 700 == 0){
            char title[16];
            sprintf(title, "%d.txt", count);
            fp = fopen(title, "w+");
            fprintf(fp,"Time =  %0.16f\n", time);
            for (i=0; i < nc_row ; i++){
   				for (j=0; j < nc_col ; j++){
                    fprintf(fp,"%0.7f   %0.7f   %0.10f  %0.10f  %0.10f  %0.10f  %0.10f\n",
                        cx[j],cy[i],r[i][j],u[i][j],v[i][j],p[i][j],q[3][i][j]);
                }
           }
           fclose(fp);
        }


    }

    /// INTEGRATION LOOP ENDS
///*******************************************************************************

    // FINAL TIMESTEP RESULTS WRITE ON FILE
    fp = fopen("testout.txt", "w+");
    //fprintf(fp,"Time =  %0.16f\n", time) ;
    for (i=0; i < nc_row ; i++){
        for (j=0; j < nc_col ; j++){
            fprintf(fp,"%0.7f   %0.7f   %0.10f  %0.10f  %0.10f  %0.10f  %0.10f\n",
                    cx[j],cy[i],r[i][j],u[i][j],v[i][j],p[i][j],q[3][i][j]);
        }
   }
   fclose(fp);

    return 0;
}


