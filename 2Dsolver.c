/*
 ============================================================================
 Name        : 2dsolver.c
 Author      : Sagar Pokharel
 Version     :
 Copyright   :
 Description : 2D Euler Solver with Energy Addition
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "2Dsolver.h"

#define PI 3.14159265358979323846

void temperature(int nc_row , int nc_col, double q[][nc_row][nc_col],double c_v, double t[][nc_col]){
/*
    Returns t as the solution for temperature
    Calculate Temperature from the conserved variables
    Use C_v, and total energy per unit volume, density and velocity
*/
	int i,j;
	double tot_energy , rho , vel_x , vel_y;
	for (i=0; i < nc_row ; i++){
		for (j=0; j < nc_col ; j++){
			tot_energy = q[3][i][j];
			rho = q[0][i][j];
			vel_x = q[1][i][j] / rho ;
			vel_y = q[2][i][j] / rho ;
			t[i][j] = (tot_energy - 0.5 * rho * (vel_x*vel_x + vel_y*vel_y))/(rho*c_v) ;

		}
	}
}


void concalculate(int nc_row, int nc_col, double r[][nc_col],double u[][nc_col],double v[][nc_col],double p[][nc_col],double gamma,double q[][nc_row][nc_col]){
    // assembles all properties in q, with conserved variables
	int i,j;

	for (i=0; i < nc_row ; i++){
			for (j=0; j < nc_col ; j++){
				q[0][i][j] = r[i][j];
				q[1][i][j] = u[i][j]*r[i][j];
				q[2][i][j] = v[i][j]*r[i][j];
				q[3][i][j] = (p[i][j]/(gamma - 1)) + 0.5 * r[i][j] * ( u[i][j]*u[i][j] + v[i][j]*v[i][j] );
			}
		}
	}


double superbee(double r,double gl,double gr)
    // Slope Limiter SUPERBEE to use for MUSCL reconstruction
{
    double et = 0.0;

    if (r<=0.0){et = 0;}
    else if (r>=0 && r<= 0.5){et = 2.0*r;}
    else if (r>0.5 && r<= 1.0){et = 1.0;}
    else if (r > 1.0){
             et = fmin(fmin(r,gr),2.0);
        }

    return et;
}


double ultrabee(double r,double gl,double gr){
    // Slope Limiter ULTRABEE to use for MUSCL reconstruction
    double et = 0.0;

    if (r<=0.0){et = 0.0;}
    else if(r>0.0){et = fmin(gl,gr);}
    return et;

}

double vanLeer(double r,double gl,double gr){
    // Slope Limiter VANLEER to use for MUSCL reconstruction
    double et = 0.0;

    if (r<=0.0){et = 0.0;}
    else if(r>0.0){
        et = fmin(((2.0*r)/(1.0+r)),gr);}
    return et;

}


double vanAlbada(double r,double gl,double gr){

    //    Purpose: to compute a VAN ALBADA type slope limiter DELTA
    double et = 0.0;

    if (r<=0.0){et = 0.0;}
    else if(r>0.0){
        et = fmin((r*(1.0+r)/(1+r*r)),gr);
    }
    return et;
}


double CFLmaintain(int nc_row, int nc_col, double r[][nc_col],double u[][nc_col],double v[][nc_col],
                    double p[][nc_col],double gamma,double CFL,double x[nc_col+1],double y[nc_row+1],int n){
//     provides dt required to maintain the provided CFL
//     n is number of time steps

    // double (*a)[y][z] = malloc(sizeof(double[x][y][z])) , (*b)[y][z] = malloc(sizeof(double[x][y][z])),
    //   (*c)[y][z] = malloc(sizeof(double[x][y][z]));
    //

	int i,j;

	// Dynamically Allocate the array : So that large arrays are not stored in stack
	double (*dx) = malloc(sizeof(double[nc_col])) ,
             (*dy) = malloc(sizeof(double[nc_row]));
	// double dx[nc_col], dy[nc_row], dt;
	double dt , spc , mspc , dts , mspcx, mspcy;

	for (i=0; i< nc_col ; i ++){
		dx[i] = x[i+1]-x[i];
	}

	for (i=0; i< nc_row ; i ++){
			dy[i] = y[i+1]-y[i];
		}
	// Largest time limit here
	dt = 1.0E12; // large number

	for (i=0; i < nc_row ; i++){
        for (j=0; j < nc_col ; j++){
            spc = sqrt(gamma*p[i][j]/r[i][j]);
            if ( !(isfinite(spc)) ){spc = 0.0 ;}
            mspc = sqrt(u[i][j]*u[i][j]+v[i][j]*v[i][j]) + spc ;
            // Individual Wave Speed Calculation
			//	mspcx = sqrt(u[i][j]*u[i][j]) + spc ;
			//	mspcy = sqrt(v[i][j]*v[i][j]) + spc ;

            dts = fmin((CFL*dx[j] / mspc),(CFL*dy[i]/ mspc));

            if ((isfinite(dts)) && dts > 0.0 ){
			dt  = fmin( dts , dt);
            }
	    }

        // Time Steps for first Initial Times
        if (dt<=0.0 || !(isfinite(dt)) || dt > 0.1*1.0e-8 || n <= 5){
            printf("TIME STEP CALCULATED = %0.12e  \n" , dt);
			dt = 1e-3 * 1.0e-7;
		}

    }


    // Free the allocated memory for dx and dy
    free(dx);
    free(dy);

    return dt;
}


void prmcalculate(int nc_row, int nc_col,double q[][nc_row][nc_col],double gamma,
                     double r[][nc_col],double u[][nc_col],double v[][nc_col],double p[][nc_col],double s[][nc_col]){

	int i,j;

	for (i=0; i < nc_row ; i++){
        for (j=0; j < nc_col ; j++){

            r[i][j] = q[0][i][j];
            u[i][j] = q[1][i][j]/r[i][j];
            v[i][j] = q[2][i][j]/r[i][j];
            p[i][j] = (gamma - 1)*(q[3][i][j] - 0.5*r[i][j]*(u[i][j]*u[i][j] + v[i][j]*v[i][j]));
            s[i][j] = sqrt(u[i][j]*u[i][j] + v[i][j]*v[i][j]);
		}
	}
}


void sing_prmcalculate(double q[],double gamma, double *r,double *u,double *v,double *p){
    *r = q[0] ;
    *u = q[1]/ (*r) ;
    *v = q[2]/ (*r) ;
    *p = (gamma - 1)*(q[3] - 0.5 * (*r) * ((*u) * (*u) + (*v) * (*v))) ;
}


void x_fun_flux(double qre[],double gamma,double ff[]){
    //    ## Define the function that calculates the flux FOR one direction - X
    //    ## At a time flux of only 1 element in each row is calculated
    //    ## a = wave speed
    //    ## Returns ff

    double rre, ure, vre, pre ;

    rre = qre[0];
    ure = qre[1]/rre;
    vre = qre[2]/rre;
    pre = (gamma - 1)*(qre[3] - 0.5*rre*(ure*ure + vre*vre));


    //    ## q[3] = E :: Totaal Energy Per Unit Volume
    //    ## FLUX FOR EACH EQUATION

    ff[0] = qre[1] ; //## Contuinity Flux
    ff[1] = ff[0]*ure + pre ; //## X-momentum Flux
    ff[2] = ff[0]*vre;
    ff[3] = (ure)*(qre[3]+pre) ;
}


void y_fun_flux(double qre[],double gamma,double ff[]){
    //    ## Define the function that calculates the flux FOR one direction - Y
    //    ## At a time flux of only 1 element in each row is calculated
    //    ## a = wave speed
    //    ## Returns ff

    double rre, ure, vre, pre ;


    rre = qre[0];
    ure = qre[1]/rre;
    vre = qre[2]/rre;
    pre = (gamma - 1)*(qre[3] - 0.5*rre*(ure*ure + vre*vre));

    //    ## q[3] = E :: Totaal Energy Per Unit Volume
    //    ## FLUX FOR EACH EQUATION

    ff[0] = qre[2] ; //## Contuinity Flux
    ff[1] = ff[0]*ure ; //## X-momentum Flux
    ff[2] = ff[0]*vre  + pre  ;
    ff[3] = (vre)*(qre[3]+pre) ;
}


void ESTIME(double DL,double UL,double PL,double DR,double UR,double PR,double GAMMA, double *SL , double *SM , double *SR){
    //#     Purpose: to compute wave speed estimates for the HLLC Riemann
    //#               solver using and adaptive approximate-state Riemann
    //#               solver including the PVRS, TRRS and TSRS solvers.
    //#               See Section 9.5, Chapter 9 of Ref. 1
    //#     Returns SL, SM , SR

	double  G1,G2,G3,G4,G5,G6,G7,QUSER;
	double CL, CR, CUP, PPV ,PMIN,PMAX,QMAX,PM,UM, PQ,PTL,PTR,GEL,GER;


    QUSER = 2.0; //##<<<<<<<<<<<< Select yourself

    CL = sqrt(GAMMA*PL/DL);
    CR = sqrt(GAMMA*PR/DR);

    G1 = (GAMMA - 1.0)/(2.0*GAMMA);
    G2 = (GAMMA + 1.0)/(2.0*GAMMA);
    G3 = (2.0 * GAMMA)/(GAMMA - 1.0);
    G4 = 2.0/(GAMMA - 1.0);
    G5 = 2.0/(GAMMA + 1.0);
    G6 = (GAMMA - 1.0)/(GAMMA + 1.0);
    G7 = (GAMMA - 1.0)/(2.0);

//#      Compute guess pressure from PVRS Riemann solver
    CUP  = 0.25*(DL + DR)*(CL + CR);
    PPV  = 0.5*(PL + PR) + 0.5*(UL - UR)*CUP;
    PPV  = fmax(0.0, PPV);
    PMIN = fmin(PL,  PR);
    PMAX = fmax(PL,  PR);
    QMAX = PMAX/PMIN;

    if (QMAX<=QUSER && (PMIN<=PPV &&  PPV<=PMAX)){

//         Select PVRS Riemann solver
         PM = PPV;
         UM = 0.5*(UL + UR) + 0.5*(PL - PR)/CUP;
    }

    else{
        if(PPV<PMIN){
            // Select Two-Rarefaction Riemann solver
            PQ  = pow((PL/PR),G1);
            UM  = (PQ*UL/CL + UR/CR + G4*(PQ - 1.0))/(PQ/CL + 1.0/CR);
            PTL = 1.0 + G7*(UL - UM)/CL;
            PTR = 1.0 + G7*(UM - UR)/CR;
            PM  = 0.5*(pow(PL*PTL,G3) + pow(PR*PTR,G3));
        }
        else{
            //   Use Two-Shock Riemann solver with PVRS as estimate
            GEL = sqrt((G5/DL)/(G6*PL + PPV));
            GER = sqrt((G5/DR)/(G6*PR + PPV));
            PM  = (GEL*PL + GER*PR - (UR - UL))/(GEL + GER);
            UM  = 0.5*(UL + UR) + 0.5*(GER*(PM - PR) - GEL*(PM - PL));
        }

    }
        // Estimate wave speeds
    if (PM<PL){
        *SL = UL - CL;
    }

    else{
        *SL = UL - CL*sqrt(1.0 + G2*(PM/PL - 1.0));
    }

    *SM = UM;

    if(PM<PR){
            *SR = UR + CR;
    }
    else {
        *SR = UR + CR*sqrt(1.0 + G2*(PM/PR - 1.0));
    }
}


void Flux_M(int nrec_row , int nrec_col , double qre[4][nrec_row][nrec_col] ,double xre[], double yre[] ,double dt,double gamma,  double accu[4][nrec_row][nrec_col]){// ## add what is needed
//    '''
//    -----------------------------------------------------------------------------------------------------
//    >>>>>>>>>>>>  d/dt () + F(net) = 0 <<<<<<<<<<<<<<<<<<<<
//    ## u = vector of the property at centroid points
//    ## cx contains the points in reference to the centroids
//    ## x locates the faces from the cells, so size of x is one more than that of cx
//    ## function 'f = u*(property)' calculates the flux.
//
//    ## Free parameter in real interval [-1,1]
//    ## for w = 0, central difference approximation is achieved
//    ## a is wavespeed
// 	accu[4][nrec_row][nrec_col] RETURNS ON THIS ARRAY
//    ------------------------------------------------------------------------------------------------'''

    //    # initialization
	int i,j,k;

    //	int nrec_row, nrec_col;
	int nre_row,nre_col;

	int nc_row,nc_col,nx,ny;



//	nrec_row = sizeof(qre[0])/sizeof(qre[0][0]);
//	nrec_col = sizeof(qre[0][0])/sizeof(qre[0][0][0]);

	nre_col = nrec_col + 1 ; // Number of x-interfaces.
	nre_row = nrec_row + 1 ; // Number of y-interfaces.



//    Including for ghost cells centroids
    nc_row = nrec_row + 4 ;
    nc_col = nrec_col + 4 ;

// Including ghost cell interfaces
    nx = nre_row + 4 ;
    ny = nre_col + 4 ;


     //// double (*a)[y][z] = malloc(sizeof(double[x][y][z])) , (*b)[y][z] = malloc(sizeof(double[x][y][z])),
    //   (*c)[y][z] = malloc(sizeof(double[x][y][z]));
    //	double (*dx) = malloc(sizeof(double[nc_col])) , (*dy) = malloc(sizeof(double[nc_row]));



    // Real cells
    // flux and accumulation initialization
    double (*Flux)[nre_row][nre_col] = malloc(sizeof(double[4][nre_row][nre_col])),
         (*xflux)[nrec_row][nre_col] = malloc(sizeof(double[4][nrec_row][nre_col])),
         (*yflux)[nre_row][nrec_col] = malloc(sizeof(double[4][nre_row][nrec_col])); // <--- allocate mem

    double (*drex)= malloc(sizeof(double[nrec_col])),
            (*drey)=malloc(sizeof(double[nrec_row])),
            (*vol_re)[nrec_col]=malloc(sizeof(double[nrec_row][nrec_col])) ; // < --- alocate mem

    double (*rre)[nrec_col]=malloc(sizeof(double[nrec_row][nrec_col])),
            (*ure)[nrec_col]=malloc(sizeof(double[nrec_row][nrec_col])),
            (*vre)[nrec_col]=malloc(sizeof(double[nrec_row][nrec_col])),
            (*pre)[nrec_col]=malloc(sizeof(double[nrec_row][nrec_col])); // < --- allocate mem

    // Initialization of values including ghost cells
    double (*q)[nc_row][nc_col] = malloc(sizeof(double[4][nc_row][nc_col])),
            (*x)=malloc(sizeof(double[nx])),
             (*y)=malloc(sizeof(double[ny]));

    double (*qirx)[nc_row][nc_col]=malloc(sizeof(double[4][nc_row][nc_col])),
            (*qilx)[nc_row][nc_col]=malloc(sizeof(double[4][nc_row][nc_col])),
            (*qiry)[nc_row][nc_col]=malloc(sizeof(double[4][nc_row][nc_col])),
            (*qily)[nc_row][nc_col]=malloc(sizeof(double[4][nc_row][nc_col]));

    double (*sx)[nc_row][nc_col] = malloc(sizeof(double[4][nc_row][nc_col])),
        (*sy)[nc_row][nc_col] = malloc(sizeof(double[4][nc_row][nc_col])),
        (*c)[nc_row][nc_col] = malloc(sizeof(double[4][nc_row][nc_col]));

    //## with ghost cells (2 left, 2 right, 2 up, 2 down)
    double (*r)[nc_col] = malloc(sizeof(double[nc_row][nc_col])),
                (*u)[nc_col] = malloc(sizeof(double[nc_row][nc_col])),
                (*v)[nc_col] = malloc(sizeof(double[nc_row][nc_col])),
                (*p)[nc_col] = malloc(sizeof(double[nc_row][nc_col])),
                (*E)[nc_col] = malloc(sizeof(double[nc_row][nc_col])),
                (*a)[nc_col] = malloc(sizeof(double[nc_row][nc_col])) ;

    double (*cflx)[nc_col] = malloc(sizeof(double[nc_row][nc_col])),
            (*cfly)[nc_col] = malloc(sizeof(double[nc_row][nc_col]));


    // For evolution of interface variable, fluxes on left and right
    double xfql[4], xfqr[4],yfql[4],yfqr[4];
    double qlx[4], qrx[4],qly[4],qry[4], HOLD;

    double wx, wy ;

    double upx, lox, upy, loy ;

    double bmx,bpx,bmy,bpy ;

    double ratiox, ratioy , glx, grx , gly , gry;
    double (*spc)[nc_col] = malloc(sizeof(double[nc_row][nc_col])) ,
         (*mspc)[nc_col] = malloc(sizeof(double[nc_row][nc_col])) ;

// s = delta for slope //
// c = wave speed //

//    ## actual cells length / volume / size
    for ( i= 0; i < nrec_col;i++){
    drex[i] = xre[i+1] - xre[i];
    }
    for ( i=0;i < nrec_row; i++){
    drey[i] = yre[i+1] - yre[i];
    }

//## [][] supply row and column number to find the AREA in 2D of the mesh
//vol_re = numpy.matmul(drey[:,None],drex[None,:])
    for ( i = 0; i< nrec_row; i++ ){
    	for (j = 0; j < nrec_col; j++){
    		vol_re[i][j] = drex[j] * drey[i] ;

    		// Primitive Variables Calculation
    		rre[i][j] = qre[0][i][j];
    		ure[i][j] = qre[1][i][j] / rre[i][j];
    		vre[i][j] = qre[2][i][j] / rre[i][j];
    		pre[i][j] = (gamma - 1)*(qre[3][i][j] - 0.5*rre[i][j]*(ure[i][j]*ure[i][j] + vre[i][j]*vre[i][j]));

    	}
    }




//    ## Create Ghost cells 2 left and 2 right, 2 up and 2 down
//    ## Apply BC, later make a function and during this phase create ghost cells and pass to flux calculation
//    ## Flux calculation will return only required values within the real cells
//    ## now apply transimissive on left and reflective on right BC's


    // Left and Right BCs
    for (i=2; i< nc_row - 2 ; i++){

//			# # LEFT TRANSMISSIVE
			u[i][0] = ure[i-2][1];
			u[i][1] = ure[i-2][0];

////			# # LEFT REFLECTIVE
//			u[i][0] = - ure[i-2][1];
//			u[i][1] = - ure[i-2][0];

			v[i][0] = vre[i-2][1];
			v[i][1] = vre[i-2][0];

			r[i][0] = rre[i-2][1];
			r[i][1] = rre[i-2][0];

			p[i][0] = pre[i-2][1];
			p[i][1] = pre[i-2][0];



//			# # RIGHT TRANSMISSIVE
			u[i][nc_col-1] = ure[i-2][nrec_col-2];
			u[i][nc_col-2] = ure[i-2][nrec_col-1];

////			# # RIGHT REFLECTIVE
			// u[i][nc_col-1] = - ure[i-2][nrec_col-2];
			// u[i][nc_col-2] = - ure[i-2][nrec_col-1];

			v[i][nc_col-1] = vre[i-2][nrec_col-2];
			v[i][nc_col-2] = vre[i-2][nrec_col-1];

			r[i][nc_col-1] = rre[i-2][nrec_col-2];
			r[i][nc_col-2] = rre[i-2][nrec_col-1];


			p[i][nc_col-1] = pre[i-2][nrec_col-2];
			p[i][nc_col-2] = pre[i-2][nrec_col-1];


    	}



	// Up and Down BCs
	for (i=2; i< nc_col - 2 ; i++){

//			# # Down TRANSMISSIVE
			v[0][i] = vre[1][i-2];
			v[1][i] = vre[0][i-2];

////			# # down REFLECTIVE
//			v[0][i] = - vre[1][i-2];
//			v[1][i] = - vre[0][i-2];

			u[0][i] = ure[1][i-2];
			u[1][i] = ure[0][i-2];

			r[0][i] = rre[1][i-2];
			r[1][i] = rre[0][i-2];

			p[0][i] = pre[1][i-2];
			p[1][i] = pre[0][i-2];



//			# # Up TRANSMISSIVE
			v[nc_row-1][i] = vre[nrec_row-2][i-2];
			v[nc_row-2][i] = vre[nrec_row-1][i-2];

////			# # UP REFLECTIVE
			// v[nc_row-1][i] = - vre[nrec_row-2][i-2];
			// v[nc_row-2][i] = - vre[nrec_row-1][i-2];

			r[nc_row-1][i] = rre[nrec_row-2][i-2];
			r[nc_row-2][i] = rre[nrec_row-1][i-2];

			u[nc_row-1][i] = ure[nrec_row-2][i-2];
			u[nc_row-2][i] = ure[nrec_row-1][i-2];


			p[nc_row-1][i] = pre[nrec_row-2][i-2];
			p[nc_row-2][i] = pre[nrec_row-1][i-2];


		}



//    ## fill the remaining slots after filling ghost cells

	for (i=2; i< nc_row - 2 ; i++){
		for (j=2; j< nc_col - 2 ; j++){

			u[i][j] = ure[i-2][j-2];
			v[i][j] = vre[i-2][j-2];
			r[i][j] = rre[i-2][j-2];
			p[i][j] = pre[i-2][j-2];

		}
	}


// // Fill the ghost cells that are at the corners as well
// // LEFT BOTTOM
u[0][0] = u[0][3];  u[0][1] = u[0][2];
u[1][0] = u[1][3];  u[1][1] = u[1][2];
v[0][0] = v[0][3];  v[0][1] = v[0][2];
v[1][0] = v[1][3];  v[1][1] = v[1][2];
r[0][0] = r[0][3];  r[0][1] = r[0][2];
r[1][0] = r[1][3];  r[1][1] = r[1][2];
p[0][0] = p[0][3];  p[0][1] = p[0][2];
p[1][0] = p[1][3];  p[1][1] = p[1][2];

// // LEFT TOP
u[nc_col-1][0] = u[nc_col-1][3];        u[nc_col-1][1] = u[nc_col-1][2];
u[nc_col-2][0] = u[nc_col-2][3];    u[nc_col-2][1] = u[nc_col-2][2];
v[nc_col-1][0] = v[nc_col-1][3];        v[nc_col-1][1] = v[nc_col-1][2];
v[nc_col-2][0] = v[nc_col-2][3];    v[nc_col-2][1] = v[nc_col-2][2];
r[nc_col-1][0] = r[nc_col-1][3];        r[nc_col-1][1] = r[nc_col-1][2];
r[nc_col-2][0] = r[nc_col-2][3];    r[nc_col-2][1] = r[nc_col-2][2];
p[nc_col-1][0] = p[nc_col-1][3];        p[nc_col-1][1] = p[nc_col-1][2];
p[nc_col-2][0] = p[nc_col-2][3];    p[nc_col-2][1] = p[nc_col-2][2];

// RIGHT BOTTOM
u[0][nc_col-1] = u[0][nc_col-4];  u[0][nc_col-2] = u[0][nc_col-3];
u[1][nc_col-1] = u[1][nc_col-4];  u[1][nc_col-2] = u[1][nc_col-3];
v[0][nc_col-1] = v[0][nc_col-4];  v[0][nc_col-2] = v[0][nc_col-3];
v[1][nc_col-1] = v[1][nc_col-4];  v[1][nc_col-2] = v[1][nc_col-3];
r[0][nc_col-1] = r[0][nc_col-4];  r[0][nc_col-2] = r[0][nc_col-3];
r[1][nc_col-1] = r[1][nc_col-4];  r[1][nc_col-2] = r[1][nc_col-3];
p[0][nc_col-1] = p[0][nc_col-4];  p[0][nc_col-2] = p[0][nc_col-3];
p[1][nc_col-1] = p[1][nc_col-4];  p[1][nc_col-2] = p[1][nc_col-3];


// RIGHT TOP
u[nc_col-1][nc_col-1] = u[nc_col-1][nc_col-4];  u[nc_col-1][nc_col-2] = u[nc_col-1][nc_col-3];
u[nc_col-2][nc_col-1] = u[nc_col-2][nc_col-4];  u[nc_col-2][nc_col-2] = u[nc_col-2][nc_col-3];
v[nc_col-1][nc_col-1] = v[nc_col-1][nc_col-4];  v[nc_col-1][nc_col-2] = v[nc_col-1][nc_col-3];
v[nc_col-2][nc_col-1] = v[nc_col-2][nc_col-4];  v[nc_col-2][nc_col-2] = v[nc_col-2][nc_col-3];
r[nc_col-1][nc_col-1] = r[nc_col-1][nc_col-4];  r[nc_col-1][nc_col-2] = r[nc_col-1][nc_col-3];
r[nc_col-2][nc_col-1] = r[nc_col-2][nc_col-4];  r[nc_col-2][nc_col-2] = p[nc_col-2][nc_col-3];
p[nc_col-1][nc_col-1] = p[nc_col-1][nc_col-4];  p[nc_col-1][nc_col-2] = p[nc_col-1][nc_col-3];
p[nc_col-2][nc_col-1] = p[nc_col-2][nc_col-4];  p[nc_col-2][nc_col-2] = p[nc_col-2][nc_col-3];

//#   ghor = [ure[:][-2],ure[:][-1]] ## reflective
//#     ghol = [ure[:][1],ure[:][0]] ## transmissive




double dxr, dxl, dyu, dyd ;


    dxr = drex[nrec_col-1];
    dxl = drex[0];

    dyu = drey[nrec_row-1];
    dyd = drey[0];

    x[0] = -dxl*2;
	x[1] = -dxl ;

	// Last = x[nx-1]
	x[nx-1] = xre[nre_col-1]+2*dxr ;
	x[nx-2] = xre[nre_col-1]+dxr ;


	  y[0] = -dyd*2 ;
	  y[1] = -dyd ;

	  	// Last = x[nx-1]
	y[ny-1] = yre[nre_row-1]+2*dyu ;
	y[ny-2] = yre[nre_row-1]+dyu ;

	for ( i=2; i < nx-2 ; i++){
		x[i] = xre[i-2];
	}

	for ( i=2; i < ny-2 ; i++){
			y[i] = yre[i-2];
		}

//    gxl = [-dxl*2, -dxl]
//    gxr = [xre[-1]+dxr, xre[-1]+2*dxr]
//
//    gyd = [-dyd*2, -dyd]
//    gyu = [yre[-1]+dyu, yre[-1]+2*dyu]
//
//    x = numpy.concatenate([gxl, xre, gxr])
//    y = numpy.concatenate([gyd, yre, gyu])

//    ## assemble all props in q, with conserved variables
	for (i=0;i< nc_row; i++){
		for ( j = 0; j < nc_col; j++){
			q[0][i][j] = r[i][j];
			q[1][i][j] = u[i][j]*r[i][j];
			q[2][i][j] = v[i][j]*r[i][j];
			q[3][i][j] = (p[i][j]/(gamma - 1)) + 0.5 * r[i][j]*(u[i][j]*u[i][j] + v[i][j]*v[i][j]);

		}
	}

	double TOL , dxs[nc_col] , centroid_x[nc_col] ,creq_x[nc_col], dys[nc_row] , centroid_y[nc_row] , creq_y[nc_row] ; //# <---- allocate mem


    TOL =  1e-6 ; // 1e-8;

//    ## Size of  each cell's length
    for (i=0 ; i< nc_col ; i++){
    	dxs[i] = x[i+1] - x[i];
    	centroid_x[i] = 0.5 * (x[i+1] + x[i]);
//# CFL number for each cell needed for sign calculation
    	creq_x[i] =  dt / (dxs[i]) ;


    }

    for(i=0 ; i < nc_row ; i++){
    	dys[i] = y[i+1] - y[i];
    	centroid_y[i] = 0.5*(y[i+1] + y[i]);
    	creq_y[i] = dt / (dys[i]);

    }



//    ## Supply vol[][], row and column to get the area of each cell
//    vol = numpy.matmul(dxs[:,None],dys[None,:])



//    cflx = numpy.zeros_like(r)
//    cfly = numpy.zeros_like(r)


//    ## As the velocity of the flow changes need to adapt the time step at each time integration keeping the CFL condition satisified
//    ## Make a new funciton that does this
//
//    ## w = 0; CD type , w = 1; Warming-beam, w = -1; Lax-Wendroff
//    # w = (2*c - np.sign(c))/3 ## Third order accurate in space and time

    

    wx = 0.0;
    wy = 0.0;

    for (i = 1;i< (nc_row-1);i++){
        for (j = 1; j < (nc_col-1); j++){
        	// CHECK THIS IF CASE // IF cases if not executed then previous values are used coz not initialized

        	spc[i][j] = sqrt(gamma*p[i][j]/r[i][j]);
   //               if ( !(isfinite(spc)) || fabs(r[i][j]) < 1.0e-6 ){spc = 0.0 ;} // when pressure and density both decrease rapidly
            mspc[i][j] = sqrt(u[i][j]*u[i][j]+v[i][j]*v[i][j]) + spc[i][j] ;
            if (fabs(u[i][j]) > 1e-4){
                cflx[i][j] = (u[i][j] + spc[i][j]*copysign(1.0,u[i][j]))*creq_x[j];
                cfly[i][j] = (v[i][j]+ spc[i][j]*copysign(1.0,v[i][j]))*creq_y[i];
            	}
            	else {
                cflx[i][j] = 0.0;
                cfly[i][j] = 0.0;

            	}


            	wx = (2.0*cflx[i][j] - copysign(1.0,cflx[i][j])) / 3;
            	wy = (2.0*cfly[i][j] - copysign(1.0,cfly[i][j])) / 3;

            	//if ( fabs(cflx[i][j]) < 1e-5){ wx = 0.0 ; }
            //	if ( fabs(cfly[i][j]) < 1e-5){ wy = 0.0 ; }

            for (k = 0 ; k < 4 ; k++){
	//            # Ratio For Limiters
				upx = q[k][i][j]-q[k][i][j-1];
				lox = q[k][i][j+1]-q[k][i][j];


//            	# Y
				upy = (q[k][i][j]-q[k][i-1][j]);
				loy = (q[k][i+1][j]-q[k][i][j]);



//            	## slope calculation on each cell
//            	## One slope calculation at a time
				sx[k][i][j] = 0.5*(1+wx)*upx+0.5*(1-wx)*lox ;

//            	            # Y
				sy[k][i][j] = 0.5*(1+wy)*upy+0.5*(1-wy)*loy;


				bmx = 1; //# 2/(1+cflx[j]) ##<<< b i-1/2 ## both bm and bp 1 is also suggested
				bpx =  1; //# 2/(1 - cflx[j]) ##<<< b i+1/2 >>>>>> CHEK ITS USAGE

				bmy = 1; //# 2/(1+cfly[i]) ##<<< b i-1/2 ## both bm and bp 1 is also suggested
				bpy =  1; //# 2/(1 - cfly[i]) ##<<< b i+1/2


////## Circumvent division by zero in later stage in calculation of r
//                if (lox[k] == 0.0){
//                     lox[k] = TOL;}
//                if (upx[k] == 0.0){
//                     upx[k] = TOL;}
//
//                if (fabs(upx[k])<TOL){
//                    upx[k] = TOL*copysign(1,upx[k]);}
//                if (fabs(lox[k])<TOL){
//                    lox[k] = TOL*copysign(1,lox[k]); }
//
//                ratiox = upx[k]/lox[k] ;
//
////                # Y
////                ## Circumvent division by zero in later stage in calculation of r
//                if (loy[k] == 0.0){
//					 loy[k] = TOL;}
//				if (upy[k] == 0.0){
//					 upy[k] = TOL;}
//
//				if (fabs(upy[k])<TOL){
//					upy[k] = TOL*copysign(1,upy[k]);}
//				if (fabs(loy[k])<TOL){
//					loy[k] = TOL*copysign(1,loy[k]);}
//
//				ratioy = upy[k]/loy[k] ;
//
//
////            ## Slope Limiters Usagae


				//## Circumvent division by zero in later stage in calculation of r
//				                if (lox == 0.0){
//				                     lox = TOL;}
//				                if (upx == 0.0){
//				                     upx = TOL;}

				                if (fabs(upx)<TOL){
				                    upx = TOL*copysign(1.0,upx);}
				                if (fabs(lox)<TOL){
				                    lox = TOL*copysign(1.0,lox); }

				                ratiox = upx/lox ;

				         //       if (!(isfinite(ratiox))){ ratiox = 1.0;}

				//                # Y
				//                ## Circumvent division by zero in later stage in calculation of r
//				                if (loy == 0.0){
//									 loy = TOL;}
//								if (upy == 0.0){
//									 upy = TOL;}

								if (fabs(upy)<TOL){
									upy = TOL*copysign(1.0,upy);}
								if (fabs(loy)<TOL){
									loy = TOL*copysign(1.0,loy);}

								ratioy = upy/loy ;

							//	if (!(isfinite(ratioy))){ ratioy = 1.0;}


				//            ## Slope Limiters Usagae
   // #----------------------------------------------------------------------------------


                glx = (2*bmx*ratiox)/(1 - wx + (1+wx)*ratiox);
                grx = (2*bmx)/(1 - wx + (1+wx)*ratiox);
                sx[k][i][j] = sx[k][i][j] * vanAlbada(ratiox,glx,grx); //## Ref: Toro's Book 13.218 -- slope limiters ## superbee ## ultrabee ## This is the  delta


                gly = (2*bmy*ratioy)/(1 - wy + (1+wy)*ratioy);
                gry = (2*bmy)/(1 - wy + (1+wy)*ratioy);
                sy[k][i][j] = sy[k][i][j] * vanAlbada(ratioy,gly,gry);


//## now calculate ul and ur
           qlx[k] = q[k][i][j] - 0.5*sx[k][i][j]; // ????
           qrx[k] = q[k][i][j] + 0.5*sx[k][i][j];

//           ## now calculate ul and ur
           qly[k] = q[k][i][j] - 0.5*sy[k][i][j];
           qry[k] = q[k][i][j] + 0.5*sy[k][i][j];



         }

//    #----------------------------------------------------------------------------------
//    ### up to here considered doing single equations at a time

            x_fun_flux(qlx, gamma, xfql);
            x_fun_flux(qrx, gamma, xfqr);

            y_fun_flux(qly, gamma, yfql);
            y_fun_flux(qry, gamma, yfqr);



            for ( k = 0; k < 4 ; k++){
            	// ## caclulate evolution of ulx and ury in dt/2
            	 HOLD = 0.5 * (dt/(x[j+1]-x[j])) * (xfql[k] - xfqr[k]) + 0.5 * (dt/(y[i+1]-y[i])) * (yfql[k] - yfqr[k]) ;
         //   HOLD = 0.5 * (dt/(x[j+1]-x[j])) * (xfql[k] - xfqr[k])  ;

            	qilx[k][i][j] = qlx[k] + HOLD ;
				qirx[k][i][j] = qrx[k] + HOLD ;


	    //	HOLD	= 0.5 * (dt/(y[i+1]-y[i])) * (yfql[k] - yfqr[k]);

//				## calculate evolution of uly and ury in dt/2
				qily[k][i][j] = qly[k] + HOLD ;
				qiry[k][i][j] = qry[k] + HOLD ;

            }

        }

    }




//        #### LOOK BELOW THIS ................>>>>>>>
//
//    ## Now using these evolved fluxes at left and right
//    ## compute the accumulated flux in each cell
//    ## divide by vol and is ready to integrate in time
//
//
//    ## Call HLLC flux calculation scheme and accumulated flux in each cell will be in fmhm
//    ## supplied qil and qir contain all the conserved variable values within the computation domain

// void HLLC(int nc_row , int nc_col , double Qil[4][nc_row][nc_col],double Qir[4][nc_row][nc_col],double GAMMA, char DIR ,int f_row , int f_col, double Flux[4][f_row][f_col]);



/// USE HLLC
//    HLLC(nc_row, nc_col, qilx,qirx,gamma,'x', nrec_row , nre_col,xflux);
//    HLLC(nc_row, nc_col ,qily,qiry,gamma,'y', nre_row, nrec_col,yflux);



 /// USE LAX-Ferdich
//    LFflux(nc_row, nc_col, qilx,qirx,gamma,'x', mspc, nrec_row , nre_col,xflux);
//    LFflux(nc_row, nc_col ,qily,qiry,gamma,'y',mspc, nre_row, nrec_col,yflux);

    /// Rusonov

    RSflux(nc_row, nc_col, qilx,qirx,gamma,'x', nrec_row , nre_col,xflux);
    RSflux(nc_row, nc_col ,qily,qiry,gamma,'y', nre_row, nrec_col,yflux);




//    ## accumulation of flux from both the sides

    for (i = 0 ; i < nrec_row ; i++){
    	for ( j = 0 ; j < nrec_col ; j++){
    		for ( k =0; k<4; k++){
    			accu[k][i][j] = (-xflux[k][i][j] + xflux[k][i][j+1])/drex[j] +  (-yflux[k][i][j] + yflux[k][i+1][j])/drey[i] ;
    		}
    	}
    }
//    accu = (-xflux[:,:,:-1] + xflux[:,:,1:])/drex[None,None,:] + (- yflux[:,:-1,:] + yflux[:,1:,:])/drey[None,:,None]



//    ## direct flux contribution in time rate of change of the property, after this integrate in time
//    dudt = accu
//    return dudt

   // double (*Flux)[nre_row][nre_col] = malloc(sizeof(double[4][nre_row][nre_col])), (*xflux)[nrec_row][nre_col] = malloc(sizeof(double[4][nrec_row][nre_col])),(*yflux)[nre_row][nrec_col] = malloc(sizeof(double[4][nre_row][nrec_col])); // <--- allocate mem
    // double (*rre)[nrec_col]=malloc(sizeof(double[nrec_row][nrec_col])),
//            (*ure)[nrec_col]=malloc(sizeof(double[nrec_row][nrec_col])),
//            (*vre)[nrec_col]=malloc(sizeof(double[nrec_row][nrec_col])),
//            (*pre)[nrec_col]=malloc(sizeof(double[nrec_row][nrec_col]));


    free(Flux); free(xflux);  free(yflux);
    free(drex); free(drey); free(vol_re);
    free(ure); free(vre); free(pre); free(rre);
    free(q); free(x); free(y);
    free(qirx); free(qilx); free(qily); free(qiry);
    free(sx); free(sy); free(c);
    free(r); free(u); free(v); free(p); free(E); free(a);
    free(cflx); free(cfly);
    free(spc); free(mspc);

}



void HLLC(int nc_row , int nc_col , double Qil[4][nc_row][nc_col],double Qir[4][nc_row][nc_col],double GAMMA, char DIR ,int f_row , int f_col, double Flux[4][f_row][f_col]){
//# *
//# C     Purpose: to compute the intercell flux FI(K, I) according
//# C              to the Godunov method with the HLLC approximate
//# C              Riemann solver. See Chap. 10 of E TORO's BOOK

	// The passed Qil and Qir are the conserved variables which contain 2 ghost cell left and 2 ghost cell right


	int  inter_row , inter_col , real_cell_row , real_cell_col , real_inter_row , real_inter_col ;//nc_row , nc_col ;

	int i, j, k, l ,m , n ;

	// number of cells after adding ghost cells
//	nc_row = sizeof(Qil[0])/sizeof(Qil[0][0]);
//	nc_col = sizeof(Qil[0][0])/sizeof(Qil[0][0][0]);

		 // Number of x-interfaces.
		 // Number of y-interfaces.


    inter_row = nc_row + 1;
    inter_col = nc_col + 1;
    real_cell_row = nc_row - 4;
    real_cell_col = nc_col - 4;
    real_inter_row = nc_row - 4 + 1 ;//## number of interfaces in real = total including ghost cells - subtract ghosts + add one
    real_inter_col = nc_col - 4 + 1 ;



    double CSR[4] = {0.0};
    double CSL[4] = {0.0};


//    double Flux[4][real_cell_row][real_inter_col] // Define this while passing from parent function

//            double qil[4][real_cell_row][real_inter_col] ;
//            double qir[4][real_cell_row][real_inter_col] ;
//
//
//            double fil[4][real_cell_row][real_inter_col] ;//= {0.0};
//            double fir[4][real_cell_row][real_inter_col] ;//= {0.0};
            double rL,uL,vL,pL,rR,uR,vR,pR;

//            double RL[],UL[],VL[],PL[],SPL[],RR[],UR[],VR[],PR[], SPR[];

            double rl,ul,vl,pl,rr,ur,vr,pr, ENEL, ENER;

            double SL , SM , SR ;

            double T_QL[4], T_QR[4] ;

    if (DIR == 'x'){

    	double (*qil)[real_cell_row][real_inter_col] = malloc(sizeof(double[4][real_cell_row][real_inter_col])) ;
		double (*qir)[real_cell_row][real_inter_col] = malloc(sizeof(double[4][real_cell_row][real_inter_col])) ;


		double (*fil)[real_cell_row][real_inter_col] = malloc(sizeof(double[4][real_cell_row][real_inter_col])) ;//= {0.0};
		double (*fir)[real_cell_row][real_inter_col] = malloc(sizeof(double[4][real_cell_row][real_inter_col])) ;//= {0.0};
		double (*RL)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col])),
                (*UL)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col])),
                (*VL)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col])),
                (*PL)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col])),
                (*SPL)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col])),
                (*RR)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col])),
                (*UR)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col])),
                (*VR)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col])),
                (*PR)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col])), 
                (*SPR)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col]));



       //  Flux calculation requires left and right values on the interface of diffrent but adjecent to the interface cells
        for (i=2; i< 2 + real_cell_row; i++){ //## change this to make consistent with y as well
            for (j = 2 ; j < inter_col-2; j++ ) {

                //printf(i,j);
                m = i-2 ;
                n = j-2 ;
                for ( k = 0 ; k < 4 ; k++){
                qil[k][m][n] = Qir[k][i][j-1] ; // for left from interface is right of upstream cell
                qir[k][m][n] = Qil[k][i][j] ;// for right from interface is left for downstream cell


                T_QL[k] = qil[k][m][n];
                 T_QR[k] = qir[k][m][n];





				}

				sing_prmcalculate(T_QL,GAMMA,&rL,&uL,&vL,&pL);
				sing_prmcalculate(T_QR,GAMMA,&rR,&uR,&vR,&pR );


                fil[0][m][n] = rL*uL;
				fil[1][m][n] = rL*uL*uL+pL ;
				fil[2][m][n] = rL*uL*vL ;
				fil[3][m][n] = uL*((qil[3][m][n])+pL) ;


				fir[0][m][n] = rR*uR;
				fir[1][m][n] = rR*uR*uR+pR ;
				fir[2][m][n] = rR*uR*vR ;
				fir[3][m][n] = uR*(qir[3][m][n]+pR) ;


                }
        }

        // prmcalculate(int nc_row, int nc_col,		double q[][nc_row][nc_col],double gamma,
        // double r[][nc_col],double u[][nc_col],double v[][nc_col],double p[][nc_col],double s[][nc_col])

        prmcalculate(real_cell_row , real_inter_col, qil,GAMMA , RL,UL,VL,PL,SPL);
        prmcalculate(real_cell_row , real_inter_col, qir,GAMMA, RR,UR,VR,PR,SPR );


        // ---------

        for (i=0; i<real_cell_row; i++){
                      for (j = 0; j < real_inter_col; j++){

                          m = i ;
                          n = j ;


                          rl = RL[i][j] ;
                          ul = UL[i][j] ;
                          vl = VL[i][j] ;
                          pl = PL[i][j] ;

                          rr = RR[i][j] ;
                          ur = UR[i][j] ;
                          vr = VR[i][j] ;
                          pr = PR[i][j] ;


                          ESTIME(rl,ul,pl,rr,ur,pr,GAMMA , &SL, &SM, &SR) ;




                          if (SL>=0.0){
                 // #          Right-going supersonic flow
                          	for ( k = 0 ; k< 4 ; k++){
                              Flux[k][m][n] = fil[k][i][j] ;
                          	}

                          }

                          if (SL<0.0 && SR>0.0){
                            //  ## SUBSONIC FLOW

                              if (SM>=0.0){
                                 // ## SUBSONIC FLOW to the RIGHT##
        //##  qsL = rL*(SL-vnL)/(SL-SM)*[1; SM*nx+uL*abs(ny); SM*ny+vL*abs(nx); qL(4)/rL + (SM-vnL)*(SM+pL/(rL*(SL-vnL)))];
                                  ENEL   = qil[3][i][j]/rl  + (SM - ul)*(SM + pl/(rl*(SL - ul))) ;
                                  CSL[0] = rl*(SL - ul)/(SL - SM) ;
                                  CSL[1] = CSL[0]*SM ;
                                  CSL[2] = CSL[0]*vl ;
                                  CSL[3] = CSL[0]*ENEL ;

                                  for ( k = 0 ; k< 4 ; k++){
                                                          Flux[k][m][n] = fil[k][i][j] + SL*(CSL[k] - qil[k][i][j]);
                                                      	}
                              }

                              else {
                        //  ## Subsonic flow to the left
                                  ENER   = qir[3][i][j]/rr  + (SM - ur)*(SM + pr/(rr*(SR - ur))) ;
                                  CSR[0] = rr*(SR - ur)/(SR - SM) ;
                                  CSR[1] = CSR[0]*SM ;
                                  CSR[2] = CSR[0]*vr ;
                                  CSR[3] = CSR[0]*ENER ;


      							for ( k = 0 ; k< 4 ; k++){
      							          Flux[k][m][n] = fir[k][i][j] + SR*(CSR[k] - qir[k][i][j]);
      							           	}
                              }
                              }
                          if (SR<=0.0){
                                //  ## Left-going supersonic flow
                          	for ( k = 0 ; k< 4 ; k++){
                          	           Flux[k][m][n] = fir[k][i][j] ;
                          	          }
                          			}


                      }

                  }

        // --------

    free(qil); free(qir);
    free(fil); free(fir);
    free(RL); free(UL); free(VL) ; free(PL) ; free(SPL);
    free(RR); free(UR); free(VR) ; free(PR) ; free(SPR);

    }



    else{

    	double (*qil)[real_inter_row][real_cell_col]= malloc(sizeof(double[4][real_inter_row][real_cell_col])) ;
		double (*qir)[real_inter_row][real_cell_col]= malloc(sizeof(double[4][real_inter_row][real_cell_col])) ;


		double (*fil)[real_inter_row][real_cell_col]= malloc(sizeof(double[4][real_inter_row][real_cell_col])) ;//= {0.0};
		double (*fir)[real_inter_row][real_cell_col]= malloc(sizeof(double[4][real_inter_row][real_cell_col])) ;//= {0.0};
		double (*RL)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col])),
                (*UL)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col])),
                (*VL)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col])),
                (*PL)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col])),
                (*SPL)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col])),
                (*RR)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col])),
                (*UR)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col])),
                (*VR)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col])),
                (*PR)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col])), 
                (*SPR)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col]));


    	            // Flux calculation requires left and right values on the interface of diffrent but adjecent to the interface cells


    	// -----------------------------
    // CHECK HERE	double Flux[4][real_cell_row][real_inter_col] // Define this while passing from parent function

//    	        double qil[4][real_cell_row][real_inter_col] = {0.0};
//    	        double qir[4][real_cell_row][real_inter_col] = {0.0};
//
//
//    	        double fil[4][real_cell_row][real_inter_col] = {0.0};
//    	        double fir[4][real_cell_row][real_inter_col] = {0.0};
//
//    	        double rL,uL,vL,pL,rR,uR,vR,pR;
//
//    	        double RL,UL,VL,PL,RR,UR,VR,PR;



    	       //  Flux calculation requires left and right values on the interface of diffrent but adjecent to the interface cells
    	        for (i=2; i< inter_row - 2; i++){ //## change this to make consistent with y as well
    	            for (j = 2 ; j < nc_col-2; j++ ) {

    	                //printf(i,j);
    	                m = i-2 ;
    	                n = j-2 ;
    	                for ( k = 0 ; k < 4 ; k++){
    	                qil[k][m][n] = Qir[k][i-1][j] ; // for left from interface is right of upstream cell
    	                qir[k][m][n] = Qil[k][i][j] ;// for right from interface is left for downstream cell

    	                T_QL[k] = qil[k][m][n];
    	                T_QR[k] = qir[k][m][n];





    	                }

    	                sing_prmcalculate(T_QL,GAMMA,&rL,&uL,&vL,&pL);
    	                sing_prmcalculate(T_QR,GAMMA,&rR,&uR,&vR,&pR );


    	                fil[0][m][n] = rL*vL;
    					fil[1][m][n] = rL*vL*uL ;
    					fil[2][m][n] = rL*vL*vL +pL ;
    					fil[3][m][n] = vL*(qil[3][m][n]+pL) ;


    					fir[0][m][n] = rR*vR;
    					fir[1][m][n] = rR*vR*uR ;
    					fir[2][m][n] = rR*vR*vR +pR ;
    					fir[3][m][n] = vR*(qir[3][m][n]+pR) ;


    	                }
    	        }
    	        prmcalculate(real_inter_row,real_cell_col,qil,GAMMA , RL,UL,VL,PL,SPL);
    	        prmcalculate(real_inter_row,real_cell_col,qir,GAMMA, RR,UR,VR,PR,SPR );


    	        // ------------

                for (i=0 ; i< real_inter_row; i++){
                        for (j = 0; j < real_cell_col; j++){

                                m = i;
                                n = j;


                                rl = RL[i][j];
                               // # ul = uL[i][j];
                               // # vl = vL[i][j];
                                pl = PL[i][j];

                                rr = RR[i][j];
                               // # ur = uR[i][j];
                               // # vr = vR[i][j];
                                pr = PR[i][j];

                                ul = VL[i][j];
                                vl = UL[i][j];

                                ur = VR[i][j];
                                vr = UR[i][j];

                                ESTIME(rl,ul,pl,rr,ur,pr,GAMMA,&SL, &SM, &SR);


                                if (SL>=0.0){
                             //   #                Right-going supersonic flow
                                	for ( k = 0 ; k< 4 ; k++){
                                			Flux[k][m][n] = fil[k][i][j];
                                	}
                                }
                                if (SL<0.0 && SR>0.0){
                                   // ## SUBSONIC FLOW

                                    if (SM>=0.0 && SL<0.0){
                                      //  ## SUBSONIC FLOW to the RIGHT##
                                        ENEL   = qil[3][i][j]/rl  + (SM - ul)*(SM + pl/(rl*(SL - ul))) ;
                                        CSL[0] = rl*(SL - ul)/(SL - SM) ;
                                        CSL[1] = CSL[0]*vl ;
                                        CSL[2] = CSL[0]*SM ;
                                        CSL[3] = CSL[0]*ENEL ;

                                        for ( k = 0 ; k< 4 ; k++){
                                        				Flux[k][m][n] = fil[k][i][j] + SL*(CSL[k] - qil[k][i][j]);
                                        							           	}

                                    }
                                    else{
                                         //   ## Subsonic flow to the left
                                        ENER   = qir[3][i][j]/rr  + (SM - ur)*(SM + pr/(rr*(SR - ur))) ;
                                        CSR[0] = rr*(SR - ur)/(SR - SM) ;
                                        CSR[1] = CSR[0]*vr ;
                                        CSR[2] = CSR[0]*SM ;
                                        CSR[3] = CSR[0]*ENER ;
                                        for ( k = 0 ; k< 4 ; k++){
                                        Flux[k][m][n] = fir[k][i][j] + SR*(CSR[k] - qir[k][i][j]);
                                     	}

                                    }


                                }

                                if(SR<=0.0){
                                       // ## Left-going supersonic flow

                                	for ( k = 0 ; k< 4 ; k++){
                                    Flux[k][m][n] = fir[k][i][j] ;
                                								           	}


                                    }

                        }

                }

    	        // ------------
                free(qil); free(qir);
                free(fil); free(fir);
                free(RL); free(UL); free(VL) ; free(PL) ; free(SPL);
                free(RR); free(UR); free(VR) ; free(PR) ; free(SPR);
                        
            
            }



//    ## Start Calculating the Fluxes
//    ## Change the looping and check the flux values used.. i vs n

//    if (DIR == 'x'){
//
//
//    }
//    else{
//
//
//
//    }


   
    

}



// Make a Structure to store Parameters for a Gaussian Profile
// Store Amplitude, x-center, y-center, sigmax, and sigmay

struct par_gauss
{
    double A, x0, y0, sigma_x , sigma_y ;

};

struct par_expo
{


};


// SOURCE FLUX
//void HLLC(int nc_row , int nc_col , double Qil[4][nc_row][nc_col],double Qir[4][nc_row][nc_col],double GAMMA, char DIR ,int f_row , int f_col, double Flux[4][f_row][f_col]);


void Source(int nc_row, int nc_col,double q[4][nc_row][nc_col], double x[], double y[], double t,double dt, double c_v, double source_accu[nc_row][nc_col]){
 /*   '''
---------------------------------------------------------------------------------------
 >>>>>>>>>>>           ### d/dt()dv = S dv ###          <<<<<<<<<<<<<<<
## Provide the source term contribution as source per unit time in every cell
 -----------------------------------------------------------------------------------------
 '''   */

//## value initialization for temperature
//    # c_v = 0.718 * 1000 #J/kg-K
//    # K = 22.22e-3  # W/m-K
//
//    # l = len(vec)
//
// # initialization
 int n_col, n_row, nx , ny ;

//    nc_row,nc_col = numpy.shape(q[0,:,:])
    n_col = nc_col + 1 ;
    n_row = nc_row + 1 ;


    nx = n_col;
    ny = n_row;

//    # Flux = numpy.zeros([n_row,n_col])

    int i, j ,k ;
    //double xflux[nc_row][n_col] , yflux[n_row][nc_col] ;
    double (*drex)=malloc(sizeof(double[nc_col])), (*drey)=malloc(sizeof(double[nc_row])) ,
     (*centroid_x)=malloc(sizeof(double[nc_col])) , (*centroid_y)=malloc(sizeof(double[nc_row])) , vol;
	 
	 double che , midx , midy , fac_tear;


    // LASER PROPERTIES DEFINITIONS
    double w0, lam , Rl_rn, f ;


    // The Dual Pulse Laser :: SECOND PULSE
    w0 = 200.0e-6 ; lam = 1064.0e-9 ; f = 300.0e-3 ;
    Rl_rn = 0.1* PI*(w0*w0) / lam ;






    // Use structure to define the gauss function paratmeters
    struct par_gauss gaus1 , gaus2 , gaus3 , gaus4 ;
    // double A , sigmax , sigmay ; 
//    source_accu = numpy.zeros([nc_row,nc_col])





     //## Centroids

         for (i=0 ; i< nc_col ; i++){
            drex[i] = (x[i+1] - x[i]);
            centroid_x[i] = 0.5 * (x[i+1] + x[i]);
//# CFL number for each cell needed for sign calculation

    }

    for(i=0 ; i < nc_row ; i++){
    	drey[i] = (y[i+1] - y[i]);
    	centroid_y[i] = 0.5*(y[i+1] + y[i]);


    }

    for ( i = 0 ; i < nc_row ; i++){
        for( j = 0 ; j< nc_col ; j++){

            source_accu[i][j] = 0.0 ;
        }
    }



/// Tear Drop Shape Energy Addition and Evolution in Time Study Case





   // ## MESH Definition Matrix For X's and Y's

   // mesh_x, mesh_y = numpy.meshgrid(centroid_x,centroid_y)

//    ## actual cells length / volume / size
//    drex = x[1:] - x[:-1]
//    drey = y[1:] - y[:-1]
//    ## [][] supply row and column number to find the AREA in 2D of the mesh
//    vol_re = numpy.matmul(drey[:,None],drex[None,:])





//    A = 0.3
//    k = -1
//    R = 0.2


   // ## The Radius Matrix
//    radius = numpy.sqrt(mesh_x**2 + mesh_y**2)

 //   ## Source present for cells below the given radius
//    has_source = numpy.less(radius,R)
//    ## source, i.e energy per unit time per unit volume
//
//    # ## Instant source at t=0
//    # if t==0.0:
//    #     S = 0.3 * has_source
//    # else:
//    #     S = 0.0
//
//    # ## Generalized Source
//    # S = A * t**k * has_source
//
//    # source_accu = S * v
//
//
//    ## Etot = 0.244816 FOR SEDOV Blast wave

    /*'''
    Ref: ValidationTestCaseSuiteforcompressiblehydrodynamicscomputationRaphael LoubereTheoreticalDivision,T-7,Los AlamosNationalLaboratory
    '''*/
//    ## instant source
//    ## first point/cell only has a source constribution
//    ## supplied value is dudt, after this integrate in time only

	fac_tear = 20000.0; // bigger means bigger width
	midx = centroid_x[(int)(nc_col/2)];
	midy = centroid_y[(int)(nc_row/2)];


// At Center
        vol = drex[(int)(nc_col/2)] * drey[(int)(nc_row/2)]; // <-- this should change with variable divisions in grid
        printf("\nVOL = %lf\n", vol);
    

    if (t<1*5.0e-9){

        for ( i = 0 ; i < nc_row ; i++){
            for( j = 0 ; j< nc_col ; j++){
        // source_accu[0,0] = 0.851072 / dt ## 0.244816/dt

        // AT CENTER --------------------------------------------------
        
        // source_accu[int(nc_row/2),int(nc_col/2)] = 1.0*0.244816/(dt*vol) ## 0.851072 / dt ## 0.244816/dt
       // source_accu[(int)(nc_row/2)][(int)(nc_col/2)] = 1.0*0.311357/(dt*vol) ; // 0.851072 / dt ## 0.244816/dt


        /// LASER ENRGY DEPOSITION
        /* 1.0E-12 J/s ; For 0.01 nanoseconds energy deposited = 10 J ; */
        // source_accu[(int)(nc_row/2)][(int)(nc_col/2)] = 50* 1.0e-3 /(dt*vol) ; // 0.851072 / dt ## 0.244816/dt

        // printf("SOURCE = %f\n", source_accu[(int)(nc_row/2)][(int)(nc_col/2)] );



//        # At Origin
//        # source_accu[0,0] = 1.0*0.0244816/(dt*vol) ## 0.851072 / dt ## 0.244816/dt
//
//        # source_accu[0,0] = 0.0*3.11357e-8/(dt*vol) ## 0.851072 / dt ## 0.244816/dt
//
//        # source_accu[0,0] = 0.25*0.311357/(dt*vol) ## 0.851072 / dt ## 0.244816/dt

/// AT CENTER END -----------------------------------------------------------

/// TEAR DROP WITH MIDDLE CIRCULAR SHAPE -------------------------------------

// // Tear Shape
//             che = pow((centroid_x[j] - midx),2) - fac_tear*pow((1.75*1.0e-3-centroid_y[i]+midy),3)*(1.75*1.0e-3+centroid_y[i]-midy);
//             vol = drex[j] * drey[i];
//             if (che <= 0.0){source_accu[i][j] = 8.0* 1.0e-3 /(dt*vol) ;}



//             fac_tear = 20000.0;
//             //more intense one as well one also tear shaped
//             che = pow((centroid_x[j] - midx),2) - fac_tear*pow((1.2*1.0e-3-centroid_y[i]+midy-0.25*1.0e-3),3)*(1.2*1.0e-3+centroid_y[i]-midy+0.25*1.0e-3);
//             if (che <= 0.0){source_accu[i][j] =source_accu[i][j] + 16* 1.0e-3 /(dt*vol) ;}
// //
// //
// //            fac_tear = 100000.0;
// //            //more intense one as well, also tear shaped
// //            che = pow((centroid_x[j] - midx),2) - fac_tear*pow((0.75*1.0e-3-centroid_y[i]+midy-0.75*1.0e-3),3)*(0.5*1.0e-3+centroid_y[i]-midy+0.75*1.0e-3);
// //            if (che <= 0.0){source_accu[i][j] =source_accu[i][j] + 20* 1.0e-3 /(dt*vol) ;}


//              /////// Circle at the middle lobe of tear
//             che = pow((centroid_x[j] - midx),2) + pow((centroid_y[i]-midy+0.85*1.0e-3),2) - pow(0.55*1.0e-3/2,2) ;
//             if (che <= 0.0){source_accu[i][j] =source_accu[i][j]+ 35* 1.0e-3 /(dt*vol) ;}
// TEAR SHAPE ENDS ---------------------------------------------------------------------------------------------





        // LASER INTENSITY USED HERE WITH GAUSSIAN DISTRIBUTION:: DATA OF SECOND PULSE
        // Source accu = Intensity ( watt per unit area :: checks out)
        // FITTED TO A FUNCTION in 2D ------------------------------------------------
        // Define the profile with these parameters
        //NOTE: xmid and ymid are not changed
        // At t = 0
         //  3.01549875e+11 -2.07422001e-01 -2.21465552e-06  9.38443224e-01  -4.69878748e-01
         //  Amp                x0                  y0             sigma_x          sigma_y

        // gaus1.A = 3.01549875e+11 ;
        // gaus1.x0 = -2.07422001e-01 ; 
        // gaus1.y0 = -2.21465552e-06 ;
        // gaus1.sigma_x = 9.38443224e-01 ;
        // gaus1.sigma_y = -4.69878748e-01 ;


        /*
        Fitted parameters: For Last TIME SNAP: t = 10ns :: 4 gaussian :: x0 , y0, sigmax, sigmay, amplitude
        [-5.09449565e-01  4.36821396e-01  7.09271093e-01  3.59628433e-01   1.99005213e+13 
        -5.09447698e-01 -4.36819324e-01  7.09275248e-01   3.59631310e-01  1.99005492e+13  
        1.24063233e+00  2.39198134e-05   8.80728819e-01  1.19575633e+00  6.58660872e+12 
        -2.05679456e+00   7.13013373e-06  6.59353934e-01  1.24931801e+00  6.10483100e+12]
        */


 // For only intensity
    //    gaus1.A = 1.99007486e+13 ; gaus1.x0 = -5.09472709e-01 ; gaus1.y0 = 4.36820539e-01; gaus1.sigma_x = 7.09271093e-01; gaus1.sigma_y = 3.59628433e-01 ;
    //    gaus2.A = 1.99007062e+13 ; gaus2.x0 = -5.09467920e-01 ; gaus2.y0 =-4.36822323e-01; gaus2.sigma_x = 7.09275248e-01 ; gaus2.sigma_y =  3.59631310e-01 ;
    //    gaus3.A =6.58660628e+12 ; gaus3.x0 = 1.24061914e+00; gaus3.y0 = 2.39198134e-05 ; gaus3.sigma_x = 8.80728819e-01; gaus3.sigma_y = 1.19575633e+00 ;
    //    gaus4.A =6.10485144e+12 ; gaus4.x0 = -2.05679838e+00; gaus4.y0 =7.13013373e-06 ; gaus4.sigma_x = 6.59353934e-01; gaus4.sigma_y = 1.24931801e+00;


         //// For I_Pulse2_0*Rho_e_0
        gaus1.x0 = -2.470534e-01 ; gaus1.y0 = 2.465843e-07 ;  gaus1.sigma_x = 5.314265e-01 ; gaus1.sigma_y = 3.952746e-01 ;  gaus1.A = 3.215602e+25 ; 
        gaus2.x0 = 7.488499e-01 ; gaus2.y0 = 1.407520e-06 ;  gaus2.sigma_x = 4.782765e-01 ; gaus2.sigma_y = 7.270592e-01 ;  gaus2.A = 8.904928e+24 ; 
        gaus3.x0 = -1.557827e+00 ; gaus3.y0 = 1.342642e-05 ;  gaus3.sigma_x = 6.963636e-01 ; gaus3.sigma_y = 9.185818e-01 ;  gaus3.A = 7.962736e+24 ; 
        gaus4.x0 = 1.913055e+00 ; gaus4.y0 = -1.765879e-05 ;  gaus4.sigma_x = -6.388076e-01 ; gaus4.sigma_y = -1.268513e+00 ;  gaus4.A = 4.262157e+24 ;


        // NOTE: The Gaussian Fit is defined for Z/Releygh_range ::> X and R/w0 ::> Y
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus1.A, gaus1.x0 , gaus1.y0, gaus1.sigma_x, gaus1.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus2.A, gaus2.x0 , gaus2.y0, gaus2.sigma_x, gaus2.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus3.A, gaus3.x0 , gaus3.y0, gaus3.sigma_x, gaus3.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus4.A, gaus4.x0 , gaus4.y0, gaus4.sigma_x, gaus4.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;

        }
    }
}else if (t<2.0*5.0e-9){
    // Different profiles in Time

        for ( i = 0 ; i < nc_row ; i++){
            for( j = 0 ; j< nc_col ; j++){
        // source_accu[0,0] = 0.851072 / dt ## 0.244816/dt

       gaus1.A = 1.99007486e+13 ; gaus1.x0 = -5.09472709e-01 ; gaus1.y0 = 4.36820539e-01; gaus1.sigma_x = 7.09271093e-01; gaus1.sigma_y = 3.59628433e-01 ;
       gaus2.A = 1.99007062e+13 ; gaus2.x0 = -5.09467920e-01 ; gaus2.y0 =-4.36822323e-01; gaus2.sigma_x = 7.09275248e-01 ; gaus2.sigma_y =  3.59631310e-01 ;
       gaus3.A =6.58660628e+12 ; gaus3.x0 = 1.24061914e+00; gaus3.y0 = 2.39198134e-05 ; gaus3.sigma_x = 8.80728819e-01; gaus3.sigma_y = 1.19575633e+00 ;
       gaus4.A =6.10485144e+12 ; gaus4.x0 = -2.05679838e+00; gaus4.y0 =7.13013373e-06 ; gaus4.sigma_x = 6.59353934e-01; gaus4.sigma_y = 1.24931801e+00;


        // NOTE: The Gaussian Fit is defined for Z/Releygh_range ::> X and R/w0 ::> Y
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus1.A, gaus1.x0 , gaus1.y0, gaus1.sigma_x, gaus1.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus2.A, gaus2.x0 , gaus2.y0, gaus2.sigma_x, gaus2.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus3.A, gaus3.x0 , gaus3.y0, gaus3.sigma_x, gaus3.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus4.A, gaus4.x0 , gaus4.y0, gaus4.sigma_x, gaus4.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;

        }
    }
}else if (t<3.0*5.0e-9){

        for ( i = 0 ; i < nc_row ; i++){
            for( j = 0 ; j< nc_col ; j++){
        // source_accu[0,0] = 0.851072 / dt ## 0.244816/dt

       gaus1.A = 1.99007486e+13 ; gaus1.x0 = -5.09472709e-01 ; gaus1.y0 = 4.36820539e-01; gaus1.sigma_x = 7.09271093e-01; gaus1.sigma_y = 3.59628433e-01 ;
       gaus2.A = 1.99007062e+13 ; gaus2.x0 = -5.09467920e-01 ; gaus2.y0 =-4.36822323e-01; gaus2.sigma_x = 7.09275248e-01 ; gaus2.sigma_y =  3.59631310e-01 ;
       gaus3.A =6.58660628e+12 ; gaus3.x0 = 1.24061914e+00; gaus3.y0 = 2.39198134e-05 ; gaus3.sigma_x = 8.80728819e-01; gaus3.sigma_y = 1.19575633e+00 ;
       gaus4.A =6.10485144e+12 ; gaus4.x0 = -2.05679838e+00; gaus4.y0 =7.13013373e-06 ; gaus4.sigma_x = 6.59353934e-01; gaus4.sigma_y = 1.24931801e+00;


        // NOTE: The Gaussian Fit is defined for Z/Releygh_range ::> X and R/w0 ::> Y
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus1.A, gaus1.x0 , gaus1.y0, gaus1.sigma_x, gaus1.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus2.A, gaus2.x0 , gaus2.y0, gaus2.sigma_x, gaus2.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus3.A, gaus3.x0 , gaus3.y0, gaus3.sigma_x, gaus3.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus4.A, gaus4.x0 , gaus4.y0, gaus4.sigma_x, gaus4.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;

        }
    }
}else if (t<4.0*5.0e-9){

// Different profiles in Time


        for ( i = 0 ; i < nc_row ; i++){
            for( j = 0 ; j< nc_col ; j++){
        // source_accu[0,0] = 0.851072 / dt ## 0.244816/dt

       gaus1.A = 1.99007486e+13 ; gaus1.x0 = -5.09472709e-01 ; gaus1.y0 = 4.36820539e-01; gaus1.sigma_x = 7.09271093e-01; gaus1.sigma_y = 3.59628433e-01 ;
       gaus2.A = 1.99007062e+13 ; gaus2.x0 = -5.09467920e-01 ; gaus2.y0 =-4.36822323e-01; gaus2.sigma_x = 7.09275248e-01 ; gaus2.sigma_y =  3.59631310e-01 ;
       gaus3.A =6.58660628e+12 ; gaus3.x0 = 1.24061914e+00; gaus3.y0 = 2.39198134e-05 ; gaus3.sigma_x = 8.80728819e-01; gaus3.sigma_y = 1.19575633e+00 ;
       gaus4.A =6.10485144e+12 ; gaus4.x0 = -2.05679838e+00; gaus4.y0 =7.13013373e-06 ; gaus4.sigma_x = 6.59353934e-01; gaus4.sigma_y = 1.24931801e+00;


        // NOTE: The Gaussian Fit is defined for Z/Releygh_range ::> X and R/w0 ::> Y
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus1.A, gaus1.x0 , gaus1.y0, gaus1.sigma_x, gaus1.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus2.A, gaus2.x0 , gaus2.y0, gaus2.sigma_x, gaus2.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus3.A, gaus3.x0 , gaus3.y0, gaus3.sigma_x, gaus3.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus4.A, gaus4.x0 , gaus4.y0, gaus4.sigma_x, gaus4.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;

        }
    }
}else if (t<5.0*5.0e-9){

// Different profiles in Time


        for ( i = 0 ; i < nc_row ; i++){
            for( j = 0 ; j< nc_col ; j++){
        // source_accu[0,0] = 0.851072 / dt ## 0.244816/dt

       gaus1.A = 1.99007486e+13 ; gaus1.x0 = -5.09472709e-01 ; gaus1.y0 = 4.36820539e-01; gaus1.sigma_x = 7.09271093e-01; gaus1.sigma_y = 3.59628433e-01 ;
       gaus2.A = 1.99007062e+13 ; gaus2.x0 = -5.09467920e-01 ; gaus2.y0 =-4.36822323e-01; gaus2.sigma_x = 7.09275248e-01 ; gaus2.sigma_y =  3.59631310e-01 ;
       gaus3.A =6.58660628e+12 ; gaus3.x0 = 1.24061914e+00; gaus3.y0 = 2.39198134e-05 ; gaus3.sigma_x = 8.80728819e-01; gaus3.sigma_y = 1.19575633e+00 ;
       gaus4.A =6.10485144e+12 ; gaus4.x0 = -2.05679838e+00; gaus4.y0 =7.13013373e-06 ; gaus4.sigma_x = 6.59353934e-01; gaus4.sigma_y = 1.24931801e+00;


        // NOTE: The Gaussian Fit is defined for Z/Releygh_range ::> X and R/w0 ::> Y
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus1.A, gaus1.x0 , gaus1.y0, gaus1.sigma_x, gaus1.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus2.A, gaus2.x0 , gaus2.y0, gaus2.sigma_x, gaus2.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus3.A, gaus3.x0 , gaus3.y0, gaus3.sigma_x, gaus3.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus4.A, gaus4.x0 , gaus4.y0, gaus4.sigma_x, gaus4.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;

        }
    }
}else if (t<6.0*5.0e-9){

// Different profiles in Time


        for ( i = 0 ; i < nc_row ; i++){
            for( j = 0 ; j< nc_col ; j++){
        // source_accu[0,0] = 0.851072 / dt ## 0.244816/dt

       gaus1.A = 1.99007486e+13 ; gaus1.x0 = -5.09472709e-01 ; gaus1.y0 = 4.36820539e-01; gaus1.sigma_x = 7.09271093e-01; gaus1.sigma_y = 3.59628433e-01 ;
       gaus2.A = 1.99007062e+13 ; gaus2.x0 = -5.09467920e-01 ; gaus2.y0 =-4.36822323e-01; gaus2.sigma_x = 7.09275248e-01 ; gaus2.sigma_y =  3.59631310e-01 ;
       gaus3.A =6.58660628e+12 ; gaus3.x0 = 1.24061914e+00; gaus3.y0 = 2.39198134e-05 ; gaus3.sigma_x = 8.80728819e-01; gaus3.sigma_y = 1.19575633e+00 ;
       gaus4.A =6.10485144e+12 ; gaus4.x0 = -2.05679838e+00; gaus4.y0 =7.13013373e-06 ; gaus4.sigma_x = 6.59353934e-01; gaus4.sigma_y = 1.24931801e+00;


        // NOTE: The Gaussian Fit is defined for Z/Releygh_range ::> X and R/w0 ::> Y
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus1.A, gaus1.x0 , gaus1.y0, gaus1.sigma_x, gaus1.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus2.A, gaus2.x0 , gaus2.y0, gaus2.sigma_x, gaus2.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus3.A, gaus3.x0 , gaus3.y0, gaus3.sigma_x, gaus3.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus4.A, gaus4.x0 , gaus4.y0, gaus4.sigma_x, gaus4.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;

        }
    }
}else if (t<7.0*5.0e-9){

// Different profiles in Time


        for ( i = 0 ; i < nc_row ; i++){
            for( j = 0 ; j< nc_col ; j++){
        // source_accu[0,0] = 0.851072 / dt ## 0.244816/dt

       gaus1.A = 1.99007486e+13 ; gaus1.x0 = -5.09472709e-01 ; gaus1.y0 = 4.36820539e-01; gaus1.sigma_x = 7.09271093e-01; gaus1.sigma_y = 3.59628433e-01 ;
       gaus2.A = 1.99007062e+13 ; gaus2.x0 = -5.09467920e-01 ; gaus2.y0 =-4.36822323e-01; gaus2.sigma_x = 7.09275248e-01 ; gaus2.sigma_y =  3.59631310e-01 ;
       gaus3.A =6.58660628e+12 ; gaus3.x0 = 1.24061914e+00; gaus3.y0 = 2.39198134e-05 ; gaus3.sigma_x = 8.80728819e-01; gaus3.sigma_y = 1.19575633e+00 ;
       gaus4.A =6.10485144e+12 ; gaus4.x0 = -2.05679838e+00; gaus4.y0 =7.13013373e-06 ; gaus4.sigma_x = 6.59353934e-01; gaus4.sigma_y = 1.24931801e+00;


        // NOTE: The Gaussian Fit is defined for Z/Releygh_range ::> X and R/w0 ::> Y
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus1.A, gaus1.x0 , gaus1.y0, gaus1.sigma_x, gaus1.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus2.A, gaus2.x0 , gaus2.y0, gaus2.sigma_x, gaus2.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus3.A, gaus3.x0 , gaus3.y0, gaus3.sigma_x, gaus3.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus4.A, gaus4.x0 , gaus4.y0, gaus4.sigma_x, gaus4.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;

        }
    }
}else if (t<8.0*5.0e-9){



        for ( i = 0 ; i < nc_row ; i++){
            for( j = 0 ; j< nc_col ; j++){
        // source_accu[0,0] = 0.851072 / dt ## 0.244816/dt

    //    gaus1.A = 1.99007486e+13 ; gaus1.x0 = -5.09472709e-01 ; gaus1.y0 = 4.36820539e-01; gaus1.sigma_x = 7.09271093e-01; gaus1.sigma_y = 3.59628433e-01 ;
    //    gaus2.A = 1.99007062e+13 ; gaus2.x0 = -5.09467920e-01 ; gaus2.y0 =-4.36822323e-01; gaus2.sigma_x = 7.09275248e-01 ; gaus2.sigma_y =  3.59631310e-01 ;
    //    gaus3.A =6.58660628e+12 ; gaus3.x0 = 1.24061914e+00; gaus3.y0 = 2.39198134e-05 ; gaus3.sigma_x = 8.80728819e-01; gaus3.sigma_y = 1.19575633e+00 ;
    //    gaus4.A =6.10485144e+12 ; gaus4.x0 = -2.05679838e+00; gaus4.y0 =7.13013373e-06 ; gaus4.sigma_x = 6.59353934e-01; gaus4.sigma_y = 1.24931801e+00;



        // NOTE: The Gaussian Fit is defined for Z/Releygh_range ::> X and R/w0 ::> Y
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus1.A, gaus1.x0 , gaus1.y0, gaus1.sigma_x, gaus1.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus2.A, gaus2.x0 , gaus2.y0, gaus2.sigma_x, gaus2.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus3.A, gaus3.x0 , gaus3.y0, gaus3.sigma_x, gaus3.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;
        source_accu[i][j] = source_accu[i][j]+gauss_2d(gaus4.A, gaus4.x0 , gaus4.y0, gaus4.sigma_x, gaus4.sigma_y, (centroid_x[j] - midx)/Rl_rn, (centroid_y[i] - midy)/w0) ;

        }
    }
}







         printf("SOURCE = %f\n", source_accu[(int)(nc_row/2)][(int)(nc_col/2)] );
        free(drex); free(drey); free(centroid_x); free(centroid_y);


}


// Function Defining the Energy addition Profile :: GAUSSIAN HERE
double gauss_2d(double A,double x0, double y0,double sigmax, double sigmay, double x, double y){
    // Remember the order of the parameters that are fitted
    double res ;

    //## Exactly Gaussian 2D
    res = A * exp(- (x - x0)*(x - x0)/(2*sigmax*sigmax) - (y - y0)*(y - y0)/(2*sigmay*sigmay) );
    return res;
}

// Lax- Fedrich



void LFflux(int nc_row , int nc_col , double Qil[4][nc_row][nc_col],double Qir[4][nc_row][nc_col],double GAMMA, char DIR ,double mspc[nc_row][nc_col], int f_row , int f_col, double Flux[4][f_row][f_col]){


// LFflux(qL,qR,gamma,normal,smax)
//    % Lax-Friedrichs flux
//
//    % Normal vectors
//    nx = normal(1);
//    ny = normal(2);
//
//    % Left state
//    rL = qL(1);
//    uL = qL(2)/rL;
//    vL = qL(3)/rL;
//    vnL= uL*nx + vL*ny;
//    pL = (gamma-1)*( qL(4) - 0.5*rL*(uL^2+vL^2) );
//    HL = ( qL(4) + pL ) / rL;
//
//    % Right state
//    rR = qR(1);
//    uR = qR(2)/rR;
//    vR = qR(3)/rR;
//    vnR= uR*nx + vR*ny;
//    pR = (gamma-1)*( qR(4) - 0.5*rR*(uR^2+vR^2) );
//    HR = ( qR(4) + pR ) / rR;
//
//    % Left and Right fluxes
//    FL=[rL*vnL; rL*vnL*uL + pL*nx; rL*vnL*vL + pL*ny; rL*vnL*HL];
//    FR=[rR*vnR; rR*vnR*uR + pR*nx; rR*vnR*vR + pR*ny; rR*vnR*HR];
//
//    % Rusanov numerical flux
//    LF = 0.5*( FR + FL + smax*(qL-qR) );
//end



	int  inter_row , inter_col , real_cell_row , real_cell_col , real_inter_row , real_inter_col ;//nc_row , nc_col ;

	int i, j, k, l ,m , n ;

	// number of cells after adding ghost cells
//	nc_row = sizeof(Qil[0])/sizeof(Qil[0][0]);
//	nc_col = sizeof(Qil[0][0])/sizeof(Qil[0][0][0]);

		 // Number of x-interfaces.
		 // Number of y-interfaces.


    inter_row = nc_row + 1;
    inter_col = nc_col + 1;
    real_cell_row = nc_row - 4;
    real_cell_col = nc_col - 4;
    real_inter_row = nc_row - 4 + 1 ;//## number of interfaces in real = total including ghost cells - subtract ghosts + add one
    real_inter_col = nc_col - 4 + 1 ;



    double CSR[4] = {0.0};
    double CSL[4] = {0.0};


//    double Flux[4][real_cell_row][real_inter_col] // Define this while passing from parent function

//            double qil[4][real_cell_row][real_inter_col] ;
//            double qir[4][real_cell_row][real_inter_col] ;
//
//
//            double fil[4][real_cell_row][real_inter_col] ;//= {0.0};
//            double fir[4][real_cell_row][real_inter_col] ;//= {0.0};
            double rL,uL,vL,pL,rR,uR,vR,pR;

//            double RL[],UL[],VL[],PL[],SPL[],RR[],UR[],VR[],PR[], SPR[];

            double rl,ul,vl,pl,rr,ur,vr,pr, ENEL, ENER;

            double SL , SM , SR ;

            double T_QL[4], T_QR[4] ;

    if (DIR == 'x'){

        double (*qil)[real_cell_row][real_inter_col] = malloc(sizeof(double[4][real_cell_row][real_inter_col])) ;//= {0.0};
		double (*qir)[real_cell_row][real_inter_col] = malloc(sizeof(double[4][real_cell_row][real_inter_col])) ;//= {0.0};

		double (*fil)[real_cell_row][real_inter_col] = malloc(sizeof(double[4][real_cell_row][real_inter_col])) ;//= {0.0};
		double (*fir)[real_cell_row][real_inter_col] = malloc(sizeof(double[4][real_cell_row][real_inter_col])) ;//= {0.0};
		double (*RL)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col])),
                (*UL)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col])),
                (*VL)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col])),
                (*PL)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col])),
                (*SPL)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col])),
                (*RR)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col])),
                (*UR)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col])),
                (*VR)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col])),
                (*PR)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col])), 
                (*SPR)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col]));

        double (*smax)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col])) ;

       //  Flux calculation requires left and right values on the interface of diffrent but adjecent to the interface cells
        for (i=2; i< 2 + real_cell_row; i++){ //## change this to make consistent with y as well
            for (j = 2 ; j < inter_col-2; j++ ) {

                //printf(i,j);
                m = i-2 ;
                n = j-2 ;
                for ( k = 0 ; k < 4 ; k++){
                qil[k][m][n] = Qir[k][i][j-1] ; // for left from interface is right of upstream cell
                qir[k][m][n] = Qil[k][i][j] ;// for right from interface is left for downstream cell


                smax[m][n] = mspc[i][j];

                T_QL[k] = qil[k][m][n];
                 T_QR[k] = qir[k][m][n];





				}

				sing_prmcalculate(T_QL,GAMMA,&rL,&uL,&vL,&pL);
				sing_prmcalculate(T_QR,GAMMA,&rR,&uR,&vR,&pR );


                fil[0][m][n] = rL*uL;
				fil[1][m][n] = rL*uL*uL+pL ;
				fil[2][m][n] = rL*uL*vL ;
				fil[3][m][n] = uL*((qil[3][m][n])+pL) ;


				fir[0][m][n] = rR*uR;
				fir[1][m][n] = rR*uR*uR+pR ;
				fir[2][m][n] = rR*uR*vR ;
				fir[3][m][n] = uR*(qir[3][m][n]+pR) ;


                }
        }

        // prmcalculate(int nc_row, int nc_col,		double q[][nc_row][nc_col],double gamma,
        // double r[][nc_col],double u[][nc_col],double v[][nc_col],double p[][nc_col],double s[][nc_col])

//        prmcalculate(real_cell_row , real_inter_col, qil,GAMMA , RL,UL,VL,PL,SPL);
//        prmcalculate(real_cell_row , real_inter_col, qir,GAMMA, RR,UR,VR,PR,SPR );


        // ---------

        for (i=0; i<real_cell_row; i++){
                      for (j = 0; j < real_inter_col; j++){


                      m = i ;
                    n = j ;


                      for ( k = 0 ; k< 4 ; k++){
                              Flux[k][m][n] = 0.5*(fil[k][i][j] + fir[k][i][j] + smax[i][j] * ( qil[k][i][j] - qir[k][i][j])) ;
                          	}

                      }
        }
    
    free(qil); free(qir);
    free(fil); free(fir);
    free(RL); free(UL); free(VL) ; free(PL) ; free(SPL);
    free(RR); free(UR); free(VR) ; free(PR) ; free(SPR);
    free(smax);
    
    }


    else{

    	double (*qil)[real_inter_row][real_cell_col]= malloc(sizeof(double[4][real_inter_row][real_cell_col])) ;
		double (*qir)[real_inter_row][real_cell_col]= malloc(sizeof(double[4][real_inter_row][real_cell_col])) ;


		double (*fil)[real_inter_row][real_cell_col]= malloc(sizeof(double[4][real_inter_row][real_cell_col])) ;//= {0.0};
		double (*fir)[real_inter_row][real_cell_col]= malloc(sizeof(double[4][real_inter_row][real_cell_col])) ;//= {0.0};
		double (*RL)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col])),
                (*UL)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col])),
                (*VL)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col])),
                (*PL)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col])),
                (*SPL)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col])),
                (*RR)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col])),
                (*UR)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col])),
                (*VR)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col])),
                (*PR)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col])), 
                (*SPR)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col]));


		double (*smax)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col])) ;

    	            // Flux calculation requires left and right values on the interface of diffrent but adjecent to the interface cells


    	// -----------------------------
    // CHECK HERE	double Flux[4][real_cell_row][real_inter_col] // Define this while passing from parent function

//    	        double qil[4][real_cell_row][real_inter_col] = {0.0};
//    	        double qir[4][real_cell_row][real_inter_col] = {0.0};
//
//
//    	        double fil[4][real_cell_row][real_inter_col] = {0.0};
//    	        double fir[4][real_cell_row][real_inter_col] = {0.0};
//
//    	        double rL,uL,vL,pL,rR,uR,vR,pR;
//
//    	        double RL,UL,VL,PL,RR,UR,VR,PR;



    	       //  Flux calculation requires left and right values on the interface of diffrent but adjecent to the interface cells
    	        for (i=2; i< inter_row - 2; i++){ //## change this to make consistent with y as well
    	            for (j = 2 ; j < nc_col-2; j++ ) {

    	                //printf(i,j);
    	                m = i-2 ;
    	                n = j-2 ;
    	                for ( k = 0 ; k < 4 ; k++){
    	                qil[k][m][n] = Qir[k][i-1][j] ; // for left from interface is right of upstream cell
    	                qir[k][m][n] = Qil[k][i][j] ;// for right from interface is left for downstream cell

    	                T_QL[k] = qil[k][m][n];
    	                T_QR[k] = qir[k][m][n];


                        smax[m][n] = mspc[i][j];


    	                }

    	                sing_prmcalculate(T_QL,GAMMA,&rL,&uL,&vL,&pL);
    	                sing_prmcalculate(T_QR,GAMMA,&rR,&uR,&vR,&pR );


    	                fil[0][m][n] = rL*vL;
    					fil[1][m][n] = rL*vL*uL ;
    					fil[2][m][n] = rL*vL*vL +pL ;
    					fil[3][m][n] = vL*(qil[3][m][n]+pL) ;


    					fir[0][m][n] = rR*vR;
    					fir[1][m][n] = rR*vR*uR ;
    					fir[2][m][n] = rR*vR*vR +pR ;
    					fir[3][m][n] = vR*(qir[3][m][n]+pR) ;


    	                }
    	        }
//    	        prmcalculate(real_inter_row,real_cell_col,qil,GAMMA , RL,UL,VL,PL,SPL);
//    	        prmcalculate(real_inter_row,real_cell_col,qir,GAMMA, RR,UR,VR,PR,SPR );


    	        // ------------

                for (i=0 ; i< real_inter_row; i++){
                        for (j = 0; j < real_cell_col; j++){

                                m = i;
                                n = j;


                                for ( k = 0 ; k< 4 ; k++){
                              Flux[k][m][n] = 0.5*(fil[k][i][j] + fir[k][i][j] + smax[i][j] * ( qil[k][i][j] - qir[k][i][j])) ;
                          	}


                }

    	        // ------------

    	    }


free(qil); free(qir);
    free(fil); free(fir);
    free(RL); free(UL); free(VL) ; free(PL) ; free(SPL);
    free(RR); free(UR); free(VR) ; free(PR) ; free(SPR);
    free(smax);

}

    


}


/// RUSONOV




void RSflux(int nc_row , int nc_col , double Qil[4][nc_row][nc_col],double Qir[4][nc_row][nc_col],double GAMMA, char DIR ,int f_row , int f_col, double Flux[4][f_row][f_col]){


	int  inter_row , inter_col , real_cell_row , real_cell_col , real_inter_row , real_inter_col ;//nc_row , nc_col ;

	int i, j, k, l ,m , n ;

	// number of cells after adding ghost cells
//	nc_row = sizeof(Qil[0])/sizeof(Qil[0][0]);
//	nc_col = sizeof(Qil[0][0])/sizeof(Qil[0][0][0]);

		 // Number of x-interfaces.
		 // Number of y-interfaces.


    inter_row = nc_row + 1;
    inter_col = nc_col + 1;
    real_cell_row = nc_row - 4;
    real_cell_col = nc_col - 4;
    real_inter_row = nc_row - 4 + 1 ;//## number of interfaces in real = total including ghost cells - subtract ghosts + add one
    real_inter_col = nc_col - 4 + 1 ;



    double CSR[4] = {0.0};
    double CSL[4] = {0.0};


//    double Flux[4][real_cell_row][real_inter_col] // Define this while passing from parent function

//            double qil[4][real_cell_row][real_inter_col] ;
//            double qir[4][real_cell_row][real_inter_col] ;
//
//
//            double fil[4][real_cell_row][real_inter_col] ;//= {0.0};
//            double fir[4][real_cell_row][real_inter_col] ;//= {0.0};
            double rL,uL,vL,pL,rR,uR,vR,pR;

//            double RL[],UL[],VL[],PL[],SPL[],RR[],UR[],VR[],PR[], SPR[];

            double rl,ul,vl,pl,rr,ur,vr,pr, ENEL, ENER;

            double SL , SM , SR ;

            double T_QL[4], T_QR[4] ;

            double HL, HR , RT , u , v , a , H ;

    if (DIR == 'x'){

        double (*qil)[real_cell_row][real_inter_col] = malloc(sizeof(double[4][real_cell_row][real_inter_col])) ;//= {0.0};
		double (*qir)[real_cell_row][real_inter_col] = malloc(sizeof(double[4][real_cell_row][real_inter_col])) ;//= {0.0};

    	double (*fil)[real_cell_row][real_inter_col] = malloc(sizeof(double[4][real_cell_row][real_inter_col])) ;//= {0.0};
		double (*fir)[real_cell_row][real_inter_col] = malloc(sizeof(double[4][real_cell_row][real_inter_col])) ;//= {0.0};
		double (*RL)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col])),
                (*UL)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col])),
                (*VL)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col])),
                (*PL)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col])),
                (*SPL)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col])),
                (*RR)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col])),
                (*UR)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col])),
                (*VR)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col])),
                (*PR)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col])), 
                (*SPR)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col]));

        double (*smax)[real_inter_col] = malloc(sizeof(double[real_cell_row][real_inter_col])) ;

       //  Flux calculation requires left and right values on the interface of diffrent but adjecent to the interface cells
        for (i=2; i< 2 + real_cell_row; i++){ //## change this to make consistent with y as well
            for (j = 2 ; j < inter_col-2; j++ ) {

                //printf(i,j);
                m = i-2 ;
                n = j-2 ;
                for ( k = 0 ; k < 4 ; k++){
                qil[k][m][n] = Qir[k][i][j-1] ; // for left from interface is right of upstream cell
                qir[k][m][n] = Qil[k][i][j] ;// for right from interface is left for downstream cell


          //      smax[m][n] = mspc[i][j];

                T_QL[k] = qil[k][m][n];
                 T_QR[k] = qir[k][m][n];





				}

				sing_prmcalculate(T_QL,GAMMA,&rL,&uL,&vL,&pL);
				sing_prmcalculate(T_QR,GAMMA,&rR,&uR,&vR,&pR );


                fil[0][m][n] = rL*uL;
				fil[1][m][n] = rL*uL*uL+pL ;
				fil[2][m][n] = rL*uL*vL ;
				fil[3][m][n] = uL*((qil[3][m][n])+pL) ;


				fir[0][m][n] = rR*uR;
				fir[1][m][n] = rR*uR*uR+pR ;
				fir[2][m][n] = rR*uR*vR ;
				fir[3][m][n] = uR*(qir[3][m][n]+pR) ;



				// Calculate Smax from Roes Average

				 // Calculate Roe Averages


                HL = ( qil[3][m][n] + pL ) / rL;
                HR = ( qir[3][m][n] + pR ) / rR;

                    RT = sqrt(rR/rL);
                    u = (uL+RT*uR)/(1+RT);
                    v = (vL+RT*vR)/(1+RT);
                    H = ( HL+RT* HR)/(1+RT);
                    a = sqrt( (GAMMA-1)*(H-(u*u+v*v)/2) );
                    smax[m][n] = fabs(sqrt(u*u + v*v)) + a ;



                }
        }

        // prmcalculate(int nc_row, int nc_col,		double q[][nc_row][nc_col],double gamma,
        // double r[][nc_col],double u[][nc_col],double v[][nc_col],double p[][nc_col],double s[][nc_col])

//        prmcalculate(real_cell_row , real_inter_col, qil,GAMMA , RL,UL,VL,PL,SPL);
//        prmcalculate(real_cell_row , real_inter_col, qir,GAMMA, RR,UR,VR,PR,SPR );


        // ---------










        for (i=0; i<real_cell_row; i++){
                      for (j = 0; j < real_inter_col; j++){


                      m = i ;
                    n = j ;





                      for ( k = 0 ; k< 4 ; k++){
                              Flux[k][m][n] = 0.5*(fil[k][i][j] + fir[k][i][j] + smax[i][j] * ( qil[k][i][j] - qir[k][i][j])) ;
                          	}

                      }
        }
    free(qil); free(qir);
    free(fil); free(fir);
    free(RL); free(UL); free(VL) ; free(PL) ; free(SPL);
    free(RR); free(UR); free(VR) ; free(PR) ; free(SPR);
    free(smax);  
    
    
    }


    else{

    	double (*qil)[real_inter_row][real_cell_col]= malloc(sizeof(double[4][real_inter_row][real_cell_col])) ;
		double (*qir)[real_inter_row][real_cell_col]= malloc(sizeof(double[4][real_inter_row][real_cell_col])) ;


		double (*fil)[real_inter_row][real_cell_col]= malloc(sizeof(double[4][real_inter_row][real_cell_col])) ;//= {0.0};
		double (*fir)[real_inter_row][real_cell_col]= malloc(sizeof(double[4][real_inter_row][real_cell_col])) ;//= {0.0};
		double (*RL)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col])),
                (*UL)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col])),
                (*VL)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col])),
                (*PL)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col])),
                (*SPL)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col])),
                (*RR)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col])),
                (*UR)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col])),
                (*VR)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col])),
                (*PR)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col])), 
                (*SPR)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col]));


		double (*smax)[real_cell_col]= malloc(sizeof(double[real_inter_row][real_cell_col])) ;

    	            // Flux calculation requires left and right values on the interface of diffrent but adjecent to the interface cells


    	// -----------------------------
    // CHECK HERE	double Flux[4][real_cell_row][real_inter_col] // Define this while passing from parent function

//    	        double qil[4][real_cell_row][real_inter_col] = {0.0};
//    	        double qir[4][real_cell_row][real_inter_col] = {0.0};
//
//
//    	        double fil[4][real_cell_row][real_inter_col] = {0.0};
//    	        double fir[4][real_cell_row][real_inter_col] = {0.0};
//
//    	        double rL,uL,vL,pL,rR,uR,vR,pR;
//
//    	        double RL,UL,VL,PL,RR,UR,VR,PR;



    	       //  Flux calculation requires left and right values on the interface of diffrent but adjecent to the interface cells
    	        for (i=2; i< inter_row - 2; i++){ //## change this to make consistent with y as well
    	            for (j = 2 ; j < nc_col-2; j++ ) {

    	                //printf(i,j);
    	                m = i-2 ;
    	                n = j-2 ;
    	                for ( k = 0 ; k < 4 ; k++){
    	                qil[k][m][n] = Qir[k][i-1][j] ; // for left from interface is right of upstream cell
    	                qir[k][m][n] = Qil[k][i][j] ;// for right from interface is left for downstream cell

    	                T_QL[k] = qil[k][m][n];
    	                T_QR[k] = qir[k][m][n];


                   //     smax[m][n] = mspc[i][j];


    	                }

    	                sing_prmcalculate(T_QL,GAMMA,&rL,&uL,&vL,&pL);
    	                sing_prmcalculate(T_QR,GAMMA,&rR,&uR,&vR,&pR );


    	                fil[0][m][n] = rL*vL;
    					fil[1][m][n] = rL*vL*uL ;
    					fil[2][m][n] = rL*vL*vL +pL ;
    					fil[3][m][n] = vL*(qil[3][m][n]+pL) ;


    					fir[0][m][n] = rR*vR;
    					fir[1][m][n] = rR*vR*uR ;
    					fir[2][m][n] = rR*vR*vR +pR ;
    					fir[3][m][n] = vR*(qir[3][m][n]+pR) ;




                HL = ( qil[3][m][n] + pL ) / rL;
                HR = ( qir[3][m][n] + pR ) / rR;

                    RT = sqrt(rR/rL);
                    u = (uL+RT*uR)/(1+RT);
                    v = (vL+RT*vR)/(1+RT);
                    H = ( HL+RT* HR)/(1+RT);
                    a = sqrt( (GAMMA-1)*(H-(u*u+v*v)/2) );
                    smax[m][n] = fabs(sqrt(u*u + v*v)) + a ;

    	                }
    	        }
//    	        prmcalculate(real_inter_row,real_cell_col,qil,GAMMA , RL,UL,VL,PL,SPL);
//    	        prmcalculate(real_inter_row,real_cell_col,qir,GAMMA, RR,UR,VR,PR,SPR );


    	        // ------------

                for (i=0 ; i< real_inter_row; i++){
                        for (j = 0; j < real_cell_col; j++){

                                m = i;
                                n = j;


                                for ( k = 0 ; k< 4 ; k++){
                              Flux[k][m][n] = 0.5*(fil[k][i][j] + fir[k][i][j] + smax[i][j] * ( qil[k][i][j] - qir[k][i][j])) ;
                          	}


                }

    	        // ------------

    	    }

    free(qil); free(qir);
    free(fil); free(fir);
    free(RL); free(UL); free(VL) ; free(PL) ; free(SPL);
    free(RR); free(UR); free(VR) ; free(PR) ; free(SPR);
    free(smax);


    }

}

void Geom_F(int nc_row, int nc_col,double q[4][nc_row][nc_col],double gamma,int NDIM, double x[], double y[], double source_geom[4][nc_row][nc_col]){
 /*   '''
---------------------------------------------------------------------------------------
 >>>>>>>>>>>           ### d/dt()dv = - N * ( G(U) / x ) dv ###          <<<<<<<<<<<<<<<
## Provide the source term contribution due to geometry change but provides only [N * ( G(U) / x )] this part; check sign while integration
N = 0 : Cartesian
N = 1 : Cylindrical
N = 2 : Spherical
`````````````````````````````````````````
Setup in such a way that; x -> r when N = 1 and y -> z
 -----------------------------------------------------------------------------------------
 '''  */

    int n_col, n_row, nx , ny ;

    n_col = nc_col + 1 ;
    n_row = nc_row + 1 ;

    nx = n_col;
    ny = n_row;
    
    int i, j ,k ;

    //double xflux[nc_row][n_col] , yflux[n_row][nc_col] ;
    double (*drex)=malloc(sizeof(double[nc_col])), (*drey)=malloc(sizeof(double[nc_row])) ,
     (*centroid_x)=malloc(sizeof(double[nc_col])) , (*centroid_y)=malloc(sizeof(double[nc_row])) , vol ,(*int_x)=malloc(sizeof(double[nc_col]));

    double (*r)[nc_col] = malloc(sizeof(double[nc_row][nc_col])),
            (*u)[nc_col] = malloc(sizeof(double[nc_row][nc_col])),
            (*v)[nc_col] = malloc(sizeof(double[nc_row][nc_col])),
            (*p)[nc_col] = malloc(sizeof(double[nc_row][nc_col])),
            (*E)[nc_col] = malloc(sizeof(double[nc_row][nc_col])),
            (*s)[nc_col] = malloc(sizeof(double[nc_row][nc_col])) ;
	 
	// double midx , midy ;

     //## Centroids
         for (i=0 ; i< nc_col ; i++){
            drex[i] = (x[i+1] - x[i]);
            centroid_x[i] = 0.5 * (x[i+1] + x[i]);
    }

    for(i=0 ; i < nc_row ; i++){
    	drey[i] = (y[i+1] - y[i]);
    	centroid_y[i] = 0.5*(y[i+1] + y[i]);

    }

    for ( i = 0 ; i < nc_row ; i++){
        for( j = 0 ; j< nc_col ; j++){
            for( k = 0 ; k< 4 ; k++){
            source_geom[k][i][j] = 0.0 ;
            }
        }
    }

	// midx = centroid_x[(int)(nc_col/2)];
	// midy = centroid_y[(int)(nc_row/2)];

    // void prmcalculate(int nc_row, int nc_col,double q[][nc_row][nc_col],double gamma,
    //                  double r[][nc_col],double u[][nc_col],double v[][nc_col],double p[][nc_col],double s[][nc_col]){
 
if(NDIM != 0){
    prmcalculate(nc_row,nc_col,q,gamma,r,u,v,p,s);

    for ( i = 0 ; i < nc_row ; i++){
        for( j = 0 ; j< nc_col ; j++){
            int_x[j] = log(x[j+1]) - log(x[j]) ;
            int_x[0] = log(centroid_x[0]);
            source_geom[0][i][j] = (r[i][j] * u[i][j])*(NDIM*int_x[j]);
            source_geom[1][i][j] = (r[i][j] * u[i][j] * u[i][j] )*(NDIM*int_x[j]);
            source_geom[2][i][j] = (r[i][j] * u[i][j] * v[i][j]) *(NDIM*int_x[j]);
            source_geom[3][i][j] = (u[i][j] * (q[3][i][j] + p[i][j]))*(NDIM*int_x[j]);
 
        }
    }
}

        free(drex); free(drey); free(centroid_x); free(centroid_y),free(int_x);
        free(r); free(u); free(v); free(p); free(E); free(s);


}