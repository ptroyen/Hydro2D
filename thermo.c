/*
 ============================================================================
 Name        : 2dsolver.c
 Author      : Sagar Pokharel
 Version     :
 Copyright   :
 Description : thermo.c
 ============================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "thermo.h"

// double eval_cp( double T);

#define PI 3.14159265358979323846
#define R0 8.3145*1000
#define MW_O2 32.0
#define MW_N2 28.0134
#define MW_M (0.80*MW_N2+0.2*MW_O2)
#define R (R0/MW_M)

#define J2KC 0.000239006
#define C2J 4.184

// int main(){


//     double cp, T ;

//     int i , j ,k ;


//     for (i=0 ; i <= 3000 ; i = i + 100){
//         T = i ;
//         cp = eval_cp(T) ;

//         printf("Tempr = %f \t Cp = %f\n", T , cp) ;
//     }



//  return 0;





// }


double eval_cp(double T){

    // Oxygen // Nitrogen

    double cp_species[2] ;
    double cp_mixture ;

    read_cp_i(T, cp_species) ;
    cp_mixture = (0.2*MW_O2 * cp_species[0] + 0.8*MW_N2 * cp_species[1])/(MW_M) ;

    assert(cp_mixture > 0.0);

    
    return cp_mixture; // Now in J/Kg-K
}  



void read_cp(double T , double cp[]){

    // Oxygen // Nitrogen

    if ( T > 6000.0){T = 6000 ; }

   if (T <= 1000.0) {
    cp[0] = 2.5983668450066881e+02 * (3.2129363999999998e+00 + T * (1.1274863499999999e-03 + T * (-5.7561504699999999e-07 + T * (1.3138772300000001e-09 + -8.7685539199999998e-13 * T))));
  } else {
    cp[0] = 2.5983668450066881e+02 * (3.6975781900000002e+00 + T * (6.1351968900000004e-04 + T * (-1.2588419900000000e-07 + T * (1.7752814800000000e-11 + -1.1364353100000000e-15 * T))));
  }

  
    if (T <= 1000.0) {
    cp[1] = 2.9680218594762238e+02 * (3.2986770000000001e+00 + T * (1.4082400000000001e-03 + T * (-3.9632219999999998e-06 + T * (5.6415150000000002e-09 + -2.4448550000000001e-12 * T))));
  } else {
    cp [1]= 2.9680218594762238e+02 * (2.9266399999999999e+00 + T * (1.4879769999999999e-03 + T * (-5.6847609999999995e-07 + T * (1.0097040000000000e-10 + -6.7533509999999998e-15 * T))));
  }

  

    // return cp;

}  // KJ/kg-K :: T in K


double iter_tempr(double e, double tempr_old){
/*
    Returns t as the solution for temperature
    Calculate Temperature from the conserved variables
    e : Internal Energy per unit mass
*/
	int i;
    double dT , h1 , h2 ,T , fn_min , slope_fn;
    double TOL , err ;
	
    T = tempr_old ;
    dT = 1.0e-12 ;
    // Tempr Iteration loop
    // for ( i = 0 ; i < 2 ; i++){

    //     fn_min = e + R*T - eval_cp(T)*T ;
    //     h2 = eval_cp(T+dT)*(T+dT); h1 = eval_cp(T-dT) * (T-dT) ;
    //     slope_fn = R - (h2 - h1) / (2*dT) ;

    //     T = T - (fn_min / slope_fn) ;

    //     // printf("Tempr in iteration = %f \n", T );

    // }

        TOL = 1.0e-3 ;
        err = 1.0 ;
        while(fabs(err) > TOL || T < 0.0 ){

        fn_min = e + R*T - eval_cp(T)*T ;
        h2 = eval_cp(T+dT)*(T+dT); h1 = eval_cp(T-dT) * (T-dT) ;
        slope_fn = R - (h2 - h1) / (2*dT) ;

        err = T ;  // old T
        T = T - (fn_min / slope_fn) ;

        err = (err - T )/ T ;

        // printf("Tempr in iteration = %f \n", T );

        if ( T < 0.0){

            printf("-ve temperature recorded in iteration \n");
        }

    }


    return T ;
}



void thermo_prmcalculate(int nc_row, int nc_col,double q[][nc_row][nc_col],double gamma,
                     double r[][nc_col],double u[][nc_col],double v[][nc_col],double p[][nc_col],double T[][nc_col],double s[][nc_col]){

	int i,j;
    double e , tempr_guess; 
	for (i=0; i < nc_row ; i++){
        for (j=0; j < nc_col ; j++){

            r[i][j] = q[0][i][j];
            u[i][j] = q[1][i][j]/r[i][j];
            v[i][j] = q[2][i][j]/r[i][j];
            s[i][j] = sqrt(u[i][j]*u[i][j] + v[i][j]*v[i][j]);

            p[i][j] = (1.4 - 1)*(q[3][i][j] - 0.5*r[i][j]*(u[i][j]*u[i][j] + v[i][j]*v[i][j])); // guess

            printf("pressure = %f\t",p[i][j]);

            tempr_guess = p[i][j] / (r[i][j]*R) ; // guess

            printf("Temperature = %f \n",tempr_guess);

            e = (q[3][i][j] / r[i][j]) - 0.5*s[i][j]*s[i][j] ; // mass specific internal energy

            T[i][j] = iter_tempr(e,tempr_guess) ;
            p[i][j] = r[i][j] * R * T[i][j] ;
            
		}
	}
}



void thermo_sing_prmcalculate(double q[],double guess_gamma, double *r,double *u,double *v,double *p, double *T){

    double e , tempr ,guess_p , vsq;

    *r = q[0] ;
    *u = q[1]/ (*r) ;
    *v = q[2]/ (*r) ;

    vsq = (*u)*(*u) + (*v)*(*v) ;
    guess_p = (guess_gamma - 1)*(q[3] - 0.5 * (*r) * ((*u) * (*u) + (*v) * (*v))) ;
    tempr = guess_p / (*r * R) ;
    e = (q[3] / (*r)) - 0.5 * vsq ;
    tempr = iter_tempr(e, tempr) ;

    *p = (*r) * R * tempr ;
    *T = tempr ;
}



void thermo_concalculate(int nc_row, int nc_col, double r[][nc_col],double u[][nc_col],double v[][nc_col],double p[][nc_col],double q[][nc_row][nc_col]){
    // assembles all properties in q, with conserved variables
	int i,j;
    double tempr , gamma, cp;

	for (i=0; i < nc_row ; i++){
			for (j=0; j < nc_col ; j++){
				q[0][i][j] = r[i][j];
				q[1][i][j] = u[i][j]*r[i][j];
				q[2][i][j] = v[i][j]*r[i][j];

        tempr = p[i][j] / (r[i][j] * R) ;
        cp = eval_cp(tempr) ;
        gamma = cp/(cp-R) ;
        // printf("Cp = %0.5f\n",cp);
        
        // printf("\ngamma = %0.4f , Tempr = %0.4f, Pressure = %0.4f , Cp = %0.4f", gamma,tempr,p[i][j], cp);

				q[3][i][j] = (p[i][j]/(gamma - 1)) + 0.5 * r[i][j] * ( u[i][j]*u[i][j] + v[i][j]*v[i][j] );

                // q[3][i][j] = (p[i][j]/(gamma - 1)) + 0.5 * r[i][j] * ( u[i][j]*u[i][j] + v[i][j]*v[i][j] ); // change gamma with temperature
			
			}
		}
	}



double thermo_CFLmaintain(int nc_row, int nc_col, double r[][nc_col],double u[][nc_col],double v[][nc_col],
                    double p[][nc_col],double CFL,double x[nc_col+1],double y[nc_row+1],int n){
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
	double dt , spc , mspc , dts , mspcx, mspcy , gamma , tempr , cp;

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

            
            tempr = p[i][j] / ( r[i][j] * R) ;

            // printf("Tempr = %f",tempr);

            
            cp = eval_cp(tempr) ;
            // printf("\ngamma = %0.4f , Tempr = %0.4f, Pressure = %0.4f , Cp = %0.4f", gamma,tempr,p[i][j], cp);

            assert(cp > 0.0); assert(tempr > 0.0);
            gamma = (cp)/(cp-R) ;

            // printf("\ngamma = %0.4f ", gamma);
            
            

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

        // dt = 1.0e-8; // FOR DMD DATA GENERATION
        // // Time Steps for first Initial Times
        // if (dt<=0.0 || !(isfinite(dt)) || dt > 0.1*1.0e-5 || n <= 5){
        //     printf("TIME STEP CALCULATED = %0.12e  \n" , dt);
		// 	dt = 1e-3 * 1.0e-5;

             if (n <= 500){
            // printf("TIME STEP CALCULATED = %0.12e  \n" , dt);
			dt = 1.0e-0 * 1.0e-9; }

            if (dt<=0.0 || !(isfinite(dt)) || dt > 1.0e-7 ){
            // printf("TIME STEP CALCULATED = %0.12e  \n" , dt);
			dt = 1.0e-8;}

    }


    // Free the allocated memory for dx and dy
    free(dx);
    free(dy);

    return dt;
}



void thermo_x_fun_flux(double qre[],double guess_gamma,double ff[]){
    //    ## Define the function that calculates the flux FOR one direction - X
    //    ## At a time flux of only 1 element in each row is calculated
    //    ## a = wave speed
    //    ## Returns ff

    double rre, ure, vre, pre , tre;

    thermo_sing_prmcalculate(qre,guess_gamma,&rre,&ure,&vre,&pre,&tre) ;

    ff[0] = qre[1] ; //## Contuinity Flux
    ff[1] = ff[0]*ure + pre ; //## X-momentum Flux
    ff[2] = ff[0]*vre;
    ff[3] = (ure)*(qre[3]+pre) ;
}


void thermo_y_fun_flux(double qre[],double guess_gamma,double ff[]){
    //    ## Define the function that calculates the flux FOR one direction - Y
    //    ## At a time flux of only 1 element in each row is calculated
    //    ## a = wave speed
    //    ## Returns ff

    double rre, ure, vre, pre ,tre;

    thermo_sing_prmcalculate(qre,guess_gamma,&rre,&ure,&vre,&pre,&tre) ;

    //    ## q[3] = E :: Totaal Energy Per Unit Volume
    //    ## FLUX FOR EACH EQUATION

    ff[0] = qre[2] ; //## Contuinity Flux
    ff[1] = ff[0]*ure ; //## X-momentum Flux
    ff[2] = ff[0]*vre  + pre  ;
    ff[3] = (vre)*(qre[3]+pre) ;
}



void thermo_RSflux(int nc_row , int nc_col , double Qil[4][nc_row][nc_col],double Qir[4][nc_row][nc_col],double guess_gamma[][nc_col], char DIR ,int f_row , int f_col, double Flux[4][f_row][f_col]){


	int  inter_row , inter_col , real_cell_row , real_cell_col , real_inter_row , real_inter_col ;//nc_row , nc_col ;

	int i, j, k, l ,m , n ;
  double GAMMA, cp , gammaL, gammaR;

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

            double rl,ul,vl,pl,rr,ur,vr,pr,tL,tR, ENEL, ENER;

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
                gammaL = guess_gamma[i][j-1];
                gammaR = guess_gamma[i][j];

          //      smax[m][n] = mspc[i][j];

                T_QL[k] = qil[k][m][n];
                T_QR[k] = qir[k][m][n];





				}

        // thermo_x_fun_flux(T_QL,guess_gamma,fil) ;
        // thermo_x_fun_flux(T_QR,guess_gamma,fir) ;

				thermo_sing_prmcalculate(T_QL,gammaL,&rL,&uL,&vL,&pL,&tL);
				thermo_sing_prmcalculate(T_QR,gammaR,&rR,&uR,&vR,&pR,&tR );

        cp = eval_cp(0.5*(tL + tR)) ;
        GAMMA = cp/(cp-R);

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
                    a = sqrt( (GAMMA-1)*(H-(u*u+v*v)/2) ); // GAMMA CHANGES FOR CALORICALLY IMPERFECT GAS
                    smax[m][n] = fabs(sqrt(u*u + v*v)) + a ;



                }
        }

        // prmcalculate(int nc_row, int nc_col,		double q[][nc_row][nc_col],double gamma,
        // double r[][nc_col],double u[][nc_col],double v[][nc_col],double p[][nc_col],double s[][nc_col])

//        prmcalculate(real_cell_row , real_inter_col, qil,GAMMA , RL,UL,VL,PL,SPL);
//        prmcalculate(real_cell_row , real_inter_col, qir,GAMMA, RR,UR,VR,PR,SPR );


        // --------

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

                      gammaL = guess_gamma[i-1][j];
                      gammaR = guess_gamma[i][j];

    	                T_QL[k] = qil[k][m][n];
    	                T_QR[k] = qir[k][m][n];

                     


                   //     smax[m][n] = mspc[i][j];


    	                }



                  

        // thermo_x_fun_flux(T_QL,guess_gamma,fil) ;
        // thermo_x_fun_flux(T_QR,guess_gamma,fir) ;

				thermo_sing_prmcalculate(T_QL,gammaL,&rL,&uL,&vL,&pL,&tL);
				thermo_sing_prmcalculate(T_QR,gammaR,&rR,&uR,&vR,&pR,&tR );

              cp = eval_cp(0.5*(tL + tR)) ;
              GAMMA = cp/(cp-R);


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


double eval_gamma_pr(double p,double r){

    double tempr, cp , GAMMA ;

    tempr = p / (r * R) ;
    cp = eval_cp(tempr);
    GAMMA = cp/(cp-R) ;

    return GAMMA ;


}


double eval_gamma_t(double tempr){

    double cp , GAMMA ;

    cp = eval_cp(tempr);
    GAMMA = cp/(cp-R) ;

    return GAMMA ;


}



void thermo_Flux_M(int nrec_row , int nrec_col, double qre[4][nrec_row][nrec_col] ,double pre_for_qre[][nrec_col],double xre[], double yre[] ,double dt,double guess_gamma,  double accu[4][nrec_row][nrec_col]){// ## add what is needed
//    '''
// For Calorically imperfect gases, So rather than iterating everytime to calculate the temperature for determining flux from cons variables pass P or T with the conservative variable as well
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
// MUSCL - HANCOCK 
//    ------------------------------------------------------------------------------------------------'''

    //    # initialization
	int i,j,k;

    //	int nrec_row, nrec_col;
	int nre_row,nre_col;

	int nc_row,nc_col,nx,ny;

  double GAMMA , cp , tempr;



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
                (*a)[nc_col] = malloc(sizeof(double[nc_row][nc_col])) ,
                (*cell_gamma)[nc_col] = malloc(sizeof(double[nc_row][nc_col]));

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

//## [][] supply row and column number to find the AREA in 2D of the mesh >>>>>>>>> USE FUNCTION
    for ( i = 0; i< nrec_row; i++ ){
    	for (j = 0; j < nrec_col; j++){
    		vol_re[i][j] = drex[j] * drey[i] ;

    		// Primitive Variables Calculation
    		rre[i][j] = qre[0][i][j];
    		ure[i][j] = qre[1][i][j] / rre[i][j];
    		vre[i][j] = qre[2][i][j] / rre[i][j];

    		// pre[i][j] = (gamma - 1)*(qre[3][i][j] - 0.5*rre[i][j]*(ure[i][j]*ure[i][j] + vre[i][j]*vre[i][j]));
        pre[i][j] = pre_for_qre[i][j] ;

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


// >>>>>>>>>>>>>>>>>>>> USE FUNCTIONS
//    ## assemble all props in q, with conserved variables
	for (i=0;i< nc_row; i++){
		for ( j = 0; j < nc_col; j++){
			q[0][i][j] = r[i][j];
			q[1][i][j] = u[i][j]*r[i][j];
			q[2][i][j] = v[i][j]*r[i][j];

        // tempr = p[i][j] / (r[i][j] * R) ;
        // cp = eval_cp(tempr);
        // GAMMA = cp/(cp-R) ;

        cell_gamma[i][j] = eval_gamma_pr(p[i][j],r[i][j]) ;
        GAMMA = cell_gamma[i][j];

			q[3][i][j] = (p[i][j]/(GAMMA - 1)) + 0.5 * r[i][j]*(u[i][j]*u[i][j] + v[i][j]*v[i][j]);

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

           GAMMA = cell_gamma[i][j] ;

        	spc[i][j] = sqrt(GAMMA*p[i][j]/r[i][j]);
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

            thermo_x_fun_flux(qlx, GAMMA, xfql);
            thermo_x_fun_flux(qrx,GAMMA, xfqr);

            thermo_y_fun_flux(qly, GAMMA, yfql);
            thermo_y_fun_flux(qry, GAMMA, yfqr);



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



// USE HLLC
//    HLLC(nc_row, nc_col, qilx,qirx,gamma,'x', nrec_row , nre_col,xflux);
//    HLLC(nc_row, nc_col ,qily,qiry,gamma,'y', nre_row, nrec_col,yflux);



 /// USE LAX-Ferdich
//    LFflux(nc_row, nc_col, qilx,qirx,gamma,'x', mspc, nrec_row , nre_col,xflux);
//    LFflux(nc_row, nc_col ,qily,qiry,gamma,'y',mspc, nre_row, nrec_col,yflux);

    /// Rusonov

    thermo_RSflux(nc_row, nc_col, qilx,qirx,cell_gamma,'x', nrec_row , nre_col,xflux);
    thermo_RSflux(nc_row, nc_col ,qily,qiry,cell_gamma,'y', nre_row, nrec_col,yflux);

    // RSflux(nc_row, nc_col, qilx,qirx,gamma,'x', nrec_row , nre_col,xflux);
    // RSflux(nc_row, nc_col ,qily,qiry,gamma,'y', nre_row, nrec_col,yflux);




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


void read_cp_i(double T , double cp[2]){


double C[2][3][11] = 
    {
            /* species O2
        pos 0: Tmin lower range limit
        pos 1: Tmax upper range limit
        pos 2-8: a1,a2,...,a7
        pos 9-10: b1,b2
        */    
        {
        {
        +200.0e0,         /* Tmin [K] */ 
        +1000.0e0,        /* Tmax [K] */
        -3.425563420e+04, /* a1 */
        +4.847000970e+02, /* a2 */
        +1.119010961e+00, /* a3 */
        +4.293889240e-03, /* a4 */
        -6.836300520e-07, /* a5 */
        -2.023372700e-09, /* a6 */
        +1.039040018e-12, /* a7 */
        -3.391454870e+03, /* b1 */
        +1.849699470e+01  /* b2 */
        },
        {
        +1000.0e0,        /* Tmin [K] */ 
        +6000.0e0,        /* Tmax [K] */
        -1.037939022e+06, /* a1 */
        +2.344830282e+03, /* a2 */
        +1.819732036e+00, /* a3 */
        +1.267847582e-03, /* a4 */
        -2.188067988e-07, /* a5 */
        +2.053719572e-11, /* a6 */
        -8.193467050e-16, /* a7 */
        -1.689010929e+04, /* b1 */
        +1.738716506e+01  /* b2 */        
        },
        {
        +6000.0e0,        /* Tmin [K] */ 
        +20000.0e0,       /* Tmax [K] */
        +4.975294300e+08, /* a1 */
        -2.866106874e+05, /* a2 */
        +6.690352250e+01, /* a3 */
        -6.169959020e-03, /* a4 */
        +3.016396027e-07, /* a5 */
        -7.421416600e-12, /* a6 */
        +7.278175770e-17, /* a7 */
        +2.293554027e+06, /* b1 */
        -5.530621610e+02  /* b2 */        
        }
        },


        // N2
        {
        {
        +200.0e0,         /* Tmin [K] */ 
        +1000.0e0,        /* Tmax [K] */
        +2.210371497e+04, /* a1 */
        -3.818461820e+02, /* a2 */
        +6.082738360e+00, /* a3 */
        -8.530914410e-03, /* a4 */
        +1.384646189e-05, /* a5 */
        -9.625793620e-09, /* a6 */
        +2.519705809e-12, /* a7 */
        +7.108460860e+02, /* b1 */
        -1.076003744e+01  /* b2 */
        },
        {
        +1000.0e0,        /* Tmin [K] */ 
        +6000.0e0,        /* Tmax [K] */
        +5.877124060e+05, /* a1 */
        -2.239249073e+03, /* a2 */
        +6.066949220e+00, /* a3 */
        -6.139685500e-04, /* a4 */
        +1.491806679e-07, /* a5 */
        -1.923105485e-11, /* a6 */
        +1.061954386e-15, /* a7 */
        +1.283210415e+04, /* b1 */
        -1.586640027e+01  /* b2 */        
        },
        {
        +6000.0e0,        /* Tmin [K] */ 
        +20000.0e0,       /* Tmax [K] */
        +8.310139160e+08, /* a1 */
        -6.420733540e+05, /* a2 */
        +2.020264635e+02, /* a3 */
        -3.065092046e-02, /* a4 */
        +2.486903333e-06, /* a5 */
        -9.705954110e-11, /* a6 */
        +1.437538881e-15, /* a7 */
        +4.938707040e+06, /* b1 */
        -1.672099740e+03  /* b2 */        
        }
        }

        

    };




double fac ;
fac = 1.0 ;


// NASA FORMAT but cp with 7 coefficients

if ( T <= 1000.0){

cp[0] = fac*(R0/MW_O2)*((C[0][0][2]/(T*T))+(C[0][0][3]/T)+C[0][0][4] + T* (C[0][0][5] + T*(C[0][0][6] + T*(C[0][0][7] + T*C[0][0][8]))));
cp[1] = fac*(R0/MW_N2)*((C[1][0][2]/(T*T))+(C[1][0][3]/T)+C[1][0][4] + T* (C[1][0][5] + T*(C[1][0][6] + T*(C[1][0][7] + T*C[1][0][8]))));

}
else if ( T <= 6000.0){

cp[0] = fac*(R0/MW_O2)*((C[0][1][2]/(T*T))+(C[0][1][3]/T)+C[0][1][4] + T* (C[0][1][5] + T*(C[0][1][6] + T*(C[0][1][7] + T*C[0][1][8]))));
cp[1] = fac*(R0/MW_N2)*((C[1][1][2]/(T*T))+(C[1][1][3]/T)+C[1][1][4] + T* (C[1][1][5] + T*(C[1][1][6] + T*(C[1][1][7] + T*C[1][1][8]))));

}else{

cp[0] = fac*(R0/MW_O2)*((C[0][2][2]/(T*T))+(C[0][2][3]/T)+C[0][2][4] + T* (C[0][2][5] + T*(C[0][2][6] + T*(C[0][2][7] + T*C[0][2][8]))));
cp[1] = fac*(R0/MW_N2)*((C[1][2][2]/(T*T))+(C[1][2][3]/T)+C[1][2][4] + T* (C[1][2][5] + T*(C[1][2][6] + T*(C[1][2][7] + T*C[1][2][8]))));



}

// printf("%f \n", R0/MW_O2);

int i , j ;

for (i=0 ; i<3 ; i++){

    // printf("\n %lf %lf %lf %lf %0.12f \n ", C[0][i][2], C[0][i][3], C[0][i][4], C[0][i][5], C[0][i][6]);
    // printf("\n %lf %lf %lf %lf %0.12f \n ", C[1][i][2], C[1][i][3], C[1][i][4], C[1][i][5], C[1][i][6]);
}



}

// From Capataile :: Air properties for ionized case
 
