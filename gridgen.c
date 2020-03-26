/*
 ============================================================================
 Name        : gridgen.c
 Author      : Sagar Pokharel
 Version     :
 Copyright   :
 Description : 2D Grid Generator : Algebraic
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "2Dsolver.h"
#include "gridgen.h"
#define PI 3.14159265358979323846

void grid(int IMAX , int JMAX ,double XL, double XR, double YL, double YU, double x[], double y[]){
/*
    Provide IMAX and JMAX to the function which are number of interfaces
    Number of interfaces are one more than the number of the cells
    Choose a clustering definition
    Write the cludtered grid specified values for x and y and pass x[IMAX] and y[JMAX]
    IMAX akin to n_interface_col or nx
    JMAX akin to n_interface_row or ny
*/


	int i,j;
    double zt[IMAX], et[JMAX] , val , m;
    FILE *fp , *fpar;

    m = 1.8 ; // more than one means more fine at the middle
		for (i=0; i < IMAX ; i++){

           // zt[i] = log(i + 1);
           // zt[i] = pow((i*1.0/(IMAX-1)),2.0);
           val = i*1.0/(IMAX-1);
           zt[i] = 3*(1-m) * val*val + 2*(m-1) * val*val*val + m * val; // fine at middle // <-- look for better function for this
            x[i] = XL + zt[i] * (XR - XL);
            printf("%f \t",x[i]);

        }
        printf("\n");
		for (j=0; j < JMAX ; j++){
			
           // et[j] = log(j + 1);
          // et[j] = pow((j*1.0/(JMAX-1)),2.0);
           // et[j] = 0.5*(1 - cos(j*1.0*PI/(JMAX-1))); // Fine at the walls
           
           val = j*1.0/(JMAX-1);
           et[j] = 3*(1-m) * val*val + 2*(m-1) * val*val*val + m * val; // fine at middle
            y[j] = YL + et[j] * (YU - YL);
            printf("%f \t",y[j]);
            

		}


// save nc_x and nc_col to file
        fpar = fopen("parameters.inp", "w+");
    //fprintf(fp,"Time =  %0.16f\n", time) ;
            fprintf(fpar,"%d \n%d\n",IMAX-1,JMAX-1);
            fclose(fpar);

    fp = fopen("mesh.txt", "w+");
    //fprintf(fp,"Time =  %0.16f\n", time) ;
    for (j=0; j < JMAX ; j++){
        for (i=0; i < IMAX ; i++){
            fprintf(fp,"%0.7f \t %0.7f\n",x[i],y[j]);
        }
   }


   fclose(fp);

   // save nc_x and nc_col to file
        fpar = fopen("output/parameters.inp", "w+");
    //fprintf(fp,"Time =  %0.16f\n", time) ;
            fprintf(fpar,"%d \n%d\n",IMAX-1,JMAX-1);
            fclose(fpar);

    fp = fopen("output/mesh.txt", "w+");
    //fprintf(fp,"Time =  %0.16f\n", time) ;
    for (j=0; j < JMAX ; j++){
        for (i=0; i < IMAX ; i++){
            fprintf(fp,"%0.7f \t %0.7f\n",x[i],y[j]);
        }
   }
fclose(fp);

} 