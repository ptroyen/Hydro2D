#ifndef SOLVER_H_
#define SOLVER_H_

void temperature(int nc_row , int nc_col, double q[][nc_row][nc_col],double c_v, double t[][nc_col]);

void concalculate(int nc_row, int nc_col, double r[][nc_col],double u[][nc_col],double v[][nc_col],double p[][nc_col],double gamma,double q[][nc_row][nc_col]);

double superbee(double r,double gl,double gr);

double vanLeer(double r,double gl,double gr);

double vanAlbada(double r,double gl,double gr);

double CFLmaintain(int nc_row, int nc_col, double r[][nc_col],double u[][nc_col],double v[][nc_col],double p[][nc_col],double gamma,double CFL,double x[nc_col+1],double y[nc_row+1],int n);

void prmcalculate(int nc_row, int nc_col,double q[][nc_row][nc_col],double gamma, double r[][nc_col],double u[][nc_col],double v[][nc_col],double p[][nc_col],double s[][nc_col]);

void sing_prmcalculate(double q[],double gamma, double *r,double *u,double *v,double *p);

void x_fun_flux(double qre[],double gamma,double ff[]);

void y_fun_flux(double qre[],double gamma,double ff[]);

void ESTIME(double DL,double UL,double PL,double DR,double UR,double PR,double GAMMA, double *SL , double *SM , double *SR);

void Flux_M(int nrec_row , int nrec_col , double qre[4][nrec_row][nrec_col] ,double xre[], double yre[] ,double dt,double gamma,  double accu[4][nrec_row][nrec_col]);

void HLLC(int nc_row , int nc_col , double Qil[4][nc_row][nc_col],double Qir[4][nc_row][nc_col],double GAMMA, char DIR ,int f_row , int f_col, double Flux[4][f_row][f_col]);

void RSflux(int nc_row , int nc_col , double Qil[4][nc_row][nc_col],double Qir[4][nc_row][nc_col],double GAMMA, char DIR ,int f_row , int f_col, double Flux[4][f_row][f_col]);

void Source(int nc_row, int nc_col,double q[4][nc_row][nc_col], double x[], double y[], double t,double dt, double c_v, double source_accu[nc_row][nc_col]);

#endif /* 2DSOLVER_H_ */
