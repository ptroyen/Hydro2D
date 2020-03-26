#ifndef THERMO_H_
#define THERMO_H_

double eval_cp( double T);

void read_cp(double T , double cp[]);

void read_cp_i(double T , double cp[2]);

double iter_tempr(double e, double tempr_old);

void thermo_prmcalculate(int nc_row, int nc_col,double q[][nc_row][nc_col],double gamma,
                     double r[][nc_col],double u[][nc_col],double v[][nc_col],double p[][nc_col],double T[][nc_col],double s[][nc_col]);

void thermo_sing_prmcalculate(double q[],double guess_gamma, double *r,double *u,double *v,double *p, double *T);

void thermo_concalculate(int nc_row, int nc_col, double r[][nc_col],double u[][nc_col],double v[][nc_col],double p[][nc_col],double q[][nc_row][nc_col]);

double thermo_CFLmaintain(int nc_row, int nc_col, double r[][nc_col],double u[][nc_col],double v[][nc_col],
                    double p[][nc_col],double CFL,double x[nc_col+1],double y[nc_row+1],int n);

void thermo_x_fun_flux(double qre[],double guess_gamma,double ff[]);

void thermo_y_fun_flux(double qre[],double guess_gamma,double ff[]);

void thermo_RSflux(int nc_row , int nc_col , double Qil[4][nc_row][nc_col],double Qir[4][nc_row][nc_col],double guess_gamma[][nc_col],
                     char DIR ,int f_row , int f_col, double Flux[4][f_row][f_col]);


double eval_gamma_pr(double p,double r);

double eval_gamma_t(double tempr);

void thermo_Flux_M(int nrec_row , int nrec_col, double qre[4][nrec_row][nrec_col] ,double pre_for_qre[][nrec_col],double xre[], double yre[] ,
                    double dt,double guess_gamma,  double accu[4][nrec_row][nrec_col]);// ## add what is needed


#endif /* THERMO_H_ */