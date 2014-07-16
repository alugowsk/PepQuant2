/*
 * ls.c                                                                      
 * ====                                                                      
 *  ____                _____                            __       ___        
 * /\  _`\             /\  __`\                         /\ \__  /'___`\      
 * \ \ \L\ \ __   _____\ \ \/\ \  __  __     __      ___\ \ ,_\/\_\ /\ \     
 *  \ \ ,__/'__`\/\ '__`\ \ \ \ \/\ \/\ \  /'__`\  /' _ `\ \ \/\/_/// /__    
 *   \ \ \/\  __/\ \ \L\ \ \ \\'\\ \ \_\ \/\ \L\.\_/\ \/\ \ \ \_  // /_\ \   
 *    \ \_\ \____\\ \ ,__/\ \___\_\ \____/\ \__/.\_\ \_\ \_\ \__\/\______/+  
 *     \/_/\/____/ \ \ \/  \/__//_/\/___/  \/__/\/_/\/_/\/_/\/__/\/_____/     
 *                  \ \_\                                                    
 *                   \/_/                                                    
 *                                                                           
 *                              Andrew Lugowski                              
 *                   Emili Lab at the University of Toronto                  
 *                                                                           
 * Version 1.00 (Aug 22, 2013)                                              
 *                                                                           
 * Functions and structs for alinging retention times based on a linear 
 *     function whose parameters are determined by linear least squares. File
 *     depends heavily on MPFIT: A MINPACK-1 Least Squares Fitting Library in
 *     C.                                                         
 *                                                                           
 * +ASCII art via Text ASCII Art Generator by Patrick Gillespie              
 */

#include "ls.h"
#include "common.h" //new2Darray
#include "global.h" //alignWindow

#include "mpfit.h" //mpfit

#include <stdlib.h> //malloc, free
#include <math.h> //fabs

/*
 * vars_struct - This is the private data structure which contains the data   
 *     points and their uncertainties. Stolen from testmpfit.c of MPFIT: A   
 *     MINPACK-1 Least Squares Fitting Library in C.    
 */
struct vars_struct {
	double *x;
	double *y;
	double *ey;
};

/*
 * fit - Fit a linear function to the data. Based heavily on the func 
 *     testlinfit from testmpfit.c of MINPACK-1 Least Squares Fitting Library.
 */
int fit(double *x, double *y, int n, double *p);

/*
 * linfunc - modified version of linfunc from testmpfit.c of MINPACK-1 Least    *
 *     Squares Fitting Library. 
 */
int linfunc(int m, int n, double *p, double *dy, double **dvec, void *vars);

///////////////////////////////////////////////////////////////////////////////
//                       END OF FUNCTION DECLARATIONS                        //
///////////////////////////////////////////////////////////////////////////////

int fit(double *x, double *y, int n, double *p){
 	int i;

	double *ey = (double *)malloc(n * sizeof(double));
	for (i=0; i<n; i++){
		ey[i] = 1;
	}  

	int status;
	struct vars_struct v;
	v.x = x;
	v.y = y;
	v.ey = ey;

	/* Call fitting function for n data points */
	status = mpfit(linfunc, n, 2, p, 0, 0, (void *) &v, 0);
 	return status;
}


int linfunc(int m, int n, double *p, double *dy, double **dvec, void *vars){
	int i;
	struct vars_struct *v = (struct vars_struct *) vars;
	double *x, *y, *ey, f;

	x = v->x;
	y = v->y;
	ey = v->ey;

	for (i=0; i<m; i++) {
		f = p[0] - p[1]*x[i];     /* Linear fit function; note f = a - b*x */
		dy[i] = (y[i] - f)/ey[i];
	}
  return 0;
}


double **leastSquares(double **rt, double *median, int peptideCount,
	int fileCount){
	double **params = new2Darray(fileCount, 2);
	int i,j;
	for(i=0; i<fileCount; ++i){
		double *x = (double*)malloc(peptideCount*sizeof(double) );
		double *y = (double*)malloc(peptideCount*sizeof(double) );
		for(j=0; j<2; ++j){
			params[i][j] = 1;
		}		
		int n = 0;
		/*Get run rt and median rt when run rt is not zero*/ 
		for(j=0; j<peptideCount; ++j){
			if(rt[j][i] != 0 && median[j] != 0
					// 2014-07-16 added below
					// check to align on well behaved peptides only
					&& fabs(rt[j][i] - median[j]) < alignWindow){
				y[n] = rt[j][i];
				x[n] = median[j];
				n++;
			}
		}
		fit(x, y, n, params[i]);
		free(x);
		free(y);
	}
	return params;
}


void align(double **ms2rt, double **ms1rt, double *ms2median, 
	double *ms1median, double **ms2params, double **ms1params,
	int peptideCount, int fileCount){

	int i, j;	

	for(i = 0; i < peptideCount; ++i){
		for(j = 0; j < fileCount; ++j){
			/*if there is no ms2 identification*/
			if(ms2rt[i][j] == 0){
				/*if median ms1 and ms2 retention times are close*/
				if(fabs(ms2median[i] - ms1median[i]) < alignWindow && 
					ms1median[i] > 0){

					ms2rt[i][j] = ms1params[j][0] - 
						ms1params[j][1]*ms1median[i];
				}else{
					ms2rt[i][j] = ms2params[j][0] -
						ms2params[j][1]*ms2median[i];
				}
			/*use ms1 over ms2 only if there are close*/
			}else if(fabs(ms2rt[i][j] - ms1rt[i][j]) < alignWindow &&
				ms1rt[i][j] > 0){

				ms2rt[i][j] = ms1rt[i][j];
			}
		}
	}
	return;
}

