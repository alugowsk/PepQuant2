/*
 * ls.h                                                                      
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
 * Header file for ls.c. Provides access to functions for calculating        
 *     alignment parameters as well as subsequent retention time alignmnet.  
 *                                                                           
 * +ASCII art via Text ASCII Art Generator by Patrick Gillespie              
 */

#ifndef LS_H
#define LS_H

/*
 * leastSquares - Given a retention time table and a vector of median        
 *     retention times fits a linear function for each median-run pair and
 *     return the regression parameters as a two dimensional array. 
 */
double **leastSquares(double **rt, double *median, int peptideCount,
	int fileCount);

/*
 * align - Given ms2 an ms1 retention time tables along with median retention
 *     time vectors and regressed parameters for linear functions align runs
 *     storing final results in the ms2 retention time table.
 */
void align(double **ms2rt, double **ms1rt, double *ms2median, 
	double *ms1median, double **ms2params, double **ms1params,
	int peptideCount, int fileCount);

#endif

