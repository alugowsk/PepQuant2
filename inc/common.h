/*
 * common.h                                                                  
 * ========                                                                  
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
 * Header file for common.c. Contains the declartions of functions for       
 *     routine or ubiquitous tasks and common constants.                      
 *                                                                           
 * +ASCII art via Text ASCII Art Generator by Patrick Gillespie              
 */

#ifndef COMMON_H
#define COMMON_H

#define MIN_CHARGE 1

/*
 * comparedouble - Return 1 if a is larger than b, -1 if b is larger than a, 
 *     or 0 if they are equal.                                   
 */
int comparedouble (const void * a, const void * b);

/*
 * median - Return the median value in a list of doubles with size elements. 
 *     If the list is comprised only of numbers equal to or greater than 0 the
 *     median will be of the non-zero elements only.         
 */
double median(double *list, int size);

/*
 * pearson - Determine the pearson correlation of list y of size STATES against
 *     a collection of equally sized list x. Store the correlations in corr and
 *     return the max correlation.            
 */
float pearson(float *y, float *x, float *corr);

/*
 * new2Darray - Allocate memory for a two dimensional double array.
 */
double **new2Darray(int x, int y);

/*
 * del2Darray - Free memory allocated for a two dimensional double array. The
 *     passed integer specifies the first dimension of the array.
 */
double **del2Darray(double **array, int x);

#endif

