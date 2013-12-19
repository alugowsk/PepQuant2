/*
 * common.c                                                                  
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
 * Function definitions for routine or ubiquitous tasks.                     
 *                                                                           
 * +ASCII art via Text ASCII Art Generator by Patrick Gillespie              
 */

#include "common.h"
#include "global.h"

#include <stdlib.h> //qsort, malloc, free
#include <string.h> //memcpy
#include <math.h> //fabs, fmax
#include <stdio.h>

int comparedouble (const void * a, const void * b){
	return ( *(double*)b - *(double*)a );
}


double median(double *list, int size){
	
	/*make a copy and sort copy in place*/
	double copyList[size];
	memcpy(copyList, list, size*sizeof(double));
	qsort(copyList, size, sizeof(double), comparedouble); 
	
	/*remove zeroes if list does not contain negative values*/
	while(size > 0 && copyList[size-1] == 0){
		size--;
	}

	/* want at least three non-zero entries to take median */
	if (size < 3) return 0;

	/* IQR should be reasonable */
	double third_quartile = (size/2)%2 ?
			(copyList[size/4]) :
			( (copyList[size/4-1] + copyList[size/4]) / 2);
	double first_quartile = (size/2)%2 ?
				(copyList[3*size/4]) :
				( (copyList[3*size/4-1] + copyList[3*size/4]) / 2);
	double iqr = third_quartile - first_quartile;
	if(iqr > 250) {
		printf("%f, %f, %f\n", first_quartile,
				(size%2)?
						(copyList[size/2]):
						( (copyList[size/2-1] + copyList[size/2]) / 2 ),
				third_quartile);
		return 0;
	}



	if(size%2){
		return (copyList[size/2]);
	}else{
		return ( (copyList[size/2-1] + copyList[size/2]) / 2 );
	}
}


float pearson(float *y, float *foundPattern, float *corr){
	int i, j;
	float ySum = 0;
	for(i = 0; i < isotopicStates; ++i){
		ySum += y[i];
	}
	float yAvg = ySum / isotopicStates;

	float maxCorr = 0;
	for(i = 0; i < maxCharge - MIN_CHARGE + 1; ++i){
		float sum = 0;
		float x[isotopicStates];
		for(j = 0; j < isotopicStates; ++j){
			x[j] = foundPattern[i * isotopicStates + j];
			sum += x[j];
		}

		float xAvg = sum / isotopicStates;
		float top = 0;
		float bottomLeft = 0;
		float bottomRight = 0;
		for(j = 0; j < isotopicStates; ++j){
			top += (x[j] - xAvg) * (y[j] - yAvg);
			bottomLeft += (x[j] - xAvg) * (x[j] - xAvg);
			bottomRight	+= (y[j] - yAvg) * (y[j] - yAvg);
		}
		float bottom = sqrt(bottomLeft * bottomRight);
		corr[i] = bottom == 0? 0 : (top / bottom);
		maxCorr = maxCorr > corr[i]? maxCorr: corr[i];
	}
	return maxCorr;	
}


double **new2Darray(int x, int y){
	double **array = (double **)malloc(x * sizeof(double*));
	int i;
	for(i = 0; i < x; ++i){
		array[i] = (double *)malloc(y * sizeof(double));
	}
	return array;
}


double **del2Darray(double **array, int x){
	int i;
	for(i = 0; i < x; ++i){
		free(array[i]);
	}
	free(array);
	array = NULL;
	return array;
}

