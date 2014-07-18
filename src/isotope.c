/*
 * isotope.c                                                                 
 * =========                                                                 
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
 * Functions for creating, merging, and deleting isotopic patterns and       
 *     collections of isotopic patterns. Also contains a modified merge sort 
 *     algorithm for sorting intensities and masses of isotopes when merging 
 *     two existing isotopic profiles.                                       
 *                                                                           
 * +ASCII art via Text ASCII Art Generator by Patrick Gillespie              
 */

#include "isotope.h"
#include "global.h" //ppmCutOff, isotopicStates

#include <stdio.h> //printf, fprintf
#include <stdlib.h> //malloc, free
#include <string.h> //memcpy
#include <math.h> //fabs
#include <ctype.h> //isupper

/*Elemental data, BLANK for initializing IsotopicPattern structs*/
//http://www.chem.ualberta.ca/~massspec/atomic_mass_abund.pdf
float BLANK[ISOTOPE_COUNT] = {1}; 
float BLANK_MASS[ISOTOPE_COUNT] = {0};
float C[ISOTOPE_COUNT] = {.9893, 0.0107, 0};
float C_MASS[ISOTOPE_COUNT] = {12.000000, 13.003355, 0};
float H[ISOTOPE_COUNT] = {.999885, .000115, 0};
float H_MASS[ISOTOPE_COUNT] = {1.007825, 2.014102, 0};
float N[ISOTOPE_COUNT] = {.99632, .00368, 0};
float N_MASS[ISOTOPE_COUNT] = {14.003074, 15.000109, 0};
float O[ISOTOPE_COUNT] = {.99757, .00038, 0.00205, 0};
float O_MASS[ISOTOPE_COUNT] = {15.994915, 16.999132, 17.999160, 0};
float P[ISOTOPE_COUNT] = {1.0, 0};
float P_MASS[ISOTOPE_COUNT] = {30.973762, 0};
float S[ISOTOPE_COUNT] = {.9493, .0076, .0429, 0.0002, 0};
float S_MASS[ISOTOPE_COUNT] = {31.972071, 32.971458, 33.967867, 35.967081, 0};
float Se[ISOTOPE_COUNT] = {.0089, .0937, .0763, .2377, .4961, .0873, .0};
float Se_MASS[ISOTOPE_COUNT] = {77.920386, 79.916378, 81.913485, 82.914136,
								83.911507, 85.910610, 0 };

float C13[ISOTOPE_COUNT] = {1, 0}; //heavy
float C13_MASS[ISOTOPE_COUNT] = {13.003355, 0}; //heavy
float N15[ISOTOPE_COUNT] = {1, 0}; //heavy
float N15_MASS[ISOTOPE_COUNT] = {15.000109, 0}; //heavy

/*
 * newIsotopicPattern - Allocate memory for a new IsotopicPattern and init with
 *     passed values. Return a pointer to new IsotopicPattern, NULL if error
 *     occured.  
 */
IsotopicPatternPointer newIsotopicPattern(float *intensity, float *mass);

/*
 * combineIsotopicPatterns - Creates a new IsotopicPattern by combining two 
 *     existing IsotopicPatterns. Note that user is responsible for deleting
 *     old IsotopicPatterns if they are no longer required.
 */
IsotopicPatternPointer combineIsotopicPatterns(IsotopicPatternPointer ipp1,
	IsotopicPatternPointer ipp2);

/*
 * makeMolecule - Create an IsotopicPattern for a molecule comprised of the 
 *     elements CHNOPS by specifying the amount of each. Return a pointer to
 *     the IsotopicPattern. Also supports C13 and N15. 
 */
IsotopicPatternPointer makeMolecule(int CCount, int HCount, int NCount,
	int OCount, int PCount, int SCount, int C13Count, int N15Count);

/*
 * merge_sort - Helper function for the slight modifcation of classic merge 
 *     sort algorithm. Sorts intensity and mass array in place.
 */
void merge_sort(float *intensity, float *mass, int left, int right);

/*
 * merge - A slight variation of the classic merge sort algorithm. Merges
 *     two halves of an array in place, doing so for two arrays simultaneously.
 *     Any rearrangement done to one must be reflected in the other. An
 *     additional modification is that if two entries share a similar mass the
 *     entries are merged.
 */
void merge(float *intensity, float *mass, int left, int right);

/*
 * swapEntries - Swap entries i and j in lists intensity and mass.
 */
void swapEntries(float *intensity, float *mass, int i, int j);


///////////////////////////////////////////////////////////////////////////////
//                       END OF FUNCTION DECLARATIONS                        //
///////////////////////////////////////////////////////////////////////////////

IsotopicPatternPointer newIsotopicPattern(float *intensity, float *mass){

	IsotopicPatternPointer ipp = (IsotopicPatternPointer)malloc(
		sizeof(IsotopicPattern));

	if(ipp == NULL){
		fprintf(stderr,
			"\nERROR: Out of memory - cannot create isotopic pattern.\n");
	}else{
		ipp->intensity = (float*)malloc(sizeof(float)*ISOTOPE_COUNT);
		ipp->mass = (float*)malloc(sizeof(float)*ISOTOPE_COUNT);
		if(ipp->intensity == NULL || ipp->mass == NULL){
		fprintf(stderr,
			"\nERROR: Out of memory - cannot create isotopic pattern.\n");
		}else{
			memcpy(ipp->intensity, intensity, ISOTOPE_COUNT*sizeof(float));
			memcpy(ipp->mass, mass, ISOTOPE_COUNT*sizeof(float));
		}
	}
	return ipp;
}


IsotopicPatternPointer combineIsotopicPatterns(IsotopicPatternPointer ipp1,
	IsotopicPatternPointer ipp2){

	if(ipp1 == NULL || ipp2 == NULL){
		printf("one of the isotopic patterns was NULL!\n"); 
	}

	int i, j;

	float *combinedIntensity = (float *)malloc(
		ISOTOPE_COUNT*ISOTOPE_COUNT*sizeof(float));
	float *combinedMass = (float *)malloc(
		ISOTOPE_COUNT*ISOTOPE_COUNT*sizeof(float));

	for(i = 0; i < ISOTOPE_COUNT; ++i){
		for(j = 0; j < ISOTOPE_COUNT; ++j){
			combinedIntensity[i*ISOTOPE_COUNT + j] = 
				ipp1->intensity[i] * ipp2->intensity[j];
			combinedMass[i*ISOTOPE_COUNT + j] =
				combinedIntensity[i*ISOTOPE_COUNT + j] == 0? 0:
				ipp1->mass[i] + ipp2->mass[j];
		}
	}
	merge_sort(combinedIntensity, combinedMass,
		0, ISOTOPE_COUNT*ISOTOPE_COUNT-1);
	
	
	IsotopicPatternPointer ipp = newIsotopicPattern(combinedIntensity,
		combinedMass);

	if(combinedIntensity != NULL){
		free(combinedIntensity);	
	}
	if(combinedMass != NULL){
		free(combinedMass);
	}

	return ipp;
}


IsotopicPatternPointer makeMolecule(int CCount, int HCount, int NCount,
	int OCount, int PCount, int SCount, int C13Count, int N15Count){

	IsotopicPatternPointer carbon = newIsotopicPattern(C, C_MASS);
	IsotopicPatternPointer hydrogen = newIsotopicPattern(H, H_MASS);
	IsotopicPatternPointer nitrogen = newIsotopicPattern(N, N_MASS);
	IsotopicPatternPointer oxygen = newIsotopicPattern(O, O_MASS);
	IsotopicPatternPointer phosphor = newIsotopicPattern(P, P_MASS);
	IsotopicPatternPointer sulfur = newIsotopicPattern(S, S_MASS);
	IsotopicPatternPointer heavyC = newIsotopicPattern(C13, C13_MASS); //heavy
	IsotopicPatternPointer heavyN = newIsotopicPattern(N15, N15_MASS); //heavy

	IsotopicPatternPointer molecule = newIsotopicPattern(BLANK, BLANK_MASS);
	int i;
	for(i = 0; i < CCount; ++i){
		IsotopicPatternPointer temp =
			combineIsotopicPatterns(molecule, carbon);
		delIsotopicPattern(molecule);
		molecule = temp;
	}
	for(i = 0; i < HCount; ++i){
		IsotopicPatternPointer temp =
			combineIsotopicPatterns(molecule, hydrogen);
		delIsotopicPattern(molecule);
		molecule = temp;
	}
	for(i = 0; i < NCount; ++i){
		IsotopicPatternPointer temp =
			combineIsotopicPatterns(molecule, nitrogen);
		delIsotopicPattern(molecule);
		molecule = temp;
	}
	for(i = 0; i < OCount; ++i){
		IsotopicPatternPointer temp =
			combineIsotopicPatterns(molecule, oxygen);
		delIsotopicPattern(molecule);
		molecule = temp;
	}
	for(i = 0; i < PCount; ++i){
		IsotopicPatternPointer temp =
			combineIsotopicPatterns(molecule, phosphor);
		delIsotopicPattern(molecule);
		molecule = temp;
	}
	for(i = 0; i < SCount; ++i){
		IsotopicPatternPointer temp =
			combineIsotopicPatterns(molecule, sulfur);
		delIsotopicPattern(molecule);
		molecule = temp;
	}
	for(i = 0; i < C13Count; ++i){
		IsotopicPatternPointer temp =
			combineIsotopicPatterns(molecule, heavyC);
		delIsotopicPattern(molecule);
		molecule = temp;
	}
	for(i = 0; i < N15Count; ++i){
		IsotopicPatternPointer temp =
			combineIsotopicPatterns(molecule, heavyN);
		delIsotopicPattern(molecule);
		molecule = temp;
	}

	delIsotopicPattern(carbon);
	delIsotopicPattern(hydrogen);
	delIsotopicPattern(nitrogen);
	delIsotopicPattern(oxygen);
	delIsotopicPattern(phosphor);
	delIsotopicPattern(sulfur);
	delIsotopicPattern(heavyC);
	delIsotopicPattern(heavyN);

	return molecule;
}



void merge_sort(float *intensity, float *mass, int left, int right){
	if(abs(right-left) <= 1){
		return;
	}else{
		merge_sort(intensity, mass, left, (left+right)/2);
		merge_sort(intensity, mass, (left+right)/2 +1, right);
	}
	merge(intensity, mass, left, right);
}


void merge(float *intensity, float *mass, int left, int right){

	int lmax = (left+right)/2;
	int rmin = lmax+1;
	int i = left;
	int j = rmin;

	/*Skip any 0 values on right side*/
	while(intensity[right] == 0){
		right--;
	}
	
	while( i <= lmax && j <= right){
		if(fabs(mass[i] - mass[j])/mass[i] < ppmCutOff){
			intensity[i] += intensity[j];
			
			mass[j] = mass[right];
			intensity[j] = intensity[right];

			mass[right] = 0;
			intensity[right] = 0;
			right--;
			
			int k;
			/*Make sure that right side maintains order*/
			for(k = right; k > j; --k){
				if(intensity[j] < intensity[k]){
					swapEntries(intensity, mass, k, j);
				}
			}
			/*Make sure that left side maintains order*/
			for(k = left; k < i; ++k){
				if(intensity[i] > intensity[k]){
					swapEntries(intensity, mass, k, i);
				}
			}
			/*reset i and j in case either left or right side were altered*/
			i = left;
			j = rmin;
		}else if(intensity[i] < intensity[j]){
			swapEntries(intensity, mass, i,  j);
			j++;
		}else{
			i++;
		}
	}
}


void swapEntries(float *intensity, float *mass, int i, int j){
	float tempInt = intensity[i];
	float tempMass = mass[i];
	intensity[i] = intensity[j];
	mass[i] = mass[j];
	intensity[j] = tempInt;
	mass[j] = tempMass;
}


IsotopicPatternPointer delIsotopicPattern(IsotopicPatternPointer ipp){

	if(ipp == NULL){
		return NULL;
	}
	if(ipp->intensity != NULL){
		free(ipp->intensity);
	}
	if(ipp->mass != NULL){
		free(ipp->mass);
	}
	free(ipp);
	ipp = NULL;
	return ipp;
}


void printIsotopicPattern(IsotopicPatternPointer ipp){
	int i;
	for(i = 0; i < ISOTOPE_COUNT; ++i){
		printf("%f\t%f\n", ipp->intensity[i], ipp->mass[i]);
	}	
}


void newIPCollection(IsotopicPatternPointer **IPCollection){

	int i;
	for(i = 0; i < COLLECTION_SIZE; ++i){
		(*IPCollection)[i] = NULL;
	}

	/*minus a water to allow for condensation reactions*/
	(*IPCollection)['A'] = makeMolecule(3, 5, 1, 1, 0, 0, 0, 0); //Ala A
	(*IPCollection)['R'] = makeMolecule(6, 12, 4, 1, 0, 0, 0, 0); //Arg R
	(*IPCollection)['N'] = makeMolecule(4, 6, 2, 2, 0, 0, 0, 0); //Asn N
	(*IPCollection)['D'] = makeMolecule(4, 5, 1, 3, 0, 0, 0, 0); //Asp D
	//(*IPCollection)['C'] = makeMolecule(3, 5, 1, 1, 0, 1, 0, 0); //Cys C
	(*IPCollection)['C'] = makeMolecule(5, 8, 2, 2, 0, 1, 0, 0);
	//Cys C + H3C2NO
	(*IPCollection)['Q'] = makeMolecule(5, 8, 2, 2, 0, 0, 0, 0); //Gln Q
	(*IPCollection)['E'] = makeMolecule(5, 7, 1, 3, 0, 0, 0, 0); //Glu E
	(*IPCollection)['G'] = makeMolecule(2, 3, 1, 1, 0, 0, 0, 0); //Gly G
	(*IPCollection)['H'] = makeMolecule(6, 7, 3, 1, 0, 0, 0, 0); //His H
	(*IPCollection)['I'] = makeMolecule(6, 11, 1, 1, 0, 0, 0, 0); //Ile I
	(*IPCollection)['L'] = makeMolecule(6, 11, 1, 1, 0, 0, 0, 0); //Leu L
	(*IPCollection)['K'] = makeMolecule(6, 12, 2, 1, 0, 0, 0, 0); //Lys K
	(*IPCollection)['M'] = makeMolecule(5, 9, 1, 1, 0, 1, 0, 0); //Met M
	(*IPCollection)['F'] = makeMolecule(9, 9, 1, 1, 0, 0, 0, 0); //Phe F
	(*IPCollection)['P'] = makeMolecule(5, 7, 1, 1, 0, 0, 0, 0); //Pro P
	(*IPCollection)['S'] = makeMolecule(3, 5, 1, 2, 0, 0, 0, 0); //Ser S
	(*IPCollection)['T'] = makeMolecule(4, 7, 1, 2, 0, 0, 0, 0); //Thr T
	(*IPCollection)['W'] = makeMolecule(11, 10, 2, 1, 0, 0, 0, 0); //Trp W
	(*IPCollection)['Y'] = makeMolecule(9, 9, 1, 2, 0, 0, 0, 0); //Tyr Y
	(*IPCollection)['V'] = makeMolecule(5, 9, 1, 1, 0, 0, 0, 0); //Val V

	//make selenocysteine
	IsotopicPatternPointer selenocysteine =
		makeMolecule(3, 5, 1, 1, 0, 0, 0, 0); //Sec U
	IsotopicPatternPointer selenium = newIsotopicPattern(Se, Se_MASS);
	IsotopicPatternPointer temp =
		combineIsotopicPatterns(selenocysteine, selenium);
	delIsotopicPattern(selenium);
	delIsotopicPattern(selenocysteine);
	(*IPCollection)['U'] = temp;

	//make heavy AA
	(*IPCollection)['r'+10] = makeMolecule(0, 12, 0, 1, 0, 0, 6, 4); //arg-10
	(*IPCollection)['r'+6] = makeMolecule(0, 12, 4, 1, 0, 0, 6, 0); //arg-6
	(*IPCollection)['r'+4] = makeMolecule(6, 12, 0, 1, 0, 0, 0, 4); //arg-4
	(*IPCollection)['k'-8] = makeMolecule(0, 12, 0, 1, 0, 0, 6, 2); //lys-8
	(*IPCollection)['k'-6] = makeMolecule(0, 12, 2, 1, 0, 0, 6, 0); //lys-6
	(*IPCollection)['k'-2] = makeMolecule(6, 12, 0, 1, 0, 0, 0, 2); //lys-2	

	//make mods and other
	(*IPCollection)[0] = makeMolecule(0, 2, 0, 1, 0, 0, 0, 0); //H2O
	(*IPCollection)[1] = makeMolecule(0, 1, 0, 3, 1, 0, 0, 0); //P03H (ph)
	(*IPCollection)[2] = makeMolecule(0, 0, 0, 1, 0, 0, 0, 0); //O (ox)
	(*IPCollection)[3] = makeMolecule(2, 2, 0, 1, 0, 0, 0, 0); //C2H2O (ac)

	return;
}


IsotopicPatternPointer* delIPCollection(IsotopicPatternPointer *IPCollection){
	int i;
	if(IPCollection == NULL){
		return NULL;
	}
	for(i = 0; i < COLLECTION_SIZE; ++i){
		if(IPCollection[i] != NULL){
			delIsotopicPattern(IPCollection[i]);
			IPCollection[i] = NULL;
		}
	}
	free(IPCollection);
	IPCollection = NULL;
	return IPCollection;
}


IsotopicPatternPointer makePeptide(
	IsotopicPatternPointer *IPCollection, char *sequence){

	int length = strlen(sequence);
	int i;
	int isHeavy = sequence[0] == '*'? 1: 0;
	
	IsotopicPatternPointer peptide = newIsotopicPattern(BLANK, BLANK_MASS);
	for(i = 0; i < length; ++i){ //might actually loop less than length times
		if(sequence[i] == '_' || sequence[i] == '.' ||  sequence[i] == '*'){
			;
		}else if(sequence[i] == '('){
			i++;
			IsotopicPatternPointer temp = NULL;
			if(sequence[i] =='a'){
				temp = combineIsotopicPatterns(peptide, IPCollection[3]);
			}else if(sequence[i] =='o'){
				temp = combineIsotopicPatterns(peptide, IPCollection[2]);
			}else if(sequence[i] =='p'){
				temp = combineIsotopicPatterns(peptide, IPCollection[1]);
			}else if(sequence[i] =='c'){
				// do nothing C is carbidomethylated by default
			}else{
				printf("Unknown mod in seq: %s\n", sequence);
			}
			delIsotopicPattern(peptide);
			peptide = temp;
			i+=2; //to skip rest of mod
		}else if(strchr(AMINO_ACIDS, sequence[i])){
			IsotopicPatternPointer temp = NULL;
			if(isHeavy && lys && sequence[i] == 'K'){
				temp = combineIsotopicPatterns(peptide, 
					IPCollection[(int)'k' - lys]); 				
			}else if(isHeavy && arg && sequence[i] == 'R'){
				temp = combineIsotopicPatterns(peptide, 
					IPCollection[(int)'r' + arg]); 
			}else{
				temp = combineIsotopicPatterns(peptide, 
					IPCollection[(int)sequence[i]]); 
			}
			delIsotopicPattern(peptide);
			peptide = temp;
		}else{
			fprintf(stderr, "\nERROR: Unknown char[%c] in sequence: %s.\n",
				sequence[i], sequence);
			return NULL;
		}
	}
	/*compensate for missing water at peptide ends*/
	IsotopicPatternPointer temp =
			combineIsotopicPatterns(peptide, IPCollection[0]);
	delIsotopicPattern(peptide);
	peptide = temp;

	/*will only be looking at top 'isotopicStates' so why not save some mem*/
	peptide->intensity = 
		realloc(peptide->intensity, isotopicStates*sizeof(float) );
	peptide->mass = 
		realloc(peptide->mass, isotopicStates*sizeof(float) );
	
	return peptide;
}

