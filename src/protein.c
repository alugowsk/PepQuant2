/*
 * protein.c                                                                 
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
 * Functions for creating, deleting and manipulating protein lists.          
 *                                                                           
 * +ASCII art via Text ASCII Art Generator by Patrick Gillespie              
 */

#include "protein.h"
#include "global.h"

#include <stdio.h> //fprintf
#include <string.h> //strlen, strcmp, strchr, strrchr
#include <stdlib.h> //malloc, free

/*
 * newProteinNode - Allocate memory for a new protein and intialize with  
 *     passed values. Return a pointer to new protein, NULL if error occured.
 */
ProteinNodePointer newProteinNode(char *id, double *intensities, 
	double *spectralCounts, int fileCount, int isHeavy);

/*
 * delProteinNode - Free all memory allocated for the ProteinNode. Return 
 *     pointer to next proteinNode in linked list if exists, otherwise return
 *     NULL.
 */
ProteinNodePointer delProteinNode(ProteinNodePointer pnp);

/*
 * addProteinNode - Add a new protein to the linked list of proteins. Proteins
 *     will be added in decreasing order (according to strcmp).
 */
ProteinNodePointer addProteinNode(ProteinNodePointer root, char *id,
	double *intensities, double *spectralCounts, int fileCount, int isHeavy);

///////////////////////////////////////////////////////////////////////////////
//                       END OF FUNCTION DECLARATIONS                        //
///////////////////////////////////////////////////////////////////////////////

ProteinNodePointer newProteinNode(char *id, double *intensities, 
	double *spectralCounts, int fileCount, int isHeavy){

	ProteinNodePointer pnp = (ProteinNodePointer)malloc(sizeof(ProteinNode));
	if(pnp == NULL){
		fprintf(stderr,
			"\nERROR: Out of memory - cannot create proteinNode!\n");
	}else{
		pnp->id = (char *)malloc( (strlen(id)+1)*sizeof(char) );
		pnp->intensities = (double *)malloc( fileCount * sizeof(double) );
		pnp->spectralCounts = (double *)malloc( fileCount * sizeof(double) );
		pnp->light = (double *)malloc( fileCount * sizeof(double) );
		pnp->heavy = (double *)malloc( fileCount * sizeof(double) );
		if(!pnp->id || !pnp->intensities || !pnp->spectralCounts ||
			!pnp->light || !pnp->heavy){
			fprintf(stderr,
				"\nERROR: Out of memory - cannot create proteinNode!\n");
		}else{
			memcpy(pnp->spectralCounts, spectralCounts,
				fileCount * sizeof(double));
			strncpy(pnp->id, id, strlen(id));
			pnp->id[strlen(id)] = '\0';
			int i;
			for(i = 0; i < fileCount; ++i){
				pnp->intensities[i] = 0;				
				pnp->light[i] = 0;
				pnp->heavy[i] = 0;
			}
			double *dest = NULL;
 			if(isHeavy < 0){
				dest = pnp->intensities;
			}else if(isHeavy == 0){
				dest = pnp->light;
			}else{
				dest = pnp->heavy;
			}
			memcpy(dest, intensities, fileCount * sizeof(double) );
			pnp->next = NULL;
		}
	}
	return pnp;
}


ProteinNodePointer delProteinNode(ProteinNodePointer pnp){
	if(pnp == NULL){
		return NULL;
	}
	if(pnp->id != NULL){
		free(pnp->id);
	}
	if(pnp->intensities != NULL){
		free(pnp->intensities);
	}
	if(pnp->spectralCounts != NULL){
		free(pnp->spectralCounts);
	}
	if(pnp->light != NULL){
		free(pnp->light);
	}
	if(pnp->heavy != NULL){
		free(pnp->heavy);
	}
	ProteinNodePointer next = pnp->next;
	pnp->next = NULL;
	free(pnp);
	pnp = NULL;
	return next;
}


ProteinNodePointer addProteinNode(ProteinNodePointer root, char *id,
	double *intensities, double *spectralCounts, int fileCount, int isHeavy){

	if(root == NULL){
		return newProteinNode(id, intensities, spectralCounts, fileCount,
			isHeavy);
	}else if(!strcmp(root->id, id) ){
		int i;
		for(i = 0; i < fileCount; ++i){
			root->spectralCounts[i] += spectralCounts[i];
		}
		if(isHeavy < 0){
			for(i = 0; i < fileCount; ++i){
				root->intensities[i] += intensities[i];
			}
		}else if(isHeavy == 0){
			for(i = 0; i < fileCount; ++i){
				root->light[i] += intensities[i];
			}
		}else{
			for(i = 0; i < fileCount; ++i){
				root->heavy[i] += intensities[i];
			}
		}
		return root;
	}else if(strcmp(root->id, id) < 0){
		ProteinNodePointer pnp = newProteinNode(id, intensities, spectralCounts,
			fileCount, isHeavy);
		pnp->next = root;
		return pnp;
	}else{
		root->next = addProteinNode(root->next, id, intensities, spectralCounts,
			fileCount, isHeavy);
		return root;
	}
}


int pep2prot(ProteinNodePointer *root, FastaPointer fasta,
	PeptidePointer *peptides, double **pepQuant, double **spectralCounts,
	int peptideCount, int fileCount){

	ProteinNodePointer pnp = NULL;
	int proteinCount = 0;

	int i;
	for(i = 0; i < peptideCount; ++i){
		int isHeavy = 0;
		if(!(arg || lys)){ //not looking for heavy stuff
			isHeavy = -1;
		}else if(peptides[i]->sequence[0] == '*'){ //is heavy
			isHeavy = 1;
		}else if(!strchr(peptides[i]->sequence, 'K') && 
			!strchr(peptides[i]->sequence, 'R')){ //not K or R
			isHeavy = -1;
		}else{ //we are looking for heavy stuff but this is light
			isHeavy = 0;
		}
		char *stripped = stripMods(peptides[i]->sequence);
		char *prots = findProtein(fasta, stripped);		
		/*check to see if only one ; in returned string, i.e. only one protein
		  in ';' seperated list.*/
		if( prots != NULL && strchr(prots, ';') == strrchr(prots, ';') ){
			trackCoverage(fasta, stripped);
			prots[strlen(prots)-1] = '\0'; //remove ';' from protein name
			pnp = addProteinNode(pnp, prots, pepQuant[i], spectralCounts[i],
				fileCount, isHeavy);
			free(prots);
		}
		free(stripped);
	}
	
	(*root) = pnp;
	while(pnp != NULL){
		proteinCount++;
		pnp = pnp->next;
	}

	return proteinCount;
}


void prot2table(ProteinNodePointer pnp, double ***protQuant,
	double ***protHQuant, double ***protLQuant, double ***protSpectralCounts,
	int proteinCount, int fileCount){

	int i, j;

	for(i = 0; i < proteinCount; ++i){
		memcpy( (*protSpectralCounts)[i], pnp->spectralCounts,
			fileCount*sizeof(double) );
		memcpy( (*protLQuant)[i], pnp->light, fileCount*sizeof(double) );
		memcpy( (*protHQuant)[i], pnp->heavy, fileCount*sizeof(double) );
		for(j = 0; j < fileCount; ++j){
			(*protQuant)[i][j] = 
				pnp->intensities[j] + pnp->light[j] + pnp->heavy[j];
		}
		pnp = pnp->next;
	}
	return;
}


ProteinNodePointer delProteinList(ProteinNodePointer pnp){
	if(pnp == NULL){
		return NULL;
	}
	while(pnp != NULL){
		pnp = delProteinNode(pnp);
	}
	return pnp;
}


void printProtTable(ProteinNodePointer pnp, SpectraFileNodePointer filelist,
	int proteinCount, int fileCount, double **table, char *filename){

	FILE *fp = fopen(filename, "w");
	if (fp == NULL)	{
   		fprintf(stderr, "Error opening file: %s !\n", filename);
	}

	while(filelist != NULL){
		fprintf(fp, "\t%s", filelist->rawFile);
		filelist = filelist->next;
	}
	fprintf(fp, "\n");

	int i;
	int j;
	for(i = 0; i < proteinCount; ++i){
		fprintf(fp, "%s", pnp->id);
		for(j = 0; j < fileCount; ++j){
			fprintf(fp, "\t%.6e", table[i][j]);
		}
		fprintf(fp, "\n");
		pnp = pnp->next;
	}
	fclose(fp);
	return;
}

