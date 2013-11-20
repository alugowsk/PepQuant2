/*
 * peptide.c                                                                 
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
 * Functions for parsing and manipulating tandem mass spectra search result. 
 *                                                                           
 * +ASCII art via Text ASCII Art Generator by Patrick Gillespie              
 */

#include "peptide.h"
#include "common.h" //MIN_CHARGE
#include "global.h" //maxCharge
#include "mzXML.h" //MZXMLPointer, readMZXML, delMZXML
#include "isotope.h" //AMINO_ACIDS

#include <stdio.h> //fprintf
#include <string.h> //strncpy, strlen, memcpy, strcmp
#include <stdlib.h> //malloc, free
#include <ctype.h> //isupper
#include <math.h> //fabs
#include <pthread.h>
#include <semaphore.h>
#include <time.h>
#include <dirent.h>
#include <stdbool.h>

#define STATQUEST_EXT ".txt"

sem_t exit_sem_t;

/*
 * Structure used to package data need for a thread to search for a pep in an
 *     mzXML.
 */
typedef struct spectraPackage {
	struct mzxml *mzXML;
	struct peptide *pep;
}SpectraPackage, *SpectraPackagePointer;

//Sentinel for red-black tree
Peptide TNILL = {NULL, NULL, NULL, NULL, BLACK, NULL, NULL, NULL};

IsotopicPatternPointer *IPCollection; //IPC collection

/*
 * newScanNode - Allocate memory for a new newScanNode and intialize with   
 *     passed values. Return a pointer to new peak, NULL if error occured.   
 */
ScanNodePointer newScanNode(int scanNum, float *intensity, float *corr);

/*
 * delScanNode - Free all memory allocated for the scanNode. Return pointer to
 *     next scanNode in linked list if exists, otherwise return NULL.
 */
ScanNodePointer delScanNode(ScanNodePointer snp);

/*
 * addScanNode - Add a scan to a linked list of ScanNodes. ScanNodes will be 
 *     added in increasing order.                                            
 */
ScanNodePointer addScanNode(ScanNodePointer root, int scanNum,
	float *intensity, float *corr);

/*
 * delScanNodeList - Free memory allocated for an entire scanNodeList.
 */
ScanNodePointer delScanNodeList(ScanNodePointer snp);

/*
 * newSpectraFileNode - Allocate memory for a new SpectraFileNode and intialize
 *     with passed values. Return a pointer to new SpectraFileNode, NULL if
 *     error occured.
 */
SpectraFileNodePointer newSpectraFileNode(char *rawFile, int scanNum,
	float *intensity, float *corr);

/*
 * addSpectraFileNode - Add a SpectraFileNode to a linked list of
 *     SpectraFileNodes. Create a scanNode for the passed scanNum and add it to
 *     the SpectraFileNode's scan list. Return a pointer to the root of the
 *     list of SpectraFileNodes.
 */
SpectraFileNodePointer addSpectraFileNode(SpectraFileNodePointer root,
	char *rawFile, int scanNum, float *intensity, float *corr);

/*
 * delSpectraFileNode - Free all memory allocated for the SpectraFileNode.  
 *     Return pointer to next SpectraFileNode in linked list if exists,     
 *     otherwise return NULL.                                                
 */
SpectraFileNodePointer delSpectraFileNode(SpectraFileNodePointer sfnp);

/*
 * rbLeftRotate - red black tree left rotate
 */
PeptidePointer rbLeftRotate(PeptidePointer root, PeptidePointer x);

/*
 * rbRightRotate - red black tree right rotate
 */
PeptidePointer rbRightRotate(PeptidePointer root, PeptidePointer y);

/*
 * rbInsertFixup - Ensures red black tree properties are maintained whilst
 *     inserting a node.
 */
PeptidePointer rbInsertFixup(PeptidePointer root, PeptidePointer z);

/*
 * rbDeleteFixup - Ensures red black tree properties are maintained whilst 
 *     removing a node.
 */
PeptidePointer rbDeleteFixup(PeptidePointer root, PeptidePointer x);

/*
 * rbDelete - Remove node z from the red black rooted at root.
 */
PeptidePointer rbDelete(PeptidePointer root, PeptidePointer z);

/*
 * rbFree - Free memory allocated for the PeptidePointer.
 */
PeptidePointer rbFree(PeptidePointer pp);

/*
 * rbTransplant - transplant function for red black trees of PeptidePointers.
 */
PeptidePointer rbTransplant(PeptidePointer root, PeptidePointer x,
	PeptidePointer y);

/*
 * treeMin - Given the root of a peptide tree, return the minimal protein node.
 */
PeptidePointer treeMin(PeptidePointer root);

/*
 * rbPrintInOrder - Given the root of a peptide tree, prints peptides in order.
 *     Used for debugging purposes.
 */
void rbPrintInOrder(PeptidePointer root);

/*
 * newPeptide - Allocate memory for a new peptide and intialize with passed  
 *     values. Return a pointer to new peptide, NULL if error occured.       
 */
PeptidePointer newPeptide(char *rawFile, int scanNum, char *sequence);

/*
 * fillPeptides - Recursive function that populates an array of PeptidePointers
 *     inOrder given a tree of PeptidePointers.F
 */
int fillPeptides(PeptidePointer root, int index, PeptidePointer *peptides);

/*
 * searchSpectra - given a scan and peptide determine whether the peptide's
 *     theoretical isotopic profile is present in the ms1 spectra and if match
 *     meets required conditions then store the hit in appropriate field in 
 *     peptide structure.
 */
void searchSpectra( void *ptr);

/*
 * checkChargeStates - search for mzs of a theoretical isotopic profile within 
 *     the ms1 spectra recording intensities for hits. Return 1 if at least one
 *     charge state showed all required peaks.
 */
int checkChargeStates(ScanPointer sp, float *foundPattern, float *isoMass);

/*
 * binSearch - a traditional binary search except rather than direct equality
 *     we require mzs be within ppmCutOff ppm of each other.
 */
float *binSearch(float mz, float *head, float *tail);

/*
 * checkHit - check if any of the charge states passed the correlation and
 *     intensity cutoffs.
 */
int checkHit(float *intensity, float *foundPattern, float *corr);

/*
 * makePeptideThreadFunc - Helper function for generating an peptide isotopic 
 *     pattern in a threaded environment.
 */
void *makePeptideThreadFunc( void *ptr );

/*
 * sendModsLeft - Rearranges the mods in an amino acid sequence so that mods 
 *     are moved to the leftmost equivalent AA. This allows peptides with the
 *     same mods but different mod localization to be treated the same. 
 */
void sendModsLeft(char *sequence);

///////////////////////////////////////////////////////////////////////////////
//                       END OF FUNCTION DECLARATIONS                        //
///////////////////////////////////////////////////////////////////////////////

ScanNodePointer newScanNode(int scanNum, float *intensity, float *corr){
	ScanNodePointer snp = (ScanNodePointer)malloc(sizeof(ScanNode));
	if(snp == NULL){
		fprintf(stderr, "\nERROR: Out of memory - cannot create scanNode!\n");
	}else{
		snp->scanNum = scanNum;
		snp->rt = 0;
		snp->next = NULL;
		snp->intensity = NULL;
		snp->corr = NULL;
	}
	if(intensity != NULL && corr != NULL){
		int size = maxCharge - MIN_CHARGE + 1;
		snp->intensity = (float*)malloc(size * sizeof(float));
		snp->corr = (float*)malloc(size * sizeof(float));
		if(snp->intensity == NULL || snp->corr == NULL){
			fprintf(stderr,
			 "\nERROR: Out of memory - cannot create scanNode!\n");
		}else{
			memcpy(snp->intensity, intensity, size * sizeof(float));
			memcpy(snp->corr, corr, size * sizeof(float));
		}
	}
	return snp;
}


ScanNodePointer delScanNode(ScanNodePointer snp){
	if(snp == NULL){
		return NULL;
	}
	if(snp->intensity != NULL){
		free(snp->intensity);
	}
	if(snp->corr != NULL){
		free(snp->corr);
	}	
	ScanNodePointer next = snp->next;
	snp->next = NULL;
	free(snp);
	snp = NULL;
	return next;
}


ScanNodePointer addScanNode(ScanNodePointer root, int scanNum,
	float *intensity, float *corr){
	
	ScanNodePointer cur;
	ScanNodePointer prev;
	for(cur = root, prev = NULL;
		cur != NULL && cur->scanNum < scanNum;
		prev = cur, cur = cur->next)
		;
	if(cur != NULL && cur->scanNum == scanNum){
		return root; //scan linked list does not change
	}
	ScanNodePointer new = newScanNode(scanNum, intensity, corr);
	new->next = cur;
	if(!prev){
		return new;
	}else{
		prev->next = new;
		return root;
	}

}


ScanNodePointer delScanNodeList(ScanNodePointer snp){
	if(snp == NULL){
		return NULL;
	}
	while(snp != NULL){
		snp = delScanNode(snp);
	}	
	return snp;
}


SpectraFileNodePointer newSpectraFileNode(char *rawFile, int scanNum,
	float *intensity, float *corr){
	SpectraFileNodePointer sfnp = 
		(SpectraFileNodePointer)malloc(sizeof(SpectraFileNode));
	if(sfnp == NULL){fprintf(stderr,
		"\nERROR: Out of memory - cannot create spectraFileNode!\n");
	}else{
		sfnp->rawFile = (char *)malloc( (strlen(rawFile) + 1)*sizeof(char));
		if(sfnp->rawFile == NULL){
			fprintf(stderr,
				"\nERROR: Out of memory - cannot create spectraFileNode!\n");
		}else{
			strncpy(sfnp->rawFile, rawFile, strlen(rawFile)+1 ); 
			//sfnp->rawFile[strlen(rawFile)] = '\0';
			sfnp->next = NULL;
			sfnp->scans = newScanNode(scanNum, intensity, corr);
		}
	}
	return sfnp;
}


SpectraFileNodePointer addSpectraFileNode(SpectraFileNodePointer root,
	char *rawFile, int scanNum, float *intensity, float *corr){

	SpectraFileNodePointer cur;
	SpectraFileNodePointer prev;
	for(cur = root, prev = NULL;
		cur != NULL && strcmp(cur->rawFile, rawFile) > 0;
		prev = cur, cur = cur->next)
		;
	if(cur != NULL && !strcmp(cur->rawFile, rawFile)){
		cur->scans = addScanNode(cur->scans, scanNum, intensity, corr);
		return root; //sfn linked list does not change
	}
	SpectraFileNodePointer new = newSpectraFileNode(rawFile, scanNum,
			intensity, corr);
	new->next = cur;
	if(!prev){
		return new;
	}else{
		prev->next = new;
		return root;
	}
}


SpectraFileNodePointer delSpectraFileNode(SpectraFileNodePointer sfnp){
	if(sfnp == NULL){
		return NULL;
	}
	SpectraFileNodePointer next =  sfnp->next;
	sfnp->next = NULL;
	if(sfnp->rawFile != NULL){
		free(sfnp->rawFile);
	}
	if(sfnp->scans != NULL){
		sfnp->scans = delScanNodeList(sfnp->scans);
	}
	free(sfnp);
	sfnp = NULL;
	return next;
}


PeptidePointer rbLeftRotate(PeptidePointer root, PeptidePointer x){
	PeptidePointer y = x->right;
	x->right = y->left;
	if(y->left != &TNILL){
		y->left->parent = x;
	}
	y->parent = x->parent;
	if(x->parent == &TNILL){
		root = y;
	}else if(x == x->parent->left){
		x->parent->left = y;
	}else{
		x->parent->right = y;
	}
	y->left = x;
	x->parent = y;
	
	return root;
}


PeptidePointer rbRightRotate(PeptidePointer root, PeptidePointer y){
	PeptidePointer x = y->left;
	y->left = x->right;
	if(x->right != &TNILL){
		x->right->parent = y;
	}
	x->parent = y->parent;
	if(y->parent == &TNILL){
		root = x;
	}else if(y == y->parent->right){
		y->parent->right = x;
	}else{
		y->parent->left = x;
	}
	x->right = y;
	y->parent = x;
	
	return root;
}


PeptidePointer rbInsertFixup(PeptidePointer root, PeptidePointer z){
	while (z->parent->colour == RED){
		if(z->parent == z->parent->parent->left){ //z->p is a left child
			PeptidePointer y = z->parent->parent->right; //sibling of z->p
			if(y->colour == RED){
				z->parent->colour = BLACK;
				y->colour = BLACK;
				z->parent->parent->colour = RED;
				z = z->parent->parent;
			}else{ //y is black
				if(z == z->parent->right){ //z is right
					z = z->parent;
					root = rbLeftRotate(root, z);
				}
				z->parent->colour = BLACK;
				z->parent->parent->colour = RED;
				root = rbRightRotate(root, z->parent->parent);
			}
		}else{ //z-> is a right child
			PeptidePointer y = z->parent->parent->left;
			if(y->colour == RED){
				z->parent->colour = BLACK;
				y->colour = BLACK;
				z->parent->parent->colour = RED;
				z = z->parent->parent;
			}else{ //y is black
				if(z == z->parent->left){
					z = z->parent;
					root = rbRightRotate(root, z);
				}
				z->parent->colour = BLACK;
				z->parent->parent->colour = RED;
				root = rbLeftRotate(root, z->parent->parent);
			}	
		}
	}
	root->colour = BLACK;
	return root;
}


PeptidePointer rbDeleteFixup(PeptidePointer root, PeptidePointer x){
	while(x != root && x->colour == BLACK){
		if(x == x->parent->left){
			PeptidePointer w = x->parent->right;
			if(w->colour == RED){
				w->colour = BLACK;
				x->parent->colour = RED;	
				root = rbLeftRotate(root, x->parent);
				w = x->parent->right;
			}
			if(w->left->colour == BLACK && w->right->colour == BLACK){
				w->colour = RED;
				x = x->parent;
			}else{
				if(w->right->colour == BLACK){
					w->left->colour = BLACK;
					w->colour = RED;
					root = rbRightRotate(root, w);
					w = x->parent->right;
				}
				w->colour = x->parent->colour;
				x->parent->colour = BLACK;
				w->right->colour = BLACK;
				root = rbLeftRotate(root, x->parent);
				x = root;
			}
		}else{
			PeptidePointer w = x->parent->left;
			if(w->colour == RED){
				w->colour = BLACK;
				x->parent->colour = RED;	
				root = rbRightRotate(root, x->parent);
				w = x->parent->left;
			}
			if(w->right->colour == BLACK && w->left->colour == BLACK){
				w->colour = RED;
				x = x->parent;
			}else{
				if(w->left->colour == BLACK){
					w->right->colour = BLACK;
					w->colour = RED;
					root = rbLeftRotate(root, w);
					w = x->parent->left;
				}
				w->colour = x->parent->colour;
				x->parent->colour = BLACK;
				w->left->colour = BLACK;
				root = rbRightRotate(root, x->parent);
				x = root;
			}
		}
	}
	x->colour = BLACK;
	return root;
}


PeptidePointer rbDelete(PeptidePointer root, PeptidePointer z){
	PeptidePointer x;
	PeptidePointer y = z;
	NodeColour origYcolour = y->colour;
	if(z->left == &TNILL){
		x = z->right;
		root = rbTransplant(root, z, z->right);	
	}else if(z->right == &TNILL){
		x = z->left;
		root = rbTransplant(root, z, z->left);	
	}else{
		y = treeMin(z->right);
		origYcolour = y->colour;
		x = y->right;
		if(y->parent == z){
			x->parent = y;
		}else{
			root = rbTransplant(root, y, y->right);	
			y->right = z->right;	
			y->right->parent = y;
		}
		root = rbTransplant(root, z, y);
		y->left = z->left;
		y->left->parent = y;
		y->colour = z->colour;
	}
	if(origYcolour == BLACK){
		root = rbDeleteFixup(root, x);
	}
	z = rbFree(z);
	return root;
}


PeptidePointer rbFree(PeptidePointer pp){
	if(pp == NULL){
		return NULL;
	}
	pp->parent = NULL;
	pp->left = NULL;
	pp->right = NULL;

	if(pp->sequence != NULL){
		free(pp->sequence);
	}
	if(pp->ip != NULL){
		pp->ip = delIsotopicPattern(pp->ip);
	}
	if(pp->spectraFiles != NULL){
		pp->spectraFiles = delSpectraFileList(pp->spectraFiles);
	}
	if(pp->ms1SpectraFiles != NULL){
		pp->ms1SpectraFiles = delSpectraFileList(pp->ms1SpectraFiles);
	}
	free(pp);
	pp = NULL;
	return pp;
}


PeptidePointer rbTransplant(PeptidePointer root, PeptidePointer x,
	PeptidePointer y){
	
	if(x->parent == &TNILL){
		root = y;
	}else if(x == x->parent->left){
		x->parent->left = y;
	}else{
		x->parent->right = y;
	}
	y->parent = x->parent;

	return root;
}


PeptidePointer treeMin(PeptidePointer root){
	while(root->left != &TNILL){
		root = root->left;
	}
	return root;
}


void rbPrintInOrder(PeptidePointer root){
	if(root != &TNILL){
		rbPrintInOrder(root->right);
		printf("[%s]%s\n", root->colour == RED? "RED" : "BLACK",
			root->sequence);
		rbPrintInOrder(root->left);
	}
}


PeptidePointer newPeptide(char *rawFile, int scanNum, char *sequence){
	PeptidePointer pp = (PeptidePointer)malloc(sizeof(Peptide));
	if(!pp){
		fprintf(stderr,
			 "\nERROR: Out of memory - cannot create Peptide!\n");
	}else{
		pp->sequence = malloc(strlen(sequence)+1);
		if(!pp->sequence){
			fprintf(stderr,
			 "\nERROR: Out of memory - cannot create Peptide!\n");
			free(pp);
		}else{
			strncpy(pp->sequence, sequence, strlen(sequence)+1);
			pp->spectraFiles = newSpectraFileNode(rawFile, scanNum, NULL, NULL);
			pp->ms1SpectraFiles = NULL;
			pp->ip = NULL;
			pp->colour = RED;
			pp->parent = &TNILL;
			pp->left = &TNILL;
			pp->right = &TNILL;
		}
	}
	return pp;
}


int fillPeptides(PeptidePointer root, int index, PeptidePointer *peptides){
	if(root != &TNILL){
		index = fillPeptides(root->left, index, peptides);
		peptides[index++] = root;
		index = fillPeptides(root->right, index, peptides);
	}
	return index;
}


void searchSpectra(void *ptr ){
	int i, j, k;

	MZXMLPointer mzXML = ((SpectraPackagePointer)ptr)->mzXML;
	PeptidePointer pp = ((SpectraPackagePointer)ptr)->pep;

	float *isoMass = (float*)malloc(
				isotopicStates*(maxCharge-MIN_CHARGE+1)*2*sizeof(float));
	float *foundPattern = isoMass + isotopicStates*(maxCharge-MIN_CHARGE+1);

	float *corr = (float*)malloc( (maxCharge-MIN_CHARGE+1)*2*sizeof(float));
	float *intensity = corr + maxCharge-MIN_CHARGE+1;

	for(k = 0; k < mzXML->scanCount; ++k){
		if(mzXML->scans[k]->msLevel == 1){
			if(mzXML->scans[k]->peaksCount >= MIN_PEAK_COUNT){
			ScanPointer sp = mzXML->scans[k];	
			for(i = MIN_CHARGE; i <= maxCharge; ++i){
				for(j = 0; j < isotopicStates; ++j){
					isoMass[(maxCharge - i)*isotopicStates+j] = 
						(pp->ip->mass[j]+(PROTON*i))/i;
					foundPattern[(maxCharge - i)*isotopicStates+j] = 0;
				}
			}

			/*search for hits*/
			int valid = checkChargeStates(sp, foundPattern, isoMass);

			/*check hit correlation and validity. record valid hits*/
			if(valid &&
				pearson(pp->ip->intensity, foundPattern, corr) >= corrCutOff &&
				checkHit(intensity, foundPattern, corr) ){

				pp->ms1SpectraFiles = addSpectraFileNode(pp->ms1SpectraFiles,
				mzXML->filename, sp->scanNum, intensity, corr);
			}
		}
		}
	}
	free(isoMass);
	free(corr);	
	sem_post (&exit_sem_t);
}


int checkChargeStates(ScanPointer sp, float *foundPattern, float *isoMass){

	float *head = sp->mzList;
	float *tail = head + sp->peaksCount - 1;
	int i , j;
	for(i = 0; i < isotopicStates*maxCharge; ++i){	
		float *index = binSearch(isoMass[i], head, tail);
		if(index != NULL){
			foundPattern[i] += sp->intList[index-head]; 
			/*check lower adjacent indices*/
			float *tempIndex = index - 1;
			while(tempIndex >= head && fabs( (*tempIndex) -
				isoMass[i])/ isoMass[i] < ppmCutOff){ 

				foundPattern[i] += sp->intList[tempIndex-head]; 
				tempIndex--;
			}
			/*check higher adjacent indices*/
			tempIndex = index + 1;
			while(tempIndex <= tail && fabs( (*tempIndex) -
				isoMass[i])/ isoMass[i] < ppmCutOff){ 

				foundPattern[i] += sp->intList[tempIndex-head]; 
				tempIndex++;
			}
		}
	}
	for(i = MIN_CHARGE-1; i < maxCharge; ++i){
		int peaks = 0;
		for(j = 0; j < isotopicStates; ++j){
			if( foundPattern[isotopicStates * i + j] > 0 ){
				peaks++;
			} 
		}
		if(peaks == isotopicStates){	
			return 1; //all isotopic states for a given charge were detected
		}
	}
	return 0;
}


float *binSearch(float mz, float *head, float *tail){
	while(head < tail && mz > *head && mz < *tail){
		float *mid = head + (tail-head)/2; //find middle
		if(fabs(*mid - mz)/mz < ppmCutOff){
			return mid;
		}else if(mz > *mid){
			head = mid+1;
		}else{
			tail= mid-1;
		}		
	}
	return NULL;
}


int checkHit(float *intensity, float *foundPattern, float *corr){
	int i, j;
	int valid = 0;
	for(i = MIN_CHARGE-1; i < maxCharge; ++i){
		intensity[i] = 0;
		if(corr[i] >= corrCutOff){
			for(j = 0; j < isotopicStates; ++j){
				intensity[i] += foundPattern[i*isotopicStates + j];
			}
			if(intensity[i] < intCutOff){
				intensity[i] = 0;
			}else{	
				valid = 1;
			}
		}else{ //prevent spurious hits from carrying over with valid hits
			corr[i] = 0;
			intensity[i] = 0;
		}
	}	
	return valid;
}


void *makePeptideThreadFunc( void *ptr ){
	PeptidePointer pp = (PeptidePointer)ptr;
	pp->ip = makePeptide(IPCollection, pp->sequence);
	sem_post(&exit_sem_t);
	
	return NULL;
}


void sendModsLeft(char *sequence){

	for(int i = strlen(sequence)-1; i >= 0; --i){
		if(sequence[i] == ')'){ // found a mod
			char moddedAA = sequence[i-4];
			for(int j = 0; j < i-4; ++j){ // for AAs upstream of modded AA
				// if an unmodified version of the same AA exists
				// special case is S and T which can sub each other
				if( (sequence[j] == moddedAA || 
					 ( strchr("ST", sequence[j]) && strchr("ST", moddedAA) )
					) && sequence[j+1] != '('){

					char mod[4];
					strncpy(mod, sequence+i-3, 4); // copy the mod string

					// shift sequence to fill removed mod
					char *p = sequence+i;
					char *s = p-4;
					while(s != sequence + j){
						*p-- = *s--;
					}
					strncpy(sequence+j+1, mod, 4); // copy mod in new loc
					break;
				}	
			}			
		}
	}
	return;
}


SpectraFileNodePointer delSpectraFileList(SpectraFileNodePointer sfnp){
	if(sfnp == NULL){
		return NULL;
	}
	while(sfnp != NULL){
		sfnp = delSpectraFileNode(sfnp);
	}	
	return sfnp;
}


void printPeptides(PeptidePointer *pp, int peptideCount){
	FILE *fp = fopen("peptides.txt", "w");
	if (fp == NULL)	{
   		fprintf(stderr, "Error opening file: %s !\n", "peptides.txt");
	}
	int i, j;
	for(i = 0; i < peptideCount; ++i){
		for(j = 0; j < isotopicStates; ++j){
			fprintf(fp, "[%8.3f, %.4f]",
				pp[i]->ip->mass[j], pp[i]->ip->intensity[j]);
		}
		fprintf(fp, " %s\n", pp[i]->sequence);
	}

	fclose(fp);
}


PeptidePointer delPeptideList(PeptidePointer root){
	while(root != &TNILL){
		root = rbDelete(root, root);
	}
	root = NULL;
	return root;
}


PeptidePointer parseMaxQuant(char *filename, PeptidePointer pp){

	if(filename == NULL){
		fprintf(stderr,
			"\nERROR: Could not parse search results. Filename was NULL.\n");
		exit(1);
	}

	/*try to open parseMaxQuant results file*/
	FILE *fp = fopen(filename, "r");
	if(fp == NULL){
		printf("\nERROR: opening %s\n", filename);
		exit(1);
	}else{
		char line[8096];
		while ( fgets (line ,8095, fp) != NULL ){
			int i;
			char *tokens = strtok(line, "\t");
			/*first line contains headers and must be skipped*/
			if(strcmp("Raw file", tokens) ){
				char *rawFile = NULL;
				int scanNum = -1;
				char *sequence = NULL;
				for(i = 0; i < 10; ++i){
					if(i == 0){
						rawFile = tokens;
					}else if(i == 2){
						scanNum = atoi(tokens);
					}else if(i == 9){
						sequence = tokens;
					}
					tokens = strtok(NULL, "\t");
				}				
				pp = addPeptide(pp, strcat(rawFile, ".mzXML"), scanNum,
					sequence);
				if( (lys && strchr(sequence, 'K')) || //if heavy K and seq has K
					(arg && strchr(sequence, 'R')) ){ //if heavy R and seq has R
					sequence[0] = '*'; //mark seq as heavy
					pp = addPeptide(pp, rawFile, scanNum, sequence);
				}
			}
		}
		fclose (fp);     
	}
	return pp;
}


PeptidePointer addPeptide(PeptidePointer root, char *rawFile, int scanNum,
	char *sequence){

	if(ignoreModSite){
		sendModsLeft(sequence);
	}

	PeptidePointer y = &TNILL;
	PeptidePointer x = root;	
	while(x != &TNILL){
		y = x;
		if(!strcmp(sequence, x->sequence)){
			x->spectraFiles = addSpectraFileNode(x->spectraFiles, rawFile,
			scanNum, NULL, NULL);
			return root;
		}else if(strcmp(sequence, x->sequence) < 0){
			x = x->left;
		}else{
			x = x->right;
		}
	}
	PeptidePointer z = newPeptide(rawFile, scanNum, sequence);
	z->parent = y;
	if(y == &TNILL){
		root = z;
	}else if(strcmp(z->sequence, y->sequence) < 0){
		y->left = z;
	}else{
		y->right = z;
	}
	root = rbInsertFixup(root, z);
	return root;
}


char *stripMods(char *moddedPeptide){
	char peptideBuffer[ strlen(moddedPeptide) ];
	int i = 0;
	int j = 0;

	while(j < strlen(moddedPeptide) ){
		if(strchr(AMINO_ACIDS, moddedPeptide[i])){
			peptideBuffer[j] = moddedPeptide[i];
			j++;
			i++;
		}else if(moddedPeptide[i] == '\0'){
			/*Pad end of the buffer with NULL to be safe*/
			peptideBuffer[j] = '\0';
			j++;
		}else{
			/*Assume anything not NULL or AA is a mod*/
			i++;
		}
	}
	
	char *peptide = (char *)malloc( (strlen(peptideBuffer)+1)*sizeof(char) );
	strncpy(peptide, peptideBuffer, strlen(peptideBuffer));
	peptide[strlen(peptideBuffer)] = '\0';

	return peptide;
}


PeptidePointer initIsotopicPatterns(PeptidePointer root,
	PeptidePointer *peptides, int peptideCount){

	int i;
	sem_unlink("exit_sem_t");
	sem_init (&exit_sem_t, 0, threadCount);
  	pthread_attr_t attr;
 	pthread_attr_init (&attr);
 	pthread_attr_setdetachstate (&attr, PTHREAD_CREATE_DETACHED);

	/*init isotopic pattern collection*/
	IPCollection = (IsotopicPatternPointer*)
		malloc(COLLECTION_SIZE * sizeof(IsotopicPatternPointer));
	newIPCollection(&IPCollection);
	/*calculate isotopic pattern for each peptide*/

	for(i = 0; i < peptideCount; i++){
		sem_wait (&exit_sem_t);
		pthread_t thread;
		pthread_create( &thread, &attr, makePeptideThreadFunc,
			(void *)peptides[i]);
	}
	for(i = 0; i < threadCount; i++){
		sem_wait (&exit_sem_t);
	}

	sem_destroy(&exit_sem_t);
	/*delete isotopic pattern collection*/
	
	for(i = 0; i < peptideCount; i++){
		if(peptides[i]->ip == NULL){/*happens if ip generation failed*/
			root = rbDelete(root, peptides[i]);
		}
	}
	IPCollection = delIPCollection(IPCollection);
		
	return root;
}


void printSearchResults(PeptidePointer *peptides, int peptideCount){

	FILE *fp = fopen("searchResults.txt", "w");
	if (fp == NULL)
		{
   		 printf("Error opening file!\n");
   		 exit(1);
	}

	FILE *fp2 = fopen("searchResultsMS2.txt", "w");
	if (fp == NULL)
		{
   		 printf("Error opening file!\n");
   		 exit(1);
	}
	int i;
	for(i = 0; i < peptideCount; ++i){
		fprintf(fp, "sequence: %s\n", peptides[i]->sequence);
		fprintf(fp2, "sequence: %s\n", peptides[i]->sequence);
		SpectraFileNodePointer sfnp = peptides[i]->ms1SpectraFiles;
		/*print ms1 info*/
		while(sfnp != NULL){
			fprintf(fp, "\t%s\n", sfnp->rawFile);
			ScanNodePointer snp = sfnp->scans;
			while(snp != NULL){
				fprintf(fp, "\t\t%d\t%.4e\t", snp->scanNum, snp->rt);
				int j;
				for(j = 0; j < maxCharge - MIN_CHARGE + 1; j++){
					fprintf(fp, "|%.2f - %.2e|\t", snp->corr[j],
						snp->intensity[j]);
				}
				fprintf(fp, "\n");
				snp = snp->next;
			}
			sfnp = sfnp->next;
		}
		/*print ms2 info*/
		sfnp = peptides[i]->spectraFiles;
		while(sfnp != NULL){
			fprintf(fp2, "\t%s\n", sfnp->rawFile);
			ScanNodePointer snp = sfnp->scans;
			while(snp != NULL){
				fprintf(fp2, "\t\t%d\t%.4e\n", snp->scanNum, snp->rt);
				snp = snp->next;
			}
			sfnp = sfnp->next;
		}
	}
	fclose(fp);
	fclose(fp2);
	return;
}


void printTable(PeptidePointer *peptides, SpectraFileNodePointer filelist,
	int peptideCount, int fileCount, double **table, char *filename){

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
	for(i = 0; i < peptideCount; ++i){
		fprintf(fp, "%s", peptides[i]->sequence);
		for(j = 0; j < fileCount; ++j){
			fprintf(fp, "\t%.6e", table[i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	return;
}


void searchMzXMLs(PeptidePointer *peptides, int peptideCount,
	SpectraFileNodePointer filelist){

	int i, j;
	/* prepare packages for threaded spectra search */
	SpectraPackagePointer *packages = malloc(peptideCount * 
		sizeof(SpectraPackagePointer));
	
	if(packages == NULL){
		fprintf(stderr, "ERROR\n");
		exit(1);
	}else{
		for(j = 0; j < peptideCount; ++j){	
			packages[j] = malloc(sizeof(SpectraPackage));
			if(packages[j] == NULL){
				fprintf(stderr,
					"ERROR: could not allocated memory for spectra search\n");
				exit(1);
			}
		}
		while(filelist != NULL){

			sem_unlink("exit_sem_t");
			sem_init (&exit_sem_t, 0, threadCount);
		  	pthread_attr_t attr;
		 	pthread_attr_init (&attr);
		 	pthread_attr_setdetachstate (&attr, PTHREAD_CREATE_DETACHED);

			time_t start, end;
			time(&start);
	
			/*read mzXML*/
			MZXMLPointer mzXML;
			
			readMZXML(filelist->rawFile, &mzXML);
			printf("\tFile: %s read %d spectra\n",
			filelist->rawFile, mzXML->scanCount);

			/*search mzXML for isotopic patterns*/
			for(j = 0; j < peptideCount; ++j){
				packages[j]->mzXML = mzXML;
				packages[j]->pep = peptides[j];	
				sem_wait (&exit_sem_t);
				pthread_t thread;
				pthread_create( &thread, &attr, (void*)&searchSpectra,
					(void*) packages[j]);
			}
	
			for(i=0; i<threadCount; i++){	
				sem_wait (&exit_sem_t);
			}

			/*get ms1 and ms2 retention time info*/
			for(j = 0; j < peptideCount; ++j){
				SpectraFileNodePointer sfnp = peptides[j]->ms1SpectraFiles;
				while(sfnp != NULL){
					if(!strcmp(sfnp->rawFile, mzXML->filename) ){
						ScanNodePointer snp = sfnp->scans;
						while(snp != NULL){
							snp->rt = 
								mzXML->scans[snp->scanNum-1]->retentionTime;
							snp = snp->next;
						}
						break;
					}
					sfnp = sfnp->next;
				}
				sfnp = peptides[j]->spectraFiles;
				while(sfnp != NULL){
					if(!strcmp(sfnp->rawFile, mzXML->filename) ){
						ScanNodePointer snp = sfnp->scans;
						while(snp != NULL){
							snp->rt = 
								mzXML->scans[snp->scanNum-1]->retentionTime;
							/* since we are not using snp->intensity for
							   ms2 scans lets use it to store TIC*/
							snp->intensity = (float*)malloc(sizeof(float));
							if(snp->intensity == NULL){
								fprintf(stderr,
								"Error allocating memory for ms2 tic info\n");
							}else{
								snp->intensity[0] =
									mzXML->scans[snp->scanNum-1]
										->totalIonCurrent;
								snp = snp->next;
							}
						}
						break;
					}
					sfnp = sfnp->next;
				}
			}	
			delMZXML(mzXML);
			time(&end);
			printf("took %f seconds\n", difftime(end, start));
			filelist = filelist->next;
			sem_destroy(&exit_sem_t);
		}
	}	
	/* recover memory allocated for packages */
	for(j = 0; j < peptideCount; ++j){	
			free(packages[j]);
			packages[j] = NULL;
	}
	free(packages);
}


int initFilelist(SpectraFileNodePointer *filelist, PeptidePointer *peptides,
	int peptideCount){

	SpectraFileNodePointer root = NULL;
	int i;
	for(i = 0; i < peptideCount; ++i){
		SpectraFileNodePointer sfnp = peptides[i]->spectraFiles;
		while(sfnp != NULL){
			root = addSpectraFileNode(root,
				sfnp->rawFile, 0, NULL, NULL);
			sfnp = sfnp->next;
		}
	}
	(*filelist) = root;
	int fileCount = 0;
	while(root != NULL){
		printf("\t%s\n", root->rawFile);
		fileCount++;
		root=root->next;
	}
	return fileCount;
}


PeptidePointer parseStatQuest(char *filename, PeptidePointer pp){

	if(filename == NULL){
		fprintf(stderr,
			"\nERROR: Could not parse search results. Filename was NULL.\n");
		exit(1);
	}

	/*try to open StatQuest file*/
	FILE *fp = fopen(filename, "r");
	if(fp == NULL){
		printf("Error opening %s\n", filename);
		exit(1);
	}else{
		char line[4096];
		while ( fgets (line ,4095, fp) != NULL ){
				int i;
				char *spectra = NULL;
				char *sequence = NULL;
				char *tokens = strtok(line, " ");
				for(i = 0; i < 13; ++i){
					if(i == 1){
						spectra = tokens;
					}else if(i == 12){
						sequence = tokens;
						if(sequence[0] == '+'){
							tokens = strtok(NULL, " ");
							sequence = tokens;
						}
						char *firstPeriod = strchr(sequence, '.');
						char *secondPeriod = strrchr(sequence, '.');
						if(firstPeriod != secondPeriod){
							*(secondPeriod+1) = '_';
						}
						*(firstPeriod-1) = '_';
						*(strchr(sequence, '\n')-1) = '\0';
					}
					tokens = strtok(NULL, " ");
				}				
				tokens = strtok(spectra, ".");
				char*rawFile;
				int scanNum = -1;
				for(i = 0; i < 2; ++i){
					if(i == 0){
						rawFile = tokens;
					}else if(i == 1){
						scanNum = atoi(tokens);
					}
					tokens = strtok(NULL, " ");
				}		

				pp = addPeptide(pp, strcat(rawFile, ".mzXML"), scanNum,
					sequence);
		}
		fclose (fp);     
	}
	return pp;
}


int getCount(PeptidePointer root){
	if(root == &TNILL){
		return 0;
	}else{
		return(getCount(root->left) + 1 + getCount(root->right));
	}	
}


PeptidePointer *inOrder(PeptidePointer root, int count){
	PeptidePointer *peptides = (PeptidePointer *)
		malloc(count * sizeof(PeptidePointer));
	fillPeptides(root, 0, peptides);
	return peptides;
}


PeptidePointer parseStatQuestdir(char *dirName, int cutoff){

	if(cutoff <= 0 || cutoff >= 100){
		fprintf(stderr, 
		"\nERROR: statQuest cutoff should be in the range [0-99]\n");
		exit(1);
	}
	
	DIR *dp;
	dp = opendir(dirName);
	PeptidePointer pp = &TNILL;	
	
	if(dp != NULL){
		struct dirent *ep = NULL;
		char suffix[strlen(STATQUEST_EXT)+2];
		sprintf(suffix, "%02d%s", cutoff, STATQUEST_EXT);
		while ( (ep = readdir (dp)) ){
			if(strstr(ep->d_name, suffix)){
				puts(ep->d_name);
				char fullPath[strlen(dirName) + strlen(ep->d_name) + 2];
				strcpy(fullPath, dirName);
				strcat(fullPath, "/");
				strcat(fullPath, ep->d_name);
				pp = parseStatQuest(fullPath, pp);
			}
		}
		(void) closedir (dp);
	}else{
		fprintf(stderr, "\nERROR: Could not open directory: %s!\n", dirName);
		exit(1);
	}
	
	return pp;
}

