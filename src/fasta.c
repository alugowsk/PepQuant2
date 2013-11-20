/*
 * fasta.c                                                                   
 * =======                                                                   
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
 * Functions for parsing and manipulating fasta files containing peptide     
 *     identifiers and sequence. Current version of file treats sequence     
 *     identifier and an description as a single id.                         
 *                                                                           
 * +ASCII art via Text ASCII Art Generator by Patrick Gillespie              
 */

#include "fasta.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * newEntry - Allocate memory for a new protein and intialize with passed id 
 *     and set the sequence to be null character. Return a pointer to new
 *     entry, NULL if error occured. 
 */
EntryPointer newEntry(char *id);

/*
 * delEntry - Free all memory allocated for the entry.
 */
void delEntry(EntryPointer ep);

/*
 * initCoverage - Initialize the coverage tracking array to zero.
 */
void initCoverage(EntryPointer ep);

/*
 * newFasta - Allocate memory for a new fasta database and intialize to default
 *     size. Return a pointer to new database, NULL if error occured. 
 */
FastaPointer newFasta();

/*
 * addEntry - Insert the entry into the fasta database expanding the database
 *     if required.                                          *
 */
void addEntry(FastaPointer fp, EntryPointer ep);

/*
 * trimFastaDB - Reallocate fasta database to memory chunk of equal size.
 *     Allows for more efficient memory usage after fasta database 
 *     initialization.
 */
void trimFastaDB(FastaPointer fp);

/*
 * appendSequence - Expand the sequence of the protein entry with the passed 
 *     string. 
 */
void appendSequence(EntryPointer ep, char *addSequence);

///////////////////////////////////////////////////////////////////////////////
//                       END OF FUNCTION DECLARATIONS                        //
///////////////////////////////////////////////////////////////////////////////

EntryPointer newEntry(char *id){
	EntryPointer ep = (EntryPointer)malloc(sizeof(Entry));
	if(ep == NULL){
		fprintf(stderr,
			"\nERROR: Out of memory - cannot create database entry.\n");
	}else{
		ep->id = (char*)malloc(sizeof(char)*(strlen(id)+1));
		ep->sequence = (char *)malloc(sizeof(char));
		ep->coverage = NULL;
		if(ep->sequence == NULL || ep->id == NULL){
		fprintf(stderr,
			"\nERROR: Out of memory - cannot create database entry.\n");
		}else{
			strncpy(ep->id, id, (strlen(id)+1)*sizeof(char));
			ep->sequence[0] = '\0';	
		}
	}
	return ep;
}


void delEntry(EntryPointer ep){
	if(ep == NULL){
		return;
	}
	if(ep->sequence != NULL){
		free(ep->sequence);
	}
	if(ep->id != NULL){
		free(ep->id);
	}
	if(ep->coverage != NULL){
		free(ep->coverage);
	}
	free(ep);
	ep = NULL;
}


void addEntry(FastaPointer fp, EntryPointer ep){
	fp->db[fp->index++] = ep;
	
	/*expand database buffer if needed*/
	if(fp->index == fp->size){
		fp->size*=2;
		Database temp = 
			(Database)realloc(fp->db, sizeof(EntryPointer)*fp->size);	
		if(temp == NULL){
			fprintf(stderr, 
				"\nERROR: Out of memory - cannot expand fasta database.\n");
		}else{
			fp->db = temp;
		}
	}
}


FastaPointer newFasta(){
	FastaPointer fp = (FastaPointer)malloc(sizeof(Fasta));
	if(fp == NULL){
		fprintf(stderr, 
			"\nERROR: Out of memory - cannot create fasta database.\n");
	}else{
		fp->size = INIT_DB_SIZE;
		fp->index = 0;
		fp->db = (Database)malloc(fp->size * sizeof(EntryPointer)); 
		if(fp->db == NULL){
			fprintf(stderr, 
				"\nERROR: Out of memory - cannot create fasta database.\n");
			free(fp);
			fp = NULL;
		}
	}
	return fp;
}


void trimFastaDB(FastaPointer fp){
	Database temp =
		(Database)realloc(fp->db, fp->index * sizeof(EntryPointer));
	if(temp == NULL){
		fprintf(stderr, 
			"\nERROR: Out of memory - cannot reallocate fasta database.\n");
	}else{
		fp->db = temp;
		fp->size = fp->index;
	}
}


void appendSequence(EntryPointer ep, char *addSequence){
	int existLength = strlen(ep->sequence);
	int addLength = strlen(addSequence);

	/*reallocate existing sequence to make room for sequence to be appended
	and then concatenate.*/
	char *temp = (char *)realloc(ep->sequence, sizeof(char) *
		 (existLength + addLength + 1) );
	if(temp == NULL){
		fprintf(stderr, 
			"\nERROR: Out of memory - cannot extend protein sequence.\n");
	}else{
		strcat(temp, addSequence);
		ep->sequence = temp;
	}
	
	return;
}


void initCoverage(EntryPointer ep){
	if(ep == NULL || ep->sequence == NULL){
		fprintf(stderr, 
			"\nERROR: Passes fasta entry was not properly initialized.\n");
		return;
	}
	ep->coverage = (int*)malloc(strlen(ep->sequence)*sizeof(int));
	if(ep->coverage != NULL){
		memset(ep->coverage, 0, strlen(ep->sequence)*sizeof(int));
	}else{
		fprintf(stderr, 
			"\nERROR: Out of memory - cannot initialize fasta coverage tracker.\n");
	}
	return;
}


int readFASTA(char *filename, FastaPointer *fasta ){
	/*try to open FASTA file*/
	FILE *fp = fopen(filename, "r");
	if(fp == NULL){
		fprintf(stderr, "\nERROR: opening %s\n", filename);
		exit(1);
	}

	(*fasta) = newFasta();

	/*read from fasta and create entry for every protein*/
	char line[MAX_LINE];
	EntryPointer ep = NULL;
	while(fgets(line,sizeof(line),fp) != NULL){

		/*remove newline*/
		line[strcspn ( line, "\n" )] = '\0';
		line[strcspn ( line, "\r" )] = '\0';

		if(line[0] == '>'){
			ep = newEntry(line);
			if(ep != NULL){
				addEntry((*fasta), ep);
			}
		}else{
			appendSequence(ep, line);
		}
	}

	trimFastaDB((*fasta));

	return 0;
}


void delFasta(FastaPointer fp){
	if(fp == NULL){
		return;
	}
	if(fp->db != NULL){
		int i;
		for(i = 0; i < fp->size; ++i){
			delEntry(fp->db[i]);
			fp->db[i] = NULL;
		}
		free(fp->db);
		fp->db = NULL;
	}
	free(fp);
	fp = NULL;
	return;
}


char *findProtein(FastaPointer fasta, char* sequence){
	int i;
	char * prots = (char*)malloc(sizeof(char) );
	prots[0] = '\0';
	for(i = 0; i < fasta->size; ++i){
		if(strstr(fasta->db[i]->sequence, sequence) != NULL){
			int combinedLength = strlen(prots) + strlen(fasta->db[i]->id) + 2;
			prots = (char *)realloc(prots, combinedLength*sizeof(char) );
			strcat(prots, fasta->db[i]->id);
			strcat(prots, ";"); //seperate proteins in list
			/*record the index of the matched protein will be used later to
			track coverage*/
			fasta->index = i; 
		}
	}
	if(strlen(prots) == 0){
		free(prots);
		return NULL;
	}else{
		return prots;	
	}
}


void trackCoverage(FastaPointer fasta, char *peptide){
	EntryPointer ep = fasta->db[fasta->index];

	if(ep->coverage == NULL){
		initCoverage(ep);
	}
	int offset = 0;
	offset = strstr(ep->sequence, peptide) - ep->sequence;
	int i;
	for(i = 0; i < strlen(peptide); ++i){
		ep->coverage[offset+i]++;
	}

	return;
}


void printCoverage(FastaPointer fasta){

	FILE *fp = fopen("coverage.txt", "w");
	if (fp == NULL){
		fprintf(stderr, "Error opening file: coverage.txt!\n");
	}
	
	int i;
	for(i = 0; i < fasta->size; ++i){
		if(fasta->db[i]->coverage != NULL){
			int j;
			int covered = 0;
			int length = strlen(fasta->db[i]->sequence);
			for(j = 0; j < length; ++j){
				if(fasta->db[i]->coverage[j] > 0){
					covered++;
				}
			}
			fprintf(fp, "%s (%d of %d)\n", fasta->db[i]->id, covered, length);
			
			for(j = 0; j < length; ++j){
				fprintf(fp, "%c\t", fasta->db[i]->sequence[j]);
			}
			fprintf(fp, "\n");
			for(j = 0; j < length; ++j){
				fprintf(fp, "%d\t", fasta->db[i]->coverage[j]);
			}
			fprintf(fp, "\n");
		}
	}	

	fclose(fp);
	return;
}

