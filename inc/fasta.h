/*
 * fasta.h                                                                   
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
 * Header file for fasta.c. Contains the definition of fasta related         
 *     structs and declaration of functions for parsing fasta files.         
 *                                                                           
 * +ASCII art via Text ASCII Art Generator by Patrick Gillespie              
 */

#ifndef FASTA_H
#define FASTA_H

#define MAX_LINE 1024
#define INIT_DB_SIZE 256

/*
 * entry - Description of a single protein. Usually in the format of one line
 *     featuring an identifier/description followed by sequence.
 */
typedef struct entry{
	char* id; 
	char* sequence;
	int* coverage;
}Entry, *EntryPointer, **Database;

/*
 * fasta - A simple container for all entries within a fasta file.
 */
typedef struct fasta{
	int size;
	int index; //can be used to maintain position within entry database
	struct entry **db;
}Fasta, *FastaPointer;

/*
 * readFASTA - Parse file with passed name and store set pointer to initialized
 *     struct.
 */
int readFASTA(char *filename, FastaPointer *fasta );

/*
 * delFasta - Free all memory allocated for the fasta database.
 */
void delFasta(FastaPointer fp);

/*
 * findProtein - Search protein sequences in FASTA database for passed
 *     sequence. Return a string of colon seperated protein ids containing the
 *     passed sequence, NULL otherwise.
 */
char *findProtein(FastaPointer fasta, char* sequence);

/*
 * trackCoverage - Given a peptide that was found to be unique to a protein in
 *     the FASTA file, record the protein sequence covered by the peptide. This
 *     only works for unique peptides and assumes the fasta index member points
 *     to the protein in question.  
 */
void trackCoverage(FastaPointer fasta, char *peptide);

/*
 * printCoverage - Print a table of proteins, their sequences, and the number
 *     of times each amino acid in the sequence was covered by a unique
 *     peptide.
 */
void printCoverage(FastaPointer fasta);

#endif

