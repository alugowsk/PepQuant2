/*
 * protein.h                                                                 
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
 * Header file for protein.c. Contains structs and functions for creating a  
 *     chain of proteins and their respective intensities.                   
 *                                                                           
 * +ASCII art via Text ASCII Art Generator by Patrick Gillespie              
 */

#ifndef PROTEIN_H
#define PROTEIN_H

#include "fasta.h" //FastaPointer
#include "peptide.h" //PeptidePointer

/*
 * proteinNode - A node for a linked list containing the ids and intensities 
 *     of proteins mapped to by identified proteins.
 */
typedef struct proteinNode {
	char *id;
	double *intensities;
	double *spectralCounts;
	double *heavy;
	double *light;
	struct proteinNode *next;
} ProteinNode, *ProteinNodePointer;

/*
 * pep2prot - Convert a table of peptide intensities into a linked list of 
 *     proteins with intensity information. Return the number of proteins in
 *     the linked list.
 */
int pep2prot(ProteinNodePointer *root, FastaPointer fasta,
	PeptidePointer *peptides, double **pepQuant, double **spectralCounts,
	int peptideCount, int fileCount, char ***proteinMap);

/*
 * prot2table - Copy the intensities values contained in a protein list to a
 *     2D array of appropriate size.
 */
void prot2table(ProteinNodePointer pnp, double ***protQuant,
	double ***protHQuant, double ***protLQuant, double ***protSpectralCounts,
	int proteinCount, int fileCount);

/*
 * delProteinList - Free all memory allocated for the protein linked list.
 *     Return pointer to linked list which should be NULL.
 */
ProteinNodePointer delProteinList(ProteinNodePointer pnp);

/*
 * printProtTable - Print a table of protein intensities for every protein, in
 *     the protein list and every file in the filelist.
 */
void printProtTable(ProteinNodePointer pnp, SpectraFileNodePointer filelist,
	int proteinCount, int fileCount, double **table, char *filename);

/*
 * delProteinMap - Free all memory allocated for the proteinMap.
 */
char **delProteinMap(char** proteinMap, int peptideCount);

#endif

