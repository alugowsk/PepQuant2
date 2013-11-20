/*
 * peptide.h                                                                 
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
 * Header file for peptide.c. Contains function declartions for parsing      
 *     search engine results and structs for storing them. A set of discovered
 *     peptides is maintained as a red-black tree. 
 *                                                                           
 * +ASCII art via Text ASCII Art Generator by Patrick Gillespie              
 */

#ifndef PEPTIDE_H
#define PEPTIDE_H

/*
 * nodeColour - possible colours for nodes of a red-black tree
 */
typedef enum nodeColour {
	RED,
	BLACK
}NodeColour;

/*
 * scanNode - A node for a linked list containing the scan numbers where a   
 *     peptide was detected for a given spectra file.                 
 */
typedef struct scanNode {
	int scanNum;
	double rt;
	float *corr; //stores correlation for different charge states for ms1
	/*stores intensity for different charge states for ms1,
	total ion current for ms2**/
	float *intensity; 
	struct scanNode *next;
} ScanNode, *ScanNodePointer;

/*
 * spectraFileNode - A node for a linked list containing the spectra files  
 *     where a given peptide was detected in at least one scan.
 */
typedef struct spectraFileNode {
	char *rawFile;
	struct scanNode *scans;
	struct spectraFileNode *next;
} SpectraFileNode, *SpectraFileNodePointer;

/*
 * peptide - A node for a red-black tree containing the peptides detected in a  
 *     tandem mass spectra search.
 */
typedef struct peptide {
	char* sequence;

	struct isotopicPattern *ip;
	struct spectraFileNode *spectraFiles;
	struct spectraFileNode *ms1SpectraFiles;

	NodeColour colour;
	struct peptide *parent;
	struct peptide *left;
	struct peptide *right;
} Peptide, *PeptidePointer;

extern Peptide TNILL; //Sentinel for red-black tree

/*
 * delSpectraFileList - Free all memory allocated for the linked list of 
 *     SpectraFileNodes by repeatedly calling delSpectraFileNode. Function 
 *     should return NULL.
 */
SpectraFileNodePointer delSpectraFileList(SpectraFileNodePointer sfnp);

/*
 * printPeptides - Print the theoretical isotopic patterns for all peptides in
 *     the passed array of PeptidePointers.
 */
void printPeptides(PeptidePointer *pp, int peptideCount);

/*
 * delPeptideList - Free all memory allocated for the peptide tree.
 *     Return pointer to tree which should be NULL.
 */
PeptidePointer delPeptideList(PeptidePointer pp);

/*
 * parseMaxQuant - Given a MaxQuant msms.txt file populate a peptide tree or 
 *     return NULL if something goes wrong.                                  
 */
PeptidePointer parseMaxQuant(char *filename, PeptidePointer pp);

/*
 * addPeptide - Add a new peptide to the red-bleack tree of peptides. Peptides
 *     will be added in increasing order (according to strcmp). The current 
 *     implementation includes modifications when sorting.
 */
PeptidePointer addPeptide(PeptidePointer root, char *rawFile, int scanNum,
	char *sequence);

/*
 * stripMods - Returns a string representation of the passed peptide with all
 *     modifications removed.
 */
char *stripMods(char *moddedPeptide);

/*
 * Generate isotopic patterns for all peptides in the array of peptidePointers.
 *     The function requires access to the orginal tree as well in case a node
 *     needs to be deleted in the case of a failure to generate an isotopic 
 *     pattern. This function will spawn multiple threads.
 */
PeptidePointer initIsotopicPatterns(PeptidePointer root, 
	PeptidePointer *peptides, int peptideCount);

/*
 * printSearchResults - Print peptide sequences and all files and scans in 
 *     which the peptides where identified. For ms1 identifications print a
 *     table of including charge and correlation information.
 */
void printSearchResults(PeptidePointer *peptides, int peptideCount);

/*
 * printTable - Print a table of peptide intensities for every peptide, in the
 *     peptide list and every file in the filelist.
 */
void printTable(PeptidePointer *peptides, SpectraFileNodePointer filelist,
	int peptideCount, int fileCount, double **table, char *filename);

/*
 * searchMzXMLs - For every mzXML file in the mzXML filelist search ms1 spectra
 *     for isotopic patterns of each peptide. Store search results within the 
 *     the corresponding peptide node.
 */
void searchMzXMLs(PeptidePointer *peptides, int peptideCount,
	SpectraFileNodePointer filelist);

/*
 * initFilelist - Using the array of peptides create a set of files from which
 *     the peptide list was derived and return the number of files in that set.
 */
int initFilelist(SpectraFileNodePointer *filelist, PeptidePointer *peptides,
	int peptideCount);

/*
 * getCount - return the number of peptides in the peptide tree
 */
int getCount(PeptidePointer root);

/*
 * inOrder - allocate space for an array of peptidePointers and populate the
 *     array in order using the Peptide tree.
 */
PeptidePointer *inOrder(PeptidePointer root, int count);

/*
 * parseStatQuestdir - Given a directory path containing StatQuest results parse
 *     those results. Only supports SEQUEST searches with no mods. The cutoff
 *      represents the confidence filter used and should be part of every 
 *      statQuest fileaname.
 */
PeptidePointer parseStatQuestdir(char *dirName, int cutOff);

#endif

