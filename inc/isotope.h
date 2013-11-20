/*
 * isotope.h                                                                 
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
 * Header file for isotope.c. Contains the definition of structs for         
 *     isotopic pattern calculation. The functions declared provide an       
 *     interface for the creation of isotopic patterns for peptide chains    
 *     supporting a limited amount of modification currently in the          
 *     MaxQuant format e.g. "(ac)". Isotopic patterns and collections of     
 *     isotopic patterns must be deleted using appropriate functions by the  
 *     user.                                                                 
 *                                                                           
 * +ASCII art via Text ASCII Art Generator by Patrick Gillespie              
 */

#ifndef ISOTOPE_H
#define ISOTOPE_H

#define PROTON 1.007276466812 //the mass of a proton in amu
#define ISOTOPE_COUNT 32 //the number of isotopes to track
#define COLLECTION_SIZE 128 //number of isotopic patterns in collection
#define AMINO_ACIDS "ARNDCQEGHILKMFPSTUWYV" //legal amino acids

/*
 * isotopicPattern - the isotopic pattern of an atom or molecule. A properly  
 *     initialized isotopicPattern will maintain a list of most intense
 *     isotopes in decreasing order of intensity.
 */
typedef struct isotopicPattern{
	float *intensity;
	float *mass;
}IsotopicPattern, *IsotopicPatternPointer;

/*
 * delIsotopicPattern - Free all memory allocated for the IsotopicPattern.
 */
IsotopicPatternPointer delIsotopicPattern(IsotopicPatternPointer ipp);

/*
 * printIsotopicPattern - Prints IsotopicPattern as two tab seperated columns.                                           *
 */
void printIsotopicPattern(IsotopicPatternPointer ipp);

/*
 * newIPCollection - Make a collection of IsotopicPatterns for amino acids and
 *     modifications. Current implementation generates IsotopicPatterns of
 *     amino acids without water so that they can be further used as building
 *     blocks for peptides without the need to remove atoms.
 */
void newIPCollection(IsotopicPatternPointer **IPCollection);

/*
 * delIPCollection - Free all memory allocated for the collection of 
 *     IsotopicPatterns.
 */
IsotopicPatternPointer* delIPCollection(IsotopicPatternPointer *IPCollection);

/*
 * makePeptide - Create an IsotopicPattern for a peptide by parsing a string
 *     representation of the peptide. Return a pointer to that peptide. Current
 *     implementation supports 21 amino acids with carbidomethyl cysteine, and
 *     (ac), (ph), (ox) MaxQuant mods.
 */
IsotopicPatternPointer makePeptide(
	IsotopicPatternPointer *IPCollection, char *sequence);

#endif

