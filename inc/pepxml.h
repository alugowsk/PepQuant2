/*
 * pepxml.h                                                                  
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
 * Functions for parsing and manipulating pepxml files containing search     
 *     results for LC-MSMS run.                                              
 *                                                                           
 * +ASCII art via Text ASCII Art Generator by Patrick Gillespie              
 */

#ifndef PEPXML_H
#define PEPXML_H

#include "peptide.h" //PeptidePointer, addPeptide

/*
 * parsePepXMLdir -  given the path of a directory of pepXML files parse all
 *     pepXML files and store the peptide identifications and related info
 *     in a Peptide linked list and return a pointer to the root node. 
 */
PeptidePointer parsePepXMLdir(char *dirName);

#endif

