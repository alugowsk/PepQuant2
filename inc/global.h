/*
 * global.h                                                                   
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
 * Access to global variables and function for parsing command line arguments.                        
 *                                                                           
 * +ASCII art via Text ASCII Art Generator by Patrick Gillespie              
 */

#ifndef GLOBAL_H
#define GLOBAL_H

#include <stdbool.h>

extern int alignWindow;
extern double corrCutOff;
extern int peakWindow;
extern char *fastaName;
extern double intCutOff;
extern int lys;
extern int maxCharge;
extern double ppmCutOff;
extern int arg;
extern int isotopicStates;
extern int threadCount;
extern int quantWindow;
extern char *statQuestdir;
extern int statQuestcutoff;
extern char *pepXMLdir;
extern char *maxQuant;
extern bool ignoreModSite;
extern char **dataList;

extern const char *gitversion;
extern const char *commit;

void parseArgs(int argc, char *argv[]);

#endif

