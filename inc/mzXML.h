/*
 * mzXML.h                                                                   
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
 * Header file for mzXML.c. Contains the definition of mzXML related         
 *     structs and declaration of functions for parsing mzXML files.         
 *                                                                           
 * +ASCII art via Text ASCII Art Generator by Patrick Gillespie              
 */

#ifndef MZXML_H
#define MZXML_H

#define MIN_PEAK_COUNT 1
#define EMPTY_PEAK_LIST -5

/*
 * scan - Description of a mass spec scan including decoded peak list.
 */
typedef struct scan {
	int scanNum;
	int msLevel;
	int peaksCount;
	float retentionTime;
	float totalIonCurrent;
	float precIntensity;
	int precCharge;
	float precMz;
	float *mzList;
	float *intList;
} Scan, *ScanPointer, **ScanDB;

/*
 * mzXML - collection of mass spec scans for single LC-MS run.
 */
typedef struct mzxml{
	char *filename;
	int scanCount;
	struct scan **scans;
}MZXML, *MZXMLPointer;

/*
 * delMZXML - Free all memory allocated for the MZXML
 */
MZXMLPointer delMZXML(MZXMLPointer mp);

/*
 * readMZXML - Parse the mzXML file identified by filename and store extracted
 *     data in the passed MZXMLPointer. Return 0 if operations completed
 *     successfully, -1 otherwise. 
 */
int readMZXML(char *filename, MZXMLPointer *mzXML );

#endif

