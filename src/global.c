/*
 * global.c                                                                   
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
 * Function for parsing command line arguments and settigns global variables.
 *                                                                           
 * +ASCII art via Text ASCII Art Generator by Patrick Gillespie              
 */

#include "global.h"

#include <stdlib.h> //atoi, atof
#include <stdio.h> //fprintf
#include <stdbool.h>

/*
 * printUsage - Print simple usage instructions.
 */
void printUsage();

/*
 * printHelp - Print more verbose help, describing options and what they do.
 */
void printHelp();

void parseArgs(int argc, char *argv[]){
	int i = 1;

	while(i < argc){
		switch( (int)argv[i][1]){
			case 'a':
			case 'A':
				alignWindow = atoi(argv[i+1]);
				i+=2;
				break;
			case 'c':
			case 'C':
				corrCutOff = atof(argv[i+1]);
				i+=2;
				break;
			case 'e':
			case 'E':
				peakWindow = atoi(argv[i+1]);
				i+=2;
				break;
			case 'f':
			case 'F':
				fastaName = argv[i+1];
				i+=2;
				break;
			case 'h':
			case 'H':
				printHelp();
				exit(1);
			case 'i':
			case 'I':
				intCutOff = atof(argv[i+1]);
				i+=2;
				break;
			case 'k':
			case 'K':
				lys = atoi(argv[i+1]);
				i+=2;
				break;
			case 'm':
			case 'M':
				maxCharge = atoi(argv[i+1]);
				i+=2;
				break;
			case 'o':
			case 'O':
				ignoreModSite = true;
				++i; // one arg, special case
				break;
			case 'p':
			case 'P':
				ppmCutOff = atof(argv[i+1]);
				i+=2;
				break;
			case 'q':
			case 'Q':
				quantWindow = atoi(argv[i+1]);
				i+=2;
				break;
			case 'r':
			case 'R':
				arg = atoi(argv[i+1]);
				i+=2;
				break;
			case 's':
			case 'S':
				isotopicStates = atoi(argv[i+1]);
				i+=2;
				break;
			case 't':
			case 'T':
				threadCount = atoi(argv[i+1]);
				i+=2;
				break;
			case 'x':
			case 'X':
				statQuestdir = argv[i+1];
				statQuestcutoff = atoi(argv[i+2]);
				i+=3; //two args special case
				break;
			case 'y':
			case 'Y':
				pepXMLdir = argv[i+1];
				i+=2;
				break;
			case 'z':
			case 'Z':
				maxQuant = argv[i+1];
				i+=2;
				break;
			default:
				printUsage();
				exit(1);
		}
	}

	// if fasta and one data source are not valid
	if( !fastaName && 
			!(maxQuant || pepXMLdir || (statQuestdir && statQuestcutoff) ) 
	){
		printUsage();
	}

	return;
}

void printUsage(){
	printf("Usage: pepquant2 -x|y|z data_source -f path_to_fasta [options]\n"
			"\tUse '-h' for help\n");
	return;
}

void printHelp(){
	printf("Usage: pepquant2 -x|y|z data_source -f path_to_fasta [options]\n"
			"Required:\n"
			"\t-f path_to_fasta_file\t\tFasta file used by search engine\n"
			"\n\tand one of:\n\n"
			"\t-x path_to_statQuest_dir int\tRun pepquant2 using statQuest\n"
			"\t\t\t\t\twith passed confidence filter results\n"
			"\t-y path_to_pepXML_dir\t\tRun pepquant2 using pepXML results\n"
			"\t-z path_to_msms.txt\t\tRun pepquant2 using MaxQuant results\n"
			"\n\n"
			"Options:\n"
			"\t-a integer\tThe permissible difference between MS2 and MS1 \n"
			"\t\t\tretention times before defaulting to MS2 over MS1\n"
			"\t\t\tDefault = 90\n"
			"\t-c float\tThe correlation cutoff for found isotopic pattern\n"
			"\t\t\tmatching to theoretical isotopic patter\n"
			"\t\t\tDefault = 0.99\n"
			"\t-e integer\tThe maximum time window within which at least one\n"
			"\t\t\tneighbor peak must be found to declare a MS1 retention\n"
			"\t\t\ttime the apical retention time.\n"
			"\t\t\tDefault = 30\n"
			"\t-i float\tThe minimal intensity of an observed isotopic\n"
			"\t\t\tpattern for a given charge state for it to be considered\n"
			"\t\t\ta valid hit.\n"
			"\t\t\tDefault = 1E6\n"
			"\t-k integer\tThe label status of lysine. A value of 6 assumes\n"
			"\t\t\tlysine is made using C13 and 8 assumes lysine is made\n"
			"\t\t\tusing both C13 and N15.\n"
			"\t\t\tDefault = 0\n"
			"\t-m integer\tThe maximum charge to be looked at by PepQuant2\n"
			"\t\t\twhen looking through MS1 spectra from isotopic patterns.\n"
			"\t\t\tDefault = 4\n"
			"\t-o\t\tWhether to ignore modification localization. If set, \n"
			"\t\t\tall mods will be pushed onto the leftmost unmodified and\n"
			"\t\t\tequivalent AA. S and T are considered to be equivalent.\n"
			"\t-p float\tThe PPM tolerance divided by 1E6. Used for\n"
			"\t\t\tgeneration of theoretical isotopic patterns and matching\n"
			"\t\t\tthese patterns to observed patterns.\n"
			"\t\t\tDefault = 0.000010\n"
			"\t-q integer\tThe time window around the calibrated retention\n"
			"\t\t\ttime used for quantification.\n"
			"\t\t\tDefault = 150\n"
			"\t-r integer\tThe label status of arginine. Value of 6 assumes\n"
			"\t\t\targinine is made using C13 and 10 assumes arginine is\n"
			"\t\t\tmade using both C13 and N15.\n"
			"\t\t\tDefault = 0\n"
			"\t-s integer\tThe number of isotopic states used in the\n"
			"\t\t\ttheoretical isotopic pattern.\n"
			"\t\t\tDefault = 4\n"
			"\t-t integer\tThe number of threads to use during isotopic\n"
			"\t\t\tpattern generation and MS1 spectra interrogation.\n"
			"\t\t\tDefault = 1\n"
			);
	return;
}

