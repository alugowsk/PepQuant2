/*
 * pepxml.c                                                                  
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

#include "pepxml.h"
#include "xml.h" //openXML, downTo, getAttribute 

#include <stdlib.h> //malloc
#include <string.h> //
#include <ctype.h> //tolower
#include <math.h> //atof, fabs
#include <stdio.h> //fprintf, printf

#include <dirent.h>

#define MSMS_RUN_SUMMARY "msms_run_summary"
#define BASE_NAME "base_name"
#define SEARCH_SUMMARY "search_summary"
#define AMINOACID_MOD "aminoacid_modification"
#define TERMINAL_MOD "terminal_modification"
#define MASS "mass"
#define DESCRIPTION "description"
#define VARIABLE "variable"
#define YES "Y"
#define SPECTRUM_QUERY "spectrum_query"
#define START_SCAN "start_scan"
#define SEARCH_RESULT "search_result"
#define SEARCH_HIT "search_hit"
#define PEPTIDE "peptide"
#define HIT_RANK "hit_rank"
#define MOD_INFO "modification_info"
#define MOD_NTERM "mod_nterm_mass"
#define MOD_CTERM "mod_cterm_mass"
#define MOD_AMINOACID "mod_aminoacid_mass"
#define POSITION "position"
#define MZXML_EXT ".mzXML"
#define PEPXML_EXT ".pepXML"
#define SEQ_MXZML_DIR ".mzXML_dta"

/*
 * mod - an amino acid modification identified either by id or mass
 */
typedef struct mod{
	char *id;
	double mass;
	struct mod *next;
}Mod, *ModPointer;

/*
 * newMod - given a modification name and mass create a new mod and return a 
 *     pointer to the given mod.
 */
ModPointer newMod(char *name, double mass);

/*
 * addMod -  add a mod to an existing linked list of mods 
 */
ModPointer addMod(ModPointer root, char *name, double mass);

/*
 * delMod - Free memory allocated for a mod. Return a pointer to the next mod
 *     in the linked list if it exists, NULL otherwise.
 */
ModPointer delMod(ModPointer mp);

/*
 * delMods - free memory allocated for the linked list of mods. 
 */
ModPointer delMods(ModPointer mp);

/*
 * getMods - Parse the searchSummary node of a pepXML file and compile a list
 *     of modifications that can be expected in the peptide spectrum matches
 *     in the msmsRunSummary section of the pepXML 
 */
ModPointer getMods(xmlNodePtr searchSummary);

/*
 * convertMod - given the mass of an amino acid modifications convert the mod
 *     to a string id using the previously compile mod linked list.
 */
char *convertMod(ModPointer mp, double mass);

/*
 * modSequence - given a hit node, convert the assigned peptide sequence and 
 *     modification information into a single character sequence in the form
 *     typical of MaxQuant.
 */
char *modSequence(xmlNodePtr hit, xmlChar *seq, ModPointer mp);

/*
 * getSearchResults - search the msmsRunSummary portion of a pepXML for hits
 *     converting pepXML format hits to MaxQuant like modded sequences. Store
 *     hits in the passed peptide linked list returning the new root of the
 *     linked list.
 */
PeptidePointer getSearchResults(PeptidePointer pp, xmlNodePtr msmsRunSummary,
	ModPointer mp, char *mzXMLname);

/*
 * Given the name of a pepXML, open and parse the pepXML storing the peptide 
 *     information in the passed peptide linked list. Return a pointer to the 
 *     new root of the peptide linked list. 
 */
PeptidePointer readPepXML(char *filename, PeptidePointer pp);

///////////////////////////////////////////////////////////////////////////////
//                       END OF FUNCTION DECLARATIONS                        //
///////////////////////////////////////////////////////////////////////////////


ModPointer newMod(char *name, double mass){
	ModPointer mp = (ModPointer)malloc( sizeof(Mod) );
	if(mp != NULL){
		mp->mass = mass;
		mp->next = NULL;
		mp->id = (char *)malloc(5 * sizeof(char) );
		if(mp->id != NULL){
			if(name != NULL){
				mp->id[0] = '(';
				mp->id[1] = tolower(name[0]);
				mp->id[2] = tolower(name[1]);
				mp->id[3] = ')';
				mp->id[4] = '\0';
			}else{
				mp->id[0] = '\0';
			}
		}else{
			fprintf(stderr, "\nERROR: Out of memory - cannot create mod!\n");
		}
	}else{
		fprintf(stderr, "\nERROR: Out of memory - cannot create mod!\n");
	}
	return mp;
}


ModPointer addMod(ModPointer root, char *name, double mass){
	if(root == NULL){
		return newMod(name, mass);
	}else if(root->next == NULL){
		root->next = newMod(name, mass);
	}else{
		addMod(root->next, name, mass);
	}
	return root;
}


ModPointer delMod(ModPointer mp){
	if(mp == NULL){
		return NULL;
	}
	ModPointer next = mp->next;
	if(mp->id != NULL){
		free(mp->id);
		mp->id = NULL;
	}
	free(mp);
	mp = NULL;
	
	return next;	
}


ModPointer delMods(ModPointer mp){
	if(mp == NULL){
		return NULL;
	}
	while(mp != NULL){
		mp = delMod(mp);
	}
	return mp;
}


ModPointer getMods(xmlNodePtr searchSummary){
	ModPointer mp = NULL;
	xmlNodePtr child;
	for(child=searchSummary->children; child; child = child->next){
		if(child->type == XML_ELEMENT_NODE && (
			!xmlStrcmp(child->name, (const xmlChar *)AMINOACID_MOD) ||
			!xmlStrcmp(child->name, (const xmlChar *)TERMINAL_MOD) ) ){
		
			xmlChar *isVariable = getAttribute( child, VARIABLE);
			xmlChar *mass = getAttribute( child, MASS);
			xmlChar *description = getAttribute( child, DESCRIPTION);
			
			if(!xmlStrcmp(isVariable, (const xmlChar *)YES)){
				mp = addMod(mp, (char*)description, atof((char*)mass));
			}else{
				mp = addMod(mp, NULL, atof((char*)mass));
			}

			xmlFree(mass);
			xmlFree(description);
			xmlFree(isVariable); 
		}
	}
	return mp;
}


char *convertMod(ModPointer mp, double mass){

	char *id = NULL;
	while(mp != NULL){
		if(fabs(mp->mass - mass) < 0.0001){ //incase header/text values differ
			id = mp->id;
			break;
		}
		mp = mp->next;
	}
	if(id == NULL){
		fprintf(stderr, "\nERROR: Unkown modification mass: %f!\n", mass);
	}	
	return id;
}


char *modSequence(xmlNodePtr hit, xmlChar *seq, ModPointer mp){
	char *mods[xmlStrlen(seq)+2];
	int i;
	for(i = 0; i < xmlStrlen(seq)+2; ++i){
		mods[i] = NULL;
	}

	xmlNodePtr modInfo = downTo(hit, MOD_INFO);
	if(modInfo){
		xmlChar *nterm = getAttribute( modInfo, MOD_NTERM);
		if(nterm){
			mods[0] = convertMod(mp, atof((char*)nterm));
		}
		xmlFree(nterm);
				
		xmlChar *cterm = getAttribute( modInfo, MOD_CTERM);
		if(cterm){
			mods[xmlStrlen(seq)+1] = convertMod(mp, atof((char*)cterm));
		}
		xmlFree(cterm);

		xmlNodePtr mod;
		for(mod=modInfo->children; mod; mod = mod->next){
		if(mod->type == XML_ELEMENT_NODE &&
			!xmlStrcmp(mod->name, (const xmlChar *)MOD_AMINOACID) ){
				xmlChar *position = getAttribute( mod, POSITION);
				xmlChar *mass = getAttribute( mod, MASS);

				mods[atoi((char*)position)] = 
					convertMod(mp, atof((char*)mass));

				xmlFree(position);
				xmlFree(mass);
			}
		}
	}

	char modSeqBuf[127] = "";
	for(i = 0; i < xmlStrlen(seq)+2; ++i){
		if(mods[i]){
			strncat(modSeqBuf, mods[i], strlen(mods[i]) );
		}
		if(i < xmlStrlen(seq)){
			strncat(modSeqBuf, (char*)seq+i, 1 );
		}
	}
	char *modSeq = (char*)malloc( (strlen(modSeqBuf)+1)*sizeof(char) );
	strncpy(modSeq, modSeqBuf, strlen(modSeqBuf) );
	modSeq[strlen(modSeqBuf)] = '\0';

	return modSeq;
}


PeptidePointer getSearchResults(PeptidePointer pp, xmlNodePtr msmsRunSummary,
	ModPointer mp, char *mzXMLname){

	xmlNodePtr child = NULL;
	for(child=msmsRunSummary->children; child; child = child->next){
		if(child->type == XML_ELEMENT_NODE &&
			!xmlStrcmp(child->name, (const xmlChar *)SPECTRUM_QUERY) ){

			xmlChar *scan = getAttribute( child, START_SCAN);
			xmlNodePtr searchResult = downTo(child, SEARCH_RESULT);			
			xmlNodePtr hit;
			for(hit=searchResult->children; hit; hit = hit->next){
				if(hit->type == XML_ELEMENT_NODE &&
					!xmlStrcmp(hit->name, (const xmlChar *)SEARCH_HIT) ){
					
					xmlChar *rank = getAttribute( hit, HIT_RANK);
					if(atoi( (char*)rank) == 1){
						xmlChar *seq = getAttribute( hit, PEPTIDE);
						char *modSeq = modSequence(hit, seq, mp);
						pp = addPeptide(pp, mzXMLname, atoi((char*)scan),
								modSeq);
						free(modSeq);
						xmlFree(seq); 
					}				
					xmlFree(rank); 
				}
			}
			xmlFree(scan); 
		}
	}

	return pp;
}  


PeptidePointer readPepXML(char *filename, PeptidePointer pp){

	xmlDocPtr doc;
	xmlNodePtr cur;

	/*try to open pepXML file*/
	openXML(&doc, filename);

	/*get properties of pepXML as a whole*/
	cur = xmlDocGetRootElement(doc); 

	/*get to msms_run_summary*/
	cur = downTo(cur, MSMS_RUN_SUMMARY);
	
	/*convert base to mzXML filename*/
	xmlChar * base = getAttribute(cur, BASE_NAME);
	char *mzXMLname = NULL;

	/*set offset to be one past the index of last path seperator*/
	char * path = (char*)base;
	int offset = 0;
    if( (offset = (strrchr(path, '\\') - path + 1) ) < 0 &&
        (offset = (strrchr(path, '/') - path + 1)) < 0 ) { 
		offset = 0; 
	}

	if(strstr((char*)base, SEQ_MXZML_DIR)){\
		/*if pepxml came from out2xml filename already has .mzXML_dta in name*/
		size_t filename_length = xmlStrlen(base) - 4 + 1 - offset;
		mzXMLname = (char*)malloc(filename_length*sizeof(char));
		strncpy(mzXMLname, (char*)base+offset, filename_length-4);
		mzXMLname[filename_length-4] = '\0';
	}else{
		size_t filename_length = xmlStrlen(base) - offset;
		mzXMLname = (char*)malloc(
			(filename_length + strlen(MZXML_EXT) + 1)*sizeof(char));
		strncpy(mzXMLname, (char*)base+offset, filename_length + 1);
		strncpy(mzXMLname+filename_length, MZXML_EXT, strlen(MZXML_EXT)+1);
	}
	xmlFree(base);

	/*get to search_summary and compile mod list*/
	xmlNodePtr searchSummary = downTo(cur, SEARCH_SUMMARY);
	ModPointer mp = getMods(searchSummary);
	
	/*populate */	
	pp = getSearchResults(pp, cur, mp, mzXMLname);
	puts(mzXMLname);
	free(mzXMLname);
	
	/*cleanup*/
    xmlFreeDoc(doc); 
	xmlCleanupParser();
	mp = delMods(mp);

	return pp;
}


PeptidePointer parsePepXMLdir(char *dirName){

	DIR *dp;
	dp = opendir(dirName);
	PeptidePointer pp = &TNILL;	
	
	if(dp != NULL){
		struct dirent *ep = NULL;
		while ( (ep = readdir (dp)) ){
			if(strstr(ep->d_name, PEPXML_EXT)){
				puts(ep->d_name);
				char fullPath[strlen(dirName) + strlen(ep->d_name) + 2];
				strcpy(fullPath, dirName);
				strcat(fullPath, "/");
				strcat(fullPath, ep->d_name);
				pp = readPepXML(fullPath, pp);
			}
		}
		(void) closedir (dp);
	}else{
		fprintf(stderr, "\nERROR: Could not open directory: %s!\n", dirName);
		exit(1);
	}
	
	return pp;
}

