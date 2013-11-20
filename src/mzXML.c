/*
 * mzXML.c                                                                   
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
 * Functions for parsing and manipulating mzXML files containing scan        
 *     information for LC-MS run. Current version of file assumes peak lists 
 *     are uncompressed, use 32 bit precision, and network byte order.       
 *                                                                           
 * +ASCII art via Text ASCII Art Generator by Patrick Gillespie              
 */

#include "mzXML.h"
#include "base64.h" //decodeQuartet
#include "xml.h" //openXML, searchForXPath

#include <malloc.h> //malloc_trim
#include <libxml/tree.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>

#include <arpa/inet.h> //ntohl

#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <stdlib.h>

/*
 * U32 - Used for converting from uint32_t to float. Idea stolen from RAMP:
 *     Random Access Minimal Parser. Code available at:
 *     http://code.google.com/p/jmzreader/wiki/MzXmlParser                 
 */
typedef union {
   uint32_t u32;
   float flt;
} U32;

/*
 * newScan - Allocate memory for a new scan and intialize with passed values.
 *     Return a pointer to new scan, NULL if error occured. 
 */
ScanPointer newScan(int scanNum, int msLevel, int peaksCount, 
	float retentionTime, float totalIonCurrent, float precIntensity,
	int precCharge, float precMz, float *mzList, float *intList);

/*
 * delScan - Free all memory allocated for the scan. 
 */
ScanPointer delScan(ScanPointer sp);

/*
 * newMZXML - Allocate memory for a new mzXML and intializes with passed 
 *     values. Memory is allocated for pointers to scans but scans themselves
 *     are not created. Return a pointer to new mzXML, NULL if error occured.
 */
MZXMLPointer newMZXML(char *filename, int scanCount);

/*
 * getProperties - Parse properties relating to entire mzXML such as scanCount.                                                *
 */
void getProperties(xmlNodePtr cur, int *scanCount);

/*
 * parseScans - Extract all required information from opened mzXML file.
 */
int parseScans(xmlDocPtr doc, int scanCount, MZXMLPointer *mzXML);

/*
 * getScanAttributes - Get attributes of the <scan> nodes.
 */
int getScanAttributes(xmlNodePtr node, int *scanNum, int* peaksCount,
	int *msLevel, float *retentionTime, float *totalIonCurrent);

/*
 * getPeaksAttributes - Get attributes of the <peak> node and also decode and
 *     convert peak list. 
 */
int getPeaksAttributes(xmlNodePtr node, int msLevel, int peaksCount,
	float *precMz, float *precIntensity, int *precCharge, float **mzList,
	float **intList);

/*
 * checkCompatability - Ensure peak list is uncompressed, 32bit and in network
 *     byte order as other format are currently unsupported.
 */
int checkCompatibility(xmlNodePtr cur);

///////////////////////////////////////////////////////////////////////////////
//                       END OF FUNCTION DECLARATIONS                        //
///////////////////////////////////////////////////////////////////////////////

ScanPointer newScan (int scanNum, int msLevel, int peaksCount, 
	float retentionTime, float totalIonCurrent, float precIntensity,
	int precCharge, float precMz, float *mzList, float *intList){

	ScanPointer sp = (ScanPointer)malloc(sizeof(Scan));
	if(sp == NULL){
		 fprintf(stderr, "\nERROR: Out of memory - cannot create scan!\n");
	}else{
		sp->scanNum = scanNum;
		sp->msLevel = msLevel;
		sp->peaksCount = peaksCount;
		sp->retentionTime = retentionTime;
		sp->totalIonCurrent = totalIonCurrent;
		sp->precIntensity = precIntensity;
		sp->precCharge = precCharge;
		sp->precMz = precMz;
		sp->mzList = mzList;
		sp->intList = intList;
	}
	return sp;
}


ScanPointer delScan(ScanPointer sp){
	if(sp == NULL){
		return NULL;
	}
	if(sp->mzList != NULL){
		free(sp->mzList);
		sp->mzList = NULL;
	}
	if(sp->intList != NULL){
		free(sp->intList);
		sp->intList = NULL;
	}
	free(sp);
	sp = NULL;
	return sp;
}


MZXMLPointer newMZXML(char *filename, int scanCount){
	MZXMLPointer mp = (MZXMLPointer)malloc(sizeof(MZXML));
	if(mp == NULL){
		fprintf(stderr, "\nERROR: Out of memory - cannot create mzXML!\n");
	}else{
		mp->filename = (char *)malloc( (strlen(filename)+1)*sizeof(char) );
		mp->scanCount = scanCount;
		mp->scans = (ScanDB)malloc(mp->scanCount * sizeof(ScanPointer));
		if(mp->scans == NULL || mp->filename == NULL){
			fprintf(stderr, "\nERROR: Out of memory - cannot create mzXML!\n");
			free(mp);
			mp = NULL;
		}else{
			strncpy(mp->filename, filename, strlen(filename));	
			mp->filename[strlen(filename)] = '\0';
		}
	}
	return mp;
}


void getProperties(xmlNodePtr cur, int *scanCount){
	for(cur=cur->children; cur; cur = cur->next){
		/*find and record scanCount attribute of <msRUN> node*/
		if(cur->type == XML_ELEMENT_NODE && 
			!strcmp( (char *)cur->name, "msRun")){
			xmlAttr* attribute = cur->properties;

			while(attribute && attribute->name && attribute->children){
				if(!strcmp((char *)attribute->name, "scanCount")){
					xmlChar* value = xmlNodeListGetString(cur->doc,
						attribute->xmlChildrenNode, 1);
					*scanCount = atoi((char *)value);
					xmlFree(value); 
				}
				attribute = attribute->next;
			}
		}
	}
	return;
}


int parseScans(xmlDocPtr doc, int scanCount, MZXMLPointer *mzXML){
	/*initialize XPath environment*/
	xmlXPathInit();
	xmlXPathContextPtr context;
	xmlXPathObjectPtr result;
	xmlNodeSetPtr nodeset = NULL;
	int status = 0;
	int i;
		
	/*run search for scan nodes*/
	status = searchForXPath(&context, &result, &doc,
		(xmlChar*) "//*[local-name()='scan']");

	/*check for consistency with mzXML <msRUN> node*/
	if(status == 0){	
		nodeset = result->nodesetval;
		if(scanCount != nodeset->nodeNr){
			fprintf(stderr,"\nERROR: Number of scans reported by mzXML(%d) ",
				scanCount);
			fprintf(stderr,"does not agree with the number found(%d).\n",
				nodeset->nodeNr);
			status = -1;
		}
	}

	for(i=0; i < nodeset->nodeNr; ++i){
		xmlNodePtr node = nodeset->nodeTab[i];
				
		int peaksCount = 0;
		int msLevel = 0;
		int scanNum = 0;
		float retentionTime = 0;
		float totalIonCurrent = 0;
		float precMz = 0;
		float precIntensity = 0;
		int precCharge = 0;
		float *mzList = NULL;
		float *intList = NULL;

		status = getScanAttributes(node, &scanNum, &peaksCount, &msLevel,
			&retentionTime, &totalIonCurrent);
		if(status == 0 || status == EMPTY_PEAK_LIST){
			status = getPeaksAttributes(node, msLevel, peaksCount, &precMz,
				&precIntensity, &precCharge, &mzList, &intList);
		}
		if(status == 0){
			(*mzXML)->scans[i] = newScan(scanNum, msLevel, peaksCount,
				retentionTime, totalIonCurrent, precIntensity, precCharge,
				precMz, mzList, intList);
		}
	}

	/*cleanup XPath Environment*/
 	xmlXPathFreeObject(result);
    xmlXPathFreeContext(context);

	return status;
}


int getScanAttributes(xmlNodePtr node, int *scanNum, int* peaksCount,
	int *msLevel, float *retentionTime, float *totalIonCurrent){
	int status = 0;
	xmlAttr *attribute = node->properties;
	/*find and record num, peaksCount, msLevel, and retentionTime attributes
	of <scan> node*/
	while(attribute && attribute->name && attribute->children){
  		xmlChar* value = xmlNodeListGetString(node->doc,
			attribute->xmlChildrenNode, 1);
		if(!strcmp((char *)attribute->name, "num")){
			*scanNum = atoi((char*)value);
		}else if(!strcmp((char *)attribute->name, "peaksCount")){
			*peaksCount = atoi((char*)value);
			if(*peaksCount < MIN_PEAK_COUNT){
				status = EMPTY_PEAK_LIST;
			}
		}else if(!strcmp((char *)attribute->name, "msLevel")){
			*msLevel = atoi((char*)value);
		}else if(!strcmp((char *)attribute->name, "retentionTime")){
			sscanf((char *)value,"PT%fS", retentionTime);
		}else if(!strcmp((char *)attribute->name, "totIonCurrent")){
			*totalIonCurrent = atof((char*)value);
		}
  		xmlFree(value); 
		attribute = attribute->next;
	}
	return status;
}


int getPeaksAttributes(xmlNodePtr node, int msLevel, int peaksCount,
	float *precMz, float *precIntensity, int *precCharge, float **mzList,
	float **intList){

	xmlNodePtr cur;
	/*find and record precursorMz, precursorCharge, and precursorIntensity of
	<precursorMz> node*/
	for(cur=node->children; cur; cur = cur->next){
		if(cur->type == XML_ELEMENT_NODE){
			if(msLevel == 2 && !strcmp( (char *)cur->name, "precursorMz") ){
				xmlAttr *attribute = cur->properties;
				while(attribute && attribute->name && attribute->children){
					xmlChar* value = xmlNodeListGetString(cur->doc,
						attribute->xmlChildrenNode, 1);
					if(!strcmp((char *)attribute->name, "precursorCharge")){
						*precCharge = atoi((char*)value);
					}else if(!strcmp((char *)attribute->name,
						"precursorIntensity")){
						*precIntensity = atof((char*)value);
					}
					xmlFree(value); 
					attribute = attribute->next;
				}
			*precMz = atof((char*)cur->children->content);
			
			}else if(!strcmp( (char *)cur->name, "peaks") &&
				peaksCount > MIN_PEAK_COUNT){

				if(checkCompatibility(cur) == 1){
					/*decode base64 encoded peak list*/
					unsigned char * encodedList =
						(unsigned char *) cur->children->content;
					int encodedLength = strlen((char*)encodedList);
					int padding =
						(encodedList[encodedLength-2] == '=')? 2 :
						(encodedList[encodedLength-1] == '=')? 1 : 0;
					int decodedLength = (encodedLength/4)*3-padding;
					//char unsigned decodedList[(encodedLength/4)*3];
					char unsigned decodedList[decodedLength];
					int i,j;
					for(i = 0, j = 0; i < encodedLength; i+=4, j+=3){
						decodeQuartet(encodedList+i, decodedList+j);
					}
					*mzList = (float *)malloc(peaksCount*sizeof(float));
					*intList = (float *)malloc(peaksCount*sizeof(float));

					/*convert from network byte order*/
					for(i = 0; i < peaksCount; ++i){
						U32 tmpA, tmpB;
						tmpA.u32 = ntohl(((uint32_t *)decodedList)[i*2]);
						tmpB.u32 = ntohl(((uint32_t *)decodedList)[i*2+1]);
						(*mzList)[i] = tmpA.flt;
						(*intList)[i] = tmpB.flt;
					}
				}
			}
		}
	}
	return 0;
}


int checkCompatibility(xmlNodePtr cur){
	xmlAttr *attribute = cur->properties;
	while(attribute && attribute->name && attribute->children){
		xmlChar* value = xmlNodeListGetString(cur->doc,
			attribute->xmlChildrenNode, 1);
		if(!strcmp((char *)attribute->name, "precision")){
			if( atoi((char *)value) != 32 ){
				printf("ERROR: Precision does not equal 32\n");
				return -1;
			}
		}else if(!strcmp((char *)attribute->name, "compressionType")){
			if(strcmp((char *)value, "none")){
				printf("ERROR: CompressionType does not equal \"none\"\n");
				return -1;
			}
		}else if(!strcmp((char *)attribute->name, "byteOrder")){
			if(strcmp((char *)value, "network")){
				printf("ERROR: ByteOrder does not equal \"network\"\n");
				return -1;
			}
		}
		xmlFree(value); 
		attribute = attribute->next;
	}
	return 1;
}


MZXMLPointer delMZXML(MZXMLPointer mp){
	if(mp == NULL){
		return NULL;
	}
	if(mp->filename != NULL){
		free(mp->filename);
	}
	if(mp->scans != NULL){
		int i;
		for(i = 0; i < mp->scanCount; ++i){
			if(mp->scans[i] != NULL){
				mp->scans[i] = delScan(mp->scans[i]);
			}
		}
		free(mp->scans);
		mp->scans = NULL;
	}
	free(mp);
	mp = NULL;
	return mp;
}


int readMZXML(char *filename, MZXMLPointer *mzXML ){
	
	xmlDocPtr doc;
	xmlNodePtr cur;
	int status; //tracks success/failure of successive function calls
	int scanCount = 0;

	/*try to open mzXML file*/
	status = openXML(&doc, filename);

	/*get properties of mzXML as a whole*/
	if(status == 0){	
		cur = xmlDocGetRootElement(doc);
		getProperties(cur, &scanCount);
		if(scanCount == 0){
			fprintf(stderr, "\nERROR: File: %s had 0 scans.\n", filename);
			status = -1;
		}
	}

	*mzXML =  newMZXML(filename, scanCount);

	/*get properties of specific runs and populate mzXML struct*/
	if(status == 0){
		status = parseScans(doc, scanCount, mzXML);
		if(status != 0){
			fprintf(stderr, "\nERROR: Parsing file: %s.\n", filename);
		}
	}
 	
    xmlFreeDoc(doc); 
	xmlCleanupParser();
	malloc_trim(0);
	
	return status;
}

