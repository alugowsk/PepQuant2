/*
 * xml.c                                                                     
 * =====                                                                     
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
 * Common functions for parsing and manipulating xml files.                  
 *                                                                           
 * +ASCII art via Text ASCII Art Generator by Patrick Gillespie              
 */

#include "xml.h"

int openXML(xmlDocPtr *doc, char *filename){
	
	if(filename == NULL){
		fprintf(stderr,
			"\nERROR: Could not parse xml. Filename was NULL.\n");
		return -1;
	}

	*doc = xmlParseFile(filename);

	if (*doc == NULL) {
		fprintf(stderr, "\nERROR: File: %s not parsed successfully.\n",
			filename);
		return -1;
	}
	return 0;
}


int searchForXPath(xmlXPathContextPtr *context, xmlXPathObjectPtr *result,
		xmlDocPtr *doc, xmlChar *xpath){
	*context = xmlXPathNewContext(*doc);
 	if(*context == NULL) {
        fprintf(stderr, "ERROR: unable to create new XPath context.\n");
        xmlFreeDoc(*doc); 
        return -1;
    }
	*result = xmlXPathEvalExpression( xpath, *context);
	if(*result == NULL) {
		fprintf(stderr, "Error: unable to evaluate xpath expression \"%s\".\n",
			xpath);
		xmlXPathFreeContext(*context); 
		xmlFreeDoc(*doc); 
		return -1;
    }
	if(xmlXPathNodeSetIsEmpty((*result)->nodesetval)){
		xmlXPathFreeObject(*result);
		fprintf(stderr, "Error: no results from searching \"%s\"\n", xpath);
		return -1;
	}
	return 0;
}


xmlNodePtr downTo (xmlNodePtr cur, char *targetNodeName){
	cur = cur->xmlChildrenNode;
	while(cur != NULL){
		if ((!xmlStrcmp(cur->name, (const xmlChar *)targetNodeName))) {
			break;
		}
		cur = cur->next;
	}
	return cur;
}


xmlChar *getAttribute(xmlNodePtr node, char *attr){
	xmlAttr* attribute = node->properties;

	xmlChar* value = NULL;

	while(attribute && attribute->name && attribute->children){
		if(!xmlStrcmp(attribute->name, (const xmlChar *)attr) ){
			value = xmlNodeListGetString(node->doc,
				attribute->xmlChildrenNode, 1);
			break;
		}
		attribute = attribute->next;
	}
	return value;
}
