/*
 * xml.h                                                                     
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

#ifndef XML_H
#define XML_H

#include <libxml/tree.h>
#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <libxml/xpathInternals.h>

/*
 * openXML - Open and parse xml file. Store pointer to file.
 */
int openXML(xmlDocPtr *doc, char *filename);

/*
 * searchForXPath - Search the parsed xml for passed XPath.
 */
int searchForXPath(xmlXPathContextPtr *context, xmlXPathObjectPtr *result,
	xmlDocPtr *doc, xmlChar *xpath);

/*
 * downTo - move a level down in the xml node heirarchy and return a xmlNodePtr
 *     to a node with targetNodeName if it exists, NULL otherwise.
 */
xmlNodePtr downTo (xmlNodePtr cur, char *targetNodeName);

/*
 * getAttribute - scan the children of passed node and look for an attribute
 *     with name attr. If found return the string value otherwise return NULL.
 */
xmlChar *getAttribute(xmlNodePtr node, char *attr);

#endif

