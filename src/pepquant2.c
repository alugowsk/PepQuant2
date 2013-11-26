/*
 * pepquant2.c                                    
 * ===========                        
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
 * The pepquant2.c main source file. Coordinates work done by all other files.
 *     Contains external variables and default settings for them. 
 *                                                                           
 * +ASCII art via Text ASCII Art Generator by Patrick Gillespie              
 */

#include "peptide.h"
#include "protein.h"
#include "fasta.h"
#include "common.h"
#include "pepxml.h"
#include "ls.h"
#include "global.h"

#include <math.h>//fabs, NAN
#include <stdlib.h> 
#include <stdio.h> //printf
#include <string.h> //strcmp
#include <stdbool.h>

int alignWindow = 90; //max discrepancy between MS1 and MS2 rt before defaulting
					  //to MS2 rt
double corrCutOff = 0.99; //minimum correlation between theoretical and observed
						  //isotopic patterns
int peakWindow = 30; //max distance for adjacent peak to accept apical MS1 rt
char *fastaName = NULL; //path of FASTA used for search
double intCutOff = 1E6; //minimum intensity for a 
						//(spectra, peptide, charge state tuple)
int maxCharge = 4; //highest charge state to look for in MS1
double ppmCutOff = 0.000010; //ppm error tolerance
int isotopicStates = 4; //number of isotopic states to look for in MS1 spectra
int threadCount = 1; //number of threads to use for multi-threaded steps
int quantWindow = 150; //half width of quantification window
char *statQuestdir = NULL; //parse statQuest results (statQuest directory)
int statQuestcutoff = 0; //confidence filter to use for statQuest results
char *pepXMLdir = NULL; //parse pepXML results (directory with .pepXML files)
char *maxQuant = NULL; //parse maxQuant results (msms.txt)
int lys = 0; //SILAC label status 
int arg = 0; //SILAC label status
bool ignoreModSite = false; //ignore modifcation localization


/*
 * genMS2rtTable - generate a ms2 retention time table. The current 
 *     implementation uses the weighted centroid (average of retention times
 *     weighted by their respective intensities.
 */
void genMS2rtTable(PeptidePointer *peptides, SpectraFileNodePointer filelist,
	int peptideCount, int fileCount, double **ms2rt);

/*
 * genMS1rtTable - generate a ms1 retention time table. The current 
 *     implementation uses the apical retention time (retention time with
 *     highest intensity).
 */
void genMS1rtTable(PeptidePointer *peptides, SpectraFileNodePointer filelist,
	int peptideCount, int fileCount, double **ms1rt);

/*
 * medianRTtimes - given both ms1 and ms2 retention time tables create vectors
 *     of median retention time for each, where the median is determined over
 *     all non-zero values. 
 */
void medianRTtimes(double *ms1median, double *ms2median, double **ms1rt,	
	double **ms2rt, int peptideCount, int fileCount, PeptidePointer *peptides);

/*
 * quant - given a retention time table, an array of peptides, and a
 *     linked list of files, create a quantity entry for every peptide, file
 *     tuple and store the value in the quantification table. 
 */
void quant(double **ms2rt, double **quantification, PeptidePointer *peptides,
	int peptideCount, SpectraFileNodePointer filelist);

/*
 * countSpectra - create a table of ms2 spectral counts at the peptide level. 
 *     Current implementation uses a table of double to easily integrate with
 *     existing function for converting peptide level quantification results to
 *     the protein level.
 */
void countSpectra(PeptidePointer *peptides, int peptideCount,
	SpectraFileNodePointer filelist, double **spectralCounts);

///////////////////////////////////////////////////////////////////////////////
//                       END OF FUNCTION DECLARATIONS                        //
///////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]){

	parseArgs(argc, argv);

	/*read FASTA file*/
	printf("Reading fasta %s\n", fastaName);
	FastaPointer fasta;
	readFASTA(fastaName, &fasta );

	/*parse input file(s) and init peptide list*/
	PeptidePointer pp = &TNILL;
	if(maxQuant){
		printf("Parsing MaxQuant file.\n");
		pp = parseMaxQuant(maxQuant, pp);
	}else if(pepXMLdir){
		printf("Parsing pepXML directory.\n");
		pp = parsePepXMLdir(pepXMLdir);
	}else if(statQuestdir){
		printf("Parsing statQuest directory.\n");
		pp = parseStatQuestdir(statQuestdir, statQuestcutoff);
	}else{
		fprintf(stderr, "\nERROR: no input file(s) specified!\n");
		exit(EXIT_FAILURE);
	}

	/*get peptides list and generate isotopic patterns*/
	int peptideCount = getCount(pp);
	printf("\tFound %d peptides!\n", peptideCount);
	PeptidePointer *peptides =  inOrder(pp, peptideCount);
	printf("Generating isotopic patterns.\n");
	pp = initIsotopicPatterns(pp, peptides, peptideCount);
	peptideCount = getCount(pp);
	printf("\tGenerated %d isotopic patterns!\n", peptideCount);
	free(peptides); //memory will be reallocated if inOrder is called again
	peptides =  inOrder(pp, peptideCount);

	printPeptides(peptides, peptideCount);

	/*get mzXML filelist*/
	printf("Generating mzXML filelist.\n");
	SpectraFileNodePointer filelist = NULL;
	int fileCount = initFilelist(&filelist, peptides, peptideCount);

	/*search mzXMLs for isotopic patterns*/
	printf("Searching ms1 spectra:\n");	
	searchMzXMLs(peptides, peptideCount, filelist);
	
	/*print results*/
	printf("Printing search results.\n");
	printSearchResults(peptides, peptideCount);

	/*calculate ms2 table*/
	printf("Calculating ms2 rt table\n");
	double **ms2rt = new2Darray(peptideCount, fileCount);
	genMS2rtTable(peptides, filelist, peptideCount, fileCount, ms2rt);
		
	/*print ms2 retention time table*/
	printf("Printing ms2 rt table.\n");
	printTable(peptides, filelist, peptideCount, fileCount, ms2rt, "ms2rt.txt");

	/*calculate ms1 table*/
	printf("Calculating ms1 rt table\n");
	double **ms1rt = new2Darray(peptideCount, fileCount);
	genMS1rtTable(peptides, filelist, peptideCount, fileCount, ms1rt);

	/*print ms1 retention time table*/
	printf("Printing ms1 rt table.\n");
	printTable(peptides, filelist, peptideCount, fileCount, ms1rt, "ms1rt.txt");
	
	/*calculate and print median ms1 and ms2 retention times*/
	printf("Calculating median retention times.\n");
	double *ms1median = (double *)malloc(peptideCount * sizeof(double) );
	double *ms2median = (double *)malloc(peptideCount * sizeof(double) );
	medianRTtimes(ms1median, ms2median, ms1rt, ms2rt, peptideCount, fileCount,
		peptides);

	/*determine ms2 retention time callibration*/
	printf("Determining ms2 retention time calibration.\n");
	double **ms2params = 
		leastSquares(ms2rt, ms2median, peptideCount, fileCount);

	/*determine ms1 retention time callibration*/
	printf("Determining ms1 retention time calibration.\n");
	double **ms1params = 
		leastSquares(ms1rt, ms1median, peptideCount, fileCount);

	/*print alignment parameters*/
	FILE *fp = fopen("ms2params.txt", "w");
	if (fp == NULL)	{
   		fprintf(stderr, "Error opening file: %s !\n", "ms2params.txt");
	}
	FILE *fp2 = fopen("ms1params.txt", "w");
	if (fp2 == NULL)	{
   		fprintf(stderr, "Error opening file: %s !\n", "ms1params.txt");
	}
	int i;
	SpectraFileNodePointer root = filelist;
	for(i=0; i<fileCount; ++i){
		fprintf(fp, "%s\t%.4e\t%.4e\n", root->rawFile, 
			ms2params[i][0], ms2params[i][1]);
		fprintf(fp2, "%s\t%.4e\t%.4e\n", root->rawFile, 
			ms1params[i][0], ms1params[i][1]);
		root = root->next;	
	}
	fclose(fp);
	fclose(fp2);

	/*align*/
	printf("Aligning.\n");
	align(ms2rt, ms1rt, ms2median, ms1median, ms2params, ms1params,
		peptideCount, fileCount);
	
	/*print callibrated rt table*/
	printf("Printing aligned rt table.\n");
	printTable(peptides, filelist, peptideCount, fileCount, ms2rt, "rt.txt");

	/*calculate intensities*/
	printf("Calculating intensities.\n");
	double **quantification = new2Darray(peptideCount, fileCount);
	quant(ms2rt, quantification, peptides, peptideCount, filelist);

	/*determine spectral counts table*/
	double **spectralCounts = new2Darray(peptideCount, fileCount);
	countSpectra(peptides, peptideCount, filelist, spectralCounts);

	/*convert protein information to proteins*/
	printf("Converting peptide information to protein level.\n");
	ProteinNodePointer pnp = NULL;

	char **proteinMap = NULL;
	int proteinCount = pep2prot(&pnp, fasta, peptides, quantification,
		spectralCounts, peptideCount, fileCount, &proteinMap);
	double **protQuant = new2Darray(proteinCount, fileCount);
	double **protLQuant = new2Darray(proteinCount, fileCount);
	double **protHQuant = new2Darray(proteinCount, fileCount);
	double **protSpectralCounts = new2Darray(proteinCount, fileCount);
	prot2table(pnp, &protQuant, &protHQuant, &protLQuant, &protSpectralCounts,
		proteinCount, fileCount);

	/*print quantitation results*/
	printf("Printing intensities and spectral counts.\n");
	printPepTable(peptides, filelist, peptideCount, fileCount, quantification,
		"quant.txt", proteinMap);
	printPepTable(peptides, filelist, peptideCount, fileCount, spectralCounts,
		"pepSpectra.txt", proteinMap);

	/*print protein quantitation results*/
	printf("Printing protein intensities.\n");
	printProtTable(pnp, filelist, proteinCount, fileCount, protQuant,
		"protQuant.txt");
	printProtTable(pnp, filelist, proteinCount, fileCount, protHQuant,
		"protHQuant.txt");
	printProtTable(pnp, filelist, proteinCount, fileCount, protLQuant,
		"protLQuant.txt");
	printProtTable(pnp, filelist, proteinCount, fileCount, protSpectralCounts,
		"protSpectra.txt");

	/* save some mem and determine H/L table in place */
	int j;
	for(i = 0; i < proteinCount; ++i){
		for(j = 0; j < fileCount; ++j){
			protHQuant[i][j] = (protLQuant[i][j] == 0)? NAN :
				protHQuant[i][j]/protLQuant[i][j];
		}
	}
	printProtTable(pnp, filelist, proteinCount, fileCount, protHQuant,
		"protHLratio.txt");

	printCoverage(fasta);
	
	/*Cleanup*/
	printf("Cleaning up.\n");
	proteinMap = delProteinMap(proteinMap, peptideCount);
	pp = delPeptideList(pp);
	free(peptides);
	pnp = delProteinList(pnp);
	filelist =  delSpectraFileList(filelist);
	delFasta(fasta);
	del2Darray(ms2rt, peptideCount);
	del2Darray(ms1rt, peptideCount);
	del2Darray(quantification, peptideCount);
	del2Darray(spectralCounts, peptideCount);
	del2Darray(protQuant, proteinCount);
	del2Darray(protLQuant, proteinCount);
	del2Darray(protHQuant, proteinCount);
	del2Darray(protSpectralCounts, proteinCount);
	free(ms2median);
	free(ms1median);

	return 0;
}


void genMS2rtTable(PeptidePointer *peptides, SpectraFileNodePointer filelist,
	int peptideCount, int fileCount, double **ms2rt){

	int i;
	int j = 0;
	SpectraFileNodePointer rootFilelist = filelist;
	for(i = 0; i < peptideCount; ++i){
		filelist = rootFilelist;
		j = 0;
		while(filelist != NULL){
			SpectraFileNodePointer sfnp = peptides[i]->spectraFiles;
			double totalIntensity = 0;
			double weightedCentroid = 0;
			while(sfnp != NULL){	
				if(!strcmp(sfnp->rawFile, filelist->rawFile) ){
					ScanNodePointer snp = sfnp->scans;
					while(snp != NULL){
						totalIntensity+=snp->intensity[0];
						snp = snp->next;
					}
					snp = sfnp->scans;
					while(snp != NULL){
						weightedCentroid+=
							(snp->intensity[0]*snp->rt)/totalIntensity;
						snp = snp->next;
					}
				}
				sfnp = sfnp->next;
			}
			ms2rt[i][j] = weightedCentroid;
			filelist = filelist->next;	
			j++;
		}
	}
	return;
}


void genMS1rtTable(PeptidePointer *peptides, SpectraFileNodePointer filelist,
	int peptideCount, int fileCount, double **ms1rt){

	int i;
	int j = 0;
	SpectraFileNodePointer rootFilelist = filelist;
	for(i = 0; i < peptideCount; ++i){
		filelist = rootFilelist;
		j = 0;
		while(filelist != NULL){
			SpectraFileNodePointer sfnp = peptides[i]->ms1SpectraFiles;
			double maxIntensity = 0;
			double maxRT = 0;
			int valid = 0;
			while(sfnp != NULL){	
				if(!strcmp(sfnp->rawFile, filelist->rawFile) ){
					ScanNodePointer snp = sfnp->scans;
					while(snp != NULL){
						int k;
						double totalIntensity = 0;
						for(k = 0; k < maxCharge - MIN_CHARGE + 1; ++k){
							totalIntensity+=snp->intensity[k];
						}
						if(totalIntensity > maxIntensity){
							maxIntensity = totalIntensity;
							maxRT = snp->rt;
						}
						snp = snp->next;
					}
					snp = sfnp->scans;
					while(snp != NULL){
						if(fabs(snp->rt - maxRT) < 
							peakWindow && snp->rt != maxRT){

							valid = 1;
						}
						snp = snp->next;
					}
				}
				sfnp = sfnp->next;
			}
			ms1rt[i][j] = (valid == 0)? 0 : maxRT;
			filelist = filelist->next;	
			j++;
		}
	}
	return;
}


void medianRTtimes(double *ms1median, double *ms2median, double **ms1rt,	
	double **ms2rt, int peptideCount, int fileCount, PeptidePointer *peptides){
	
	/*Calculate ms1 and ms2 retention times*/
	int i;
	for(i = 0; i < peptideCount; ++i){
		ms1median[i] = median(ms1rt[i], fileCount);
		ms2median[i] = median(ms2rt[i], fileCount);
	}

	/*print ms1 and ms2 median retention times*/
	FILE *fp = fopen("median.txt", "w");
	if (fp == NULL){
		printf("Error opening file!\n");
	}

	for(i = 0; i < peptideCount; ++i){
		fprintf(fp, "%s\t%8.4f\t%8.4f\n",peptides[i]->sequence, ms1median[i],
			ms2median[i]);
	}
	fclose(fp);

	return;
}

void quant(double **ms2rt, double **quantification, PeptidePointer *peptides,
	int peptideCount, SpectraFileNodePointer filelist){

	SpectraFileNodePointer rootFilelist = filelist;
	int i;
	int j = 0;

	for(i = 0; i < peptideCount; ++i){
		filelist = rootFilelist;
		j = 0;
		while(filelist != NULL){
			SpectraFileNodePointer sfnp = peptides[i]->ms1SpectraFiles;
			double totalIntensity = 0;
			while(sfnp != NULL){	
				if(!strcmp(sfnp->rawFile, filelist->rawFile) ){
					ScanNodePointer snp = sfnp->scans;
					while(snp != NULL){
						if(fabs(snp->rt - ms2rt[i][j]) < quantWindow){
							int k;
							int valid = 0;
							double total = 0;
							for(k = 0; k < maxCharge - MIN_CHARGE + 1; k++){
								if(snp->corr[k] >= corrCutOff){
									valid = 1;
								}
								total += snp->intensity[k];
							}
							totalIntensity+= valid == 1? total : 0;
						}
						snp = snp->next;
					}
					break;
				}
				sfnp = sfnp->next;
			}
			quantification[i][j] = totalIntensity;
			filelist = filelist->next;	
			j++;
		}
	}
	return;
}

void countSpectra(PeptidePointer *peptides, int peptideCount,
	SpectraFileNodePointer filelist, double **spectralCounts){

	SpectraFileNodePointer rootFilelist = filelist;
	int i;
	int j = 0;

	for(i = 0; i < peptideCount; ++i){
		filelist = rootFilelist;
		j = 0;
		while(filelist != NULL){
			SpectraFileNodePointer sfnp = peptides[i]->spectraFiles;
			double count = 0;
			while(sfnp != NULL){	
				if(!strcmp(sfnp->rawFile, filelist->rawFile) ){
					ScanNodePointer snp = sfnp->scans;
					while(snp != NULL){
						count+=1;
						snp = snp->next;
					}
					break;
				}
				sfnp = sfnp->next;
			}
			spectralCounts[i][j] = count;
			filelist = filelist->next;	
			j++;
		}
	}
	return;
}
