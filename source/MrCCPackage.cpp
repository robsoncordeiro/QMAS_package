/**********************************************************************
* GBDI Arboretum - Copyright (c) 2002-2009 GBDI-ICMC-USP
*
*                           Homepage: http://gbdi.icmc.usp.br/arboretum
**********************************************************************/
/* ====================================================================
 * The GBDI-ICMC-USP Software License Version 1.0
 *
 * Copyright (c) 2009 Grupo de Bases de Dados e Imagens, Instituto de
 * Ciências Matemáticas e de Computação, University of São Paulo -
 * Brazil (the Databases and Image Group - Intitute of Matematical and
 * Computer Sciences).  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 *
 * 3. The end-user documentation included with the redistribution,
 *    if any, must include the following acknowledgment:
 *       "This product includes software developed by Grupo de Bases
 *        de Dados e Imagens, Instituto de Ciências Matemáticas e de
 *        Computação, University of São Paulo - Brazil (the Databases
 *        and Image Group - Intitute of Matematical and Computer
 *        Sciences)"
 *
 *    Alternately, this acknowledgment may appear in the software itself,
 *    if and wherever such third-party acknowledgments normally appear.
 *
 * 4. The names of the research group, institute, university, authors
 *    and collaborators must not be used to endorse or promote products
 *    derived from this software without prior written permission.
 *
 * 5. The names of products derived from this software may not contain
 *    the name of research group, institute or university, without prior
 *    written permission of the authors of this software.
 *
 * THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED.  IN NO EVENT SHALL THE AUTHORS OF THIS SOFTWARE OR
 * ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
 * OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 * ====================================================================
 *                                            http://gbdi.icmc.usp.br/
 */
/**
* @file
* This file demonstrates the usage of the stCorrelationClustering class.
*
* @version 1.0
* @author Robson Leonardo Ferreira Cordeiro (robson@icmc.usp.br)
* @author Agma Juci Machado Traina (agma@icmc.usp.br)
* @author Christos Faloutsos (christos@cs.cmu.edu)
* @author Caetano Traina Jr (caetano@icmc.usp.br)
*/
// Copyright (c) 2002-2009 GBDI-ICMC-USP

//------------------------------------------------------------------------------
// MrCCProject
//------------------------------------------------------------------------------

#include <string.h>
#include <time.h>

#ifndef __GNUG__
	#include "stCorrelationClustering.h"
#endif //__GNUG__

#ifdef __GNUG__
	#include "stCorrelationClustering.cpp"
#endif //__GNUG__


// default values
#define TAM_LINE 500
#define NORMALIZE_FACTOR 0 // Independent

// global variables
clock_t startTime;

/**
* Initiates the measurement of run time.
*/
void initClock() {
	startTime=clock();
}//end initClock

/**
* Prints the elapsed time.
*/
void printElapsed() { 
	printf("Elapsed time: %0.3lf sec.\n",(double)(clock()-startTime)/CLOCKS_PER_SEC);
}//end printElapsed

/**
* Initiates the clustering process.
*
* @param argc Number of arguments received.
* @param argv Array with the arguments received.
* @return success (0), error (1).
*/
int main(int argc, char **argv) {
	initClock(); // initiates the meassurement of run time

	// first validations
	if (argc != 8) {
		printf("Usage: MrCCProject <pThreshold> <H> <hardClustering> <database_file> <result_file> <number_of_dimensions> <number_of_objects>");
		return 1; //error
	}//end if

	if (atoi(argv[2]) < 2) {
		printf("MrCC needs at least two resolution levels (H >= 2) to perform the clustering process.");
		return 1; //error
	}//end if
	
	// opens/creates the used files
	FILE *database, *result;
    database=fopen(argv[4], "r");
	result=fopen(argv[5], "w");
	if (!(database&&result)) {
		printf("MrCC could not open the database file or create the result file.");
		return 1; //error
	}//end if
	
	// reads objects from the source database
	int classId, numberOfDimensions=atoi(argv[6]), numberOfObjects=atoi(argv[7]);
	double **objectsArray = new double*[numberOfObjects];
	for (int i=0; i<numberOfObjects; i++) {
      objectsArray[i] = new double[numberOfDimensions];
	}//end for
	for (int i=0; i<numberOfObjects; i++) {
		for (int j=0; j<numberOfDimensions; j++) {
		   fscanf(database,"%lf",&objectsArray[i][j]);
		}//end for
		fscanf(database,"%d",&classId); // discarts the classes
    }//end for
	fclose(database); // the database file will not be used anymore

    // creates an object of the class stCorrelationClustering
    stCorrelationClustering *sCluster = new stCorrelationClustering(objectsArray, numberOfDimensions, NORMALIZE_FACTOR, numberOfObjects, 
																	(2*numberOfDimensions), -1, atof(argv[1]), atoi(argv[2]), atoi(argv[3]) );

    // looks for correlation clusters
    sCluster->findCorrelationClusters();

	// mounts the result file
	int numBetaClusters = sCluster->getNumBetaClusters(), numCorrelationClusters = sCluster->getNumCorrelationClusters(), betaCluster,
		*correlationClustersBelongings = sCluster->getCorrelationClustersBelongings();
	char **dimCorrelationClusters = sCluster->getDimCorrelationClusters(), line[TAM_LINE], belongsTo, str[20];
	double **minBetaClusters = sCluster->getMinBetaClusters(), **maxBetaClusters = sCluster->getMaxBetaClusters(),
		   *normalizeSlope = sCluster->getCalcTree()->getNormalizeSlope(),
	       *normalizeYInc = sCluster->getCalcTree()->getNormalizeYInc();

	// axes relevant to the found clusters
	for (int i=0; i<numCorrelationClusters; i++) {
		strcpy(line,""); // empty line
		strcat(line,"ClusterResult");
		sprintf(str,"%d",i+1);		
		strcat(line,str);		
		for (int j=0; j<numberOfDimensions; j++) {
			(dimCorrelationClusters[i][j]) ? strcat(line," 1") : strcat(line," 0");
		}//end for
		strcat(line,"\n");
		// fputs(line,result); // writes the relevant axes to the current cluster in the result file
	}//end for

	// labels each point after the clusters found
	// fputs("LABELING\n",result);
	int outlier;
	for (int point=0; point < numberOfObjects; point++) {
		outlier = 1;
		for (betaCluster=0; betaCluster < numBetaClusters; betaCluster++) {
			belongsTo=1;
			// undoes the normalization and verify if the current point belongs to the current beta-cluster
			for (int dim=0; belongsTo && dim<numberOfDimensions; dim++) {				
				if (! (objectsArray[point][dim] >= ((minBetaClusters[betaCluster][dim]*normalizeSlope[dim])+normalizeYInc[dim]) && 
					   objectsArray[point][dim] <= ((maxBetaClusters[betaCluster][dim]*normalizeSlope[dim])+normalizeYInc[dim])) ) {
					belongsTo=0; // this point does not belong to the current beta-cluster
				}//end if
			} // end for
			if (belongsTo) {
				outlier = 0;
				sprintf(str, "%d,%d\n", point+1, correlationClustersBelongings[betaCluster]+1);
				fputs(str, result);
			}//end if
		}//end for
		if (outlier) {
			sprintf(str, "%d,%d\n", point+1, 0);
                        fputs(str, result);
		}
	}//end for
	fclose(result); // the result file will not be used anymore

    // disposes objectsArray
	for (int i=0;i<numberOfObjects;i++) {
		delete [] objectsArray[i];
	}//end for
    delete [] objectsArray;

	// disposes the used structures
	delete sCluster;

	printElapsed(); // prints the elapsed time
	//getc(stdin);
	return 0; // success
}//end main
