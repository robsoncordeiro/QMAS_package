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
* This file implements the class stCorrelationClustering.
*
* @version 1.0
* @author Robson Leonardo Ferreira Cordeiro (robson@icmc.usp.br)
* @author Agma Juci Machado Traina (agma@icmc.usp.br)
* @author Christos Faloutsos (christos@cs.cmu.edu)
* @author Caetano Traina Jr (caetano@icmc.usp.br)
*/
// Copyright (c) 2002-2009 GBDI-ICMC-USP

//------------------------------------------------------------------------------
// class stCorrelationClustering
//------------------------------------------------------------------------------

#include "stCorrelationClustering.h"

// default values
#define initialLevel 2

stCorrelationClustering::stCorrelationClustering (
                            double ** objectsArray, int numberOfDimensions,
                            int normalizeFactor, int numberOfObjects,
                            int centralConvolutionValue, int neighbourhoodConvolutionValue,
                            double pThreshold, int H, int hardClustering) {

    // stores the number of dimensions (database), the number of grid levels (counting tree) 
	// and the kind of clustering (soft or hard)
    this->numberOfDimensions = numberOfDimensions;
    this->H = H;
	this->hardClustering = hardClustering;

    // builds the counting tree and inserts objects on it
    calcTree = new stCountingTree(numberOfDimensions, H);
    fastDistExponent(numberOfObjects, objectsArray, normalizeFactor);

    // builds vectors to describe the positions of a cluster center in the data space
    minBetaClusterCenter = new double[numberOfDimensions];
    maxBetaClusterCenter = new double[numberOfDimensions];

    // builds vectors to the parents of a cluster center and of a neighbour
    betaClusterCenterParents = new stPairCell*[H];
    neighbourParents = new stPairCell*[H];

    // builds auxiliar vectors used to search for the relevant attributes
    attributesRelevance = new double[numberOfDimensions];

    // builds auxiliar vector to describe which neighbours belong to a found cluster
    neighbourhood = new char[numberOfDimensions];

    // stores the convolution matrix (center and direct neighbours)
    this->centralConvolutionValue=centralConvolutionValue;
    this->neighbourhoodConvolutionValue=neighbourhoodConvolutionValue;

    // stores the pThreshold
    this->pThreshold = pThreshold;

    // initiates the number of found clusters
    numBetaClusters=numCorrelationClusters=0;

    // defines the maximum number of beta-clusters
    maxNumBetaClusters=2000;

    // builds pointers to describe the found clusters
    minBetaClusters = new double*[maxNumBetaClusters];
    maxBetaClusters = new double*[maxNumBetaClusters];
    dimBetaClusters = new char*[maxNumBetaClusters];
	dimCorrelationClusters = new char*[maxNumBetaClusters];
	correlationClustersBelongings = new int[maxNumBetaClusters];
    for (int i=0; i<maxNumBetaClusters; i++) { // initiation
      minBetaClusters[i]=maxBetaClusters[i]=0;
      dimBetaClusters[i]=dimCorrelationClusters[i]=0;
	  correlationClustersBelongings[i]=-1;
    }//end for

}//end stCorrelationClustering::stCorrelationClustering

//---------------------------------------------------------------------------
stCorrelationClustering::~stCorrelationClustering() {

    // disposes the used structures
    delete [] neighbourhood;    
    delete [] attributesRelevance;
    delete [] minBetaClusterCenter;
    delete [] maxBetaClusterCenter;
	delete [] correlationClustersBelongings;
    delete [] betaClusterCenterParents;
    delete [] neighbourParents;
    for (int i=0; i<numBetaClusters; i++) {
      delete [] minBetaClusters[i];
      delete [] maxBetaClusters[i];
      delete [] dimBetaClusters[i];
    }//end for
	for (int i=0; i<numCorrelationClusters; i++) {
	  delete [] dimCorrelationClusters[i];	
	}//end for
	delete [] dimCorrelationClusters;
    delete [] minBetaClusters;
    delete [] maxBetaClusters;
    delete [] dimBetaClusters;
    delete calcTree;

}//end stCorrelationClustering::~stCorrelationClustering

//---------------------------------------------------------------------------
void stCorrelationClustering::findCorrelationClusters() {

    // defines when a new cluster is found
    int ok, center, total;
    do { // looks for a cluster in each loop
      ok=0; // no new cluster was found
      // defines the initial grid level to analyze
      int level=initialLevel;
      do { // analyzes each level until a new cluster is found
        // apply the convolution matrix to each grid cell in the current level
        walkThroughConvolution(level); // finds the cell with the bigest convolution value		
        
		if (betaClusterCenter) {          
          betaClusterCenter->usedCell = 1; // visited cell		  
		  if (level) {
			// pointer to a neighbour of the father
			stPairCell *fatherNeighbour;			
			for (int i=0; i<numberOfDimensions; i++) {			   
	  			// initiates total with the number of points in the father
				total=betaClusterCenterParents[level-1]->sumOfPoints;
				// discovers the number of points in the center
				if (betaClusterCenter->index->getBitValue(i,numberOfDimensions)) {
  				  center = total - betaClusterCenterParents[level-1]->P[i];
				} else {
				  center = betaClusterCenterParents[level-1]->P[i];
				}//end if
				// looks for the points in the direct neighbours of the father
				fatherNeighbour = internalNeighbour(i,betaClusterCenterParents[level-1],(level-1)?betaClusterCenterParents[level-2]:0);
				if (fatherNeighbour) {
				  total += fatherNeighbour->sumOfPoints;
				}//end if
				fatherNeighbour = externalNeighbour(i,betaClusterCenterParents[level-1],betaClusterCenterParents,0,level-1);
				if (fatherNeighbour) {
				  total += fatherNeighbour->sumOfPoints;
				}//end if
  			    // percentual of points in the center related to the average
                attributesRelevance[i] = (100*center)/((double)total/6);				
				// right critical value for the statistical test
				int criticalValue = GetCriticalValueBinomialRight2(total, (double)1/6, pThreshold);
				if (center > criticalValue) {
				  ok=1; // new cluster found
				}//end if
			}//end for
		  } else { // analyzes each dimension based on the points distribution of the entire database
				// initiate the total of points
				total=calcTree->getSumOfPoints();			
				for (int i=0; i<numberOfDimensions; i++) {			   
					// discovers the number of points in the center
					if (betaClusterCenter->index->getBitValue(i,numberOfDimensions)) {
						center = total - (calcTree->getP())[i];
					} else {
						center = (calcTree->getP())[i];
    				}//end if			
  			        // percentual of points in the center related to the average
                    attributesRelevance[i] = (100*center)/((double)total/2);					
					// right critical value for the statistical test
					int criticalValue = GetCriticalValueBinomialRight2(total, (double)1/2, pThreshold);
					if (center > criticalValue) {
					  ok=1; // new cluster found
					}//end if
				}//end for
		  }//end if		
		}//end if

        if (!ok) {
          level++; // next level to be analyzed
        }//end if
      } while (!ok && level < H);//end do while

      if (ok) { // if a new cluster was found...	    

        // discovers the cThreshold based on the minimum description length method
        double cThreshold = calcCThreshold(attributesRelevance);

        // new cluster found
        numBetaClusters++;
		printf("a beta-cluster was found at the Counting-tree level %d.\n",level); // prints the level in which a new beta-cluster was found
        minBetaClusters[numBetaClusters-1] = new double[numberOfDimensions];
        maxBetaClusters[numBetaClusters-1] = new double[numberOfDimensions];
        dimBetaClusters[numBetaClusters-1] = new char[numberOfDimensions];

        // important dimensions
        for (int i=0;i<numberOfDimensions;i++) {
          dimBetaClusters[numBetaClusters-1][i]=(attributesRelevance[i] >= cThreshold);
        }//end for

        // analyzes neighbours in important dimensions to verify which of them also belong to the found cluster
        for (int i=0;i<numberOfDimensions;i++) {
          neighbourhood[i]='N'; // no direct neighbour belongs to the cluster
        }//end for
        // center's position in the data space
        cellPosition(betaClusterCenter,betaClusterCenterParents,minBetaClusterCenter,maxBetaClusterCenter,level);
        // for each important dimension, analyzes internal and external neighbours to decide if they also
        // belong to the cluster
        for (int i=0;i<numberOfDimensions;i++) {
          if (dimBetaClusters[numBetaClusters-1][i]) {
            // looks for the internal neighbour
            neighbour = internalNeighbour(i,betaClusterCenter,(level)?betaClusterCenterParents[level-1]:0);
            if (neighbour) { // internal neighbour in important dimension always belongs to the cluster
              // neighbour's position in the data space
              cellPositionDimensionE_j(neighbour,betaClusterCenterParents,&minNeighbour,&maxNeighbour,level,i);
              if (maxBetaClusterCenter[i] > maxNeighbour) {
                neighbourhood[i]='I'; // inferior neighbour in i belongs to the cluster
              } else {
                neighbourhood[i]='S'; // superior neighbour in i belongs to the cluster
              }//end if
            }//end if

            // looks for the external neighbour
            for (int j=0;j<level;j++) {
              neighbourParents[j]=betaClusterCenterParents[j];
            }//end for
            neighbour = externalNeighbour(i,betaClusterCenter,betaClusterCenterParents,neighbourParents,level);
            if (neighbour) { // analyzes external neighbour to decide if it belongs to the cluster
              if (neighbourhood[i] == 'N') {
                // neighbour's position in the data space
                cellPositionDimensionE_j(neighbour,neighbourParents,&minNeighbour,&maxNeighbour,level,i);
                if (maxBetaClusterCenter[i] > maxNeighbour) {
                  neighbourhood[i]='I'; // inferior neighbour in i belongs to the cluster
                } else {
                  neighbourhood[i]='S'; // superior neighbour in i belongs to the cluster
                }//end if
              } else {
                neighbourhood[i]='B'; // both inferior and superior neighbours in i belong to the cluster
              }//end if
            }//end if
          }//end if
        }//end for        

        // stores the description of the found cluster
        double length = maxBetaClusterCenter[0]-minBetaClusterCenter[0];
        for (int i=0;i<numberOfDimensions;i++) {
          if (dimBetaClusters[numBetaClusters-1][i]) { // dimension important to the cluster
            // analyzes if the neighbours in i also belong to the cluster
            switch (neighbourhood[i]) {
              case 'B': // both inferior and superior neighbours in i belong to the cluster
                  minBetaClusterCenter[i]-=length;
                  maxBetaClusterCenter[i]+=length;
                  break;
              case 'S': // superior neighbour in i belongs to the cluster
                  maxBetaClusterCenter[i]+=length;
                  break;
              case 'I': // inferior neighbour in i belongs to the cluster
                  minBetaClusterCenter[i]-=length;
            }//end switch
            // new cluster description - relevant dimension
            minBetaClusters[numBetaClusters-1][i] = minBetaClusterCenter[i];
            maxBetaClusters[numBetaClusters-1][i] = maxBetaClusterCenter[i];
          } else {
            // new cluster description - irrelevant dimension
            minBetaClusters[numBetaClusters-1][i] = 0;
            maxBetaClusters[numBetaClusters-1][i] = 1;
          }//end if
        }//end for
      }//end if
      // stops when no new cluster is found
    } while (ok);//end do while

	printf("\n%d beta-clusters were found.\n",numBetaClusters); // prints the number of beta-clusters found
	mergeBetaClusters(); // merges clusters that share some database space
	printf("\n%d correlation clusters left after the merging fase.\n",numCorrelationClusters); // prints the number of correlation clusters found

}//end stCorrelationClustering::findCorrelationClusters

//---------------------------------------------------------------------------
double stCorrelationClustering::calcCThreshold(double *attributesRelevance) {

  double *sortedRelevance = new double[numberOfDimensions];
  for (int i=0;i<numberOfDimensions;i++) {
    sortedRelevance[i]=attributesRelevance[i];
  }//end for  
  qsort(sortedRelevance,numberOfDimensions,sizeof(double),compare); // sorts the relevances vector
  double cThreshold = sortedRelevance[minimumDescriptionLength(sortedRelevance)];
  delete [] sortedRelevance;
  return cThreshold;

}//end stCorrelationClustering::calcCThreshold

//---------------------------------------------------------------------------
int stCorrelationClustering::minimumDescriptionLength(double *sortedRelevance) {

  int cutPoint=-1;
  double preAverage, postAverage, descriptionLength, minimumDescriptionLength;
  for (int i=0;i<numberOfDimensions;i++) {
    descriptionLength=0;
    // calculates the average of both sets
    preAverage=0;
    for (int j=0;j<i;j++) {
      preAverage+=sortedRelevance[j];
    }//end for
    if (i) {
      preAverage/=i;      
	  descriptionLength += (ceil(preAverage)) ? (log10(ceil(preAverage))/log10((double)2)) : 0;	// changes the log base from 10 to 2
    }//end if
    postAverage=0;
    for (int j=i;j<numberOfDimensions;j++) {
      postAverage+=sortedRelevance[j];
    }//end for
    if (numberOfDimensions-i) {
      postAverage/=(numberOfDimensions-i);      
	  descriptionLength += (ceil(postAverage)) ? (log10(ceil(postAverage))/log10((double)2)) : 0;	// changes the log base from 10 to 2
    }//end if
    // calculates the description length
    for (int j=0;j<i;j++) {      
	  descriptionLength += (ceil(fabs(preAverage-sortedRelevance[j]))) ? (log10(ceil(fabs(preAverage-sortedRelevance[j])))/log10((double)2)) : 0;	// changes the log base from 10 to 2
    }//end for
    for (int j=i;j<numberOfDimensions;j++) {      
	  descriptionLength += (ceil(fabs(postAverage-sortedRelevance[j]))) ? (log10(ceil(fabs(postAverage-sortedRelevance[j])))/log10((double)2)) : 0;	// changes the log base from 10 to 2
    }//end for
    // verify if this is the best cut point
    if (cutPoint==-1 || descriptionLength < minimumDescriptionLength) {
      cutPoint=i;
      minimumDescriptionLength = descriptionLength;
    }//end if
  }//end for
  return cutPoint;

}//end stCorrelationClustering::minimumDescriptionLength

//---------------------------------------------------------------------------
void stCorrelationClustering::walkThroughConvolution(int level) {

  // pointers to the parents of a cell
  stPairCell **parentsVector = new stPairCell*[level];
  for (int i=0;i<level;i++) {
    parentsVector[i]=0;
  }//end for
  int biggestConvolutionValue = -MAXINT;
  betaClusterCenter = 0;
  // iniciate the recursive process to analyze a level
  walkThroughConvolutionRecursive(calcTree->getFirstLevel(),parentsVector,level,0,&biggestConvolutionValue);
  // disposes parentsVector
  delete [] parentsVector;

}//end stCorrelationClustering::walkThroughConvolution

//---------------------------------------------------------------------------
void stCorrelationClustering::walkThroughConvolutionRecursive(stPairCell *cell,
                               stPairCell **cellParents, int level, int currentLevel,
                               int *biggestConvolutionValue) {

  if (cell) {
    if (currentLevel < level) { // goes to the next level storing the pointer to the father
      while (cell) {
        cellParents[currentLevel] = cell;
        walkThroughConvolutionRecursive(cell->nextLevel,cellParents,level,(currentLevel+1),biggestConvolutionValue);
        cell=cell->nextCell; // next cell in the same node
      }//end while
    } else { // this is the level to be analyzed
      int newConvolutionValue;
      char clusterFoundBefore;
      double *maxCell = new double[numberOfDimensions];
      double *minCell = new double[numberOfDimensions];
      while (cell) {
        // Does not analyze cells analyzed before and cells that can't be the biggest convolution center.
        // It speeds up the algorithm, specially when neighbourhoodConvolutionValue <= 0
        if ( (!cell->usedCell) && ((neighbourhoodConvolutionValue > 0) || ((cell->sumOfPoints*centralConvolutionValue) > *biggestConvolutionValue)) ) {
          // discovers the position of cell in the data space
          cellPosition(cell,cellParents,minCell,maxCell,level);
          // verifies if this cell belongs to a cluster found before
          clusterFoundBefore=0;
          for (int i=0;!clusterFoundBefore && i<numBetaClusters;i++) {
            clusterFoundBefore = 1;
            for (int j=0; clusterFoundBefore && j<numberOfDimensions; j++) {
              // Does not cut off cells in a level upper than the level where a cluster was found
              if (!(maxCell[j] <= maxBetaClusters[i][j] && minCell[j] >= minBetaClusters[i][j])) {
                clusterFoundBefore = 0;
              }//end if
            }//end for
          }//end for

          if (!clusterFoundBefore) { // the cell doesn't belong to any found cluster
            // applies the convolution matrix to cell
            if (neighbourhoodConvolutionValue) {
              newConvolutionValue=applyConvolution(cell,cellParents,currentLevel);
            } else {
              newConvolutionValue=centralConvolutionValue*cell->sumOfPoints; // when the neighbourhood weight is zero
            }//end if
            if (newConvolutionValue > *biggestConvolutionValue) {
              // until now, cell is the biggest convolution value
              *biggestConvolutionValue = newConvolutionValue;
              betaClusterCenter = cell;
              for (int j=0 ; j<currentLevel ; j++) {
                betaClusterCenterParents[j] = cellParents[j];
              }//end for
            }//end if
          }//end if
        }//end if
        cell=cell->nextCell; // next cell in the same node
      }//end while
      // disposes minCell and maxCell
      delete [] minCell;
      delete [] maxCell;
    }//end if
  }//end if

}//end stCorrelationClustering::walkThroughConvolutionRecursive

//---------------------------------------------------------------------------
int stCorrelationClustering::applyConvolution(stPairCell *cell,
                                stPairCell **cellParents, int level) {

  stPairCell *neighbour;
  int newValue = cell->sumOfPoints*centralConvolutionValue;
  // looks for the neighbours
  for (int k=0;k<numberOfDimensions;k++) {
    neighbour = internalNeighbour(k,cell,((level)?cellParents[level-1]:0));
    if (neighbour) {
      newValue+=(neighbour->sumOfPoints*neighbourhoodConvolutionValue);
    }//end if
    neighbour = externalNeighbour(k,cell,cellParents,0,level);
    if (neighbour) {
      newValue+=(neighbour->sumOfPoints*neighbourhoodConvolutionValue);
    }//end if
  }//end for
  // return the cell value after applying the convolution matrix
  return newValue;

}//end stCorrelationClustering::applyConvolution

//---------------------------------------------------------------------------
char stCorrelationClustering::cellPosition(stPairCell *cell,
                                stPairCell **cellParents, double *min, double *max, int level) {

  if (cell) {
    if (level) {
      cellPosition(cellParents[level-1],cellParents,min,max,level-1);
      for (int i=0; i<numberOfDimensions; i++) {
        if (cell->index->getBitValue(i,numberOfDimensions)) { // bit in the position i is 1
          min[i] += ((max[i]-min[i])/2);
        } else { // bit in the position i is 0
          max[i] -= ((max[i]-min[i])/2);
        }//end if
      }//end for
    } else { // level zero
      for (int i=0; i<numberOfDimensions; i++) {
        if (cell->index->getBitValue(i,numberOfDimensions)) { // bit in the position i is 1
          min[i] = 0.5;
          max[i] = 1;
        } else { // bit in the position i is 0
          min[i] = 0;
          max[i] = 0.5;
        }//end if
      }//end for
    }//end if
    return 1;
  }//end if
  return 0;

}//end stCorrelationClustering::cellPosition

//---------------------------------------------------------------------------
char stCorrelationClustering::cellPositionDimensionE_j(stPairCell *cell,
                                stPairCell **cellParents, double *min, double *max, int level, int j) {

  if (cell) {
    if (level) {
      cellPositionDimensionE_j(cellParents[level-1],cellParents,min,max,level-1,j);      
      if (cell->index->getBitValue(j,numberOfDimensions)) { // bit in the position j is 1
        *min += ((*max-*min)/2);
      } else { // bit in the position j is 0
        *max -= ((*max-*min)/2);
      }//end if      
    } else { // level zero
      if (cell->index->getBitValue(j,numberOfDimensions)) { // bit in the position j is 1
        *min = 0.5;
        *max = 1;
      } else { // bit in the position j is 0
        *min = 0;
        *max = 0.5;
      }//end if
    }//end if
    return 1;
  }//end if
  return 0;

}//end stCorelationClustering::cellPositionDimensionE_j

//---------------------------------------------------------------------------
stPairCell *stCorrelationClustering::externalNeighbour(int dimIndex, stPairCell *cell,
                                                        stPairCell **cellParents,
                                                        stPairCell **neighbourParents, int level) {

   if (level) {
     stPairCell *father = cellParents[level-1], *fathersNeighbour;
     if (cell->index->getBitValue(dimIndex,numberOfDimensions) ^ father->index->getBitValue(dimIndex,numberOfDimensions)) { // XOR - different bit values -> starts going down in the tree
       // creates the index that the father's neighbour should have
       stCellId *fathersNeighbourId = new stCellId(numberOfDimensions);
       *fathersNeighbourId=*father->index;
       fathersNeighbourId->invertBit(dimIndex,numberOfDimensions);
       // looks for the father's neighbour
       if ((*fathersNeighbourId == *father->index) < 0) {
         father = (level-1) ? cellParents[level-2]->nextLevel : calcTree->getFirstLevel(); // through the previous level, goes to the beginning of the list
       }//end if
       while (father && (*fathersNeighbourId == *father->index) > 0) {
         father = father->nextCell;
       }//end while
       fathersNeighbour = (father && *fathersNeighbourId == *father->index) ? 0 : father;
       delete fathersNeighbourId;
     } else {  // equal bit values -> continues going up in the tree
       fathersNeighbour = externalNeighbour(dimIndex,father,cellParents,neighbourParents,level-1); // recursively, finds the external neighbour of the father
     }//end if
     if (fathersNeighbour) { // father's neighbour was found
       // creates the index that the son's neighbour should have
       stCellId *sonsNeighbourId = new stCellId(numberOfDimensions);
       *sonsNeighbourId=*cell->index;
       sonsNeighbourId->invertBit(dimIndex,numberOfDimensions);
       // looks for the son's neighbour
       cell = fathersNeighbour->nextLevel;
       while (cell && (*sonsNeighbourId == *cell->index) > 0) {
         cell=cell->nextCell;
       }//end while
       if (cell && *sonsNeighbourId == *cell->index) {
         cell = 0;
       }//end if
       if (neighbourParents) {
         neighbourParents[level-1] = fathersNeighbour; // if it existis, stores the pointer to the father's neighbour
       }//end if
       delete sonsNeighbourId;
       return cell; // return the external son's neighbour in the dimension i
     }//end if
     return 0; // there is no father's neighbour
   }//end if
   return 0; // a cell in level zero never has an external neighbour

}//end stCorrelationClustering::externalNeighbour

//---------------------------------------------------------------------------
stPairCell *stCorrelationClustering::internalNeighbour(int dimIndex,
                                                        stPairCell *cell,
                                                        stPairCell *father) {

  // creates the id that the neighbour should have
  stCellId *neighboursId = new stCellId(numberOfDimensions);
  *neighboursId=*cell->index;
  neighboursId->invertBit(dimIndex,numberOfDimensions);
  // looks for the neighbour
  if ((*cell->index == *neighboursId) > 0) {
    cell = (father) ? father->nextLevel : calcTree->getFirstLevel(); // through the previous level, goes to the beginning of the list
  }//end if
  while (cell && (*cell->index == *neighboursId) < 0) {
    cell = cell->nextCell;
  }//end while
  if (cell && *cell->index == *neighboursId) {
    cell = 0;
  }//end if
  delete neighboursId;
  return cell;

}//end stCorrelationClustering::internalNeighbour

//---------------------------------------------------------------------------
void stCorrelationClustering::fastDistExponent(int numberOfObjects,
                                double **objectsArray, int normalizeFactor) {

   double *minD, *maxD, biggest;
   double *onePoint, *resultPoint, *a, *b; // y=Ax+B to normalize each dataset.
   double normalizationFactor = 1.0;

   minD = (double *) calloc ((1+numberOfDimensions),sizeof(double));
   maxD = (double *) calloc ((1+numberOfDimensions),sizeof(double));
   a = (double *) calloc(numberOfDimensions,sizeof(double));
   b = (double *) calloc(numberOfDimensions,sizeof(double));
   onePoint = (double *) calloc(numberOfDimensions,sizeof(double));
   resultPoint = (double *) calloc(numberOfDimensions,sizeof(double));

   // normalizes the data
   minMax(numberOfObjects, objectsArray, minD, maxD);
   biggest = 0; // for Normalize==0, 1 or 3
   // normalize=0->Independent, =1->mantain proportion, =2->Clip
   //          =3->Geo Referenced
   if (normalizeFactor == 2) {
     biggest = MAXDOUBLE;
   }//end if

   for (int i=0; i<numberOfDimensions; i++) {
     a[i] = (maxD[i] - minD[i]) * normalizationFactor; //a takes the range of each dimension
     b[i] = minD[i];
     if (a[i] == 0) {
       a[i] = 1;
     }//end if
   }//end for

   for (int i=0; i<numberOfDimensions; i++) {
     if ((normalizeFactor < 2 || normalizeFactor == 3) && biggest < a[i]) {
       biggest = a[i];
     }//end if
     if (normalizeFactor == 2 && biggest > a[i]) {
       biggest = a[i];
     }//end if
   }//end for

   if (normalizeFactor != 0) {
     for (int i=0; i<numberOfDimensions; i++) {
       a[i] = biggest; // normalized keeping proportion
     }//end for
     /* when we have the proportional normalization, every A[i] are gonna receive
     the biggest range.*/
   }//end if

   if (normalizeFactor >= 0) {
     calcTree->normalize(a,b); // if there is some normalization
   }//end if

   //process each point of objectArray
   for (int i=0; i<numberOfObjects; i++) {
     for (int j=0; j<numberOfDimensions; j++) {
       onePoint[j] = objectsArray[i][j];
     }//end for
     calcTree->point(onePoint,resultPoint); //add to the grid structure
   }//end for

   // disposes used memory
   delete[] onePoint;
   delete[] resultPoint;
   delete[] a;
   delete[] b;
   delete[] minD;
   delete[] maxD;

}//end stCorrelationClustering::FastDistExponent

//---------------------------------------------------------------------------
void stCorrelationClustering::minMax(int numberOfObjects, double**objectsArray,
                                double * min, double * max) {

  for (int j=0; j<numberOfDimensions; j++){ // sets the values to the minimum/maximum possible here
    min[j] = MAXDOUBLE;
    max[j] = -MAXDOUBLE;
  }// end for
  // looking for the minimum and maximum values
  for (int i=0; i<numberOfObjects; i++) {
    for (int j=0; j<numberOfDimensions; j++) {
      if (objectsArray[i][j] < min[j]) {
        min[j] = objectsArray[i][j];
      }//end if
      if (objectsArray[i][j] > max[j]) {
        max[j] = objectsArray[i][j];
      }//end if
    }//end for
  }//end for

}//end stCorrelationClustering::MinMax

//---------------------------------------------------------------------------
void stCorrelationClustering::mergeBetaClusters() {

	int shareSpace, i=0, j, k, aux;
	// merges beta-clusters
	while (i<numBetaClusters) {
		j=i+1;		
		while (j<numBetaClusters) {
			// discovers if beta-cluster i shares database space with beta-cluster j
			shareSpace=1;
			for(k=0; shareSpace && k<numberOfDimensions; k++) {
				if (!(maxBetaClusters[i][k] > minBetaClusters[j][k] && minBetaClusters[i][k] < maxBetaClusters[j][k])) {
					shareSpace=0; // beta-clusters i and j do not share database space
				}//end if
			}//end for
			if (hardClustering && shareSpace) { // merges both beta-clusters
				if (correlationClustersBelongings[i]==-1 && correlationClustersBelongings[j]==-1) { // both clusters belong to no merged cluster
					correlationClustersBelongings[i]=correlationClustersBelongings[j]=numCorrelationClusters++; // new merged cluster
				} else {
					if (correlationClustersBelongings[i]!=-1 && correlationClustersBelongings[j]!=-1) { // both clusters belong to some merged cluster(s)
						if (correlationClustersBelongings[i]!=correlationClustersBelongings[j]) { // both clusters belong to different merged clusters
							numCorrelationClusters--;
							for (k=0; k<numBetaClusters; k++) {
								if (k!=i && k!=j && (correlationClustersBelongings[k]==correlationClustersBelongings[i] || correlationClustersBelongings[k]==correlationClustersBelongings[j])) {
									correlationClustersBelongings[k]=(correlationClustersBelongings[i]>correlationClustersBelongings[j]) ? correlationClustersBelongings[j] : correlationClustersBelongings[i];
								}//end if
							}//end for
							if (correlationClustersBelongings[i]>correlationClustersBelongings[j]) {
								aux = correlationClustersBelongings[i]; // deleted cluster 
								correlationClustersBelongings[i]=correlationClustersBelongings[j];								 
							} else { 
								aux = correlationClustersBelongings[j]; // deleted cluster 
								correlationClustersBelongings[j]=correlationClustersBelongings[i];
							}//end if
							for (k=0; k<numBetaClusters; k++) {
								if (correlationClustersBelongings[k] > aux) {
									correlationClustersBelongings[k]--;
								}//end if
							}//end for
						}//end if
					} else { // only one of the beta-clusters belongs to some merged cluster
						(correlationClustersBelongings[i]==-1) ? correlationClustersBelongings[i]=correlationClustersBelongings[j] : correlationClustersBelongings[j]=correlationClustersBelongings[i];
					}//end if
				}//end if
			}//end if
			j++; // next beta-cluster
		}//end while
		if (correlationClustersBelongings[i] == -1) {
			correlationClustersBelongings[i] = numCorrelationClusters++; // new merged cluster
		}//end if
		i++; // next cluster
	}//end while

	// important dimensions to the merged clusters
	for (i=0; i<numCorrelationClusters; i++) {
		dimCorrelationClusters[i] = new char[numberOfDimensions];
		for (j=0; j<numberOfDimensions;j++) {
			dimCorrelationClusters[i][j]=0;
		}//end for
	}//end for
	for (i=0; i<numBetaClusters; i++) {
		for (j=0; j<numberOfDimensions; j++) {
			if (!dimCorrelationClusters[correlationClustersBelongings[i]][j]) {
				dimCorrelationClusters[correlationClustersBelongings[i]][j] = dimBetaClusters[i][j];
			}//end if
		}//end for
	}//end for

}//end stCorrelationClustering::mergeBetaClusters
