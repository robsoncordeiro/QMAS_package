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
* This file implements the classes stCountingTree and stPairCell.
*
* @version 3.0
* @author Caetano Traina Jr (caetano@icmc.usp.br)
* @author Agma Juci Machado Traina (agma@icmc.usp.br)
* @author Christos Faloutsos (christos@cs.cmu.edu)
* @author Elaine Parros Machado de Sousa (parros@icmc.usp.br)
* @author Ana Carolina Riekstin (anacarol@grad.icmc.usp.br)
* @author Robson Leonardo Ferreira Cordeiro (robson@icmc.usp.br)
*/
// Copyright (c) 2002-2009 GBDI-ICMC-USP

#include "stCountingTree.h"

//------------------------------------------------------------------------------
// class stCountingTree
//------------------------------------------------------------------------------
stCountingTree::stCountingTree (int dimensions, int H) {
   //initializing
   numberOfDimensions = dimensions;
   this->H = H;
   firstLevel = 0;
   sumOfPoints = 0;

   normalizeSlope = new double[numberOfDimensions];
   normalizeYInc = new double[numberOfDimensions];
   P = new int[numberOfDimensions];
   for(int i=0; i<numberOfDimensions; i++) {
      normalizeSlope[i] = 1;
      normalizeYInc[i] = 0;
	  P[i] = 0;
   }//end for

}//end stCountingTree::stCountingTree

//---------------------------------------------------------------------------
stCountingTree::~stCountingTree() {
   // recursively deletes the tree from the root
   reset(firstLevel); 
   // deletes the other arrays
   delete [] normalizeSlope;
   delete [] normalizeYInc;
   delete [] P;
   
}//end stCountingTree::~stCountingTree

//---------------------------------------------------------------------------
void stCountingTree::reset(stPairCell * cell) {
   stPairCell * aux;
   while (cell) { // while there are cells in this node
     reset(cell->nextLevel); // resets the subtree
     aux = cell; // holds the cell
     cell = cell->nextCell; // gets next cell in the node
     delete aux; // deletes the holded cell
   }//end while
}//end stCountingTree::reset

//---------------------------------------------------------------------------
void stCountingTree::normalize(double * slope, double * yInc) {
   for (int i=0; i<numberOfDimensions; i++) {
      normalizeSlope[i] = slope[i];
      normalizeYInc[i] = yInc[i];
   }//end for
}//end stCountingTree::normalize

//---------------------------------------------------------------------------
char stCountingTree::point(double * point, double * normPoint) {
   double * min = new double[numberOfDimensions];
   double * max = new double[numberOfDimensions];
   char out = 0;

   // normalizes the point
   for (int i=0; i<numberOfDimensions; i++) {
      normPoint[i] = (point[i]-normalizeYInc[i])/normalizeSlope[i];
      if (normPoint[i] < 0 || normPoint[i] > 1) {
         out = 1; // invalid point
      }//end if
      min[i] = 0;
      max[i] = 1;
   }//end for
      
   if (!out) { // if the point is valid
	   sumOfPoints++; // counts this point in level 0
	   // recursively counts this point in all levels deeper than level 0
	   cellProcess(1, &firstLevel, min, max, normPoint, 0);
   }//end if
   
   delete [] min;
   delete [] max;
   
   return !out; // return the "validity" of the point
}//end stCountingTree::point

//---------------------------------------------------------------------------
void stCountingTree::cellProcess(int currentLevel, stPairCell ** pPC, double * min,
							     double * max, double * normPoint, stPairCell * parentCell) {

   int comparing; // the result when comparing a pair of ids
   stCellId * cellId; // id of the cell that covers normPoint in currentLevel
   double middle;
   stPairCell * pairCell = *pPC; // pointer in the previous level that links the 
								 //	current node to the tree
								 // if currentLevel == 1, pairCell == firstLevel
				                 // otherwise, pairCell == parentCell->nextLevel

   if (currentLevel < H) { // H determines when the recursive insertion
						   // of normPoint must stop
	   
	   // mounts the id of the cell that covers normPoint in the current level
	   // and stores in min/max this cell's lower and upper bounds
	   cellId = new stCellId(numberOfDimensions);
	   for (int i=0; i<numberOfDimensions; i++) {
		  middle = (max[i]-min[i])/2 + min[i];
		  if (normPoint[i] > middle) { // if normPoint[i] == middle, considers that the point
			                           // belongs to the cell in the lower half regarding i
    		 cellId->invertBit(i,numberOfDimensions); // turns the bit i to one
			 min[i] = middle;
		  } else {
			 // the bit i remains zero
			 max[i] = middle;
		  }//end if
	   }//end for

	   // looks for a cell in the current node whose index is equal to
	   // the computed one (cellId) and, if it does not exist,
	   // creates a cell with the computed id and inserts it in the
	   // sorted list (based on the id) that represents the current node
	   if (!pairCell) { // the new cell will be the first one in the 
		                // current node (empty list)
		  pairCell = new stPairCell(numberOfDimensions); // creates a new cell
		  *pPC = pairCell; // links the new cell to its parent
		  *(pairCell->index) = *cellId; // copies the computed id to the pairCell's id
	   } else { // this node already have cell(s) and thus, we follow the list...
		  char done = 0;
		  stPairCell * auxCell = 0;
		  while (pairCell && !done) { // searches the list looking for a cell
			                          // with an id equal to the computed one
			 comparing = (*(pairCell->index) == *cellId); // compares the ids
			 if (!comparing) { // pairCell is the cell we are looking for
				done = 1;
				pairCell->sumOfPoints++; // one more point covered by pairCell
			 } else { // pairCell is not the cell we are looking for, so,
				      // continue the walk in the list...                  
				if (comparing < 0) { // the id in pairCell is smaller than
									 // the computed one, so, we must keep looking
				   auxCell = pairCell; // holds the previous cell
				   pairCell = pairCell->nextCell; // next cell
				} else { // the id in pairCell is bigger than the computed one,
					     // so, a cell with the computed id does not exist
					     // and we must create and insert it in the list 
					     // before pairCell
				   if (!auxCell) { // inserts a new cell in the beginning of the list
					  auxCell = new stPairCell(numberOfDimensions); // creates a new cell
					  *pPC = auxCell; // links the new cell to its parent
					  auxCell->nextCell = pairCell; // links the existing list 
					                                // to the new cell
					  pairCell = auxCell;
				   } else { // inserts a new cell in the middle of the list
					  pairCell = new stPairCell(numberOfDimensions); // creates a new cell
					  // inserts pairCell in the list, after auxCell
					  pairCell->nextCell = auxCell->nextCell;
					  auxCell->nextCell = pairCell;
				   }//end if
				   *(pairCell->index) = *cellId; // copies the computed id to the new cell
				   done = 1;
				}//end if
			 }//end if
		  }// end while - searching list
		  if (!done) { // the whole list was travelled without finding in the desired cell
					   // thus, insert a new cell in the end of the list
			 pairCell = new stPairCell(numberOfDimensions); // creates a new cell
			 auxCell->nextCell = pairCell; // inserts pairCell in the end of the list
			 *(pairCell->index) = *cellId; // copies the computed id to the new cell
		  }//end if
	   }//end if

	   // updates the half space counts in level currentLevel-1
	   if (parentCell) { // updates in the parent cell (currentLevel-1 > 0)
		 for (int i=0; i<numberOfDimensions; i++) {
		   if (!cellId->getBitValue(i,numberOfDimensions)) {
			 parentCell->P[i]++; // counts normPoint, since it is in the 
			                     // parent's lower half regarding i
		   }//end if
		 }//end for
	   } else { // updates in the root (currentLevel-1 == 0)
		   for (int i=0; i<numberOfDimensions; i++) {
		     if (!cellId->getBitValue(i,numberOfDimensions)) {
  			   P[i]++; // counts normPoint, since it is in the 
			           // root's lower half regarding i
		     }//end if
		   }//end for		 
	   }//end if	  
	   
	   // deletes the computed id
	   delete cellId;

	   // counts normPoint in level currentLevel + 1
	   cellProcess(currentLevel+1,&(pairCell->nextLevel), min, max, normPoint, pairCell);
   }//end if

}//end stCountingTree::cellProcess
