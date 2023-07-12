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
* This file defines the classes stCountingTree and stPairCell.
*
* @version 3.0
* @author Caetano Traina Jr (caetano@icmc.usp.br)
* @author Agma Juci Machado Traina (agma@icmc.usp.br)
* @author Christos Faloutsos (christos@cs.cmu.edu)
* @author Elaine Parros Machado de Sousa (parros@icmc.usp.br)
* @author Ana Carolina Riekstin (anacarol@grad.icmc.usp.br)
* @author Robson Leonardo Ferreira Cordeiro (robson@icmc.usp.br)
*
*/
// Copyright (c) 2002-2009 GBDI-ICMC-USP

#ifndef __STCOUNTINGTREE_H
#define __STCOUNTINGTREE_H

#include "stCellId.h"

//----------------------------------------------------------------------------
// class stPairCell
//----------------------------------------------------------------------------
/**
* This class represents the structure of a Counting-tree cell.
*
* @version 3.0
* @author Caetano Traina Jr (caetano@icmc.usp.br)
* @author Agma Juci Machado Traina (agma@icmc.usp.br)
* @author Christos Faloutsos (christos@cs.cmu.edu)
* @author Elaine Parros Machado de Sousa (parros@icmc.usp.br)
* @author Ana Carolina Riekstin (anacarol@grad.icmc.usp.br)
* @author Robson Leonardo Ferreira Cordeiro (robson@icmc.usp.br)
*
*/
//---------------------------------------------------------------------------
class stPairCell {

   public:

      /**
      * Creates a Counting-tree cell.
      *
      * @param dimension The number of dimensions in the dataset.
      */
      stPairCell(int dimension) {
         index = new stCellId(dimension);
         sumOfPoints = 1; // if the cell exists, 
						  // then it contains at least one point
		 
		 // null pointers
         nextLevel = 0; 
         nextCell = 0;
         
         usedCell = 0;
         P = new int[dimension];
         for (int i=0;i<dimension;i++) {
           P[i]=0;
         }//end for

      }//end stPairCell

      /**
      * Disposes this structure.
      */
      ~stPairCell() {         
		 delete index;
         delete [] P;

      }//end ~stPairCell

      /**
      * Flag to control the usage of the cell.
      */
      char usedCell;

      /**
      * Vector with the half space counts of the cell 
	  * regarding each dimension (lower half).
      */
      int *P;

      /**
      * The spatial position of the cell inside its parent cell.
      */
      stCellId * index; // loc in the paper

      /**
      * Sum of points hitting the cell.
      */
      int sumOfPoints; // n in the paper

      /**
      * Pointer to the next Counting-tree level.
      */
      stPairCell * nextLevel; // ptr in the paper

      /**
      * Pointer to the next cell in the same tree node (linked list).
      */
      stPairCell * nextCell;

};//end stPairCell
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Class stCountingTree
//----------------------------------------------------------------------------
/**
* This class implements a structure to gather the points of a dataset and
* manipulates them.
*
* @version 3.0
* @author Caetano Traina Jr (caetano@icmc.usp.br)
* @author Agma Juci Machado Traina (agma@icmc.usp.br)
* @author Christos Faloutsos (christos@cs.cmu.edu)
* @author Elaine Parros Machado de Sousa (parros@icmc.usp.br)
* @author Ana Carolina Riekstin (anacarol@grad.icmc.usp.br)
* @author Robson Leonardo Ferreira Cordeiro (robson@icmc.usp.br)
* 
*/
//----------------------------------------------------------------------------
class stCountingTree{

   public:

      /**
      * Constructor method to initialize the memory structure.
      *
      * @param dimensions The number of dimensions in the dataset.
      * @param H The number of levels where the points are gathered in.
      */
      stCountingTree (int dimensions, int H);

      /**
      * Disposes this structure.
      */
      ~stCountingTree();

      /**
      * Gets the pointer to the first cell in the tree level one.
	  *
	  * @return firstLevel.
      */
      stPairCell * getFirstLevel() {
         return firstLevel;
      }//end getFirstLevel

      /**
      * Gets The normalizeSlope array.
	  *
	  * @return normalizeSlope.
      */
      double * getNormalizeSlope(){
         return normalizeSlope;
      }//end getNormalizeSlope

      /**
      * Gets The normalizeYInc array.
	  *
	  * @return normalizeYInc.
      */
      double * getNormalizeYInc(){
         return normalizeYInc;
      }//end getNormalizeYInc

      /**
      * Gets the half space counts of the root.
	  *
	  * @return P.
      */
      int * getP(){
         return P;
      }//end getP

      /**
      * Gets the count of points in the root.
	  *
	  * @return sumOfPoints.
      */
      int getSumOfPoints(){
         return sumOfPoints;
      }//end getSumOfPoints

      /**
      * Sets the array of normalization constants to be applied to each point,
      * when inserting it in the memory structure.
      *
      * @param slope Array with the values of X points.
      * @param yInc  Array with the values of Y points.
      */
      void normalize(double * slope, double * yInc);

      /**
      * Inserts one point into the memory structure.
      *
      * @param point The point.
      * @param normPoint Array to receive the point after 
	  *                  normalization.
      * @return A value that determines if the point was 
	  *         correctly inserted (1) or not (0).
      */
      char point(double * point, double * normPoint);

   private:

      /**
      * Number of dimensions in the dataset.
      */
      int numberOfDimensions;

      /**
      * Number of different hyper-square sizes (resolution levels).
      */
      int H;

      /**
      * Pointer to the node in level one of the tree.
      */
      stPairCell * firstLevel;

      /**
      * Vector with the half space counts of the root
	  * regarding each dimension (lower half).
      */
      int *P;

      /**
      * Count of points in the root.
      */
      int sumOfPoints;

      /**
      * Vectors used to store the normalization constants to be
      * applied to each point when inserting it in the memory structure.
      */
      double * normalizeSlope;
      double * normalizeYInc;

      /**
      * Recursively deletes the tree from the root.
	  *
	  * @param cell Pointer to the current subtree to be deleted.
      */
      void reset(stPairCell * cell);

      /**
      * Internal recursive procedure used to insert one point into one
      * node in one tree level.
	  *
	  * @param currentLevel The tree level to insert the point.
	  * @param pPC If currentLevel == 1, a pointer to the pointer firstLevel;
				   Otherwise, a pointer to the pointer parentCell->nextLevel.
	  * @param min The minimum bound of the current node's parent
	  *            in each dimension.
	  * @param max The maximum bound of the current node's parent
	  *            in each dimension.
	  * @param normPoint The point (previously normalized) to be inserted.
	  * @param parentCell Pointer to the current node's parent.
      */
      void cellProcess(int currentLevel, stPairCell ** pPC,
                       double * min, double * max, 
					   double * normPoint, stPairCell * parentCell);

};//end stCountingTree

#endif //__STSTORAGEGRID_H
