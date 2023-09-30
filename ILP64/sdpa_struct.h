/* -------------------------------------------------------------

This file is a component of SDPA
Copyright (C) 2004 SDPA Project

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

------------------------------------------------------------- */

// printing presicion of vectors and matrices
#define P_FORMAT ((char*)"%+8.3e")
#define NO_P_FORMAT "NOPRINT"


#ifndef __sdpa_struct_h__
#define __sdpa_struct_h__

#include "sdpa_include.h"
#include "sdpa_block.h"

#define DATA_CAPSULE 1
  // DATA_CAPSULE 0 : Three Arrays (row,column,sp_ele)
  // DATA_CAPSULE 1 : Capsuled data storage
  
namespace sdpa {

class Vector
{
public:
  SDPA_INT nDim;
  double* ele;

  Vector();
  Vector(SDPA_INT nDim, double value = 0.0);
  ~Vector();

  void initialize(SDPA_INT nDim, double value = 0.0);
  void initialize(double value);
  void terminate();

  void setZero();
  void display(FILE* fpout = stdout, char* printFormat = P_FORMAT);
  void display(FILE* fpout,double scalar, char* printFormat = P_FORMAT);
  bool copyFrom(Vector& other);
};

class BlockVector
{
public:
  SDPA_INT  nBlock;
  SDPA_INT* blockStruct;

  Vector* ele;
  
  BlockVector();
  BlockVector(BlockStruct& bs, double value = 0.0);
  BlockVector(SDPA_INT nBlock, SDPA_INT* blockStruct, double value = 0.0);
  ~BlockVector();
  
  void initialize(BlockStruct& bs, double value = 0.0);
  void initialize(SDPA_INT nBlock, SDPA_INT* blockStruct, double value = 0.0);
  void initialize(double value);
  void terminate();

  void setZero();
  void display(FILE* fpout = stdout, char* printFormat = P_FORMAT);
  bool copyFrom(BlockVector& other);
};

class SparseMatrix
{
public:
  SDPA_INT nRow, nCol;

  enum Type { SPARSE, DENSE};
  Type type;
  
  SDPA_INT NonZeroNumber;
  // for memory
  SDPA_INT NonZeroCount;
  // currentry stored
  SDPA_INT NonZeroEffect;
  // use for calculation of F1,F2,F3 

  // for Dense
  double* de_ele;

  // for Sparse ; 0:sparse 1:dense
  enum dsType {DSarrays, DScapsule};
  dsType DataStruct;

  // for Sparse Data1 // dsArrays
  SDPA_INT*    row_index;
  SDPA_INT*    column_index;
  double* sp_ele;

  // for Sparse Data2 // dsCapsule
  typedef struct{
    SDPA_INT vRow;
    SDPA_INT vCol;
    double vEle;
  } SparseElement; // __attribute__( (aligned (16)));

  SparseElement* DataS;

  SparseMatrix();
  SparseMatrix(SDPA_INT nRow,SDPA_INT nCol, Type type, SDPA_INT NonZeroNumber);
  ~SparseMatrix();

  #if DATA_CAPSULE 
  void initialize(SDPA_INT nRow,SDPA_INT nCol, Type type, SDPA_INT NonZeroNumber,
		  dsType DataStruct = DScapsule);
  #else
  void initialize(SDPA_INT nRow,SDPA_INT nCol, Type type, SDPA_INT NonZeroNumber,
		  dsType DataStruct = DSarrays);
  #endif
  void terminate();

  void display(FILE* fpout = stdout, char* printFormat = P_FORMAT);
  bool copyFrom(SparseMatrix& other);

  void changeToDense(bool forceChange = false);
  void setZero();
  void setIdentity(double scalar = 1.0);

  bool sortSparseIndex(SDPA_INT&i, SDPA_INT& j);
};

class DenseMatrix
{
public:
  SDPA_INT nRow, nCol;

  enum Type { DENSE, COMPLETION};
  Type type;
  
  double* de_ele;

  DenseMatrix();
  DenseMatrix(SDPA_INT nRow,SDPA_INT nCol, Type type);
  ~DenseMatrix();

  void initialize(SDPA_INT nRow,SDPA_INT nCol, Type type);
  void terminate();
  
  void display(FILE* fpout = stdout, char* printFormat = P_FORMAT);
  bool copyFrom(DenseMatrix& other);
  bool copyFrom(SparseMatrix& other);

  void setZero();
  void setIdentity(double scalar = 1.0);
};

class SparseLinearSpace
{
public:
  SDPA_INT  SDP_sp_nBlock;
  SDPA_INT  SOCP_sp_nBlock;
  SDPA_INT  LP_sp_nBlock;

  SDPA_INT*  SDP_sp_index;
  SDPA_INT*  SOCP_sp_index;
  SDPA_INT*  LP_sp_index;

  SparseMatrix* SDP_sp_block;
  SparseMatrix* SOCP_sp_block;
  double* LP_sp_block;
  
  SparseLinearSpace();
  SparseLinearSpace(SDPA_INT SDP_nBlock, SDPA_INT* SDP_blockStruct, 
		    SDPA_INT* SDP_NonZeroNumber,
		    SDPA_INT SOCP_nBlock, SDPA_INT* SOCP_blockStruct,
		    SDPA_INT* SOCP_NonZeroNumber,
		    SDPA_INT LP_nBlock, bool* LP_NonZeroNumber);
  SparseLinearSpace(SDPA_INT SDP_sp_nBlock, 
                    SDPA_INT* SDP_sp_index,
                    SDPA_INT* SDP_sp_blockStruct, 
                    SDPA_INT* SDP_sp_NonZeroNumber,
                    SDPA_INT SOCP_sp_nBlock, 
                    SDPA_INT* SOCP_sp_index,
                    SDPA_INT* SOCP_sp_blockStruct,
                    SDPA_INT* SOCP_sp_NonZeroNumber,
                    SDPA_INT LP_sp_nBlock, 
                    SDPA_INT* LP_sp_index);
  ~SparseLinearSpace();

  // dense form of block index
  void initialize(SDPA_INT SDP_nBlock, SDPA_INT* SDP_blockStruct, 
		    SDPA_INT* SDP_NonZeroNumber,
		    SDPA_INT SOCP_nBlock, SDPA_INT* SOCP_blockStruct,
		    SDPA_INT* SOCP_NonZeroNumber,
		    SDPA_INT LP_nBlock, bool* LP_NonZeroNumber);
  // sparse form of block index      2008/02/27 kazuhide nakata
  void initialize(SDPA_INT SDP_sp_nBlock, 
                  SDPA_INT* SDP_sp_index,
                  SDPA_INT* SDP_sp_blockStruct, 
                  SDPA_INT* SDP_sp_NonZeroNumber,
                  SDPA_INT SOCP_sp_nBlock, 
                  SDPA_INT* SOCP_sp_index,
                  SDPA_INT* SOCP_sp_blockStruct,
                  SDPA_INT* SOCP_sp_NonZeroNumber,
                  SDPA_INT LP_sp_nBlock, 
                  SDPA_INT* LP_sp_index);
  void terminate();
  
  void changeToDense(bool forceChange=false);
  void display(FILE* fpout = stdout, char* printFormat = P_FORMAT);
  bool copyFrom(SparseLinearSpace& other);
  
  void setElement_SDP(SDPA_INT block, SDPA_INT nCol, SDPA_INT nRow, double ele);
  void setElement_SOCP(SDPA_INT block, SDPA_INT nCol, SDPA_INT nRow, double ele);
  void setElement_LP(SDPA_INT block, double ele);

  void setZero();
  void setIdentity(double scalar = 1.0);
  // no check
  bool sortSparseIndex(SDPA_INT&l , SDPA_INT& i, SDPA_INT& j);
};

class DenseLinearSpace
{
 public:
  SDPA_INT  SDP_nBlock;
  SDPA_INT  SOCP_nBlock;
  SDPA_INT  LP_nBlock;

  DenseMatrix* SDP_block;
  DenseMatrix* SOCP_block;
  double* LP_block;

  DenseLinearSpace();
  DenseLinearSpace(BlockStruct& bs);
  ~DenseLinearSpace();
  void initialize(BlockStruct& bs);
  void terminate();

  void display(FILE* fpout = stdout, char* printFormat = P_FORMAT);
  void displaySolution(BlockStruct& bs, FILE* fpout = stdout,
		       char* printFormat = P_FORMAT);
  bool copyFrom(DenseLinearSpace& other);
  void setElement_SDP(SDPA_INT block, SDPA_INT nCol, SDPA_INT nRow, double ele);
  void setElement_SOCP(SDPA_INT block, SDPA_INT nCol, SDPA_INT nRow, double ele);
  void setElement_LP(SDPA_INT block, double ele);
  void setZero();
  void setIdentity(double scalar = 1.0);
};

} // end of namespace 'sdpa'

#endif // __sdpa_struct_h__
