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

#include <sched.h>

#include "sdpa_newton.h"
#include "sdpa_parts.h"
#include "sdpa_jordan.h"
#include "sdpa_linear.h"  
#include "sdpa_algebra.h"
#include "ilasm.h"

namespace sdpa {

SDPA_INT  Newton::Column_Number = 0;
SDPA_INT  Newton::Calc_F1       = 0;
SDPA_INT  Newton::Calc_F2       = 0;

static int iteration;

Newton::Newton()
{
  useFormula = NULL;
  FormulaEstimate = NULL;
  useThread = NULL;
  bMat_type  = DENSE;

  // Caution: if SDPA doesn't use sparse bMat, 
  //          following variables are indefinite.
  this->SDP_nBlock = -1;
  SDP_number = NULL;  SDP_location_sparse_bMat = NULL;
  SDP_constraint1 = NULL;  SDP_constraint2 = NULL;
  SDP_blockIndex1 = NULL;  SDP_blockIndex2 = NULL;
  this->SOCP_nBlock = -1;
  SOCP_number = NULL;  SOCP_location_sparse_bMat = NULL;
  SOCP_constraint1 = NULL;  SOCP_constraint2 = NULL;
  SOCP_blockIndex1 = NULL;  SOCP_blockIndex2 = NULL;
  this->LP_nBlock = -1;
  LP_number = NULL;  LP_location_sparse_bMat = NULL;
  LP_constraint1 = NULL;  LP_constraint2 = NULL;
  LP_blockIndex1 = NULL;  LP_blockIndex2 = NULL;

  diagonalIndex = NULL;
  NUM_THREADS  = 1;
}

Newton::Newton(SDPA_INT m, BlockStruct& bs)
{
  initialize(m, bs);
}

Newton::~Newton()
{
  terminate();
}

void Newton::initialize(SDPA_INT m, BlockStruct& bs)
{
  gVec.initialize(m);

  SDP_nBlock  = bs.SDP_nBlock;
  SOCP_nBlock = bs.SOCP_nBlock;
  LP_nBlock   = bs.LP_nBlock;
  
  DxMat.initialize(bs);
  DyVec.initialize(m);
  DzMat.initialize(bs);
  r_zinvMat.initialize(bs);
  x_rd_zinvMat.initialize(bs);

  // Memory allocation of useFormula is moved to computeFormula_SDP
  // NewArray(useFormula,FormulaType,m*SDP_nBlock);


  bMat_type = DENSE;
  // Caution: if SDPA doesn't use sparse bMat, 
  //          following variables are indefinite.
  this->SDP_nBlock = -1;
  SDP_number = NULL;  SDP_location_sparse_bMat = NULL;
  SDP_constraint1 = NULL;  SDP_constraint2 = NULL;
  SDP_blockIndex1 = NULL;  SDP_blockIndex2 = NULL;
  this->SOCP_nBlock = -1;
  SOCP_number = NULL;  SOCP_location_sparse_bMat = NULL;
  SOCP_constraint1 = NULL;  SOCP_constraint2 = NULL;
  SOCP_blockIndex1 = NULL;  SOCP_blockIndex2 = NULL;
  this->LP_nBlock = -1;
  LP_number = NULL;  LP_location_sparse_bMat = NULL;
  LP_constraint1 = NULL;  LP_constraint2 = NULL;
  LP_blockIndex1 = NULL;  LP_blockIndex2 = NULL;

  diagonalIndex = NULL;
}

void Newton::terminate()
{

  if (bMat_type == SPARSE){

    if (SDP_location_sparse_bMat && SDP_constraint1 && SDP_constraint2
	&& SDP_blockIndex1 && SDP_blockIndex2) {
      for (SDPA_INT l=0; l<SDP_nBlock; ++l) {
	DeleteArray(SDP_location_sparse_bMat[l]);
	DeleteArray(SDP_constraint1[l]);
	DeleteArray(SDP_constraint2[l]);
	DeleteArray(SDP_blockIndex1[l]);
	DeleteArray(SDP_blockIndex2[l]);
      }
      DeleteArray(SDP_number);
      DeleteArray(SDP_location_sparse_bMat);
      DeleteArray(SDP_constraint1);
      DeleteArray(SDP_constraint2);
      DeleteArray(SDP_blockIndex1);
      DeleteArray(SDP_blockIndex2);
    }
#if 0
    if (SOCP_location_sparse_bMat && SOCP_constraint1 && SOCP_constraint2
	&& SOCP_blockIndex1 && SOCP_blockIndex2) {
      for (SDPA_INT l=0; l<SOCP_nBlock; ++l) {
	DeleteArray(SOCP_location_sparse_bMat[l]);
	DeleteArray(SOCP_constraint1[l]);
	DeleteArray(SOCP_constraint2[l]);
	DeleteArray(SOCP_blockIndex1[l]);
	DeleteArray(SOCP_blockIndex2[l]);
      }
      DeleteArray(SOCP_number);
      DeleteArray(SOCP_location_sparse_bMat);
      DeleteArray(SOCP_constraint1);
      DeleteArray(SOCP_constraint2);
      DeleteArray(SOCP_blockIndex1);
      DeleteArray(SOCP_blockIndex2);
    }
#endif
    if (LP_location_sparse_bMat && LP_constraint1 && LP_constraint2
	&& LP_blockIndex1 && LP_blockIndex2) {
      for (SDPA_INT l=0; l<LP_nBlock; ++l) {
	DeleteArray(LP_location_sparse_bMat[l]);
	DeleteArray(LP_constraint1[l]);
	DeleteArray(LP_constraint2[l]);
	DeleteArray(LP_blockIndex1[l]);
	DeleteArray(LP_blockIndex2[l]);
      }
      DeleteArray(LP_number);
      DeleteArray(LP_location_sparse_bMat);
      DeleteArray(LP_constraint1);
      DeleteArray(LP_constraint2);
      DeleteArray(LP_blockIndex1);
      DeleteArray(LP_blockIndex2);
    }

    DeleteArray(diagonalIndex);
    sparse_bMat.terminate();

  } else { // bMat_type == DENSE
    bMat.terminate();
  }

  SDPA_INT m = gVec.nDim;
  gVec.terminate();
  DxMat.terminate();
  DyVec.terminate();
  DzMat.terminate();
  r_zinvMat.terminate();
  x_rd_zinvMat.terminate();

  if (useFormula) {
    for (SDPA_INT j=0; j<m; ++j) {
      DeleteArray(useFormula[j]);
    }
    DeleteArray(useFormula);
  }

  if (FormulaEstimate) {
    for (SDPA_INT j=0; j<m; ++j) {
      DeleteArray(FormulaEstimate[j]);
    }
    DeleteArray(FormulaEstimate);
  }
  
  if (useThread) {
    for (SDPA_INT j=0; j<m; ++j) {
      DeleteArray(useThread[j]);
    }
    DeleteArray(useThread);
  }

}

void Newton::initialize_dense_bMat(SDPA_INT m)
{
  //  bMat_type = DENSE;
  //  printf("DENSE computations\n");
  bMat.initialize(m,m,DenseMatrix::DENSE);
}

  // 2008/03/12 kazuhide nakata
void Newton::initialize_sparse_bMat(SDPA_INT m)
{

  //  bMat_type = SPARSE;
  //  printf("SPARSE computation\n");

  // initialize sparse_bMat by Chordal::makeGraph
  //  sparse_bMat.display();

  bool isEmptyMatrix = false;
  // make index of diagonalIndex
  NewArray(diagonalIndex,SDPA_INT,m+1);
  SDPA_INT k=0;
  for (SDPA_INT index=0; index<sparse_bMat.NonZeroCount; index++){
    if (sparse_bMat.row_index[index] == sparse_bMat.column_index[index]) {
      diagonalIndex[k] = index;
      if (sparse_bMat.row_index[index] != k+1) {
	rMessage("The matrix [" << (sparse_bMat.row_index[index]-1)
		 << "] is empty");
	isEmptyMatrix = true;
	diagonalIndex[k+1] = diagonalIndex[k];
	k++;
      }
      k++;
    }
  }
  if (isEmptyMatrix) {
    rMessage("Input Data Error :: Some Input Matricies are Empty");
  }
    
  diagonalIndex[m] = sparse_bMat.NonZeroCount;
  #if 0
  rMessage("diagonalIndex = ");
  for (SDPA_INT index=0; index <m; ++index) {
    printf(" [%ld:%ld]",index,diagonalIndex[index]);
  }
  printf("\n");
  #endif
}

  // 2008/03/12 kazuhide nakata
void Newton::initialize_bMat(SDPA_INT m, Chordal& chordal,
			     InputData& inputData,
			     FILE* Display,
                             FILE* fpOut)
{
  /* Create clique tree */

  switch (chordal.best) {
  case SELECT_DENSE: 
    bMat_type = DENSE;
    if (Display) {
      fprintf(Display,"Schur computation : DENSE \n");
    }
    if (fpOut) {
      fprintf(fpOut,"Schur computation : DENSE \n");
    }
    initialize_dense_bMat(m);
    // Here, we release MUMPS and sparse_bMat
    chordal.terminate();
    break;
  case SELECT_MUMPS_BEST: 
    bMat_type = SPARSE;
    if (Display) {
      fprintf(Display,"Schur computation : SPARSE \n");
    }
    if (fpOut) {
      fprintf(fpOut,"Schur computation : SPARSE \n");
    }
    initialize_sparse_bMat(m);
    make_aggrigateIndex(inputData);
    break;
  default: 
    rError("Wrong Ordering Obtained");
    break;
  }

}


SDPA_INT Newton::binarySearchIndex(SDPA_INT i, SDPA_INT j)
{
  // binary search for index of sparse_bMat 
  SDPA_INT t = -1;
  // We store only the lower triangular
  SDPA_INT ii = i, jj = j;
  if (i<j) {
    jj = i;
    ii = j;
  }
  SDPA_INT begin  = diagonalIndex[jj]; 
  SDPA_INT end    = diagonalIndex[jj+1]-1;
  SDPA_INT target = (begin + end) / 2;
  while (end - begin > 1){
    if (sparse_bMat.row_index[target] < ii+1){
      begin = target;
      target = (begin + end) / 2;
    } else if (sparse_bMat.row_index[target] > ii+1){
      end = target;
      target = (begin + end) / 2;
    } else if (sparse_bMat.row_index[target] == ii+1) {
      t = target;
      break;
    }
  }
  if (t == -1){
    if (sparse_bMat.row_index[begin] == ii+1){
      t = begin;
    } else if (sparse_bMat.row_index[end] == ii+1){
      t = end;
    } else {
      #if 0
      SDPA_INT m = sparse_bMat.nRow;
      rMessage("Trouble ii = " << ii << " jj = " << j <<  " m = " << m);
      for (SDPA_INT k = 0; k<sparse_bMat.NonZeroCount; ++k) {
	fprintf(stdout,"[%04d,%04d] at %04d\n",
		sparse_bMat.row_index[k],
		sparse_bMat.column_index[k],
		k);
      }
      #endif
      // rError("Newton::make_aggrigateIndex program bug");
    }
  } 
  return t;
}
  

void Newton::make_aggrigateIndex_SDP(InputData& inputData)
{
  SDP_nBlock = inputData.SDP_nBlock;
  NewArray(SDP_number,SDPA_INT,SDP_nBlock);

  // memory allocate for aggrigateIndex
  NewArray(SDP_constraint1,SDPA_INT*,SDP_nBlock);
  NewArray(SDP_constraint2,SDPA_INT*,SDP_nBlock);
  NewArray(SDP_blockIndex1,SDPA_INT*,SDP_nBlock);
  NewArray(SDP_blockIndex2,SDPA_INT*,SDP_nBlock);
  NewArray(SDP_location_sparse_bMat,SDPA_INT*,SDP_nBlock);

  for (SDPA_INT l=0; l<SDP_nBlock; l++){
    const SDPA_INT size = (inputData.SDP_nConstraint[l] + 1) 
      * inputData.SDP_nConstraint[l] / 2;
    SDP_number[l] = size;
    NewArray(SDP_constraint1[l],SDPA_INT,size);
    NewArray(SDP_constraint2[l],SDPA_INT,size);
    NewArray(SDP_blockIndex1[l],SDPA_INT,size);
    NewArray(SDP_blockIndex2[l],SDPA_INT,size);
    NewArray(SDP_location_sparse_bMat[l],SDPA_INT,size);
  }

  for (SDPA_INT l = 0; l<SDP_nBlock; l++){
    SDPA_INT NonZeroCount = 0;

    for (SDPA_INT k1=0; k1<inputData.SDP_nConstraint[l]; k1++){
      SDPA_INT j = inputData.SDP_constraint[l][k1];
      SDPA_INT jb = inputData.SDP_blockIndex[l][k1];
      SDPA_INT jnz = inputData.A[j].SDP_sp_block[jb].NonZeroEffect;

      for (SDPA_INT k2=0; k2<inputData.SDP_nConstraint[l]; k2++){
	SDPA_INT i = inputData.SDP_constraint[l][k2];
	SDPA_INT ib = inputData.SDP_blockIndex[l][k2];
	SDPA_INT inz = inputData.A[i].SDP_sp_block[ib].NonZeroEffect;

	if ((jnz < inz) || ((inz == jnz) && (i < j))){
	  continue;
	}
	// set index which A_i and A_j are not zero matrix
	SDPA_INT target = binarySearchIndex(i,j);
	if (target == -1) {
	  rMessage("("<<(i+1)<<","<<(j+1)<<") might have index");
	  SDP_number[l]--;
	  continue;
	}
	SDP_constraint1[l][NonZeroCount] = i;
	SDP_constraint2[l][NonZeroCount] = j;
	SDP_blockIndex1[l][NonZeroCount] = ib;
	SDP_blockIndex2[l][NonZeroCount] = jb;
	SDP_location_sparse_bMat[l][NonZeroCount] = target;
	NonZeroCount++;
      }
    } // for k1
  } //for l  lth block
}


void Newton::make_aggrigateIndex_SOCP(InputData& inputData)
{
  SOCP_nBlock = inputData.SOCP_nBlock;
  NewArray(SOCP_number,SDPA_INT,SOCP_nBlock);
  if (SOCP_number == NULL) {
    rError("Newton::make_aggrigateIndex_SOCP memory exhausted ");
  }

  NewArray(SOCP_constraint1,SDPA_INT*,SOCP_nBlock);
  NewArray(SOCP_constraint2,SDPA_INT*,SOCP_nBlock);
  NewArray(SOCP_blockIndex1,SDPA_INT*,SOCP_nBlock);
  NewArray(SOCP_blockIndex2,SDPA_INT*,SOCP_nBlock);
  NewArray(SOCP_location_sparse_bMat,SDPA_INT*,SOCP_nBlock);

  for (SDPA_INT l=0; l<SOCP_nBlock; l++){
    SDPA_INT size = (inputData.SOCP_nConstraint[l]+1)
      * inputData.SOCP_nConstraint[l]/2;
    SOCP_number[l] = size;
    NewArray(SOCP_constraint1[l],SDPA_INT,size);
    NewArray(SOCP_constraint2[l],SDPA_INT,size);
    NewArray(SOCP_blockIndex1[l],SDPA_INT,size);
    NewArray(SOCP_blockIndex2[l],SDPA_INT,size);
    NewArray(SOCP_location_sparse_bMat[l],SDPA_INT,size);
  }

  for (SDPA_INT l = 0; l<SOCP_nBlock; l++){
    SDPA_INT NonZeroCount = 0;

    for (SDPA_INT k1=0; k1<inputData.SOCP_nConstraint[l]; k1++){
      SDPA_INT j = inputData.SOCP_constraint[l][k1];
      SDPA_INT jb = inputData.SOCP_blockIndex[l][k1];
      SDPA_INT jnz = inputData.A[j].SOCP_sp_block[jb].NonZeroEffect;

      for (SDPA_INT k2=0; k2<inputData.SOCP_nConstraint[l]; k2++){
	SDPA_INT i = inputData.SOCP_constraint[l][k2];
	SDPA_INT ib = inputData.SOCP_blockIndex[l][k2];
	SDPA_INT inz = inputData.A[i].SOCP_sp_block[ib].NonZeroEffect;

	if ((jnz < inz) || ((inz == jnz) && (i < j))){
	  continue;
	}

	// set index which A_i and A_j are not zero matrix
	SOCP_constraint1[l][NonZeroCount] = i;
	SOCP_constraint2[l][NonZeroCount] = j;
	SOCP_blockIndex1[l][NonZeroCount] = ib;
	SOCP_blockIndex2[l][NonZeroCount] = jb;

	SDPA_INT target = binarySearchIndex(i,j);
	SOCP_location_sparse_bMat[l][NonZeroCount] = target;
	NonZeroCount++;
      }
    } // for k1
  } //for l  lth block
}

void Newton::make_aggrigateIndex_LP(InputData& inputData)
{
  LP_nBlock = inputData.LP_nBlock;
  NewArray(LP_number,SDPA_INT,LP_nBlock);

  // memory allocate for aggrigateIndex
  NewArray(LP_constraint1,SDPA_INT*,LP_nBlock);
  NewArray(LP_constraint2,SDPA_INT*,LP_nBlock);
  NewArray(LP_blockIndex1,SDPA_INT*,LP_nBlock);
  NewArray(LP_blockIndex2,SDPA_INT*,LP_nBlock);
  NewArray(LP_location_sparse_bMat,SDPA_INT*,LP_nBlock);

  for (SDPA_INT l=0; l<LP_nBlock; l++){
    SDPA_INT size = (inputData.LP_nConstraint[l]+1)
      * inputData.LP_nConstraint[l]/2;
    LP_number[l] = size;
    NewArray(LP_constraint1[l],SDPA_INT,size);
    NewArray(LP_constraint2[l],SDPA_INT,size);
    NewArray(LP_blockIndex1[l],SDPA_INT,size);
    NewArray(LP_blockIndex2[l],SDPA_INT,size);
    NewArray(LP_location_sparse_bMat[l],SDPA_INT,size);
  }

  for (SDPA_INT l = 0; l<LP_nBlock; l++){
    SDPA_INT NonZeroCount = 0;

    for (SDPA_INT k1=0; k1<inputData.LP_nConstraint[l]; k1++){
      SDPA_INT j = inputData.LP_constraint[l][k1];
      SDPA_INT jb = inputData.LP_blockIndex[l][k1];

      for (SDPA_INT k2=0; k2<inputData.LP_nConstraint[l]; k2++){
	SDPA_INT i = inputData.LP_constraint[l][k2];
	SDPA_INT ib = inputData.LP_blockIndex[l][k2];

	if (i < j){
	  continue;
	}

	// set index which A_i and A_j are not zero matrix
	LP_constraint1[l][NonZeroCount] = i;
	LP_constraint2[l][NonZeroCount] = j;
	LP_blockIndex1[l][NonZeroCount] = ib;
	LP_blockIndex2[l][NonZeroCount] = jb;
	
	SDPA_INT target = binarySearchIndex(i,j);
	LP_location_sparse_bMat[l][NonZeroCount] = target;
	NonZeroCount++;
      }
    } // for k1
  } //for l  lth block
}

void Newton::make_aggrigateIndex(InputData& inputData)
{
  make_aggrigateIndex_SDP(inputData);
  //  make_aggrigateIndex_SOCP(inputData);
  make_aggrigateIndex_LP(inputData);
}

void Newton::computeFormula_SDP(InputData& inputData,
				double DenseRatio, double Kappa)
{
  SDPA_INT m = inputData.b.nDim;
  SDPA_INT SDP_nBlock = inputData.SDP_nBlock;

  NewArray(useFormula, FormulaType*, m);
  NewArray(FormulaEstimate, double*, m);
  NewArray(useThread, SDPA_INT*, m);

  for (SDPA_INT j=0; j<m; ++j) {
    NewArray(useFormula[j], FormulaType, inputData.A[j].SDP_sp_nBlock);
    //    printf("inputData.A[%ld].SDP_sp_nBlock = %ld, SDP_nBlock = %ld\n", j, inputData.A[j].SDP_sp_nBlock, SDP_nBlock);
  }

  for (SDPA_INT j=0; j<m; ++j) {
    NewArray(FormulaEstimate[j], double, inputData.A[j].SDP_sp_nBlock);
  }
  for (SDPA_INT j=0; j<m; ++j) {
    for (SDPA_INT k=0; k<inputData.A[j].SDP_sp_nBlock; ++k) {
      FormulaEstimate[j][k] = 0.0;
    }
  }

  for (SDPA_INT i=0; i<m; ++i) {
    NewArray(useThread[i], SDPA_INT, SDP_nBlock);
  }

  // need not to initialize useFormula

  
  SDPA_INT** upNonZeroCount;
  NewArray(upNonZeroCount, SDPA_INT*, m);
  for (SDPA_INT j=0; j<m; ++j) {
    NewArray(upNonZeroCount[j], SDPA_INT, inputData.A[j].SDP_sp_nBlock);
  }
  for (SDPA_INT j=0; j<m; ++j) {
    for (SDPA_INT jb=0; jb < inputData.A[j].SDP_sp_nBlock; ++jb) {
      upNonZeroCount[j][jb] = 0;
    }
  }

  SparseLinearSpace* A = inputData.A;

  #if 0
  for (SDPA_INT k=0; k<m; ++k) {
    for (SDPA_INT l=0; l<inputData.A[0].nBlock; ++l) {
      rMessage("A[" << k << "].ele[" << l << "] ="
	       << inputData.A[k].ele[l].NonZeroEffect);
    }
  }
  #endif

  // Count sum of number of elements
  // that each number of elements are less than own.

  for (SDPA_INT l=0; l<SDP_nBlock; ++l) {
    for (SDPA_INT k1=0; k1 < inputData.SDP_nConstraint[l];k1++){
      SDPA_INT j = inputData.SDP_constraint[l][k1];
      SDPA_INT jb = inputData.SDP_blockIndex[l][k1];
      SDPA_INT jnz = A[j].SDP_sp_block[jb].NonZeroEffect;
      SDPA_INT up = jnz;
      // rMessage("up = " << up);

      for (SDPA_INT k2=0; k2 < inputData.SDP_nConstraint[l];k2++){
	SDPA_INT i = inputData.SDP_constraint[l][k2];
	SDPA_INT ib = inputData.SDP_blockIndex[l][k2];
	SDPA_INT inz = A[i].SDP_sp_block[ib].NonZeroEffect;
	//		printf("%ld %ld %ld %ld %ld %ld\n",i,ib,inz, j, jb,jnz);
	if (inz < jnz) {
	  up += inz;
	}
#if 1
	else if ((jnz == inz) && i>j ) {
	  up += inz;
	}
#endif
      }
      upNonZeroCount[j][jb] = up;
      // rMessage("up = " << up);
    }
  }

  // Determine which formula
  double ff=0, ff1=0, ff2=0, ff3=0;
  Calc_F1 = Calc_F2 = 0;
  for (SDPA_INT l=0; l<SDP_nBlock; ++l) {
    SDPA_INT countf1,countf2,countf3;
    countf1 = countf2 = countf3 = 0;
    for (SDPA_INT k=0; k < inputData.SDP_nConstraint[l]; k++){
      //      printf("inputData.SDP_nConstraint[%ld] = %ld\n", l, inputData.SDP_nConstraint[l]);
      SDPA_INT j =  inputData.SDP_constraint[l][k];
      SDPA_INT jb =  inputData.SDP_blockIndex[l][k];
      double jnz = inputData.A[j].SDP_sp_block[jb].NonZeroEffect;

      double f1,f2,f3;
      double n       = inputData.A[j].SDP_sp_block[jb].nRow;
      double up      = upNonZeroCount[j][jb];

#if 1
      f1 = Kappa*n*jnz + n*n*n + Kappa*up;
      f2 = Kappa*n*jnz + Kappa*(n+1)*up;
      #if 0
      f3 = Kappa*(2*Kappa*jnz+1)*up/Kappa;
      #else
      f3 = Kappa*(2*Kappa*jnz+1)*up;
      #endif

#endif
      ff1 += f1;
      ff2 += f2;
      ff3 += f3;
      //rMessage("up = " << up << " nonzero = " << nonzero);
      //rMessage("f1=" << f1 << " f2=" << f2 << " f3=" << f3);
      //      printf("%ld %ld %lf %lf %lf\n",j,jb,f1,f2,f3);
      //printf("%ld %ld %le %le %le\n",k,l,f1,f2,f3);
      if (inputData.A[j].SDP_sp_block[jb].type == SparseMatrix::DENSE) {
	// if DENSE, we use only F1 or F2,
	// that is we don't use F3
	if (f1<f2) {
	  useFormula[j][jb] = F1;
	  //	  FormulaEstimate[j][jb] = f1;
	  countf1++;
	  ff += f1;
	} else {
	  useFormula[j][jb] = F2;
	  //	  FormulaEstimate[j][jb] = f2;
	  countf2++;
	  ff += f2;
	}
      } else {
	//printf("n = %ld, inz = %ld, n * n = %ld, inz / (n * n) = %lf\n", (int)n, (int)inz, (int)(n*n), (double)inz/(n*n));
	// this case is SPARSE
	if (f1<f2 && f1<f3) {
	  if ((n <= 200) && (2.0 * n >= jnz)) {
	    //rMessage("line " << k << " is F3:1");
	    useFormula[j][jb] = F3;
	    //	    FormulaEstimate[j][jb] = f3;
	    countf3++;
	    ff += f3;
	  }
	  else {
	    //rMessage("line " << k << " is F1");
	    useFormula[j][jb] = F1;
	    //	    FormulaEstimate[j][jb] = f1;
	    countf1++;
	    ff += f1;
	  }
	} else if (f2<f3) {
	  //rMessage("line " << k << " is F2");
	  useFormula[j][jb] = F2;
	  //	  FormulaEstimate[j][jb] = f2;
	  countf2++;
	  ff += f2;
	} else {
	  //rMessage("line " << k << " is F3:2");
	  useFormula[j][jb] = F3;
	  //	  FormulaEstimate[j][jb] = f3;
	  countf3++;
	  ff += f3;
	}
      }
    }

    Calc_F1 += countf1;
    Calc_F2 += countf2;
    // rMessage("Kappa = " << Kappa);
    #if 0
    rMessage("count f1 = " << countf1
	     << ":: count f2 = " << countf2
	     << ":: count f3 = " << countf3);
    #endif
  } // end of 'for (SDPA_INT l)'

  for (SDPA_INT j=0; j<m; ++j) {
    DeleteArray(upNonZeroCount[j]);
  }
  DeleteArray(upNonZeroCount);

  return;
}

void Newton::computeFormula_SDP_new(InputData& inputData)
{
  SDPA_INT m = inputData.b.nDim;
  SDPA_INT SDP_nBlock = inputData.SDP_nBlock;

  SDPA_INT RowNum, ThreadNum;
  SDPA_INT k1, j, jb;
  double *TotalTime;
  double Sum, TotalSum = 0.0;
  double *TmpEstimate;
  double **TmpEstimate2;

  NewArray(TotalTime, double, SDP_nBlock);
  NewArray(TmpEstimate, double, NUM_THREADS);
  NewArray(TmpEstimate2, double*, m);
  for (SDPA_INT i=0; i<m; ++i) {
    NewArray(TmpEstimate2[i], double, NUM_THREADS);
  }

  for (SDPA_INT i=0; i<m; ++i) {
    for (SDPA_INT k=0; k<NUM_THREADS; ++k) {
      TmpEstimate2[i][k] = 0.0;
    }
  }

  for (SDPA_INT i=0; i<m; ++i) {
    for (SDPA_INT k=0; k<SDP_nBlock; ++k) {
      useThread[i][k] = 0;
    }
  }

  //  printf("NUM_THREADS=%ld, SDP_nBlock = %ld, m = %ld\n", NUM_THREADS, SDP_nBlock, m);
  for (SDPA_INT l=0; l<SDP_nBlock; ++l) {

    for (SDPA_INT i=0; i<NUM_THREADS; ++i) {
      TmpEstimate[i] = 0.0;
    }

    TotalTime[l] = 0.0;
    for (k1 = 0; k1 < inputData.SDP_nConstraint[l]; ++k1) {
      j = inputData.SDP_constraint[l][k1];
      jb = inputData.SDP_blockIndex[l][k1];

      TotalTime[l] += FormulaEstimate[j][jb];
    }

    for (k1 = 0, ThreadNum = 0; k1 < inputData.SDP_nConstraint[l]; ++k1) {

      j = inputData.SDP_constraint[l][k1];
      jb = inputData.SDP_blockIndex[l][k1];
#if 0
      if (useFormula[j][jb] != F3){
	useThread[j][jb] = -1;
	continue;
      }
#endif
      //      printf("l = %ld, k1 = %ld, j = %ld, jb = %ld\n", l, k1, j, jb);
      //      printf("ThreadNum = %ld\n", ThreadNum);
      TmpEstimate[ThreadNum] += FormulaEstimate[j][jb];
      TmpEstimate2[j][ThreadNum] += FormulaEstimate[j][jb];
      //      printf("TotalTime[%ld] / NUM_THREADS = %lf, TimeEstimate[%ld] = %lf\n", l, TotalTime[l] / NUM_THREADS, ThreadNum, TmpEstimate[ThreadNum]);
      useThread[j][jb] = ThreadNum;
      //            printf("useThread[%ld][%ld] = %ld\n", j, jb, useThread[j][jb]);

      if (TmpEstimate[ThreadNum] > (TotalTime[l] / NUM_THREADS)) {
	ThreadNum++;
      }
    }

#if 0
    for (SDPA_INT i=0; i<m; ++i) {
      for (SDPA_INT k=0; k<NUM_THREADS; ++k) {
	//	printf("Block = %ld, TmpEstimate2[%ld][%ld] = %lf\n", l, i, k, TmpEstimate2[i][k]);
	printf("useThread[%ld][%ld] = %ld\n", i, k, useThread[i][k]);
      }
    }
#endif
  }

  if (TotalTime)
    DeleteArray(TotalTime);

  if (TmpEstimate) 
    DeleteArray(TmpEstimate);

  if (TmpEstimate2) {
    for (SDPA_INT j=0; j<m; ++j) {
      DeleteArray(TmpEstimate2[j]);
    }
    DeleteArray(TmpEstimate2);
  }

  return;
}


void Newton::compute_rMat(Newton::WHICH_DIRECTION direction,
			  AverageComplementarity& mu,
			  DirectionParameter& beta,
			  Solutions& currentPt,
			  WorkVariables& work)
{

  //     CORRECTOR ::  r_zinv = (-XZ -dXdZ + mu I)Z^{-1}
  // not CORRECTOR ::  r_zinv = (-XZ + mu I)Z^{-1}
  double target = beta.value*mu.current;
  Lal::let(r_zinvMat,'=',currentPt.invzMat,'*',&target);
  Lal::let(r_zinvMat,'=',r_zinvMat,'+',currentPt.xMat,&DMONE);

  if (direction == CORRECTOR) {
    // work.DLS1 = Dx Dz Z^{-1}
    Jal::ns_jordan_triple_product(work.DLS1,DxMat,DzMat,
				  currentPt.invzMat,work.DLS2);
    Lal::let(r_zinvMat,'=',r_zinvMat,'+',work.DLS1,&DMONE);
  }

  //  rMessage("r_zinvMat = ");
  //  r_zinvMat.display();
}

void Newton::Make_gVec(Newton::WHICH_DIRECTION direction,
		       InputData& inputData,
		       Solutions& currentPt,
		       Residuals& currentRes,
		       AverageComplementarity& mu,
		       DirectionParameter& beta,
		       Phase& phase,
		       WorkVariables& work,
		       ComputeTime& com)
{
  TimeStart(START1);
  // rMessage("mu = " << mu.current);
  // rMessage("beta = " << beta.value);
  compute_rMat(direction,mu,beta,currentPt,work);

  TimeEnd(END1);

  com.makerMat += TimeCal(START1,END1);

  TimeStart(START2);
  TimeStart(START_GVEC_MUL);

  // work.DLS1 = R Z^{-1} - X D Z^{-1} = r_zinv - X D Z^{-1}
  if (phase.value == SolveInfo:: pFEAS
      || phase.value == SolveInfo::noINFO) {

    if (direction == CORRECTOR) {
      // x_rd_zinvMat is computed in PREDICTOR step
      Lal::let(work.DLS1,'=',r_zinvMat,'+',x_rd_zinvMat,&DMONE);
    } else {
      // currentPt is infeasilbe, that is the residual
      // dualMat is not 0.
      //      x_rd_zinvMat = X D Z^{-1}
      Jal::ns_jordan_triple_product(x_rd_zinvMat,currentPt.xMat,
				    currentRes.dualMat,currentPt.invzMat,
				    work.DLS2);
      Lal::let(work.DLS1,'=',r_zinvMat,'+',x_rd_zinvMat,&DMONE);
    } // if (direction == CORRECTOR)

  } else {
    // dualMat == 0
    work.DLS1.copyFrom(r_zinvMat);
  }
  
  //  rMessage("work.DLS1");
  //  work.DLS1.display();

  TimeEnd(END_GVEC_MUL);
  com.makegVecMul += TimeCal(START_GVEC_MUL,END_GVEC_MUL);
    
  inputData.multi_InnerProductToA(work.DLS1,gVec);
  Lal::let(gVec,'=',gVec,'*',&DMONE);
  // rMessage("gVec =  ");
  // gVec.display();

  #if 0
  if (phase.value == SolveInfo:: dFEAS
      || phase.value == SolveInfo::noINFO) {
  #endif
    Lal::let(gVec,'=',gVec,'+',currentRes.primalVec);
  #if 0
  }
  #endif
  
  TimeEnd(END2);
  com.makegVec += TimeCal(START2,END2);
}

void Newton::calF1(double& ret, DenseMatrix& G,
		    SparseMatrix& Ai)
{
  Lal::let(ret,'=',Ai,'.',G);
}

void Newton::calF2(double& ret,
		   DenseMatrix& F, DenseMatrix& G,
		   DenseMatrix& invZ, SparseMatrix& Ai,
		    bool& hasF2Gcal)
{
  SDPA_INT alpha,beta;
  double value1,value2;

  SDPA_INT n    = Ai.nRow;
  // rMessage(" using F2 ");
  switch (Ai.type) {
  case SparseMatrix::SPARSE:
    // rMessage("F2::SPARSE  " << Ai.NonZeroCount);
    ret = 0.0;
    for (SDPA_INT index = 0; index < Ai.NonZeroCount; ++index) {
#if DATA_CAPSULE
      alpha  = Ai.DataS[index].vRow;
      beta   = Ai.DataS[index].vCol;
      value1 = Ai.DataS[index].vEle;
#else
      alpha  = Ai.row_index[index];
      beta   = Ai.column_index[index];
      value1 = Ai.sp_ele[index];
#endif      
      value2 = ddot_f77(&n, &F.de_ele[alpha+n*0], &n,
			&invZ.de_ele[0+n*beta], &IONE);
      ret += value1*value2;
      if (alpha!=beta) {
	value2 = ddot_f77(&n, &F.de_ele[beta+n*0], &n,
			  &invZ.de_ele[0+n*alpha], &IONE);
	ret += value1*value2;
      }
    }
    break;
  case SparseMatrix::DENSE:
    // G is temporary matrix
    // rMessage("F2::DENSE");
    if (hasF2Gcal == false) {
      // rMessage(" using F2 changing to F1");
      Lal::let(G,'=',F,'*',invZ);
      hasF2Gcal = true;
    }
    Lal::let(ret,'=',Ai,'.',G);
    break;
  } // end of switch
}

inline void Newton::calF3(double& ret,
		    DenseMatrix& X, DenseMatrix& invZ,
		    SparseMatrix& Ai, SparseMatrix& Aj)
{
  // Ai and Aj are SPARSE
  ret = 0.0;
  double sum;
  const SDPA_INT nCol = X.nCol;
  // rMessage("Aj.NonZeroCount = " << Aj.NonZeroCount);
  for (SDPA_INT index1=0; index1<Aj.NonZeroCount; ++index1) {
#if DATA_CAPSULE
    SDPA_INT alpha = Aj.DataS[index1].vRow;
    SDPA_INT beta  = Aj.DataS[index1].vCol;
    double value1 = Aj.DataS[index1].vEle;
#else
    SDPA_INT alpha = Aj.row_index[index1];
    SDPA_INT beta  = Aj.column_index[index1];
    double value1 = Aj.sp_ele[index1];
#endif
    sum = 0.0;
    for (SDPA_INT index2=0; index2<Ai.NonZeroCount; ++index2) {
#if DATA_CAPSULE
      SDPA_INT gamma = Ai.DataS[index2].vRow;
      SDPA_INT delta  = Ai.DataS[index2].vCol;
      double value2 = Ai.DataS[index2].vEle;
#else
      SDPA_INT gamma = Ai.row_index[index2];
      SDPA_INT delta  = Ai.column_index[index2];
      double value2 = Ai.sp_ele[index2];
#endif
      double plu = value2*invZ.de_ele[delta+nCol*beta]
        * X.de_ele[gamma+nCol*alpha];
      sum += plu;
      if (gamma!=delta) {
        double plu2 = value2*invZ.de_ele[gamma+nCol*beta]
          * X.de_ele[delta+nCol*alpha];
        sum += plu2;
      }
    }
    ret += value1*sum;
    if (alpha==beta) {
      continue;
    }
    sum = 0.0;
    for (SDPA_INT index2=0; index2<Ai.NonZeroCount; ++index2) {
#if DATA_CAPSULE
      SDPA_INT gamma = Ai.DataS[index2].vRow;
      SDPA_INT delta  = Ai.DataS[index2].vCol;
      double value2 = Ai.DataS[index2].vEle;
#else
      SDPA_INT gamma = Ai.row_index[index2];
      SDPA_INT delta  = Ai.column_index[index2];
      double value2 = Ai.sp_ele[index2];
#endif
      double plu = value2*invZ.de_ele[delta+nCol*alpha]
        * X.de_ele[gamma+nCol*beta];
      sum += plu;
      if (gamma!=delta) {
        double plu2 = value2*invZ.de_ele[gamma+nCol*alpha]
          * X.de_ele[delta+nCol*beta];
        sum += plu2;
      }
    }
    ret += value1*sum;
  } // end of 'for (index1)'
  return;
}

inline void Newton::calF1_thread(double& ret, DenseMatrix& G,
				 SparseMatrix& Ai)
{
  Lal::let(ret,'=',Ai,'.',G);
}

void Newton::calF2_thread(double& ret,
			  DenseMatrix& F, DenseMatrix& G,
			  DenseMatrix& invZ, SparseMatrix& Ai,
			  bool& hasF2Gcal)
{
  SDPA_INT alpha,beta;
  double value1,value2;

  SDPA_INT n    = Ai.nRow;
  // rMessage(" using F2 ");
  switch (Ai.type) {
  case SparseMatrix::SPARSE:
    // rMessage("F2::SPARSE  " << Aj.NonZeroCount);
    ret = 0.0;
    for (SDPA_INT index = 0; index < Ai.NonZeroCount; ++index) {
#if DATA_CAPSULE
      alpha  = Ai.DataS[index].vRow;
      beta   = Ai.DataS[index].vCol;
      value1 = Ai.DataS[index].vEle;
#else
      alpha  = Ai.row_index[index];
      beta   = Ai.column_index[index];
      value1 = Ai.sp_ele[index];
#endif
      value2 = ddot_f77(&n, &F.de_ele[alpha+n*0], &n,
			&invZ.de_ele[0+n*beta], &IONE);
      ret += value1*value2;
      if (alpha!=beta) {
	value2 = ddot_f77(&n, &F.de_ele[beta+n*0], &n,
			  &invZ.de_ele[0+n*alpha], &IONE);
	ret += value1*value2;
      }
    }
    break;
  case SparseMatrix::DENSE:
    // G is temporary matrix
    // rMessage("F2::DENSE");
    if (hasF2Gcal == false) {
      // rMessage(" using F2 changing to F1");
      Lal::let(G,'=',F,'*',invZ);
      hasF2Gcal = true;
    }
    Lal::let(ret,'=',Ai,'.',G);
    break;
  } // end of switch
}

#if 0
inline void Newton::calF3_thread(double& ret,
				 DenseMatrix& X, DenseMatrix& invZ,
				 SparseMatrix& Ai, SparseMatrix& Aj)
{
  // Ai and Aj are SPARSE
  ret = 0.0;
  double sum;
  const SDPA_INT nCol = X.nCol;
  // rMessage("Aj.NonZeroCount = " << Aj.NonZeroCount);
  for (SDPA_INT index1=0; index1<Aj.NonZeroCount; ++index1) {
#if DATA_CAPSULE
    SDPA_INT alpha = Aj.DataS[index1].vRow;
    SDPA_INT beta  = Aj.DataS[index1].vCol;
    double value1 = Aj.DataS[index1].vEle;
#else
    SDPA_INT alpha = Aj.row_index[index1];
    SDPA_INT beta  = Aj.column_index[index1];
    double value1 = Aj.sp_ele[index1];
#endif
    sum = 0.0;
    
    const SDPA_INT nCol = X.nCol;
    double *Xalpha, *Xbeta, *invZalpha, *invZbeta;
    Xalpha = &X.de_ele[nCol*alpha];
    Xbeta  = &X.de_ele[nCol*beta];
    invZalpha = &invZ.de_ele[nCol*alpha];
    invZbeta  = &invZ.de_ele[nCol*beta];
    
    for (SDPA_INT index2=0; index2<Ai.NonZeroCount; ++index2) {
#if DATA_CAPSULE
      SDPA_INT gamma = Ai.DataS[index2].vRow;
      SDPA_INT delta  = Ai.DataS[index2].vCol;
      double value2 = Ai.DataS[index2].vEle;
#else
      SDPA_INT gamma = Ai.row_index[index2];
      SDPA_INT delta  = Ai.column_index[index2];
      double value2 = Ai.sp_ele[index2];
#endif
      double plu = value2*invZbeta[delta] * Xalpha[gamma];

      sum += plu;
      if (gamma!=delta) {
        double plu2 = value2*invZbeta[gamma] * Xalpha[delta];
        sum += plu2;
      }
    }
    ret += value1*sum;
    if (alpha==beta) {
      continue;
    }
    sum = 0.0;
    for (SDPA_INT index2=0; index2<Ai.NonZeroCount; ++index2) {
#if DATA_CAPSULE
      SDPA_INT gamma = Ai.DataS[index2].vRow;
      SDPA_INT delta  = Ai.DataS[index2].vCol;
      double value2 = Ai.DataS[index2].vEle;
#else
      SDPA_INT gamma = Ai.row_index[index2];
      SDPA_INT delta  = Ai.column_index[index2];
      double value2 = Ai.sp_ele[index2];
#endif      
      double plu = value2*invZalpha[delta] * Xbeta[gamma];
      sum += plu;
      if (gamma!=delta) {
        double plu2 = value2*invZalpha[gamma] * Xbeta[delta];
        sum += plu2;
      }
    }
    ret += value1*sum;
  } // end of 'for (index1)'
  return;
}
#endif

inline void Newton::calF3_thread_1x1(double& ret,
				     DenseMatrix& X, DenseMatrix& invZ,
				     SparseMatrix& Ai, SparseMatrix& Aj)
{
  // Ai and Aj are SPARSE
  ret = 0.0;
  double sum = 0.0;
  SDPA_INT index1 = 0, index2 = 0;

  // rMessage("Aj.NonZeroCount = " << Aj.NonZeroCount);
#if DATA_CAPSULE
  SDPA_INT alpha = Aj.DataS[index1].vRow;
  SDPA_INT beta  = Aj.DataS[index1].vCol;
  double value1 = Aj.DataS[index1].vEle;
  SDPA_INT gamma = Ai.DataS[index2].vRow;
  SDPA_INT delta  = Ai.DataS[index2].vCol;
  double value2 = Ai.DataS[index2].vEle;
#else
  SDPA_INT alpha = Aj.row_index[index1];
  SDPA_INT beta  = Aj.column_index[index1];
  double value1 = Aj.sp_ele[index1];
  SDPA_INT gamma = Ai.row_index[index2];
  SDPA_INT delta  = Ai.column_index[index2];
  double value2 = Ai.sp_ele[index2];
#endif

  const SDPA_INT nCol = X.nCol;
  double *Xalpha, *Xbeta, *invZalpha, *invZbeta;
  Xalpha = &X.de_ele[nCol*alpha];
  Xbeta  = &X.de_ele[nCol*beta];
  invZalpha = &invZ.de_ele[nCol*alpha];
  invZbeta  = &invZ.de_ele[nCol*beta];
  
  double plu = value2*invZbeta[delta] * Xalpha[gamma];
  sum += plu;
  if (gamma!=delta) {
    double plu2 = value2*invZbeta[gamma] * Xalpha[delta];
    sum += plu2;
  }
  ret += value1*sum;
  if (alpha==beta) {
    return;
  }
  sum = 0.0;
  plu = value2*invZalpha[delta] * Xbeta[gamma];
  sum += plu;
  if (gamma!=delta) {
    double plu2 = value2*invZalpha[gamma] * Xbeta[delta];
    sum += plu2;
  }
  ret += value1*sum;

  return;
}

#if 0
inline void Newton::calF3_thread(double& ret,
				 DenseMatrix& X, DenseMatrix& invZ,
				 SparseMatrix& Ai, SparseMatrix& Aj)
{
  // Ai and Aj are SPARSE
  ret = 0.0;
  SDPA_INT index1, index2;
  SDPA_INT alpha, beta, gamma, delta;
  double value1, value2, plu, plu2;
  double sum;
  const SDPA_INT nCol = X.nCol;
  double *Xalpha, *Xbeta, *invZalpha, *invZbeta;
  for (index1=0; index1<Aj.NonZeroCount; ++index1) {
    alpha = Aj.DataS[index1].vRow;
    beta  = Aj.DataS[index1].vCol;
    value1 = Aj.DataS[index1].vEle;
    sum = 0.0;
    
    Xalpha = &X.de_ele[nCol*alpha];
    Xbeta  = &X.de_ele[nCol*beta];
    invZalpha = &invZ.de_ele[nCol*alpha];
    invZbeta  = &invZ.de_ele[nCol*beta];
    
    for (index2=0; index2<Ai.NonZeroCount; ++index2) {
      #if DATA_CAPSULE
      gamma = Ai.DataS[index2].vRow;
      delta  = Ai.DataS[index2].vCol;
      value2 = Ai.DataS[index2].vEle;
      #else
      gamma = Ai.row_index[index2];
      delta  = Ai.column_index[index2];
      value2 = Ai.sp_ele[index2];
      #endif
      plu = value2*invZbeta[delta] * Xalpha[gamma];
      sum += plu;
      if (gamma!=delta) {
        plu2 = value2*invZbeta[gamma] * Xalpha[delta];
        sum += plu2;
      }

      if (alpha==beta) {
	continue;
      }

      plu = value2*invZalpha[delta] * Xbeta[gamma];
      sum += plu;
      if (gamma!=delta) {
        plu2 = value2*invZalpha[gamma] * Xbeta[delta];
        sum += plu2;
      }
    }
    ret += value1*sum;
  } // end of 'for (index1)'
  return;
}
#endif

#if 1
inline void Newton::calF3_thread(double& ret,
				 DenseMatrix& X, DenseMatrix& invZ,
				 SparseMatrix& Ai, SparseMatrix& Aj)
{
  // Ai and Aj are SPARSE
  ret = 0.0;
  double sum;
  double *Xalpha, *Xbeta, *invZalpha, *invZbeta;
  const SDPA_INT nCol = X.nCol;
  //  rMessage("Aj.NonZeroCount = " << Aj.NonZeroCount);
  for (SDPA_INT index1=0; index1<Aj.NonZeroCount; ++index1) {
    #if DATA_CAPSULE
    SDPA_INT alpha = Aj.DataS[index1].vRow; 
    SDPA_INT beta  = Aj.DataS[index1].vCol;
    double value1 = Aj.DataS[index1].vEle;
    #else
    SDPA_INT alpha = Aj.row_index[index1]; 
    SDPA_INT beta  = Aj.column_index[index1];
    double value1 = Aj.sp_ele[index1];
    #endif
    Xalpha =    &X.de_ele[nCol*alpha];
    Xbeta  =    &X.de_ele[nCol*beta];
    invZalpha = &invZ.de_ele[nCol*alpha];
    invZbeta  = &invZ.de_ele[nCol*beta];
    sum = 0.0;
    //    rMessage("Ai.NonZeroCount = " << Ai.NonZeroCount);
    for (SDPA_INT index2=0; index2<Ai.NonZeroCount; ++index2) {
      #if DATA_CAPSULE
      SDPA_INT gamma = Ai.DataS[index2].vRow;
      SDPA_INT delta  = Ai.DataS[index2].vCol;
      double value2 = Ai.DataS[index2].vEle;
      #else
      SDPA_INT gamma = Ai.row_index[index2]; 
      SDPA_INT delta = Ai.column_index[index2];
      double value2 = Ai.sp_ele[index2];
      #endif
      double plu = value2 * invZbeta[delta] * Xalpha[gamma];
      sum += plu;
      if (gamma!=delta) {
        double plu2 = value2 * invZbeta[gamma] * Xalpha[delta];
        sum += plu2;
      }
    }
    ret += value1*sum;
    if (alpha==beta) {
      continue;
    }
    sum = 0.0;
    for (SDPA_INT index2=0; index2<Ai.NonZeroCount; ++index2) {
      #if DATA_CAPSULE
      SDPA_INT gamma = Ai.DataS[index2].vRow;
      SDPA_INT delta  = Ai.DataS[index2].vCol;
      double value2 = Ai.DataS[index2].vEle;
      #else
      SDPA_INT gamma = Ai.row_index[index2]; 
      SDPA_INT delta  = Ai.column_index[index2];
      double value2 = Ai.sp_ele[index2];
      #endif
      double plu = value2 * invZalpha[delta] * Xbeta[gamma];
      sum += plu;
      if (gamma!=delta) {
        double plu2 = value2 * invZalpha[gamma] * Xbeta[delta];
        sum += plu2;
      }
    }
    ret += value1*sum;
  } // end of 'for (index1)'
  return;
}
#endif

inline void Newton::calF3_thread_2(double& ret,
				   DenseMatrix& X, DenseMatrix& invZ,
				   SparseMatrix& Ai, SparseMatrix& Aj)
{
  // Ai and Aj are SPARSE
  ret = 0.0;
  double sum = 0.0;
  SDPA_INT index1, index2, alpha, beta, gamma, delta;
  double value1, plu;

  double *Xalpha, *Xbeta, *invZalpha, *invZbeta;
  const SDPA_INT nCol = X.nCol;

  // rMessage("Aj.NonZeroCount = " << Aj.NonZeroCount);
  for (index1=0; index1<Aj.NonZeroCount; ++index1) {
    #if DATA_CAPSULE
    alpha  = Aj.DataS[index1].vRow;
    beta   = Aj.DataS[index1].vCol;
    value1 = Aj.DataS[index1].vEle;
    #else
    alpha = Aj.row_index[index1]; 
    beta  = Aj.column_index[index1];
    value1 = Aj.sp_ele[index1];
    #endif

    Xalpha    = &X.de_ele[nCol*alpha];
    Xbeta     = &X.de_ele[nCol*beta];
    invZalpha = &invZ.de_ele[nCol*alpha];
    invZbeta  = &invZ.de_ele[nCol*beta];

    if (alpha != beta) {
      sum = 0.0;
      for (index2=0; index2 < Ai.NonZeroCount; index2++) {
	#if DATA_CAPSULE
	gamma = Ai.DataS[index2].vRow;
	delta  = Ai.DataS[index2].vCol;
	#else
	gamma = Ai.row_index[index2];
	delta  = Ai.column_index[index2];
	#endif
	
	if (gamma != delta) {
	  plu = invZbeta [delta] * Xalpha[gamma]
	      + invZbeta [gamma] * Xalpha[delta]
	      + invZalpha[delta] * Xbeta [gamma]
	      + invZalpha[gamma] * Xbeta [delta];
	} else {
	  plu = invZbeta [delta] * Xalpha[gamma]
	      + invZalpha[delta] * Xbeta [gamma];
	}
	#if DATA_CAPSULE
	sum += Ai.DataS[index2].vEle * plu;
	#else
	sum += Ai.sp_ele[index2] * plu;
	#endif
      }
      ret += value1 * sum;
    } else {
      sum = 0.0;
      for (index2=0; index2 < Ai.NonZeroCount; index2++) {
	#if DATA_CAPSULE
	gamma = Ai.DataS[index2].vRow;
	delta  = Ai.DataS[index2].vCol;
	#else
	gamma = Ai.row_index[index2];
	delta  = Ai.column_index[index2];
	#endif
	
	if (gamma != delta) {
	  plu = invZbeta[delta] * Xalpha[gamma]
	      + invZbeta[gamma] * Xalpha[delta];
	} else {
	  plu = invZbeta[delta] * Xalpha[gamma];
	}
	#if DATA_CAPSULE
	sum += Ai.DataS[index2].vEle * plu;
	#else
	sum += Ai.sp_ele[index2] * plu;
	#endif	
      }
      ret += value1 * sum;
    }
  }
  return;
}

void Newton::compute_bMat_dense_SDP_thread(InputData& inputData,
                                    Solutions& currentPt,
                                    WorkVariables& work,
                                    ComputeTime& com)
{
  // rMessage("NUM_THREADS = " << NUM_THREADS);
  pthread_t*  handle;
  NewArray(handle,pthread_t,NUM_THREADS);
  thread_arg_t* targ;
  NewArray(targ,thread_arg_t,NUM_THREADS);

  SDPA_INT m = currentPt.mDim;
  SDPA_INT SDP_nBlock = inputData.SDP_nBlock;

  iteration++;
  if (iteration >= 2)
    computeFormula_SDP_new(inputData);

  for (SDPA_INT k=0; k<NUM_THREADS; k++) {
    targ[k].Num_of_Threads  = NUM_THREADS;
    targ[k].mDIM            = m;
    targ[k].SDP_nBlock      = SDP_nBlock;
    targ[k].bMat            = &bMat;
    targ[k].useFormula      = useFormula;
    targ[k].useThread       = useThread;
    targ[k].FormulaEstimate = FormulaEstimate;
    targ[k].inputData       = &inputData;
    targ[k].currentPt       = &currentPt;
    targ[k].work            = &work;
    targ[k].com             = &com;
  }

  //  fprintf(stderr, "hoge  .. %ld\n", NUM_THREADS);

  for (SDPA_INT l=0; l<SDP_nBlock; l++) {
    Column_Number = 0;

    for (SDPA_INT k=0; k<NUM_THREADS; k++) {
      targ[k].Block_Number = l;
      targ[k].thread_num   = k;
    }
    
    if (iteration >= 2) {
#pragma omp parallel for schedule(dynamic)
      for (SDPA_INT k=0; k<NUM_THREADS; k++) {
	compute_bMat_dense_SDP_thread_func_new(&targ[k]);
      } 
    }
    else {
#pragma omp parallel for schedule(dynamic)
      for (SDPA_INT k=0; k<NUM_THREADS; k++) {
	compute_bMat_dense_SDP_thread_func(&targ[k]);
      }
    }
  }

  DeleteArray(handle);
  DeleteArray(targ);
}

void* Newton::compute_bMat_dense_SDP_thread_func(void *arg)
{
  SDPA_INT l, m, n;
  SDPA_INT k, k1;
  SDPA_INT SDP_nBlock;
  DenseMatrix work1, work2;
  static struct timeval B_NDIAG_START1, B_NDIAG_END1;

  FILE *fp;
  char filename[100];
  SDPA_INT t1, t2, ts;

  thread_arg_t *targ = (thread_arg_t *)arg;
 
  l = targ->Block_Number;
  k = targ->thread_num;
  m = targ->mDIM;
  n = targ->Num_of_Threads;
  SDP_nBlock = targ->SDP_nBlock;

  //  sprintf(filename, "output_thread_%ld.txt", k);
  //  if ((fp=fopen(filename,"a")) == NULL) {
  //    rError("Cannot Open Output File " << filename);
  //  }
  //  fprintf(fp, "iteration = %d, thread_num = %ld, get_apicid = %d\n\n", iteration, k, get_apicid());

  //  printf("targ-> Block_Number = %ld\n", targ-> Block_Number); 

  // DenseMatrix& xMat = targ->currentPt->xMat.SDP_block[l];
  // DenseMatrix& invzMat = targ->currentPt->invzMat.SDP_block[l];
  work1.initialize(targ->work->DLS1.SDP_block[l].nRow,
		   targ->work->DLS1.SDP_block[l].nCol,
		   DenseMatrix::DENSE);
  work2.initialize(targ->work->DLS2.SDP_block[l].nRow,
		   targ->work->DLS2.SDP_block[l].nCol,
		   DenseMatrix::DENSE);

  //  k1 = k;
  while(1) {
    t1 = rdtsc();
#if 1
#pragma omp critical
    k1 = Column_Number++;
#endif

    //    printf("inputData->SDP_nConstraint[%ld] = %ld\n", l, targ->inputData->SDP_nConstraint[l]);
    if (k1 >= targ->inputData->SDP_nConstraint[l])
      break;

    SDPA_INT j = targ->inputData->SDP_constraint[l][k1];
    SDPA_INT jb = targ->inputData->SDP_blockIndex[l][k1];

    FormulaType formula = targ->useFormula[j][jb];
    if ((formula==F1) || (formula==F2)) {
      //      k1 += n;
      continue;
    }

    SDPA_INT jnz = targ->inputData->A[j].SDP_sp_block[jb].NonZeroEffect;
    SparseMatrix& Aj = targ->inputData->A[j].SDP_sp_block[jb];

    TimeStart(B_NDIAG_START1);
#if 0
    TimeStart(B_NDIAG_START2);

    bool hasF2Gcal = false;

    if (formula==F1) {
#pragma omp critical
      {
	Lal::let(work1,'=',targ->currentPt->xMat.SDP_block[l],'*',Aj);
	Lal::let(work2,'=',work1,'*',targ->currentPt->invzMat.SDP_block[l]);
      }
    } else if (formula==F2) {
#pragma omp critical
      {
      // Lal::let(work1,'=',Ai,'*',targ->currentPt->invzMat.SDP_block[l]);
	Lal::let(work1,'=',targ->currentPt->xMat.SDP_block[l],'*',Aj);
      }
      hasF2Gcal = false;
    }
    TimeEnd(B_NDIAG_END2);
    targ->com->B_PRE += TimeCal(B_NDIAG_START2,B_NDIAG_END2);
#endif

    for (SDPA_INT k2=targ->inputData->SDP_nConstraint[l]-1; k2 >= 0; k2--) {
      SDPA_INT i = targ->inputData->SDP_constraint[l][k2];
      SDPA_INT ib = targ->inputData->SDP_blockIndex[l][k2];
      SDPA_INT inz = targ->inputData->A[i].SDP_sp_block[ib].NonZeroEffect;
      SparseMatrix& Ai = targ->inputData->A[i].SDP_sp_block[ib];
           
      if ((jnz < inz) || ( (inz == jnz) && (i<j))) {
	continue;
      }
      double value;
      if ((Ai.NonZeroCount == 1) && (Aj.NonZeroCount == 1)) {
	calF3_thread_1x1(value,
			   targ->currentPt->xMat.SDP_block[l],
			   targ->currentPt->invzMat.SDP_block[l],
			   Ai,Aj);
      }
      else {
	calF3_thread_2(value,
			 targ->currentPt->xMat.SDP_block[l],
			 targ->currentPt->invzMat.SDP_block[l],
			 Ai,Aj);
      }

#if 0
      switch (formula) {
      case F3:
	// rMessage("calF3");
#if DEBUG
	printf("F3 in %ld\n", targ->thread_num);
#endif
	if ((Ai.NonZeroCount == 1) && (Aj.NonZeroCount == 1))
	  calF3_thread_1x1(value,
			   targ->currentPt->xMat.SDP_block[l],
			   targ->currentPt->invzMat.SDP_block[l],
			   Ai,Aj);
	else
	  calF3_thread_2(value,
			 targ->currentPt->xMat.SDP_block[l],
			 targ->currentPt->invzMat.SDP_block[l],
			 Ai,Aj);
	break;
      case F1:
	// rMessage("calF1");
#if DEBUG
	printf("F1 in %ld\n", targ->thread_num);
#endif
	calF1_thread(value,work2,Ai);
	break;
      case F2:
	// rMessage("calF2 ");
#if DEBUG
	printf("F2 in %ld\n", targ->thread_num);
#endif
	calF2_thread(value,work1,work2,
		     targ->currentPt->invzMat.SDP_block[l],Ai,hasF2Gcal);
	break;
      } // end of switch
#endif
      if (i!=j) {
	targ->bMat->de_ele[i+m*j] += value;
	targ->bMat->de_ele[j+m*i] += value;
      } else {
	targ->bMat->de_ele[i+m*i] += value;
      }
    } // end of 'for (SDPA_INT j)'

    TimeEnd(B_NDIAG_END1);
    double t = TimeCal(B_NDIAG_START1,B_NDIAG_END1);
    //    printf("targ->com->B_F3 = %lf, t = %lf\n", targ->com->B_F3, t);
#pragma omp critical
    targ->com->B_F3 += t;
#if 0
    switch (formula) {
    case F1: targ->com->B_F1 += t; break;
    case F2: targ->com->B_F2 += t; break;
    case F3: targ->com->B_F3 += t; break;
    }
#endif

    t2 = rdtsc();
    ts = t2 - t1;

    targ->FormulaEstimate[j][jb] += ts;

    //    fprintf(stdout, "j = %ld, jb = %ld, Estimate = %lf, Time = %ld\n", j, jb, targ->FormulaEstimate[j][jb], ts);
    //    k1 += n;
  }

  //  fprintf(fp, "\n");
  //  fclose(fp);

  work1.terminate();
  work2.terminate();

  return NULL;
}

void* Newton::compute_bMat_dense_SDP_thread_func_new(void *arg)
{
  SDPA_INT l, m, n;
  SDPA_INT k, k1;
  SDPA_INT SDP_nBlock;
  DenseMatrix work1, work2;
  static struct timeval B_NDIAG_START1, B_NDIAG_END1;

  //  FILE *fp;
  //  char filename[100];
  SDPA_INT t1, t2, ts;

  thread_arg_t *targ = (thread_arg_t *)arg;
 
  l = targ->Block_Number;
  k = targ->thread_num;
  m = targ->mDIM;
  n = targ->Num_of_Threads;
  SDP_nBlock = targ->SDP_nBlock;

  //  sprintf(filename, "output_thread_%ld.txt", k);
  //  if ((fp=fopen(filename,"a")) == NULL) {
  //    rError("Cannot Open Output File " << filename);
  //  }
  // DenseMatrix& xMat = targ->currentPt->xMat.SDP_block[l];
  // DenseMatrix& invzMat = targ->currentPt->invzMat.SDP_block[l];
  work1.initialize(targ->work->DLS1.SDP_block[l].nRow,
		   targ->work->DLS1.SDP_block[l].nCol,
		   DenseMatrix::DENSE);
  work2.initialize(targ->work->DLS2.SDP_block[l].nRow,
		   targ->work->DLS2.SDP_block[l].nCol,
		   DenseMatrix::DENSE);

  k1 = -1;
  while(1) {
    t1 = rdtsc();
    k1++;

    //    printf("inputData->SDP_nConstraint[%ld] = %ld\n", l, targ->inputData->SDP_nConstraint[l]);
    if (k1 >= targ->inputData->SDP_nConstraint[l])
      break;

    SDPA_INT j = targ->inputData->SDP_constraint[l][k1];
    SDPA_INT jb = targ->inputData->SDP_blockIndex[l][k1];

    FormulaType formula = targ->useFormula[j][jb];
    if ((formula==F1) || (formula==F2) || targ->useThread[j][jb] != k) {
      //      k1 += n;
      continue;
    }

    SDPA_INT jnz = targ->inputData->A[j].SDP_sp_block[jb].NonZeroEffect;
    SparseMatrix& Aj = targ->inputData->A[j].SDP_sp_block[jb];

    TimeStart(B_NDIAG_START1);
#if 0
    TimeStart(B_NDIAG_START2);

    bool hasF2Gcal = false;

    if (formula==F1) {
#pragma omp critical
      {
	Lal::let(work1,'=',targ->currentPt->xMat.SDP_block[l],'*',Aj);
	Lal::let(work2,'=',work1,'*',targ->currentPt->invzMat.SDP_block[l]);
      }
    } else if (formula==F2) {
#pragma omp critical
      {
      // Lal::let(work1,'=',Ai,'*',targ->currentPt->invzMat.SDP_block[l]);
	Lal::let(work1,'=',targ->currentPt->xMat.SDP_block[l],'*',Aj);
      }
      hasF2Gcal = false;
    }
    TimeEnd(B_NDIAG_END2);
    targ->com->B_PRE += TimeCal(B_NDIAG_START2,B_NDIAG_END2);
#endif

    for (SDPA_INT k2=targ->inputData->SDP_nConstraint[l]-1; k2 >= 0; k2--) {
      SDPA_INT i = targ->inputData->SDP_constraint[l][k2];
      SDPA_INT ib = targ->inputData->SDP_blockIndex[l][k2];
      SDPA_INT inz = targ->inputData->A[i].SDP_sp_block[ib].NonZeroEffect;
      SparseMatrix& Ai = targ->inputData->A[i].SDP_sp_block[ib];
           
      if ((jnz < inz) || ( (inz == jnz) && (i<j))) {
	continue;
      }
      double value;
#if 0
      if ((Ai.NonZeroCount == 1) && (Aj.NonZeroCount == 1)) {
	calF3_thread_1x1(value,
			   targ->currentPt->xMat.SDP_block[l],
			   targ->currentPt->invzMat.SDP_block[l],
			   Ai,Aj);
      }
      else {
	calF3_thread_2(value,
			 targ->currentPt->xMat.SDP_block[l],
			 targ->currentPt->invzMat.SDP_block[l],
			 Ai,Aj);
      }
#endif
      calF3_thread_2(value,
		     targ->currentPt->xMat.SDP_block[l],
		     targ->currentPt->invzMat.SDP_block[l],
		     Ai,Aj);

#if 0
      switch (formula) {
      case F3:
	// rMessage("calF3");
#if DEBUG
	printf("F3 in %ld\n", targ->thread_num);
#endif
	if ((Ai.NonZeroCount == 1) && (Aj.NonZeroCount == 1))
	  calF3_thread_1x1(value,
			   targ->currentPt->xMat.SDP_block[l],
			   targ->currentPt->invzMat.SDP_block[l],
			   Ai,Aj);
	else
	  calF3_thread_2(value,
			 targ->currentPt->xMat.SDP_block[l],
			 targ->currentPt->invzMat.SDP_block[l],
			 Ai,Aj);
	break;
      case F1:
	// rMessage("calF1");
#if DEBUG
	printf("F1 in %ld\n", targ->thread_num);
#endif
	calF1_thread(value,work2,Ai);
	break;
      case F2:
	// rMessage("calF2 ");
#if DEBUG
	printf("F2 in %ld\n", targ->thread_num);
#endif
	calF2_thread(value,work1,work2,
		     targ->currentPt->invzMat.SDP_block[l],Ai,hasF2Gcal);
	break;
      } // end of switch
#endif
      if (i!=j) {
	targ->bMat->de_ele[i+m*j] += value;
	targ->bMat->de_ele[j+m*i] += value;
      } else {
	targ->bMat->de_ele[i+m*i] += value;
      }
    } // end of 'for (SDPA_INT j)'

    TimeEnd(B_NDIAG_END1);
    double t = TimeCal(B_NDIAG_START1,B_NDIAG_END1);
    //    printf("targ->com->B_F3 = %lf, t = %lf\n", targ->com->B_F3, t);
#pragma omp critical
    targ->com->B_F3 += t;
#if 0
    switch (formula) {
    case F1: targ->com->B_F1 += t; break;
    case F2: targ->com->B_F2 += t; break;
    case F3: targ->com->B_F3 += t; break;
    }
#endif

    t2 = rdtsc();
    ts = t2 - t1;

    targ->FormulaEstimate[j][jb] += ts;

    //    fprintf(stdout, "j = %ld, jb = %ld, Estimate = %lf, Time = %ld\n", j, jb, targ->FormulaEstimate[j][jb], ts);
  }

  //  fprintf(fp, "\n");
  //  fclose(fp);

  work1.terminate();
  work2.terminate();

  return NULL;
}

void Newton::compute_bMat_dense_SDP(InputData& inputData,
                                    Solutions& currentPt,
                                    WorkVariables& work,
                                    ComputeTime& com)
{
  SDPA_INT m = currentPt.mDim;
  SDPA_INT SDP_nBlock = inputData.SDP_nBlock;

  for (SDPA_INT l=0; l<SDP_nBlock; ++l) {
    DenseMatrix& xMat = currentPt.xMat.SDP_block[l];
    DenseMatrix& invzMat = currentPt.invzMat.SDP_block[l];
    DenseMatrix& work1 = work.DLS1.SDP_block[l];
    DenseMatrix& work2 = work.DLS2.SDP_block[l];

    for (SDPA_INT k1=0; k1<inputData.SDP_nConstraint[l]; k1++) {
      SDPA_INT j = inputData.SDP_constraint[l][k1];
      SDPA_INT jb = inputData.SDP_blockIndex[l][k1];

      FormulaType formula = useFormula[j][jb];
      if (formula==F3)
	continue;

      SDPA_INT jnz = inputData.A[j].SDP_sp_block[jb].NonZeroEffect;
      SparseMatrix& Aj = inputData.A[j].SDP_sp_block[jb];

      // ---------------------------------------------------
      // formula = F3; // this is force change
      // ---------------------------------------------------
      TimeStart(B_NDIAG_START1);
      TimeStart(B_NDIAG_START2);

      bool hasF2Gcal = false;
      if (formula==F1) {
	// Lal::let(work1,'=',Ai,'*',invzMat);
	// Lal::let(work2,'=',xMat,'*',work1);
	Lal::let(work1,'=',xMat,'*',Aj);
	Lal::let(work2,'=',work1,'*',invzMat);
      } else if (formula==F2) {
	// Lal::let(work1,'=',Ai,'*',invzMat);
	Lal::let(work1,'=',xMat,'*',Aj);
	hasF2Gcal = false;
	// Lal::let(gMat.ele[l],'=',xMat.ele[l],'*',fMat.ele[l]);
      }
      TimeEnd(B_NDIAG_END2);
      com.B_PRE += TimeCal(B_NDIAG_START2,B_NDIAG_END2);

#pragma omp parallel for schedule(dynamic)
      for (SDPA_INT k2=0; k2<inputData.SDP_nConstraint[l]; k2++) {
	SDPA_INT i = inputData.SDP_constraint[l][k2];
	SDPA_INT ib = inputData.SDP_blockIndex[l][k2];
	SDPA_INT inz = inputData.A[i].SDP_sp_block[ib].NonZeroEffect;
	SparseMatrix& Ai = inputData.A[i].SDP_sp_block[ib];
          
	// Select the formula A[i] or the formula A[j].
	// Use formula that has more NonZeroEffects than others.
	// We must calculate i==j.
          
	if ((jnz < inz) || ( (inz == jnz) && (i<j))) {
	  continue;
	}

	double value;
	switch (formula) {
	case F3:
	  // rMessage("calF3");
	  calF3(value,xMat,invzMat,Ai,Aj);
	  break;
	case F1:
	  // rMessage("calF1");
	  calF1(value,work2,Ai);
	  break;
	case F2:
	  // rMessage("calF2 ");
	  calF2(value,work1,work2,invzMat,Ai,hasF2Gcal);
	  // calF1(value2,gMat.ele[l],A[j].ele[l]);
	  // rMessage("calF2:  " << (value-value2));
	  break;
	} // end of switch
	if (i!=j) {
	  bMat.de_ele[i+m*j] += value;
	  bMat.de_ele[j+m*i] += value;
	} else {
	  bMat.de_ele[i+m*i] += value;
	}
      } // end of 'for (SDPA_INT j)'

      TimeEnd(B_NDIAG_END1);
      double t = TimeCal(B_NDIAG_START1,B_NDIAG_END1);
      switch (formula) {
      case F1: com.B_F1 += t; break;
      case F2: com.B_F2 += t; break;
      case F3: com.B_F3 += t; break;
      }
    } // end of 'for (SDPA_INT i)'
  } // end of 'for (SDPA_INT l)'
}

void Newton::compute_bMat_dense_SDP123(InputData& inputData,
                                    Solutions& currentPt,
                                    WorkVariables& work,
                                    ComputeTime& com)
{
  SDPA_INT m = currentPt.mDim;
  SDPA_INT SDP_nBlock = inputData.SDP_nBlock;

  for (SDPA_INT l=0; l<SDP_nBlock; ++l) {
    DenseMatrix& xMat = currentPt.xMat.SDP_block[l];
    DenseMatrix& invzMat = currentPt.invzMat.SDP_block[l];
    DenseMatrix& work1 = work.DLS1.SDP_block[l];
    DenseMatrix& work2 = work.DLS2.SDP_block[l];

    for (SDPA_INT k1=0; k1<inputData.SDP_nConstraint[l]; k1++) {
      SDPA_INT j = inputData.SDP_constraint[l][k1];
      SDPA_INT jb = inputData.SDP_blockIndex[l][k1];

      FormulaType formula = useFormula[j][jb];

      SDPA_INT jnz = inputData.A[j].SDP_sp_block[jb].NonZeroEffect;
      SparseMatrix& Aj = inputData.A[j].SDP_sp_block[jb];

      // ---------------------------------------------------
      // formula = F3; // this is force change
      // ---------------------------------------------------
      TimeStart(B_NDIAG_START1);
      TimeStart(B_NDIAG_START2);

      bool hasF2Gcal = false;
      if (formula==F1) {
	// Lal::let(work1,'=',Ai,'*',invzMat);
	// Lal::let(work2,'=',xMat,'*',work1);
	Lal::let(work1,'=',xMat,'*',Aj);
	Lal::let(work2,'=',work1,'*',invzMat);
      } else if (formula==F2) {
	// Lal::let(work1,'=',Ai,'*',invzMat);
	Lal::let(work1,'=',xMat,'*',Aj);
	hasF2Gcal = false;
	// Lal::let(gMat.ele[l],'=',xMat.ele[l],'*',fMat.ele[l]);
      }
      TimeEnd(B_NDIAG_END2);
      com.B_PRE += TimeCal(B_NDIAG_START2,B_NDIAG_END2);

#pragma omp parallel for schedule(dynamic)
      for (SDPA_INT k2=0; k2<inputData.SDP_nConstraint[l]; k2++) {
	SDPA_INT i = inputData.SDP_constraint[l][k2];
	SDPA_INT ib = inputData.SDP_blockIndex[l][k2];
	SDPA_INT inz = inputData.A[i].SDP_sp_block[ib].NonZeroEffect;
	SparseMatrix& Ai = inputData.A[i].SDP_sp_block[ib];
          
	// Select the formula A[i] or the formula A[j].
	// Use formula that has more NonZeroEffects than others.
	// We must calculate i==j.
          
	if ((jnz < inz) || ( (inz == jnz) && (i<j))) {
	  continue;
	}

	double value;
	switch (formula) {
	case F3:
	  // rMessage("calF3");
	  calF3(value,xMat,invzMat,Ai,Aj);
	  break;
	case F1:
	  // rMessage("calF1");
	  calF1(value,work2,Ai);
	  break;
	case F2:
	  // rMessage("calF2 ");
	  calF2(value,work1,work2,invzMat,Ai,hasF2Gcal);
	  // calF1(value2,gMat.ele[l],A[j].ele[l]);
	  // rMessage("calF2:  " << (value-value2));
	  break;
	} // end of switch
	if (i!=j) {
	  bMat.de_ele[i+m*j] += value;
	  bMat.de_ele[j+m*i] += value;
	} else {
	  bMat.de_ele[i+m*i] += value;
	}
      } // end of 'for (SDPA_INT j)'

      TimeEnd(B_NDIAG_END1);
      double t = TimeCal(B_NDIAG_START1,B_NDIAG_END1);
      switch (formula) {
      case F1: com.B_F1 += t; break;
      case F2: com.B_F2 += t; break;
      case F3: com.B_F3 += t; break;
      }
    } // end of 'for (SDPA_INT i)'
  } // end of 'for (SDPA_INT l)'
}

void Newton::compute_bMat_dense_SDP3(InputData& inputData,
                                    Solutions& currentPt,
                                    WorkVariables& work,
                                    ComputeTime& com)
{
  SDPA_INT m = currentPt.mDim;
  SDPA_INT SDP_nBlock = inputData.SDP_nBlock;

  for (SDPA_INT l=0; l<SDP_nBlock; ++l) {
    DenseMatrix& xMat = currentPt.xMat.SDP_block[l];
    DenseMatrix& invzMat = currentPt.invzMat.SDP_block[l];
    DenseMatrix& work1 = work.DLS1.SDP_block[l];
    DenseMatrix& work2 = work.DLS2.SDP_block[l];

#pragma omp parallel for schedule(dynamic)
    for (SDPA_INT k1=0; k1<inputData.SDP_nConstraint[l]; k1++) {
      SDPA_INT j = inputData.SDP_constraint[l][k1];
      SDPA_INT jb = inputData.SDP_blockIndex[l][k1];

      FormulaType formula = useFormula[j][jb];

      SDPA_INT jnz = inputData.A[j].SDP_sp_block[jb].NonZeroEffect;
      SparseMatrix& Aj = inputData.A[j].SDP_sp_block[jb];

      // ---------------------------------------------------
      // formula = F3; // this is force change
      // ---------------------------------------------------
      TimeStart(B_NDIAG_START1);
      if (formula != F3)
	continue;

      bool hasF2Gcal = false;
#if 0
      TimeStart(B_NDIAG_START2);

      if (formula==F1) {
	// Lal::let(work1,'=',Ai,'*',invzMat);
	// Lal::let(work2,'=',xMat,'*',work1);
	Lal::let(work1,'=',xMat,'*',Aj);
	Lal::let(work2,'=',work1,'*',invzMat);
      } else if (formula==F2) {
	// Lal::let(work1,'=',Ai,'*',invzMat);
	Lal::let(work1,'=',xMat,'*',Aj);
	hasF2Gcal = false;
	// Lal::let(gMat.ele[l],'=',xMat.ele[l],'*',fMat.ele[l]);
      }
      TimeEnd(B_NDIAG_END2);
      com.B_PRE += TimeCal(B_NDIAG_START2,B_NDIAG_END2);
#endif

      for (SDPA_INT k2=0; k2<inputData.SDP_nConstraint[l]; k2++) {
	SDPA_INT i = inputData.SDP_constraint[l][k2];
	SDPA_INT ib = inputData.SDP_blockIndex[l][k2];
	SDPA_INT inz = inputData.A[i].SDP_sp_block[ib].NonZeroEffect;
	SparseMatrix& Ai = inputData.A[i].SDP_sp_block[ib];
          
	// Select the formula A[i] or the formula A[j].
	// Use formula that has more NonZeroEffects than others.
	// We must calculate i==j.
          
	if ((jnz < inz) || ( (inz == jnz) && (i<j))) {
	  continue;
	}

	double value;
	switch (formula) {
	case F3:
	  // rMessage("calF3");
#if 0
	  if ((Ai.NonZeroCount == 1) && (Aj.NonZeroCount == 1)) {
	    calF3_thread_2(value,xMat,invzMat,Ai,Aj);
	  }
	  else {
	    calF3_thread_2(value,xMat,invzMat,Ai,Aj);
	  }
#endif
	  calF3_thread_2(value,xMat,invzMat,Ai,Aj);
	  break;
	case F1:
	  // rMessage("calF1");
	  calF1(value,work2,Ai);
	  break;
	case F2:
	  // rMessage("calF2 ");
	  calF2(value,work1,work2,invzMat,Ai,hasF2Gcal);
	  // calF1(value2,gMat.ele[l],A[j].ele[l]);
	  // rMessage("calF2:  " << (value-value2));
	  break;
	} // end of switch
	if (i!=j) {
	  bMat.de_ele[i+m*j] += value;
	  bMat.de_ele[j+m*i] += value;
	} else {
	  bMat.de_ele[i+m*i] += value;
	}
      } // end of 'for (SDPA_INT j)'

      TimeEnd(B_NDIAG_END1);
      double t = TimeCal(B_NDIAG_START1,B_NDIAG_END1);
      switch (formula) {
      case F1: com.B_F1 += t; break;
      case F2: com.B_F2 += t; break;
      case F3: com.B_F3 += t; break;
      }
    } // end of 'for (SDPA_INT i)'
  } // end of 'for (SDPA_INT l)'
}

void Newton::setNumThreads(FILE* Display, FILE* fpOut, SDPA_INT NumThreads)
{
  if (NumThreads == 0) { // Automatic from OMP_NUM_THREADS
    char* env1      = NULL;
    env1 = getenv("OMP_NUM_THREADS");
    if (env1 != NULL) {
      NUM_THREADS = atoi(env1);
    }
    else {
      NUM_THREADS = 1;
    }
  }
  else {
    NUM_THREADS = NumThreads;
  }
  if (Display) {
    fprintf(Display,"NumThreads  is set as %d\n", NUM_THREADS);
  }
  if (fpOut) {
    fprintf(fpOut,  "NumThreads  is set as %d\n", NUM_THREADS);
  }
  
}

void Newton::compute_bMat_sparse_SDP_thread(InputData& inputData,
                                    Solutions& currentPt,
                                    WorkVariables& work,
                                    ComputeTime& com)
{
  pthread_t*  handle;
  NewArray(handle,pthread_t,NUM_THREADS);
  thread_arg_t* targ;
  NewArray(targ,thread_arg_t,NUM_THREADS);

  SDPA_INT m = currentPt.mDim;
  SDPA_INT SDP_nBlock = inputData.SDP_nBlock;

  iteration++;
  if (iteration >= 2)
    computeFormula_SDP_new(inputData);

  for (SDPA_INT k=0; k<NUM_THREADS; k++) {
    targ[k].mDIM            = m;
    targ[k].SDP_nBlock      = SDP_nBlock;
    targ[k].SDP_number      = SDP_number;
    targ[k].sparse_bMat     = &sparse_bMat;
    targ[k].SDP_constraint1 = SDP_constraint1;
    targ[k].SDP_constraint2 = SDP_constraint2;
    targ[k].SDP_blockIndex1 = SDP_blockIndex1;
    targ[k].SDP_blockIndex2 = SDP_blockIndex2;
    targ[k].SDP_location_sparse_bMat = SDP_location_sparse_bMat;    
    targ[k].useFormula      = useFormula;
    targ[k].useThread       = useThread;
    targ[k].FormulaEstimate = FormulaEstimate;
    targ[k].inputData       = &inputData;
    targ[k].currentPt       = &currentPt;
    targ[k].work            = &work;
    targ[k].com             = &com;
   }

  for (SDPA_INT l=0; l<SDP_nBlock; l++) {
    Column_Number = 0;

    for (SDPA_INT k=0; k<NUM_THREADS; k++) {
      targ[k].Block_Number = l;
      targ[k].thread_num   = k;
    }

    if (iteration >= 2) {
#pragma omp parallel for schedule(dynamic)
      for (SDPA_INT k=0; k<NUM_THREADS; k++) {
	compute_bMat_sparse_SDP_thread_func_new(&targ[k]);
      }
    }
    else {
#pragma omp parallel for schedule(dynamic)
      for (SDPA_INT k=0; k<NUM_THREADS; k++) {
	compute_bMat_sparse_SDP_thread_func(&targ[k]);
      }
    }
  }

  DeleteArray(handle);
  DeleteArray(targ);
}

void* Newton::compute_bMat_sparse_SDP_thread_func(void *arg)
{
  SDPA_INT l;
  SDPA_INT iter;
  SDPA_INT SDP_nBlock;
  DenseMatrix work1, work2;
  SDPA_INT previous_j = -1;

  thread_arg_t *targ = (thread_arg_t *)arg;
 
  SDPA_INT t1, t2, ts;

  l          = targ->Block_Number;
  SDP_nBlock = targ->SDP_nBlock;

  //  printf("targ-> Block_Number = %ld\n", targ-> Block_Number); 

  // DenseMatrix& xMat = targ->currentPt->xMat.SDP_block[l];
  // DenseMatrix& invzMat = targ->currentPt->invzMat.SDP_block[l];
  work1.initialize(targ->work->DLS1.SDP_block[l].nRow,
		   targ->work->DLS1.SDP_block[l].nCol,
		   DenseMatrix::DENSE);
  work2.initialize(targ->work->DLS2.SDP_block[l].nRow,
		   targ->work->DLS2.SDP_block[l].nCol,
		   DenseMatrix::DENSE);

  TimeStart(B_NDIAG_START1);

  while(1) {
    t1 = rdtsc();
#pragma omp critical
    iter = Column_Number++;

    if (iter >= targ->SDP_number[l])
      break;

    SDPA_INT j = targ->SDP_constraint2[l][iter];
    SDPA_INT jb = targ->SDP_blockIndex2[l][iter];
    //    printf("j = %ld, jb = %ld\n", j, jb);
    SparseMatrix& Aj = targ->inputData->A[j].SDP_sp_block[jb];
    
    FormulaType formula = targ->useFormula[j][jb];
    // rMessage("FormulaType = " << formula);

    if (j != previous_j){
      TimeStart(B_NDIAG_START2);

      if (formula==F1) {
#pragma omp critical
	{
	Lal::let(work1,'=',targ->currentPt->xMat.SDP_block[l],'*',Aj);
	Lal::let(work2,'=',work1,'*',targ->currentPt->invzMat.SDP_block[l]);
	}
      } else if (formula==F2) {
#pragma omp critical
	{
	// Lal::let(work1,'=',Ai,'*',targ->currentPt->invzMat.SDP_block[l]);
	Lal::let(work1,'=',targ->currentPt->xMat.SDP_block[l],'*',Aj);
	}
      }
      TimeEnd(B_NDIAG_END2);
      targ->com->B_PRE += TimeCal(B_NDIAG_START2,B_NDIAG_END2);
    }
    
    SDPA_INT i = targ->SDP_constraint1[l][iter];
    SDPA_INT ib = targ->SDP_blockIndex1[l][iter];
    //    printf("i = %ld, ib = %ld\n", i, ib);
    SparseMatrix& Ai = targ->inputData->A[i].SDP_sp_block[ib];

    double value;
    bool dummyHasF2Gcal = true;
    switch (formula) {
    case F3:
      //      rMessage("calF3");
      if ((Ai.NonZeroCount == 1) && (Aj.NonZeroCount == 1))
	calF3_thread_1x1(value,
			 targ->currentPt->xMat.SDP_block[l],
			 targ->currentPt->invzMat.SDP_block[l],
			 Ai,Aj);
      else
	calF3_thread_2(value,
		       targ->currentPt->xMat.SDP_block[l],
		       targ->currentPt->invzMat.SDP_block[l],
		       Ai,Aj);
      break;
    case F1:
      //      rMessage("calF1");
      calF1_thread(value,work2,Ai);
      break;
    case F2:
      //      rMessage("calF2 ");
      calF2_thread(value,work1,work2,
		   targ->currentPt->invzMat.SDP_block[l],Ai,dummyHasF2Gcal);
      break;
    } // end of switch
    targ->sparse_bMat->sp_ele[targ->SDP_location_sparse_bMat[l][iter]] += value;
    previous_j = j;

    t2 = rdtsc();
    ts = t2 - t1;

    targ->FormulaEstimate[j][jb] += ts;

    //    fprintf(stdout, "j = %ld, jb = %ld, Estimate = %lf, Time = %ld\n", j, jb, targ->FormulaEstimate[j][jb], ts);    
  }
#if 0
  TimeEnd(B_NDIAG_END1);
  double t = TimeCal(B_NDIAG_START1,B_NDIAG_END1);
  switch (formula) {
  case F1: targ->com->B_F1 += t; break;
  case F2: targ->com->B_F2 += t; break;
  case F3: targ->com->B_F3 += t; break;
  }
#endif

  work1.terminate();
  work2.terminate();

  return NULL;
}

void* Newton::compute_bMat_sparse_SDP_thread_func_new(void *arg)
{
  SDPA_INT l, k;
  SDPA_INT iter;
  SDPA_INT SDP_nBlock;
  DenseMatrix work1, work2;
  SDPA_INT previous_j = -1;
  thread_arg_t *targ = (thread_arg_t *)arg;

  SDPA_INT t1, t2, ts;
 
  l          = targ->Block_Number;
  k          = targ->thread_num;
  SDP_nBlock = targ->SDP_nBlock;

  //  printf("targ-> Block_Number = %ld\n", targ-> Block_Number); 

  // DenseMatrix& xMat = targ->currentPt->xMat.SDP_block[l];
  // DenseMatrix& invzMat = targ->currentPt->invzMat.SDP_block[l];
  work1.initialize(targ->work->DLS1.SDP_block[l].nRow,
		   targ->work->DLS1.SDP_block[l].nCol,
		   DenseMatrix::DENSE);
  work2.initialize(targ->work->DLS2.SDP_block[l].nRow,
		   targ->work->DLS2.SDP_block[l].nCol,
		   DenseMatrix::DENSE);

  TimeStart(B_NDIAG_START1);

  iter = -1;
  while(1) {
    t1 = rdtsc();
    iter++;

    if (iter >= targ->SDP_number[l])
      break;

    SDPA_INT j = targ->SDP_constraint2[l][iter];
    SDPA_INT jb = targ->SDP_blockIndex2[l][iter];
    if (targ->useThread[j][jb] != k) {
      continue;
    }

    SparseMatrix& Aj = targ->inputData->A[j].SDP_sp_block[jb];
    FormulaType formula = targ->useFormula[j][jb];
    // rMessage("FormulaType = " << formula);

    if (j != previous_j){
      TimeStart(B_NDIAG_START2);

      if (formula==F1) {
#pragma omp critical
	{
	Lal::let(work1,'=',targ->currentPt->xMat.SDP_block[l],'*',Aj);
	Lal::let(work2,'=',work1,'*',targ->currentPt->invzMat.SDP_block[l]);
	}
      } else if (formula==F2) {
#pragma omp critical
	{
	// Lal::let(work1,'=',Ai,'*',targ->currentPt->invzMat.SDP_block[l]);
	Lal::let(work1,'=',targ->currentPt->xMat.SDP_block[l],'*',Aj);
	}
      }
      TimeEnd(B_NDIAG_END2);
      targ->com->B_PRE += TimeCal(B_NDIAG_START2,B_NDIAG_END2);
    }
    
    SDPA_INT i = targ->SDP_constraint1[l][iter];
    SDPA_INT ib = targ->SDP_blockIndex1[l][iter];
    SparseMatrix& Ai = targ->inputData->A[i].SDP_sp_block[ib];

    double value;
    bool dummyHasF2Gcal = true;
    switch (formula) {
    case F3:
      //      rMessage("calF3");
      if ((Ai.NonZeroCount == 1) && (Aj.NonZeroCount == 1))
	calF3_thread_1x1(value,
			 targ->currentPt->xMat.SDP_block[l],
			 targ->currentPt->invzMat.SDP_block[l],
			 Ai,Aj);
      else
	calF3_thread_2(value,
		       targ->currentPt->xMat.SDP_block[l],
		       targ->currentPt->invzMat.SDP_block[l],
		       Ai,Aj);
      break;
    case F1:
      //      rMessage("calF1");
      calF1_thread(value,work2,Ai);
      break;
    case F2:
      //      rMessage("calF2 ");
      calF2_thread(value,work1,work2,
		   targ->currentPt->invzMat.SDP_block[l],Ai,dummyHasF2Gcal);
      break;
    } // end of switch
    targ->sparse_bMat->sp_ele[targ->SDP_location_sparse_bMat[l][iter]] += value;
    previous_j = j;

    t2 = rdtsc();
    ts = t2 - t1;

    targ->FormulaEstimate[j][jb] += ts;

    //    fprintf(stdout, "j = %ld, jb = %ld, Estimate = %lf, Time = %ld\n", j, jb, targ->FormulaEstimate[j][jb], ts);       
  }
#if 0
  TimeEnd(B_NDIAG_END1);
  double t = TimeCal(B_NDIAG_START1,B_NDIAG_END1);
  switch (formula) {
  case F1: targ->com->B_F1 += t; break;
  case F2: targ->com->B_F2 += t; break;
  case F3: targ->com->B_F3 += t; break;
  }
#endif

  work1.terminate();
  work2.terminate();

  return NULL;
}

void Newton::compute_bMat_sparse_SDP(InputData& inputData,
				     Solutions& currentPt,
				     WorkVariables& work,
				     ComputeTime& com)
{
  TimeStart(B_NDIAG_START1);
  TimeStart(B_NDIAG_START2);

  for (SDPA_INT l=0; l<SDP_nBlock; ++l) {
    DenseMatrix& xMat = currentPt.xMat.SDP_block[l];
    DenseMatrix& invzMat = currentPt.invzMat.SDP_block[l];
    DenseMatrix& work1 = work.DLS1.SDP_block[l];
    DenseMatrix& work2 = work.DLS2.SDP_block[l];
    SDPA_INT previous_j = -1;

    for (SDPA_INT iter = 0; iter < SDP_number[l]; iter++){
      //      TimeStart(B_NDIAG_START1);
      SDPA_INT j = SDP_constraint2[l][iter];
      SDPA_INT jb = SDP_blockIndex2[l][iter];
      SparseMatrix& Aj = inputData.A[j].SDP_sp_block[jb];
      FormulaType formula = useFormula[j][jb];
      
      if (j != previous_j){
	// ---------------------------------------------------
	// formula = F3; // this is force change
	// ---------------------------------------------------
	TimeStart(B_NDIAG_START2);
	
	if (formula==F1) {
	  // Lal::let(work1,'=',Ai,'*',invzMat);
	  // Lal::let(work2,'=',xMat,'*',work1);
	  Lal::let(work1,'=',xMat,'*',Aj);
	  Lal::let(work2,'=',work1,'*',invzMat);
	} else if (formula==F2) {
	  // Lal::let(work1,'=',Ai,'*',invzMat);
	  Lal::let(work1,'=',xMat,'*',Aj);
	  // Lal::let(gMat.ele[l],'=',xMat.ele[l],'*',fMat.ele[l]);
	}
	TimeEnd(B_NDIAG_END2);
	com.B_PRE += TimeCal(B_NDIAG_START2,B_NDIAG_END2);
      }
      
      SDPA_INT i = SDP_constraint1[l][iter];
      SDPA_INT ib = SDP_blockIndex1[l][iter];
      SparseMatrix& Ai = inputData.A[i].SDP_sp_block[ib];

      #if 0
      printf("iter = %ld, loc = %ld, i = %ld, ib = %ld, j = %ld, jb = %ld\n",
	     iter, SDP_location_sparse_bMat[l][iter], i, ib, j, jb);
      #endif
      
      double value;
      bool dummyHasF2Gcal = true;
      switch (formula) {
      case F3:
	// rMessage("calF3");
	if ((Ai.NonZeroCount == 1) && (Aj.NonZeroCount == 1))
	  calF3_thread_1x1(value,xMat,invzMat,Ai,Aj);
	else
	  calF3_thread(value,xMat,invzMat,Ai,Aj);
        // calF3(value,work1,work2,xMat,invzMat,Ai,Aj);
	break;
      case F1:
	// rMessage("calF1");
	calF1(value,work2,Ai);
	break;
      case F2:
	// rMessage("calF2 ");
	calF2(value,work1,work2,invzMat,Ai,dummyHasF2Gcal);
	// calF1(value2,gMat.ele[l],A[j].ele[l]);
	// rMessage("calF2:  " << (value-value2));
	break;
      } // end of switch
      sparse_bMat.sp_ele[SDP_location_sparse_bMat[l][iter]] += value;
      previous_j = j;
    } // end of 'for (SDPA_INT index)'
#if 0
    TimeEnd(B_NDIAG_END1);
    double t = TimeCal(B_NDIAG_START1,B_NDIAG_END1);
    switch (formula) {
    case F1: com.B_F1 += t; break;
    case F2: com.B_F2 += t; break;
    case F3: com.B_F3 += t; break;
    }
#endif
  } // end of 'for (SDPA_INT l)'
}

#if 0
void Newton::compute_bMat_dense_SCOP(InputData& inputData,
				     Solutions& currentPt,
				     WorkVariables& work,
				     ComputeTime& com)
{
    rError("current version does not support SOCP");
}

void Newton::compute_bMat_sparse_SOCP(InputData& inputData,
				      Solutions& currentPt,
				      WorkVariables& work,
				      ComputeTime& com)
{
    rError("current version does not support SOCP");
}
#endif

void Newton::compute_bMat_dense_LP(InputData& inputData,
				   Solutions& currentPt,
				   WorkVariables& work,
				   ComputeTime& com)
{
  SDPA_INT m = currentPt.mDim;
  SDPA_INT LP_nBlock = inputData.LP_nBlock;

  TimeEnd(B_DIAG_START1);
  for (SDPA_INT l=0; l<LP_nBlock; ++l) {
    double xMat = currentPt.xMat.LP_block[l];
    double invzMat = currentPt.invzMat.LP_block[l];

#pragma omp parallel for schedule(dynamic)
      for (SDPA_INT k1=0; k1<inputData.LP_nConstraint[l]; k1++) {
	SDPA_INT j = inputData.LP_constraint[l][k1];
	SDPA_INT jb = inputData.LP_blockIndex[l][k1];
	//	SDPA_INT inz = inputData.A[i].LP_sp_block[ib].NonZeroEffect;
	double Aj = inputData.A[j].LP_sp_block[jb];

	for (SDPA_INT k2=k1; k2<inputData.LP_nConstraint[l]; k2++) {
	  SDPA_INT i = inputData.LP_constraint[l][k2];
	  SDPA_INT ib = inputData.LP_blockIndex[l][k2];
	  //	  SDPA_INT jnz = inputData.A[j].LP_sp_block[jb].NonZeroEffect;
	  double Ai = inputData.A[i].LP_sp_block[ib];

	  double value;
	  value = xMat * invzMat * Ai * Aj;

	  if (i!=j) {
	    bMat.de_ele[i+m*j] += value;
	    bMat.de_ele[j+m*i] += value;
	  } else {
	    bMat.de_ele[i+m*i] += value;
	  }
	} // end of 'for (SDPA_INT j)'
      } // end of 'for (SDPA_INT i)'
  } // end of 'for (SDPA_INT l)'
  TimeEnd(B_DIAG_END1);
  com.B_DIAG += TimeCal(B_DIAG_START1,B_DIAG_END1);
}

void Newton::compute_bMat_sparse_LP(InputData& inputData,
				    Solutions& currentPt,
				    WorkVariables& work,
				    ComputeTime& com)
{
  TimeEnd(B_DIAG_START1);
  for (SDPA_INT l=0; l<LP_nBlock; ++l) {
    double xMat = currentPt.xMat.LP_block[l];
    double invzMat = currentPt.invzMat.LP_block[l];
    
    for (SDPA_INT iter = 0; iter < LP_number[l]; iter++){
      SDPA_INT j = LP_constraint2[l][iter];
      SDPA_INT jb = LP_blockIndex2[l][iter];
      double Aj = inputData.A[j].LP_sp_block[jb];

      SDPA_INT i = LP_constraint1[l][iter];
      SDPA_INT ib = LP_blockIndex1[l][iter];
      double Ai = inputData.A[i].LP_sp_block[ib];
      
      double value;
      value = xMat * invzMat * Ai * Aj;
      sparse_bMat.sp_ele[LP_location_sparse_bMat[l][iter]] += value;
    } // end of 'for (SDPA_INT iter)
  } // end of 'for (SDPA_INT l)'
  TimeEnd(B_DIAG_END1);
  com.B_DIAG += TimeCal(B_DIAG_START1,B_DIAG_END1);
}



void Newton::Make_bMat(InputData& inputData,
		       Solutions& currentPt,
		       WorkVariables& work,
		       ComputeTime& com)
{
  TimeStart(START3);
  if (bMat_type == SPARSE){
    // set sparse_bMat zero
    sdpa_dset(sparse_bMat.NonZeroCount, 0.0, sparse_bMat.sp_ele, 1);
    //    compute_bMat_sparse_SDP(inputData,currentPt,work,com);
    compute_bMat_sparse_SDP_thread(inputData,currentPt,work,com);
    //   compute_bMat_sparse_SOCP(inputData,currentPt,work,com);
    compute_bMat_sparse_LP(inputData,currentPt,work,com);
  } else {
    bMat.setZero();
#if 1
    if ((Calc_F1 > 0) || (Calc_F2 > 0))
      compute_bMat_dense_SDP(inputData,currentPt,work,com);
#endif
    //    compute_bMat_dense_SDP123(inputData,currentPt,work,com);
    compute_bMat_dense_SDP3(inputData,currentPt,work,com);
    //compute_bMat_dense_SDP_thread(inputData,currentPt,work,com);
    //    compute_bMat_dense_SOCP(inputData,currentPt,work,com);
    compute_bMat_dense_LP(inputData,currentPt,work,com);
  }
  #if 0
  rMessage("bMat =  ");
  if (bMat_type == SPARSE) {
    bMat.display();
  }
  else {
    sparse_bMat.display();
  }
  #endif
  TimeEnd(END3);
  com.makebMat += TimeCal(START3,END3);
}

bool Newton::compute_DyVec(Newton::WHICH_DIRECTION direction,
			   SDPA_INT m,
			   InputData& inputData,
			   Chordal& chordal,
			   Solutions& currentPt,
			   WorkVariables& work,
			   ComputeTime& com,
			   FILE* Display, FILE* fpOut)
{
  if (direction == PREDICTOR) {
    TimeStart(START3_2);
    
    if (bMat_type == SPARSE){
      bool ret = chordal.factorizeSchur(m,diagonalIndex, Display, fpOut);
      if (ret == SDPA_FAILURE) {
	return SDPA_FAILURE;
      }
    } else {
      bool ret = Lal::choleskyFactorWithAdjust(bMat);
      if (ret == SDPA_FAILURE) {
	return SDPA_FAILURE;
      }
    }
    // rMessage("Cholesky of bMat =  ");
    // bMat.display();
    // sparse_bMat.display();
    TimeEnd(END3_2);
    com.choleskybMat += TimeCal(START3_2,END3_2);
  }
  // bMat is already cholesky factorized.


  TimeStart(START4);
  if (bMat_type == SPARSE){
    DyVec.copyFrom(gVec);
    chordal.solveSchur(DyVec);
  } else {
    Lal::let(DyVec,'=',bMat,'/',gVec);
  }
  TimeEnd(END4);
  com.solve += TimeCal(START4,END4);
  // rMessage("DyVec =  ");
  // DyVec.display();
  return SDPA_SUCCESS;
}

void Newton::compute_DzMat(InputData& inputData,
			   Residuals& currentRes,
			   Phase& phase,
			   ComputeTime& com)
{
  TimeStart(START_SUMDZ);
  inputData.multi_plusToA(DyVec, DzMat);
  Lal::let(DzMat,'=',DzMat,'*',&DMONE);
  if (phase.value == SolveInfo:: pFEAS
      || phase.value == SolveInfo::noINFO) {
    Lal::let(DzMat,'=',DzMat,'+',currentRes.dualMat);
  }
  TimeEnd(END_SUMDZ);
  com.sumDz += TimeCal(START_SUMDZ,END_SUMDZ);
}

void Newton::compute_DxMat(Solutions& currentPt,
			   WorkVariables& work,
			   ComputeTime& com)
{
  TimeStart(START_DX);
  // work.DLS1 = dX dZ Z^{-1}
  Jal::ns_jordan_triple_product(work.DLS1,currentPt.xMat,DzMat,
				currentPt.invzMat,work.DLS2);
  // dX = R Z^{-1} - dX dZ Z^{-1}
  Lal::let(DxMat,'=',r_zinvMat,'+',work.DLS1,&DMONE);
  TimeEnd(END_DX);
  TimeStart(START_SYMM);
  Lal::getSymmetrize(DxMat);
  TimeEnd(END_SYMM);
  // rMessage("DxMat =  ");
  // DxMat.display();
  com.makedX += TimeCal(START_DX,END_DX);
  com.symmetriseDx += TimeCal(START_SYMM,END_SYMM);
}


bool Newton::Mehrotra(Newton::WHICH_DIRECTION direction,
		      SDPA_INT m,
		      InputData& inputData,
		      Chordal& chordal,
		      Solutions& currentPt,
		      Residuals& currentRes,
		      AverageComplementarity& mu,
		      DirectionParameter& beta,
		      Switch& reduction,
		      Phase& phase,
		      WorkVariables& work,
		      ComputeTime& com,
		      FILE* Display, FILE* fpOut)
{
  //   rMessage("xMat, yVec, zMat =  ");
  //   currentPt.xMat.display();
  //   currentPt.yVec.display();
  //   currentPt.zMat.display();

  Make_gVec(direction, inputData, currentPt, currentRes,
	    mu, beta, phase, work, com);
  //  gVec.display();

  if (direction == PREDICTOR) {
    Make_bMat(inputData, currentPt, work, com);
    // bMat.display();
    // display_sparse_bMat();
    // display_index();
  }


  bool ret = compute_DyVec(direction,
			   m, inputData, chordal,
			   currentPt, work, com, Display, fpOut);
  if (ret == SDPA_FAILURE) {
    return SDPA_FAILURE;
  }
  //  rMessage("cholesky factorization =  ");
  //  sparse_bMat.display();

  TimeStart(START5);

  compute_DzMat(inputData, currentRes, phase, com);
  compute_DxMat(currentPt, work, com);

  TimeEnd(END5);
  com.makedXdZ += TimeCal(START5,END5);

  // rMessage("DxMat, DyVec, DzMat =  ");
  //   DxMat.display();
  //   DyVec.display();
  //   DzMat.display();

  return true;
}

void Newton::display(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }

  fprintf(fpout,"rNewton.DxMat = \n");
  DxMat.display(fpout);
  fprintf(fpout,"rNewton.DyVec = \n");
  DyVec.display(fpout);
  fprintf(fpout,"rNewton.DzMat = \n");
  DzMat.display(fpout);
}

void Newton::display_index(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }
  printf("display_index: %d %d %d\n",SDP_nBlock,SOCP_nBlock,LP_nBlock);

  for (SDPA_INT l=0; l<SDP_nBlock; l++){
    printf("SDP:%dth block\n",l);
    for (SDPA_INT k=0; k<SDP_number[l]; k++){
      printf("SDP(i=%d,ib=%d,j=%d,jb=%d) for target = %d\n",
	     SDP_constraint1[l][k],SDP_blockIndex1[l][k],
	     SDP_constraint2[l][k],SDP_blockIndex2[l][k], 
	     SDP_location_sparse_bMat[l][k]);
    }
  }

  for (SDPA_INT l=0; l<SOCP_nBlock; l++){
    printf("SOCP:%dth block\n",l);
    for (SDPA_INT k=0; k<SOCP_number[l]; k++){
      printf("SOCP(i=%d,ib=%d,j=%d,jb=%d) for target = %d\n",
	     SOCP_constraint1[l][k],SOCP_blockIndex1[l][k],
	     SOCP_constraint2[l][k],SOCP_blockIndex2[l][k], 
	     SOCP_location_sparse_bMat[l][k]);
    }
  }

  for (SDPA_INT l=0; l<LP_nBlock; l++){
    printf("LP:%dth block\n",l);
    for (SDPA_INT k=0; k<LP_number[l]; k++){
      printf("LP(i=%d,ib=%d; j=%d,jb=%d) for target = %d\n",
	     LP_constraint1[l][k],LP_blockIndex1[l][k],
	     LP_constraint2[l][k],LP_blockIndex2[l][k], 
	     LP_location_sparse_bMat[l][k]);
    }

  }

}

void Newton::display_sparse_bMat(FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }
  fprintf(fpout,"{\n");
  for (SDPA_INT index=0; index<sparse_bMat.NonZeroCount; ++index) {
    SDPA_INT i        = sparse_bMat.row_index[index];
    SDPA_INT j        = sparse_bMat.column_index[index];
    double value = sparse_bMat.sp_ele[index];
    fprintf(fpout,"val[%d,%d] = %e\n", i,j,value);
  }
  fprintf(fpout,"}\n");
}

} // end of namespace 'sdpa'

