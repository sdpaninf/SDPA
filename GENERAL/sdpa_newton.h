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

#ifndef __sdpa_newton_h__
#define __sdpa_newton_h__

#include "sdpa_chordal.h"

#define SparseCholesky 1

#ifdef GOTO_BLAS
extern "C" {
  void goto_set_num_threads(int);
};
#endif

namespace sdpa {

class Newton;

class Solutions;
class InputData;
class Residuals;
class WorkVariables;

class ComputeTime;
class Parameter;
class StepLength;
class DirectionParameter;
class Switch;
class RatioInitResCurrentRes;
class SolveInfo;
class Phase;
class AverageComplementarity;



class Newton
{
public:
  enum bMat_Sp_De {SPARSE, DENSE};
  bMat_Sp_De bMat_type;

  SparseMatrix sparse_bMat;
  DenseMatrix bMat; // the coefficent of Schur complement
  Vector      gVec; // the right hand side of Schur complement
  DenseLinearSpace DxMat;
  Vector           DyVec;
  DenseLinearSpace DzMat;

  DenseLinearSpace r_zinvMat;
  DenseLinearSpace x_rd_zinvMat;
  
  enum FormulaType {F1,F2,F3};
  FormulaType** useFormula;
  double**      FormulaEstimate;
  SDPA_INT**    useThread;

  // Caution: 
  // if SDPA doesn't use sparse bMat, following variables are indefinite.
  //
  // nBLock : number of block
  // nConstraint[k]: number of combination of nonzero matrices in k-th block
  // when A[k].block[i] and A[k].block[j] are nonzero matrices, 
  //     i             <-> constraint1[k][t]
  //     j             <-> constraint2[k][t]
  //     A[k].block[i] <-> A[k].sp_block[blockIndex1[k][t]]
  //     A[k].block[j] <-> A[k].sp_block[blockIndex2[k][t]]
  //     B_{ij}        <-> sparse_bMat.sp_ele[location_sparse_bMat[k][t]]
  SDPA_INT SDP_nBlock;  SDPA_INT* SDP_number;  
  SDPA_INT** SDP_constraint1;  SDPA_INT** SDP_constraint2;
  SDPA_INT** SDP_blockIndex1;  SDPA_INT** SDP_blockIndex2;
  SDPA_INT** SDP_location_sparse_bMat;
  SDPA_INT SOCP_nBlock;  SDPA_INT* SOCP_number;  
  SDPA_INT** SOCP_constraint1;  SDPA_INT** SOCP_constraint2;
  SDPA_INT** SOCP_blockIndex1;  SDPA_INT** SOCP_blockIndex2;
  SDPA_INT** SOCP_location_sparse_bMat;
  SDPA_INT LP_nBlock;  SDPA_INT* LP_number;  
  SDPA_INT** LP_constraint1;  SDPA_INT** LP_constraint2;
  SDPA_INT** LP_blockIndex1;  SDPA_INT** LP_blockIndex2;
  SDPA_INT** LP_location_sparse_bMat;

  // from index of aggrigate sparsity pattern to index of sparse_bMat
  // B_{ii} <-> sparse_bMat[diagonalIndex[i]]
  SDPA_INT* diagonalIndex;
  // B_{ij} for all i is between diagonalIndex[j] and rowStartIndex[j+1]

  Newton();
  Newton(SDPA_INT m, BlockStruct& bs);
  ~Newton();
  
  void initialize(SDPA_INT m, BlockStruct& bs);

  void terminate();

  void initialize_dense_bMat(SDPA_INT m);
  // 2008/03/12 kazuhide nakata
  void initialize_sparse_bMat(SDPA_INT m);
  // 2008/03/12 kazuhide nakata
  void initialize_bMat(SDPA_INT m, Chordal& chordal, InputData& inputData, 
                       FILE* Display, FILE* fpOut);

  SDPA_INT binarySearchIndex(SDPA_INT i, SDPA_INT j);
  void make_aggrigateIndex_SDP(InputData& inputData);
  void make_aggrigateIndex_SOCP(InputData& inputData);
  void make_aggrigateIndex_LP(InputData& inputData);
  void make_aggrigateIndex(InputData& inputData);

  void computeFormula_SDP(InputData& inputData,
			  double DenseRatio,double Kappa);

  void computeFormula_SDP_new(InputData& inputData);

  enum WHICH_DIRECTION {PREDICTOR, CORRECTOR};
  void compute_rMat(WHICH_DIRECTION direction,
		    AverageComplementarity& mu,
		    DirectionParameter& beta,
		    Solutions& cuurentPt,
		    WorkVariables& work);

  void Make_gVec(Newton::WHICH_DIRECTION direction,
	       InputData& inputData,
	       Solutions& currentPt,
	       Residuals& currentRes,
	       AverageComplementarity& mu,
	       DirectionParameter& beta,
	       Phase& phase,
	       WorkVariables& work,
	       ComputeTime& com);

  void calF1(double& ret, DenseMatrix& G,
	     SparseMatrix& Ai);
  void calF2(double& ret, DenseMatrix& F, DenseMatrix& G,
	     DenseMatrix& invZ, SparseMatrix& Ai, bool& hasF2Gcal);
  void calF3(double& ret,
	     DenseMatrix& X, DenseMatrix& invZ,
	     SparseMatrix& Ai, SparseMatrix& Aj);

  // B_{i,j} = (X A_i Z^{-1}) \bullet A_j
  void compute_bMat_dense_SDP(InputData& inputData,
			      Solutions& currentPt,
			      WorkVariables& work,
			      ComputeTime& com);

  void compute_bMat_dense_SDP123(InputData& inputData,
			      Solutions& currentPt,
			      WorkVariables& work,
			      ComputeTime& com);

  void compute_bMat_dense_SDP3(InputData& inputData,
			      Solutions& currentPt,
			      WorkVariables& work,
			      ComputeTime& com);

  void compute_bMat_dense_SDP_thread(InputData& inputData,
			      Solutions& currentPt,
			      WorkVariables& work,
			      ComputeTime& com);

  // B_{i,j} = (X A_i Z^{-1}) \bullet A_j 
  void compute_bMat_sparse_SDP(InputData& inputData,
			       Solutions& currentPt,
			       WorkVariables& work,
			       ComputeTime& com);

  void compute_bMat_sparse_SDP_thread(InputData& inputData,
			       Solutions& currentPt,
			       WorkVariables& work,
			       ComputeTime& com);

  void compute_bMat_dense_SOCP(InputData& inputData,
			       Solutions& currentPt,
			       WorkVariables& work,
			       ComputeTime& com);

  void compute_bMat_sparse_SOCP(InputData& inputData,
			       Solutions& currentPt,
			       WorkVariables& work,
			       ComputeTime& com);

  void compute_bMat_dense_LP(InputData& inputData,
			     Solutions& currentPt,
			     WorkVariables& work,
			     ComputeTime& com);

  void compute_bMat_sparse_LP(InputData& inputData,
			      Solutions& currentPt,
			      WorkVariables& work,
			      ComputeTime& com);

  void Make_bMat(InputData& inputData,
		 Solutions& currentPt,
		 WorkVariables& work,
		 ComputeTime& com);

  bool compute_DyVec(Newton::WHICH_DIRECTION direction,
		     SDPA_INT m,
		     InputData& inputData,
		     Chordal& chordal,
		     Solutions& currentPt,
		     WorkVariables& work,
		     ComputeTime& com,
		     FILE* Display, FILE* fpOut);

  void compute_DzMat(InputData& inputData,
		     Residuals& currentRes,
		     Phase& phase,
		     ComputeTime& com);
  
  void compute_DxMat(Solutions& currentPt,
		     WorkVariables& work,
		     ComputeTime& com);
  
  
  bool Mehrotra(WHICH_DIRECTION direction,
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
		FILE* Display, FILE* fpOut);
  
  void display(FILE* fpout=stdout);
  void display_index(FILE* fpout=stdout);
  void display_sparse_bMat(FILE* fpout=stdout);

  static pthread_mutex_t job_mutex;
  static pthread_cond_t  job_cond;
  static SDPA_INT  Column_Number;
  static bool mutex_flag;
  static SDPA_INT  Calc_F1;
  static SDPA_INT  Calc_F2;

  static void calF1_thread(double& ret, DenseMatrix& G,
			   SparseMatrix& Aj);
  static void calF2_thread(double& ret, DenseMatrix& F, DenseMatrix& G,
			   DenseMatrix& X, SparseMatrix& Aj,
			   bool& hasF2Gcal);
  static void calF3_thread(double& ret,
			   DenseMatrix& X, DenseMatrix& invZ,
			   SparseMatrix& Ai, SparseMatrix& Aj);
  static void calF3_thread_1x1(double& ret,
			       DenseMatrix& X, DenseMatrix& invZ,
			       SparseMatrix& Ai, SparseMatrix& Aj);
  static void calF3_thread_2(double& ret,
			     DenseMatrix& X, DenseMatrix& invZ,
			     SparseMatrix& Ai, SparseMatrix& Aj);

  static  void* compute_bMat_dense_SDP_thread_func(void *arg);
  static  void* compute_bMat_dense_SDP_thread_func_new(void *arg);
  static  void* compute_bMat_sparse_SDP_thread_func(void *arg);
  static  void* compute_bMat_sparse_SDP_thread_func_new(void *arg);

  SDPA_INT NUM_THREADS;
  SDPA_INT NUM_GOTOBLAS;
  void setNumThreads(FILE* Display, FILE* fpOut, SDPA_INT NumThreads=0);
};

 typedef struct _thread_arg {
   SDPA_INT Block_Number;
   SDPA_INT thread_num;
   SDPA_INT Num_of_Threads;
   SDPA_INT mDIM;
   SDPA_INT SDP_nBlock;
   SDPA_INT *SDP_number;
   SDPA_INT **SDP_constraint1;
   SDPA_INT **SDP_constraint2;
   SDPA_INT **SDP_blockIndex1;
   SDPA_INT **SDP_blockIndex2;
   SDPA_INT **SDP_location_sparse_bMat;
   DenseMatrix* bMat;
   SparseMatrix* sparse_bMat;
   Newton::FormulaType** useFormula;
   SDPA_INT **useThread;
   double** FormulaEstimate;
   InputData* inputData;
   Solutions* currentPt;
   WorkVariables* work;
   ComputeTime* com;
 } thread_arg_t;

} // end of namespace 'sdpa'

#endif // __sdpa_newton_h__
