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

#define DIMACS_PRLONG INT 0
#define DIMACS_ERROR_NUMBER 6

#include "sdpa_io.h"
#include "sdpa_linear.h"
#include "sdpa_jordan.h"
#include <vector>
#include <algorithm>

namespace sdpa {

void IO::read(FILE* fpData, FILE* fpout, SDPA_INT& m, char* str)
{
  while (true) {
    volatile SDPA_INT dummy=0; dummy++;//for gcc-3.3 bug
    fgets(str,lengthOfString,fpData);
    if (str[0]=='*' || str[0]=='"') {
      fprintf(fpout,"%s",str);
    } else {
      sscanf(str,"%ld",&m);
      break;
    }
  }
}

void IO::read(FILE* fpData, SDPA_INT & nBlock)
{
  fscanf(fpData,"%d",&nBlock);
}

void IO::read(FILE* fpData, BlockStruct& bs)
{
  for (SDPA_INT l=0; l<bs.nBlock; ++l) {
    fscanf(fpData,"%*[^0-9+-]%d",&bs.blockStruct[l]);
  }
  // only for SDP and LP
  for (SDPA_INT l=0; l<bs.nBlock; ++l) {
    if (bs.blockStruct[l] > 0 ) {
      bs.blockType[l] = BlockStruct::btSDP;
    }
    if (bs.blockStruct[l] < 0 ) {
      bs.blockType[l] = BlockStruct::btLP;
    }
  }
}

void IO::read(FILE* fpData, Vector& b)
{
  for (SDPA_INT k=0; k<b.nDim; ++k) {
    fscanf(fpData,"%*[^0-9+-]%lf",&b.ele[k]);
  }
}

void IO::read(FILE* fpData, DenseLinearSpace& xMat,
	      Vector& yVec, DenseLinearSpace& zMat,
	      BlockStruct& bs, bool inputSparse)
{
  // yVec is opposite sign
  SDPA_INT k=0;
  double tmp;
  if (fscanf(fpData,"%lf",&tmp) > 0) {
    // if y[0] locates the first charcter in fpData
    // then we need the following line
    yVec.ele[k] = -tmp;
    // rMessage("yVec.ele[" << k << "] = " << -tmp);
    k++;
  }
  for (; k<yVec.nDim; ++k) {
    fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
    yVec.ele[k] = -tmp;
    // rMessage("yVec.ele[" << k << "] = " << -tmp);
  }

  if (inputSparse) {
    // sparse case , zMat , xMat in this order
    SDPA_INT i,j,l,target;
    double value;
    while (true) {
      if (fscanf(fpData,"%*[^0-9+-]%d",&target)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&l)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&i)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&j)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%lf",&value)<=0) {
	break;
      }

      if (value == 0.0)
        continue;
      
#if 0
      rMessage("target = " << target
	       << ": l " << l
	       << ": i " << i
	       << ": j " << j
	       << ": value " <<value);
      #endif

      if (bs.blockType[l-1] == BlockStruct::btSDP) {
	SDPA_INT l2 = bs.blockNumber[l-1];
	if (target==1) {
	  zMat.setElement_SDP(l2,i-1,j-1,value);
	} else {
	  xMat.setElement_SDP(l2,i-1,j-1,value);
	}
      } else if (bs.blockType[l-1] == BlockStruct::btSOCP) {
	rError("io:: current version does not support SOCP");
	SDPA_INT l2 = bs.blockNumber[l-1];
	if (target==1) {
	  zMat.setElement_SOCP(l2,i-1,j-1,value);
	} else {
	  xMat.setElement_SOCP(l2,i-1,j-1,value);
	}
      } else if (bs.blockType[l-1] == BlockStruct::btLP) {
	if (i != j){
	  rError("io:: LP part  3rd element != 4th element\n"
		 "column should be the same as row in LP part.");
	}
	#if 0
	rMessage("l = " << l
		 << ": blockNumber[l-1] = " << bs.blockNumber[l-1]
		 << ": index = " << bs.blockNumber[l-1]+i-1
		 << ": i = " << i);
	#endif
	if (target==1) {
	  zMat.setElement_LP(bs.blockNumber[l-1]+i-1,value);
	} else {
	  xMat.setElement_LP(bs.blockNumber[l-1]+i-1,value);
	}
      }
    } // end of 'while (true)'
  } else {
    // dense case , zMat , xMat in this order
    // for SDP
    for (SDPA_INT l=0; l<bs.nBlock; ++l) {
      if (bs.blockType[l] == BlockStruct::btSDP) {
	SDPA_INT l2 =   bs.blockNumber[l];
	SDPA_INT size = bs.blockStruct[l];
	for (SDPA_INT i=0; i<size; ++i) {
	  for (SDPA_INT j=0; j<size; ++j) {
	    double tmp;
	    fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	    if (i<=j && tmp!=0.0) {
	      zMat.setElement_SDP(l2,i,j,tmp);
	    }
	  }
	}
      }
      else if (bs.blockType[l] == BlockStruct::btSOCP) {
	rError("io:: current version does not support SOCP");
      }
      else if (bs.blockType[l] == BlockStruct::btLP) {
	SDPA_INT size  = bs.blockStruct[l];
	SDPA_INT index = bs.blockNumber[l];
	for (SDPA_INT j=0; j<size; ++j) {
	  double tmp;
	  fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	  if (tmp!=0.0) {
	    zMat.setElement_LP(index,tmp);
	  }
	  index++;
	}
      }
    }
    
    for (SDPA_INT l=0; l<bs.nBlock; ++l) {
      if (bs.blockType[l] == BlockStruct::btSDP) {
	SDPA_INT l2   = bs.blockNumber[l];
	SDPA_INT size = bs.blockStruct[l];
	for (SDPA_INT i=0; i<size; ++i) {
	  for (SDPA_INT j=0; j<size; ++j) {
	    double tmp;
	    fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	    if (i<=j && tmp!=0.0) {
	      xMat.setElement_SDP(l2,i,j,tmp);
	    }
	  }
	}
      }
      else if (bs.blockType[l] == BlockStruct::btSOCP) {
	rError("io:: current version does not support SOCP");
      }
      else if (bs.blockType[l] == BlockStruct::btLP) {
	SDPA_INT size  = bs.blockStruct[l];
	SDPA_INT index = bs.blockNumber[l];
	for (SDPA_INT j=0; j<size; ++j) {
	  double tmp;
	  fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	  if (tmp!=0.0) {
	    xMat.setElement_LP(index,tmp);
	  }
	  index++;
	}
      }
    }
  } // end of 'if (inputSparse)'
}

// 2008/02/27 kazuhide nakata
// without LP_ANonZeroCount
#if 1
void IO::read(FILE* fpData, SDPA_INT m,
	      BlockStruct& bs,
	      InputData& inputData, bool isDataSparse)
{
  inputData.initialize_bVec(m);
  read(fpData,inputData.b);
  long position = ftell(fpData);

  // C,A must be accessed "double".

  //   initialize block struct of C and A
  setBlockStruct(fpData, inputData, m, bs,
                 position, isDataSparse);
  //   rMessage(" C and A initialize over");
    
  setElement(fpData, inputData, m, bs,
             position, isDataSparse);
  //   rMessage(" C and A have been read");
}
#endif

// 2008/02/27 kazuhide nakata   
// without LP_ANonZeroCount
void IO::setBlockStruct(FILE* fpData, InputData& inputData, SDPA_INT m,
			BlockStruct& bs,
                        long position, bool isDataSparse)
{
  // seed the positon of C in the fpData
  fseek(fpData, position, 0);

  vector<int>* SDP_index;
  NewArray(SDP_index,vector<int>,m+1);
  vector<int>* SOCP_index;
  NewArray(SOCP_index,vector<int>,m+1);
  vector<int>* LP_index;
  NewArray(LP_index,vector<int>,m+1);

  // for SDP
  SDPA_INT SDP_sp_nBlock;
  SDPA_INT* SDP_sp_index;
  SDPA_INT* SDP_sp_blockStruct;
  SDPA_INT* SDP_sp_NonZeroNumber;
  NewArray(SDP_sp_index,SDPA_INT,bs.SDP_nBlock);
  NewArray(SDP_sp_blockStruct,SDPA_INT,bs.SDP_nBlock);
  NewArray(SDP_sp_NonZeroNumber,SDPA_INT,bs.SDP_nBlock);
  // for SOCP
  SDPA_INT SOCP_sp_nBlock;
  SDPA_INT* SOCP_sp_blockStruct;
  SDPA_INT* SOCP_sp_index;
  SDPA_INT* SOCP_sp_NonZeroNumber;
  // for LP
  SDPA_INT LP_sp_nBlock;
  SDPA_INT* LP_sp_index;
  NewArray(LP_sp_index,SDPA_INT,bs.LP_nBlock);

  if (isDataSparse) {
    SDPA_INT i,j,k,l;
    double value;
    while (true) {
      if (fscanf(fpData,"%*[^0-9+-]%d",&k)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&l)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&i)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&j)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%lf",&value)<=0) {
	break;
      }
      
      if (bs.blockType[l-1] == BlockStruct::btSDP) {
        SDPA_INT l2 = bs.blockNumber[l-1];
        SDP_index[k].push_back(l2);
      } else if (bs.blockType[l-1] == BlockStruct::btSOCP) {
        rError("io:: current version does not support SOCP");
        SDPA_INT l2 = bs.blockNumber[l-1];;
        SOCP_index[k].push_back(l2);
      } else if (bs.blockType[l-1] == BlockStruct::btLP) {
        if (i!=j){
          printf("invalid data file k:%d, l:%d, i:%d, j:%d, value:%lf\n",
                 k,l,i,j,value);
          rError("IO::initializeLinearSpace");
        }
        SDPA_INT l2 = bs.blockNumber[l-1];
        LP_index[k].push_back(l2+i-1);
      } else {
        rError("io::read not valid blockType");
      }
    }// end of 'while (true)'
    
  } else { // isDataSparse == false
    
    // constant matrix
    for (SDPA_INT l=0; l<bs.nBlock; ++l){
      if (bs.blockType[l] == BlockStruct::btSDP) {
        SDPA_INT l2   = bs.blockNumber[l];
        SDPA_INT size = bs.SDP_blockStruct[l2];
        for (SDPA_INT i=0; i<size; ++i) {
          for (SDPA_INT j=0; j<size; ++j) {
            double tmp;
            fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
            if (i<=j && tmp!=0.0) {
              SDP_index[0].push_back(l2);
            }
          }
        }
      } else if (bs.blockType[l] == BlockStruct::btSOCP) {
        rError("io:: current version does not support SOCP");
      } else if (bs.blockType[l] == BlockStruct::btLP) { // LP part
	SDPA_INT start = bs.blockNumber[l];
	SDPA_INT size  = bs.blockStruct[l];
        for (SDPA_INT j=0; j<size; ++j) {
          double tmp;
          fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
          if (tmp!=0.0) {
              LP_index[0].push_back(start+j);
          }
        }
      } else {
        rError("io::read not valid blockType");
      }
    }
    // data matrices
    for (SDPA_INT k=0; k<m; ++k) {
      for (SDPA_INT l=0; l<bs.nBlock; ++l){
        if (bs.blockType[l] == BlockStruct::btSDP) {
          SDPA_INT l2 =   bs.blockNumber[l];
          SDPA_INT size = bs.SDP_blockStruct[l2];
          for (SDPA_INT i=0; i<size; ++i) {
            for (SDPA_INT j=0; j<size; ++j) {
              double tmp;
              fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
              if (i<=j && tmp!=0.0) {
                SDP_index[k+1].push_back(l2);
              }
            }
          }
        } else if (bs.blockType[l] == BlockStruct::btSOCP) {
          rError("io:: current version does not support SOCP");
        } else if (bs.blockType[l] == BlockStruct::btLP) {
	  SDPA_INT start = bs.blockNumber[l];
	  SDPA_INT size  = bs.blockStruct[l];
          for (SDPA_INT j=0; j<size; ++j) {
            double tmp;
            fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
            if (tmp!=0.0) {
              LP_index[k+1].push_back(start+j);
            }
          }
        } else {
          rError("io::read not valid blockType");
        }
      }
    }
    
  } // end of 'if (isDataSparse)'

  NewArray(inputData.A,SparseLinearSpace,m);
  for (SDPA_INT k=0 ; k<m+1; k++){
    sort(SDP_index[k].begin(),SDP_index[k].end());
    SDP_sp_nBlock = 0;
    SDPA_INT previous_index = -1;
    SDPA_INT index;
    for (SDPA_INT i=0; i<SDP_index[k].size(); i++){
      index = SDP_index[k][i];
      if (previous_index != index){
        SDP_sp_index[SDP_sp_nBlock] = index;
        SDP_sp_blockStruct[SDP_sp_nBlock] = bs.SDP_blockStruct[index];
        SDP_sp_NonZeroNumber[SDP_sp_nBlock] = 1;
        previous_index = index;
        SDP_sp_nBlock++;
      } else {
        SDP_sp_NonZeroNumber[SDP_sp_nBlock-1]++;
      }
    }

    // dummy initialization to surpress compiler warning
    SOCP_sp_nBlock        = 0;
    SOCP_sp_blockStruct   = NULL;
    SOCP_sp_index         = NULL;
    SOCP_sp_NonZeroNumber = NULL;
    
    sort(LP_index[k].begin(),LP_index[k].end());
    LP_sp_nBlock=0;
    previous_index = -1;
    for (SDPA_INT i=0; i<LP_index[k].size(); i++){
      index = LP_index[k][i];
      if (previous_index != index){
        LP_sp_index[LP_sp_nBlock] = index;
        previous_index = index;
        LP_sp_nBlock++;
      }
    }

    if (k==0){
      inputData.C.initialize(SDP_sp_nBlock, 
                             SDP_sp_index,
                             SDP_sp_blockStruct, 
                             SDP_sp_NonZeroNumber,
                             SOCP_sp_nBlock, 
                             SOCP_sp_blockStruct, 
                             SOCP_sp_index,
                             SOCP_sp_NonZeroNumber,
                             LP_sp_nBlock, 
                             LP_sp_index);
    } else {
      inputData.A[k-1].initialize(SDP_sp_nBlock, 
                                  SDP_sp_index,
                                  SDP_sp_blockStruct, 
                                  SDP_sp_NonZeroNumber,
                                  SOCP_sp_nBlock, 
                                  SOCP_sp_blockStruct, 
                                  SOCP_sp_index,
                                  SOCP_sp_NonZeroNumber,
                                  LP_sp_nBlock, 
                                  LP_sp_index);
    }
  }

  DeleteArray(SDP_index);
  DeleteArray(SOCP_index);
  DeleteArray(LP_index);

  DeleteArray(SDP_sp_index);
  DeleteArray(SDP_sp_blockStruct);
  DeleteArray(SDP_sp_NonZeroNumber);
  DeleteArray(SDP_sp_NonZeroNumber);
#if 0
  DeleteArray(SOCP_sp_index);
  DeleteArray(SOCP_sp_blockStruct);
  DeleteArray(SOCP_sp_NonZeroNumber);
#endif
  DeleteArray(LP_sp_index);
}


// 2008/02/27 kazuhide nakata   
// without LP_ANonZeroCount
void IO::setElement(FILE* fpData, InputData& inputData, SDPA_INT m,
		    BlockStruct& bs,
                    long position, bool isDataSparse)
{
  // in Sparse, read C,A[k]

  // seed the positon of C in the fpData
  fseek(fpData, position, 0);

  if (isDataSparse) {
    SDPA_INT i,j,k,l;
    double value;
    while (true) {
      if (fscanf(fpData,"%*[^0-9+-]%d",&k)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&l)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&i)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%d",&j)<=0) {
	break;
      }
      if (fscanf(fpData,"%*[^0-9+-]%lf",&value)<=0) {
	break;
      }
#if 0
      rMessage("input k:" << k <<
	       " l:" << l <<
	       " i:" << i <<
	       " j:" << j);
#endif     

      if (bs.blockType[l-1] == BlockStruct::btSDP) {
	SDPA_INT l2 = bs.blockNumber[l-1];
	if (k==0) {
	  inputData.C.setElement_SDP(l2,i-1,j-1,-value);
	} else {
	  inputData.A[k-1].setElement_SDP(l2,i-1,j-1,value);
	}
      } else if (bs.blockType[l-1] == BlockStruct::btSOCP) {
	rError("io:: current version does not support SOCP");
	SDPA_INT l2 = bs.blockNumber[l-1];
	if (k==0) {
	  inputData.C.setElement_SOCP(l2,i-1,j-1,-value);
	} else {
	  inputData.A[k-1].setElement_SOCP(l2,i-1,j-1,value);
	}
      } else if (bs.blockType[l-1] == BlockStruct::btLP) {
	if (i != j){
	  rError("io:: LP part  3rd element != 4th element\n"
		 "column should be same as row in LP part.");
	}
	if (k==0) {
	  inputData.C.setElement_LP(bs.blockNumber[l-1]+i-1,-value);
	} else {
	  inputData.A[k-1].setElement_LP(bs.blockNumber[l-1]+i-1,value);
	}
      } else {
	rError("io::read not valid blockType");
      }
    } 
  } else {  // dense

    // constant matrix
    for (SDPA_INT l=0; l<bs.nBlock; ++l){
      if (bs.blockType[l] == BlockStruct::btSDP) {
	SDPA_INT l2   = bs.blockNumber[l];
	SDPA_INT size = bs.SDP_blockStruct[l2];
	for (SDPA_INT i=0; i<size; ++i) {
	  for (SDPA_INT j=0; j<size; ++j) {
	    double tmp;
	    fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	    if (i<=j && tmp!=0.0) {
	      inputData.C.setElement_SDP(l2,i,j,-tmp);
	    }
	  }
	}
      } else if (bs.blockType[l] == BlockStruct::btSOCP) {
	rError("io:: current version does not support SOCP");
      } else if (bs.blockType[l] == BlockStruct::btLP) {
	SDPA_INT start = bs.blockNumber[l];
	SDPA_INT size  = bs.blockStruct[l];
	for (SDPA_INT j=0; j<size; ++j) {
	  double tmp;
	  fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	  if (tmp!=0.0) {
	    inputData.C.setElement_LP(start+j,-tmp);
	  }
	}
      } else {
	rError("io::read not valid blockType");
      }
    }

    // data matrices
    for (SDPA_INT k=0; k<m; ++k) {
      
      for (SDPA_INT l=0; l<bs.nBlock; ++l){
	if (bs.blockType[l] == BlockStruct::btSDP) {
	  SDPA_INT l2   = bs.blockNumber[l];
	  SDPA_INT size = bs.SDP_blockStruct[l2];
	  for (SDPA_INT i=0; i<size; ++i) {
	    for (SDPA_INT j=0; j<size; ++j) {
	      double tmp;
	      fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	      if (i<=j && tmp!=0.0) {
		inputData.A[k].setElement_SDP(l2,i,j,tmp);
	      }
	    }
	  }
	} else if (bs.blockType[l] == BlockStruct::btSOCP) {
	  rError("io:: current version does not support SOCP");
	} else if (bs.blockType[l] == BlockStruct::btLP) {
	  SDPA_INT start = bs.blockNumber[l];
	  SDPA_INT size  = bs.blockStruct[l];
	  for (SDPA_INT j=0; j<size; ++j) {
	    double tmp;
	    fscanf(fpData,"%*[^0-9+-]%lf",&tmp);
	    if (tmp!=0.0) {
	      inputData.A[k].setElement_LP(start+j,tmp);
	    }
	  }
	} else {
	  rError("io::read not valid blockType");
	}
      }
    } // for k

  } // end of 'if (isDataSparse)'

}

void IO::printHeader(FILE* fpout, FILE* Display)
{
  if (fpout) {
    fprintf(fpout,"   mu      thetaP  thetaD  objP      objD "
	    "     alphaP  alphaD  beta \n");
    fflush(fpout);
  }
  if (Display) {
    fprintf(Display,"   mu      thetaP  thetaD  objP      objD "
	    "     alphaP  alphaD  beta \n");
    fflush(Display);
  }
}

void IO::printOneIteration(SDPA_INT pIteration,
			    AverageComplementarity& mu,
			    RatioInitResCurrentRes& theta,
			    SolveInfo& solveInfo,
			    StepLength& alpha,
			    DirectionParameter& beta,
			    FILE* fpout,
			    FILE* Display)
{
  FILE* fp = NULL;
  for (SDPA_INT fp_index=0; fp_index<2; ++fp_index) {
    if (fp_index == 0) {
      fp = fpout;
    }
    else {
      fp = Display;
    }
    if (fp == NULL) {
      continue;
    }
  
    #if REVERSE_PRIMAL_DUAL
    fprintf(fp,"%2d %4.1e %4.1e %4.1e %+7.2e %+7.2e"
	    " %4.1e %4.1e %4.2e\n", pIteration, mu.current,
	    theta.dual, theta.primal,
	    -solveInfo.objValDual,-solveInfo.objValPrimal,
	    alpha.dual, alpha.primal, beta.value);

    #else
    fprintf(fp,"%2d %4.1e %4.1e %4.1e %+7.2e %+7.2e"
	    " %4.1e %4.1e %4.2e\n", pIteration, mu.current,
	    theta.primal, theta.dual,
	    solveInfo.objValPrimal, solveInfo.objValDual,
	    alpha.primal, alpha.dual, beta.value);
    #endif
  }
}


void IO::printLastInfo(SDPA_INT pIteration,
		       AverageComplementarity& mu,
		       RatioInitResCurrentRes& theta,
		       SolveInfo& solveInfo,
		       StepLength& alpha,
		       DirectionParameter& beta,
		       Residuals& currentRes,
		       Phase & phase,
		       Solutions& currentPt,
		       InputData& inputData,
                       WorkVariables& work,
		       double cputime,
		       ComputeTime& com,
		       Parameter& param,
		       FILE* fpout,
		       FILE* Display,
		       bool printTime)
{
  // SDPA_INT nDim = currentPt.nDim;

  printOneIteration(pIteration,mu,theta,solveInfo,alpha,
		    beta, fpout, Display);

  double mean = (fabs(solveInfo.objValPrimal)
		 + fabs(solveInfo.objValDual)) / 2.0;
  double PDgap = fabs(solveInfo.objValPrimal
		      - solveInfo.objValDual);
  // double dominator;
  double relgap;
  if (mean < 1.0) {
    relgap = PDgap;
  } else {
    relgap = PDgap/mean;
  }

  // double gap    = mu.current*nDim;
  double gap = solveInfo.objValPrimal - solveInfo.objValDual;
  double digits = 1000; // 1000 means infinity in this case
  digits = -log10(fabs(PDgap/mean));

  FILE* fp = NULL;
  for (SDPA_INT fp_index = 0; fp_index < 2; fp_index++) {
    if (fp_index == 0) {
      fp = Display;
    }
    else {
      fp = fpout;
    }
    if (fp == NULL) {
      continue;
    }
    fprintf(fp, "\n");
    phase.display(fp);
    fprintf(fp, "   Iteration = %d\n",  pIteration);
    fprintf(fp, "          mu = ");
    fprintf(fp, param.infPrint, mu.current);
    fprintf(fp, "\n");
    fprintf(fp, "relative gap = ");
    fprintf(fp, param.infPrint, relgap);
    fprintf(fp, "\n");
    fprintf(fp, "        gap  = ");
    fprintf(fp, param.infPrint, gap);
    fprintf(fp, "\n");
    fprintf(fp, "     digits  = ");
    fprintf(fp, param.infPrint, digits);
    fprintf(fp, "\n");
#if REVERSE_PRIMAL_DUAL
    fprintf(fp, "objValPrimal = ");
    fprintf(fp, param.infPrint, -solveInfo.objValDual);
    fprintf(fp, "\n");
    fprintf(fp, "objValDual   = ");
    fprintf(fp, param.infPrint, -solveInfo.objValPrimal);
    fprintf(fp, "\n");
    fprintf(fp, "p.feas.error = ");
    fprintf(fp, param.infPrint, currentRes.normDualMat);
    fprintf(fp, "\n");
    fprintf(fp, "d.feas.error = ");
    fprintf(fp, param.infPrint, currentRes.normPrimalVec);
    fprintf(fp, "\n");
#else
    fprintf(fp, "objValPrimal = ");
    fprintf(fp, param.infPrint, solveInfo.objValPrimal);
    fprintf(fp, "\n");
    fprintf(fp, "objValDual   = ");
    fprintf(fp, param.infPrint, solveInfo.objValDual);
    fprintf(fp, "\n");
    fprintf(fp, "p.feas.error = ");
    fprintf(fp, param.infPrint, currentRes.normPrimalVec);
    fprintf(fp, "\n");
    fprintf(fp, "d.feas.error = ");
    fprintf(fp, param.infPrint, currentRes.normDualMat);
    fprintf(fp, "\n");
#endif
    if (printTime == true) {
      fprintf(fp, "total time   = %.6f\n",cputime);
    }

  }

  if (fpout) {
    param.display(fpout,param.infPrint);
    com.display(fpout);
  }

  #if DIMACS_PRINT
  double dimacs_error[DIMACS_ERROR_NUMBER+1];
  computeDimacs(dimacs_error,solveInfo, currentRes, currentPt,
		inputData, work);
  printDimacs(dimacs_error, param.xPrint, fpout);
  printDimacs(dimacs_error, param.xPrint, Display);
  #endif

}

void IO::computeDimacs(double* dimacs_error,
		       SolveInfo& solveInfo,
		       Residuals& currentRes,
		       Solutions& currentPt,
		       InputData& inputData,
		       WorkVariables& work)
{
  double b1     = Lal::getOneNorm(inputData.b);
  double c1     = Lal::getOneNorm(inputData.C);
  double p_norm = sqrt(Lal::getTwoNorm(currentRes.primalVec));
  double d_norm = sqrt(Lal::getTwoNorm(currentRes.dualMat));
  double x_min  = Jal::getMinEigen(currentPt.xMat,work);
  double z_min  = Jal::getMinEigen(currentPt.zMat,work);

  #if 0
  printf("b1:%e\n",b1);
  printf("c1:%e\n",c1);
  printf("p_norm:%e\n",p_norm);
  printf("d_norm:%e\n",d_norm);
  printf("x_min:%e\n",x_min);
  printf("z_min:%e\n",z_min);
  #endif
  
  double ctx = solveInfo.objValPrimal;
  double bty = solveInfo.objValDual;
  double xtz = 0.0;
  Lal::let(xtz,'=',currentPt.xMat,'.',currentPt.zMat);

  for (SDPA_INT i=0; i<=6; ++i) {
    dimacs_error[i] = 0.0;
  }

  dimacs_error[1] = p_norm / (1+b1);
  dimacs_error[2] = max( 0.0, - x_min / (1+b1));
  dimacs_error[3] = d_norm / (1+c1);
  dimacs_error[4] = max( 0.0, - z_min / (1+c1));
  dimacs_error[5] = (ctx - bty) / (1 + fabs(ctx) + fabs(bty));
  dimacs_error[6] = xtz / (1 + fabs(ctx) + fabs(bty));
}

void IO::printDimacs(double* DimacsError,char* printFormat,
		     FILE* fpout)
{
  if (fpout == NULL) {
    return;
  }
  fprintf(fpout,  "\n");
  fprintf(fpout,  "* DIMACS_ERRORS * \n");
  fprintf(fpout,  "err1 = ");
  fprintf(fpout,  printFormat, DimacsError[1]);
  fprintf(fpout, "  [||Ax-b|| / (1+||b||_1)]\n");
  fprintf(fpout,  "err2 = ");
  fprintf(fpout,  printFormat, DimacsError[2]);
  fprintf(fpout, "  [max(0, -lambda(x)/(1+||b||_1))]\n");
  fprintf(fpout,  "err3 = ");
  fprintf(fpout,  printFormat, DimacsError[3]);
  fprintf(fpout, "  [||A^Ty + z - c || / (1+||c||_1)]\n");
  fprintf(fpout,  "err4 = ");
  fprintf(fpout,  printFormat, DimacsError[4]);
  fprintf(fpout, "  [max(0, -lambda(z)/(1+||c||_1))]\n");
  fprintf(fpout,  "err5 = ");
  fprintf(fpout,  printFormat, DimacsError[5]);
  fprintf(fpout, "  [(<c,x> - <b,y>) / (1 + |<c,x>| + |<b,y>|)]\n");
  fprintf(fpout,  "err6 = ");
  fprintf(fpout,  printFormat, DimacsError[6]);
  fprintf(fpout, "  [<x,z> / (1 + |<c,x>| + |<b,y>|)]\n");
  fprintf(fpout,  "\n");
}

void IO::printSolution(BlockStruct& bs, Solutions& currentPt,
		       Parameter& param, FILE* fpout)
{
  if (fpout != NULL) {
    #if 1
    #if REVERSE_PRIMAL_DUAL
    fprintf(fpout,"xVec = \n");
    currentPt.yVec.display(fpout,1.0,param.xPrint);
    fprintf(fpout,"xMat = \n");
    currentPt.zMat.displaySolution(bs,fpout,param.XPrint);
    fprintf(fpout,"yMat = \n");
    currentPt.xMat.displaySolution(bs,fpout,param.YPrint);
    #else
    fprintf(fpout,"xMat = \n");
    currentPt.xMat.displaySolution(bs,fpout,param.XPrint);
    fprintf(fpout,"yVec = \n");
    currentPt.yVec.display(fpout,1.0,param.xPrint);
    fprintf(fpout,"zMat = \n");
    currentPt.zMat.displaySolution(bs,fpout,param.YPrint);
    #endif
    #endif
  }
}

} // end of namespace 'sdpa'
