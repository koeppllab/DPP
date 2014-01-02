#include "mex.h"
#include "math.h"



int CountmRNA(double *rPath, double *deltamRNA, 
        int offIdx, int synthIdx, mwSize numEvents)
{
  
    int i, counter;
    double geneVal;
    int rValue, mRNACount;
    
    counter = 0;
    mRNACount = 0;
    
    //mexPrintf("Idx: %d %d\n", geneIdx, mRNAIdx);
    
    for (i=0; i<numEvents; i++)
    {
        rValue = (int)rPath[i];
        
        if (rValue == synthIdx)
        {
            mRNACount = mRNACount + 1;
        }
        else if (rValue == offIdx)
        {
            //mexPrintf("Gene turned off (mRNA produced: %d)\n", mRNACount);
            deltamRNA[counter] = (double)mRNACount;
            mRNACount = 0;
            counter++;
        }
        //mexPrintf("r: %d  %d\n", mRNACount, rValue);

    }
    
    return counter;
  
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  double *rPath, *xPath, *deltamRNA;
  int offIdx, synthIdx, numGeneOns, numEvents;

 

  
  /* The input must be a noncomplex scalar double.*/
  numEvents = mxGetN(prhs[0]);
  
  /* Create matrix for the return argument. */
  plhs[0] = mxCreateDoubleMatrix(1,numEvents, mxREAL);
 
   
  /* Assign pointers to each input and output. */
  rPath = mxGetPr(prhs[0]);
  offIdx = mxGetScalar(prhs[1]);
  synthIdx = mxGetScalar(prhs[2]);

  deltamRNA = mxGetPr(plhs[0]);
  
  /* Call the timestwo subroutine. */
  numGeneOns = CountmRNA(rPath, deltamRNA, offIdx, synthIdx, numEvents);
  mxSetN(plhs[0], numGeneOns);
  
}
