#include "mex.h"
#include "math.h"
#include "DPP.h"



void SamplePath(double *xSampled, double *xPath, double *tPath, double *grid,
        mwSize numStates, mwSize numEvents, mwSize gridLength)
{
  
    int i, counter, tPathIdx;
    double currTime;
    double currX[numStates];
    //currX[0] = 10;
    
    tPathIdx = 0;
    
    for (i=0; i<gridLength; i++)
    {
        currTime = grid[i];
        //mexPrintf("currTime=%f tPathTime=%f\n", currTime, tPath[tPathIdx]);
        
        while (currTime >= tPath[tPathIdx])
        {
            tPathIdx++;
            
            if (tPathIdx > (numEvents - 1))
            {       
                break;
            }
        }
        
        tPathIdx--;
        
        GetValue(currX, xPath, tPathIdx, numStates);
        SetValue(xSampled, currX, i, numStates);
        
        //mexPrintf("gridIdx=%f\n", currX[0]);
    }  
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  double *xPath, *tPath, *grid, *xSampled;
  int numStates, numEvents, gridLength;

 

  
  /* The input must be a noncomplex scalar double.*/
  numEvents = mxGetN(prhs[0]);
  numStates = mxGetM(prhs[0]);
  gridLength = mxGetN(prhs[2]);
  
  /* Create matrix for the return argument. */
  plhs[0] = mxCreateDoubleMatrix(numStates,gridLength, mxREAL);
  xSampled = mxGetPr(plhs[0]);
   
  /* Assign pointers to each input and output. */
  xPath = mxGetPr(prhs[0]);
  tPath = mxGetPr(prhs[1]);
  grid = mxGetPr(prhs[2]);
  
  /* Call the subroutine. */
  SamplePath(xSampled, xPath, tPath, grid, numStates, numEvents, gridLength);
  
}
