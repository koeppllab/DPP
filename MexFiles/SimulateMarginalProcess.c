#include "mex.h"
#include "math.h"
#include "DPP.h"



void ExecuteReaction(double *x, double *pre, double *post, int reactionIndex, int numReactions, int numSpecies)
{
    
    int i, idx;
    double change;
    
    for (i=0; i < numSpecies; i++)
    {
        idx = numReactions * i + reactionIndex;
        change = post[idx] - pre[idx];
        
        x[i] = x[i] + change;
    }

}



int SimulatePath(double *x, double *t, double *react, double *G, 
         double *x0, double *c, double *pre, double *post, int n, 
         double maxTime, int tvRateIndex, double *inputTime, double *inputLevels,
         double *syncTime, double startTime, double *populationRateIdx, 
         double *a, double *b, double *initialA, double *initialB, 
         mwSize numReactions,
         mwSize numSpecies,
         mwSize numPopulationRates,
         mwSize numInputLevels,
         mwSize numSyncTimes)
{
    
    int i, k, j, arrayIdx, nextReaction,
            tvRateHeterogeneous, syncCounter, effNumSyncTimes;
    double h[numReactions], g[numReactions], currentX[numReactions], 
           alpha, beta, IntG[numPopulationRates], intG,
           R[numPopulationRates], r;
    double tau, minTau, currentTime, t0, inputValue;
    
    double cTmp[numReactions];
    
    int sync;
    
    SetValue(cTmp, c, 0, numReactions);

    SetValue(x, x0, 0, numSpecies);
    SetValue(currentX, x0, 0, numSpecies);
    
    t0 = startTime;
    SetValue(t, &t0, 0, 1);
    currentTime = startTime;
    
    SetValue(R, initialA, 0, numPopulationRates);
    SetValue(IntG, initialB, 0, numPopulationRates);
    
    //mexPrintf("%d\n", numPopulationRates);
    
//     for (k=0; k<numPopulationRates; k++)
//     {
//         mexPrintf("%f \n", R[k]);
//     }
//     mexPrintf("\n");

    
    
    cTmp[tvRateIndex - 1] = GetPiecewiseConstantInput(inputTime, inputLevels,
            currentTime, numInputLevels);
    
    tvRateHeterogeneous = GetIndex(populationRateIdx, numPopulationRates, tvRateIndex);
    
    syncCounter = 0;
    
    //if initial time is not zero, initialize synccounter
    
    
    //mexPrintf("%d %f %f\n", numSyncTimes, syncTime[1], startTime);
    
    //Check for the maximum sync time smaller than maxTime
    for (i=0; i<numSyncTimes; i++)
    {
        if (syncTime[i] > maxTime)
        {
            effNumSyncTimes = i;
            break;
        }
    }
    //mexPrintf("Num SyncTimes %d\n", numSyncTimes);
    
    for (i=0; i<effNumSyncTimes; i++)
    {
        //mexPrintf("Effective syncTime %f\n", syncTime[i]);
        if (syncTime[i] >= t0)
        {
            syncCounter = i;
            break;
        }
    }
    
    if (syncTime[syncCounter] < t0)
    {
        syncCounter = effNumSyncTimes;
    }
    
   

    for (i=0; i<n; i++)
    {

        CalculatePropensities(h, g, currentX, cTmp, pre, numSpecies, numReactions, 1);
        //cVarying[i] = cTmp[tvRateIndex - 1];
               
        //Set varying+heterogeneous propensities
        if (tvRateHeterogeneous > -1)
        {
            g[tvRateIndex - 1] = cTmp[tvRateIndex - 1] * g[tvRateIndex - 1];
        }
        
        SetValue(G, g, i, numReactions);
        //SetValue(H, h, i, numReactions);
        
        for (k=0; k<numReactions; k++)
        {
            //mexPrintf("%f ", h[k]);
            arrayIdx = GetIndex(populationRateIdx, numPopulationRates, k+1);
            
            if (arrayIdx > -1)
            {
                alpha = a[arrayIdx];
                beta = b[arrayIdx];
                intG = IntG[arrayIdx];
                r = R[arrayIdx];
                
                tau = DrawLomax(alpha, beta, intG, r, g[k]);
            } 
            else
            {
                tau = DrawExponential(1/h[k]);
            }
            
            if ((tau < minTau) || (k == 0 ))
            {
                minTau = tau;
                nextReaction = k;
            }
        }

        //mexPrintf(" minTau=%f \n", minTau);
        
        if (mxIsInf(minTau))
        {
            
            //mexPrintf("No Reaction to fire! Stop...\n");
            SetValue(t, &maxTime, i+1, 1);
            SetValue(x, currentX, i+1, numSpecies); 
            return i+1;
            //break;
        }
        
        currentTime = currentTime + minTau;
        
        //mexPrintf("Maximum time reached (%f s, %f s)! Stop...\n", maxTime, minTau);
        
        if ((currentTime > maxTime) && (syncCounter == effNumSyncTimes))
        {
            //mexPrintf("Maximum time reached (%f s, %f s)! Stop...\n", maxTime, minTau);
            SetValue(t, &maxTime, i+1, 1);
            SetValue(x, currentX, i+1, numSpecies);
            return i+1;
        }
        
        // Do not execute reaction if we need to sync!!
        if (currentTime > syncTime[syncCounter])
        {
            cTmp[tvRateIndex - 1] = 
                    GetPiecewiseConstantInput(inputTime, inputLevels,
                    syncTime[syncCounter], numInputLevels);

            minTau = syncTime[syncCounter] - currentTime + minTau ;
            currentTime = syncTime[syncCounter];
            syncCounter++;
            nextReaction = -1;
        }
        else
        {
            //mexPrintf("Next reaction: %d\n", nextReaction+1);
            ExecuteReaction(currentX, pre, post, nextReaction, numReactions, numSpecies);
        }
        
        //Set current state and time
        SetValue(t, &currentTime, i+1, 1);
        SetValue(x, currentX, i+1, numSpecies); 
        
        react[i] = (double) nextReaction + 1;
        
        for (k=0; k<numPopulationRates; k++)
        {
            // -1 because of MATLAB indexing!!
            IntG[k] = IntG[k] + g[(int)populationRateIdx[k] - 1] * minTau;
            
        }
        
        arrayIdx = GetIndex(populationRateIdx, numPopulationRates, nextReaction+1);

        //mexPrintf("%d %d\n", arrayIdx, nextReaction);
        
        if (arrayIdx > -1)
        {
            R[arrayIdx] = R[arrayIdx] + 1;
            //mexPrintf("%f %d %d %d\n", IntG[arrayIdx], R[arrayIdx], arrayIdx, numPopulationRates);
            //mexPrintf("%f %f %d %d\n", IntG[arrayIdx], R[arrayIdx], numPopulationRates, arrayIdx);
        }
        
     

    }
    
    
    return i;
}






        
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  double *x0, *pre, *post, *a, *b, *initialA, 
         *initialB, *syncTime, *inputTime, *inputLevels,
         *x, *t, *r, *G, *c, *populationRateIdx,
         *cVarying;
    
  int n, tvRateIndex, numPopulationRates, numEvents, numInputLevels;
  
  int i;
  double *Tau;
  
  double maxTime, startTime;
  
  mwSize numSpecies, numReactions, numStates,
         numSyncTimes;
  
  /* Check for proper number of arguments. */
  //if(nrhs!=3) {
  //  mexErrMsgTxt("Three inputs required.");
  //} else if(nlhs>2) {
  //  mexErrMsgTxt("Too many output arguments.");
  //}
  
  /* The input must be a noncomplex scalar double.*/
  numSpecies = mxGetM(prhs[0]);
  numReactions = mxGetM(prhs[1]);
  numPopulationRates = mxGetM(prhs[11]);
  numInputLevels = mxGetN(prhs[7]);
  numSyncTimes = mxGetN(prhs[9]);

  x0 = mxGetPr(prhs[0]);
  c = mxGetPr(prhs[1]);
  pre = mxGetPr(prhs[2]);
  post = mxGetPr(prhs[3]);
  n = mxGetScalar(prhs[4]);
  maxTime = mxGetScalar(prhs[5]);
  tvRateIndex = mxGetScalar(prhs[6]);
  inputTime = mxGetPr(prhs[7]);
  inputLevels = mxGetPr(prhs[8]);
  syncTime = mxGetPr(prhs[9]);
  startTime = mxGetScalar(prhs[10]);
  populationRateIdx = mxGetPr(prhs[11]); 
  a = mxGetPr(prhs[12]);
  b = mxGetPr(prhs[13]);
  initialA = mxGetPr(prhs[14]);
  initialB = mxGetPr(prhs[15]);
  

  /* Create matrix for the return argument. */
  plhs[0] = mxCreateDoubleMatrix(numSpecies, n, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1, n, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(1, n, mxREAL);
  plhs[3] = mxCreateDoubleMatrix(numReactions, n, mxREAL);
  //plhs[4] = mxCreateDoubleMatrix(numReactions, n, mxREAL);
  //plhs[5] = mxCreateDoubleMatrix(1, n, mxREAL);
  
  //plhs[6] = mxCreateDoubleMatrix(1, 100000, mxREAL);
  
  /* Assign pointers to each input and output. */
  x = mxGetPr(plhs[0]);
  t = mxGetPr(plhs[1]);
  r = mxGetPr(plhs[2]);
  //H = mxGetPr(plhs[3]);
  G = mxGetPr(plhs[3]);
  //cVarying = mxGetPr(plhs[5]);

  //for (i=0;i<100000; i++)
  //{     
  //  Tau[i] = DrawLomax(3, 50, 0, 0, 1);
  
  //}
  
  numEvents = SimulatePath(x, t, r, G, x0, c, pre, post, n,
          maxTime, tvRateIndex, inputTime, inputLevels,
          syncTime, startTime, populationRateIdx, a, b, initialA, initialB, 
          numReactions, numSpecies, numPopulationRates, numInputLevels,
          numSyncTimes);
  
  mxSetN(plhs[0], numEvents + 1);
  mxSetN(plhs[1], numEvents + 1);
  mxSetN(plhs[2], numEvents);
  mxSetN(plhs[3], numEvents);
  //mxSetN(plhs[4], numEvents); 
  //mxSetN(plhs[5], numEvents);
 
  //CalculatePropensities(h, g, x, c, pre, numSpecies, numReactions, numStates);
}
