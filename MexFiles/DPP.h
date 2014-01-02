/* -----------------------------------------------------
 * Objective: Commonly used helper functions
 *------------------------------------------------------
 * Author: Christoph Zechner
 *-----------------------------------------------------*/

#include "mex.h"
#include "math.h"

/* -----------------------------------------------------
 * Name: Factorial
 *------------------------------------------------------
 * Purpose: Calculates factorial (used for propensities).
 *-----------------------------------------------------*/
int Factorial(int x)
{
    if (x == 0)
    {
        return 1;
    }
    
    return x * Factorial(x - 1);
}


/* -----------------------------------------------------
 * Name: NChooseK
 *------------------------------------------------------
 * Purpose: Calculates binomial coefficient 
 *          ("n choose k").
 *-----------------------------------------------------*/
int NChooseK(int n, int k)
{
    /* return Factorial(n) / (Factorial(k) * Factorial(n - k)) */
    
    int i, prod;
    
    prod = 1;
    
    for (i=1; i<=k; i++)
    {
        prod = prod * (n - k + i) / i;
    }
    
    return prod;
}

/* -----------------------------------------------------
 * Name: CalculatePropensities
 *------------------------------------------------------
 * Purpose: Calculates mass-action propensities 
 *          (i.e., "c*g(x)").
 *-----------------------------------------------------*/
void CalculatePropensities(double *h, double *g, double *x,
        double *c, double *pre, mwSize numSpecies, mwSize numReactions,
        mwSize numStates)
{
    
    mwSize i;
    mwSize k;
    mwSize l;
    mwSize propIdx;
    
    int stochProduct;
    int preVal;
    int stateValue;

    for (l=0; l<numStates; l++)
    {
        for (i=0; i<numReactions; i++)
        {
            stochProduct = 1;
            
            for (k=0; k<numSpecies; k++)
            {
                preVal = (int)pre[numReactions * k + i];
                stateValue = (int)x[l*numSpecies + k];
                
                
                if (preVal > 0)
                {
                    if (stateValue >= preVal)
                    {
                        stochProduct = stochProduct * NChooseK(stateValue, preVal);
                    }
                    else
                    {
                        stochProduct = 0;
                        break;
                    }
                }
            }

            propIdx = l*numReactions + i;
            g[propIdx] = stochProduct;
            h[propIdx] = c[i] * stochProduct;
            
            //printf("State %d: %d (%d) \n", k, stateValue, stochProduct);
        }
    } 
}

/* -----------------------------------------------------
 * Name: SetValue
 *------------------------------------------------------
 * Purpose: Sets the row of a matrix (at a given index).
 *-----------------------------------------------------*/
void SetValue(double *X, double *x, int index, mwSize numRows)
{
    int i;
    
    for (i = 0; i<numRows; i++)
    {
        X[numRows * index + i] = x[i];
    }
        
}

/* -----------------------------------------------------
 * Name: GetValue
 *------------------------------------------------------
 * Purpose: Gets the row of a matrix (at a given index).
 *-----------------------------------------------------*/
void GetValue(double *x, double *X, int index, mwSize numRows)
{
    int i;
    
    for (i = 0; i<numRows; i++)
    {
        x[i] = X[numRows * index + i];
    }

}

/* -----------------------------------------------------
 * Name: DrawExponential
 *------------------------------------------------------
 * Purpose: Draws an exponential waiting time with mean
 *          mu.
 *-----------------------------------------------------*/
double DrawExponential(double mu)
{
 
    double u;
    u = (double) rand() / RAND_MAX;
    
    return mu * log(1/u);
    
}


/* -----------------------------------------------------
 * Name: DrawLomax
 *------------------------------------------------------
 * Purpose: Draws a Lomax-distributed waiting time.
 *-----------------------------------------------------*/
double DrawLomax(double a, double b, double intG, double r, double g)
{
    double u, n;
    n = r;
 
    u = (double) rand() / RAND_MAX;
    
    return (double) -1*(intG + b - exp((log(1/u) + 
            log(intG + b)*(a + n))/(a + n))) / g;
}


/* -----------------------------------------------------
 * Name: GetIndex
 *------------------------------------------------------
 * Purpose: Find index of a given element in a vector.
 *-----------------------------------------------------*/
int GetIndex(double *array, int size, int value)
{
    int i;
    
    for (i=0; i<size; i++)
    {

        //Add 1 since MATLAB indexing starts at 1!!!
        if ((int)array[i] == value)
        {
            return i;
        }
    }
    
    return -1;
        
}

/* -----------------------------------------------------
 * Name: GetPiecewiseConstantInput
 *------------------------------------------------------
 * Purpose: Returns value of piece-wise-constant 
 *          function at a given time.
 *-----------------------------------------------------*/
double GetPiecewiseConstantInput(double *inputTimes, double *inputLevels, 
        double t, int numInputLevels)
{

    int i, counter;
    
    counter = 0;
    
    if (t < inputTimes[0])
    {
        //mexPrintf("Smaller: %f %f %d\n", t, inputTimes[0], numInputLevels);
        return inputLevels[0];
    }
    
    if (t > inputTimes[numInputLevels - 1])
    {
        //mexPrintf("Bigger: %f %f %d\n", t, inputTimes[0], numInputLevels);
        return inputLevels[numInputLevels - 1];
    }
    
    for (i=0; i<numInputLevels; i++)
    {
        if (t >= inputTimes[i])
        {
            counter++;
            //mexPrintf("NORMAL: %f %f %d\n", t, inputTimes[i], numInputLevels);
        }
    }
    
    return inputLevels[counter];
}



