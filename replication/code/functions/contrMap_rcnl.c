/*==========================================================
 * c_contrMap_rcnl.c - modifies Nevo (2001) contraction mapping
 * to accommodate an "inside good" nest.  
 *
 * Inputs include:
 *      tol = scalar indicating the contraction mapping tolerance
 *      maxiter = scalar indicating maximimum number of iterations
 *      expmvalold = vector as in Nevo (2001) Matlab code
 *      expmu = matrix as in Nevo (2001) Matlab code
 *      cdindex = vector as in Nevo (2001) Matlab code
 *      rho = nested logit parameter (rho=0 collapses nests)
 *
 * Output is a vector, named expmval internally, that is the exponentiated 
 * mean consumer valuations that align observed and predicted shares given 
 * the nesed logit parameter and the candidate nonlinear parameters that 
 * give rise to the expmu matrix.  Auxiliary output vectors are the market-
 * specific tolerances and iterations.
 *
 * The calling syntax is:
 *
 * [expmval,mkttol,mktiter] = c_contMap_rcnl(tol,maxiter,expmvalold,expmu,cdindex,rho)
 *
 *========================================================*/
/* $Written by Nathan Miller, January 30, 2016$ */

#include "mex.h"
#include "math.h"


//mexPrintf("3 to the power of 12 equals: %f\n",pow(3,12));

/*==========================================================
 * The printArray function makes it simpler to print arrays,
 * thereby faciliating the debugging process. 
 *========================================================*/
void printArray(int size, char* name, double *array) {
    int j;
    for (j=0; j<size; j++){
        mexPrintf("%s[%d]=%f\n", name, j, array[j]);
    }
}


/*==========================================================
 * The contrMap function contains the computational routine.
 *========================================================*/
void contrMap(double tol, double maxiter, 
        double *expmvalold, double *expmval, double *expmu,
        mwSize nJ, mwSize nI, mwSize nM, 
        double *cdindex, double *shares, double rho, 
        double *mktnorm, double *mktiter)
{
    /* Assigning variable names */
    mwSize j;
    mwSize i;
    mwSize mkt;
    mwSize ii;
    mwSize k;
    mwSize nJM; 
    int mStart;
    int mEnd;
    double invnI;  
    double maxabsdev;
    double absdev;
    double norm;      
    int iter;
    mxArray *expmvaloldMx;
    double  *expmvaloldM; 
    mxArray *sharesMx;
    double  *sharesM;     
    mxArray *expmuMx;
    double  *expmuM;     
    double   incVal;    
    mxArray *exputilMx;
    double  *exputilM;    
    mxArray *pSharesiMx;
    double  *pSharesiM;     
    mxArray *pSharesMx;
    double  *pSharesM; 
    mxArray *expmvalMx;
    double  *expmvalM;    
    double  oneLessRho;
    double  invOneLessRho;

    mxArray *mvaloldMx;
    double  *mvaloldM; 
    mxArray *mvalMx;
    double  *mvalM;     
    
    // Moves division outside the for loops //
    invnI = 1.0 / (double)nI;
    oneLessRho = (1-rho); 
    invOneLessRho = 1 / (double)oneLessRho;
    
    //mexPrintf("invnI: %f \n",invnI);
    //mexPrintf("nI: %d \n",nI);
    //printArray(nM, "cdindex", cdindex);    
    //mexPrintf("oneLessRho: %f \n",oneLessRho);
    //mexPrintf("invOneLessRho: %f \n",invOneLessRho);
    
    /*==========================================================
     * Looping through markets creates time savings because the iterations
     * for each market end when the convergence criterion is met within the
     * market, rather than when it is met across all markets.
     *========================================================*/    
    for (mkt=0; mkt<nM; mkt++) {
        
        /* Start/stop variables to help create market-specific objects */
        if (mkt==0){
            mStart = 0;
            nJM = cdindex[mkt];
        }
        else {
            mStart = cdindex[mkt-1];
            nJM = cdindex[mkt] - cdindex[mkt-1];
        }           
        mEnd = cdindex[mkt];
            
        /* Create market-specific objects: expmvalold, expmu, shares*/
        expmvaloldMx = mxCreateDoubleMatrix((mwSize)nJM,1,mxREAL); 
        expmvaloldM  = mxGetPr(expmvaloldMx);

        sharesMx = mxCreateDoubleMatrix((mwSize)nJM,1,mxREAL); 
        sharesM  = mxGetPr(sharesMx);
            
        expmuMx = mxCreateDoubleMatrix(nJM*nI,1,mxREAL); 
        expmuM  = mxGetPr(expmuMx);            
                 
        /* Fill in market-specific objects with input data*/            
        for (j=mStart; j<mEnd; j++) {            
            expmvaloldM[j-mStart] = expmvalold[j];  
            //mexPrintf("expmvaloldM[%d]: %f \n",j-mStart,expmvaloldM[j-mStart]);
            sharesM[j-mStart] = shares[j];  
            for (i=0; i<nI; i++){
                expmuM[nJM*i+j-mStart] = expmu[nJ*i+j];  
                //mexPrintf("expmuM[%d]: %f \n",nJM*i+j-mStart,expmuM[nJM*i+j-mStart]);
            }
        }
        
        /* Market-specific intermediate data sets */
        exputilMx = mxCreateDoubleMatrix(nJM*nI,1,mxREAL); 
        exputilM  = mxGetPr(exputilMx);  
        
        pSharesiMx = mxCreateDoubleMatrix(nJM*nI,1,mxREAL); 
        pSharesiM  = mxGetPr(pSharesiMx);  

        pSharesMx = mxCreateDoubleMatrix(nJM,1,mxREAL); 
        pSharesM  = mxGetPr(pSharesMx);    

        expmvalMx = mxCreateDoubleMatrix(nJM,1,mxREAL); 
        expmvalM  = mxGetPr(expmvalMx);            
        
        mvalMx = mxCreateDoubleMatrix(nJM,1,mxREAL); 
        mvalM  = mxGetPr(mvalMx); 
        
        mvaloldMx = mxCreateDoubleMatrix(nJM,1,mxREAL); 
        mvaloldM  = mxGetPr(mvaloldMx);       
        
        
         /*==========================================================
         * Conducting the market-specific contraction mapping.
         *========================================================*/  
        iter = 0;
        norm = 1.0;       
        while (norm>tol && iter<maxiter) {        
            maxabsdev=0;
            
            /* Calculating market-individual-specific shares */
            for (ii=0; ii<nI; ii++) {
                incVal = 0;
                for (k=0; k<nJM; k++) {
                    exputilM[(ii*nJM)+k] = pow(expmvaloldM[k]*expmuM[(ii*nJM)+k],invOneLessRho);
                    incVal = incVal + exputilM[(ii*nJM)+k];
                    //mexPrintf("incVal in loop: %f \n",incVal);
                }
                incVal = oneLessRho*log(incVal); 
                //mexPrintf("incVal out loop: %f \n",incVal);   
                //mexPrintf("incVal2 out loop: %f \n",log(1+exp(incVal)));   
                for (k=0; k<nJM; k++) {
                    pSharesiM[(ii*nJM)+k] = exputilM[(ii*nJM) + k] * exp(incVal) / exp(incVal*invOneLessRho) / exp(log(1+exp(incVal)));
                    //mexPrintf("pSharesiM: %f \n", exputilM[(ii*nJM) + k] * exp(incVal) / exp(incVal*invOneLessRho) / exp(log(1+exp(incVal))));   
                }
            }
            /*========================================================= 
            * Aggregates to obtain predicted market-specific shares,
            * calculates expmval and deviations, and updates expmvalold
            *=========================================================*/
            for (k=0; k<nJM; k++) {
                pSharesM[k]=0;
                for (ii=0; ii<nI; ii++) {
                    pSharesM[k] = pSharesM[k]+(pSharesiM[(ii*nJM)+k] * invnI);
                }
                //mexPrintf("pSharesM: %f \n", pSharesM[k]);
                
                expmvalM[k] = expmvaloldM[k] * pow(sharesM[k]/pSharesM[k],oneLessRho);
                //expmvalM[k] = expmvaloldM[k] * sharesM[k] / pSharesM[k];
                
                //mvaloldM[k]=log(expmvaloldM[k]);
                //mvalM[k]= mvaloldM[k] + oneLessRho*(log(sharesM[k])-log(pSharesM[k]));
                //expmvalM[k] = exp(mvalM[k]);
                
                absdev = fabs(expmvalM[k]-expmvaloldM[k]);
                //mexPrintf("absdev: %f \n",absdev);
                if (absdev>maxabsdev){
                    maxabsdev=absdev;
                }
                //mexPrintf("maxabsdev: %f \n",maxabsdev);               
                expmvaloldM[k]=expmvalM[k];
            }
            //mexPrintf("maxabsdev: %f \n",maxabsdev);
            //printArray(nJM,"sharesM", sharesM);
            //printArray(nJM,"pSharesM", pSharesM);
            iter = iter+1;
            norm = maxabsdev;
        }

        /*=========================================================
         * Slots market-specific expmvalM into expmval after convergence.
         * The vector exmpval is what is returned to the gateway function
         * and also to Matlab.
         *=========================================================*/
        for (j=mStart; j<mEnd; j++) {            
            expmval[j] = expmvalM[j-mStart];  
        }
        mktnorm[mkt]=norm;
        mktiter[mkt]=iter;
        
        //mexPrintf("nJM: %d \n",nJM);
        //mexPrintf("mStart: %d \n",mStart);
        //mexPrintf("mEnd: %d \n",mEnd);
        //printArray(nJM,"expmvaloldM", expmvaloldM);
        //printArray(nJM,"sharesM", sharesM);
        //printArray(nJM*nI,"expmuM", expmuM);
        //printArray(nJM*nI,"exputilM", exputilM);
        //printArray(nJM*nI,"pSharesiM", pSharesiM);
        //printArray(nJM,"pSharesM", pSharesM);        
        //printArray(nJM,"expmvalM", expmvalM);   
        
        /* Destroying market-specific data to prevent memory leaks */
        mxDestroyArray(expmvaloldMx);
        mxDestroyArray(expmuMx);
        mxDestroyArray(sharesMx);
        mxDestroyArray(exputilMx);
        mxDestroyArray(pSharesiMx);        
        mxDestroyArray(pSharesMx);
        mxDestroyArray(expmvalMx);  
        mxDestroyArray(mvalMx);  
        mxDestroyArray(mvaloldMx);
        
        /* *** MARKET LOOP ENDS *** */
    }
   
   //mexPrintf("Hello!\n");
   //mexPrintf("norm %f \n",norm);
   //mexPrintf("iter %d \n",iter);
   //printArray(nM, "cdindex", cdindex);    
   //printArray(nJ, "expmvalold", expmvalold);
   //printArray(nJ*nI, "expmu", expmu);
   //printArray(nJ, "expmval", expmval);

/* *** END OF contrMap FUNCTION *** */
}
    



/*====================================================================
 * mexFunction is the gateway beween Matlab and the contrMap function.
 *===================================================================*/
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double tol;          /* contraction mapping tolerance */
    double maxiter;      /* max iterations in contraction mapping */    
    double *expmvalold;  /* initial vector of mean valuations */
    double *expmu;       /* matrix of individual-specific deviations */   
    double *cdindex;     /* vector identifies last obs of each market */
    double *shares;      /* vector of observed market shares */
    double rho ;         /* nested logit parameter */

    size_t nJ ;          /* number of product-mkt observations */
    size_t nI ;          /* number of individual draws  */
    size_t nM ;          /* number of markets  */    
    double *expmval;     /* output vector of updated mean valuations */    
    double *mktnorm;     /* output vector of market tolerances */    
    double *mktiter;     /* output vector of market # iterations */  
    
    /* get the value of tol  */
    tol = mxGetScalar(prhs[0]);

    /* get the value of maxiter  */
    maxiter = mxGetScalar(prhs[1]);    
    
    /* create a pointer to the real exmvalold data  */
    expmvalold = mxGetPr(prhs[2]);

    /* get dimensions of the inputs */
    nJ = mxGetM(prhs[3]);
    nI = mxGetN(prhs[3]);    
    nM = mxGetM(prhs[4]);
    
    // mexPrintf("nJ=%d; nI=%d; nM=%d\n", nJ, nI, nM);
    
    /* create a pointer to the real expmu data  */
    expmu = mxGetPr(prhs[3]);

    /* create a pointer to the real cdindex data */
    cdindex = mxGetPr(prhs[4]);    

    /* create a pointer to the real share data  */
    shares = mxGetPr(prhs[5]);      

    /* create a pointer to the nested logit parameter */
    rho = mxGetScalar(prhs[6]);
        
    /* create the output data 1 */
    plhs[0] = mxCreateDoubleMatrix((mwSize)nJ,1,mxREAL);

    /* create the output data 2 */
    plhs[1] = mxCreateDoubleMatrix((mwSize)nM,1,mxREAL);    
 
    /* create the output data 2 */
    plhs[2] = mxCreateDoubleMatrix((mwSize)nM,1,mxREAL);      
    
    /* get a pointer to the real data in the output matrix */
    expmval = mxGetPr(plhs[0]);

    /* get a pointer to the realized tolerance */
    mktnorm = mxGetPr(plhs[1]);

    /* get a pointer to the realized tolerance */
    mktiter = mxGetPr(plhs[2]);    
    
    /* call the computational routine */
    contrMap(tol,maxiter,expmvalold,expmval,expmu,
            (mwSize)nJ,(mwSize)nI,(mwSize)nM,
            cdindex,shares,rho,mktnorm,mktiter);    

    
    
  //  printArray(nJ, "final expmval", expmval);
  //  printArray(nJ, "final exputil", exputil);
  //  mexPrintf("Done!\n");   
    
    
}
  
    
  