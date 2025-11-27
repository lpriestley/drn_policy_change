#include "mex.h"
#include <math.h>
#include <limits.h>
#include "df_RLR_MAP_routines.c"

#ifndef MAX
#define MAX(A,B) (((A) > (B)) ? (A) : (B))
#endif

#ifndef MIN
#define MIN(A,B) (((A) > (B)) ? (B) : (A))
#endif

void mexFunction(int             nlhs,      /* No. of output arguments */
                 mxArray         *plhs[],   /* Output arguments. */ 
                 int             nrhs,      /* No. of input arguments. */
                 const mxArray   *prhs[])   /* Input arguments. */
{
   int       nacq = 0;   /* # of acquisistions. */
   int       nvox = 0;   /* # of voxels. */
   int       np = 0;     /* # of priors. */
   int       nm = 0;     /* Noise model. */
   int       m = 0;      /* Scratch pad. */
   int       n = 0;      /* Scratch pad. */
   int       i = 0;
   int       pi[8];      /* Indicies of priors. */
   double    *sv = NULL;
   double    *y = NULL;
   double    *b = NULL;
   double    *g = NULL;
   double    *pmu = NULL;
   double    *pvar = NULL;
   double    *theta = NULL;
   double    *tmp = NULL;
   double    *ei = NULL;
   double    *hist = NULL;
   

   if (nrhs == 0) mexErrMsgTxt("usage: theta = df_MAP(theta,y,b,g,nm,pr,pri)");
   if ((nrhs!=4) && (nrhs!=5) && (nrhs!=7)) mexErrMsgTxt("df_MAP: 4, 5 or 7 input arguments required");
   if ((nlhs!=1) && (nlhs!=2) && (nlhs!=3)) mexErrMsgTxt("df_MAP: 1, 2 or 3 output arguments required");

   /* 
      Get second argument, which should be an mxn matrix, where
      m is number of acquisitions and n is number of voxels.
   */ 

   if (!mxIsNumeric(prhs[1]) || mxIsComplex(prhs[1]) || mxIsSparse(prhs[1]) || !mxIsDouble(prhs[1]))
   {
      mexErrMsgTxt("df_MAP: Y must be numeric, real, full and double");
   }
   if ((nacq = mxGetM(prhs[1])) > MLEN)
   {
      mexErrMsgTxt("df_MAP: Maximum # of acquisitons exceeded"); /* To spot transposed Y. */
   }
   if ((nvox = mxGetN(prhs[1])) > MVOX)
   {
      mexErrMsgTxt("df_MAP: Mismatch between # of start vectors and data");
   }
   y = mxGetPr(prhs[1]);

   /*
      Get first argument, which should be an (1+7n)x1 vector, where
      n is number of voxels. These are the starting estimates of
      the parameters. 
   */
   if (!mxIsNumeric(prhs[0]) || mxIsComplex(prhs[0]) || mxIsSparse(prhs[0]) || !mxIsDouble(prhs[0]))
   {
      mexErrMsgTxt("df_MAP: theta must be numeric, real, full and double");
   }
   m = mxGetM(prhs[0]); n = mxGetN(prhs[0]);
   if ((m != 1+7*nvox) || (n != 1))
   {
      mexErrMsgTxt("df_MAP: Wrong format of theta, must be (1+7*nvox)*n");
   }
   sv = mxGetPr(prhs[0]);

   /*
      Get third argument. Should be an mx1 vector with b-values.
   */
   
   if (!mxIsNumeric(prhs[2]) || mxIsComplex(prhs[2]) || mxIsSparse(prhs[2]) || !mxIsDouble(prhs[2]))
   {
      mexErrMsgTxt("df_MAP: b must be numeric, real, full and double");
   }

   if (MIN(mxGetM(prhs[2]),mxGetN(prhs[2])) != 1)
   {
      mexErrMsgTxt("df_MAP: b must be mx1 vector");
   }
   if (MAX(mxGetM(prhs[2]),mxGetN(prhs[2])) != nacq)
   {
      mexErrMsgTxt("df_MAP: b must be mx1 vector");
   }
   b = mxGetPr(prhs[2]);

   /*
      Get fourth argument. Should be an mx3 vector with gradients.
   */
   
   if (!mxIsNumeric(prhs[3]) || mxIsComplex(prhs[3]) || mxIsSparse(prhs[3]) || !mxIsDouble(prhs[3]))
   {
      mexErrMsgTxt("df_MAP: g must be numeric, real, full and double");
   }
   if ((mxGetM(prhs[3]) != nacq) || (mxGetN(prhs[3]) != 3))
   {
      mexErrMsgTxt("df_MAP: g must be mx3 matrix");
   }
   g = mxGetPr(prhs[3]);

   /*
      Get fifth argument. Should be a scalar with the value 1 or 2.
   */

   if (nrhs > 4)
   {
      if (!mxIsNumeric(prhs[4]) || mxIsComplex(prhs[4]) || mxIsSparse(prhs[4]) || !mxIsDouble(prhs[4]))
      {
         mexErrMsgTxt("df_MAP: nm must be numeric, real, full and double");
      }
      if ((MAX(mxGetM(prhs[4]),mxGetN(prhs[4])) != 1) || 
          ((mxGetScalar(prhs[4]) != 2.0) && (mxGetScalar(prhs[4]) != 1.0)))
      {
         mexErrMsgTxt("df_MAP: nm should be a scalar 1 or 2");
      }
      nm = ((int) mxGetScalar(prhs[4]));
   }
   else
   {
      nm = 2; /* Rician noise defualt. */ 
   }
   
   if (nrhs > 5)
   {
      /* 
         Get sixth argument. Should be a kx2 matrix,
         where k is the # of parameters with priors.
      */
      if (!mxIsNumeric(prhs[5]) || mxIsComplex(prhs[5]) || mxIsSparse(prhs[5]) || !mxIsDouble(prhs[5]))
      {
         mexErrMsgTxt("df_MAP: pr must be numeric, real, full and double");
      }
      if ((np = mxGetM(prhs[5])) > 8 || mxGetN(prhs[5]) != 2) 
      {
         mexErrMsgTxt("df_MAP: pr must be a kx2 matrix where k<=8");
      }
      pmu = mxGetPr(prhs[5]);
      pvar = &(pmu[np]);

      /* 
         Get seventh argument. Should be a kx1 matrix,
         where k is the # of parameters with priors.
      */
      if (!mxIsNumeric(prhs[6]) || mxIsComplex(prhs[6]) || mxIsSparse(prhs[6]) || !mxIsDouble(prhs[6]))
      {
         mexErrMsgTxt("df_MAP: pi must be numeric, real, full and double");
      }
      if (MAX(mxGetM(prhs[6]),mxGetN(prhs[6])) != np || MIN(mxGetM(prhs[6]),mxGetN(prhs[6])) != 1) 
      {
         mexErrMsgTxt("df_MAP: mismatch between pr and pi");
      }
      tmp = mxGetPr(prhs[6]);
      for (i=0; i<np; i++) {pi[i] = ((int) (tmp[i]-1.0+0.2));}
   }
   else
   {
      np = 0; pmu = NULL; pvar = NULL;
   }

   /*
      Allocate memory for output.
   */

   plhs[0] = mxCreateDoubleMatrix(1+7*nvox,1,mxREAL);
   theta = mxGetPr(plhs[0]);
   for (i=0; i<(1+7*nvox); i++)
   {
      theta[i] = sv[i];
   }
   if (nlhs > 1) 
   {
      plhs[1] = mxCreateDoubleMatrix(nacq,nvox,mxREAL);
      ei = mxGetPr(plhs[1]);
   }
   else {ei = (double *) mxCalloc(nacq*nvox,sizeof(double));}
   if (nlhs > 2)  
   {
      plhs[2] = mxCreateDoubleMatrix(MAXITER,1,mxREAL);
      hist = mxGetPr(plhs[2]);
   }
   else {hist = (double *) mxCalloc(MAXITER,sizeof(double));}

   /*
      Do the fitting.
   */
    
   if (fit_RLR(g,b,y,nacq,nvox,pmu,pvar,pi,np,nm,theta,ei,hist) < 0)
   {
      for (i=0; i<(1+7*nvox); i++)
      {
         theta[i] = mxGetNaN();
      }
   } 

   if (nlhs < 3) {mxFree(hist);}
   if (nlhs < 2) {mxFree(ei);}

   return;
}
