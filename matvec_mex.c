#include <omp.h>
#include "matrix.h"
#include "mex.h"

static void 
Matvec(int n, const mwIndex *ia, const mwIndex *ja, const double *a,
       const double *x, double *y)
{
    int i, j;
    double t;

#pragma omp parallel for private(i,j,t)
    for (i=0; i<n; i++)
    {
        t = 0.;
        for (j=ia[i]; j<ia[i+1]; j++)
            t += a[j]*x[ja[j]];
        y[i] = t;
    }
}

// y = matvec(a, x);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int n;

    if (nrhs != 2)
        mexErrMsgTxt("matvec mex function called with bad number of arguments");

    n = mxGetN(prhs[0]);
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL); // create solution vector

    Matvec(n, mxGetJc(prhs[0]), mxGetIr(prhs[0]), mxGetPr(prhs[0]),
        mxGetPr(prhs[1]), mxGetPr(plhs[0]));
}
