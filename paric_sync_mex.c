#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "mex.h"

#define WALLCLOCK(time) do {                                 \
      unsigned long val;                                       \
      volatile unsigned int a, d;                              \
      __asm__ __volatile__("rdtsc" : "=a" (a), "=d" (d) : );   \
      val = ((unsigned long) a)|(((unsigned long)d)<<32);      \
      (time) = val / 3330000000.;                              \
    } while(0)


int check_sorted(int n, const mwIndex *ia, const mwIndex *ja)
{
    int i, j;

    for (i=0; i<n; i++)
    {
        for (j=ia[i]; j<ia[i+1]-1; j++)
        {
            if (ja[j] >= ja[j+1])
            {
                printf("%d not sorted %d %d\n", i, ja[j], ja[j+1]);
                return -1;
            }
        }
    }

    return 0;
}

// compute upper triangular cholesky factor
// this version is synchronous
// warning: values in the initial guess will be overwritten
void parilu_sym_sync(int n, int nnz,
        const int *rowind, const int *colind, const double *val,
        const mwIndex *iau, const mwIndex *jau, double *au, double *ag,
        int numiter, int numthreads)
{
    int iter, k;
    int i, j;
    double s;
    int il, iu, jl, ju;
    int failed = 0; // shared by all threads

        for (iter=1; iter<=numiter; iter++)
        {
#pragma omp parallel for num_threads(numthreads) private(k,i,j,s,il,iu,jl,ju) schedule(dynamic,4096)
            for (k=0; k<nnz; k++)
            {
                i = rowind[k];
                j = colind[k];
                s = val[k];

                il = iau[i];
                iu = iau[j];
                while (il < iau[i+1] && iu < iau[j+1])
                {
                    jl = jau[il];
                    ju = jau[iu];

                    if (jl < ju)
                        il++;
                    else if (ju < jl)
                        iu++;
                    else
                    {
                        // we are going to modify this u entry
                        s -= ag[il] * ag[iu];
                        il++;
                        iu++;
                    }
                }

                // undo the last operation (it must be the last)
                s += ag[il-1]*ag[iu-1];

                // modify u entry
                if (i == j)
                {
                    if (s <= 0.) failed = 1;
                    au[iu-1] = sqrt(s);
                }
                else
                    au[iu-1] = s / ag[il-1];
            }
            // end omp loop

            // copy computed values into guess
            for (k=0; k<iau[n]; k++)
                ag[k] = au[k];

            if (failed)
                mexErrMsgTxt("negative or zero pivot");
        }
}

// matlab sparse matrices are stored by columns
// row indices are indexed with base 0
// ia indices use base 0

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int *rowind, *colind; // must be int32
    const double *val;
    mwIndex *iau, *jau, *iag, *jag;
    double *au, *ag;
    int numiter;
    int numthreads;
    int n, nnz;

    // u = paric_mex(rowind, colind, val, u, numiter, numthreads)
    if (nrhs != 6)
        mexErrMsgTxt("mex function called with bad number of input arguments");
    if (nlhs != 1)
        mexErrMsgTxt("mex function called with bad number of output arguments");

    rowind     = (int *) mxGetData(prhs[0]);
    colind     = (int *) mxGetData(prhs[1]);
    val        =         mxGetPr(prhs[2]);
    numiter    = (int)  *mxGetPr(prhs[4]);
    numthreads = (int)  *mxGetPr(prhs[5]);

    // initial guess u, which is called g
      n     = mxGetM(prhs[3]);
    iag     = mxGetJc(prhs[3]);
    jag     = mxGetIr(prhs[3]);
     ag     = mxGetPr(prhs[3]);
    if (check_sorted(n, iag, jag))
        mexErrMsgTxt("initial guess not sorted");

    // copy initial guess u and use as output matrix (leave input untouched)
    plhs[0] = mxDuplicateArray(prhs[3]);
    iau     = mxGetJc(plhs[0]);
    jau     = mxGetIr(plhs[0]);
     au     = mxGetPr(plhs[0]);
    if (check_sorted(n, iau, jau))
        mexErrMsgTxt("duplicated U not sorted");

    nnz = mxGetM(prhs[0]);

    parilu_sym_sync(n, nnz, rowind, colind, val, iau, jau, au, ag, numiter, numthreads);
}