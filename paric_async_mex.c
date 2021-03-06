#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "mex.h"
#include "check_sorted.h"

void paric_async(int n, int nnz,
        const int *rowind, const int *colind, const double *val,
        const mwIndex *iau, const mwIndex *jau, double *au, 
        int numiter, int numthreads)
{
    int iter, k;
    int i, j;
    double s;
    int il, iu, jl, ju;
    int failed = 0; // shared by all threads

#pragma omp parallel num_threads(numthreads) private(iter,k,i,j,s,il,iu,jl,ju)
    {
        for (iter=1; iter<=numiter; iter++)
        {
#pragma omp for schedule(dynamic,4096) nowait
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
                        s -= au[il] * au[iu];
                        il++;
                        iu++;
                    }
                }

                // undo the last operation (it must be the last)
                s += au[il-1]*au[iu-1];

                // modify u entry
                if (i == j)
                {
                    if (s <= 0.) failed = 1;
                    au[iu-1] = sqrt(s);
                }
                else
                    au[iu-1] = s / au[il-1];
            }
            // end omp loop

            if (failed)
                mexErrMsgTxt("negative or zero pivot");
        }
    }
}

// matlab sparse matrices are stored by columns
// row indices are indexed with base 0
// ia indices use base 0

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int *rowind, *colind; // must be int32
    const double *val;
    mwIndex *iau, *jau;
    double *au;
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

    // copy initial guess u and use as output matrix (leave input untouched)
    plhs[0] = mxDuplicateArray(prhs[3]);
    iau     = mxGetJc(plhs[0]);
    jau     = mxGetIr(plhs[0]);
     au     = mxGetPr(plhs[0]);

    nnz = mxGetM(prhs[0]);
    n   = mxGetM(prhs[3]);

    if (check_sorted(n, iau, jau))
        mexErrMsgTxt("matrix U not sorted");

    paric_async(n, nnz, rowind, colind, val, iau, jau, au, numiter, numthreads);
}
