#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "mex.h"
#include "check_sorted.h"

void pariutdu(int n, int nnz,
        const int *rowind, const int *colind, const double *val,
        const mwIndex *iau, const mwIndex *jau, double *au, 
        double *ad, // diagonal entries
        int numiter, int numthreads)
{
    int iter, k;
    int i, j;
    double s;
    int il, iu, jl, ju;

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
                        s -= au[il] * au[iu] * ad[jl];
                        il++;
                        iu++;
                    }
                }

                // undo the last operation (it must be the last)
                s += au[il-1]*au[iu-1]*ad[i];

                // modify u entry
                if (i == j)
                    ad[i] = s;
                else
                    au[iu-1] = s / ad[i];
            }
            // end omp loop
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
    mwIndex *iad, *jad;
    double *ad;
    int numiter;
    int numthreads;
    int n, nnz;

    // [u d] = pariutdu_mex(rowind, colind, val, u, d, numiter, numthreads)
    if (nrhs != 7)
        mexErrMsgTxt("mex function called with bad number of input arguments");
    if (nlhs != 2)
        mexErrMsgTxt("mex function called with bad number of output arguments");

    rowind     = (int *) mxGetData(prhs[0]);
    colind     = (int *) mxGetData(prhs[1]);
    val        =         mxGetPr(prhs[2]);
    numiter    = (int)  *mxGetPr(prhs[5]);
    numthreads = (int)  *mxGetPr(prhs[6]);

    // copy initial guess u and use as output matrix (leave input untouched)
    plhs[0] = mxDuplicateArray(prhs[3]);
    iau     = mxGetJc(plhs[0]);
    jau     = mxGetIr(plhs[0]);
     au     = mxGetPr(plhs[0]);
    plhs[1] = mxDuplicateArray(prhs[4]);
    iad     = mxGetJc(plhs[1]);
    jad     = mxGetIr(plhs[1]);
     ad     = mxGetPr(plhs[1]);

    nnz = mxGetM(prhs[0]);
    n   = mxGetM(prhs[3]);

    if (check_sorted(n, iau, jau))
        mexErrMsgTxt("matrix U not sorted");

    // check that d is diagonal
    if (iad[n] != n)
        mexErrMsgTxt("input D not diagonal");

    pariutdu(n, nnz, rowind, colind, val, iau, jau, au, ad, numiter, numthreads);
}
