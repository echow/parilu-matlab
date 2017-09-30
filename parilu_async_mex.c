#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "mex.h"
#include "check_sorted.h"

void parilu(int n, int nnz, 
        const int *rowind, const int *colind, const double *val,
        const mwIndex *ial, const mwIndex *jal, double *al,
        const mwIndex *iau, const mwIndex *jau, double *au, 
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

                il = ial[i];
                iu = iau[j];
                while (il < ial[i+1] && iu < iau[j+1])
                {
                    jl = jal[il];
                    ju = jau[iu];

                    if (jl < ju)
                        il++;
                    else if (ju < jl)
                        iu++;
                    else
                    {
                        s -= al[il] * au[iu];
                        il++;
                        iu++;
                    }
                }

                s += al[il-1]*au[iu-1];

                if (i>j)
                    al[il-1] = s / au[iu-1];
                else
                    au[iu-1] = s;
            }
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
    mwIndex *ial, *jal;
    double *al;
    mwIndex *iau, *jau;
    double *au;
    int numiter;
    int numthreads;
    int n, nnz;

    // [l u] = paric_mex(rowind, colind, val, l, u, numiter, numthreads)
    // l must be stored in transposed form (CSC format)
    if (nrhs != 7)
        mexErrMsgTxt("mex function called with bad number of input arguments");
    if (nlhs != 2)
        mexErrMsgTxt("mex function called with bad number of output arguments");

    rowind     = (int *) mxGetData(prhs[0]);
    colind     = (int *) mxGetData(prhs[1]);
    val        =         mxGetPr(prhs[2]);
    numiter    = (int)  *mxGetPr(prhs[5]);
    numthreads = (int)  *mxGetPr(prhs[6]);

    // copy initial l and u and use as output matrix (leave input untouched)
    plhs[0] = mxDuplicateArray(prhs[3]);
    ial     = mxGetJc(plhs[0]);
    jal     = mxGetIr(plhs[0]);
     al     = mxGetPr(plhs[0]);
    plhs[1] = mxDuplicateArray(prhs[4]);
    iau     = mxGetJc(plhs[1]);
    jau     = mxGetIr(plhs[1]);
     au     = mxGetPr(plhs[1]);

    nnz = mxGetM(prhs[0]);
    n   = mxGetM(prhs[3]);

    if (check_sorted(n, ial, jal))
        mexErrMsgTxt("matrix L not sorted");
    if (check_sorted(n, iau, jau))
        mexErrMsgTxt("matrix U not sorted");

    parilu(n, nnz, rowind, colind, val, ial, jal, al, 
                                        iau, jau, au, numiter, numthreads);
}
