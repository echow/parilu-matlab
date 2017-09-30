#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "mex.h"
#include "check_sorted.h"

void pariutdu_sync(int n, int nnz,
        const int *rowind, const int *colind, const double *val,
        const mwIndex *iau, const mwIndex *jau, double *au, 
        double *ad, // diagonal entries
        double *ag, double *agd, // frozen factors
        int numiter, int numthreads)
{
    int iter, k;
    int i, j;
    double s;
    int il, iu, jl, ju;

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
                        s -= ag[il] * ag[iu] * agd[jl];
                        il++;
                        iu++;
                    }
                }

                // undo the last operation (it must be the last)
                s += ag[il-1]*ag[iu-1]*agd[i];

                // modify u entry
                if (i == j)
                    ad[i] = s;
                else
                    au[iu-1] = s / agd[i];
            }
            // end omp loop

            // copy computed values into guess
            // but only if there is another iteration
            if (iter < numiter)
            {
                for (k=0; k<iau[n]; k++)
                    ag[k] = au[k];
                for (k=0; k<n; k++)
                    agd[k] = ad[k];
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
    double *ag, *agd; // frozen factors
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

    // copy initial guess u,d and use as output matrix (leave input untouched)
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

    // allocate temp arrays if more than one sweep
    mxArray *utemp, *dtemp;
    if (numiter == 1)
    {
        utemp = (mxArray *) prhs[3]; // discard const
        dtemp = (mxArray *) prhs[4]; // discard const
    }
    else
    {
        utemp = mxDuplicateArray(prhs[3]);
        dtemp = mxDuplicateArray(prhs[4]);
    }

    // get pointers to frozen factors
    ag  = mxGetPr(utemp);
    agd = mxGetPr(dtemp);

    pariutdu_sync(n, nnz, rowind, colind, val, iau, jau, au, ad, 
        ag, agd, numiter, numthreads);
}
