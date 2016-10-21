/* The matrix is saved in the column-major order.
 * Instead of sorting Ir array after computation, I build up a supplementary
 * array, and rearrange Ir array with one pass O(nrow), but this involves
 * allocating an array of size nrow. This should be fast if (A_deriv*b_val)
 * is not too sparse, but if it is really sparse, it might be faster to
 * sort instead, which will give approximately
 * O(nderiv * nnz_per_column * log(nnz_per_column)). */

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{

    /* Read In Value Vector */
    mwSize nrow;
    nrow = mxGetM(prhs[0]);

    /* Read In Derivative Matrix */
    mwSize nderiv, nnz;
    mwIndex *irsA,*jcsA;
    double *prA;
    nderiv  = mxGetN(prhs[1]);
    irsA    = mxGetIr(prhs[1]);
    jcsA    = mxGetJc(prhs[1]);
    prA     = mxGetPr(prhs[1]);
    nnz     = jcsA[nderiv];
    
    /* Read in Vector V */
    mwSize counter=0;

    /* Allocate Data Matrix for Output*/
    double *srA;
    mwIndex *lirs,*ljcs;
    lirs    = mxMalloc( nnz * sizeof(*lirs));
    ljcs    = mxMalloc( (nderiv+1) * sizeof(*ljcs));
    srA     = mxMalloc( nnz * sizeof(*srA));
    
    /* Logical Vector */
    if (mxIsLogical(prhs[0])) {
        bool *prV;
        prV = mxGetLogicals(prhs[0]);
    
        /* For a full vector */
        if (!mxIsSparse(prhs[0])) {
            mwIndex k,j,i;
            for (j=0 ; j < nderiv ; ++j) {
                ljcs[j] = counter;
    
                for (i=jcsA[j]; i<jcsA[j+1] ; ++i) {
                    k = irsA[i];
                    if (prV[k]) {
                        lirs[counter]=k;
                        srA[counter] = prA[i] * prV[k];
                        ++counter;
                    }
                }
            }
            ljcs[nderiv]=counter;
    
            mxRealloc(lirs, counter * sizeof(*lirs));
            mxRealloc(srA, counter * sizeof(*srA));
        }
        /* For a Sparse Vector */
        else {
            mwIndex *irsV,*jcsV,*pos,old,i;
            irsV    = mxGetIr(prhs[0]);
            jcsV    = mxGetJc(prhs[0]);
            pos     = mxMalloc(nrow * sizeof(*pos));
    
            old=0;
            for (i=0; i<jcsV[1]; ++i) {
                for (; old<irsV[i]; ++old) {
                    pos[old] = -1;
                }
                pos[old] = i;
                ++old;
            }
            for (i=old; i<nrow; ++i) {
                pos[i]=-1;
            }
    
            /* Multiply */
            mwIndex k,j,tmp;
            for (j=0 ; j < nderiv ; ++j) {
                ljcs[j] = counter;
    
                for (i=jcsA[j]; i<jcsA[j+1] ; ++i) {
                    k = irsA[i];
                    if (pos[k] != (-1)) {
                        lirs[counter]= k;
                        srA[counter] = prA[i] * prV[pos[k]];
                        ++counter;
                    }
                }
            }
            ljcs[nderiv]=counter;
    
            mxRealloc(lirs, counter * sizeof(*lirs));
            mxRealloc(srA, counter * sizeof(*srA));
        }
    }
    /* Real Vector */
    else {
        double *prV;
        prV  = mxGetPr(prhs[0]);
    
        /* For a full vector */
        if (!mxIsSparse(prhs[0])) {
            mwIndex k,j,i;
            for (j=0 ; j < nderiv ; ++j) {
                ljcs[j] = counter;
    
                for (i=jcsA[j]; i<jcsA[j+1] ; ++i) {
                    k = irsA[i];
                    if (prV[k]!=0.0) {
                        lirs[counter]=k;
                        srA[counter] = prA[i] * prV[k];
                        ++counter;
                    }
                }
            }
            ljcs[nderiv]=counter;
    
            /* Reallocate Output */
            mxRealloc(lirs, counter * sizeof(*lirs));
            mxRealloc(srA, counter * sizeof(*srA));
        }
        /* For a Sparse Vector */
        else {
            mwIndex *irsV,*jcsV,*pos,old,i;
            irsV    = mxGetIr(prhs[0]);
            jcsV    = mxGetJc(prhs[0]);
            pos     = mxMalloc(nrow * sizeof(*pos));
    
            old=0;
            for (i=0; i<jcsV[1]; ++i) {
                for (; old<irsV[i]; ++old) {
                    pos[old] = -1;
                }
                pos[old] = i;
                ++old;
            }
            for (i=old; i<nrow; ++i) {
                pos[i]=-1;
            }
    
            /* Multiply */
            mwIndex k,j,tmp;
            for (j=0 ; j < nderiv ; ++j) {
                ljcs[j] = counter;
    
                for (i=jcsA[j]; i<jcsA[j+1] ; ++i) {
                    k = irsA[i];
                    if (pos[k] != (-1)) {
                        lirs[counter]= k;
                        srA[counter] = prA[i] * prV[pos[k]];
                        ++counter;
                    }
                }
            }
            ljcs[nderiv]=counter;
            
            /* Reallocate Output */
            mxRealloc(lirs, counter * sizeof(*lirs));
            mxRealloc(srA, counter * sizeof(*srA));
        }
    }
    
    /* Set Output */
    plhs[0] = mxCreateSparse(nrow,nderiv,counter,mxREAL);
    if (counter > 0) {
        mxSetIr(plhs[0],lirs);
        mxSetJc(plhs[0],ljcs);
        mxSetPr(plhs[0],srA);
    }
}
