/* This fuMatrix is saved in column-major order.
 * Instead of sorting Ir array after computation, I believe up a supplement
 * array, and rearrange Ir array with one pass O(nrow), but this involves
 * allocating an array of size nrow. This should be fast if (A_deriv*b_val)
 * is not too sparse, but if it is really sparse, it might be faster to
 * sort instead, which will give approximately
 * O(nderiv * nnz_per_column * log(nnz_per_column)). */

#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{

    /* Read In Matrix A */
    double *prA;
	double *srA;
	double *prV,tmpval;
	double old_deriv;
	mwSize mA, nderiv, nnz;
    mwIndex *irsA,*jcsA;
	mwIndex *lirs,*ljcs;
	mwSize ncol,counter=0,nrow;
	mwIndex *pointer, rownnz, pointed,tmpind;
	mwIndex k,j,i;
	/* mwIndex *irsV,*jcsV,*pos,old,i; */
	mwIndex *irsV,*jcsV,*pos,old;
	/* mwIndex k,j,tmp; */
	mwIndex tmp;
    mA      = mxGetM(prhs[0]);
    nderiv  = mxGetN(prhs[0]);
    irsA    = mxGetIr(prhs[0]);
    jcsA    = mxGetJc(prhs[0]);
    prA     = mxGetPr(prhs[0]);
    nnz     = jcsA[nderiv];
    
    /* Allocate Data Matrix for Output*/
	/* double *srA; */
    /* mwIndex *lirs,*ljcs; */
    lirs    = mxMalloc( nnz * sizeof(*lirs));
    ljcs    = mxMalloc( (nderiv+1) * sizeof(*ljcs));
    srA     = mxMalloc( nnz * sizeof(*srA));
    
    /* Read in Vector V */
    /* mwSize ncol,counter=0,nrow; */
    /* mwIndex *pointer, rownnz, pointed,tmpind; */
    /* double *prV,tmpval */;
    ncol    = mxGetM(prhs[1]);
    prV     = mxGetPr(prhs[1]);
    nrow    = mA/ncol;
    pointer = mxMalloc( nrow * sizeof(*pointer));
    
    /* For a full vector */
    if (!mxIsSparse(prhs[1])) {
        /* mwIndex k,j,i, */ 
		old_deriv=nderiv+1;
        for (j=0 ; j < nderiv ; ++j) {
            ljcs[j] = counter;
            
            /* Reset pointer Array */
            rownnz=0;
            for (i = 0; i<nrow;++i) {
                pointer[i]=-1;
            }
            
            for (i=jcsA[j]; i<jcsA[j+1] ; ++i) {
                k = irsA[i];
                pointed = pointer[k%nrow];
                if (prV[k / nrow]!=0.0) {
                    if (pointed != -1) {
                        srA[pointed] += prA[i] * prV[k / nrow];
                    }
                    else {
                        pointer[k%nrow] = counter+rownnz;
                        srA[counter+rownnz] = prA[i] * prV[k / nrow];
                        lirs[counter+rownnz]= k%nrow;
                        ++rownnz;
                    }
                }
            }
            
            /* Post Sort Ir Array */
            for (i=0;i<nrow;++i){
                pointed=pointer[i];
                if (pointed!=-1) {
                    tmpind=lirs[counter];
                    tmpval=srA[counter];
                    lirs[counter]=i;
                    srA[counter]=srA[pointed];
                    srA[pointed]=tmpval;
                    lirs[pointed]=tmpind;
                    pointer[tmpind]=pointed;
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
        /* mwIndex *irsV,*jcsV,*pos,old,i; */
        irsV    = mxGetIr(prhs[1]);
        jcsV    = mxGetJc(prhs[1]);
        pos     = mxMalloc(ncol * sizeof(*pos));
    
        old=0;
        for (i=0; i<jcsV[1]; ++i) {
            for (; old<irsV[i]; ++old) {
                pos[old] = -1;
            }
            pos[old] = i;
            ++old;
        }
        for (i=old; i<ncol; ++i) {
            pos[i]=-1;
        }
       
        /* Multiply */
        /* mwIndex k,j,tmp */;
        for (j=0 ; j < nderiv ; ++j) {
            ljcs[j] = counter;
            
            /* Reset pointer Array */
            rownnz=0;
            for (i = 0; i<nrow; ++i) {
                pointer[i]=-1;
            }
            
            for (i=jcsA[j]; i<jcsA[j+1] ; ++i) {
                k = irsA[i];
                tmp = k / nrow;
                if (pos[tmp] != (-1)) {
                    pointed = pointer[k%nrow];
                    if (pointed != -1) {
                        srA[pointed] += prA[i] * prV[pos[tmp]];
                        /* mexPrintf("%f\n",srA[pointed]); */
                    }
                    else {
                        pointer[k%nrow] = counter+rownnz;
                        srA[counter+rownnz] = prA[i] * prV[pos[tmp]];
                        lirs[counter+rownnz]= k%nrow;
                        ++rownnz;
                        /* mexPrintf("%f\n",srA[counter+rownnz-1]); */
                    }
                }
            }
                        
            /* Post Sort Ir Array */
            for (i=0;i<nrow;++i){
                pointed=pointer[i];
                if (pointed!=-1) {
                    tmpind=lirs[counter];
                    tmpval=srA[counter];
                    lirs[counter]=i;
                    srA[counter]=srA[pointed];
                    srA[pointed]=tmpval;
                    lirs[pointed]=tmpind;
                    pointer[tmpind]=pointed;
                    ++counter;
                    /* mexPrintf("%d\n",counter); */
                }
            }
        }
        ljcs[nderiv]=counter;
        /* mexPrintf("%d here\n",counter); */
    
        mxRealloc(lirs, counter * sizeof(*lirs));
        mxRealloc(srA, counter * sizeof(*srA));
    }
    
    /* Set Output */
    plhs[0] = mxCreateSparse(nrow,nderiv,counter,mxREAL);
	if (counter!=0) {
		mxSetIr(plhs[0],lirs);
		mxSetJc(plhs[0],ljcs);
		mxSetPr(plhs[0],srA);
	}
    mxFree(pointer);
    
}