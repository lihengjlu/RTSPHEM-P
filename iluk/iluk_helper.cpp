#include "iluk_helper.h"
#include <vector>

using std::vector;

// This function changes the maximum number of nonzeros in the mxArray 
// arr to nzmax_new. It is assumed that arr has only real-valued entries.
void reallocNzmax(mxArray *arr, double** A, mwIndex** prowA, mwSize nzmax_new)
{	
	if(arr == NULL) {
		mexErrMsgTxt("mxArray is NULL");
	}
	
	if(nzmax_new < 0 || nzmax_new > mxGetN(arr)*mxGetM(arr)) {
		mexErrMsgTxt("Invalid input for new number of nonzeros");
	}
	
	void *newptr;
	newptr = mxRealloc(mxGetPr(arr), nzmax_new*sizeof(double));
	if(newptr == NULL) mexErrMsgTxt("Insufficient memory for reallocation");
	mxSetPr(arr, static_cast<double*>(newptr));
	newptr = mxRealloc(mxGetIr(arr), nzmax_new*sizeof(mwIndex));
	if(newptr == NULL) mexErrMsgTxt("Insufficient memory for reallocation");
	mxSetIr(arr, static_cast<mwIndex*>(newptr));
	mxSetNzmax(arr, nzmax_new);
	*A = mxGetPr(arr);
	*prowA = mxGetIr(arr);
}

// This function takes a sparse N x N matrix in CSC format and 
// overwrites it with its transpose.
void sparse_sqr_transpose(double *A, mwIndex *prowA, mwIndex *pcolindexA, mwSize N)
{
	vector<mwIndex> temp(N + 1, 0);
	vector<double> Atmp(pcolindexA[N]);
	vector<mwIndex> prowAtmp(pcolindexA[N]);
	vector<mwIndex> pcolindexAtmp(N + 1, 0);
	
	for(mwIndex i = 0; i != N; ++i) {
        for(mwIndex j = pcolindexA[i]; j != pcolindexA[i + 1]; ++j)
            ++pcolindexAtmp[prowA[j] + 1];
    }

	temp[1] = pcolindexAtmp[1];
	for(mwIndex i = 2; i <= N; ++i) {
        pcolindexAtmp[i] += pcolindexAtmp[i - 1];
        temp[i] = pcolindexAtmp[i];
    }
    
    for(mwIndex i = 0; i != N; i++) {
        for(mwIndex j = pcolindexA[i]; j != pcolindexA[i + 1]; ++j) {
            prowAtmp[temp[prowA[j]]] = i;
            Atmp[temp[prowA[j]]++] = A[j];
        }
    }

	// overwrite the input arrays
	pcolindexA[0] = 0;
	for(mwIndex i = 1; i <= N; ++i) {
		pcolindexA[i] = pcolindexAtmp[i];
	}
	
	for(mwIndex i = 0; i != pcolindexAtmp[N]; ++i) {
		prowA[i] = prowAtmp[i];
		A[i] = Atmp[i];
	}
}