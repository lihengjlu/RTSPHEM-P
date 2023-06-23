#include "mex.h"
#include "iluk_helper.h"
#include <algorithm>
#include <ctime>
#include <utility>
#include <vector>

#define CHECK_INPUT
#undef  CHECK_INPUT
#define ENABLE_PROFILE
#undef ENABLE_PROFILE

using std::inplace_merge;
using std::max;
using std::min;
using std::pair;
using std::upper_bound;
using std::vector;

bool comp_first_greater(mwIndex, const pair<mwIndex, int>&);
bool b_search(const vector<pair<mwIndex,int> >&, mwIndex);
/*************************************************************************************************
 * File: iluk_Cmex.c
 * Author: Killian Miller
 * Copyright (c) 2014 Killian Miller. All rights reserved.
 *
 * [L, U] = iluk_Cmex(A, lfil);
 * [L, U, levs] = iluk_Cmex(A, lfil);
 *
 * Input: 
 *         A     Sparse N by N matrix
 *         lfil  Level of fill (lfil >= 0)
 * Output:
 *         L     Sparse N by N lower triangular factor 
 *         U     Sparse N by N upper triangular factor
 *         levs  Sparse N by N matrix with levels of permitted fill-ins
 *               Note that levs is the transposed levels matrix.
 *
 * Computes an Incomplete LU factorization of the matrix A with level of
 * fill given by lfil. Note that A must be the transpose of the matrix for 
 * which you want the ILU(k) factorization. The algorithm proceeds in two 
 * phases. In the first phase locations of permitted fill-in entries are 
 * determined based on the level of fill parameter lfil and are stored in  
 * the data structure level. In the second phase the numeric values of  
 * permitted fill-in entries are computed and the L and U factors are 
 * constructed.
 *
 * Memory usage:
 *
 * L factor: nnz(L)*(sizeof(double) + sizeof(mwIndex)) + N*sizeof(mwIndex)
 * U factor: nnz(U)*(sizeof(double) + sizeof(mwIndex)) + N*sizeof(mwIndex)
 * level: N*sizeof(vector<pair<mwIndex, int> >) + (nnz(L) + nnz(U) - N)*sizeof(pair<mwIndex, int>)
 * The symbolic and numeric phases each use two temporary work arrays of length N that are
 * destroyed at the end of the respective phases.
 *************************************************************************************************/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{     
	/*
	 *  declare work variables
	 */
	int lfil;
	mwSize N, nnzU, nnzL, nnzUold, nnzLold;
	static const int realloc_thresh = 0;
	typedef vector<pair<mwIndex,int> >::iterator vitr;
	typedef vector<pair<mwIndex,int> >::const_iterator c_vitr;
	
	// lfil is the level of fill parameter
	// N is the order of the matrix A
	// realloc_thresh (>= 0) is a constant that causes reallocation of
	// of the L and U factors if nzmax(.) - nnz(.) > realloc_thresh
	
	/*
	 *  declare variables for sparse matrix A
	 */
	double *A;
	mwIndex *prowA, *pcolindexA;
    
    /*
	 *  declare variables for output
	 */
    mwIndex *prowL, *pcolindexL, *prowU, *pcolindexU;
	double *L, *U;
    
    // check input and output
	#ifdef CHECK_INPUT
	if (nrhs != 2) {
	   mexErrMsgTxt("Two input arguments required");
	}
	if (nlhs < 2 || nlhs > 3) {
	   mexErrMsgTxt("Number of outputs must be between two and three");
	}
	#endif
    
    // sparse matrix A
	#ifdef CHECK_INPUT
	if (!(mxIsSparse(prhs[0])) || mxIsComplex(prhs[0]) || mxIsEmpty(prhs[0])) {
        mexErrMsgTxt("Input matrix must be nonempty, real, and sparse");
	}
	#endif
	A = mxGetPr(prhs[0]);		   // real vector for A
	prowA = mxGetIr(prhs[0]);	   // row indices for elements of A
	pcolindexA = mxGetJc(prhs[0]); // index into the columns
	N = mxGetN(prhs[0]);           // number of columns of A
	#ifdef CHECK_INPUT
	if(N != mxGetM(prhs[0])) {
		mexErrMsgTxt("Input matrix must be square");
	}
	#endif
	
	// fill level lfil
	lfil = static_cast<int>(mxGetScalar(prhs[1]));
	++lfil; // use a shifted level of fill
	#ifdef CHECK_INPUT
	if(lfil < 0) {
		mexErrMsgTxt("Fill level must be at least 0");
	}
	#endif
	
	// level data structure: 0 <= level[i][j] <= lfil is the level of fill for 
	// permitted fill-in entries
	// initially level[i][j] = 1 for all (i,j) such that A(i,j) ~= 0
	#ifdef ENABLE_PROFILE
	clock_t t_sym = clock();
	#endif
	vector<vector<pair<mwIndex,int> > > level(N);
	level[0].reserve(pcolindexA[1]);
	for(mwIndex k = 0; k != pcolindexA[1]; ++k) {
		level[0].push_back(pair<mwIndex,int>(prowA[k], 1));
	}
	// symbolic phase
	// determine where fill-ins are permitted to occur based on the level of fill
	// count the number of nonzeros in L and U
	nnzU = pcolindexA[1];
	nnzL = 0;
	{
	mwSize nzmax = N;
	// optimize memory usage in the case of no fill-ins, i.e., ILU(0)
	if(lfil == 1) {
		nzmax = pcolindexA[1];
		for(mwIndex j = 1; j != N; ++j) {
			nzmax = max(nzmax, pcolindexA[j + 1] - pcolindexA[j]);
		}
	}
	vector<int> work_lev(N, 0); // work vector for level values
	vector<mwIndex> work_idx(nzmax); // index array for work_lev
	for(mwIndex j = 1; j != N; ++j) {
		mwSize nzmax = 0;
		// initialize level[j] work array
		for(mwIndex k = pcolindexA[j]; k != pcolindexA[j + 1]; ++k) {
			work_lev[prowA[k]] = 1;
			work_idx[nzmax++] = prowA[k];
		}
		mwIndex ii = 0;
		while(ii < nzmax && work_idx[ii] < j) {
			mwIndex i = work_idx[ii];
			++nnzL;
			// itrk->first points to the first column index in row i that is > i
			c_vitr itrk = upper_bound(level[i].begin(), level[i].end(), i, comp_first_greater);
			mwIndex nzmax_old = nzmax;
			while(itrk != level[i].end()) {
				int weight = work_lev[i] + itrk->second;			
				// is (i, itrk->first) already a fill-in? 
				if(work_lev[itrk->first] != 0) {
					work_lev[itrk->first] = min(weight, work_lev[itrk->first]);
				}
				else if(weight <= lfil) {
					work_lev[itrk->first] = weight;
					work_idx[nzmax++] = itrk->first;
				}
				++itrk;
			}
			// merge sorted sublists [0, nzmax_old) and [nzmax_old, nzmax) in work_idx
			// note that elements in positions 0,...,ii are already in their sorted
			// position
			if(nzmax > nzmax_old) {
				inplace_merge(work_idx.begin() + ii + 1, work_idx.begin() + nzmax_old, work_idx.begin() + nzmax);
			}
			++ii;
		}	
		nnzU += nzmax;
		// fill level[j] of data structure
		level[j].reserve(nzmax);
		for(mwIndex k = 0; k != nzmax; ++k) {
			level[j].push_back(pair<mwIndex,int>(work_idx[k], work_lev[work_idx[k]]));
			work_lev[work_idx[k]] = 0;
		}
	}
	} // artificial scope to force destruction of work objects
	nnzU -= nnzL;
	nnzL += N; // include unit diagonal elements
	#ifdef ENABLE_PROFILE
	t_sym = clock() - t_sym;
	printf("Time for symbolic factorization (sec): %f\n",static_cast<double>(t_sym)/CLOCKS_PER_SEC);
	#endif

	// create output matrices
	plhs[0] = mxCreateSparse(N, N, nnzL, mxREAL);
	L = mxGetPr(plhs[0]);
	prowL = mxGetIr(plhs[0]);
	pcolindexL = mxGetJc(plhs[0]);

	plhs[1] = mxCreateSparse(N, N, nnzU, mxREAL);
	U = mxGetPr(plhs[1]);
	prowU = mxGetIr(plhs[1]);
	pcolindexU = mxGetJc(plhs[1]);

	// first rows of L and U are known
	pcolindexU[0] = 0;
	pcolindexU[1] = pcolindexA[1];
	for(mwIndex j = pcolindexA[0]; j != pcolindexA[1]; ++j) {
		prowU[j] = prowA[j];
		U[j] = A[j];
	}
	// check that first nonzero in first row of U is a diagonal element
	if(prowU[0] != 0) {
		mexErrMsgTxt("U has zero diagonal element");
	}
	pcolindexL[0] = 0;
	pcolindexL[1] = 1;
	prowL[0] = 0;
	L[0] = 1.0;

	// numeric phase
	// compute the values of permitted fill-in entries and build L and U
	#ifdef ENABLE_PROFILE
	clock_t t_num = clock();
	#endif
	nnzLold = nnzL; nnzUold = nnzU;
	nnzL = 1; nnzU = pcolindexU[1];
	{
	vector<double> work(N, 0.0); // work array for rows of A
	vector<char> work_lev(N, 0); // work array for levels
	for(mwIndex j = 1; j != N; ++j) {
		// unpack jth row of A into the work array
		for(mwIndex k = pcolindexA[j]; k != pcolindexA[j + 1]; ++k) {
			work[prowA[k]] = A[k];
		}
		// unpack jth row of level into a work array
		for(c_vitr itr = level[j].begin(); itr != level[j].end(); ++itr) {
			work_lev[itr->first] = 1;
		}
		c_vitr itr = level[j].begin();
		while(itr != level[j].end() && itr->first < j) {
			double alpha = work[itr->first] / U[pcolindexU[itr->first]];
			work[itr->first] = alpha;
			for(mwIndex k = pcolindexU[itr->first] + 1; k != pcolindexU[itr->first + 1]; ++k) {
				// is prowU[k] a permitted fill-in?
				if(work_lev[prowU[k]] == 1) {
					work[prowU[k]] -= alpha * U[k];
				}
				// if memory is an issue you can avoid the work_lev array
				// by using a binary search, however, doing so may
				// significantly increase computing time, especially when
				// the level of fill is large
				// if(b_search(level[j], prowU[k])) {
				//	work[prowU[k]] -= alpha * U[k];
				// }
			}
			++itr;
		}
		// update row j of L and U and zero the work array
		itr = level[j].begin();
		while(itr != level[j].end()) {
			if(work[itr->first] != 0.0) {
				if(itr->first < j) {
					L[nnzL] = work[itr->first];
					prowL[nnzL++] = itr->first;
				}
				else {
					U[nnzU] = work[itr->first];
					prowU[nnzU++] = itr->first;
				}
				work[itr->first] = 0.0;
			}
			work_lev[itr->first] = 0;
			++itr;
		}
		// check that first nonzero in row j of U is a diagonal element
		if(prowU[pcolindexU[j]] != j) {
			mexErrMsgTxt("U has a zero diagonal element");
		}
		L[nnzL] = 1.0;
		prowL[nnzL++] = j;
		pcolindexL[j + 1] = nnzL;
		pcolindexU[j + 1] = nnzU;
	}
	} // artificial scope to force destruction of work object
	#ifdef ENABLE_PROFILE
	t_num = clock() - t_num;
	printf("Time for numeric factorization (sec): %f\n",static_cast<double>(t_num)/CLOCKS_PER_SEC);
	#endif

	// decrease the maximum number of nonzeros in L and U to nnzL and nnzU if necessary
	if(nnzLold - nnzL > realloc_thresh) {
		reallocNzmax(plhs[0], &L, &prowL, nnzL); 
	}
	if(nnzUold - nnzU > realloc_thresh) {
		reallocNzmax(plhs[1], &U, &prowU, nnzU);
	}
	
	// transpose L and U to obtain CSC format for output
	sparse_sqr_transpose(L, prowL, pcolindexL, N);
	sparse_sqr_transpose(U, prowU, pcolindexU, N);
	
	// optionally output the data structure of levels as a sparse matrix
	if(nlhs == 3) {
		plhs[2] = mxCreateSparse(N, N, nnzUold + nnzLold - N, mxREAL);
		double *levs = mxGetPr(plhs[2]);
		mwIndex *prowlevs = mxGetIr(plhs[2]);
		mwIndex *pcolindexlevs = mxGetJc(plhs[2]);
		
		pcolindexlevs[0] = 0;
		mwIndex nnzlevs = 0;
		for(int j = 0; j != N; ++j) {
			for(c_vitr itr = level[j].begin(); itr != level[j].end(); ++itr) {
				prowlevs[nnzlevs] = itr->first;
				levs[nnzlevs++] = itr->second;
			}
			pcolindexlevs[j + 1] = nnzlevs;
		}
	}
}	

// binary search that searches the vector arr for elem with comparisons
// based on the first element in each pair
bool b_search(const vector<pair<mwIndex,int> >& arr, mwIndex elem)
{
	if(arr.size() == 0)
		return false;	
	vector<pair<mwIndex,int> >::size_type imin = 0, imax = arr.size() - 1;
	while(imin <= imax) {
		vector<pair<mwIndex,int> >::size_type imid = imin + (imax - imin) / 2;
		if(arr[imid].first == elem)
			return true;
		else if(arr[imid].first < elem)
			imin = imid + 1;
		else
			imax = imid - 1;
	}
	return false;
}

// comparison function for std::upper_bound
bool comp_first_greater(mwIndex val, const pair<mwIndex, int>& p)
{
	return (val < p.first);
}