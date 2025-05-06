#include "mex.h"
// #include <math.h>

/*
INPUT: U
OUTPUT: Ux, Uy
// Jan Michalek, 6.6.2014
Equivalent to the following code:
Ux = [diff(U,1,2), U(:,1000) - U(:,1000)];
Uy = [diff(U,1,1); U(1000,:)-U(1000,:)];
*/
#ifdef _DEBUG
#undef _DEBUG
#include <omp.h>
#define _DEBUG
#else
#include <omp.h>
#endif

#define	U_IN            prhs[0]
#define	nr_threads_IN   prhs[1]
#define Dv_OUT          plhs[0]
/* The gateway routine. */
void mexFunction( int nlhs, mxArray *plhs[],
								 int nrhs, const mxArray *prhs[] )
{
	double *U, *Dv;
	int row, col, slc, vol, delta_vol;
	int i_vol_next,i;

	int R, C, S, V;
	const mwSize *dims;
	int nr_threads;

	/* check for the proper number of arguments */
	if(nrhs != 2)
		mexErrMsgTxt("2 input arguments are required.");
	if(nlhs != 1)
		mexErrMsgTxt("Need 1 output argument.");

	dims = mxGetDimensions(U_IN);
	R = dims[0];//number of rows
	C = dims[1];//number of columns
	S = dims[2];//number of slices
	V = dims[3];//number of volumes

	/* get pointer to the input */
	U = mxGetPr(U_IN);
	//Ui = mxGetPi(U_IN);
	nr_threads= (int)*mxGetPr(nr_threads_IN);
	omp_set_num_threads(nr_threads);

	/* create Dv */
	//Dv_OUT = mxCreateDoubleMatrix(R, C, mxREAL);
	Dv_OUT = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);
	Dv = mxGetPr(Dv_OUT);  // Vr = mxGetPr(plhs[2]);

	/* MAIN COMPUTATION */
#pragma omp parallel for \
	default(none) \
	shared (U,Dv,R,C,S,V)\
	private (row,col,slc,vol,delta_vol,i_vol_next,i)
	for (vol = 1; vol <= V; vol++)
	{
		if(vol<V)
		{
			delta_vol=R*C*S;//linear distance of the next neighbour in vol-direction,(row,col,slc,vol+1)
		}
		else
		{//Neumann boundary condition
			delta_vol=0;
		}
		for (slc = 1; slc <= S; slc++)
		{
			for (col = 1; col <= C; col++)
			{
				for (row = 1; row <= R; row++) 
				{
					i=row+(col-1)*R+(slc-1)*R*C+(vol-1)*R*C*S-1;//linear index of the array element (row,col,slc,vol), must start with zero
					i_vol_next=i+delta_vol;//linear index of the next neighbour in vol-direction,(row,col,slc,vol+1)
					Dv[i] = U[i_vol_next] - U[i];
				}//for (row = 1; row <= R; row++)
			}//for (col = 1; col <= C; col++)
		}//for (slc = 1; slc <= S; slc++)
	}//for (vol = 1; vol <= V; vol++)

	return;
}




