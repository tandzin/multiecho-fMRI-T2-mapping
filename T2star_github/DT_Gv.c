#include "mex.h"
// #include <math.h>
//Jan Michalek, 6.6.2014
/*
INPUT: Gv
OUTPUT: Ua

Equivalent to the following code:
Ux = [diff(Gr,1,2), Gr(:,1) - Gr(:,i)]; 
Uy = [diff(Gr,1,1); Gr(1,:)-Gr(m,:)];

*/
#ifdef _DEBUG
#undef _DEBUG
#include <omp.h>
#define _DEBUG
#else
#include <omp.h>
#endif

//gradient components in the respective directions
#define	Gv_IN            prhs[0]
#define	nr_threads_IN    prhs[1]
#define Ua_OUT           plhs[0]
/* The gateway routine. */
void mexFunction( int nlhs, mxArray *plhs[],
								 int nrhs, const mxArray *prhs[] )
{
	double *Gv, *Ua;
	int row, col, slc, vol, delta_vol, i_vol_prev;
	int R, C, S, V;
	int i;
	const mwSize *dims;
	int nr_threads;

	/* check for the proper number of arguments */
	if(nrhs != 2)
		mexErrMsgTxt("2 input arguments are required.");
	if(nlhs != 1)
		mexErrMsgTxt("Need 1 output argument.");
	//R = mxGetM(Gr_IN);
	//C = mxGetN(Gr_IN);
	dims = mxGetDimensions(Gv_IN);
	R = dims[0];//number of rows
	C = dims[1];//number of columns
	S = dims[2];//number of slices
	V = dims[3];//number of volumes

	/* get pointers to the inputs */
	Gv = mxGetPr(Gv_IN);
	nr_threads= (int)*mxGetPr(nr_threads_IN);
	omp_set_num_threads(nr_threads);

	/* create the output */
	Ua_OUT = mxCreateNumericArray(4, dims, mxDOUBLE_CLASS, mxREAL);
	//Ua_OUT = *mxCreateDoubleMatrix(R, C, S, V, mxREAL);
	Ua = mxGetPr(Ua_OUT); 

	/* MAIN COMPUTATION */
#pragma omp parallel for \
	default(none) \
	shared (Gv,Ua,V,S,C,R)\
	private (row,col,slc,vol,delta_vol,i_vol_prev,i)

	for (vol = 1; vol <= V; vol++)
	{
		if(vol>1)
		{
			delta_vol=R*C*S;//linear distance of the previous neighbour in vol-direction,(row,col,slc,vol-1)
		}
		else
		{
			delta_vol=0;//Neumann boundary condition
			//delta_vol=(1-C)*R;//periodic boundary condition
		}
		for (slc = 1; slc <= S; slc++)
		{
			for (col = 1; col <= C; col++)
			{
				for (row = 1; row <= R; row++) 
				{
					//N= R*C*S*V – 1, indexing starts at zero
					//i corresponds to the 1st dimension in the transposed Di matrix, Di', which has dimensions [N,4]
					i=row+(col-1)*R+(slc-1)*R*C+(vol-1)*R*C*S-1;//linear index of the array element (row,col,slc,vol), must start with zero
					i_vol_prev=i-delta_vol;//linear index of the previous neighbour in vol-direction,(row,col,slc,vol-1)
					Ua[i] = Gv[i_vol_prev] - Gv[i];
				}//for (row = 1; row <= R; row++)
			}//for (col = 1; col <= C; col++)
		}//for (slc = 1; slc <= S; slc++)
	}//for (vol = 1; vol <= V; vol++)

	return;
}




