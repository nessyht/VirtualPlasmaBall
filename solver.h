#include "misc.h"
#include "cublas_wrapper.h"
#include "tools/globals.h"

/*
 *      solver.h
 *
 *      This file contains various iterative solvers for the given problem.
 *
 *      @author Simon Schoelly
 *
 */

template<class FT>
__global__ void d2d(const int m, FT const * const x, FT const * const mask, FT const dxyz, FT * const xprime, FT * const yprime)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	if (mask[tid] == GAS)
	{
		xprime[tid] = (x[tid - 1] - x[tid + 1])/dxyz;
		yprime[tid] = (x[tid - m] - x[tid + m])/dxyz;
	}
}

template<class FT>

__global__ void d3d(const int m, FT const * const x, FT const * const mask, FT const alpha, FT * const xprime, FT * const yprime, FT * const zprime)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;

	int mm = m * m;

	if (mask[tid] == GAS)
	{
		xprime[tid] = x[tid - 1] - x[tid + 1];
		yprime[tid] = x[tid - m] - x[tid + m];
		zprime[tid] = x[tid - mm] - x[tid + mm];
	}
}

 /*
  *      multiplies the vector x with the matrix A for the 2D problem
  *
  *      @param FT Field Type - Either float or double
  *      @param m >= 1 - grid size m
  *      @param alpha constant alpha > 0
  *      @param x != NULL - input vector x of length m*m
  *      @param b != NULL - output vector of length m*m
  *
  *      @return A*x
  */
template<class FT>
__global__ void multiply_by_A(int const m, FT const alpha, FT const dxyz2, FT const * const mask, FT const * const x, FT * const b) {
	int n = m * m;
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= n) {
		return;
	}
	//int neu = 0;
	//double neum = 0;
	FT value = (FT(4))*x[tid];
	//if (tid % m != 0) {
	//	if (mask[tid - 1] == BORDER)
	//	{
	//		neum += x[tid - 1];
	//		++neum;
	//	}
	//	value -= x[tid - 1];
	//}
	//if ((tid + 1) % m != 0) {
	//	if (mask[tid + 1] == BORDER)
	//	{
	//		neum += x[tid + 1];
	//		++neum;
	//	}
	//	
	//	value -= x[tid + 1];
	//}
	//if (tid + m < n) {
	//	if (mask[tid - m] == BORDER)
	//	{
	//		neum += x[tid - m];
	//		++neum;
	//	}
	//	
	//	value -= x[tid - m];
	//}
	//if (tid - m >= 0) {
	//	if (mask[tid + m] == BORDER)
	//	{
	//		neum += x[tid + m];
	//		++neum;
	//	}
	//	
	//	value -= x[tid + m];
	//}
	//b[tid] = value;
	////if (neu > 0) 
	////	b[tid] = neum + FT(neu)*x[tid];

	if (mask[tid] == ELECTRODE || mask[tid] == PLASMA || mask[tid] == FAKEELECTRODE 
		|| mask[tid] == BORDER || mask[tid] == CONDUCTOR)
	{
		b[tid] = x[tid] / dxyz2;
		return;
	}

	//// neumann boundary
	////FT off_diags = 0;
	////int bound = 0;
	////if (mask[tid - 1] == BORDER)
	////{
	////	++bound;
	////	off_diags += x[tid + 1];
	////}
	////if (mask[tid + 1] == BORDER)
	////{
	////	++bound;
	////	off_diags += x[tid - 1];
	////}
	////if (mask[tid + m] == BORDER)
	////{
	////	++bound;
	////	off_diags += x[tid - m];
	////}
	////if (mask[tid - m] == BORDER)
	////{
	////	++bound;
	////	off_diags += x[tid + m];
	////}
	////if (bound > 0)
	////{
	////	value = x[tid] * bound;
	////	value -= off_diags;
	////	return;
	////}
	//// 
	value -= x[tid - 1];
	value -= x[tid + 1];
	value -= x[tid + m];
	value -= x[tid - m];

	b[tid] = value / dxyz2;



	//FT current = b[tid];
	//FT value = -(FT(4))*x[tid];

	//if (mask[tid] == GAS || mask[tid] == AIR)
	//{
	//	value += x[tid - 1];
	//	value += x[tid + 1];
	//	value += x[tid - m];
	//	value += x[tid + m];
	//}
	//else
	//	value = x[tid];

	//b[tid] = value;


	//if (tid % m != 0)
	//{
	//	switch ((int)mask[tid - 1])
	//	{
	//	case ELECTRODE:
	//		value -= -10000;
	//		break;

	//	case PLASMA:
	//		value -= -10000;
	//		break;

	//	case BORDER:
	//		value -= 0;
	//		break;

	//	default:
	//		value -= x[tid - 1];
	//	}
	//	//value -= x[tid - 1];
	//}
	//if ((tid + 1) % m != 0) 
	//{
	//	switch ((int)mask[tid + 1])
	//	{
	//	case ELECTRODE:
	//		value -= -10000;
	//		break;

	//	case PLASMA:
	//		value -= -10000;
	//		break;

	//	case BORDER:
	//		value -= 0;
	//		break;

	//	default:
	//		value -= x[tid + 1];

	//	}
	//	//value -= x[tid + 1];
	//}
	//if (tid + m < n) 
	//{
	//	switch ((int)mask[tid + m])
	//	{
	//	case ELECTRODE:
	//		value -= -10000;
	//		break;

	//	case PLASMA:
	//		value -= -10000;
	//		break;

	//	case BORDER:
	//		value -= 0;
	//		break;

	//	default:
	//		value -= x[tid + m];
	//		break;

	//	}
	//	//value -= x[tid + m];
	//}
	//if (tid - m >= 0) 
	//{
	//	switch ((int)mask[tid - m])
	//	{
	//	case ELECTRODE:
	//		value -= -10000;
	//		break;
	//	case PLASMA:
	//		value -= -10000;
	//		break;

	//	case BORDER:
	//		value -= 0;
	//		break;

	//	default:
	//		value -= x[tid - m];
	//	}
	//	//value -= x[tid - m];
	//}
	
	//if (mask[tid] == ELECTRODE) // || mask[tid] == PLASMA)
	//	b[tid] = x[tid];	
	////	
	//if (mask[tid] == BORDER || mask[tid] == GLASS)
	//	b[tid] = x[tid];
	////if (mask[tid] == PLASMA)
	//	b[tid] = 0; // out: 0, x[tid], -x[tid]

}

/*
 *      multiplies the vector x with the matrix A for the 3D problem
 *
 *      @param FT Field Type - Either float or double
 *      @param m >= 1 - grid size m
 *      @param alpha constant alpha > 0
 *      @param x != NULL - input vector x of length m*m*m
 *      @param b != NULL - output vector of length m*m*m
 *
 *      @return A*x
 */
template<class FT>
__global__ void multiply_by_A3D(int const m, FT const alpha, FT const * const mask, FT const * const x, FT * const b) {
	int n = m * m * m;
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= n) {
		return;
	}

	FT value = (alpha + FT(6))*x[tid];

	if (tid + (m*m) < n) {
		value -= x[tid + m * m];
	}
	if (tid - (m*m) >= 0) {
		value -= x[tid - m * m];
	}

	if (tid % m != 0) {
		value -= x[tid - 1];
	}
	if ((tid + 1) % m != 0) {
		value -= x[tid + 1];
	}

	if (tid / (m*m) == (tid + m) / (m*m)) {
		value -= x[tid + m];
	}

	if ((tid + m * m) / (m*m) == (tid + m * m - m) / (m*m)) {
		value -= x[tid - m];
	}

	b[tid] = value;

	if (mask[tid] == ELECTRODE)
		b[tid] = x[tid] - (1.0 / alpha);


	if (x[tid] == 0)
		b[tid] = x[tid];
}

/*
*      chebyshev iteration solver for the 2D problem
*
*      @param FT Field Type - Either float or double
*      @param cublas_handle must be initalized with cublasCreate
*      @param cusparse_handle must be initalized with cusparseCreate
*      @param 4096 >= m >= 1 - grid size m
*      @param alpha constant alpha > 0
*      @param b != NULL - input vector x of length m*m
*      @param x != NULL - output vector of length m*m
*      @param max_iter >= 0 the maximum number ot iterations before the algorithm halts
*      @param tolerance >= 0 the algorithm halts when the norm of the residum has shrunk more than tolerance
*      @param preconditiner != NUll - a preconditioner for A
*
*      @return A\b
*      @return number of iterations until convergence
*/

template<class FT>
int solve_with_chebyshev_iteration(cublasHandle_t const             cublas_handle,
	cusparseHandle_t const           cusparse_handle,
	int const                        m,
	FT const                         alpha,
	FT const * const                 b,
	FT * const                       x,
	int const                        max_iter,
	FT const                         tolerance,
	Preconditioner<FT> * const       preconditioner)
{

	FT const FT_ONE(1);
	FT const FT_MINUS_ONE(-1);

	int n = m * m;

	FT lambda_max = FT(1);
	FT lambda_min = (FT(8) + alpha) / (FT(8) + alpha + FT(16) / alpha);



	FT d = (lambda_max + lambda_min) / FT(2);
	FT c = (lambda_max - lambda_min) / FT(2);

	FT *Ax, *Ax2, *r, *z, *p;
	cudaMalloc((void **)&Ax, n * sizeof(FT));
	cudaMalloc((void **)&Ax2, n * sizeof(FT));
	cudaMalloc((void **)&r, n * sizeof(FT));
	cudaMalloc((void **)&z, n * sizeof(FT));
	cudaMalloc((void **)&p, n * sizeof(FT));

	if (preconditioner != NULL) {
		preconditioner->init(m, alpha, cublas_handle, cusparse_handle);
	}


	// x = 0
	device_memset<FT>(x, FT(0), n);

	// Ax = A*x
	multiply_by_A<FT> << <divide_and_round_up(n, 1024), 1024 >> > (m, alpha, x, Ax);

	// r = b
	cublas_copy(cublas_handle, n, b, r);
	// r = r - Ax = b - Ax
	cublas_axpy(cublas_handle, n, &FT_MINUS_ONE, Ax, r);

	FT norm0;
	cublas_nrm2(cublas_handle, n, r, &norm0);

	int num_iter;
	for (num_iter = 1; num_iter <= max_iter; ++num_iter) {
		FT fa, fb;
		// z = M \ r
		preconditioner->run(r, z);
		if (num_iter == 1) {
			// p = z
			cublas_copy(cublas_handle, n, z, p);

			fa = FT(1) / d;
		}
		else {
			fb = (c*fa / FT(2));
			fb = fb * fb;
			fa = FT(1) / (d - fb / fa);

			// p = fb * p
			cublas_scal(cublas_handle, n, &fb, p);

			// p = p + z
			cublas_axpy(cublas_handle, n, &FT_ONE, z, p);
		}

		// x = x + fa*p
		cublas_axpy(cublas_handle, n, &fa, p, x);

		// Ax = A*x
		multiply_by_A<FT> << <divide_and_round_up(n, 1024), 1024 >> > (m, alpha, x, Ax);

		// r = b
		cublas_copy(cublas_handle, n, b, r);
		// r = r - Ax = b - Ax
		cublas_axpy(cublas_handle, n, &FT_MINUS_ONE, Ax, r);

		FT norm;
		cublas_nrm2(cublas_handle, n, r, &norm);
		if (norm <= tolerance * norm0) {
			break;
		}
	}

	return num_iter;
}



/*
*      pcg solver for the 2D problem
*
*      @param FT Field Type - Either float or double
*      @param cublas_handle must be initalized with cublasCreate
*      @param cusparse_handle must be initalized with cusparseCreate
*      @param 4096 >= m >= 1 - grid size m
*      @param alpha constant alpha > 0
*      @param b != NULL - input vector x of length m*m
*      @param x != NULL - output vector of length m*m
*      @param max_iter >= 0 the maximum number ot iterations before the algorithm halts
*      @param tolerance >= 0 the algorithm halts when the norm of the residum has shrunk more than tolerance
*      @param preconditiner if NULL then then cg algorithm is used
*
*      @return A\b
*      @return number of iterations until convergence
*/
template<class FT>
int solve_with_conjugate_gradient(cublasHandle_t   const cublas_handle,
	cusparseHandle_t const cusparse_handle,
	int const				m,
	FT const				alpha,
	FT const				dxyz,
	FT const * const		b,
	FT const * const		mask,
	FT * const				x,
	FT * const				xprime,
	FT * const				yprime,
	int						max_iter,
	FT						tolerance,
	Preconditioner<FT> *	preconditioner = NULL)
{
	int const n = m * m;
	float const dxyz2 = dxyz * dxyz;

	if (preconditioner != NULL) {
		preconditioner->init(m, alpha, cublas_handle, cusparse_handle);
	}

	device_memset<FT>(x, FT(0), n); // x = 0
	FT *p, *r, *h, *Ax, *q;
	cudaMalloc((void **)&p, n * sizeof(FT));
	cudaMalloc((void **)&r, n * sizeof(FT));
	cudaMalloc((void **)&h, n * sizeof(FT));
	cudaMalloc((void **)&Ax, n * sizeof(FT));
	cudaMalloc((void **)&q, n * sizeof(FT));
	FT beta, c, ph;

	FT const FT_ONE(1);
	FT const FT_MINUS_ONE(-1);

	// Ax = A*x
	multiply_by_A<FT> << <divide_and_round_up(n, 1024), 1024 >> > (m, alpha, dxyz2, mask, x, Ax);
	// r = b
	cublas_copy(cublas_handle, n, b, r);
	// r = r - Ax = b - Ax
	cublas_axpy(cublas_handle, n, &FT_MINUS_ONE, Ax, r);

	FT norm0;
	cublas_nrm2(cublas_handle, n, r, &norm0);

	if (preconditioner == NULL) {
		// p = r
		cublas_copy(cublas_handle, n, r, p);
	}
	else {
		// p = M \ r
		preconditioner->run(r, p);
		// q = p
		cublas_copy(cublas_handle, n, p, q);
	}
	int iter_num;
	for (iter_num = 1; iter_num <= max_iter; ++iter_num) {
		if (preconditioner == NULL) {
			// beta = <r, r>
			cublas_dot(cublas_handle, n, r, r, &beta);
		}
		else {
			// beta = <r, q>
			cublas_dot(cublas_handle, n, r, q, &beta);
		}

		// h = Ap
		multiply_by_A<FT> << <divide_and_round_up(n, 1024), 1024 >> > (m, alpha, dxyz2, mask, p, h);

		// ph = <p, h>
		cublas_dot(cublas_handle, n, p, h, &ph);

		c = beta / ph;

		//  x = x + c*p
		cublas_axpy(cublas_handle, n, &c, p, x);

		// r = r - c*h
		FT minus_c = -c;
		cublas_axpy(cublas_handle, n, &minus_c, h, r);

		FT norm;
		cublas_nrm2(cublas_handle, n, r, &norm);
		if (norm <= tolerance * norm0) {
			break;
		}

		if (preconditioner != NULL) {
			// q = B \ r
			preconditioner->run(r, q);
		}

		if (preconditioner == NULL) {
			// rr = <r, r>
			FT rr;
			cublas_dot(cublas_handle, n, r, r, &rr);

			beta = rr / beta;
		}
		else {
			// rq = <r, q>
			FT rq;
			cublas_dot(cublas_handle, n, r, q, &rq);

			beta = rq / beta;
		}

		// p = beta * p
		cublas_scal(cublas_handle, n, &beta, p);

		if (preconditioner == NULL) {
			// p = r + p
			cublas_axpy(cublas_handle, n, &FT_ONE, r, p);
		}
		else {
			// p = q + p
			cublas_axpy(cublas_handle, n, &FT_ONE, q, p);
		}
	}

	cudaFree(p);
	cudaFree(r);
	cudaFree(h);
	cudaFree(Ax);
	cudaFree(q);

	d2d << <divide_and_round_up(n, 1024), 1024 >> > (m, x, mask, dxyz, xprime, yprime);

	return iter_num;
}

/*
*      pcg solver for the 3D problem
*
*      @param FT Field Type - Either float or double
*      @param cublas_handle must be initalized with cublasCreate
*      @param cusparse_handle must be initalized with cusparseCreate
*      @param 256 >= m >= 1 - grid size m
*      @param alpha constant alpha > 0
*      @param b != NULL - input vector x of length m*m*m
*      @param x != NULL - output vector of length m*m*m
*      @param max_iter >= 0 the maximum number ot iterations before the algorithm halts
*      @param tolerance >= 0 the algorithm halts when the norm of the residum has shrunk more than tolerance
*      @param preconditiner if NULL then then cg algorithm is used
*
*      @return A\b
*      @return number of iterations until convergence
*/
template<class FT>
int solve_with_conjugate_gradient3D(cublasHandle_t   const cublas_handle,
	cusparseHandle_t const cusparse_handle,
	int const				m,
	FT const				alpha,
	FT const * const		b,
	FT const * const		mask,
	FT * const				x,
	FT * const				xprime,
	FT * const				yprime,
	FT * const				zprime,
	int						max_iter,
	FT						tolerance,
	Preconditioner<FT> *	preconditioner = NULL)
{
	int const n = m * m*m;

	if (preconditioner != NULL) {
		preconditioner->init(m, alpha, cublas_handle, cusparse_handle);
	}

	device_memset<FT>(x, FT(0), n); // x = 0
	FT *p, *r, *h, *Ax, *q;
	cudaMalloc((void **)&p, n * sizeof(FT));
	cudaMalloc((void **)&r, n * sizeof(FT));
	cudaMalloc((void **)&h, n * sizeof(FT));
	cudaMalloc((void **)&Ax, n * sizeof(FT));
	cudaMalloc((void **)&q, n * sizeof(FT));
	FT beta, c, ph;

	FT const FT_ONE(1);
	FT const FT_MINUS_ONE(-1);

	// Ax = A*x
	multiply_by_A3D<FT> << <divide_and_round_up(n, 1024), 1024 >> > (m, alpha, mask, x, Ax);
	// r = b
	cublas_copy(cublas_handle, n, b, r);
	// r = r - Ax = b - Ax
	cublas_axpy(cublas_handle, n, &FT_MINUS_ONE, Ax, r);

	FT norm0;
	cublas_nrm2(cublas_handle, n, r, &norm0);

	if (preconditioner == NULL) {
		// p = r
		cublas_copy(cublas_handle, n, r, p);
	}
	else {
		// p = M \ r
		preconditioner->run(r, p);
		// q = p
		cublas_copy(cublas_handle, n, p, q);
	}
	int iter_num;
	for (iter_num = 1; iter_num <= max_iter; ++iter_num) {
		if (preconditioner == NULL) {
			// beta = <r, r>
			cublas_dot(cublas_handle, n, r, r, &beta);
		}
		else {
			// beta = <r, q>
			cublas_dot(cublas_handle, n, r, q, &beta);
		}

		// h = Ap
		multiply_by_A3D<FT> << <divide_and_round_up(n, 1024), 1024 >> > (m, alpha, mask, p, h);

		// ph = <p, h>
		cublas_dot(cublas_handle, n, p, h, &ph);

		c = beta / ph;

		//  x = x + c*p
		cublas_axpy(cublas_handle, n, &c, p, x);

		// r = r - c*h
		FT minus_c = -c;
		cublas_axpy(cublas_handle, n, &minus_c, h, r);

		FT norm;
		cublas_nrm2(cublas_handle, n, r, &norm);
		if (norm <= tolerance * norm0) {
			break;
		}

		if (preconditioner != NULL) {
			// q = B \ r
			preconditioner->run(r, q);
		}

		if (preconditioner == NULL) {
			// rr = <r, r>
			FT rr;
			cublas_dot(cublas_handle, n, r, r, &rr);

			beta = rr / beta;
		}
		else {
			// rq = <r, q>
			FT rq;
			cublas_dot(cublas_handle, n, r, q, &rq);

			beta = rq / beta;
		}

		// p = beta * p
		cublas_scal(cublas_handle, n, &beta, p);

		if (preconditioner == NULL) {
			// p = r + p
			cublas_axpy(cublas_handle, n, &FT_ONE, r, p);
		}
		else {
			// p = q + p
			cublas_axpy(cublas_handle, n, &FT_ONE, q, p);
		}
	}

	cudaFree(p);
	cudaFree(r);
	cudaFree(h);
	cudaFree(Ax);
	cudaFree(q);

	d3d << <divide_and_round_up(n, 1024), 1024 >> > (m, x, mask, alpha, xprime, yprime, zprime);


	return iter_num;
}