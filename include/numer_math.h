#ifndef NUMERICAL_MATHEMATICS_H_
#define NUMERICAL_MATHEMATICS_H_

void print_matrix(const char* info, int m, int n, double* M, bool column_major);

void translate_to_origin_2d(int n, double* p, double* xcenter);
void scale_pts_2d(int n, double* p, double k);
void normalize_pts_2d(int n, double* p, double averMag, double* T_colmaj);
void symm_biliear_form_3d_2_rowvec(double* v, double* w, double a);

double bilinear_form_RP2(double* F, int n, double* u, double* v);

/**
RQ-decomposition based on the Gram-Schmidt process
@param M    input matrix(column-major)
@param R    output upper triangular matrix
@param Q    output orthononal matrix
@return     decomposition error. ||M-RQ||
*/
double RQ_decomp_GS_3d(double* M, double* R, double* Q);

/**
RQ-decomposition based on the Gram-Schmidt process
@param M    input matrix(column-major)
@param R    output upper triangular matrix
@param Q    output orthononal matrix
@return     decomposition error. ||M-RQ||
*/
float  RQ_decomp_GS_3d(float* M, float* R, float* Q);

/**
RQ-decomposition based on the Householder mapping
@param M    input matrix(column-major)
@param R    output upper triangular matrix
@param Q    output orthononal matrix
@return     decomposition error. ||M-RQ||
*/
double RQ_decomp_Householder_3d(double* tM, double* tR, double* tQ);

/**
RQ-decomposition based on the Householder mapping
@param M    input matrix(column-major)
@param R    output upper triangular matrix
@param Q    output orthononal matrix
@return     decomposition error. ||M-RQ||
*/
float RQ_decomp_Householder_3d(float* tM, float* tR, float* tQ);

/** Singular value decomposition in single precision
@param m    number of rows of the input matrix
@param n    number of columns of the input matrix
@param A    input matrix(column-major order)
@param sv   the pointer to the array of singular values(the array size will be min(m, n))
@param U    output matrix of left singular vectors(column-major order)
@param V    output matrix of right singular vectors(column-major order)
@param thin compute option. If true, the columns of U and V will shrink to min(m, n)
*/
void svd_sp(int m, int n, float* A, float* sv, float* U, float* V, bool thin=true);

/** Singular value decomposition in double precision
@param m    number of rows of the input matrix
@param n    number of columns of the input matrix
@param A    input matrix(column-major order)
@param sv   the pointer to the array of singular values(the array size will be min(m, n))
@param U    output matrix of left singular vectors(column-major order)
@param V    output matrix of right singular vectors(column-major order)
@param thin compute option. If true, the columns of U and V will shrink to min(m, n)
*/
void svd_dp(int m, int n, double* A, double* sv, double* U, double* V, bool thin=true);

/** Singular value decomposition in double precision
@param A    The pointer of the input matrix(Eigen::MatrixXd)
@param sv   the pointer to the array of singular values(Eigen::MatrixXd). 
@param U    output matrix of left singular vectors((Eigen::MatrixXd)
@param V    output matrix of right singular vectors((Eigen::MatrixXd)
@param thin compute option. If true, the columns of U and V will shrink to min(m, n)
*/
void svd_dp(void* A, void* sv, void* U, void* V, bool thin=true);

/** Find a nonzero kernel vector v such that Av = 0.

This uses the singular value decomposition. 
If A is under-determined, tha is, the columns is larger than rows, 
the kernel vector becomes the last column vector of the right singular vector matrix V and
the return value is zero.
In the other cases, the return value is the smallest singular value which may be nonzero, 
and the minimizer becomes the corresponding singular vector.

@param m   number of rows of the given matrix
@param n   number of columns of the given matrix
@param A   given matrix(column-major order)
@param v   pointer to the vector minimizing Ax
@return    the minimum value of Ax
*/
double find_nonzero_kernel_vector(int m, int n, double* A, double* v);

void quaternionToRotation(double* quat, double* R);
void compute_covariance_matrix(int n, double* p, double* center, double* C);

#endif