#ifndef KABSCH_RMSD
#define KABSCH_RMSD

#include <cblas.h>
#include <lapacke.h>

#include <iterator>

#include <iostream>
#include <fstream>

#include <tuple>

namespace kabsch
{

/**

	Calculates the RMSD 

	@param P - coordinates of molecule P
	@param Q - coordinates of molecule Q
	@param N - number of atoms

	@return double - the RMSD

**/
template <class M>
double rmsd(
	M *P,
	M *Q,
	const unsigned int N)
{
    double rmsd {0.0};
    const unsigned int D {3};
	const unsigned int size = N*D;

    for(unsigned int i = 0; i < size; ++i) {
        rmsd += (P[i] - Q[i])*(P[i] - Q[i]);
    }

    return sqrt(rmsd/N);
}


template <class M>
M* centroid(
	M *coordinates,
	unsigned int n_atoms)
{
    double x {0};
    double y {0};
    double z {0};
    // unsigned int size = sizeof(coordinates);
    // unsigned int n_atoms = size / 3;

	const unsigned int size = n_atoms*3;

    unsigned int i = 0;
    while(i<size)
    {
        x += coordinates[i++];
        y += coordinates[i++];
        z += coordinates[i++];
    }

    x /= n_atoms;
    y /= n_atoms;
    z /= n_atoms;

	M *centroid = new M[3];
	centroid[0] = x;
	centroid[1] = y;
	centroid[2] = z;

    return centroid;
}


template <class Matrix>
Matrix* multiply(Matrix *A, Matrix *B,
    const int M,
    const int N,
    const int K)
{
    double one = 1.0;

	// TODO Matrix C[M*N] = {0}
	// Compilers doesnt like the following syntax
	Matrix *C = new Matrix[M*N] {0};

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
        M, N, K, one,
        A, K,
        B, N, one,
        C, N);

    return C;
}


template <class Matrix>
Matrix* transpose_multiply(
	Matrix *A,
	Matrix *B,
    const int M,
    const int N,
    const int K)
{
    double one = 1.0;

    Matrix *C = new Matrix[M*N] {0};

    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
        M, N, K, one,
        A, M,
        B, N, one,
        C, N);

    return C;
}


template <class Matrix>
std::tuple<Matrix*, Matrix*, Matrix*> matrix_svd(
	Matrix *A,
	int rows,
	int cols)
{
    // lapack_int LAPACKE_dgesvd( int matrix_layout, char jobu, char jobvt,
    //     lapack_int m, lapack_int n,
    //     double* a, lapack_int lda,
    //     double* s, double* u, lapack_int ldu,
    //     double* vt, lapack_int ldvt,
    //     double* superb );

    Matrix *U = new Matrix[cols*rows] {0};
    Matrix *S = new Matrix[rows] {0};
    Matrix *VT = new Matrix[cols*rows] {0};
    Matrix *superb = new Matrix[cols*rows] {0};
    int info;

    info = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A',
        rows, cols,
        A, rows,
        S,
        U, rows,
        VT, rows,
        superb);

    if(info > 0)
    {
        // TODO What if failed
    }

	delete[] superb;

    return std::make_tuple(U, S, VT);
}


template <class M>
double determinant3x3(M A)
{
    // determinant of a square 3x3 matrix
    double det = A[0]*A[4]*A[8]
        +A[1]*A[5]*A[6]
        +A[2]*A[3]*A[7]
        -A[2]*A[4]*A[6]
        -A[1]*A[3]*A[8]
        -A[0]*A[5]*A[7];

    return det;
}


template <class M, class T>
M* kabsch(
	M *P,
	M *Q,
	const T n_atoms)
{
	// M *U = new M[3*3] {0};
	// U[0] = 1.0;
	// U[4] = 1.0;
	// U[8] = 1.0;
	// return U;

    M *U;
    M *S;
    M *V;

    M *C = transpose_multiply(P, Q, 3, 3, n_atoms);

	std::tie(U, S, V) = matrix_svd(C, 3, 3);

    // Getting the sign of the det(U)*(V) to decide whether we need to correct
    // our rotation matrix to ensure a right-handed coordinate system.
    if(determinant3x3(U)*determinant3x3(V) < 0.0)
    {
		U[3*0+2] = -U[3*0+2];
		U[3*2+2] = -U[3*1+2];
		U[3*1+2] = -U[3*2+2];
    }


    M *rotation = multiply(U, V, 3, 3, 3);

	delete[] C;
	delete[] U;
	delete[] S;
	delete[] V;

    return rotation;
}


template <class M, class T>
M* kabsch_rotate(
	M *P,
	M *Q,
	T n_atoms)
{
    M *U = kabsch(P, Q, n_atoms);
    M *product = multiply(P, U, n_atoms, 3, 3);
	delete[] U;
    return product;
}


template <class M, class T>
double kabsch_rmsd(
	M *P,
	M *Q,
	T n_atoms)
{
    M *P_rotated = kabsch_rotate(P, Q, n_atoms);
	double rmsdval = rmsd(P_rotated, Q, n_atoms);
	delete[] P_rotated;
	return rmsdval;
}

} // namespace rmsd

#endif
