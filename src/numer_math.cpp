#include "numer_math.h"
#include <iostream>
#include <vector>

#include <Eigen/Dense>

void print_matrix(const char* info, int m, int n, double* M, bool column_major)
{
    std::cout << info << " : " << std::endl;

    if (column_major)
    {
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
            {

                std::cout << M[m*j + i] << " ";
            }
            std::cout << std::endl;
        }
    }
    else
    {
        for (int i = 0; i < m; ++i)
        {
            for (int j = 0; j < n; ++j)
            {
                std::cout << M[n*i + j] << " ";
            }
            std::cout << std::endl;
        }
    }
}

void translate_to_origin_2d(int n, double* p, double* xcenter)
{
    xcenter[0] = xcenter[1] = 0.0;

    for (int i = 0; i < n; ++i)
    {
        xcenter[0] += p[2 * i + 0];
        xcenter[1] += p[2 * i + 1];
    }

    if (n > 0)
    {
        xcenter[0] /= (float)n;
        xcenter[1] /= (float)n;
    }

    for (int i = 0; i < n; ++i)
    {
        p[2 * i + 0] -= xcenter[0];
        p[2 * i + 1] -= xcenter[1];
    }
}

double average_magnitude_2d(int n, double*p)
{
    double mag = 0.0;
    for (int i = 0; i < n; ++i)
    {
        Eigen::Vector2d tp(p[2 * i + 0], p[2 * i + 1]);

        mag += tp.norm();
    }

    if (n > 0)
    {
        mag /= (double)n;
    }

    return mag;
}

void scale_pts_2d(int n, double* p, double k)
{
    for (int i = 0; i < n; ++i)
    {
        p[2 * i + 0] *= k;
        p[2 * i + 1] *= k;
    }
}

void normalize_pts_2d(int n, double* p, double averMag, double* T_colmaj)
{
    Eigen::Vector2d xcenter;
    translate_to_origin_2d(n, p, xcenter.data());
    double k = averMag / average_magnitude_2d(n, p);
    scale_pts_2d(n, p, k);

    Eigen::Map<Eigen::Matrix3d> T(T_colmaj);
    T.setIdentity();
    T.block<2, 2>(0, 0) *= k;
    T(0, 2) = -k*xcenter(0);
    T(1, 2) = -k*xcenter(1);
}

void symm_biliear_form_3d_2_rowvec(double* v, double* w, double* a)
{
    a[0] = v[0]*w[0];
    a[1] = v[0]*w[1] + v[1]*w[0];
    a[2] = v[0]*w[2] + v[2]*w[0];
    a[3] = v[1]*w[1];   
    a[4] = v[1]*w[2] + v[2]*w[1];
    a[5] = v[2]*w[2];
}

double bilinear_form_RP2(double* tF, int n, double* tu, double* tv)
{
    Eigen::Map<Eigen::Matrix3d> F(tF);

    double result = 0.0;

    for (int i = 0; i < n; ++i)
    {
        Eigen::Vector3d u(tu[2*i+0], tu[2*i+1], 1.0);
        Eigen::Vector3d v(tv[2*i+0], tv[2*i+1], 1.0);

        result += u.dot(F*v);
    }

    return result;
}

Eigen::Vector3d householder(Eigen::Vector3d v, Eigen::Vector3d x)
{
    return x - 2.0*v.dot(x) / (v.dot(v))* v;
}

Eigen::Vector2d householder(Eigen::Vector2d v, Eigen::Vector2d x)
{
    return x - 2.0*v.dot(x) / (v.dot(v))* v;
}

double find_householder_matrix_3(Eigen::Vector3d& x, Eigen::Vector3d& e, Eigen::Matrix3d& tQ)
{
    double alpha = x.norm();
    //std::cout << "alpha=" << alpha << std::endl;

    double sign = 1.0;
    double xe = x.dot(e);
    if (xe < 0.0)   sign = -1.0;
    //std::cout << "sign=" << sign << std::endl;

    Eigen::Vector3d te[3];
    te[0] = Eigen::Vector3d(1.0, 0.0, 0.0);
    te[1] = Eigen::Vector3d(0.0, 1.0, 0.0);
    te[2] = Eigen::Vector3d(0.0, 0.0, 1.0);

    Eigen::Vector3d v = x - sign*alpha*e;
    //std::cout << "v=" << v << std::endl;

    tQ.block<3, 1>(0, 0) = householder(v, te[0]);
    tQ.block<3, 1>(0, 1) = householder(v, te[1]);
    tQ.block<3, 1>(0, 2) = householder(v, te[2]);

    return sign*alpha; // NOTE : this return value may not be necessary for later steps.
}

double find_householder_matrix_2(Eigen::Vector2d& x, Eigen::Vector2d& e, Eigen::Matrix2d& tQ)
{
    double alpha = x.norm();

    double sign = 1.0;
    double xe = x.dot(e);
    if (xe < 0.0)   sign = -1.0;

    Eigen::Vector2d te[2];
    te[0] = Eigen::Vector2d(1.0, 0.0);
    te[1] = Eigen::Vector2d(0.0, 1.0);


    Eigen::Vector2d v = x - sign*alpha*e;

    tQ.block<2, 1>(0, 0) = householder(v, te[0]);
    tQ.block<2, 1>(0, 1) = householder(v, te[1]);

    return sign*alpha; // NOTE : this return value may not be necessary for later steps.
}

double RQ_decomp_Householder_3d(double* tM, double* tR, double* tQ)
{
    Eigen::Map<Eigen::Matrix3d> M(tM);
    Eigen::Map<Eigen::Matrix3d> R(tR);
    Eigen::Map<Eigen::Matrix3d> Q(tQ);
    //std::cout << "M=\n" << M << std::endl;

    Eigen::Vector3d x = M.block<1, 3>(2, 0);
    //std::cout << "x=" << x.transpose() << std::endl;

    Eigen::Vector3d e(0.0, 0.0, 1.0);
    Eigen::Matrix3d tQ1;
    double pivot[3];
    pivot[2] = find_householder_matrix_3(x, e, tQ1);
    //std::cout << "tQ1=" << tQ1 << std::endl;

    Eigen::Matrix3d Q1 = tQ1.transpose();
    Eigen::Matrix3d M1 = M*Q1;
    //std::cout << "M1=\n" << M1 << std::endl;

    Eigen::Vector2d y = M1.block<1, 2>(1, 0);
    Eigen::Vector2d e2(0.0, 1.0);
    Eigen::Matrix2d tQ2;
    pivot[1] = find_householder_matrix_2(y, e2, tQ2);

    Eigen::Matrix3d Q2 = Eigen::Matrix3d::Identity();
    Q2.block<2, 2>(0, 0) = tQ2.transpose();

    R = M1*Q2;

    Q = Q2.transpose()*Q1.transpose();

    Eigen::Matrix3d errMat = M - R*Q;
    double err = errMat.norm();
    //std::cout << "error=" << err << std::endl;
    if (err > 1.0e-6)
    {
        std::cerr << "WARNING! Decomposition error = " << errMat.norm() << std::endl;
    }

    R /= R(2, 2);

    return err;
}

Eigen::Vector3f householder(Eigen::Vector3f v, Eigen::Vector3f x)
{
    return x - 2.0f*v.dot(x) / (v.dot(v))* v;
}

Eigen::Vector2f householder(Eigen::Vector2f v, Eigen::Vector2f x)
{
    return x - 2.0f*v.dot(x) / (v.dot(v))* v;
}

float find_householder_matrix_3(Eigen::Vector3f& x, Eigen::Vector3f& e, Eigen::Matrix3f& tQ)
{
    float alpha = x.norm();
    //std::cout << "alpha=" << alpha << std::endl;

    float sign = 1.0f;
    float xe = x.dot(e);
    if (xe < 0.0f)   sign = -1.0f;
    //std::cout << "sign=" << sign << std::endl;

    Eigen::Vector3f te[3];
    te[0] = Eigen::Vector3f(1.0, 0.0, 0.0);
    te[1] = Eigen::Vector3f(0.0, 1.0, 0.0);
    te[2] = Eigen::Vector3f(0.0, 0.0, 1.0);

    Eigen::Vector3f v = x - sign*alpha*e;
    //std::cout << "v=" << v << std::endl;

    tQ.block<3, 1>(0, 0) = householder(v, te[0]);
    tQ.block<3, 1>(0, 1) = householder(v, te[1]);
    tQ.block<3, 1>(0, 2) = householder(v, te[2]);

    return sign*alpha; // NOTE : this return value may not be necessary for later steps.
}

float find_householder_matrix_2(Eigen::Vector2f& x, Eigen::Vector2f& e, Eigen::Matrix2f& tQ)
{
    float alpha = x.norm();

    float sign = 1.0;
    float xe = x.dot(e);
    if (xe < 0.0)   sign = -1.0;

    Eigen::Vector2f te[2];
    te[0] = Eigen::Vector2f(1.0, 0.0);
    te[1] = Eigen::Vector2f(0.0, 1.0);


    Eigen::Vector2f v = x - sign*alpha*e;

    tQ.block<2, 1>(0, 0) = householder(v, te[0]);
    tQ.block<2, 1>(0, 1) = householder(v, te[1]);

    return sign*alpha; // NOTE : this return value may not be necessary for later steps.
}

float RQ_decomp_Householder_3d(float* tM, float* tR, float* tQ)
{
    Eigen::Map<Eigen::Matrix3f> M(tM);
    Eigen::Map<Eigen::Matrix3f> R(tR);
    Eigen::Map<Eigen::Matrix3f> Q(tQ);
    //std::cout << "M=\n" << M << std::endl;

    Eigen::Vector3f x = M.block<1, 3>(2, 0);
    //std::cout << "x=" << x.transpose() << std::endl;

    Eigen::Vector3f e(0.0, 0.0, 1.0);
    Eigen::Matrix3f tQ1;
    float pivot[3];
    pivot[2] = find_householder_matrix_3(x, e, tQ1);
    //std::cout << "tQ1=" << tQ1 << std::endl;

    Eigen::Matrix3f Q1 = tQ1.transpose();
    Eigen::Matrix3f M1 = M*Q1;
    //std::cout << "M1=\n" << M1 << std::endl;

    Eigen::Vector2f y = M1.block<1, 2>(1, 0);
    Eigen::Vector2f e2(0.0, 1.0);
    Eigen::Matrix2f tQ2;
    pivot[1] = find_householder_matrix_2(y, e2, tQ2);

    Eigen::Matrix3f Q2 = Eigen::Matrix3f::Identity();
    Q2.block<2, 2>(0, 0) = tQ2.transpose();

    R = M1*Q2;

    Q = Q2.transpose()*Q1.transpose();

    Eigen::Matrix3f errMat = M - R*Q;
    float err = errMat.norm();
    //std::cout << "error=" << err << std::endl;
    if (err > 1.0e-6)
    {
        std::cerr << "WARNING! Decomposition error = " << errMat.norm() << std::endl;
    }

    R /= R(2, 2);

    return err;
}

double RQ_decomp_GS_3d(double* tM, double* tR, double* tQ)
{
    Eigen::Map<Eigen::Matrix3d> M(tM);
    Eigen::Map<Eigen::Matrix3d> R(tR);
    Eigen::Map<Eigen::Matrix3d> Q(tQ);

    M *= 1.0 / M.row(2).norm();

    R.setConstant(0.0);

    // 3rd row
    R(2, 2) = M.row(2).norm();
    Q.row(2) = M.row(2);
    Q.row(2).normalize();

    // 2nd row
    R(1, 2) = M.row(1).dot(Q.row(2));
    Q.row(1) = M.row(1) - R(1, 2)*Q.row(2);
    R(1, 1) = Q.row(1).norm();
    Q.row(1).normalize();

    // 1st row
    R(0, 1) = M.row(0).dot(Q.row(1));
    R(0, 2) = M.row(0).dot(Q.row(2));
    Q.row(0) = M.row(0) - R(0, 1)*Q.row(1) - R(0, 2)*Q.row(2);
    R(0, 0) = Q.row(0).norm();
    Q.row(0).normalize();

    double det = Q.row(0).cross(Q.row(1)).dot(Q.row(2));

    double k = 1.0;

    if (det < 0.0)
    {
        k = -1.0f;
        Q *= k;

        //std::cout << "k=" << k << std::endl;
    }

    det = Q.row(0).cross(Q.row(1)).dot(Q.row(2));

    if (det < 0.0f)
    {
        std::cerr << "WARNING! Orientation reversing in rotation matrix" << std::endl;
    }

    Eigen::Matrix3d errMat = k*M - R*Q;
    double error = errMat.norm();
    if (error > 1.0e-6)
    {
        std::cerr << "WARNING! Decomposition error = " << errMat.norm() << std::endl;
    }

    return error;

    //Eigen::Vector3d p(pts3d[0].x, pts3d[0].y, pts3d[0].z);
    //std::cout << Q*p << std::endl;
    //std::cout << "det=" << Q.row(2).cross(Q.row(1)).dot(Q.row(0)) << std::endl;

    //std::cout << "R=\n" << R << "\n"
    //    << "Q=\n" << Q << std::endl;

    //std::cout << "RQ=\n" << R*Q << std::endl;
    //std::cout << "M-RQ=" << M - R*Q << std::endl;
    //std::cout << "QTQ=\n" << Q.transpose()*Q << std::endl;

}

float  RQ_decomp_GS_3d(float* tM, float* tR, float* tQ)
{
    Eigen::Map<Eigen::Matrix3f> M(tM);
    Eigen::Map<Eigen::Matrix3f> R(tR);
    Eigen::Map<Eigen::Matrix3f> Q(tQ);

    M *= 1.0f / M.row(2).norm();

    R.setConstant(0.0);

    // 3rd row
    R(2, 2) = M.row(2).norm();
    Q.row(2) = M.row(2);
    Q.row(2).normalize();

    // 2nd row
    R(1, 2) = M.row(1).dot(Q.row(2));
    Q.row(1) = M.row(1) - R(1, 2)*Q.row(2);
    R(1, 1) = Q.row(1).norm();
    Q.row(1).normalize();

    // 1st row
    R(0, 1) = M.row(0).dot(Q.row(1));
    R(0, 2) = M.row(0).dot(Q.row(2));
    Q.row(0) = M.row(0) - R(0, 1)*Q.row(1) - R(0, 2)*Q.row(2);
    R(0, 0) = Q.row(0).norm();
    Q.row(0).normalize();

    float det = Q.row(0).cross(Q.row(1)).dot(Q.row(2));

    float k = 1.0;

    if (det < 0.0)
    {
        k = -1.0f;
        Q *= k;

        //std::cout << "k=" << k << std::endl;
    }

    det = Q.row(0).cross(Q.row(1)).dot(Q.row(2));

    if (det < 0.0f)
    {
        std::cerr << "WARNING! Orientation reversing in rotation matrix" << std::endl;
    }

    Eigen::Matrix3f errMat = k*M - R*Q;
    float error = errMat.norm();
    if (error > 1.0e-6)
    {
        std::cerr << "WARNING! Decomposition error = " << errMat.norm() << std::endl;
    }

    return error;

    //Eigen::Vector3d p(pts3d[0].x, pts3d[0].y, pts3d[0].z);
    //std::cout << Q*p << std::endl;
    //std::cout << "det=" << Q.row(2).cross(Q.row(1)).dot(Q.row(0)) << std::endl;

    //std::cout << "R=\n" << R << "\n"
    //    << "Q=\n" << Q << std::endl;

    //std::cout << "RQ=\n" << R*Q << std::endl;
    //std::cout << "M-RQ=" << M - R*Q << std::endl;
    //std::cout << "QTQ=\n" << Q.transpose()*Q << std::endl;
}

void svd_dp(int m, int n, double* tA, double* tsv, double* tU, double* tV, bool thin)
{
    Eigen::Map<Eigen::MatrixXd> A(tA, m, n);

    //std::cout << A << std::endl;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd;

    if (thin)
        svd = Eigen::JacobiSVD<Eigen::MatrixXd>(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    else
        svd = Eigen::JacobiSVD<Eigen::MatrixXd>(A, Eigen::ComputeFullU | Eigen::ComputeFullV);

    int p = m < n ? m : n;
    Eigen::Map<Eigen::VectorXd> sv(tsv, p);
    sv = svd.singularValues();
    //std::cout << "singular values=" << sv.transpose() << std::endl;

    Eigen::Map<Eigen::MatrixXd> U(tU, m, m), V(tV, n, n);
    U = svd.matrixU();
    V = svd.matrixV();

    //std::cout << "U=\n" << U << std::endl;
    //std::cout << "V=\n" << V << std::endl;


}

//void svd_dp(const Eigen::MatrixXd& A, Eigen::VectorXd& sv, Eigen::MatrixXd& U, Eigen::MatrixXd& V, bool thin)
//{
//    Eigen::JacobiSVD<Eigen::MatrixXd> svd;
//    if(thin)
//        svd = Eigen::JacobiSVD<Eigen::MatrixXd>(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
//    else
//        svd = Eigen::JacobiSVD<Eigen::MatrixXd>(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
//
//    sv = svd.singularValues();
//
//    U = svd.matrixU();
//    V = svd.matrixV();
//}

void svd_dp(void* A, void* tsv, void* tU, void* tV, bool thin)
{

    Eigen::JacobiSVD<Eigen::MatrixXd> svd;
    if(thin)
        svd = Eigen::JacobiSVD<Eigen::MatrixXd>(*(Eigen::MatrixXd*)A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    else
        svd = Eigen::JacobiSVD<Eigen::MatrixXd>(*(Eigen::MatrixXd*)A, Eigen::ComputeFullU | Eigen::ComputeFullV);

    Eigen::VectorXd* sv = (Eigen::VectorXd*)tsv;
    *sv = svd.singularValues();
    //std::cout << svd.singularValues() << std::endl;

    Eigen::MatrixXd* U = (Eigen::MatrixXd*)tU;
    *U = svd.matrixU();

    Eigen::MatrixXd* V = (Eigen::MatrixXd*)tV;
    *V = svd.matrixV();

    //std::cout << "V=" << svd.matrixV() << std::endl;
}

double find_nonzero_kernel_vector(int m, int n, double* tA, double* tv)
{
    Eigen::Map<Eigen::MatrixXd> A(tA, m, n);
    Eigen::Map<Eigen::VectorXd> v(tv, n);
    //std::cout << "mxn=" << m << "x" << n << std::endl;

    //std::cout << "A=\n" << A << std::endl;
    //std::cout << "v=" << v.transpose() << std::endl;

    Eigen::VectorXd sv;
    Eigen::MatrixXd U, V;

    svd_dp(&A, &sv, &U, &V, false);

    v = V.col(n - 1);

    //std::cout << V.col(n - 1) << std::endl;

    if (m < n)
        return 0.0;
    else
        return sv(n - 1);
}

void quaternionToRotation(double* quat, double* tR)
{
    Eigen::Vector4d q(quat);
    Eigen::Map<Eigen::Matrix3d> R(tR);

    double s = q(0);

    double x = q(1);
    double y = q(2);
    double z = q(3);

    R(0, 0) = 1.0 - 2.0*y*y - 2.0*z*z;
    R(1, 0) = 2.0*x*y + 2.0*s*z;
    R(2, 0) = 2.0*x*z - 2.0*s*y;

    R(0, 1) = 2.0*x*y - 2.0*s*z;
    R(1, 1) = 1.0 - 2.0*x*x - 2.0*z*z;
    R(2, 1) = 2.0*y*z + 2.0*s*x;

    R(0, 2) = 2.0*x*z + 2.0*s*y;
    R(1, 2) = 2.0*y*z - 2.0*s*x;
    R(2, 2) = 1.0 - 2.0*x*x - 2.0*y*y;
}

void compute_covariance_matrix(int n, double* tp, double* tcenter, double* tC)
{
    Eigen::Vector3d* p = (Eigen::Vector3d*)tp;

    Eigen::Vector3d c(tcenter);

    Eigen::Map<Eigen::Matrix3d> C(tC);

    Eigen::MatrixXd C1(3, n), C2(n, 3);
    for (int i = 0; i < n; ++i)
    {
        C1.block<3, 1>(0, i) = (p[i]-c);
        C2.block<1, 3>(i, 0) = (p[i]-c).transpose();
    }

    C = C1*C2;
}