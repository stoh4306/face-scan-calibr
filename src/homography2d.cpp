#include "homography2d.h"

#include <vector>

void Homography2D::symm_bilinear_form_2_rowvec(Eigen::Vector3d v, Eigen::Vector3d w, Eigen::VectorXd& a)
{
    a(0) = v(0)*w(0);
    a(1) = v(0)*w(1) + v(1)*w(0);
    a(2) = v(0)*w(2) + v(2)*w(0);
    a(3) = v(1)*w(1);
    a(4) = v(1)*w(2) + v(2)*w(1);
    a(5) = v(2)*w(2);
}

void Homography2D::find_image_absolute_conic_equations(Eigen::VectorXd& v, Eigen::VectorXd& w)
{
    Eigen::Vector3d h[3];
    h[0] = H_.block<3, 1>(0, 0);
    h[1] = H_.block<3, 1>(0, 1);
    h[2] = H_.block<3, 1>(0, 2);

    // Equation 1 : h1^T w h2 = 0
    symm_bilinear_form_2_rowvec(h[0], h[1], v);

    // Equation 2 : h1^T w h1 = h2^T w h2
    symm_bilinear_form_2_rowvec(h[0], h[0], w);

    Eigen::VectorXd tw(6);
    symm_bilinear_form_2_rowvec(h[1], h[1], tw);

    w = w - tw;
}

mg::Vector2d Homography2D::translate_to_origin(int n, mg::Vector2d* p)
{
    mg::Vector2d c(0.0f, 0.0f);

    for (int i = 0; i < n; ++i)
    {
        c += p[i];
    }

    if (n > 0)
    {
        c /= (float)n;
    }

    for (int i = 0; i < n; ++i)
    {
        p[i] -= c;
    }

    return c;
}

double Homography2D::average_magnitude(int n, mg::Vector2d* p)
{
    double mag = 0.0;
    for (int i = 0; i < n; ++i)
    {
        mag += p[i].norm();
    }

    if (n > 0)
    {
        mag /= (float)n;
    }

    return mag;
}

void Homography2D::scaleBy(int n, mg::Vector2d* p, double k)
{
    for (int i = 0; i < n; ++i)
    {
        p[i] *= k;
    }
}

void Homography2D::find_normalize_transform(int n, mg::Vector2d* p, double expectedMag, Eigen::Matrix3d& T)
{
    mg::Vector2d center = translate_to_origin(n, p);
    double k = expectedMag / average_magnitude(n, p);
    scaleBy(n, p, k);

    T = Eigen::Matrix3d::Identity();
    T.block<2, 2>(0, 0) *= k;
    T(0, 2) = -k*center.x;
    T(1, 2) = -k*center.y;
}

double Homography2D::forward_projection_error()
{
    int n = numPoints_;

    double error = 0.0;

    for (int i = 0; i < n; ++i)
    {
        Eigen::Vector3d p(srcPts_[i].x, srcPts_[i].y, 1.0f);
        
        Eigen::Vector3d ttq = H_*p;

        Eigen::Vector2d tq(ttq(0) / ttq(2), ttq(1) / ttq(2));

        Eigen::Vector2d q(tgtPts_[i].data());

        error += (q - tq).norm();//(q - tq).squaredNorm();
    }

    if (n > 0)
    {
        error /= (double)n;
    }

    return error;
}

double Homography2D::forward_projection_error(int n, mg::Vector2d* src, mg::Vector2d* tgt, Eigen::Matrix3d& H)
{
    double error = 0.0;

    for (int i = 0; i < n; ++i)
    {
        Eigen::Vector3d p(src[i].x, src[i].y, 1.0f);

        Eigen::Vector3d ttq = H*p;

        Eigen::Vector2d tq(ttq(0) / ttq(2), ttq(1) / ttq(2));

        Eigen::Vector2d q(tgt[i].data());

        error += (q - tq).squaredNorm();
    }

    if (n > 0)
    {
        error /= (double)n;
    }

    return error;
}

Eigen::Matrix3d Homography2D::computeHomography(bool normalized)
{
    std::vector< mg::Vector2d > tempSrcPts(numPoints_);
    std::vector< mg::Vector2d > tempTgtPts(numPoints_);

    memcpy(tempSrcPts.data(), srcPts_, sizeof(mg::Vector2d)*numPoints_);
    memcpy(tempTgtPts.data(), tgtPts_, sizeof(mg::Vector2d)*numPoints_);

    Eigen::Matrix3d srcTrans, tgtTrans; // Normalization transform
    srcTrans.setIdentity();
    tgtTrans.setIdentity();

    std::vector< mg::Vector2d > & p = tempSrcPts;
    std::vector< mg::Vector2d > & q = tempTgtPts;

    if (normalized)
    {
        const double SQRT_TWO = pow(2.0, 0.5);

        find_normalize_transform(numPoints_, p.data(), SQRT_TWO, srcTrans);
        find_normalize_transform(numPoints_, q.data(), SQRT_TWO, tgtTrans);

        //for (int i = 0; i < numPoints_; ++i)
        //{
        //    std::cout << p[i] << "->" << q[i] << "\n";
        //}

        //std::cout << "src transform : \n" << srcTrans << std::endl;
        //std::cout << "tgt transform : \n" << tgtTrans << std::endl;
    }

    // Construct a matrix A to find the homography.
    // Note that the homography is determined by solving Ah = 0.
    const int n = numPoints_;
    Eigen::MatrixXd A(2 * n, 9);

    A.setZero();

    for (int i = 0; i < n; ++i)
    {
        Eigen::Vector3d v = Eigen::Vector3d(p[i].x, p[i].y, 1.0f);
        A.block<1, 3>(2 * i + 0, 3) = -v;
        A.block<1, 3>(2 * i + 0, 6) = q[i].y * v;
        A.block<1, 3>(2 * i + 1, 0) = v;
        A.block<1, 3>(2 * i + 1, 6) = -q[i].x * v;
    }

    Eigen::JacobiSVD<Eigen::MatrixXd> svd;
    if(A.rows() < A.cols())
        svd = Eigen::JacobiSVD<Eigen::MatrixXd>(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
    else
        svd = Eigen::JacobiSVD<Eigen::MatrixXd>(A, Eigen::ComputeThinU | Eigen::ComputeThinV);

    //Eigen::VectorXd sv = svd.singularValues();
    //std::cout << "singular values : " << sv.transpose() << std::endl;
    //int numEigenvalues = sv.size();
    //std::cout << " Smallest singular value : " << sv(numEigenvalues - 1) << std::endl;
    Eigen::VectorXd V = svd.matrixV().col(A.cols() - 1);//svd.matrixV().col(numEigenvalues - 1);

    //std::cout << "mag(V)=" << V.norm() << std::endl;

    H_(0, 0) = V(0);
    H_(0, 1) = V(1);
    H_(0, 2) = V(2);

    H_(1, 0) = V(3);
    H_(1, 1) = V(4);
    H_(1, 2) = V(5);

    H_(2, 0) = V(6);
    H_(2, 1) = V(7);
    H_(2, 2) = V(8);
    
    //float error = forward_projection_error(n, p.data(), q.data(), H_);

    H_ = tgtTrans.inverse()*(H_*srcTrans);

    //float error = forward_projection_error();
    //std::cout << "error(forward)=" << error << std::endl;

    return H_;
}