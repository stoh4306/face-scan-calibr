#ifndef HOMOGRAPHY_2D_H_
#define HOMOGRAPHY_2D_H_

#include "vec.h"

#include <Eigen/Dense>

class Homography2D
{
public:
    Homography2D() : numPoints_(0), srcPts_(NULL), tgtPts_(NULL) {}

    void setPoints(int n, mg::Vector2d* from, mg::Vector2d* to)
    {
        numPoints_ = n;

        srcPts_ = from;
        tgtPts_ = to;
    }

    void setTransform(Eigen::Matrix3d H)
    {
        H_ = H;
    }

    Eigen::Matrix3d computeHomography(bool normalized = true);

    void find_normalize_transform(int n, mg::Vector2d* p, double expectedMag, Eigen::Matrix3d& T);
    mg::Vector2d translate_to_origin(int N, mg::Vector2d* p);
    double average_magnitude(int N, mg::Vector2d* p);
    void scaleBy(int N, mg::Vector2d* p, double k);

    // Absolute conic
    void find_image_absolute_conic_equations(Eigen::VectorXd& v, Eigen::VectorXd& w);
    void symm_bilinear_form_2_rowvec(Eigen::Vector3d v, Eigen::Vector3d w, Eigen::VectorXd& a);

    // statistics
    double forward_projection_error();
    double forward_projection_error(int n, mg::Vector2d* src, mg::Vector2d* tgt, Eigen::Matrix3d& H);

public:
    int numPoints_;

    mg::Vector2d* srcPts_;
    mg::Vector2d* tgtPts_;

    Eigen::Matrix3d H_;
};

#endif