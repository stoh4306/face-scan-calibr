#include "calibr.h"
#include "reconst.h"

#include "numer_math.h"
#include "homography2d.h"
#include "vec.h"

#include <vector>
#include <fstream>

#include <opencv2/opencv.hpp>

bool read_calibrations(const char* K_file, const char* CP_file, float* tK, float* tCP)
{
    std::ifstream inFile(K_file);
    if (!inFile.is_open())  return false;

    Eigen::Map<Eigen::Matrix3f> K(tK);
    inFile >> K(0, 0) >> K(0, 1) >> K(0, 2)
        >> K(1, 0) >> K(1, 1) >> K(1, 2)
        >> K(2, 0) >> K(2, 1) >> K(2, 2);
    inFile.close();

    std::ifstream inFile2(CP_file);
    if (!inFile2.is_open())  return false;

    Eigen::Map<Eigen::Matrix4f> CP(tCP);
    inFile2
        >> CP(0, 0) >> CP(0, 1) >> CP(0, 2) >> CP(0, 3)
        >> CP(1, 0) >> CP(1, 1) >> CP(1, 2) >> CP(1, 3)
        >> CP(2, 0) >> CP(2, 1) >> CP(2, 2) >> CP(2, 3)
        >> CP(3, 0) >> CP(3, 1) >> CP(3, 2) >> CP(3, 3);
    inFile2.close();

    return true;
}

bool find_check_corners(int w, int h, unsigned char* color_check_img, int board_width, int board_height, double* d_corner_pts, bool showImage)
{
    cv::Mat img(h, w, CV_8UC3, color_check_img);
    //cv::imshow("input image", img);
    //cv::waitKey(0);

    std::vector< cv::Point2f > corner_pts;

    cv::Size boardSize(board_width, board_height);
    bool found = cv::findChessboardCorners(img, boardSize, corner_pts);

    //for (int i = 0; i < corner_pts.size(); ++i)
    //{
    //    std::cout << corner_pts[i] << std::endl;
    //}

    //cv::Mat tempImg;
    //img.copyTo(tempImg);
    //cv::drawChessboardCorners(tempImg, boardSize, cv::Mat(corner_pts), true);

    //cv::imshow("check_corners", tempImg);
    //cv::waitKey(0);

    if (found)
    {
        cv::Mat imgGray;
        cv::cvtColor(img, imgGray, cv::COLOR_BGR2GRAY);
        cv::cornerSubPix(imgGray, corner_pts, cv::Size(11, 11),
            cv::Size(-1, -1), cv::TermCriteria( cv::TermCriteria::EPS + cv::TermCriteria::MAX_ITER, 30, 0.1));

        for (int i = 0; i < (int)corner_pts.size(); ++i)
        {
            d_corner_pts[2 * i + 0] = corner_pts[i].x;
            d_corner_pts[2 * i + 1] = corner_pts[i].y;
        }

        //std::cout << "SUB-PIXEL ACCURACY" << std::endl;
        //for (int i = 0; i < corner_pts.size(); ++i)
        //{
        //    std::cout << corner_pts[i] << std::endl;
        //}

        if (showImage)
        {
            cv::drawChessboardCorners(img, boardSize, cv::Mat(corner_pts), true);
            cv::imshow("check_corners", img);
            cv::waitKey(0);
        }

        return true;
    }
    else
    {
        std::cout << "Error, can't find corner points in the input check image" << std::endl;
        return false;
    }
}

void find_homography_plane_to_image(int n, double* plane_pts, double* image_pts, double* H_colmaj)
{
    Homography2D hom2d;

    hom2d.setPoints(n, (mg::Vector2d*)plane_pts, (mg::Vector2d*)image_pts);
    
    Eigen::Map<Eigen::Matrix3d> H(H_colmaj);
    H = hom2d.computeHomography(true);

    double error = hom2d.forward_projection_error();

    std::cout << "- Homography computed for plane and image : error=" << (float)error << "\n"
        << "H = \n" << H << std::endl;

}

void compute_iac_from_plane2image_homography(int numHom, double* tH, double* tomega)
{
    //-----------------------------------------------------------------------------------------
    // Circular points in the plane of infinity : (1, i, 0), (1,-i, 0)
    // The image of circular points for a homography is lying on the image of absolute conic w.
    // H(1,i,0)=h1+i*h2, H(1,-i,0)=h1-i*h2.
    // For each homography H= [h1, h2, h3], two equations are given:
    //   h1^T w h2 = 0  and  h1^T w h1 = h2^T w h2
    //-----------------------------------------------------------------------------------------
    std::cout << "- Compute the image of the absolute conic : " << std::endl;
    Eigen::MatrixXd A(2 * numHom, 6);

    for (int i = 0; i < numHom; ++i)
    {
        Eigen::Map<Eigen::Matrix3d> H(tH + 9 * i);

        Homography2D hom;
        hom.setTransform(H);

        Eigen::VectorXd v(6), w(6);
        hom.find_image_absolute_conic_equations(v, w);

        A.block<1, 6>(2 * i + 0, 0) = v;
        A.block<1, 6>(2 * i + 1, 0) = w;
    }

    // Compute a unit kernel vector x minimizing ||Ax||
    Eigen::VectorXd V(6);
    find_nonzero_kernel_vector(2 * numHom, 6, A.data(), V.data());

    //Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
    //Eigen::VectorXd sv = svd.singularValues();
    //int numEigenvalues = sv.size();
    ////std::cout << " Smallest singular value : " << sv(numEigenvalues - 1) << std::endl;
    //Eigen::VectorXd V = svd.matrixV().col(numEigenvalues - 1);

    //std::cout << "mag(V)=" << V.norm() << std::endl;
    Eigen::Map<Eigen::Matrix3d> omega(tomega);
    omega(0, 0) = V(0);
    omega(0, 1) = V(1);
    omega(0, 2) = V(2);

    omega(1, 1) = V(3);
    omega(1, 2) = V(4);

    omega(2, 2) = V(5);

    omega(1, 0) = omega(0, 1);
    omega(2, 0) = omega(0, 2);
    omega(2, 1) = omega(1, 2);

    std::cout << ". omega(IAC) : \n" << omega << "\n" << std::endl;

    //Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(omega);
    //std::cout << "eigen values : " << eigensolver.eigenvalues().transpose() << std::endl;
}

void find_calibration_matrix_from_iac(double* tw, double* tK, int method)
{
    std::cout << "- Compute K from the image of absolute conic" << std::endl;

    Eigen::Map<Eigen::Matrix3d> w(tw), K(tK);

    Eigen::Matrix3d inv_w = w.inverse();
    //std::cout << "\nomega^(-1)=\n" << inv_w << std::endl;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(inv_w, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::Vector3d sv = svd.singularValues();
    //std::cout << "singular values = " << sv.transpose() << std::endl;

    Eigen::Matrix3d V = svd.matrixV();
    //std::cout << "det=" << V.determinant() << "\n"
    //    << "VV^T=\n" << V*V.transpose() << std::endl;

    Eigen::Matrix3d D = Eigen::Matrix3d::Identity();
    D(0, 0) = sqrt(sv(0));
    D(1, 1) = sqrt(sv(1));
    D(2, 2) = sqrt(sv(2));

    Eigen::Matrix3d M = V*D;
    //std::cout << "M=\n" << M << std::endl;
    //std::cout << "M*M^T=\n" << M*M.transpose() << std::endl;
    //std::cout << "(M*M^T)^(-1)=\n" << (M*M.transpose()).inverse() << std::endl;
    double cholesky_error = (inv_w - (M*M.transpose())).norm() / 9.0;
    std::cout << ". Cholesky decomp. error=" << (float)cholesky_error << std::endl;

    Eigen::Matrix3d R;

    double err = 1.0;
    if (method == 0) // Gramschmidt method
    {
        double err = RQ_decomp_GS_3d(M.data(), K.data(), R.data());

        if (err < 1.0e-6)
        {
            double k = M(2, 2) / R(2, 2);

            std::cout << ". Camera matrix(K) : \n" << K << std::endl;
            //std::cout << "\nR=\n" << R << std::endl;
            //std::cout << "k*K*R=" << k*K*R << std::endl;
        }
        else // Householder approach taken for error cases
        {
            // The householder method is known to be more stable than the Gramschmidt process.
            std::cout << ". Householder method for RQ-decomposition:" << std::endl;
            RQ_decomp_Householder_3d(M.data(), K.data(), R.data());
            std::cout << ". Camera matrix(K) : \n" << K << std::endl;
        }
    }

    else if (method == 1) // Householder method
    {
        // The householder method is known to be more stable than the Gramschmidt process.
        std::cout << ". Householder method for RQ-decomposition:" << std::endl;
        RQ_decomp_Householder_3d(M.data(), K.data(), R.data());
        std::cout << ". Camera matrix(K) : \n" << K << std::endl;
    }
}

void transform_image_pts_to_normalized_coord_frame(int n, double* img_pts, double* K, double* n_img_pts)
{
    Eigen::Map<Eigen::Matrix3d> tK(K);
    Eigen::Matrix3d invK = tK.inverse();
    
    for (int i = 0; i < n; ++i)
    {
        Eigen::Vector3d tx(img_pts[2*i+0], img_pts[2*i+1], 1.0);

        Eigen::Vector3d nx = invK*tx;

        n_img_pts[2 * i + 0] = nx(0);
        n_img_pts[2 * i + 1] = nx(1);
    }
}

void find_fundamental_matrix_normalized_linear(int n, double* left_img_pts,
    double* right_img_pts, double* tF)
{
    std::cout << "- FINDING FUNDAMENTAL MATRIX : " << "\n";

    // Normalize points
    std::vector< mg::Vector2d > tempSrcPts(n);
    std::vector< mg::Vector2d > tempTgtPts(n);
    for (int i = 0; i < n; ++i)
    {
        tempSrcPts[i] = mg::Vector2d(left_img_pts[2 * i + 0], left_img_pts[2 * i + 1]);
        tempTgtPts[i] = mg::Vector2d(right_img_pts[2 * i + 0], right_img_pts[2 * i + 1]);
    }
    std::cout << tempSrcPts[0] << "-" << tempTgtPts[0] << std::endl;

    Eigen::Matrix3d srcTrans, tgtTrans; // Normalization transform
    srcTrans.setIdentity();
    tgtTrans.setIdentity();

    std::vector< mg::Vector2d > p = tempSrcPts;
    std::vector< mg::Vector2d > q = tempTgtPts;

    const double SQRT_TWO = pow(2.0, 0.5);

    normalize_pts_2d(n, (double*)p.data(), SQRT_TWO, srcTrans.data());
    normalize_pts_2d(n, (double*)q.data(), SQRT_TWO, tgtTrans.data());

    //for (int i = 0; i < numPoints_; ++i)
    //{
    //    std::cout << p[i] << "->" << q[i] << "\n";
    //}
    std::cout << "src transform : \n" << srcTrans << std::endl;
    std::cout << "tgt transform : \n" << tgtTrans << std::endl;
   
    Eigen::MatrixXd A(n, 9);

    for (int i = 0; i < n; ++i)
    {
        Eigen::VectorXd tv(9);

        mg::Vector3d tp, tq;  // q^T F p = 0
        tp = mg::Vector3d(p[i].x, p[i].y, 1.0);
        tq = mg::Vector3d(q[i].x, q[i].y, 1.0);

        tv(0) = tq.x*tp.x, tv(1) = tq.x*tp.y, tv(2) = tq.x*tp.z;
        tv(3) = tq.y*tp.x, tv(4) = tq.y*tp.y, tv(5) = tq.y*tp.z;
        tv(6) = tq.z*tp.x, tv(7) = tq.z*tp.y, tv(8) = tq.z*tp.z;

        A.block<1, 9>(i, 0) = tv;
    }
    //std::cout << rowIndex << " " << totalNumPts << std::endl;

    //double ssv;
    Eigen::VectorXd v(9);
    find_nonzero_kernel_vector(n, 9, A.data(), v.data());

    Eigen::Map<Eigen::Matrix3d> F(tF);
    F(0, 0) = v(0);
    F(0, 1) = v(1);
    F(0, 2) = v(2);

    F(1, 0) = v(3);
    F(1, 1) = v(4);
    F(1, 2) = v(5);

    F(2, 0) = v(6);
    F(2, 1) = v(7);
    F(2, 2) = v(8);
    std::cout << "F=" << F << std::endl;

    //
    //double error = bilinear_form(F, totalNumPts, q.data(), p.data());
    //std::cout << "error=" << error << std::endl;
    //std::cout << "det(F)=" << F.determinant() << std::endl;

    // Singular value decomposition
    Eigen::Vector3d sv;
    Eigen::Matrix3d U, V;
    svd_dp(3,3, F.data(), sv.data(), U.data(), V.data());
    std::cout << "singular values=" << sv.transpose() << std::endl;

    // Make F as rank(F)=2
    Eigen::Matrix3d D = Eigen::Matrix3d::Identity();
    D(0, 0) = sv(0);
    D(1, 1) = sv(1);
    D(2, 2) = 0.0;

    F = U*(D*V.transpose());

    svd_dp(3, 3, F.data(), sv.data(), U.data(), V.data());

    //error = bilinear_form(F, totalNumPts, q.data(), p.data());
    //std::cout << "-> error=" << error << "(rank(F)=2)" << std::endl;

    F = tgtTrans.transpose()*(F*srcTrans);

    std::cout << ". F = \n" << F << std::endl;

    //std::cout << p[0] << "->" << tempSrcPts[0] << std::endl;
    double error = bilinear_form_RP2(tF, n, (double*)tempTgtPts.data(), (double*)tempSrcPts.data());
    std::cout << "error = " << error << std::endl;// " (after denormalization)" << std::endl;
    std::cout << "det(F)=" << F.determinant() << std::endl;

}

void print_camera_pose_center(const char* info, Eigen::Matrix3d R, Eigen::Vector3d t)
{
    Eigen::Matrix3d tR = R.transpose();

    std::cout
        << "\n"
        << "- " << info << " : \n"
        << ". cam pose :\n" << tR << "\n"
        << ". cam center : " << (-tR*t).transpose() << std::endl;
}

double rotational_differ(Eigen::Matrix3d R)
{
    return abs(1.0 - R(0, 0)) + abs(1.0 - R(1, 1)) + abs(1.0 - R(2, 2));
}

void compute_camera_matrix_from_essential_matrix(int n, double* np, double* nq,
    double* K1, double* K2, double* E, double* tP1, double* tP2)
{
    //---------------------------
    // NOTE : x2^T E x1 = 0
    //---------------------------

    // Scale the essential matrix
    Eigen::Vector3d lambda;
    Eigen::Matrix3d U, V;
    svd_dp(3,3, E, lambda.data(), U.data(), V.data());//singular_value_decomp(E, lambda, U, V);

    std::cout << ". singular values of E = " << lambda.transpose() << std::endl;
    std::cout << "U=\n" << U << ",\n" << "V=\n" << V << std::endl;

    Eigen::Matrix3d Z, W;
    Z.setZero(), Z(0, 1) = 1.0, Z(1, 0) = -1.0;
    W.setZero(), W(0, 1) = -1.0, W(1, 0) = 1.0, W(2, 2) = 1.0;

    Eigen::Vector3d u3 = U.col(2);
    std::cout << "u3=" << u3.transpose() << std::endl;

    // Investigate four cases
    Eigen::Matrix3d R;
    Eigen::Vector3d t;

    // Case 1
    R = U*(W*V.transpose());
    t = u3;

    if (R.determinant() < 0.0)
    {
        R *= -1.0;
        t *= -1.0;
    }

    Eigen::Matrix3d tR = R.transpose();

    print_camera_pose_center("case 1", R, t);

    double rot_diff = rotational_differ(tR);
    Eigen::Vector3d cam_center = -tR*t;
    std::cout << "rot_diff=" << rot_diff << ", center_x=" << cam_center(0) << std::endl;

    Eigen::Map<Eigen::MatrixXd> P2(tP2, 3, 4);
    if (rot_diff < 0.5 && cam_center(0) > 0.0) // Right camera is on the right side of the left camera
    {
        P2.block<3, 3>(0, 0) = R;
        P2.block<3, 1>(0, 3) = t;

        std::cout << "***** CASE 1 selected" << std::endl;
        return;
    }

    // Case 2
    t *= -1.0;
    print_camera_pose_center("case 2", R, t);

    rot_diff = rotational_differ(tR);
    cam_center = -tR*t;
    std::cout << "rot_diff=" << rot_diff << ", center_x=" << cam_center(0) << std::endl;
    if (rot_diff < 0.5 && cam_center(0) > 0.0) // Right camera is on the right side of the left camera
    {
        P2.block<3, 3>(0, 0) = R;
        P2.block<3, 1>(0, 3) = t;

        std::cout << "***** CASE 2 selected" << std::endl;

        return;
    }

    // Case 3
    R = U*(W.transpose()*V.transpose());
    t = u3;

    if (R.determinant() < 0.0)
    {
        R *= -1.0;
        t *= -1.0;
    }

    tR = R.transpose();

    print_camera_pose_center("case 3", R, t);

    rot_diff = rotational_differ(tR);
    cam_center = -tR*t;
    std::cout << "rot_diff=" << rot_diff << ", center_x=" << cam_center(0) << std::endl;
    if (rot_diff < 0.5 && cam_center(0) > 0.0) // Right camera is on the right side of the left camera
    {
        P2.block<3, 3>(0, 0) = R;
        P2.block<3, 1>(0, 3) = t;

        std::cout << "***** CASE 3 selected" << std::endl;

        return;
    }

    // Case 4
    t *= -1.0;
    print_camera_pose_center("case 4", R, t);

    rot_diff = rotational_differ(tR);
    cam_center = -tR*t;
    std::cout << "rot_diff=" << rot_diff << ", center_x=" << cam_center(0) << std::endl;
    if (rot_diff < 0.5 && cam_center(0) > 0.0) // Right camera is on the right side of the left camera
    {
        P2.block<3, 3>(0, 0) = R;
        P2.block<3, 1>(0, 3) = t;

        std::cout << "***** CASE 4 selected" << std::endl;

        return;
    }

    std::cout << "***** ERROR****** Camera matrix cannot be selected" << std::endl;

    return;
}

int SingleCameraCalibr::computeTotalNumPoints()
{
    total_num_pts_ = 0;
    for (int i = 0; i < num_images_; ++i)
        total_num_pts_ += num_pts_[i];
    return total_num_pts_;
}

void SingleCameraCalibr::convert_image_pts_to_normalized_coord()
{
    computeTotalNumPoints();

    if (normalized_image_pts_)  delete[] normalized_image_pts_;
    normalized_image_pts_ = new double[2 * total_num_pts_];

    transform_image_pts_to_normalized_coord_frame(total_num_pts_, image_pts_, K_, normalized_image_pts_);
}

void SingleCameraCalibr::find_intrinsic_matrix()
{
    int method = 0; // 0: Gram-Schmidt, 1: House-holder 
    find_calibration_matrix_from_iac(Omega_, K_, method);
}

void SingleCameraCalibr::compute_image_absolute_conics()
{
    compute_iac_from_plane2image_homography(num_images_, H_.data(), Omega_);
}

void SingleCameraCalibr::compute_homography_planes_to_images()
{
    H_.resize(num_images_ * 9);

    int start_index_pt = 0;
    for (int i = 0; i < num_images_; ++i)
    {
        Eigen::Map<Eigen::Matrix4d> H(H_.data() + 9 * i);
        find_homography_plane_to_image(num_pts_[i], plane_pts_ + 2 * start_index_pt, image_pts_ + 2 * start_index_pt, H.data());

        start_index_pt += num_pts_[i];
    }
}

SingleCameraCalibr::~SingleCameraCalibr()
{
    if (normalized_image_pts_)
    {
        delete[] normalized_image_pts_;
        normalized_image_pts_ = NULL;
    }
}

CameraCalibrWithCheck::~CameraCalibrWithCheck()
{
    if (num_pts_)
    {
        delete[] num_pts_;
        num_pts_ = NULL;
    }

    if (plane_pts_)
    {
        delete[] plane_pts_;
        plane_pts_ = NULL;
    }

    //if (image_pts_)
    //{
    //    delete[] image_pts_;
    //    image_pts_ = NULL;
    //}
}

int CameraCalibrWithCheck::find_corner_points(const char* filename, bool vis_corner_pts)
{
    corner_pts_.reserve(2*num_images_*board_width_*board_height_);

    if (num_pts_)   delete[] num_pts_;
    num_pts_ = new int[num_images_];

    int total_num_pts = 0;

    //std::cout << "w=" << w_ << ", h=" << h_ << std::endl; 

    //std::cout << "Finding corner points..." << std::endl;

    for (int i = 0; i < num_images_; ++i)
    {
        //std::cout << "i=" << i << std::endl;
        std::vector<double> corners(2*board_width_*board_height_);
        bool found = find_check_corners(w_, h_, 
            check_img_data_ + i * 3 * w_*h_, board_width_, board_height_, corners.data(), vis_corner_pts);

        if (found)
        {
            corner_pts_.insert(corner_pts_.end(), corners.begin(), corners.end());
            num_pts_[i] = (int)corners.size() / 2;
        }
        else
        {
            num_pts_[i] = 0;
            std::cout << "No corner points found for image[" << i << "]" << std::endl;
			return i;
        }

        total_num_pts += num_pts_[i];
    }

    total_num_pts_ = total_num_pts;

    image_pts_ = corner_pts_.data();
    std::cout << "image point pointer : " << image_pts_ << std::endl;

    std::ofstream outFile(filename);
    outFile
        << "NUM_PTS 	= " << total_num_pts << "\n"
        << "OFFSET      = 0" << "\n"
        << "PTS_DATA    = " << "\n";

    for (int i = 0; i < total_num_pts_; ++i)
    {
        outFile << corner_pts_[2 * i + 0] << " " << corner_pts_[2 * i + 1] << "\n";
    }
    outFile.close();

    if (plane_pts_)     delete[] plane_pts_;
    plane_pts_ = new double[2 * total_num_pts];

    const double a = board_square_size_;
    const double xmax = (board_width_ - 1)*a;
    
    std::vector<double> tp(2 * board_width_*board_height_);
    int index = 0;
    for (int j = 0; j < board_height_; ++j)
    {
        double y = j*a;
        for (int i = 0; i < board_width_; ++i)
        {
            //double x = i*a;
            double x = xmax - i*a;
            
            tp[2 * index + 0] = x;
            tp[2 * index + 1] = y;

            index++;
        }
    }
    std::cout << mg::Vector2d(image_pts_) << ", " << mg::Vector2d(image_pts_ + 2 * (total_num_pts-1)) << std::endl;
    std::cout << mg::Vector2d(corner_pts_[0], corner_pts_[1]) << ", " << mg::Vector2d(corner_pts_[2 * (total_num_pts - 1) + 0], corner_pts_[2 * (total_num_pts - 1) + 1]) << std::endl;
    std::cout << mg::Vector2d(tp.data()) << ", " << mg::Vector2d(tp.data() + 2 * (board_width_*board_height_ - 1)) << std::endl;

    double* p = plane_pts_;

    for (int i = 0; i < num_images_; ++i)
    {
        if (num_pts_[i] > 0)
        {
            memcpy(p, tp.data(), sizeof(double) * 2 * board_width_*board_height_);
            p += 2 * num_pts_[i];
        }
    }
	return -1;
}
int CameraCalibrWithCheck::find_corner_points(std::vector<Eigen::Vector2d>* corner_pts, bool vis_corner_pts)
{

	corner_pts_.reserve(2 * num_images_*board_width_*board_height_);

	if (num_pts_)   delete[] num_pts_;
	num_pts_ = new int[num_images_];

	int total_num_pts = 0;
	for (int i = 0; i < num_images_; ++i)
	{
		std::vector<double> corners(2 * board_width_*board_height_);
		bool found = find_check_corners(w_, h_,
			check_img_data_ + i * 3 * w_*h_, board_width_, board_height_, corners.data(), vis_corner_pts);

		if (found)
		{
			corner_pts_.insert(corner_pts_.end(), corners.begin(), corners.end());
			num_pts_[i] = (int)corners.size() / 2;
		}
		else
		{
			num_pts_[i] = 0;
			std::cout << "No corner points found for image[" << i << "]" << std::endl;
			return i;
		}

		total_num_pts += num_pts_[i];
	}


	total_num_pts_ = total_num_pts;

	image_pts_ = corner_pts_.data();
	std::cout << "image point pointer : " << image_pts_ << std::endl;

	corner_pts->clear();
	for (int i = 0; i < total_num_pts_; ++i)
	{
		corner_pts->insert(corner_pts->end(), Eigen::Vector2d(corner_pts_[2 * i + 0], corner_pts_[2 * i + 1]));
	}

	if (plane_pts_)     delete[] plane_pts_;
	plane_pts_ = new double[2 * total_num_pts];

	const double a = board_square_size_;
	const double xmax = (board_width_ - 1)*a;
	

	std::vector<double> tp(2 * board_width_*board_height_);
	int index = 0;
	for (int j = 0; j < board_height_; ++j)
	{
		double y = j * a;
		for (int i = 0; i < board_width_; ++i)
		{
			double x = xmax - i * a;

			tp[2 * index + 0] = x;
			tp[2 * index + 1] = y;

			index++;
		}
	}


	std::cout << mg::Vector2d(image_pts_) << ", " << mg::Vector2d(image_pts_ + 2 * (total_num_pts - 1)) << std::endl;
	std::cout << mg::Vector2d(corner_pts_[0], corner_pts_[1]) << ", " << mg::Vector2d(corner_pts_[2 * (total_num_pts - 1) + 0], corner_pts_[2 * (total_num_pts - 1) + 1]) << std::endl;
	std::cout << mg::Vector2d(tp.data()) << ", " << mg::Vector2d(tp.data() + 2 * (board_width_*board_height_ - 1)) << std::endl;

	double* p = plane_pts_;

	for (int i = 0; i < num_images_; ++i)
	{
		if (num_pts_[i] > 0)
		{
			memcpy(p, tp.data(), sizeof(double) * 2 * board_width_*board_height_);
			p += 2 * num_pts_[i];
		}
	}
	return -1;

}

void StereoCamCalibr::reconstruct_3d_points(std::string cornerPointfileName)
{
    int N = leftCamPtr_->totalNumPoints();
    double* np = leftCamPtr_->normalized_image_pts_;
    double* nq = rightCamPtr_->normalized_image_pts_;

    std::vector<Eigen::Vector3d> X(N);
    for (int i = 0; i < N; ++i)
    {
        linear_triangulation(np+2*i, nq+2*i, P_l_, P_r_, X[i].data());
    }

    write_point_cloud_in_ply((outFolder+"\\reconst.ply").c_str(), (int)X.size(), X.data()->data());

    std::ofstream outFile(outFolder + cornerPointfileName);
    outFile << N << "\n";
    for (int i = 0; i < N; ++i)
    {
        outFile << X[i](0) << " " << X[i](1) << " " << X[i](2) << "\n";
    }
    outFile.close();
}

void StereoCamCalibr::reconstruct_3d_points(std::vector<Eigen::Vector3d>* corner_pts)
{
	int N = leftCamPtr_->totalNumPoints();
	double* np = leftCamPtr_->normalized_image_pts_;
	double* nq = rightCamPtr_->normalized_image_pts_;

	std::vector<Eigen::Vector3d> X(N);
	for (int i = 0; i < N; ++i)
	{
		linear_triangulation(np + 2 * i, nq + 2 * i, P_l_, P_r_, X[i].data());
		corner_pts->insert(corner_pts->end(), X[i]);
	}
}

void StereoCamCalibr::reconstruct_3d_points_geom_corr()
{
    int N = leftCamPtr_->totalNumPoints();
    double* np = leftCamPtr_->normalized_image_pts_;
    double* nq = rightCamPtr_->normalized_image_pts_;

    std::vector<Eigen::Vector3d> X(N);
    for (int i = 0; i < N; ++i)
    {
        linear_triangulation(np + 2 * i, nq + 2 * i, P_l_, P_r_, X[i].data());
    }

    write_point_cloud_in_ply("reconst.ply", (int)X.size(), X.data()->data());
}

void StereoCamCalibr::find_fundamental_matrix()
{
    int totalNumPts = leftCamPtr_->totalNumPoints();
    std::cout << "total num pts=" << totalNumPts << std::endl;

    std::cout << "pointer : " << leftCamPtr_->image_pts_ << ", " << rightCamPtr_->image_pts_ << std::endl;
    find_fundamental_matrix_normalized_linear(totalNumPts, leftCamPtr_->image_pts_, rightCamPtr_->image_pts_, F_);
}

void StereoCamCalibr::find_essential_matrix()
{
    // Get the essential matrix from intrinsic parameters
    //Eigen::Matrix3d tK[2]; // 0: left, 1: right
    //read_intrinsic_param(left_K_file, tK[0], "LEFT_K");
    //read_intrinsic_param(right_K_file, tK[1], "RIGHT_K");
    Eigen::Map<Eigen::Matrix3d> K_l(leftCamPtr_->K_);
    Eigen::Map<Eigen::Matrix3d> K_r(rightCamPtr_->K_);

    // Find the essential matrix
    Eigen::Map<Eigen::Matrix3d> F(F_);

    Eigen::Map<Eigen::Matrix3d> E(E_);
    E = K_r.transpose()*(F*K_l);
    std::cout << "\n-Essential matrix(E) : \n" << E << std::endl;
}
//"left_k.txt","right_k.txt""left_cam_pose.txt","right_cam_pose.txt"
void StereoCamCalibr::find_projection_marix(std::string leftKfileName,std::string rightKfileName,std::string leftCamPosefileName,std::string rightCamPosefileName)
{
    // Image points in terms of normalized coordinates
    int numPts = leftCamPtr_->totalNumPoints();
    leftCamPtr_->convert_image_pts_to_normalized_coord();
    rightCamPtr_->convert_image_pts_to_normalized_coord();

    double* np = leftCamPtr_->normalized_image_pts_;
    double* nq = rightCamPtr_->normalized_image_pts_;

    // Camera matrix computed
    Eigen::Map<Eigen::Matrix<double, 3, 4>> P_l(P_l_);
    Eigen::Map<Eigen::Matrix<double, 3, 4>> P_r(P_r_);

    //NOTE : It is always possible to set P_l = [ I | 0 ]
    P_l.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity();
    P_l.block<3, 1>(0, 3) = Eigen::Vector3d::Zero();

    double* K_l = leftCamPtr_->K_;
    double* K_r = rightCamPtr_->K_;
    compute_camera_matrix_from_essential_matrix(numPts, np, nq,
        K_l, K_r, E_, P_l.data(), P_r.data());

    // Compute camera pose and center
    Eigen::Matrix3d R_l, R_r;
    R_l = P_l.block<3, 3>(0, 0);
    R_r = P_r.block<3, 3>(0, 0);

    Eigen::Map<Eigen::Matrix3d> pose_l(campose_l_), pose_r(campose_r_);
    pose_l = R_l.transpose();
    pose_r = R_r.transpose();

    Eigen::Map<Eigen::Vector3d> c_l(camcenter_l_), c_r(camcenter_r_);
    c_l = -pose_l*P_l.block<3, 1>(0, 3);
    c_r = -pose_r*P_r.block<3, 1>(0, 3);

    //
    std::cout << "- CAMERA POSE :\n"
        << "  left :\n" << pose_l << "\n"
        << "  right : \n" << pose_r << "\n"
        << "- CAMERA CENTER :\n"
        << "  left : " << c_l.transpose() << "\n"
        << "  right : " << c_r.transpose() << "\n" << std::endl;

    // Export camera parameter matrices to disk
    Eigen::Matrix3d K1(K_l), K2(K_r);
	std::ofstream outFile_1(outFolder + leftKfileName), outFile_2(outFolder + rightKfileName);
    outFile_1 << K1(0, 0) << " " << K1(0, 1) << " " << K1(0, 2) << "\n"
        << K1(1, 0) << " " << K1(1, 1) << " " << K1(1, 2) << "\n"
        << K1(2, 0) << " " << K1(2, 1) << " " << K1(2, 2) << std::endl;

    outFile_2
        << K2(0, 0) << " " << K2(0, 1) << " " << K2(0, 2) << "\n"
        << K2(1, 0) << " " << K2(1, 1) << " " << K2(1, 2) << "\n"
        << K2(2, 0) << " " << K2(2, 1) << " " << K2(2, 2) << std::endl;

    outFile_1.close();
    outFile_2.close();

    Eigen::Matrix4d CP1, CP2;
    CP1.setIdentity();
    CP2.setIdentity();
    CP1.block<3, 3>(0, 0) = pose_l;
    CP2.block<3, 3>(0, 0) = pose_r;
    CP1.block<3, 1>(0, 3) = c_l;
    CP2.block<3, 1>(0, 3) = c_r;

	std::ofstream cpFile_1(outFolder + leftCamPosefileName), cpFile_2(outFolder + rightCamPosefileName);
    cpFile_1 
        << CP1(0, 0) << " " << CP1(0, 1) << " " << CP1(0, 2) << " " << CP1(0, 3) << "\n"
        << CP1(1, 0) << " " << CP1(1, 1) << " " << CP1(1, 2) << " " << CP1(1, 3) << "\n"
        << CP1(2, 0) << " " << CP1(2, 1) << " " << CP1(2, 2) << " " << CP1(2, 3) << "\n"
        << CP1(3, 0) << " " << CP1(3, 1) << " " << CP1(3, 2) << " " << CP1(3, 3) << "\n"
        << std::endl;

    cpFile_2
        << CP2(0, 0) << " " << CP2(0, 1) << " " << CP2(0, 2) << " " << CP2(0, 3) << "\n"
        << CP2(1, 0) << " " << CP2(1, 1) << " " << CP2(1, 2) << " " << CP2(1, 3) << "\n"
        << CP2(2, 0) << " " << CP2(2, 1) << " " << CP2(2, 2) << " " << CP2(2, 3) << "\n"
        << CP2(3, 0) << " " << CP2(3, 1) << " " << CP2(3, 2) << " " << CP2(3, 3) << "\n"
        << std::endl;

    cpFile_1.close();
    cpFile_2.close();
}
void StereoCamCalibr::fine_projection_matrix(Eigen::Matrix3d* k_left, Eigen::Matrix3d* k_right, Eigen::Matrix4d* cam_pos_left, Eigen::Matrix4d* cam_pos_right)
{
	// Image points in terms of normalized coordinates
	int numPts = leftCamPtr_->totalNumPoints();
	leftCamPtr_->convert_image_pts_to_normalized_coord();
	rightCamPtr_->convert_image_pts_to_normalized_coord();

	double* np = leftCamPtr_->normalized_image_pts_;
	double* nq = rightCamPtr_->normalized_image_pts_;

	// Camera matrix computed
	Eigen::Map<Eigen::Matrix<double, 3, 4>> P_l(P_l_);
	Eigen::Map<Eigen::Matrix<double, 3, 4>> P_r(P_r_);

	//NOTE : It is always possible to set P_l = [ I | 0 ]
	P_l.block<3, 3>(0, 0) = Eigen::Matrix3d::Identity();
	P_l.block<3, 1>(0, 3) = Eigen::Vector3d::Zero();

	double* K_l = leftCamPtr_->K_;
	double* K_r = rightCamPtr_->K_;
	compute_camera_matrix_from_essential_matrix(numPts, np, nq,
		K_l, K_r, E_, P_l.data(), P_r.data());

	// Compute camera pose and center
	Eigen::Matrix3d R_l, R_r;
	R_l = P_l.block<3, 3>(0, 0);
	R_r = P_r.block<3, 3>(0, 0);

	Eigen::Map<Eigen::Matrix3d> pose_l(campose_l_), pose_r(campose_r_);
	pose_l = R_l.transpose();
	pose_r = R_r.transpose();

	Eigen::Map<Eigen::Vector3d> c_l(camcenter_l_), c_r(camcenter_r_);
	c_l = -pose_l * P_l.block<3, 1>(0, 3);
	c_r = -pose_r * P_r.block<3, 1>(0, 3);

	//
	std::cout << "- CAMERA POSE :\n"
		<< "  left :\n" << pose_l << "\n"
		<< "  right : \n" << pose_r << "\n"
		<< "- CAMERA CENTER :\n"
		<< "  left : " << c_l.transpose() << "\n"
		<< "  right : " << c_r.transpose() << "\n" << std::endl;

	// Export camera parameter matrices to disk
	

	Eigen::Matrix3d K1(K_l), K2(K_r);
	for (int i = 0; i < 9; ++i)
	{
		(*k_left)(i) = K1(i);
		(*k_right)(i) = K2(i);
	}

	cam_pos_left->setIdentity();
	cam_pos_right->setIdentity();

	cam_pos_left->block<3, 3>(0, 0) = pose_l;
	cam_pos_right->block<3, 3>(0, 0) = pose_r;
	cam_pos_left->block<3, 1>(0, 3) = c_l;
	cam_pos_right->block<3, 1>(0, 3) = c_r;

}

void compute_rigid_transform(int numPts, double* src, double* tgt, double* Trans)
{
	Eigen::Vector3d* srcPt = (Eigen::Vector3d*)src;
	Eigen::Vector3d* tgtPt = (Eigen::Vector3d*)tgt;
	Eigen::Map<Eigen::Matrix4d> T(Trans);

	if (numPts <= 0)  return;

	// Find centroids for source and target poits
	Eigen::Vector3d tp(0.0, 0.0, 0.0), tq(0.0, 0.0, 0.0);

	for (int i = 0; i < numPts; ++i)
	{
		tp += srcPt[i];
		tq += tgtPt[i];
	}

	tp /= (double)numPts;
	tq /= (double)numPts;
	//std::cout << "Centroids of points : \n" << "src :\n" << tp << "\n" << "tgt :\n" << tq << std::endl;

	// Prepare correspondences with centroids at the origin
	std::vector< Eigen::Vector3d > p(numPts), q(numPts);
	for (int i = 0; i < numPts; ++i)
	{
		p[i] = srcPt[i] - tp;
		q[i] = tgtPt[i] - tq;
	}

	//-------------------------------------------------------------------------------
	// Maximize q^T(N)q  which is equivalent to finding the maximum eigenvector of N
	//-------------------------------------------------------------------------------
	double S[3][3];
	std::memset(S, 0, sizeof(double) * 9);

	for (int y = 0; y < 3; ++y)
	{
		for (int x = 0; x < 3; ++x)
		{
			for (int i = 0; i < numPts; ++i)
			{
				S[x][y] += p[i](x) * q[i](y);
			}
		}
	}

	Eigen::Matrix4d N;

	N(0, 0) = S[0][0] + S[1][1] + S[2][2];
	N(0, 1) = S[1][2] - S[2][1];
	N(0, 2) = S[2][0] - S[0][2];
	N(0, 3) = S[0][1] - S[1][0];

	N(1, 1) = S[0][0] - S[1][1] - S[2][2];
	N(1, 2) = S[0][1] + S[1][0];
	N(1, 3) = S[2][0] + S[0][2];

	N(2, 2) = -S[0][0] + S[1][1] - S[2][2];
	N(2, 3) = S[1][2] + S[2][1];

	N(3, 3) = -S[0][0] - S[1][1] + S[2][2];

	N(1, 0) = N(0, 1);
	N(2, 0) = N(0, 2);
	N(2, 1) = N(1, 2);
	N(3, 0) = N(0, 3);
	N(3, 1) = N(1, 3);
	N(3, 2) = N(2, 3);

	// Compute eigenvalues and eigenvectors of N
	// NOTE : N has four real eigenvalues since N is symmetric.
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix4d> eigenSolver(N);
	//std::cout << "Eigenvalues :\n" << eigenSolver.eigenvalues() << std::endl;
	//std::cout << "Eigenvectors : \n" << eigenSolver.eigenvectors() << std::endl;

	// Choose the largest eigenvector(last index)
	Eigen::Vector4d quat = eigenSolver.eigenvectors().col(3);
	quat.normalize();

	Eigen::Matrix3d R;
	quaternionToRotation(quat.data(), R.data());
	//std::cout << "R=\n" << R << std::endl;

	Eigen::Vector3d t;
	t = tq - R * tp;
	//std::cout << "t=\n" << t << std::endl;

	T.setIdentity();
	T.block<3, 3>(0, 0) = R;
	T.block<3, 1>(0, 3) = t;
}

void find_2D_affine_transform(int n, double* p, double* q, double* A, double* b)
{
	if (n <= 0)
	{
		std::cout << "Error, the number of points should be positive. Nothings done!" << std::endl;
		return;
	}
	//std::cout << "n=" << n << std::endl;

	std::vector<Eigen::Vector2d> srcPt(n), tgtPt(n);
	memcpy(srcPt.data(), p, sizeof(Eigen::Vector2d)*n);
	memcpy(tgtPt.data(), q, sizeof(Eigen::Vector2d)*n);

	Eigen::Vector2d src_center(srcPt[0]), tgt_center(tgtPt[0]);

	for (int i = 0; i < n; ++i)
	{
		srcPt[i] -= src_center;
		tgtPt[i] -= tgt_center;
	}

	Eigen::MatrixXd P(n, 2), Q(n, 2);

	for (int i = 0; i < n; ++i)
	{
		P(i, 0) = srcPt[i](0);
		P(i, 1) = srcPt[i](1);

		Q(i, 0) = tgtPt[i](0);
		Q(i, 1) = tgtPt[i](1);
	}

	Eigen::MatrixXd tP = P.transpose();
	Eigen::Matrix2d M = tP * P;
	if (fabs(M.determinant()) < 1.0e-12)
	{
		std::cout << "CAUTION! Almost zero determinant" << std::endl;
		for (int i = 0; i < n; ++i)
		{
			std::cout << srcPt[i].transpose() << " - " << tgtPt[i].transpose() << "\n";
		}
		std::cout << std::endl;
	}

	Eigen::Matrix2d invM = M.inverse();

	Eigen::MatrixXd B = invM * tP;

	Eigen::Map<Eigen::Matrix2d> tA(A);
	tA = (B*Q).transpose();

	Eigen::Map<Eigen::Vector2d> tb(b);
	tb = -tA * src_center + tgt_center;
}

void find_affine_transform(int n, double* p, double* q, double* A, double* b)
{
	if (n <= 0)
	{
		std::cout << "Error, the number of points should be positive. Nothings done!" << std::endl;
		return;
	}
	//std::cout << "n=" << n << std::endl;

	std::vector<Eigen::Vector3d> srcPt(n), tgtPt(n);
	memcpy(srcPt.data(), p, sizeof(Eigen::Vector3d)*n);
	memcpy(tgtPt.data(), q, sizeof(Eigen::Vector3d)*n);

	Eigen::Vector3d src_center(srcPt[0]), tgt_center(tgtPt[0]);

	for (int i = 0; i < n; ++i)
	{
		srcPt[i] -= src_center;
		tgtPt[i] -= tgt_center;
	}

	Eigen::MatrixXd P(n, 3), Q(n, 3);

	for (int i = 0; i < n; ++i)
	{
		P(i, 0) = srcPt[i](0);
		P(i, 1) = srcPt[i](1);
		P(i, 2) = srcPt[i](2);

		Q(i, 0) = tgtPt[i](0);
		Q(i, 1) = tgtPt[i](1);
		Q(i, 2) = tgtPt[i](2);
	}

	Eigen::MatrixXd tP = P.transpose();
	Eigen::Matrix3d M = tP*P;
	if (fabs(M.determinant()) < 1.0e-12)
	{
		std::cout << "CAUTION! Almost zero determinant" << std::endl;
		for (int i = 0; i < n; ++i)
		{
			std::cout << srcPt[i].transpose() << " - " << tgtPt[i].transpose() << "\n";
		}
		std::cout << std::endl;
	}
	
	Eigen::Matrix3d invM = M.inverse();

	Eigen::MatrixXd B = invM * tP;

	Eigen::Map<Eigen::Matrix3d> tA(A);
	tA = (B*Q).transpose();

	Eigen::Map<Eigen::Vector3d> tb(b);
	tb = -tA * src_center + tgt_center;
}

void find_similarity_transform(int n, double* src_pt, double* tgt_pt, double& scale, double* scale_trans, double* rigid_trans,
    double* similarity_trans)
{
    if (n <= 0)
    {
        std::cout << "Error, the number of points should be positive. Nothings done!" << std::endl;
        return;
    }
    std::cout << "n=" << n << std::endl;

    Eigen::Vector3d* srcPt = (Eigen::Vector3d*)src_pt;
    Eigen::Vector3d* tgtPt = (Eigen::Vector3d*)tgt_pt;

    Eigen::Vector3d src_center(0.0, 0.0, 0.0), tgt_center(0.0, 0.0, 0.0);

    for (int i = 0; i < n; ++i)
    {
        src_center += srcPt[i];
        tgt_center += tgtPt[i];
    }

    src_center /= (double)n;
    tgt_center /= (double)n;
    std::cout << "center[0]=" << src_center.transpose() << std::endl;
    std::cout << "center[1]=" << tgt_center.transpose() << std::endl;

    // Compute the covariance matrices
    Eigen::Matrix3d C[2];
    compute_covariance_matrix(n, (double*)srcPt, src_center.data(), C[0].data());
    compute_covariance_matrix(n, (double*)tgtPt, tgt_center.data(), C[1].data());

    // Find the scaling factor from covariance matrices
    Eigen::JacobiSVD<Eigen::MatrixXd> svd1(C[0], Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::JacobiSVD<Eigen::MatrixXd> svd2(C[1], Eigen::ComputeFullU | Eigen::ComputeFullV);

    Eigen::Vector3d sv[2];
    sv[0] = svd1.singularValues();
    sv[1] = svd2.singularValues();
    std::cout << "sv[0] = " << sv[0].transpose() << std::endl;
    std::cout << "sv[1] = " << sv[1].transpose() << std::endl;

    scale = sqrt( sv[0].dot(sv[1]) / sv[0].squaredNorm() );
    std::cout << "scale factor = " << scale << std::endl;

    Eigen::Map<Eigen::Matrix4d> S(scale_trans);
    S.setIdentity();
    S.block<3, 3>(0, 0) *= scale;
    //S.block<3, 1>(0, 3) *= (1.0 - scale);

    std::vector<Eigen::Vector3d> scaled_src_pt(n);
    for (int i = 0; i < n; ++i)
    {
        scaled_src_pt[i] = scale*srcPt[i];//scale*(srcPt[i] - src_center)+src_center;
        std::cout << i << " : " << scaled_src_pt[i].transpose() << std::endl;
    }

    Eigen::Map<Eigen::Matrix4d> T(rigid_trans);
    compute_rigid_transform(n, scaled_src_pt.data()->data(), (double*)tgtPt, T.data());

    Eigen::Map<Eigen::Matrix4d> finalTrans(similarity_trans);
    finalTrans = T*S;
}//*/

//Umeyama approach
void find_similarity_transform_umeyama(int n, double* src_pt, double* tgt_pt, double& scale, double* rigid_trans,
    double* similarity_trans)
{

    if (n <= 0)
    {
        std::cout << "Error, the number of points should be positive. Nothings done!" << std::endl;
        return;
    }
    //std::cout << "n=" << n << std::endl;

    Eigen::Vector3d* srcPt = (Eigen::Vector3d*)src_pt;
    Eigen::Vector3d* tgtPt = (Eigen::Vector3d*)tgt_pt;

    Eigen::Vector3d src_center(0.0, 0.0, 0.0), tgt_center(0.0, 0.0, 0.0);

    for (int i = 0; i < n; ++i)
    {
        src_center += srcPt[i];
        tgt_center += tgtPt[i];
    }

    src_center /= (double)n;
    tgt_center /= (double)n;
    //std::cout << "center[0]=" << src_center.transpose() << std::endl;
    //std::cout << "center[1]=" << tgt_center.transpose() << std::endl;

    double sigma2[2] = { 0.0, 0.0 };

    Eigen::Matrix3d Sigma_xy;
    Sigma_xy.setZero();

    for (int i = 0; i < n; ++i)
    {
        sigma2[0] += (srcPt[i] - src_center).squaredNorm();
        sigma2[1] += (tgtPt[i] - tgt_center).squaredNorm();

        Sigma_xy += (tgtPt[i] - tgt_center)*(srcPt[i] - src_center).transpose();
    }

    sigma2[0] /= (double)n;
    sigma2[1] /= (double)n;

    Sigma_xy /= (double)n;

    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Sigma_xy, Eigen::ComputeFullU | Eigen::ComputeFullV);

    Eigen::Matrix3d U = svd.matrixU();
    Eigen::Matrix3d V = svd.matrixV();

    Eigen::Vector3d sv = svd.singularValues();
    //std::cout << "singular values = " << sv.transpose() << std::endl;
    //std::cout << "det(Sigma_xy) = " << Sigma_xy.determinant() << std::endl;
    //std::cout << "sigma2_x=" << sigma2[0] << std::endl;

    Eigen::Matrix3d D;
    D.setZero();
    D(0, 0) = sv(0), D(1, 1) = sv(1), D(2, 2) = sv(2);

    int rank = 3;
    double eps = 1.0E-9;

    for (int i = 0; i < 3; ++i)
    {
        if (sv(i) < eps)
            rank--;
    }

    if (rank < 2)
    {
        std::cout << "Error, the covariance matrix has rank <=1." << "\n"
            << "Impossible to determine a proper similarity transform."
            << std::endl;
        return;
    }

    Eigen::Matrix3d S;
    S.setIdentity();

    double det_sigma_xy = Sigma_xy.determinant();
    if (det_sigma_xy < 0.0)
    {
        S(2, 2) = -1.0;
    }

    if (rank == 2)
    {
        S.setIdentity();

        if (U.determinant()*V.determinant() < 0.0)
            S(2, 2) = -1.0;
    }

    Eigen::Matrix3d R;
    R = U*S*V.transpose();

    double k = (D*S).trace() / sigma2[0];
    scale = k;

    Eigen::Vector3d t;
    t = tgt_center - k*R*src_center;

    Eigen::Map<Eigen::Matrix4d> rigidTrans(rigid_trans);
    rigidTrans.setIdentity();
    rigidTrans.block<3, 3>(0, 0) = R;
    rigidTrans.block<3, 1>(0, 3) = t;

    Eigen::Map<Eigen::Matrix4d> SimTrans(similarity_trans);
    SimTrans.setIdentity();
    SimTrans.block<3, 3>(0, 0) = k*R;
    SimTrans.block<3, 1>(0, 3) = t;

    double err = 0.0;
    for (int i = 0; i < n; ++i)
    {
        err += (tgtPt[i] - (k*R*srcPt[i] + t)).squaredNorm();
    }
    err /= n;
    std::cout << "error=" << err << std::endl;
}

void find_similarity_transform_3(int n, double* src_pt, double* tgt_pt, double& scale, double* scale_trans, double* rigid_trans,
    double* similarity_trans)
{
    if (n <= 0)
    {
        std::cout << "Error, the number of points should be positive. Nothings done!" << std::endl;
        return;
    }
    std::cout << "n=" << n << std::endl;

    Eigen::Vector3d* srcPt = (Eigen::Vector3d*)src_pt;
    Eigen::Vector3d* tgtPt = (Eigen::Vector3d*)tgt_pt;

    Eigen::Vector3d src_center(0.0, 0.0, 0.0), tgt_center(0.0, 0.0, 0.0);

    for (int i = 0; i < n; ++i)
    {
        src_center += srcPt[i];
        tgt_center += tgtPt[i];
    }

    src_center /= (double)n;
    tgt_center /= (double)n;
    std::cout << "center[0]=" << src_center.transpose() << std::endl;
    std::cout << "center[1]=" << tgt_center.transpose() << std::endl;

    double totalSumEdgeLength[2] = { 0.0, 0.0 };
    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            totalSumEdgeLength[0] += (srcPt[i] - srcPt[j]).squaredNorm();
            totalSumEdgeLength[1] += (tgtPt[i] - tgtPt[j]).squaredNorm();
        }
    }

    scale = sqrt(totalSumEdgeLength[1] / totalSumEdgeLength[0]);

    Eigen::Map<Eigen::Matrix4d> S(scale_trans);
    S.setIdentity();
    S.block<3, 3>(0, 0) *= scale;
    //S.block<3, 1>(0, 3) *= (1.0 - scale);

    std::vector<Eigen::Vector3d> scaled_src_pt(n);
    for (int i = 0; i < n; ++i)
    {
        scaled_src_pt[i] = scale*srcPt[i];//scale*(srcPt[i] - src_center)+src_center;
        std::cout << i << " : " << scaled_src_pt[i].transpose() << std::endl;
    }

    Eigen::Map<Eigen::Matrix4d> T(rigid_trans);
    compute_rigid_transform(n, scaled_src_pt.data()->data(), (double*)tgtPt, T.data());

    Eigen::Map<Eigen::Matrix4d> finalTrans(similarity_trans);
    finalTrans = T*S;
}