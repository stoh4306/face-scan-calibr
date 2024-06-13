#ifndef MG_CALIBRATION_H_
#define MG_CALIBRATION_H_

#include <vector>
#include <string>
#include <Eigen/Dense>

bool find_check_corners(int w, int h, unsigned char* color_check_img, int board_width, int board_height, double* corner_pts, bool showImage=false);
void find_homography_plane_to_image(int n, double* plane_pts, double* image_pts, double* H_colmaj);
void compute_iac_from_plane2image_homography(int numHom, double* H, double* omega);
void find_calibration_matrix_from_iac(double* tw, double* tK, int method = 0);

void find_fundamental_matrix_normalized_linear(int n, double* left_img_pts, double* right_img_pts, double* F);
void transform_image_pts_to_normalized_coord_frame(int n, double* img_pts, double* K, double* n_img_pts);
void compute_camera_matrix_from_essential_matrix(int n, double* np, double* nq,
    double* K1, double* K2, double* tE, double* tP1, double* tP2);

void compute_rigid_transform(int numPts, double* srcPt, double* tgtPt, double* T);
void find_affine_transform(int n, double* p, double* q, double* A, double* b);
void find_2D_affine_transform(int n, double* p, double* q, double* A, double* b);

void find_similarity_transform(int n, double* src_pt, double* tgt_pt, double& scale, double* scale_trans, double* rigid_trans,
    double* similarity_trans);
void find_similarity_transform_umeyama(int n, double* src_pt, double* tgt_pt, double& scale, double* rigid_trans,
    double* similarity_trans);
void find_similarity_transform_3(int n, double* src_pt, double* tgt_pt, double& scale, double* scale_trans, double* rigid_trans,
    double* similarity_trans);
// Intrinsic parameter extraction from check images
//void find_intrinsic_param_from_check_images(int num_images, int w, int h, unsigned char* color_check_img,
//    int board_width, int board_height, double square_size, double* K);

bool read_calibrations(const char* K_file, const char* CP_file, float* K, float* CP);

class SingleCameraCalibr
{
public:
    SingleCameraCalibr() : num_images_(0), total_num_pts_(0), num_pts_(NULL), plane_pts_(NULL), image_pts_(NULL), 
        normalized_image_pts_(NULL), calibrated_(false) {}
    ~SingleCameraCalibr();

    void setNumImages(int n)    { num_images_ = n; }
    void setNumPointsInImage(int* npts) { num_pts_ = npts; }
    void setImagePoints(double* p)  { image_pts_ = p; }
    void set3DPoints(double* p)     { plane_pts_ = p; }
    void setPoints(int n, int* npts, double* p3d, double* image_pts) 
    { 
        num_images_ = n, num_pts_ = npts, plane_pts_ = p3d, image_pts_ = image_pts;
    }
    int computeTotalNumPoints();
    int totalNumPoints()    { return total_num_pts_; }

    void compute_homography_planes_to_images();
    void compute_image_absolute_conics();
    void find_intrinsic_matrix();
    void set_calibration_status(bool status)  { calibrated_ = status; }

    void convert_image_pts_to_normalized_coord();

public:
    int         num_images_;
    int         total_num_pts_;
    int*        num_pts_;   // num_pts[i] : number of points in image[i]
    double*     plane_pts_;    // Points should be lie on a plane and they have 2d coordinates
    double*     image_pts_;

    // The following array will be allocated within this class if necessary. 
    // So, a memory free is needed at the destruction
    double*     normalized_image_pts_; 

    std::vector<double> H_;     // homography from plane to image
    double Omega_[9]; // image of absolute conic
    double K_[9];               // intrinsic matrix

    bool calibrated_;
};

class CameraCalibrWithCheck : public SingleCameraCalibr
{
public:
    CameraCalibrWithCheck() : SingleCameraCalibr(),
        w_(0), h_(0), check_img_data_(NULL), board_width_(0), board_height_(0) {}

    ~CameraCalibrWithCheck();

    void setImages(int n, int w, int h, unsigned char* check_img_data)
    {
        num_images_ = n;

        w_ = w;
        h_ = h;

        check_img_data_ = check_img_data; 
    }

    void setCheckBoard(int board_width, int board_height, double square_size)
    {
        board_width_ = board_width;
        board_height_ = board_height;
        board_square_size_ = square_size;
    }

    int find_corner_points(const char* filename, bool vis_corner_pts=false);
	int find_corner_points(std::vector<Eigen::Vector2d>* corner_pts, bool vis_corner_pts = false);

public:
    int w_, h_;
    unsigned char* check_img_data_;     // color image
    int board_width_, board_height_;
    double board_square_size_;

    std::vector<double> corner_pts_;
};

class StereoCamCalibr
{
    typedef CameraCalibrWithCheck* CamCalibrPtr;

public:
    StereoCamCalibr() : leftCamPtr_(NULL), rightCamPtr_(NULL),outFolder("") {}

    void setCam(CamCalibrPtr left, CamCalibrPtr right)  { leftCamPtr_ = left, rightCamPtr_ = right; }
    void setLeftCam(CamCalibrPtr left)    { leftCamPtr_ = left; }
    void setRightCam(CamCalibrPtr right)  { rightCamPtr_ = right; }

    void find_fundamental_matrix();
    void find_essential_matrix();
    //void find_projection_marix();
	void find_projection_marix(std::string leftKfileName = "left_k.txt"
		, std::string rightKfileName = "right_k.txt",
		std::string leftCamPosefileName = "left_cam_pose.txt",
		std::string rightCamPosefileName = "right_cam_pose.txt");
	void fine_projection_matrix(Eigen::Matrix3d* k_left, Eigen::Matrix3d* k_right, Eigen::Matrix4d* cam_pos_left, Eigen::Matrix4d* cam_pos_right);

    void reconstruct_3d_points(std::string cornerPointfileName = "corner_pts_0.txt");
	void reconstruct_3d_points(std::vector<Eigen::Vector3d>* corner_pts);
    void reconstruct_3d_points_geom_corr();
    void set_output_folder(std::string folder) { outFolder = folder; }
public:
    CamCalibrPtr leftCamPtr_, rightCamPtr_;
    std::string outFolder;

    double F_[9];
    double E_[9];
    double P_l_[12], P_r_[12];
    double campose_l_[9], campose_r_[9];
    double camcenter_l_[3], camcenter_r_[3];
};
#endif