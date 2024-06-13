#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "numer_math.h"
#include "calibr.h"
#include "reconst.h"
#include "vec.h"

#include <Eigen/Dense>
#include <opencv2/opencv.hpp>

typedef std::vector< mg::Vector2d > Point2dArray;

bool read_check_corner_points(const char* inputFile, int& numImages, int& width, int& height,
    int& board_width, int& board_height, float& d, std::vector< Point2dArray >& projCorner)
{
    std::ifstream inFile(inputFile);

    if (!inFile.is_open())
    {
        std::cout << "Error, can't open the input data : " << inputFile << std::endl;
        return false;
    }

    std::string temp;
    inFile >> temp;

    do
    {
        if (temp == "NUM_IMAGES")
        {
            inFile >> temp >> numImages;
            std::getline(inFile, temp, '\n');
        }

        else if (temp == "IMAGE_WIDTH")
        {
            inFile >> temp >> width;
            std::getline(inFile, temp, '\n');
        }

        else if (temp == "IMAGE_HEIGHT")
        {
            inFile >> temp >> height;
            std::getline(inFile, temp, '\n');
        }

        else if (temp == "BOARD_WIDTH")
        {
            inFile >> temp >> board_width;
            std::getline(inFile, temp, '\n');
        }

        else if (temp == "BOARD_HEIGHT")
        {
            inFile >> temp >> board_height;
            std::getline(inFile, temp, '\n');
        }

        else if (temp == "SQUARE_SIZE")
        {
            inFile >> temp >> d;
            std::getline(inFile, temp, '\n');
        }

        else if (temp == "CORNER_DATA")
        {
            inFile >> temp;

            projCorner.resize(numImages);

            const int numCorners = board_width*board_height;
            for (int i = 0; i < numImages; ++i)
            {
                projCorner[i].resize(numCorners);

                for (int j = 0; j < numCorners; ++j)
                {
                    inFile >> projCorner[i][j].x >> projCorner[i][j].y;
                }
            }//*/
        }
        inFile >> temp;
    } while (!inFile.eof());

    inFile.close();

    // print 
    std::cout << ". NUM_IMAGES = " << numImages << "\n"
        << ". IMAGE_SIZE = " << width << " x " << height << "\n"
        << ". BOARD_SIZE = " << board_width << " x " << board_height << "\n"
        << ". " << projCorner[0][0] << ", " << projCorner[numImages - 1][board_width*board_height - 1]
        << std::endl;

    return true;
}

void set_corner_coordinates_on_board(int board_width, int board_height, float square_size, Point2dArray& corner)
{
    corner.resize(board_width*board_height);

    float y = 0.0f;
    float dy = square_size;

    for (int j = 0; j < board_height; ++j)
    {
        float x = 0.0f;
        float dx = square_size;

        for (int i = 0; i < board_width; ++i)
        {
            corner[j*board_width + i] = mg::Vector2d(x, y);
            x += dx;
        }

        y += dy;
    }

    //std::cout << corner[0] << ", " << corner[board_width*board_height - 1] << std::endl;
    //int index = 0;
    //for (int j = 0; j < board_height; ++j)
    //{
    //    for (int i = 0; i < board_width; ++i)
    //    {
    //        std::cout << corner[index] << " ";
    //        index++;
    //    }
    //    std::cout << "\n";
    //}
    //std::cout << std::endl;
}

void read_intrinsic_param(const char* filename, Eigen::Matrix3d& K, const char* mess)
{
    std::ifstream inFile(filename);
    if (!inFile.is_open())
    {
        std::cout << "Error, can't open the intrinsic parameter data file : " << filename << std::endl;
        return;
    }

    inFile >> K(0, 0) >> K(0, 1) >> K(0, 2)
        >> K(1, 0) >> K(1, 1) >> K(1, 2)
        >> K(2, 0) >> K(2, 1) >> K(2, 2);

    inFile.close();

    std::cout << mess << " : " << "\n" << K << std::endl;

    return;
}

void convert_image_pts_to_normalized_coord(std::vector< Point2dArray >& x, Eigen::Matrix3d K,
    std::vector< Eigen::Vector3d >& nx)
{
    int totalNumPts = 0;
    int numImages = (int)x.size();
    for (int i = 0; i < numImages; ++i)
    {
        totalNumPts += (int)x[i].size();
    }

    nx.resize(totalNumPts);

    Eigen::Matrix3d invK = K.inverse();

    int index = 0;
    for (int i = 0; i < numImages; ++i)
    {
        for (int j = 0; j < (int)x[i].size(); ++j)
        {
            Eigen::Vector3d tx(x[i][j].x, x[i][j].y, 1.0);

            nx[index] = invK*tx;

            index++;
        }
    }

    /*if (x[0][0].x*nx[0](0) < 0.0)
    {
    for (int i = 0; i < totalNumPts; ++i)
    {
    nx[i] *= -1.0;
    }
    }

    std::cout << x[0][0] << "->" << nx[0].transpose() << std::endl;*/
}

int main_1(int argc, char** argv)
{
    //=====================================================
    // Singular value decomposition test
    //=====================================================
    std::cout
        << "************************************************" << "\n"
        << " SINGULAR VALUE DECOMPOSITION" << "\n"
        << "************************************************" << std::endl;

    std::cout << "1. Array inputs" << std::endl;
    double A[9]; // column-major
    A[0] = 0.00412804, A[3] = 0.537757,  A[6] = -0.0535764;
    A[1] = 0.476805,   A[4] = 0.0213649, A[7] = -11.4263;
    A[2] = 0.0389476,  A[5] = 11.3104,   A[8] = 0.020639;

    Eigen::Map<Eigen::Matrix3d> AA(A);
    std::cout << "- Input : A = \n" << AA << std::endl;

    double sv[3];
    double U[9];
    double V[9];

    svd_dp(3, 3, A, sv, U, V);

    bool col_major = true;
    std::cout << "- Output : S, U, V" << std::endl;
    print_matrix("S=", 1, 3, sv, col_major);
    print_matrix("U=", 3, 3, U, col_major);
    print_matrix("V=", 3, 3, V, col_major);

    std::cout << "\n2. Matrix(Eigen) inputs : A" << std::endl;
    Eigen::Map<Eigen::MatrixXd> tA(A, 3, 3);
    Eigen::VectorXd tsv;
    Eigen::MatrixXd tU, tV;
    svd_dp(&tA, &tsv, &tU, &tV);
    std::cout << "tU=\n" << tU << std::endl;
    std::cout << "tV=\n" << tV << std::endl;

    std::cout << "\n3. General rectangular matrix input" << std::endl;
    Eigen::MatrixXd B(4, 5);
    B.setZero();
    B(0, 0) = 1.0;
    B(0, 4) = 2.0;
    B(1, 2) = 3.0;
    B(3, 1) = 4.0;
    std::cout << "- Input : B=\n" << B << std::endl;

    bool thin = false;
    svd_dp(&B, &tsv, &tU, &tV, thin);
    std::cout << "- Output : " << "\n";
    std::cout << "tU=\n" << tU << std::endl;
    std::cout << "tV=\n" << tV << std::endl;

    std::cout << "\n4. Find the approximate solution for Ax=0" << std::endl;
    std::cout << "- Input : B" << std::endl;

    Eigen::VectorXd x(5);
    double minVal = find_nonzero_kernel_vector(4, 5, B.data(), x.data());
    std::cout << "- Output : \n"
        << "min of Bx = " << minVal << ", x=" << x.transpose() << std::endl;

    std::cout << "- Input : A" << std::endl;
    minVal = find_nonzero_kernel_vector(3, 3, A, x.data());
    std::cout << "- Output : \n"
        << "min of Ax = " << minVal << ", x=" << x.transpose() << std::endl;

    //=====================================================
    // RQ-decomposition test
    //=====================================================
    std::cout
        << "************************************************" << "\n"
        << " RQ DECOMPOSITION" << "\n"
        << "************************************************" << std::endl;

    A[0] = 5.82819e-007, A[3] = 2.89432e-008, A[6] = -0.000315824;
    A[1] = 2.89432e-008, A[4] = 5.96646e-007, A[7] = -0.000247933;
    A[2] = -0.000315824, A[5] = -0.000247933, A[8] =  1.0;

    Eigen::Map<Eigen::Matrix3d> M(A);
    std::cout << "- Input : A=" << "\n" << M << std::endl;

    Eigen::Matrix3d R, Q;
    RQ_decomp_GS_3d(A, R.data(), Q.data());
    std::cout << "- Output : R,Q" << "\n"
        << "R=\n" << R << "\n"
        << "Q=\n" << Q << "\n"
        << "err=" << (M - R*Q).norm() << std::endl;

    //=====================================================
    // Homography, image of absolute conic and intrinsic camera matrix
    //=====================================================
    std::cout
        << "************************************************" << "\n"
        << " HOMOGRAPHY, IMAGE OF ABSOLUTE CONIC" << "\n"
        << " AND INTRINSIC CAMERA MATRIX" << "\n"
        << "************************************************" << std::endl;
    std::vector<double> img_pts(3 * 8);
    img_pts[0] = 155, img_pts[1] = 151;
    img_pts[2] = 482, img_pts[3] = 78;
    img_pts[4] = 217, img_pts[5] = 411;
    img_pts[6] = 490, img_pts[7] = 331;
    img_pts[8] = 595, img_pts[9] = 87;
    img_pts[10] = 895, img_pts[11] = 195;
    img_pts[12] = 596, img_pts[13] = 336;
    img_pts[14] = 837, img_pts[15] = 459;
    img_pts[16] = 490, img_pts[17] = 387;
    img_pts[18] = 778, img_pts[19] = 465;
    img_pts[20] = 343, img_pts[21] = 601;
    img_pts[22] = 689, img_pts[23] = 720;

    std::vector<double> plane_pts(8);
    plane_pts[0] = 0, plane_pts[1] = 0;
    plane_pts[2] = 1, plane_pts[3] = 0;
    plane_pts[4] = 0, plane_pts[5] = 1;
    plane_pts[6] = 1, plane_pts[7] = 1;

    std::vector<Eigen::Matrix3d> H(3);
    find_homography_plane_to_image(4, plane_pts.data(), img_pts.data(), H[0].data());
    find_homography_plane_to_image(4, plane_pts.data(), img_pts.data()+8, H[1].data());
    find_homography_plane_to_image(4, plane_pts.data(), img_pts.data()+16, H[2].data());

    Eigen::Matrix3d omega;
    compute_iac_from_plane2image_homography(3, H.data()->data(), omega.data());
    //std::cout << "iac:\n" << omega << std::endl;

    Eigen::Matrix3d K;
    find_calibration_matrix_from_iac(omega.data(), K.data(), 1);

    //=====================================================
    // Fundamental matrix, essential matrix and 3d reconstruction
    //=====================================================
    if (argc <= 4)
    {
        std::cout << "***** Test has been stopped.\n"
            << "***** At least 4 input arguments are neeeded to proceed" << "\n"
            << "***** Usage : " << argv[0] << " arg1 arg2 arg3 arg4" << "\n"
            << "- arg1 : left image points" << "\n"
            << "- arg2 : right image points" << "\n"
            << "- arg3 : intrinsic camera matrix for left" << "\n"
            << "- arg4 : intrinsic camera matrix for right" << "\n"
            << std::endl;
        return 0;
    }

    std::cout
        << "************************************************" << "\n"
        << " FUNDAMENTAL MATRIX, ESSENTIAL MATRIX" << "\n"
        << " AND 3D RECONSTRUCTION" << "\n"
        << "************************************************" << std::endl;

    const char* left_input_file = argv[1];
    const char* right_input_file = argv[2];
    const char* left_K_file = argv[3];
    const char* right_K_file = argv[4];

    // Load input data
    int numImages = 0;
    int width = 0, height = 0;
    int board_width = 0, board_height = 0;
    float square_size = 0.0f;

    // Read projected corner data
    std::vector< Point2dArray > projCorner[2];
    std::cout
        << "--------------------------------------------------\n"
        << " LEFT IMAGE:\n"
        << "--------------------------------------------------"
        << std::endl;
    read_check_corner_points(left_input_file, numImages, width, height,
        board_width, board_height, square_size, projCorner[0]);

    std::cout
        << "--------------------------------------------------\n"
        << " RIGHT IMAGE:\n"
        << "--------------------------------------------------"
        << std::endl;
    read_check_corner_points(right_input_file, numImages, width, height,
        board_width, board_height, square_size, projCorner[1]);

    // Set the coordinates of the corner points on the board
    Point2dArray corner;
    set_corner_coordinates_on_board(board_width, board_height, square_size, corner);

    // Find the fundamental matrix
    int totalNumPts = 0;
    std::vector< Point2dArray >& left_proj_corner = projCorner[0];
    std::vector< Point2dArray >& right_proj_corner = projCorner[1];
    Point2dArray left_pts, right_pts;
    for (int i = 0; i < left_proj_corner.size(); ++i)
    {
        for (int j = 0; j < left_proj_corner[i].size(); ++j)
        {
            left_pts.push_back(left_proj_corner[i][j]);
            right_pts.push_back(right_proj_corner[i][j]);
            totalNumPts++;
        }
    }

    Eigen::Matrix3d F;
    find_fundamental_matrix_normalized_linear(totalNumPts, (double*)left_pts.data(), (double*)right_pts.data(), F.data());

    // Get the essential matrix from intrinsic parameters
    Eigen::Matrix3d tK[2]; // 0: left, 1: right
    read_intrinsic_param(left_K_file, tK[0], "LEFT_K");
    read_intrinsic_param(right_K_file, tK[1], "RIGHT_K");

    // Image points in terms of normalized coordinates
    std::vector< Eigen::Vector3d > np, nq; // np: left, nq: right
    convert_image_pts_to_normalized_coord(projCorner[0], tK[0], np);
    convert_image_pts_to_normalized_coord(projCorner[1], tK[1], nq);

    Eigen::Matrix3d E;
    E = tK[1].transpose()*(F*tK[0]);
    std::cout << "\n-Essential matrix(E) : \n" << E << std::endl;

    // Camera matrix computed
    Eigen::Matrix<double, 3, 4> P[2];
    P[0].block<3, 3>(0, 0) = Eigen::Matrix3d::Identity();
    P[0].block<3, 1>(0, 3) = Eigen::Vector3d::Zero();

    compute_camera_matrix_from_essential_matrix(totalNumPts, (double*)np.data(), (double*)nq.data(), 
        tK[0].data(), tK[1].data(), E.data(), P[0].data(), P[1].data());

    std::vector<Eigen::Vector3d> X(np.size());
    for (int i = 0; i < np.size(); ++i)
    {
        linear_triangulation(np[i].data(), nq[i].data(), P[0].data(), P[1].data(), X[i].data());
    }

    write_point_cloud_in_ply("reconst.ply", (int)X.size(), X.data()->data());

    return 0;
}

// Single camera calibration(intrinsic parameter extraction from check images)
int main_2(int argc, char** argv)
{
    int numImages = atoi(argv[1]);
    const char* imgPrefix = argv[2];
    int paddings = atoi(argv[3]);
    const char* ext = argv[4];

    int bw = atoi(argv[5]); // board width
    int bh = atoi(argv[6]); // board height
    double bsize = atof(argv[7]); // square size in board

    std::cout << "- LOADING IMAGES : " << imgPrefix << "(" << numImages << ")...";

    std::vector<unsigned char> imgData;
 
    int w, h;
    for (int i = 0; i < numImages; ++i)
    {
        std::ostringstream filename;
        filename << imgPrefix << std::setw(paddings) << std::setfill('0') << i << "." << ext;

        cv::Mat tempImg = cv::imread(filename.str().c_str());

        w = tempImg.cols;
        h = tempImg.rows;

        imgData.insert(imgData.end(), tempImg.data, tempImg.data + 3 * w*h);

        //cv::Mat testImg(h, w, CV_8UC3, imgData.data() + i*(3 * w*h));
        //cv::imshow("check image", testImg);
        //cv::waitKey(0);
    }
    std::cout << "(" << w << "x" << h << ")" << std::endl;

    CameraCalibrWithCheck camCalibr;
    camCalibr.setImages(numImages, w, h, imgData.data());
    camCalibr.setCheckBoard(bw, bh, bsize);

    camCalibr.find_corner_points("corner_img_pts.txt");

    camCalibr.compute_homography_planes_to_images();
    camCalibr.compute_image_absolute_conics();
    camCalibr.find_intrinsic_matrix();

    return 0;
}

// Stereo calibration with check pattern and check corner points reconstruction
int main(int argc, char** argv)
{
    if (argc < 10)
    {
        std::cout << "Usage : calibr_test_1.exe arg[1] ... arg[10]" << "\n"
            << " - arg[1] : number of images" << "\n"
            << " - arg[2] : image prefix" << "\n"
            << " - arg[3] : number of paddings" << "\n"
            << " - arg[4] : image extension" << "\n"
            << " - arg[5] : start frame" << "\n"
            << " - arg[6] : board width" << "\n"
            << " - arg[7] : board height" << "\n"
            << " - arg[8] : board size" << "\n"
            << " - arg[9] : left image folder" << "\n"
            << " - arg[10] : right image folder" << std::endl;
        return -1;
    }
    int numImages = atoi(argv[1]);
    const char* imgPrefix = argv[2];
    int paddings = atoi(argv[3]);
    const char* ext = argv[4];
    int start_frame = atoi(argv[5]);

    int bw = atoi(argv[6]); // board width
    int bh = atoi(argv[7]); // board height
    double bsize = atof(argv[8]); // square size in board

    const char* left_folder = argv[9];
    const char* right_folder = argv[10];

    std::cout << "- LOADING IMAGES : " << imgPrefix << "(" << numImages << ")...";

    std::vector<unsigned char> imgData[2];

    int w, h;
    for (int i = 0; i < numImages; ++i)
    {
        int frame = i + start_frame;
        for (int j = 0; j < 2; ++j)
        {
            std::ostringstream filename;
            if (j==0)
                filename << left_folder << "/" << imgPrefix << std::setw(paddings) << std::setfill('0') << frame << "." << ext;
            else
                filename << right_folder << "/" << imgPrefix << std::setw(paddings) << std::setfill('0') << frame << "." << ext;

            std::cout << filename.str() << std::endl;

            cv::Mat tempImg;
            tempImg = cv::imread(filename.str().c_str());
            if (tempImg.empty())
            {
                std::cout << "Error, can't load image from " << filename.str() << std::endl;
                return 0;
            }

            w = tempImg.cols;
            h = tempImg.rows;

            imgData[j].insert(imgData[j].end(), tempImg.data, tempImg.data + 3 * w*h);

            //cv::Mat testImg(h, w, CV_8UC3, imgData.data() + i*(3 * w*h));
            //cv::imshow("check image", testImg);
            //cv::waitKey(0);
        }
    }
    std::cout << "(" << w << "x" << h << ")" << std::endl;


    // Find intrinsic parameters -> K 
    CameraCalibrWithCheck camCalibr[2];
    std::string corner_pts_filename[2];
    corner_pts_filename[0] = "left_pts_0.txt";
    corner_pts_filename[1] = "right_pts_0.txt";
    bool vis_result = false;
    for (int i = 0; i < 2; ++i)
    {
        //std::cout << "Set images" << std::endl;
        camCalibr[i].setImages(numImages, w, h, imgData[i].data());
        //std::cout << "Set check board" << std::endl;
        camCalibr[i].setCheckBoard(bw, bh, bsize);
        std::cout << "-Finding corner points for camera[" << i << "]" << std::endl;
        camCalibr[i].find_corner_points(corner_pts_filename[i].c_str(), vis_result);
        camCalibr[i].compute_homography_planes_to_images();
        camCalibr[i].compute_image_absolute_conics();
        camCalibr[i].find_intrinsic_matrix();
    }

    // Find the Fundamental matrix for image correspondences
    StereoCamCalibr stereoCalibr;
    stereoCalibr.setCam(&camCalibr[0], &camCalibr[1]);
    stereoCalibr.find_fundamental_matrix();
    stereoCalibr.find_essential_matrix();
    stereoCalibr.find_projection_marix();

    stereoCalibr.reconstruct_3d_points();

    return 0;
}

// Resizing
int main_resize_images(int argc, char** argv)
{
    if (argc < 6)
    {
        std::cout << "Usage : resize_img_seq.exe num_images img_prefix paddings ext start_frame" << std::endl;
        return -1;
    }

    int numImages = atoi(argv[1]);
    const char* imgPrefix = argv[2];
    int paddings = atoi(argv[3]);
    const char* ext = argv[4];
    int start_frame = atoi(argv[5]);

    int w, h;
    for (int i = 0; i < numImages; ++i)
    {
        int frame = i + start_frame;

        std::ostringstream infilename, outfilename;
        infilename << imgPrefix << std::setw(paddings) << std::setfill('0') << frame << "." << ext;
        outfilename << "resized_" << std::setw(paddings) << std::setfill('0') << frame << "." << ext;

        cv::Mat tempImg;
        tempImg = cv::imread(infilename.str().c_str());
        if (tempImg.empty())
        {
            std::cout << "Error, can't load image from " << infilename.str() << std::endl;
            return 0;
        }

        w = tempImg.cols;
        h = tempImg.rows;

        cv::Mat resizedImg;
        cv::resize(tempImg, resizedImg, cv::Size(w/2, h/2));

        cv::imwrite(outfilename.str().c_str(), resizedImg);
        std::cout << infilename.str() << "->" << outfilename.str() << std::endl;
    }

    return 0;
}
