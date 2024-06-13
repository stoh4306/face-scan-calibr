#include "reconst.h"
#include "calibr.h"
#include "numer_math.h"

#include "knn_search.h"

#include <Eigen/Dense>

#include <iostream>
#include <fstream>

#include <opencv2/opencv.hpp>

void warpImage(int w, int h, int type, void* srcImg, void* tgtImg, double* pH)
{
    Eigen::Map<Eigen::Matrix3d> H(pH);
    Eigen::Matrix3f tH = H.cast<float>().transpose();

    cv::Mat cvH(3, 3, CV_32FC1, tH.data());
    //memcpy(cvH.data, tH.data(), sizeof(float) * 9);
    std::cout << "cvH=" << cvH << std::endl;

    cv::Mat src(h, w, type, srcImg);
    cv::Mat dst(h, w, type, tgtImg);
    cv::warpPerspective(src, dst, cvH, src.size());
}

void draw_epipolar_lines(cv::Mat img, Eigen::Vector3d e, cv::Scalar color)
{
    int n = 5;

    int delta_y = img.rows / (n + 1);

    std::cout << "e=" << e.transpose() << std::endl;
    if (abs(e(2) - 1.0) < 1.0e-9)
    {
        //std::cout << "e=" << e.transpose() << std::endl;
        for (int i = 0; i < n; ++i)
        {
            int x0 = 0, y0 = delta_y*(i + 1);

            Eigen::Vector3d p(x0, y0, 1.0);

            double t = (img.cols - 1 - p(0)) / (e(0) - p(0));

            Eigen::Vector3d q = (1.0 - t)*p + t*e;

            std::cout << "line[" << i << "] : " << p.transpose() << "->" << q.transpose() << std::endl;

            cv::line(img, cv::Point((int)p(0), (int)p(1)), cv::Point((int)q(0), (int)q(1)), color);
        }
    }
    else // epipole at inifity
    {
        for (int i = 0; i < n; ++i)
        {
            int x0 = 0, y0 = delta_y*(i + 1);

            Eigen::Vector3d p(x0, y0, 1.0);

            double t = (img.cols - 1) / (e(0) - p(0) + img.cols);

            Eigen::Vector3d q = (1.0 - t)*p + t*e;
            q /= q(2);

            std::cout << "line[" << i << "] : " << p.transpose() << "->" << q.transpose() << std::endl;

            cv::line(img, cv::Point((int)p(0), (int)p(1)), cv::Point((int)q(0), (int)q(1)), color);
        }
    }

    //cv::imshow("epipolar lines", img);
    //cv::waitKey(0);
}

void draw_horizontal_lines(cv::Mat img, int numLines, cv::Scalar color)
{
    int n = numLines;

    int delta_y = img.rows / (n + 1);

    for (int i = 1; i < n+1; ++i)
    {
        cv::line(img, cv::Point(0, i*delta_y), cv::Point(img.cols - 1, i*delta_y), color);
    }
    //cv::imshow("epipolar lines", img);
    //cv::waitKey(0);
}

void rectify(double* tK1, double* tP1, double B, int numPts, double* x2, double* X, 
    double* tK2, double* tP2,
    double* H, int w, int h, int type, void* left_img, void* right_img, void* rect_right_img, std::string folder)
{
    //--------------------------------------------------------
    // Idea and procedure
    //--------------------------------------------------------
    // 1. New projection for P2'=K[I|t'] with t'=[-|t|, 0, 0]
    // 2. New projection pts x2' of X for P2, i.e. x2'=P2*X
    // 3. Compute homography H to send x1 tp x2'
    // 4. Warp the right image with H. 
    //    The warped image must be the rectified image with the left image

    Eigen::Map<Eigen::Matrix3d> K1(tK1);
    Eigen::Map<Eigen::MatrixXd> P1(tP1, 3, 4);

    // 1. Caculate P2
    Eigen::Map<Eigen::Matrix3d> K2(tK2);
    K2 = K1;
    Eigen::Map<Eigen::MatrixXd> P2(tP2, 3, 4);
    P2.block<3, 3>(0, 0).setIdentity();
    Eigen::Vector3d t = Eigen::Vector3d(B, 0, 0);
    P2.block<3, 1>(0, 3) = -t;
    P2 = K2*P2;

    // 2. Project reconstructed point in 3D by P2
    std::vector< Eigen::Vector2d > x3(numPts);
    for (int i = 0; i < numPts; ++i)
    {
        Eigen::Vector4d tp(X[3*i+0], X[3*i+1], X[3*i+2], 1.0);
        Eigen::Vector3d tq;
        tq = P2*tp;

        tq /= tq(2);

        x3[i](0) = tq(0);
        x3[i](1) = tq(1);

        Eigen::Vector2d tx2(x2[2 * i + 0], x2[2 * i + 1]);
        //std::cout << i << " : " << tx2.transpose() << "->" << x3[i].transpose() << std::endl;
    }

    // 3. Compute homography to send x2 and x3
    find_homography_plane_to_image(numPts, x2, x3.data()->data(), H);

    // 4. warp image
    warpImage(w, h, type, right_img, rect_right_img, H);

    cv::Mat rectImg(h, w, type, rect_right_img);
    cv::imwrite(folder+"rect_right.bmp", rectImg);

    
    cv::Mat tempImg[2];

    cv::Mat leftImg(h, w, type, left_img);
    leftImg.copyTo(tempImg[0]);;
    rectImg.copyTo(tempImg[1]);
    // Draw horizontal lines
    draw_horizontal_lines(tempImg[0], 10, cv::Scalar(0, 0, 255));
    draw_horizontal_lines(tempImg[1], 10, cv::Scalar(0, 0, 255));
    cv::imwrite(folder + "left_rect_lines.bmp", tempImg[0]);
	cv::imwrite(folder+ "right_rect_lines.bmp", tempImg[1]);

    // Draw epipolar lines
    //Eigen::Vector3d inf_e1 = P1*Eigen::Vector4d(B, 0, 0, 1);
    //Eigen::Vector3d inf_e2 = P2*Eigen::Vector4d(0, 0, 0, 1);
    //draw_epipolar_lines(tempImg[0], inf_e1, cv::Scalar(0, 0, 255));
    //draw_epipolar_lines(tempImg[1], inf_e2, cv::Scalar(0, 0, 255));
    //cv::imwrite("left_rect_epi_line.bmp", tempImg[0]);
    //cv::imwrite("right_rect_epi_line.bmp", tempImg[1]);

}

void rectify2(double* tK1, double* tP1, double B, int numPts, double* x2,
    double* tK2, double* camPose2, double* tP2,
    double* H, int w, int h, int type, void* left_img, void* right_img, void* rect_right_img,std::string folder)
{
    //--------------------------------------------------------
    // Idea and procedure
    //--------------------------------------------------------
    // 1. New projection for P2'=K[I|t'] with t'=[-|t|, 0, 0]
    // 2. New projection pts x2' of X for P2, i.e. x2'=P2*X
    // 3. Compute homography H to send x1 tp x2'
    // 4. Warp the right image with H. 
    //    The warped image must be the rectified image with the left image

    Eigen::Map<Eigen::Matrix3d> K1(tK1);
    Eigen::Map<Eigen::MatrixXd> P1(tP1, 3, 4);

    // 1. Caculate P2
    Eigen::Map<Eigen::Matrix3d> K2(tK2);
    Eigen::Matrix3d origK2 = K2;
    K2 = K1;
    Eigen::Map<Eigen::MatrixXd> P2(tP2, 3, 4);
    P2.block<3, 3>(0, 0).setIdentity();
    Eigen::Vector3d t = Eigen::Vector3d(B, 0, 0);
    P2.block<3, 1>(0, 3) = -t;
    P2 = K2*P2;

    //// 2. Project reconstructed point in 3D by P2
    //std::vector< Eigen::Vector2d > x3(numPts);
    //for (int i = 0; i < numPts; ++i)
    //{
    //    Eigen::Vector4d tp(X[3 * i + 0], X[3 * i + 1], X[3 * i + 2], 1.0);
    //    Eigen::Vector3d tq;
    //    tq = P2*tp;

    //    tq /= tq(2);

    //    x3[i](0) = tq(0);
    //    x3[i](1) = tq(1);

    //    Eigen::Vector2d tx2(x2[2 * i + 0], x2[2 * i + 1]);
    //    //std::cout << i << " : " << tx2.transpose() << "->" << x3[i].transpose() << std::endl;
    //}

    //// 3. Compute homography to send x2 and x3
    //find_homography_plane_to_image(numPts, x2, x3.data()->data(), H);
    Eigen::Map<Eigen::Matrix3d> tH(H);
    Eigen::Map<Eigen::Matrix3d> R(camPose2);
    tH = K1*R*origK2.inverse();

    // 4. warp image
    warpImage(w, h, type, right_img, rect_right_img, H);

    cv::Mat rectImg(h, w, type, rect_right_img);
    cv::imwrite(folder + "rect_right.bmp", rectImg);

    // Drar epipolar lines
    cv::Mat tempImg[2];

    cv::Mat leftImg(h, w, type, left_img);
    leftImg.copyTo(tempImg[0]);;
    rectImg.copyTo(tempImg[1]);
    Eigen::Vector3d inf_e1 = P1*Eigen::Vector4d(B, 0, 0, 1);
    Eigen::Vector3d inf_e2 = P2*Eigen::Vector4d(0, 0, 0, 1);
    draw_epipolar_lines(tempImg[0], inf_e1, cv::Scalar(0, 0, 255));
    draw_epipolar_lines(tempImg[1], inf_e2, cv::Scalar(0, 0, 255));
    cv::imwrite(folder + "left_rect_epi_line.bmp", tempImg[0]);
    cv::imwrite(folder + "right_rect_epi_line.bmp", tempImg[1]);

}

void stereo_matching(int w, int h, int data_type, void* left_img_data, void* right_img_data,
    int blockSize, int maxDisparity, float* disp, const char* out_disparity_file)
{
    cv::Mat leftImg(h, w, data_type, left_img_data);
    cv::Mat rightImg(h, w, data_type, right_img_data);

	
    // Compute disparities from block matching algorithm
    std::cout << "- Computing correspondences and disparities..." << std::endl;
    cv::Mat dispImg;

    int numberOfDisparities = maxDisparity;
    numberOfDisparities = numberOfDisparities > 0 ? numberOfDisparities : ((w / 8) + 15) & -16;
    int SADWindowSize = blockSize;
    SADWindowSize = SADWindowSize > 0 ? SADWindowSize : 3;
    int speckleWindowSize = 100;
    int speckleRange = 32;
    int cn = leftImg.channels();

#if CV_MAJOR_VERSION < 3
    cv::StereoSGBM sgbm;
    sgbm.preFilterCap = 63;
    //sgbm.SADWindowSize = SADWindowSize > 0 ? SADWindowSize : 3;
    sgbm.SADWindowSize = SADWindowSize;
    sgbm.P1 = 8 * cn*sgbm.SADWindowSize*sgbm.SADWindowSize;
    sgbm.P2 = 32 * cn*sgbm.SADWindowSize*sgbm.SADWindowSize;
    sgbm.minDisparity = 0;// -numberOfDisparities;
    sgbm.numberOfDisparities = numberOfDisparities;
    sgbm.uniquenessRatio = 10;
    sgbm.speckleWindowSize = speckleWindowSize;//bm.state->speckleWindowSize;
    sgbm.speckleRange = speckleRange; //bm.state->speckleRange;
    sgbm.disp12MaxDiff = 1;
    sgbm.fullDP = true; // alg == STEREO_HH;
    sgbm(leftImg, rightImg, dispImg);
#else
    cv::Ptr<cv::StereoSGBM> sgbm = cv::StereoSGBM::create(0, numberOfDisparities, 
        SADWindowSize,
		//3,
		8 * cn*SADWindowSize*SADWindowSize, 32 * cn*SADWindowSize*SADWindowSize,
        1, 63, 10, speckleWindowSize, speckleRange, cv::StereoSGBM::MODE_SGBM);
    sgbm->compute(leftImg, rightImg, dispImg);
#endif

    short* dispImgData = (short*)dispImg.data;
    for (int i = 0; i < w*h; ++i)
    {
        disp[i] = dispImgData[i] / 16.0f;
    }

    cv::Mat disp8;
    dispImg.convertTo(disp8, CV_8U, 255 / (numberOfDisparities*16.));
    cv::imshow("disparity", disp8);
    cv::waitKey(0);

    std::ofstream outFile(out_disparity_file);
    outFile << w << " " << h << "\n";

    for (int i = 0; i < h; ++i)
    {
        for (int j = 0; j < w; ++j)
        {
            float d = disp[i*w + j];
            outFile << d << " ";
        }
        outFile << "\n";
    }
    outFile.close();
}

void stereo_matching(int w, int h, int data_type, void* left_img_data, void* right_img_data,
    int blockSize, int maxDisparity, double* disp, const char* out_disparity_file)
{
    cv::Mat leftImg(h, w, data_type, left_img_data);
    cv::Mat rightImg(h, w, data_type, right_img_data);
	//cv::imshow("left", leftImg);
	//cv::waitKey(0);

	//cv::imshow("right", rightImg);
	//cv::waitKey(0);
    // Compute disparities from block matching algorithm
    std::cout << "- Computing disparities(block matching)..." << std::endl;
	cv::Mat dispImg;

	int numberOfDisparities = maxDisparity;
	numberOfDisparities = numberOfDisparities > 0 ? numberOfDisparities : ((w / 8) + 15) & -16;
	
	int SADWindowSize = blockSize;
	SADWindowSize = SADWindowSize > 0 ? SADWindowSize : 3;
	int speckleWindowSize = 100;
	int speckleRange = 32;
	int cn = leftImg.channels();

	//numberOfDisparities = MAX(maxDisparity,0);

#if CV_MAJOR_VERSION < 3
	cv::StereoSGBM sgbm;
	sgbm.preFilterCap = 63;
	//sgbm.SADWindowSize = SADWindowSize > 0 ? SADWindowSize : 3;
	sgbm.SADWindowSize = SADWindowSize;
	sgbm.P1 = 8 * cn*sgbm.SADWindowSize*sgbm.SADWindowSize;
	sgbm.P2 = 32 * cn*sgbm.SADWindowSize*sgbm.SADWindowSize;
	sgbm.minDisparity = 0;// -numberOfDisparities;
	sgbm.numberOfDisparities = numberOfDisparities;
	sgbm.uniquenessRatio = 10;
	sgbm.speckleWindowSize = speckleWindowSize;//bm.state->speckleWindowSize;
	sgbm.speckleRange = speckleRange; //bm.state->speckleRange;
	sgbm.disp12MaxDiff = 1;
	sgbm.fullDP = 0; // alg == STEREO_HH;
	sgbm(leftImg, rightImg, dispImg);
#else
	int minDisparity = 0;// -maxDisparity / 5;
	//numberOfDisparities = maxDisparity;// minDisparity;

	cv::Ptr<cv::StereoSGBM> sgbm = cv::StereoSGBM::create(
		minDisparity, //0,//-numberOfDisparities/4, //minDisparity
		numberOfDisparities,//numDisparities
		SADWindowSize, //blockSize
		8 * cn*SADWindowSize*SADWindowSize, //P1
		32 * cn*SADWindowSize*SADWindowSize,//P2
		1, //disp12MaxDiff
		127, //63, //preFilterCap
		10,//40, //uniquenessRatio
		speckleWindowSize, //speckleWindowSize
		speckleRange, //speckleRange
		cv::StereoSGBM::MODE_SGBM);//cv::StereoSGBM::MODE_SGBM);//mode

	sgbm->compute(leftImg, rightImg, dispImg);
#endif

	//cv::Mat dispImg2;
	//dispImg.convertTo(dispImg2, CV_8UC1, 1 / 16.0);
	//cv::imshow("disp", dispImg2);
	//cv::waitKey(0);

    short* dispImgData = (short*)dispImg.data;
    for (int i = 0; i < w*h; ++i)
    {
        disp[i] = dispImgData[i] / 16.0;
    }
	
	//cv::waitKey(0);
    //cv::Mat disp8, resizedDisp8;
    //dispImg.convertTo(disp8, CV_8U, 255 / (numberOfDisparities*16.));
	//cv::imshow("disparity", disp8);
	//cv::resize(disp8, resizedDisp8, cv::Size(disp8.cols / 2, disp8.rows / 2));
    //cv::imshow("disparity", resizedDisp8);
    //cv::waitKey(0);//*/

    /*std::ofstream outFile(out_disparity_file);
    outFile << w << " " << h << "\n";

    for (int i = 0; i < h; ++i)
    {
        for (int j = 0; j < w; ++j)
        {
            int d = disp[i*w + j];
            outFile << d << " ";
        }
        outFile << "\n";
    }
    outFile.close();*/
}

void linear_triangulation(double* tx1, double* tx2, double* tP1, double* tP2, double* tX)
{
    Eigen::Vector3d x1(tx1[0], tx1[1], 1.0), x2(tx2[0], tx2[1], 1.0);
    Eigen::Map<Eigen::MatrixXd> P1(tP1, 3, 4), P2(tP2, 3, 4);
    Eigen::Map<Eigen::Vector3d> X(tX);

    Eigen::MatrixXd A(4, 4);

    A.block<1, 4>(0, 0) = x1(0)*P1.row(2) - P1.row(0);
    A.block<1, 4>(1, 0) = x1(1)*P1.row(2) - P1.row(1);

    A.block<1, 4>(2, 0) = x2(0)*P2.row(2) - P2.row(0);
    A.block<1, 4>(3, 0) = x2(1)*P2.row(2) - P2.row(1);

    double ssv;
    Eigen::VectorXd V(4);

    ssv = find_nonzero_kernel_vector(4, 4, A.data(), V.data());//compute_smallest_singular_vector(A, ssv, V);

    X(0) = V(0) / V(3);
    X(1) = V(1) / V(3);
    X(2) = V(2) / V(3);

    Eigen::Vector3d pv1 = P1*V;
    Eigen::Vector3d pv2 = P2*V;

    pv1 /= pv1(2);
    pv2 /= pv2(2);

    /*std::cout << "TRIANGULATION : X = " << X.transpose() 
        << " (proj. err=" << (x1-pv1).norm() << ", " << (x2-pv2).norm() << ") " << std::endl;
    std::cout
        << "x1=" << x1.transpose() << ", P1*V=" << pv1.transpose() << "\n"
        << "x2=" << x2.transpose() << ", P2*V=" << pv2.transpose() << "\n"
        << std::endl;//*/
}

void point2rawdepth(int n, float* p, float* K_data, float* T_data, int width, int height, float* rawDepth, bool holeFilling, const char* filename)
{
	Eigen::Matrix3f K(K_data);
	Eigen::Matrix<float, 3, 4> T = Eigen::Matrix4f(T_data).block<3, 4>(0, 0);

	Eigen::Matrix<float, 3, 4> P = K * T;
    //std::cout << "- Projection : " << P << std::endl;

	std::fill(rawDepth, rawDepth + width * height, 0.0f);

	float INF = std::numeric_limits<float>::infinity();

	cv::Mat depthMap(height, width, CV_16UC1);
	unsigned short * depthData = (unsigned short*)depthMap.data;

	//float zmin = INF, zmax = -INF;

	for (int i = 0; i < n; ++i)
	{
		Eigen::Vector4f x(p[3*i+0], p[3*i+1], p[3*i+2], 1.0f);
		Eigen::Vector3f y = P * x;

		Eigen::Vector2i q;
		q(0) = (int)floor(y(0) / y(2));
		q(1) = (int)floor(y(1) / y(2));

		if( y(2) >= 0.0f && y(2) < INF && q(0) >= 0 && q(0) < width && q(1) >= 0 && q(1) < height)
		{
			rawDepth[q(1)*width + q(0)] = y(2);
		}

		//if (y(2) < zmin)	zmin = y(2);
		//if (y(2) > zmax)	zmax = y(2);
	}

	//std::cout << "z : " << zmin << " ~ " << zmax << std::endl;

	if (holeFilling)
	{
		std::vector<float> tempRawDepth(width*height);
		std::copy(rawDepth, rawDepth + width * height, tempRawDepth.data());

		for (int j = 1; j < height - 1; ++j)
		{
			for (int i = 1; i < width - 1; ++i)
			{
				int cid = j * width + i;
				if (rawDepth[cid] <= 0.0f)
				{
					float* nf[8];
					nf[0] = rawDepth + cid - 1 - width;
					nf[1] = rawDepth + cid - width;
					nf[2] = rawDepth + cid + 1 - width;

					nf[3] = rawDepth + cid - 1;
					nf[4] = rawDepth + cid + 1;

					nf[5] = rawDepth + cid - 1 + width;
					nf[6] = rawDepth + cid + width;
					nf[7] = rawDepth + cid + 1 + width;

					float depth = 0.0f;
					float weight = 0.0f;
					for (int nid = 0; nid < 8; ++nid)
					{
						if (*nf[nid] > 0.0f)
						{
							depth += *nf[nid];
							weight += 1.0f;
						}
					}
					if (weight > 0.0f)	depth /= weight;

					tempRawDepth[cid] = depth;
				}
			}
		}
		std::copy(tempRawDepth.begin(), tempRawDepth.end(), rawDepth);
	}//*/


	for (int i = 0; i < width*height; ++i)
	{
		if (rawDepth[i] <= 0.0f)	depthData[i] = 0;
		else
		{
			float depthInMili = rawDepth[i] * 1000.0f;
			depthData[i] = (unsigned short)(depthInMili);
		}
	}
	cv::imwrite("depth.png", depthMap);

    if (std::string(filename).size() > 0)
	{
		std::ofstream outFile(filename, std::ios::binary);
		//outFile << width << " " << height << "\n";
		outFile.write((char*)rawDepth, sizeof(float)*width*height);
		outFile.close();
	}
}

double W2d(double r)
{
	if (r < 1.0)
		return 1.0 - 1.5*r*r*(1.0 - 0.5*r);
	else if (r < 2.0)
		return 0.25*(2.0 - r)*(2.0 - r)*(2.0 - r);
	else
		return 0.0;
}

void remedy_error_region_by_feature_tracking(int num_corr, double* p, double* q,
	int w, int h, double* tK, float* base_depth, float* input_depth,
	unsigned char* error_region_mask, float* output_depth,
	bool updateColor, unsigned char* base_color, unsigned char* input_color, unsigned char* colorMask,
	float alpha, float beta)
{
	//------------------------------------------------------------
	// ASSUME :
	//------------------------------------------------------------
	//  . {p} and {q} are points in 2D image.
	//  . {q} is not in the error region.
	//------------------------------------------------------------
	
	//------------------------------------------------------------
	// Find 3D points for reference and target frames
	//------------------------------------------------------------
	std::vector<Eigen::Vector3d> x(num_corr), y(num_corr);

	Eigen::Matrix3d K(tK);
	Eigen::Matrix3d invK = K.inverse();

	for (int i = 0; i < num_corr; ++i)
	{
		int ti = (int)floor(p[2 * i + 0]);
		int tj = (int)floor(p[2 * i + 1]);

		int pindex = tj * w + ti;

		double z = base_depth[pindex];

		x[i] = z * invK * Eigen::Vector3d(p[2 * i + 0], p[2 * i + 1], 1.0);

		ti = (int)floor(q[2 * i + 0]);
		tj = (int)floor(q[2 * i + 1]);

		pindex = tj * w + ti;

		z = input_depth[pindex];
		y[i] = z * invK * Eigen::Vector3d(q[2 * i + 0], q[2 * i + 1], 1.0);

		//if (i == 1)
		//{
		//	std::cout << i << " : " << x[i].transpose() << "..." << y[i].transpose() << std::endl;
		//}
	}

	// Convert float vector array for knn search
	std::vector<float>	tp(num_corr * 2), query(num_corr*2);
	for (int i = 0; i < 2*num_corr; ++i)
	{
		tp[i] = (float)p[i];
		query[i] = tp[i];
	}

	// Search neighbor points
	KnnSearch knnSearch(num_corr, tp.data(), 2);

	float averSpacing = knnSearch.find_aver_spacing();
	//std::cout << "- average spacing = " << averSpacing << std::endl;

	const int max_nn = 10;

	std::vector<int> indices(num_corr * max_nn);
	std::vector<float> dists(num_corr * max_nn);

	knnSearch.nearestSearch(num_corr, query.data(), max_nn, indices.data(), dists.data());
	//for (int i = 0; i < num_corr; ++i)
	//{
	//	std::cout << i << " (" << query[2*i+0] << "," << query[2 * i + 0] << ") : ";
	//	for (int j = 0; j < max_nn; ++j)
	//	{
	//		std::cout << indices[i*max_nn + j] << " ";
	//	}
	//	std::cout << std::endl;
	//}

	//----------------------------------------------------------------
	// Find local affine maps for all p's : y = Ax+b
	//----------------------------------------------------------------
	std::cout << "- Finding local affine maps with correspondences..." << std::flush;

	std::vector<Eigen::Matrix3d> A(num_corr);
	std::vector<Eigen::Vector3d> b(num_corr);
	std::vector<Eigen::Matrix4d> invAffine(num_corr);

	std::vector<int> nid(max_nn), nd(max_nn);

	int nn = 0;
	
	for (int i = 0; i < num_corr; ++i)
	{
		A[i].setIdentity();

		nn = 0;
		for (int j = 0; j < max_nn; ++j)
		{
			if (j > 0 && dists[i*max_nn + j] < 1.0e-9)	break;

			nid[nn] = indices[i*max_nn + j];
			nd[nn] = dists[i*max_nn + j];
			nn++;
		}

		nid.resize(nn);
		nd.resize(nn);

		/*if (nn < max_nn)
		{
			std::cout << i << " : ";
			for (int j = 0; j < max_nn; ++j)
			{
				std::cout << indices[i*max_nn + j] << " ";
			}
			std::cout << std::endl;
		}*/

		if (nn < 4) continue;

		std::vector<Eigen::Vector3d> np(nn), nq(nn);

		for (int j = 0; j < nn; ++j)
		{
			np[j] = x[nid[j]];
			nq[j] = y[nid[j]];
		}

		find_affine_transform(nn, (double*)np.data(), (double*)nq.data(), A[i].data(), b[i].data());

		Eigen::Matrix4d T;
		T.setIdentity();
		T.block<3, 3>(0, 0) = A[i];
		T.block<3, 1>(0, 3) = b[i];
		invAffine[i] = T.inverse();

		//std::cout << "3D Affine [" << i << "] :\n" << T << "\n"
		//	<< " det=" << T.block<3,3>(0,0).determinant() << std::endl;

		/*if (i == 1)
		{
			std::cout << i << "(" << nn << ") : ";
			for (int j = 0; j < nn; ++j)
			{
				std::cout << nid[j] << " ";
			}
			std::cout << std::endl;

			std::cout << "correspondences" << std::endl;
			for (int j = 0; j < nn; ++j)
			{
				std::cout << nid[j] << " :\n"
					<< x[nid[j]].transpose() << "(" << p[2 * nid[j] + 0] << "," << p[2 * nid[j] + 1] << ")\n"
					<< y[nid[j]].transpose() << "(" << q[2 * nid[j] + 0] << "," << q[2 * nid[j] + 1] << ")\n"
					<< (A[i]*x[nid[j]]+b[i]).transpose() << std::endl;
			}
		}//*/
	}
	//std::cout << "done" << std::endl;

	// Construct another 2D search structure on the target image
	std::vector<float> tq(num_corr * 2);
	for (int i = 0; i < tq.size(); ++i)
	{
		tq[i] = q[i];
	}

	KnnSearch knnSearch_target(num_corr, tq.data(), 2);
	knnSearch_target.nearestSearch(num_corr, (float*)tq.data(), max_nn, indices.data(), dists.data());

	//-----------------------------------------------------------------------
	// Find 2D affine maps from the target image to the source image
	// Note that this is necessary because we need find neighbor points 
	// for an error point on the target image.
	//-----------------------------------------------------------------------
	//std::cout << "- Finding 2D local affine maps with correspondences..." << std::flush;

	std::vector<Eigen::Matrix2d> A2(num_corr);
	std::vector<Eigen::Vector2d> b2(num_corr);
	//std::vector<Eigen::Matrix3d> TA(num_corr);

	nid.resize(max_nn), nd.resize(max_nn);

	nn = 0;

	for (int i = 0; i < num_corr; ++i)
	{
		A2[i].setIdentity();

		nn = 0;
		for (int j = 0; j < max_nn; ++j)
		{
			if (j > 0 && dists[i*max_nn + j] < 1.0e-9)	break;

			nid[nn] = indices[i*max_nn + j];
			nd[nn] = dists[i*max_nn + j];
			nn++;
		}

		nid.resize(nn);
		nd.resize(nn);

		//if (nn < max_nn)
		//{
		//	std::cout << i << " : ";
		//	for (int j = 0; j < max_nn; ++j)
		//	{
		//		std::cout << indices[i*max_nn + j] << " ";
		//	}
		//	std::cout << std::endl;
		//}

		if (nn < 3) continue;

		std::vector<Eigen::Vector2d> np(nn), nq(nn);

		for (int j = 0; j < nn; ++j)
		{
			nq[j] = Eigen::Vector2d(q[2 * nid[j] + 0], q[2 * nid[j] + 1]);
			np[j] = Eigen::Vector2d(p[2 * nid[j] + 0], p[2 * nid[j] + 1]);
		}

		find_2D_affine_transform(nn, (double*)nq.data(), (double*)np.data(), A2[i].data(), b2[i].data());

		//TA[i].setIdentity();
		//TA[i].block<2, 2>(0, 0) = A2[i];
		//TA[i].block<2, 1>(0, 2) = b2[i];

		Eigen::Matrix3d T;
		T.setIdentity();
		T.block<2, 2>(0, 0) = A2[i];
		T.block<2, 1>(0, 2) = b2[i];
		//std::cout << "2D Affine [" << i << "] :\n" << T << "\n"
		//	<< " det=" << T.block<2, 2>(0, 0).determinant() << std::endl;

		/*if (i == 88)
		{
			std::cout << i << "(" << nn << ") : ";
			for (int j = 0; j < nn; ++j)
			{
				std::cout << nid[j] << "(" << nd[j] << ") ";
			}
			std::cout << std::endl;

			std::cout << "correspondences" << std::endl;
			for (int j = 0; j < nn; ++j)
			{
				std::cout << nid[j] << " :\n"
					<< "(" << np[j].transpose() << ")\n"
					<< "(" << nq[j].transpose() << ")\n"
					<< (A2[i] * np[j] + b2[i]).transpose() << std::endl;
			}
		}//*/
	}
	std::cout << "done" << std::endl;

	//--------------------------------------------------------------------------------
	// Find the points in the error region as a new set of query points
	//--------------------------------------------------------------------------------
	std::vector<Eigen::Vector2f> errImgPt(w*h);
	std::vector<Eigen::Vector2i> errPixel(w*h);
	int nep = 0;
	for (int j = 0; j < h; ++j)
	{
		for (int i = 0; i < w; ++i)
		{
			int pid = j * w + i;
			if (error_region_mask[pid] >= 5)
			{
				errImgPt[nep] = Eigen::Vector2f(i, j);
				errPixel[nep] = Eigen::Vector2i(i, j);

				//if (i == 935 && j == 1229)
				//	std::cout << "errPt[" << nep << "]=" << errPixel[nep].transpose() << std::endl;

				nep++;
			}
		}// end for(i)
	}// end for(j)

	errImgPt.resize(nep);
	errPixel.resize(nep);
	//std::cout << "- # of error points = " << nep << std::endl;

	// Find neighbors for error points
	//std::cout << "- Finding neighbors of error points..." << std::flush;
	indices.resize(nep*max_nn);
	dists.resize(nep*max_nn);
	knnSearch_target.nearestSearch(nep, (float*)errImgPt.data(), max_nn, indices.data(), dists.data());
	//std::cout << "done" << std::endl;

	//-------------------------------------------------------------------------------------
	// Find points corresponding to error points on the source image and in 3D space
	//-------------------------------------------------------------------------------------
	std::vector<Eigen::Vector2f> corrPt(nep);
	std::vector<Eigen::Vector3d> corrPt3d(nep);

	nid.resize(max_nn);
	nd.resize(max_nn);

	for (int i = 0; i < nep; ++i)
	{
		int nn = 0;
		for (int j = 0; j < max_nn; ++j)
		{
			if (j > 0 && dists[i*max_nn + j] < 1.0e-9)	break;

			nid[nn] = indices[i*max_nn + j];
			nd[nn] = dists[i*max_nn + j];
			nn++;
		}
		nid.resize(nn);
		nd.resize(nn);

		double supp = 0.6 * sqrt(nd[nn - 1]);

		Eigen::Vector2d errPtOnSrc(0.0, 0.0);

		double sumW = 0.0;
		for (int j = 0; j < nn; ++j)
		{
			double w = W2d(sqrt(nd[j]) / supp);

			Eigen::Vector2d  tempPt = A2[nid[j]] * Eigen::Vector2d(errImgPt[i](0), errImgPt[i](1)) + b2[nid[j]];
			errPtOnSrc += w * tempPt;

			//if (i == 37727)
			//{
			//	std::cout << "w=" << w << ", " << errImgPt[i].transpose() << "->" << tempPt.transpose() << std::endl;
			//}

			sumW += w;
		}

		if(sumW > 0.0)	errPtOnSrc /= sumW;

		int ti = (int)floor(errPtOnSrc(0));
		int tj = (int)floor(errPtOnSrc(1));

		int pindex = tj * w + ti;

		double z = base_depth[pindex];

		corrPt[i] = Eigen::Vector2f(errPtOnSrc(0), errPtOnSrc(1));
		corrPt3d[i] = z * invK * Eigen::Vector3d(errPtOnSrc(0), errPtOnSrc(1), 1.0);

		//if (i == 37727)
		//{
		//	std::cout << "corrPtOnSrc[" << i << "] = " << corrPt[i].transpose() << "\n"
		//		<< "corrPt3dOnSrc[" << i << "] = " << corrPt3d[i].transpose() << std::endl;
		//}
	}

	//--------------------------------------------------------------------------
	// Map corresponding 3D points to the targe space 
	// and project to the target image to update the depth map 
	//--------------------------------------------------------------------------
	std::cout << "- Mapping 3D points in the reference space to the target space..." << std::flush;

	std::vector<Eigen::Vector3d> corrected_err_pt(nep);

	// Search again for corresponding points on the source image
	knnSearch.nearestSearch(nep, (float*)corrPt.data(), max_nn, indices.data(), dists.data());
	
	memcpy(output_depth, input_depth, sizeof(float)*w*h);

	for (int i = 0; i < nep; ++i)
	{
		int nn = 0;
		for (int j = 0; j < max_nn; ++j)
		{
			if (j > 0 && dists[i*max_nn + j] < 1.0e-9)	break;

			nid[nn] = indices[i*max_nn + j];
			nd[nn] = dists[i*max_nn + j];
			nn++;
		}

		Eigen::Vector3d pt3dInTgt(0.0, 0.0, 0.0);

		double supp = 0.6 * sqrt(nd[nn - 1]);
		double sumW = 0.0;

		for (int j = 0; j < nn; ++j)
		{
			double w = W2d(sqrt(nd[j]) / supp);

			Eigen::Vector3d tempPt = A[nid[j]] * corrPt3d[i] + b[nid[j]];
			pt3dInTgt += w*tempPt;

			//if (i == 37727)
			//{
			//	std::cout << "w=" << w << ", " << corrPt3d[i].transpose() << "->" << tempPt.transpose() << std::endl;
			//}

			sumW += w;
		}

		pt3dInTgt /= sumW;

		corrected_err_pt[i] = pt3dInTgt;

		float corrected_z = (K * corrected_err_pt[i])(2);

		Eigen::Vector3d prPt = K * pt3dInTgt;
		prPt /= prPt(2);
		/*if (i == 37727)
		{
			std::cout << "Mapped 3D point in target space : " << pt3dInTgt.transpose() << std::endl;
			std::cout << "projected point in target image = " << prPt.transpose() << std::endl;
			std::cout << "(" << errPixel[i].transpose() << ")" << std::endl;
		}//*/

		output_depth[errPixel[i](1)*w +errPixel[i](0)] = corrected_z;
	}
	std::cout << std::endl;

	if (updateColor)
	{
		for (int i = 0; i < nep; ++i)
		{
			int pi = (int)corrPt[i](0);
			int pj = (int)corrPt[i](1);

			int bgr[3];
			bgr[0] = alpha * base_color[3 * (pj*w + pi) + 0] + beta;
			bgr[1] = alpha * base_color[3 * (pj*w + pi) + 1] + beta;
			bgr[2] = alpha * base_color[3 * (pj*w + pi) + 2] + beta;

			if (bgr[0] > 255)	bgr[0] = 255;
			if (bgr[1] > 255)	bgr[1] = 255;
			if (bgr[2] > 255)	bgr[2] = 255;

			if (bgr[0] < 0)		bgr[0] = 0;
			if (bgr[1] < 0)		bgr[1] = 0;
			if (bgr[2] < 0)		bgr[2] = 0;

			input_color[3 * (errPixel[i](1)*w + errPixel[i](0)) + 0] = bgr[0];
			input_color[3 * (errPixel[i](1)*w + errPixel[i](0)) + 1] = bgr[1];
			input_color[3 * (errPixel[i](1)*w + errPixel[i](0)) + 2] = bgr[2];

			colorMask[3 * (errPixel[i](1)*w + errPixel[i](0)) + 0] = bgr[0];
			colorMask[3 * (errPixel[i](1)*w + errPixel[i](0)) + 1] = bgr[1];
			colorMask[3 * (errPixel[i](1)*w + errPixel[i](0)) + 2] = bgr[2];
		}
	}
}

void write_point_cloud_in_ply(const char* filename, int n, double* p)
{
    std::ofstream outFile(filename);

    outFile << "ply" << "\n"
        << "format ascii 1.0" << "\n"
        << "comment" << "\n"
        << "comment" << "\n"
        << "element vertex " << n << "\n"
        << "property float x" << "\n"
        << "property float y" << "\n"
        << "property float z" << "\n"
        << "end_header" << "\n";

    for (int i = 0; i < n; ++i)
    {
        Eigen::Vector3d tp(p + 3 * i);
        outFile << tp.transpose() << "\n";
    }

    outFile.close();
}

void write_point_cloud_in_ply(const char* filename, int n, double* p, unsigned char* c)
{
    std::ofstream outFile(filename);

    outFile << "ply" << "\n"
        << "format ascii 1.0" << "\n"
        << "comment" << "\n"
        << "comment" << "\n"
        << "element vertex " << n << "\n"
        << "property float x" << "\n"
        << "property float y" << "\n"
        << "property float z" << "\n"
        << "property uchar red " << "\n"
        << "property uchar green" << "\n"
        << "property uchar blue" << "\n"
        << "end_header" << "\n";

    for (int i = 0; i < n; ++i)
    {
        Eigen::Vector3d tp(p + 3 * i);
        outFile << tp.transpose() << " ";

        // BGR
        outFile 
            << (int)c[3 * i + 2] << " " 
            << (int)c[3 * i + 1] << " " 
            << (int)c[3 * i + 0] << "\n";
    }

    outFile.close();
}

void write_point_cloud_in_pcd(const char* filename, int n, double* p)
{
	std::ofstream outFile(filename, std::ios::binary);

	if (!outFile.is_open())
	{
		std::cerr << "Error! Can't open " << filename << " for writing" << std::endl;
		return;
	}

	std::ostringstream ss;

	// header with number of points
	ss << "# .PCD v.7 - Point Cloud Data file format\n"
		<< "VERSION .7\n"
		<< "FIELDS x y z rgb\n"
		<< "SIZE 4 4 4 4\n"
		<< "TYPE F F F F\n"
		<< "COUNT 1 1 1 1\n"
		<< "WIDTH " << n << "\n"
		<< "HEIGHT 1\n"
		<< "VIEWPOINT 0 0 0 1 0 0 0\n"
		<< "POINTS " << n << "\n"
		//<< "DATA ascii\n";
		<< "DATA binary\n";

	outFile.write((char*)ss.str().c_str(), ss.str().size());

	int r, g, b;
	//unsigned int rgb;

	float tp[3];
	for (int i = 0; i < n; ++i)
	{
		// position
		tp[0] = (float)p[3 * i + 0];
		tp[1] = (float)p[3 * i + 1];
		tp[2] = (float)p[3 * i + 2];
		outFile.write((char*)tp, sizeof(float) * 3);
	}

	outFile.close();
}
void write_point_cloud_in_pcd(const char* filename, int n, double* p, unsigned char* c)
{
    std::ofstream outFile(filename, std::ios::binary);

    if (!outFile.is_open())
    {
        std::cerr << "Error! Can't open " << filename << " for writing" << std::endl;
        return;
    }

    std::ostringstream ss;

    // header with number of points
    ss << "# .PCD v.7 - Point Cloud Data file format\n"
        << "VERSION .7\n"
        << "FIELDS x y z rgb\n"
        << "SIZE 4 4 4 4\n"
        << "TYPE F F F F\n"
        << "COUNT 1 1 1 1\n"
        << "WIDTH " << n << "\n"
        << "HEIGHT 1\n"
        << "VIEWPOINT 0 0 0 1 0 0 0\n"
        << "POINTS " << n << "\n"
        //<< "DATA ascii\n";
        << "DATA binary\n";

    outFile.write((char*)ss.str().c_str(), ss.str().size());

    int r, g, b;
    //unsigned int rgb;

    float tp[3];
    for (int i = 0; i < n; ++i)
    {
        // position
        tp[0] = (float)p[3 * i + 0];
        tp[1] = (float)p[3 * i + 1];
        tp[2] = (float)p[3 * i + 2];
        outFile.write((char*)tp, sizeof(float) * 3);

        // color(BGR)
        b = c[3 * i + 0];
        g = c[3 * i + 1];
        r = c[3 * i + 2];

        unsigned int rgb = ((r << 16) | (g << 8) | (b << 0));

        outFile.write((char*)&rgb, sizeof(unsigned int));
    }

    outFile.close();
}

void write_point_cloud_in_pcd(const char* filename, int n, double* p, unsigned char* c, double* normal)
{
	std::ofstream outFile(filename, std::ios::binary);

	if (!outFile.is_open())
	{
		std::cerr << "Error! Can't open " << filename << " for writing" << std::endl;
		return;
	}

	std::ostringstream ss;

	// header with number of points
	ss << "# .PCD v.7 - Point Cloud Data file format\n"
		<< "VERSION .7\n"
		<< "FIELDS x y z rgb normal_x normal_y normal_z\n"
		<< "SIZE 4 4 4 4 4 4 4\n"
		<< "TYPE F F F F F F F\n"
		<< "COUNT 1 1 1 1 1 1 1\n"
		<< "WIDTH " << n << "\n"
		<< "HEIGHT 1\n"
		<< "VIEWPOINT 0 0 0 1 0 0 0\n"
		<< "POINTS " << n << "\n"
		//<< "DATA ascii\n";
		<< "DATA binary\n";

	outFile.write((char*)ss.str().c_str(), ss.str().size());

	int r, g, b;
	//unsigned int rgb;

	float tp[3];
	for (int i = 0; i < n; ++i)
	{
		// position
		tp[0] = (float)p[3 * i + 0];
		tp[1] = (float)p[3 * i + 1];
		tp[2] = (float)p[3 * i + 2];
		outFile.write((char*)tp, sizeof(float) * 3);

		// color(BGR)
		b = c[3 * i + 0];
		g = c[3 * i + 1];
		r = c[3 * i + 2];

		unsigned int rgb = ((r << 16) | (g << 8) | (b << 0));

		outFile.write((char*)&rgb, sizeof(unsigned int));

		tp[0] = (float)normal[3 * i + 0];
		tp[1] = (float)normal[3 * i + 1];
		tp[2] = (float)normal[3 * i + 2];
		outFile.write((char*)tp, sizeof(float) * 3);
	}

	outFile.close();
}