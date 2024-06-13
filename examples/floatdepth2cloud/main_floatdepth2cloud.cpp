#include <iostream>
#include <opencv2/opencv.hpp>
#include <Eigen/Dense>

#include "calibr.h"
#include "cloud2.h"

int main(int argc, char** argv)
{
	std::cout << "**** Convert 32-bit float depth data stored in 4-channel image to cloud and disparity map" << std::endl;
	if(argc < 5)
	{
		std::cout << "- Usage : " << argv[0] << " <depth_img> <K_file> <Right_CamPose_file> <texture_img>" << std::endl;
		return -1;
	}
	const char* imgFile = argv[1];
	const char* K_file = argv[2];
	const char* RCP_file = argv[3]; // NOTE : Camera pose for the right camera
	const char* tex_file = argv[4]; 

	std::cout << "- Reading depth image file : " << imgFile << std::flush;
	cv::Mat depthImg = cv::imread(argv[1], -1);
	std::cout << std::endl;

	int width = depthImg.cols;
	int height = depthImg.rows;

	float* depthData = (float*)depthImg.data;

	auto result = std::minmax_element(depthData, depthData + width * height);
	//std::cout << *result.first << ", " << *result.second << std::endl;

	Eigen::Matrix3f K;
	Eigen::Matrix4f RCP;

	std::cout << "- Reading calibration : " << K_file << ", " << RCP_file << std::flush;
	read_calibrations(K_file, RCP_file, K.data(), RCP.data());
	std::cout << std::endl;

	Eigen::Matrix4f LCP = Eigen::Matrix4f::Identity();
	
	Eigen::Matrix3f invK = K.inverse();
	//std::cout << "- invK = \n" << invK << std::endl;

	PointCloud2 cloud;
	cloud.setSize(width*height);

	cv::Mat dispImg(height, width, CV_8UC4);
	float* dispData = (float*)dispImg.data;
	std::fill(dispData, dispData + width * height, 0.0f);
	Eigen::Matrix4f world2local = RCP.inverse();

	//Texture file
	cv::Mat texImg = cv::imread(tex_file, -1);
	int tex_width = texImg.cols;
	int tex_height = texImg.rows; //std::cout << tex_width << " " << tex_height << std::endl;
	int ch = texImg.channels(); //std::cout << "channels=" << ch << std::endl;
	unsigned char* texData = texImg.data;
	
	unsigned int count = 0;

	float INF = std::numeric_limits<float>::infinity();

	for (int j = 0; j < height; ++j)
	{
		for (int i = 0; i < width; ++i)
		{
			const float z = depthData[j*width + i];
			if (z > 0.0f && z < INF)
			{
				Eigen::Vector3f q(i, j, 1.0f);
				Eigen::Vector3f p = z * (invK*q);

				Eigen::Vector4f x = Eigen::Vector4f(p(0), -p(1), -p(2), 1.0f);
				cloud.pos_[count] = mg::Vector3f(x.data());
				//cloud.col_[count] = mg::Vector3f(255.0f, 255.0f, 255.0f);
				cloud.col_[count] = mg::Vector3f(texData[ch*(j*tex_width + i) + 0], 
												 texData[ch*(j*tex_width + i) + 1],
												 texData[ch*(j*tex_width + i) + 2]);
				count++;

				// disparity
				Eigen::Vector4f y = world2local * Eigen::Vector4f(p(0), p(1), p(2), 1.0f);
				Eigen::Vector3f tq = K*Eigen::Vector3f(y.data());
				
				float u = tq(0) / tq(2);
				
				//float v = tq(1) / tq(2);
				//if (fabs(v - j) > 0.001)
				//	std::cout << "(u,v)=" << u << ", " << v << std::endl;

				if ( u >= 0.0 && u < width)
				{
					dispData[j*width + i] = i - u;
				}
			}
		}
	}

	cloud.pos_.resize(count);

	bool colorFlag = false;
	if (!texImg.empty()) colorFlag = true;
	bool normalFlag = false;
	const char* cloudFile = "cloud.ply";
	std::cout << "- Exporting cloud and disparity : " << cloudFile << ", " << std::flush;
	cloud.exportToPly(cloudFile, colorFlag, normalFlag);
	
	const char* dispFile = "disparity.png";
	std::cout << dispFile << "..." << std::flush;
	cv::imwrite("disparity.png", dispImg);
	std::cout << "done" << std::endl;

	//cv::Mat dispImg2(height, width, CV_8UC1);
	//unsigned char* disp8data = dispImg2.data;
	//for (int i = 0; i < width*height; ++i)
	//{
	//	disp8data[i] = (unsigned int)(dispData[i] / width*255.9999);
	//}
	//cv::imshow("disp8", dispImg2);
	//cv::waitKey(0);

	return 0;
}