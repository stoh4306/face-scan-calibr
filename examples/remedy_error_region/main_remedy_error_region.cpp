#include "calibr.h"
#include "reconst.h"

#include <vector>
#include <iostream>
#include <fstream>
#include <string>

#include <opencv2/opencv.hpp>
#include <Eigen/Dense>

bool readTrackedPoints(const char* filename, std::vector<Eigen::Vector2d>& p)
{
	std::cout << " . Reading tracking points : " << filename << "..." << std::flush;

	std::ifstream inFile(filename);
	if (!inFile.is_open())
	{
		std::cout << "Error, can't open " << filename << " for reading!" << std::endl;
		return false;
	}

	int numPts = 0;

	inFile >> numPts;

	std::string temp;
	std::getline(inFile, temp, '\n');

	p.resize(numPts);

	for (int i = 0; i < numPts; ++i)
	{
		inFile >> p[i](0) >> p[i](1);
		std::getline(inFile, temp, '\n');
	}

	inFile.close();

	//std::cout << p[0].transpose() << "-" << p[numPts - 1].transpose() << std::endl;

	std::cout << "done" << std::endl;

	return true;
}

bool readDepth32File(const char* filename, std::vector<float>& depth)
{
	cv::Mat img = cv::imread(filename, -1);

	if (img.empty())
	{
		std::cout << "Error, can't open image : " << filename << std::endl;
		return false;
	}

	if (img.channels() != 4)
	{
		std::cout << "Error, improper depth file. The file should be in a 32-bit image format." << std::endl;
		return false;
	}

	int w = img.cols;
	int h = img.rows;

	depth.resize(w*h);
	memcpy(depth.data(), img.data, sizeof(float)*w*h);

	return true;
}

bool readDepth16File(const char* filename, std::vector<float>& depth)
{
	std::cout << " . Reading depth image : " << filename << "..." << std::flush;

	cv::Mat depthImg = cv::imread(filename, -1);

	if (depthImg.empty())
	{
		std::cout << "ERROR" << std::endl;
		return false;
	}

	int w = depthImg.cols;
	int h = depthImg.rows;

	depth.resize(w*h);

	unsigned short* depthData = (unsigned short*)depthImg.data;

	for (int i = 0; i < w*h; ++i)
	{
		depth[i] = depthData[i] * 0.001f;
	}

	//std::cout << filename << " : " << depth[0] << "-" << depth[w*h - 1] << std::endl;
	//auto result = std::minmax_element(depth.begin(), depth.end());
	//std::cout << filename << " : min = " << *result.first << ", max=" << *result.second << std::endl;
	std::cout << "done" << std::endl;

	return true;
}

void remove_points_inside_mask(std::vector<Eigen::Vector2d>& srcPt,
	std::vector<Eigen::Vector2d>& tgtPt,
	cv::Mat& maskImg)
{
	std::cout << "- Removing points inside mask..." << std::flush;

	int n = srcPt.size();
	std::vector<Eigen::Vector2d> tempSrcPt(n), tempTgtPt(n);

	int width = maskImg.cols;

	unsigned char* maskData = maskImg.data;

	int id = 0;

	std::vector<int> removedPtId(n);
	int idrp = 0;
	for (int i = 0; i < tgtPt.size(); ++i)
	{
		int u = (int)tgtPt[i](0);
		int v = (int)tgtPt[i](1);

		int pid = v * width + u;

		if (maskData[pid] < 50)
		{
			tempSrcPt[id] = srcPt[i];
			tempTgtPt[id] = tgtPt[i];

			id++;
		}
		else
		{
			removedPtId[idrp] = i;
			idrp++;
		}
	}

	tempSrcPt.resize(id);
	tempTgtPt.resize(id);
	removedPtId.resize(idrp);
	std::cout << n << "->" << id << std::flush;
	std::cout << "(removed pts : ";
	for (int i = 0; i < idrp; ++i)
		std::cout << removedPtId[i] << " ";
	std::cout << ")" << std::endl;

	std::swap(srcPt, tempSrcPt);
	std::swap(tgtPt, tempTgtPt);
}

bool readIntrinsicMatrix(const char* filename, Eigen::Matrix3d& K)
{
	std::cout << " . Reading K file : " << filename << "..." << std::flush;

	std::ifstream inFile(filename);
	if (!inFile.is_open())
	{
		std::cout << "ERROR" << std::endl;
		return false;
	}
	std::cout << std::endl;

	inFile >> K(0, 0) >> K(0, 1) >> K(0, 2)
		>> K(1, 0) >> K(1, 1) >> K(1, 2)
		>> K(2, 0) >> K(2, 1) >> K(2, 2);
	inFile.close();

	//std::cout << "K=\n" << K << std::endl;

	return true;
}

int main(int argc, char** argv)
{
	std::cout << "=================================================" << "\n"
		<< " Remedy error region" << "\n"
		<< "=================================================" << std::endl;

	if (argc < 9)
	{
		std::cout << "- Usage : remedy_error_region.exe "
			<< "<src_pts> <tgt_pts> <src_depth> <tgt_depth> <mask> <K_file> <update_color> <src_color> <tgt_color> <alpha> <beta>"
			<< std::endl;
		return 0;
	}

	const char* srcPtFile = argv[1];
	const char* tgtPtFile = argv[2];
	const char* srcDepthFile = argv[3];
	const char* tgtDepthFile = argv[4];
	const char* maskFile = argv[5];
	const char* calib_K_file = argv[6];
	const bool  updateColor = atoi(argv[7]) > 0 ? true : false;
	const char* srcColorFile = argv[8];
	const char* tgtColorFile = argv[9];
	const float alpha = atof(argv[10]);
	const float beta = atof(argv[11]);

	// Read tracked points
	std::vector<Eigen::Vector2d> srcPt, tgtPt;
	readTrackedPoints(srcPtFile, srcPt);
	readTrackedPoints(tgtPtFile, tgtPt);

	// Read depth files
	std::vector<float> srcDepth, tgtDepth;
	readDepth16File(srcDepthFile, srcDepth);
	readDepth16File(tgtDepthFile, tgtDepth);

	// Read mask file
	std::cout << " . Reading mask image : " << maskFile << "..." << std::flush;

	cv::Mat maskImg = cv::imread(maskFile, -1);
	if (maskImg.empty())
	{
		std::cout << "ERROR" << std::endl;
		std::cout << "*** Program has stopped abnormally" << std::endl;
		return 0;
	}
	else
		std::cout << "done" << std::endl;

	int width = maskImg.cols;
	int height = maskImg.rows;

	// Read K file
	Eigen::Matrix3d K;
	readIntrinsicMatrix(calib_K_file, K);

	// Read RGB images
	cv::Mat srcColor, tgtColor, colorMask;
	if (updateColor)
	{
		std::cout << " . Reading RGB images : " << srcColorFile << ", " << tgtColorFile << "..." << std::flush;

		srcColor = cv::imread(srcColorFile);
		tgtColor = cv::imread(tgtColorFile);

		if (srcColor.empty() || tgtColor.empty())
		{
			std::cout << "ERROR" << std::endl;
			std::cout << "*** Program has stopped abnormally" << std::endl;
			return 0;
		}

		colorMask = cv::Mat::zeros(srcColor.size(), srcColor.type());
	}

	// Refine tracked points by removing points inside the mask
	remove_points_inside_mask(srcPt, tgtPt, maskImg);

	// 
	std::vector<float> outDepth(width*height);
	remedy_error_region_by_feature_tracking(srcPt.size(), (double*)srcPt.data(), (double*)tgtPt.data(),
		width, height, K.data(), srcDepth.data(), tgtDepth.data(), maskImg.data, outDepth.data(),
		updateColor, srcColor.data, tgtColor.data, colorMask.data, alpha, beta);

	cv::Mat outDepthImg(height, width, CV_8UC4);
	memcpy(outDepthImg.data, outDepth.data(), sizeof(float)*width*height);

	const char* outFile = "outDepth.png";
	cv::imwrite(outFile, outDepthImg);

	std::cout << "- Output depth file exported : " << outFile << std::endl;

	if (updateColor)
	{
		const char* outColorFile = "outColor.png";
		const char* outColorMaskFile = "outColorMask.png";
		cv::imwrite(outColorFile, tgtColor);
		cv::imwrite(outColorMaskFile, colorMask);
		std::cout << "- Output color file exported : " << outColorFile << ", " << outColorMaskFile << std::endl;
	}

	return 0;
}