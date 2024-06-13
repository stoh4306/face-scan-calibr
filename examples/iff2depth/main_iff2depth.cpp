
#include<opencv2/opencv.hpp>
#include<iostream>
#include<string>
#include"renew_iff.h"


void decompose_image_file_name(std::string imgFileName, std::string& prefix, int& paddings, int& startFrame, std::string& ext)
{
	size_t found = imgFileName.find_last_of('.');
	ext = imgFileName.substr(found + 1);

	int count = 0;
	for (size_t i = found - 1; i >= 0; --i)
	{
		int ascii = (int)imgFileName[i];

		if (ascii >= 48 && ascii <= 57)
			count++;
		else
			break;
	}

	prefix = imgFileName.substr(0, imgFileName.size() - (count + 1 + ext.size()));
	paddings = count;

	std::string numStr = imgFileName.substr(prefix.size(), count);

	startFrame = atoi(numStr.c_str());

	std::cout << prefix << " " << count << "(" << startFrame << ") " << ext << std::endl;

}


int main(int argc, char* argv[]) {
	//std::string fileName = "iff8_0019.iff";
	if (argc < 2) {
		std::cout << "iffOpen.exe startFile numOfFrame" << std::endl;
		return -1;
	}
	std::string prefix, ext;
	int paddings, startFrame;
	int numOfFrame = 0;
	numOfFrame = atoi(argv[2]);
	decompose_image_file_name(argv[1], prefix, paddings, startFrame, ext);
	if (ext != "iff")
		return -2;
	for (int i = 0; i < numOfFrame; ++i) {
		std::stringstream sstream;
		sstream << prefix << std::setw(paddings) << std::setfill('0') << i + startFrame << "." << ext;

		std::cout << sstream.str() << " load..." << std::endl;
		loadIff(sstream.str());
	}


	//std::string fileName = loadIff(fileName);
	//cv::Mat mat;
	//mat = cv::imread("1.exr", cv::IMREAD_ANYCOLOR | cv::IMREAD_ANYDEPTH);
	//cv::namedWindow("1");
	//cv::imshow("1",mat);
	//cv::waitKey();
	return 0;
}

