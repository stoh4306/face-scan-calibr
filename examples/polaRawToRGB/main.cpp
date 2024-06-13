//#include<Windows.h>
#include<iostream>
#include<vector>
#include<opencv2/opencv.hpp>
#include<string>
#include<sstream>




void decompose_image_file_name(std::string imgFileName, std::string& prefix, int& paddings, int& startFrame, std::string& ext);

std::vector<std::wstring> getFileInFolderW(std::wstring folder, std::wstring ext) {
	std::vector<std::wstring> names;
	wchar_t search_path[255];
	::wsprintfW(search_path, L"%s/%s", folder.c_str(), ext.c_str());
	WIN32_FIND_DATAW fd;
	HANDLE hFind = ::FindFirstFileW(search_path, &fd);

	if (hFind != INVALID_HANDLE_VALUE) {
		do {
			if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
				names.push_back(fd.cFileName);
			}
		} while (::FindNextFileW(hFind, &fd));
		::FindClose(hFind);
	}
	return names;
}

std::vector<std::string> getFileInFolderA(std::string folder, std::string ext) {
	std::vector<std::string> names;
	char search_path[255];
	::wsprintf(search_path, "%s\\*.%s", folder.c_str(), ext.c_str());
	WIN32_FIND_DATA fd;
	HANDLE hFind = ::FindFirstFile(search_path, &fd);

	if (hFind != INVALID_HANDLE_VALUE) {
		do {
			if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
				names.push_back(fd.cFileName);
			}
		} while (::FindNextFile(hFind, &fd));
		::FindClose(hFind);
	}
	return names;
}


int main(int argc, char* argv[]) {
	if (argc < 4) {
		std::cout << "plImgeLoad.exe folder ext outPrefix( + num_) paddings" << std::endl;
		getchar();
		return -1;
	}
		
	int paddings = atoi(argv[4]);
	std::string outPreFix,ext,outFolder;
	std::string inFolder = argv[1];
	//int paddings;
	//int startIndex;
	////decompose_image_file_name(argv[1], outPreFix, paddings, startIndex, ext);
	//int end = startIndex + atoi(argv[2]);
	std::vector<std::string> files = getFileInFolderA(argv[1], argv[2]);
	//std::wstring tmpWstring = argv[3];
	outPreFix = argv[3];
	//outPreFix.assign(tmpWstring.begin(), tmpWstring.end());

	outFolder = argv[1];
	//tmpWstring = argv[1];
	//outFolder.assign(tmpWstring.begin(), tmpWstring.end());

	ext = argv[2];

	for (int i = 0; i < files.size(); ++i) {
		std::string fileName;
		fileName = inFolder + "\\" + files[i];
		cv::Mat mat;
		mat = cv::imread(fileName, cv::IMREAD_UNCHANGED);
		int w = mat.cols;
		int h = mat.rows;
		cv::Mat mat1, mat2, mat3, mat4;
		mat1 = cv::Mat(h / 2, w / 2, mat.type());
		mat2 = cv::Mat(h / 2, w / 2, mat.type());
		mat3 = cv::Mat(h / 2, w / 2, mat.type());
		mat4 = cv::Mat(h / 2, w / 2, mat.type());
		for (int j = 0; j < h - 1; j += 2) {
			for (int i = 0; i < w - 1; i += 2) {
				mat1.at<cv::uint8_t>(j / 2, i / 2) = mat.at<cv::uint8_t>(j, i);
				mat2.at<cv::uint8_t>(j / 2, i / 2) = mat.at<cv::uint8_t>(j, i + 1);
				mat3.at<cv::uint8_t>(j / 2, i / 2) = mat.at<cv::uint8_t>(j + 1, i);
				mat4.at<cv::uint8_t>(j / 2, i / 2) = mat.at<cv::uint8_t>(j + 1, i + 1);
			}
		}
		cv::Mat mat11, mat22, mat33, mat44;
		cv::cvtColor(mat1, mat11, cv::COLOR_BayerRG2RGB);
		cv::cvtColor(mat2, mat22, cv::COLOR_BayerRG2RGB);
		cv::cvtColor(mat3, mat33, cv::COLOR_BayerRG2RGB);
		cv::cvtColor(mat4, mat44, cv::COLOR_BayerRG2RGB);

		//std::string fileName;
		std::string name1 = outFolder + "\\" + outPreFix + "1_";
		std::string name2 = outFolder + "\\" + outPreFix + "2_";
		std::string name3 = outFolder + "\\" + outPreFix + "3_";
		std::string name4 = outFolder + "\\" + outPreFix + "4_";
		//cv::imwrite("mat1.png", mat1);
		//cv::imwrite("mat2.png", mat2);
		//cv::imwrite("mat3.png", mat3);
		//cv::imwrite("mat4.png", mat4);
		//fileName = argv[1];
		//fileName = fileName.replace(fileName.size() - 4, fileName.size() - 4, name1);
		//int lastFolder = fileName.find_last_of("\\") + 1;
		//if (lastFolder < 0) {
		//	lastFolder = fileName.find_last_of("/") + 1;
		//}
		std::ostringstream outFileStream;
		outFileStream << name1 << std::setw(paddings) << std::setfill('0') << i << "." << ext;
		cv::imwrite(outFileStream.str(), mat11);
		std::cout << outFileStream.str() << " Save." << std::endl;
		outFileStream.str("");

		outFileStream << name2 << std::setw(paddings) << std::setfill('0') << i << "." << ext;
		cv::imwrite(outFileStream.str(), mat22);
		std::cout << outFileStream.str() << " Save." << std::endl;
		outFileStream.str("");
		outFileStream << name3 << std::setw(paddings) << std::setfill('0') << i << "." << ext;
		cv::imwrite(outFileStream.str(), mat33);
		std::cout << outFileStream.str() << " Save." << std::endl;
		outFileStream.str("");
		outFileStream << name4 << std::setw(paddings) << std::setfill('0') << i << "." << ext;
		cv::imwrite(outFileStream.str(), mat44);
		std::cout << outFileStream.str() << " Save." << std::endl;
		outFileStream.str("");
		std::cout << "==========================" << std::endl;
	}


	//cv::Mat mat;
	//mat = cv::imread(argv[1], cv::IMREAD_UNCHANGED);
	//int w = mat.cols;
	//int h = mat.rows;
	//cv::Mat mat1, mat2, mat3, mat4;
	//mat1 = cv::Mat(h / 2, w / 2, mat.type());
	//mat2 = cv::Mat(h / 2, w / 2, mat.type());
	//mat3 = cv::Mat(h / 2, w / 2, mat.type());
	//mat4 = cv::Mat(h / 2, w / 2, mat.type());
	//for (int j = 0; j < h - 1; j += 2) {
	//	for (int i = 0; i < w - 1; i += 2) {
	//		mat1.at<cv::uint8_t>(j / 2, i / 2) = mat.at<cv::uint8_t>(j, i);
	//		mat2.at<cv::uint8_t>(j / 2, i / 2) = mat.at<cv::uint8_t>(j, i + 1);
	//		mat3.at<cv::uint8_t>(j / 2, i / 2) = mat.at<cv::uint8_t>(j + 1, i);
	//		mat4.at<cv::uint8_t>(j / 2, i / 2) = mat.at<cv::uint8_t>(j + 1, i + 1);
	//	}
	//}
	////for (int j = 0; j < h - 1; j += 4) {
	////	for (int i = 0; i < w - 1; i += 4) {
	////		mat1.at<cv::uint8_t>(j / 2, i / 2) = mat.at<cv::uint8_t>(j, i);
	////		mat1.at<cv::uint8_t>(j / 2, i / 2 + 1) = mat.at<cv::uint8_t>(j, i + 1);
	////		mat1.at<cv::uint8_t>(j / 2 + 1, i / 2) = mat.at<cv::uint8_t>(j + 1, i);
	////		mat1.at<cv::uint8_t>(j / 2 + 1, i / 2 + 1) = mat.at<cv::uint8_t>(j + 1, i + 1);
	////					
	////		mat2.at<cv::uint8_t>(j / 2, i / 2) = mat.at<cv::uint8_t>(j, i + 2);
	////		mat2.at<cv::uint8_t>(j / 2, i / 2 + 1) = mat.at<cv::uint8_t>(j, i + 3);
	////		mat2.at<cv::uint8_t>(j / 2 + 1, i / 2) = mat.at<cv::uint8_t>(j + 1, i + 2);
	////		mat2.at<cv::uint8_t>(j / 2 + 1, i / 2 + 1) = mat.at<cv::uint8_t>(j + 1, i + 3);
	////					
	////		mat3.at<cv::uint8_t>(j / 2, i / 2) = mat.at<cv::uint8_t>(j + 2, i);
	////		mat3.at<cv::uint8_t>(j / 2, i / 2 + 1) = mat.at<cv::uint8_t>(j + 2, i + 1);
	////		mat3.at<cv::uint8_t>(j / 2 + 1, i / 2) = mat.at<cv::uint8_t>(j + 3, i);
	////		mat3.at<cv::uint8_t>(j / 2 + 1, i / 2 + 1) = mat.at<cv::uint8_t>(j + 3, i + 1);
	////					
	////		mat4.at<cv::uint8_t>(j / 2, i / 2) = mat.at<cv::uint8_t>(j + 2, i + 2);
	////		mat4.at<cv::uint8_t>(j / 2, i / 2 + 1) = mat.at<cv::uint8_t>(j + 2, i + 3);
	////		mat4.at<cv::uint8_t>(j / 2 + 1, i / 2) = mat.at<cv::uint8_t>(j + 3, i + 2);
	////		mat4.at<cv::uint8_t>(j / 2 + 1, i / 2 + 1) = mat.at<cv::uint8_t>(j + 3, i + 3);
	////	}
	////}
	////for (int j = 0; j < h-1; j+=4) {
	////	for (int i = 0; i < w-1; i+=4) {
	////		mat1.at<cv::Vec3b>(j / 2, i / 2) = mat.at<cv::Vec3b>(j, i);
	////		mat1.at<cv::Vec3b>(j / 2, i / 2 + 1) = mat.at<cv::Vec3b>(j, i + 1);
	////		mat1.at<cv::Vec3b>(j / 2 + 1, i / 2) = mat.at<cv::Vec3b>(j + 1, i);
	////		mat1.at<cv::Vec3b>(j / 2 + 1, i / 2 + 1) = mat.at<cv::Vec3b>(j + 1, i + 1);
	////		mat2.at<cv::Vec3b>(j / 2, i / 2) = mat.at<cv::Vec3b>(j , i + 2);
	////		mat2.at<cv::Vec3b>(j / 2, i / 2 + 1) = mat.at<cv::Vec3b>(j , i + 3);
	////		mat2.at<cv::Vec3b>(j / 2 + 1, i / 2) = mat.at<cv::Vec3b>(j + 1 , i + 2);
	////		mat2.at<cv::Vec3b>(j / 2 + 1, i / 2 + 1) = mat.at<cv::Vec3b>(j + 1 , i + 3);
	////		mat3.at<cv::Vec3b>(j / 2, i / 2) = mat.at<cv::Vec3b>(j + 2, i);
	////		mat3.at<cv::Vec3b>(j / 2, i / 2 + 1) = mat.at<cv::Vec3b>(j + 2, i + 1);
	////		mat3.at<cv::Vec3b>(j / 2 + 1, i / 2) = mat.at<cv::Vec3b>(j + 3, i);
	////		mat3.at<cv::Vec3b>(j / 2 + 1, i / 2 + 1) = mat.at<cv::Vec3b>(j + 3, i + 1);
	////		mat4.at<cv::Vec3b>(j / 2, i / 2) = mat.at<cv::Vec3b>(j + 2, i + 2);
	////		mat4.at<cv::Vec3b>(j / 2, i / 2 + 1) = mat.at<cv::Vec3b>(j + 2, i + 3);
	////		mat4.at<cv::Vec3b>(j / 2 + 1, i / 2) = mat.at<cv::Vec3b>(j + 3, i + 2);
	////		mat4.at<cv::Vec3b>(j / 2 + 1, i / 2 + 1) = mat.at<cv::Vec3b>(j + 3, i + 3);
	////	}
	////}
	//cv::Mat mat11, mat22, mat33, mat44;
	//cv::cvtColor(mat1, mat11, cv::COLOR_BayerRG2RGB);
	//cv::cvtColor(mat2, mat22, cv::COLOR_BayerRG2RGB);
	//cv::cvtColor(mat3, mat33, cv::COLOR_BayerRG2RGB);
	//cv::cvtColor(mat4, mat44, cv::COLOR_BayerRG2RGB);
	//std::string fileName;
	//std::string name1 = "Covert1";
	//std::string name2 = "Covert2";
	//std::string name3 = "Covert3";
	//std::string name4 = "Covert4";
	//cv::imwrite("mat1.png", mat1);
	//cv::imwrite("mat2.png", mat2);
	//cv::imwrite("mat3.png", mat3);
	//cv::imwrite("mat4.png", mat4);
	//fileName = argv[1];
	////fileName = fileName.replace(fileName.size() - 4, fileName.size() - 4, name1);
	//int lastFolder = fileName.find_last_of("\\") + 1;
	//if (lastFolder < 0) {
	//	lastFolder = fileName.find_last_of("/") + 1;
	//}
	//fileName.insert(lastFolder, name1);
	//cv::imwrite(fileName, mat11);
	//std::cout << fileName << " Save." << std::endl;
	//fileName = argv[1];
	//cv::imwrite(fileName.insert(lastFolder, name2), mat22);
	//std::cout << fileName << " Save." << std::endl;
	//fileName = argv[1];
	//cv::imwrite(fileName.insert(lastFolder, name3), mat33);
	//std::cout << fileName << " Save." << std::endl;
	//fileName = argv[1];
	//cv::imwrite(fileName.insert(lastFolder, name4), mat44);
	//std::cout << fileName << " Save." << std::endl;
	////std::vector<cv::Mat> mat2(16);
	////for(int i=0;i<16;++i)
	////	mat2[i] = cv::Mat(h / 4, w / 4, mat.type());
	////for (int j = 0; j < h; ++j) {
	////	for (int i = 0; i < w; ++i) {
	////		int modi = i % 4;
	////		int modj = j % 4;
	////		mat2[modj * 4 + modi].at<cv::Vec3b>(j / 4, i / 4) = mat.at<cv::Vec3b>(j, i);
	////	}
	////}
	////for (int i = 0; i < 16; ++i) {
	////	std::stringstream ss;
	////	ss << "mat" << i << ".png";
	////	cv::imwrite(ss.str(), mat2[i]);
	////}
	return 0;
}


void decompose_image_file_name(std::string imgFileName, std::string& prefix, int& paddings, int& startFrame, std::string& ext)
{
	size_t found = imgFileName.find_last_of('.');
	ext = imgFileName.substr(found + 1);

	int count = 0;
	for (int i = found - 1; i >= 0; --i)
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