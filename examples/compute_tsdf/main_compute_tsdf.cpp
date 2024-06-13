#include <iostream>
#include <vector>
#include <fstream>

#include "cloud2.h"
#include "calibr.h"

#include <Eigen/Dense>

void write_volume_data(const char* filename, float* center,
	unsigned int nx, unsigned int ny, unsigned int nz, float vsize,
	std::vector<float>& f)
{
	std::ofstream outFile(filename, std::ios::binary);

	float bmin[3], bmax[3];
	bmin[0] = bmin[1] = bmin[2] = -0.5f*nx*vsize;
	bmin[0] += center[0];
	bmin[1] += center[1];
	bmin[2] += center[2];

	bmax[0] = bmax[1] = bmax[2] = 0.5f*nx*vsize;
	bmax[0] += center[0];
	bmax[1] += center[1];
	bmax[2] += center[2];

	outFile.write((char*)bmin, sizeof(float) * 3);
	outFile.write((char*)bmax, sizeof(float) * 3);

	outFile.write((char*)&nx, sizeof(unsigned int));
	outFile.write((char*)&ny, sizeof(unsigned int));
	outFile.write((char*)&nz, sizeof(unsigned int));

	outFile.write((char*)f.data(), sizeof(float)*nx*ny*nz);

	outFile.close();
}


int main(int argc, char** argv)
{
	std::cout << "***************************************************" << "\n"
		<< " COMPUTE TRUNCATED SIGNED DISTANCE FUNCION" << "\n"
		<< "***************************************************" << std::endl;

	if (argc <= 8)
	{
		std::cout << "- Usage : " << argv[0] << " arg1 arg2 ... arg7" << "\n"
			<< " . arg1 : point cloud filename" << "\n"
			<< " . arg2 : camera matrix K filename " << "\n"
			<< " . arg3 : camera pose filename" << "\n"
			<< " . arg4 : width(depthmap)" << "\n"
			<< " . arg5 : height(depthmap)" << "\n"
			<< " . arg6 : number of voxels in each direction" << "\n"
			<< " . arg7 : truncated distance size in number of voxels" << "\n"
			<< " . arg8 : depthmap hole-filling(0: no, 1: yes)"
			<< std::endl;

		return 0;
	}
	const char* cloudFile	= argv[1];
	const char* K_file		= argv[2];
	const char* CP_file		= argv[3];
	const int	width		= atoi(argv[4]);
	const int	height		= atoi(argv[5]);
	const int	blockSize	= atoi(argv[6]);
	const int	truncSize	= atoi(argv[7]);
	const bool  holeFill	= atoi(argv[8]) > 0 ? true : false;
	
	PointCloud2 cloud;
	cloud.loadData(argv[1]);

	if (cloud.getSize() == 0)
	{
		std::cout << "- No points in the cloud file -> Nothings done!" << std::endl;
		return 0;
	}

	float K[9], CP[16];
	read_calibrations(argv[2], argv[3], K, CP);
	
	mg::Vector3f bmin, bmax;
	cloud.computeBoundingBox(bmin, bmax);
	if (bmin.z < 0.0)
	{
		std::cout << "- Campose is multiplied by diag(1, -1, -1, 1):\n";
		Eigen::Map<Eigen::Matrix4f> camPose(CP);
		//std::cout << camPose << "->";
		camPose.col(1) *= -1.0f;
		camPose.col(2) *= -1.0f;
		//std::cout << camPose << std::endl;
	}

	std::vector<float> rawDepth(width*height);
	std::cout << "- Calculating raw depth data...";
	cloud.convertToRawDepth(width, height, rawDepth.data(), K, CP, holeFill, "raw_depth.data");
	std::cout << "done" << std::endl;

	//test
	//DepthCloud depthCloud;
	//depthCloud.loadData("depth.png");
	//depthCloud.loadCalibrData(K_file);
	//depthCloud.loadCamPoseData(CP_file);
	//depthCloud.depth2cloud(false);
	//depthCloud.exportToPly("depth2cloud.ply");

	float xform_data[16];
	Eigen::Map<Eigen::Matrix4f> T(xform_data);
	Eigen::Matrix4f camPose(CP);
	T = camPose.inverse();
	//std::cout << "T=\n" << T << std::endl;

	float blockCenter[3];
	float voxelSize;
	std::vector<float> tsdf(blockSize*blockSize*blockSize);
	//std::cout << "- Computing tsdf...";
	cloud.computeTSDF(width, height, rawDepth.data(), K, T.data(), 
					  blockSize, blockCenter, voxelSize, truncSize, tsdf.data());
	//std::cout << "done" << std::endl;

	const char* tsdfFile = "tsdf.vox";
	std::cout << "- Writing TSDF..." << tsdfFile << "...";
	write_volume_data("tsdf.vox", blockCenter, blockSize, blockSize, blockSize, voxelSize, tsdf);
	std::cout << "done" << std::endl;

	return 0;
}