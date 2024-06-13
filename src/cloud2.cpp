#include "cloud2.h"

#include "reconst.h"

#include <fstream>
#include <map>
#include <cstdlib>
#include <ctime>
#include <utility>

#include <opencv2/opencv.hpp>

#include <Eigen/Dense>

float bilinear_interpolation(float* f, int index, int w, float s, float t)
{
	return f[index] * (1.0 - s)*(1.0 - t) + f[index + 1] * s*(1.0 - t)
		+ f[index + w] * (1.0 - s)*t + f[index + w + 1] * s*t;
}

void PointCloud2::computeTSDF(int width, int height, float* rawDepth, float* KData, float* TData,
							  int blockSize, float* blockCenterData, float& voxelSize, 
							  int truncSize, float* tsdf)
{
	// Compute the bounding box for the point cloud
	// and set the block information
	mg::Vector3f tbmin, tbmax;
	computeBoundingBox(tbmin, tbmax);
	Eigen::Vector3f bmin(tbmin.data()), bmax(tbmax.data());

	Eigen::Map<Eigen::Vector3f> blockCenter(blockCenterData);
	blockCenter = 0.5f*(bmin + bmax);

	// Extend the bounding box by some (e.g. 10%)
	Eigen::Vector3f delta = bmax - bmin;
	bmin = blockCenter - 0.55*delta;
	bmax = blockCenter + 0.55*delta;
	std::cout << "- Bounding box : (" << bmin.transpose() << ")~(" << bmax.transpose() << ")" << std::endl;

	delta *= 1.1f;
	float blockLength = std::max(std::max(delta.x(), delta.y()), delta.z());
	voxelSize = blockLength / blockSize;
	delta = Eigen::Vector3f(blockLength, blockLength, blockLength);

	bmin = blockCenter - 0.5f*delta;
	bmax = blockCenter + 0.5f*delta;

	std::cout << "- Block Info :\n"
		<< " . center : (" << blockCenter.transpose() << ")\n"
		<< " . size : " << blockSize << "x" << blockSize << "x" << blockSize << "(length=" << blockLength << ")\n"
		<< " . voxel size : " << voxelSize << std::endl;

	// Set calibration and projection matrices
	Eigen::Matrix3f K(KData);
	Eigen::Matrix4f T(TData);

	Eigen::Matrix<float, 3, 4> Proj = K * T.block<3, 4>(0, 0);
	std::cout << "- Projection :\n" << Proj << std::endl;

	// Compute the truncated signed distance funtions on the voxel centers
	float INF = std::numeric_limits<float>::infinity();
	float eps = std::numeric_limits<float>::epsilon();

	float truncDist = voxelSize * truncSize;

	std::fill(tsdf, tsdf + blockSize * blockSize*blockSize, -1.0f);

	Eigen::Vector3f p;
	Eigen::Vector3f voxel(voxelSize, voxelSize, voxelSize);
	unsigned int count = 0, count1 = 0, count2 = 0;
	for (int k = 0; k < blockSize; ++k)
	{
		for (int j = 0; j < blockSize; ++j)
		{
			for (int i = 0; i < blockSize; ++i)
			{
				p = voxelSize * Eigen::Vector3f(i, j, k) + bmin + 0.5f*voxel;

				Eigen::Vector3f q = Proj * Eigen::Vector4f(p(0), p(1), p(2), 1.0f);

				float z = q(2);

				q /= z;

				int u, v;
				u = (int)floor(q(0));
				v = (int)floor(q(1));

				if (   z < INF && z > -INF && u > 0 && u < width-1 && v > 0 && v < height-1 
					&& rawDepth[v * width + u] > 0.0f && rawDepth[v * width + u] < INF)
				{
					// NOTE : We are using bilinear interpolation for smoother depthmap
					float du = q(0) - (float)u;
					float dv = q(1) - (float)v;

					float depth = bilinear_interpolation(rawDepth, v*width + u, width, du, dv);

					if (fabs(depth - z) < truncDist)
					{
						float f = (depth - z) / truncDist; // This guarantees |f|<=1

						tsdf[k*blockSize*blockSize + j * blockSize + i] = f;
					}
					else
					{
						if (depth-z > truncDist)
							tsdf[k*blockSize*blockSize + j * blockSize + i] = 1.0f;
					}
				}
				//else
				//{
				//	count++;
				//}
			}// end for(i)
		}// end for(j)
	}// end for(k)

	//std::cout << "- # of voxels with  improper truncated signed distances =" << count << " " << count1 << std::endl;
}

void PointCloud2::convertToRawDepth(int width, int height, float* rawDepth, float* K, float* camPose, bool holeFill, const char* filename)
{
	Eigen::Matrix4f CP(camPose);
	Eigen::Matrix4f T = CP.inverse();

	point2rawdepth(pos_.size(), pos_.data()->data(), K, T.data(), width, height, rawDepth, holeFill, filename);
}

void DepthCloud::computeCloudNormalInLocal(unsigned int numPts, 
									Point3fArray& p_local,
									std::map<unsigned int, unsigned int>& pixId2ptId)
{
	const int w = width_;
	const int h = height_;

	normal_.resize(numPts);
	std::fill(normal_.begin(), normal_.end(), mg::Vector3f(0, 0, -1));

	std::map<unsigned int, unsigned int>::iterator it;
	for (it = pixId2ptId.begin(); it != pixId2ptId.end(); ++it)
	{
		unsigned int center = it->first;

		unsigned int vid[9];
		// 0 - 1 - 2
		// 3 - 4 - 5
		// 6 - 7 - 8
		vid[0] = center + w - 1;
		vid[1] = center + w;
		vid[2] = center + w + 1;

		vid[3] = center - 1;
		vid[4] = center;
		vid[5] = center + 1;

		vid[6] = center - w - 1;
		vid[7] = center - w;
		vid[8] = center - w + 1;

		std::vector<int> nid(9);
		int count = 0;
		for (int k = 0; k < 9; ++k)
		{
			if (pixId2ptId.find(vid[k]) != pixId2ptId.end())
			{
				nid[count] = pixId2ptId[vid[k]];
				count++;
			}
		}
		nid.resize(count);

		if (count > 4)
		{

			Eigen::Vector3f centroid(0.0f, 0.0f, 0.0f);
			for (int k = 0; k < count; ++k)
			{
				centroid += Eigen::Vector3f(p_local[nid[k]].data());
			}
			centroid /= (float)count;

			Eigen::Matrix3f A = Eigen::Matrix3f::Zero();
			for (int k = 0; k < count; ++k)
			{
				Eigen::Vector3f dp = Eigen::Vector3f(p_local[nid[k]].data()) - centroid;

				A += dp * dp.transpose();
			}

			A /= (float)count;

			Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigenSolve(A);
			Eigen::Vector3f normal = eigenSolve.eigenvectors().col(0);

			//std::cout << it->first << "," << it->second << " :" << eigenSolve.eigenvalues().transpose() << std::endl;

			if (normal(2) > 0.0)
			{
				normal *= -1.0f;
			}

			normal_[it->second] = mg::Vector3f(normal.data());
		}
		/*else
		{
			std::cout << "Neighbor count=" << count << std::endl;
		}*/
	}
}

void DepthCloud::computeCloudInLocal(unsigned int& numPts,
									std::map<unsigned int, unsigned int>& pixId2ptId)
{
	const unsigned int w = width_;
	const unsigned int h = height_;

	pos_.resize(w*h);

	unsigned short * data = (unsigned short*)depth_;

	Eigen::Matrix3f invK(invK_);

	for (unsigned int j = 0; j < h; ++j)
	{
		for (unsigned int i = 0; i < w; ++i)
		{
			unsigned int pixId = j * w + i;
			unsigned short d = data[pixId];

			if (d > 0)
			{
				float z = 0.001f* data[pixId];

				Eigen::Vector3f p((float)i, (float)j, 1.0f);

				Eigen::Vector3f x = z * (invK*p);

				pos_[numPts] = mg::Vector3f(x.data());

				pixId2ptId[pixId] = numPts;

				numPts++;
			}

		}// end for(i)
	}// end for(j)

	pos_.resize(numPts);
}

void DepthCloud::setTextureToCloud(unsigned int numPts, std::map<unsigned int, unsigned int>& pixId2ptId)
{
	if (texture_ == NULL)
	{
		col_.resize(numPts);
		std::fill(col_.begin(), col_.end(), mg::Vector3f(255, 255, 255));
	}
	else
	{
		col_.resize(numPts);

		unsigned char* texData = (unsigned char*)texture_;
		std::map<unsigned int, unsigned int>::iterator it;
		for (it = pixId2ptId.begin(); it != pixId2ptId.end(); ++it)
		{
			col_[it->second].x = texData[3 * it->first + 0];
			col_[it->second].y = texData[3 * it->first + 1];
			col_[it->second].z = texData[3 * it->first + 2];
		}
	}
}

void DepthCloud::setTextureToCloud(unsigned int numPts, std::map<unsigned int, unsigned int>& pixId2ptId,
	unsigned char* textureData, int channels)
{
	if (textureData == NULL)
	{
		col_.resize(numPts);
		std::fill(col_.begin(), col_.end(), mg::Vector3f(255, 255, 255));
	}
	else
	{
		col_.resize(numPts);

		if (channels >= 3)
		{
			std::map<unsigned int, unsigned int>::iterator it;
			for (it = pixId2ptId.begin(); it != pixId2ptId.end(); ++it)
			{
				col_[it->second].x = textureData[channels*it->first + 0];
				col_[it->second].y = textureData[channels*it->first + 1];
				col_[it->second].z = textureData[channels*it->first + 2];
			}
		}
		else if (channels == 1)
		{
			std::map<unsigned int, unsigned int>::iterator it;
			for (it = pixId2ptId.begin(); it != pixId2ptId.end(); ++it)
			{
				col_[it->second].x = textureData[it->first];
				col_[it->second].y = textureData[it->first];
				col_[it->second].z = textureData[it->first];
			}
		}
		else
		{
			std::fill(col_.begin(), col_.end(), mg::Vector3f(255, 255, 255));
		}
	}
}

void DepthCloud::transform4x4ToCloud(float* trans)
{
	unsigned int numPts = (unsigned int)pos_.size();

	Eigen::Matrix4f T(trans);
	Eigen::Matrix3f R = T.block<3, 3>(0, 0);

	for (unsigned int i = 0; i < numPts; ++i)
	{
		Eigen::Vector4f x = T * Eigen::Vector4f(pos_[i].x, pos_[i].y, pos_[i].z, 1.0f);
		pos_[i] = mg::Vector3f(x.data());

		Eigen::Vector3f n_w = R * Eigen::Vector3f(normal_[i].data());
		normal_[i] = mg::Vector3f(n_w.data());
	}
}

void DepthCloud::depth2cloud(bool inWorldCoor)
{
	const int w = width_;
	const int h = height_;

	if (dataType_ == CV_8U)
	{

	}
	else if (dataType_ == CV_16U)
	{
		unsigned int numPts = 0;
		std::map<unsigned int, unsigned int> pixId2ptId;

		// Compute point cloud from depth data
		computeCloudInLocal(numPts, pixId2ptId);

		// Set texture to cloud
		setTextureToCloud(numPts, pixId2ptId);

		// Compute normals
		computeCloudNormalInLocal(numPts, pos_, pixId2ptId);

		if (inWorldCoor)
		{
			Eigen::Matrix<float, 3, 4> cp(cam_pose_);
			Eigen::Matrix4f local2World = Eigen::Matrix4f::Identity();
			local2World.block<3, 4>(0, 0) = cp;

			Eigen::Matrix4f T = Eigen::Matrix4f::Identity();
			T(1, 1) = T(2, 2) = -1.0;

			local2World *= T;

			transform4x4ToCloud(local2World.data());
		}

	}
}
void DepthCloud::depth2cloud(unsigned char* textureData, int channels, bool inWorldCoor)
{
	const int w = width_;
	const int h = height_;

	if (dataType_ == CV_8U)
	{

	}
	else if (dataType_ == CV_16U)
	{
		unsigned int numPts = 0;
		std::map<unsigned int, unsigned int> pixId2ptId;

		// Compute point cloud from depth data
		computeCloudInLocal(numPts, pixId2ptId);

		// Set texture to cloud
		setTextureToCloud(numPts, pixId2ptId, textureData, channels);

		// Compute normals
		computeCloudNormalInLocal(numPts, pos_, pixId2ptId);

		if (inWorldCoor)
		{
			Eigen::Matrix<float, 3, 4> cp(cam_pose_);
			Eigen::Matrix4f local2World = Eigen::Matrix4f::Identity();
			local2World.block<3, 4>(0, 0) = cp;

			Eigen::Matrix4f T = Eigen::Matrix4f::Identity();
			T(1, 1) = T(2, 2) = -1.0;

			local2World *= T;

			transform4x4ToCloud(local2World.data());
		}
	}
}

bool DepthCloud::loadTextureImage(const char* filename)
{
	cv::Mat texImg = cv::imread(filename);

	if (texImg.empty())
	{
		std::cout << "Error, can't open texture image : " << filename << std::endl;
		return false;
	}

	freeTextureData();

	const unsigned int tw = texImg.cols;
	const unsigned int th = texImg.rows;

	texture_ = malloc(sizeof(unsigned char) * 3 * tw*th);

	int channels = texImg.channels();

	if (channels == 3)
	{
		memcpy(texture_, texImg.data, sizeof(unsigned char) * 3 * tw*th);
	}
	else if (channels == 4)
	{
		unsigned char* texData = (unsigned char*)texture_;
		for (unsigned int i = 0; i < tw*th; ++i)
		{
			texData[3 * i + 0] = texImg.data[channels*i + 0];
			texData[3 * i + 1] = texImg.data[channels*i + 1];
			texData[3 * i + 2] = texImg.data[channels*i + 2];
		}
	}

	return true;
}

bool DepthCloud::loadCamPoseData(const char* filename)
{
	std::ifstream inFile(filename);
	if (!inFile.is_open())
	{
		std::cout << "Error, can't open calibration data : " << filename << std::endl;
		return false;
	}

	// Column major order!
	inFile
		>> cam_pose_[0] >> cam_pose_[3] >> cam_pose_[6] >> cam_pose_[9]
		>> cam_pose_[1] >> cam_pose_[4] >> cam_pose_[7] >> cam_pose_[10]
		>> cam_pose_[2] >> cam_pose_[5] >> cam_pose_[8] >> cam_pose_[11];

	inFile.close();

	Eigen::Matrix<float, 3, 4> camPose(cam_pose_);
	//std::cout << "Camera Pose : " << "\n" << camPose << std::endl;

	Eigen::Matrix4f extCamPose = Eigen::Matrix4f::Identity();
	extCamPose.block<3, 4>(0, 0) = camPose;

	Eigen::Matrix4f invExtCamPose = extCamPose.inverse();
	Eigen::Map<Eigen::Matrix<float, 3, 4>> invCamPose(inv_cam_pose_);
	invCamPose = invExtCamPose.block<3, 4>(0, 0);

	//std::cout << "Inverse Camera Pose = \n" << invCamPose << std::endl;

	inFile.close();

	return true;
}

bool DepthCloud::loadCalibrData(const char* filename)
{
	std::ifstream inFile(filename);
	if (!inFile.is_open())
	{
		std::cout << "Error, can't open calibration data : " << filename << std::endl;
		return false;
	}

	// Column major order!
	inFile 
		>> K_[0] >> K_[3] >> K_[6]
		>> K_[1] >> K_[4] >> K_[7]
		>> K_[2] >> K_[5] >> K_[8];

	inFile.close();

	Eigen::Matrix3f K(K_);
	//std::cout << "K : \n" << K << std::endl;

	Eigen::Map<Eigen::Matrix3f> invK(invK_);
	invK = K.inverse();

	//std::cout << "Inverse K = \n" << invK << std::endl;

	inFile.close();

	Eigen::Map<Eigen::Matrix<float, 3, 4>> camPose(cam_pose_);
	Eigen::Matrix4f I = Eigen::Matrix4f::Identity();
	camPose = I.block<3, 4>(0, 0);

	return true;
}

void DepthCloud::setDepthImage(unsigned int w, unsigned int h, void* depth, int type)
{
	width_ = w;
	height_ = h;
	dataType_ = type;

	if (type == CV_8U)
	{
		freeDepthData();
		depth_ = malloc(sizeof(unsigned char)*w*h);
		memcpy(depth_, depth, sizeof(unsigned char)*w*h);
	}
	else if (dataType_ == CV_16U)
	{
		freeDepthData();
		depth_ = malloc(sizeof(unsigned short)*w*h);
		memcpy(depth_, depth, sizeof(unsigned short)*w*h);
	}
	else
	{
		std::cout << "Error, unknown data type : " << type << std::endl;
	}
}

bool DepthCloud::loadDepthImage(const char* filename)
{
	cv::Mat depthImg = cv::imread(filename, -1);
	if (depthImg.empty())
	{
		std::cout << "Error, can't load depth image : " << filename << std::endl;
		return false;
	}

	dataType_ = depthImg.depth();
	unsigned int w = depthImg.cols;
	unsigned int h = depthImg.rows;

	if (dataType_ == CV_8U)
	{
		freeDepthData();
		depth_ = malloc(sizeof(unsigned char)*w*h);
		memcpy(depth_, depthImg.data, sizeof(unsigned char)*w*h);
	}
	else if (dataType_ == CV_16U)
	{
		freeDepthData();
		depth_ = malloc(sizeof(unsigned short)*w*h);
		memcpy(depth_, depthImg.data, sizeof(unsigned short)*w*h);
	}
	else
	{
		std::cout << "Error, depth image should be in 8 or 16 bits." << std::endl;
		return false;
	}

	width_ = w;
	height_ = h;

	return true;
}

bool DepthCloud::loadData(const char* filename)
{
	return loadDepthImage(filename);
}

void DepthCloud::freeDepthData()
{
	if (depth_)
	{
		free(depth_);
		depth_ = NULL;
	}
}

void DepthCloud::freeTextureData()
{
	if (texture_)
	{
		free(texture_);
		texture_ = NULL;
	}
}

//void DepthCloud::setDataType(int type)
//{
//	dataType_ = type; 
//}

//int DepthCloud::getDataType()
//{
//	return dataType_;
//}

//unsigned int DepthCloud::width()
//{
//	return width_;
//}

//unsigned int DepthCloud::height()
//{
//	return height_;
//}

//void*  DepthCloud::getDepthData()
//{
//	return depth_;
//}
DepthCloud::~DepthCloud()
{
	freeDepthData();
	freeTextureData();
}

DepthCloud::DepthCloud() : depth_(NULL), width_(0), height_(0),
	texture_(NULL), tex_width_(0), tex_height_(0)
{

}

void PointCloud2::initDiffuse()
{
	unsigned int n = (unsigned int)pos_.size();

	diffuse_.resize(n);

	std::fill(diffuse_.begin(), diffuse_.end(), mg::Vector3f(0.0f, 0.0f, 0.0f));
}

void PointCloud2::applyDirLightToDiffuse(mg::Vector3f lv, float mag)
{
	unsigned int n = (unsigned int)pos_.size();

	diffuse_.resize(n);

	for (unsigned int i = 0; i < n; ++i)
	{
		float lvn = lv.dotProduct(normal_[i]);

		if (lvn > 0.0f)
		{
			diffuse_[i] += (mag * lvn)*col_[i];
		}
	}
}

void  PointCloud2::applyDirLightToSpecular(mg::Vector3f lv, float mag)
{
	unsigned int n = (unsigned int)pos_.size();

	specular_[0].resize(n);
	specular_[1].resize(n);
	specular_[2].resize(n);

	for (unsigned int i = 0; i < n; ++i)
	{
		float lvn = lv.dotProduct(normal_[i]);

		specular_[0][i] = mg::Vector3f(0.0f, 0.0f, 0.0f);
		specular_[1][i] = mg::Vector3f(0.0f, 0.0f, 0.0f);
		specular_[2][i] = mg::Vector3f(0.0f, 0.0f, 0.0f);

		if (lvn > 0.0f)
		{
			specular_[0][i] = mag * col_[i].x*(-lv + 2.0f*lvn*normal_[i]);
			specular_[1][i] = mag * col_[i].y*(-lv + 2.0f*lvn*normal_[i]);
			specular_[2][i] = mag * col_[i].z*(-lv + 2.0f*lvn*normal_[i]);
		}
	}
}

void  PointCloud2::setAmbient(float mag)
{
	unsigned int n = (unsigned int)pos_.size();

	ambient_.resize(n);

	for (unsigned int i = 0; i < n; ++i)
	{
		ambient_[i] = mag * col_[i];
	}
}

void  PointCloud2::setPhongShadeFactor(float diff, float spec, float ambient)
{
	kd_ = diff;
	ks_ = spec;
	ka_ = ambient;
}

void PointCloud2::setSpecularExponent(float exp)
{
	alpha_ = exp;
}

bool PointCloud2::exportToPly(const char* filename, bool color, bool normal)
{
	std::ofstream outFile(filename, std::ios::binary);

	if (!outFile.is_open())
	{
		std::cerr << "Error, can't open " << filename << " for writing" << std::endl;
		return false;
	}

	outFile << "ply" << "\n"
		<< "format binary_little_endian 1.0" << "\n"
		<< "comment" << "\n"
		<< "comment" << "\n"
		<< "element vertex " << pos_.size() << "\n"
		<< "property float x" << "\n"
		<< "property float y" << "\n"
		<< "property float z" << "\n"
		<< "property uchar red " << "\n"
		<< "property uchar green" << "\n"
		<< "property uchar blue" << "\n";

	if (normal)
	{
		outFile
			<< "property float nx" << "\n"
			<< "property float ny" << "\n"
			<< "property float nz" << "\n";
	}

	outFile << "end_header" << "\n";

	unsigned char c[3];

	if (color && normal)
	{
		for (int i = 0; i < pos_.size(); ++i)
		{
			outFile.write((char*)&pos_[i], sizeof(mg::Vector3f));

			c[0] = (unsigned char)col_[i].z;
			c[1] = (unsigned char)col_[i].y;
			c[2] = (unsigned char)col_[i].x;

			outFile.write((char*)c, sizeof(unsigned char) * 3);

			outFile.write((char*)&normal_[i], sizeof(mg::Vector3f));
		}
	}
	else if (color)
	{
		for (int i = 0; i < pos_.size(); ++i)
		{
			outFile.write((char*)&pos_[i], sizeof(mg::Vector3f));

			c[0] = (unsigned char)col_[i].z;
			c[1] = (unsigned char)col_[i].y;
			c[2] = (unsigned char)col_[i].x;

			outFile.write((char*)c, sizeof(unsigned char) * 3);
		}
	}
	else if(normal)
	{
		for (int i = 0; i < pos_.size(); ++i)
		{
			outFile.write((char*)&pos_[i], sizeof(mg::Vector3f));

			c[0] = 255.0f;
			c[1] = 255.0f;
			c[2] = 255.0f;

			outFile.write((char*)c, sizeof(unsigned char) * 3);

			outFile.write((char*)&normal_[i], sizeof(mg::Vector3f));
		}
	}
	else
	{
		for (int i = 0; i < pos_.size(); ++i)
		{
			outFile.write((char*)&pos_[i], sizeof(mg::Vector3f));

			c[0] = 255.0f;
			c[1] = 255.0f;
			c[2] = 255.0f;

			outFile.write((char*)c, sizeof(unsigned char) * 3);
		}
	}

	outFile.close();

	return true;
}

void PointCloud2::reflectionX(float x0)
{
	unsigned int N = getSize();
	for (unsigned int i = 0; i < N; ++i)
	{
		pos_[i].x = 2.0f*x0 - pos_[i].x;
		normal_[i].x *= -1.0f;
	}
}

void PointCloud2::reflectionY(float y0)
{
	unsigned int N = getSize();
	for (unsigned int i = 0; i < N; ++i)
	{
		pos_[i].y = 2.0f*y0 - pos_[i].y;
		normal_[i].y *= -1.0f;
	}
}

void PointCloud2::reflectionZ(float z0)
{
	unsigned int N = getSize();
	for (unsigned int i = 0; i < N; ++i)
	{
		pos_[i].z = 2.0f*z0 - pos_[i].z;
		normal_[i].z *= -1.0f;
	}
}

//void PointCloud2::perturb_z_depths(int pert_method, float lambda)
//{
//	if (pert_method == MULTIPLE_WAVELENGTH)
//	{
//		unsigned int N = getSize();
//		for (unsigned int i = 0; i < N; ++i)
//		{
//			pos_[i].z = (int)floor(pos_[i].z / lambda - 0.5f)*lambda;
//		}
//	}
//	else if (pert_method == RANDOM_NOISES)
//	{
//		unsigned int N = getSize();
//
//		srand((unsigned int)time(NULL));
//
//		for (unsigned int i = 0; i < N; ++i)
//		{
//			pos_[i].z = pos_[i].z + (float)rand() / RAND_MAX * lambda;
//		}
//	}
//	else
//	{
//		std::cout << "Error, unknown perturbation method" << std::endl;
//		return;
//	}
//}

mg::Vector3f PointCloud2::computeCenter()
{
	mg::Vector3f c(0.0, 0.0, 0.0);

	unsigned int n = pos_.size();

	for (unsigned int i = 0; i < n; ++i)
	{
		c += pos_[i];
	}

	c /= (float)n;

	return c;
}

void PointCloud2::computeBoundingBox(mg::Vector3f& bmin, mg::Vector3f& bmax)
{
	float INF = std::numeric_limits<float>::infinity();
	bmin = mg::Vector3f(INF, INF, INF);
	bmax = -bmin;

	unsigned int numPts = (unsigned int)pos_.size();

	for (unsigned int i = 0; i < numPts; ++i)
	{
		if (pos_[i].x < bmin.x)	bmin.x = pos_[i].x;
		if (pos_[i].y < bmin.y)	bmin.y = pos_[i].y;
		if (pos_[i].z < bmin.z)	bmin.z = pos_[i].z;

		if (pos_[i].x > bmax.x)	bmax.x = pos_[i].x;
		if (pos_[i].y > bmax.y)	bmax.y = pos_[i].y;
		if (pos_[i].z > bmax.z)	bmax.z = pos_[i].z;
	}
}

void  PointCloud2::setSize(unsigned int n)
{
	pos_.resize(n);
	normal_.resize(n);
	col_.resize(n);

	ambient_.resize(n);
	diffuse_.resize(n);
	specular_[0].resize(n);
	specular_[1].resize(n);
	specular_[2].resize(n);
}

unsigned int PointCloud2::getSize()
{
	return (unsigned int)pos_.size();
}


void PointCloud2::setPoints(unsigned int n, mg::Vector3f* p)
{
	setSize(n);
	memcpy(pos_.data(), p, sizeof(mg::Vector3f)*n);
}

void PointCloud2::setNormals(unsigned int n, mg::Vector3f* v)
{
	setSize(n);
	memcpy(normal_.data(), v, sizeof(mg::Vector3f)*n);
}

void PointCloud2::setColors(unsigned int n, mg::Vector3f* c)
{
	setSize(n);
	memcpy(col_.data(), c, sizeof(mg::Vector3f)*n);
}

void PointCloud2::setConstantColor(BGR c)
{
	unsigned int n = (unsigned int)pos_.size();

	for (unsigned int i = 0; i < n; ++i)
	{
		col_[i] = mg::Vector3f(c.b, c.g, c.r);
	}
}

void PointCloud2::rescaleColor(float factor)
{
	int n = (int)col_.size();
	for (int i = 0; i < n; ++i)
	{
		col_[i] *= factor;
	}
}

bool PointCloud2::loadDataInPly(const char* filename)
{
	std::cout << "- Reading " << filename << "...";

	std::ifstream inFile(filename, std::ios::binary);

	if (!inFile.is_open())
	{
		std::cerr << "Error, can't open PLY file for reading : " << filename << std::endl;
		return false;
	}

	std::stringstream buffer;
	buffer << inFile.rdbuf();

	std::string temp;

	int num_pts = 0;

	bool isBinary = true;
	bool hasPoint = false, hasColor = false, hasNormal = false;

	// Header
	std::vector<std::pair<std::string, std::string>>	vertexProperty;

	do
	{
		if (temp == "format")
		{
			buffer >> temp;
			if (temp == "ascii") isBinary = false;
		}
		else if (temp == "vertex")
		{
			buffer >> num_pts;

			setSize(num_pts);
		}
		else if (temp == "property")
		{
			if (num_pts > 0)
			{
				std::string type, name;
				buffer >> type;
				buffer >> name;

				vertexProperty.push_back(std::make_pair(type, name));

				if (name == "x")
				{
					hasPoint = true;
				}
				else if (name == "red")
				{
					hasColor = true;
				}
				else if (name == "nx")
				{
					hasNormal = true;
				}
			}
		}

		buffer >> temp;

	} while (temp != "end_header");

	std::cout << "num pts=" << num_pts << "...";
	std::cout << "position=" << hasPoint << ", color=" << hasColor << ", normal=" << hasNormal << std::endl;

	if (isBinary)
	{
		std::getline(buffer, temp, '\n');

		unsigned char cdata[3];
		for (int i = 0; i < num_pts; ++i)
		{
			if (hasPoint)
				buffer.read((char*)&pos_[i], sizeof(mg::Vector3f));
			if (hasColor)
			{
				buffer.read((char*)cdata, sizeof(unsigned char) * 3);
				col_[i] = mg::Vector3f(cdata[2], cdata[1], cdata[0]);
			}
			if (hasNormal)
				buffer.read((char*)&normal_[i], sizeof(mg::Vector3f));

			if (vertexProperty.size() > 9)
			{
				float curv;
				buffer.read((char*)&curv, sizeof(float));
			}
		}
	}
	else // ASCII
	{
		for (int i = 0; i < num_pts; ++i)
		{
			if (hasPoint)
				buffer >> pos_[i].x >> pos_[i].y >> pos_[i].z;
			if (hasColor)
				buffer >> col_[i].z >> col_[i].y >> col_[i].x;
			if (hasNormal)
				buffer >> normal_[i].x >> normal_[i].y >> normal_[i].z;

			if (vertexProperty.size() > 9)
			{
				float curv;
				buffer >> curv;
			}
		}
	}

	//std::cout << pos_[0] << " " << col_[0] << " " << normal_[0] << std::endl;
	//std::cout << pos_[1] << " " << col_[1] << " " << normal_[1] << std::endl;
	//std::cout << pos_[num_pts - 1] << " " << col_[num_pts - 1] << " " << normal_[num_pts - 1] << std::endl;

	inFile.close();

	return true;
}

bool PointCloud2::loadData(const char* filename)
{
	return loadDataInPly(filename);
}

void PointCloud2::setData(int n, double* p, unsigned char* c)
{
	setSize(n);

	for (int i = 0; i < n; ++i)
	{
		if (p != nullptr)
		{
			Eigen::Vector3d tp(p + 3 * i);
			pos_[i].x = tp.x();
			pos_[i].y = tp.y();
			pos_[i].z = tp.z();
		}

		if (c != nullptr)
		{
			col_[i].x = c[3 * i + 0];
			col_[i].y = c[3 * i + 1];
			col_[i].z = c[3 * i + 2];
		}
	}
}

PointCloud2::PointCloud2() : kd_(0.0f), ks_(0.0f), ka_(1.0f), alpha_(1.0f) {}