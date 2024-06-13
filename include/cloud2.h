#pragma once

#include "vec.h"
#include <vector>
#include <map>

typedef std::vector<mg::Vector3f> Point3fArray;
typedef std::vector<mg::Vector3d> Point3dArray;
typedef std::vector<mg::Vector3i> Point3iArray;

struct BGR { unsigned char b; unsigned char g; unsigned char r; };

class PointCloud2
{
public:
	PointCloud2();

	//! Set the number of points
	void	setSize(unsigned int n);

	//! Get the number of points
	unsigned int getSize();

	//! Set the points
	void	setPoints(unsigned int n, mg::Vector3f* p);

	//! Set the normals
	void	setNormals(unsigned int n, mg::Vector3f* v);

	//! Set the colors
	void	setColors(unsigned int n, mg::Vector3f* c);

	//! Set the color as a constant value
	void	setConstantColor(BGR c);

	//! Rescale the color values
	void	rescaleColor(float factor);

	//! Compute the bounding box of the points
	void	computeBoundingBox(mg::Vector3f& bmin, mg::Vector3f& bmax);

	//! Compute the center of mass
	mg::Vector3f computeCenter();

	//! Z-depths perturbation
	//void	perturb_z_depths(int pert_method, float lambda);

	//! Reflection 
	void	reflectionX(float x0);
	void	reflectionY(float y0);
	void	reflectionZ(float z0);

	virtual bool loadData(const char* filename);

	//! Apply the directional light with the diffusiness
	void	applyDirLightToDiffuse(mg::Vector3f lv, float mag);

	void	initDiffuse();

	void	applyDirLightToSpecular(mg::Vector3f lv, float mag);

	void	setAmbient(float mag);

	//! Set the phong shading factors
	void	setPhongShadeFactor(float diff, float spec, float ambient);

	void	setSpecularExponent(float exp);

	//! Export the point cloud data in PLY
	bool	exportToPly(const char* filename, bool color=true, bool normal=true);

	void 	convertToRawDepth(int width, int height, float* rawDepth, float* K, float* camPose,  bool holeFill, const char* filename="");

	void	computeTSDF(int width, int height, float* rawDepth, float* KData, float* TData, 
						int blockSize, float* blockCenter, float& voxelSize, int truncSize, float* tsdf);

public: // JK
	void setData(int n, double* p, unsigned char* c);

protected:
	bool loadDataInPly(const char* filename);

public:
	Point3fArray pos_;
	Point3fArray normal_;
	Point3fArray col_;

	Point3fArray ambient_;
	Point3fArray diffuse_;
	Point3fArray specular_[3];

	float kd_;		// diffuse 
	float ks_;		// specular
	float ka_;		// ambient
	float alpha_;	// specular exponent
};

class DepthCloud : public PointCloud2
{
public:
	DepthCloud();
	~DepthCloud();

	virtual bool loadData(const char* filename);
	bool	loadDepthImage(const char* filename);
	bool	loadCalibrData(const char* filename); 
	bool	loadCamPoseData(const char* filename);
	void	setDepthImage(unsigned int w, unsigned int h, void* depth, int type);
	void	freeDepthData();
	bool	loadTextureImage(const char* filename);
	void	freeTextureData();

	void	depth2cloud(unsigned char* textureData, int channels, bool inWorldCoor = false);
	void	depth2cloud(bool inWorldCoor = false);

	void	computeCloudInLocal(unsigned int& numPts, std::map<unsigned int, unsigned int>& pixId2ptId);
	void	computeCloudNormalInLocal(unsigned int numPts, 
					Point3fArray& p_local, std::map<unsigned int, unsigned int>& pixId2ptId);
	void	setTextureToCloud(unsigned int numPts, std::map<unsigned int, unsigned int>& pixId2ptId);
	void	setTextureToCloud(unsigned int numPts, std::map<unsigned int, unsigned int>& pixId2ptId,
								unsigned char* textureData, int channels);
	void	transform4x4ToCloud(float* trans);

protected:
		
public:
	void* depth_;
	int	dataType_; // CV_8U or CV_16U
	unsigned int width_, height_;

	void* texture_;
	int textureDataType_;
	unsigned int tex_width_, tex_height_;

	float K_[9], invK_[9];
	float cam_pose_[12], inv_cam_pose_[12]; // orientation, center
};