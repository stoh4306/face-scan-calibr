#ifndef MG_3D_RECONST_H_
#define MG_3D_RECONST_H_

#include<string>
void warpImage(int w, int h, int type, void* srcImg, void* tgtImg, double* pH);
void rectify(double* tK1, double* tP1, double B, int numPts, double* x2, double* X,
    double* tK2, double* tP2,
    double* H, int w, int h, int type, void* left_img, void* right_img, void* rect_right_img,std::string folder = "");
void rectify2(double* tK1, double* tP1, double B, int numPts, double* x2,
    double* tK2, double* camPose2, double* tP2,
    double* H, int w, int h, int type, void* left_img, void* right_img, void* rect_right_img,std::string folder = "");
void stereo_matching(int w, int h, int type, void* left_img_data, void* right_img_data,
    int blockSize, int maxDisparity, float* disp, const char* out_disparity_file);
void stereo_matching(int w, int h, int type, void* left_img_data, void* right_img_data,
    int blockSize, int maxDisparity, double* disp, const char* out_disparity_file);
void linear_triangulation(double* x1, double* x2, double* P1, double* P2, double* X);

void point2rawdepth(int n, float* p, float* K_data, float* T_data, int width, int height, float* rawDepth, bool holeFilling = true, const char* filename="");

void remedy_error_region_by_feature_tracking(int num_corr, double* p, double* q,
	int w, int h, double* K, float* base_depth, float* input_depth, 
	unsigned char* error_region_mask, float* output_depth,
	bool updateColor, unsigned char* base_color, unsigned char* input_color, unsigned char* colorMask,
	float alpha, float beta);

void write_point_cloud_in_ply(const char* filename, int n, double* p);
void write_point_cloud_in_ply(const char* filename, int n, double* p, unsigned char* c);
void write_point_cloud_in_pcd(const char* filename, int n, double* p);
void write_point_cloud_in_pcd(const char* filename, int n, double* p, unsigned char* c);
void write_point_cloud_in_pcd(const char* filename, int n, double* p, unsigned char* c, double* normal);
#endif