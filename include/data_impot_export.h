#pragma once

#include <vector>
#include <string>
#include <Eigen/Dense>


bool import_matrix(const char* filename, int rows, int cols, double* data);
void export_matrix(const char* filename, int rows, int cols, double* data);

bool import_2d_points(const char* filename, std::vector<Eigen::Vector2d>* pts);
void export_2d_points(const char* filename, std::vector<Eigen::Vector2d>* pts);
bool import_3d_points(const char* filename, std::vector<Eigen::Vector3d>* pts);
void export_3d_points(const char* filename, std::vector<Eigen::Vector3d>* pts);

void export_volume_data(const char* filename, float* center, unsigned int nx, unsigned int ny, unsigned int nz, float vsize, std::vector<float>& f);
