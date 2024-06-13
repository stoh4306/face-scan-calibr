#include "data_impot_export.h"

#include <fstream>

bool import_matrix(const char* filename, int rows, int cols, double* data)
{
	if (filename == nullptr)
	{
		// Error Message
		return false;
	}

	std::ifstream inFile(filename);
	if (!inFile.is_open())
	{
		// Error Message
		return false;
	}

	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols; ++j)
		{
			inFile >> data[i + (j * rows)];
		}
	}
	inFile.close();

	return true;
}
void export_matrix(const char* filename, int rows, int cols, double* data)
{
	if (filename == nullptr)
	{
		// Error Message
		return;
	}

	std::ofstream outFile(filename);

	for (int i = 0; i < rows; ++i)
	{
		for (int j = 0; j < cols; ++j)
		{
			outFile << data[i + (rows * j)] << " ";
		}
		outFile << std::endl;
	}
}

bool import_2d_points(const char* filename, std::vector<Eigen::Vector2d>* pts)
{
	if (filename == nullptr)
		return false;

	std::ifstream inFile(filename);
	if (!inFile.is_open())
		return false;

	int pts_cnt;
	inFile >> pts_cnt;

	pts->clear();
	for (int i = 0; i < pts_cnt; ++i)
	{
		Eigen::Vector2d point;
		inFile >> point(0) >> point(1);
		pts->push_back(point);
	}
	inFile.close();
	return true;
}
void export_2d_points(const char* filename, std::vector<Eigen::Vector2d>* pts)
{
	if (filename == nullptr || pts == nullptr)
		return;

	size_t pts_cnt = pts->size();

	std::ofstream outFile(filename);
	outFile
		<< pts_cnt << std::endl;

	for (size_t i = 0; i < pts_cnt; ++i)
	{
		outFile << (*pts)[i].x() << " " << (*pts)[i].y() << std::endl;
	}
	outFile.close();

}
bool import_3d_points(const char* filename, std::vector<Eigen::Vector3d>* pts)
{
	if (filename == nullptr)
		return false;

	std::ifstream inFile(filename);
	if (!inFile.is_open())
		return false;

	int pts_cnt;
	inFile >> pts_cnt;

	pts->clear();
	for (int i = 0; i < pts_cnt; ++i)
	{
		Eigen::Vector3d point;
		inFile >> point(0) >> point(1) >> point(2);
		pts->push_back(point);
	}
	inFile.close();
	return true;
}
void export_3d_points(const char* filename, std::vector<Eigen::Vector3d>* pts)
{
	if (filename == nullptr || pts == nullptr)
		return;

	size_t pts_cnt = pts->size();

	std::ofstream outFile(filename);
	outFile << pts_cnt << std::endl;
	for (size_t i = 0; i < pts_cnt; ++i)
	{
		outFile << (*pts)[i].x() << " " << (*pts)[i].y() << " " << (*pts)[i].z() << std::endl;
	}
	outFile.close();
}

void export_volume_data(const char* filename, float* center, unsigned int nx, unsigned int ny, unsigned int nz, float vsize, std::vector<float>& f)
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