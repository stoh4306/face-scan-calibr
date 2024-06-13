#pragma once

#include <vector>

class KnnSearch
{
public:
	KnnSearch();
	KnnSearch(int npts, float* pts_data, int dim);
	~KnnSearch();

	void rebuild(int npts, float* pts_data);
	void nearestSearch(int n, float* query_pts, int max_nn, int* indices, float* dists);
	void radiusSearch(int n, float* query_pts, float radius,
		std::vector<std::vector<int>>& indices, std::vector<std::vector<float>>& dists);

	float find_aver_spacing();

public:
	int dim_;
	void * flann_index_;
	std::vector<float> data_;
};
