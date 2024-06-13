#include "knn_search.h"

#include <flann/flann.h>

float KnnSearch::find_aver_spacing()
{
	int max_nn = 2;

	int npts = data_.size() / dim_;
	flann::Matrix<float> query(new float[npts * dim_], npts, dim_);
	memcpy(query.ptr(), data_.data(), sizeof(float) * dim_ * npts);

	flann::Matrix<int> indices(new int[npts*max_nn], npts, max_nn);
	flann::Matrix<float> dists(new float[npts*max_nn], npts, max_nn);

	((flann::Index<flann::L2<float> >*)flann_index_)->knnSearch(query, indices, dists, max_nn, flann::SearchParams());

	double psep = 0.0;
	for (int i = 0; i < npts; ++i)
	{
		psep += sqrt(dists[i][1]);
	}

	if (npts > 0)
	{
		psep /= npts;
	}

	delete[] indices.ptr();
	delete[] dists.ptr();

	return psep;
}

void KnnSearch::radiusSearch(int n, float* query_pts, float radius,
	std::vector<std::vector<int>>& indices, std::vector<std::vector<float>>& dists)
{
	flann::Matrix<float> query(query_pts, n, dim_);

	((flann::Index<flann::L2<float> >*)flann_index_)->radiusSearch(query, indices, dists, radius, flann::SearchParams());
}

void KnnSearch::nearestSearch(int n, float* query_pts, int max_nn, int* t_indices, float* t_dists)
{
	flann::Matrix<float> query(query_pts, n, dim_);


	flann::Matrix<int> indices(t_indices, n, max_nn);
	flann::Matrix<float> dists(t_dists, n, max_nn);

	((flann::Index<flann::L2<float> >*)flann_index_)->knnSearch(query, indices, dists, max_nn, flann::SearchParams());
}

void KnnSearch::rebuild(int npts, float* pts_data)
{
	//delete[] data_.ptr();
	data_.resize(npts * dim_);
	memcpy(data_.data(), pts_data, sizeof(float) * dim_ * npts);

	flann::Matrix<float> data(data_.data(), npts, dim_);

	if (flann_index_)
	{
		flann::Index<flann::L2<float> >* temp_ptr = (flann::Index<flann::L2<float> >*)flann_index_;
		delete temp_ptr;
		flann_index_ = NULL;
	}	

	flann::Index<flann::L2<float> >* index_ptr = new flann::Index<flann::L2<float> >(data, flann::KDTreeIndexParams());
	flann_index_ = index_ptr;

	//std::cout << " . Building linear index..." << std::flush;
	index_ptr->buildIndex();
	//std::cout << "done" << std::endl;
}

KnnSearch::~KnnSearch()
{
	if(flann_index_)
	{
		flann::Index<flann::L2<float> >* temp_ptr = (flann::Index<flann::L2<float> >*)flann_index_;
		delete temp_ptr;
		flann_index_ = NULL;
	}
}

KnnSearch::KnnSearch(int npts, float* pts_data, int dim) : flann_index_(NULL), dim_(dim)
{
	rebuild(npts, pts_data);
}

KnnSearch::KnnSearch() : flann_index_(NULL), dim_(3)
{

}