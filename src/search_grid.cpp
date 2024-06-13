#include "search_grid.h"
#include <string.h>

#include <iostream>

void SearchGrid2D::find_neighbors(float* query_point, float r, std::vector<int>& neighbor_ind)
{
    neighbor_ind.clear();

    // Find the cell containing the query point
    int ci = (int)((query_point[0] - bmin_.x) / h_);
    int cj = (int)((query_point[1] - bmin_.y) / h_);
    //std::cout << "(ci, cj)=" << ci << ", " << cj << std::endl;

    // Find the cell radius
    int b = (int)(r / h_) + 1;
    //std::cout << "b=" << std::endl;

    int j0 = cj - b; if (j0 < 0)  j0 = 0;
    int j1 = cj + b; if (j1 >= ny_)  j1 = ny_-1;

    int i0 = ci - b; if (i0 < 0)     i0 = 0;
    int i1 = ci + b; if (i1 >= nx_)  i1 = nx_ - 1;
                         
    for (int j = j0; j <= j1; ++j)
    {
        for (int i = i0; i <= i1; ++i)
        {
            int cindex = j*nx_ + i;

            neighbor_ind.insert(neighbor_ind.end(), pts_of_cell_[cindex].begin(), pts_of_cell_[cindex].end());
        }
    }

    //return cindex;
}

void SearchGrid2D::build()
{
    const int numCells = nx_*ny_;

    // Initialize
    pts_of_cell_.resize(numCells);
    for (int i = 0; i < numCells; ++i)
    {
        pts_of_cell_[i].clear();
    }

    // Build search structure
    const int numPts = (int)p_.size();
    for (int i = 0; i < numPts; ++i)
    {
        int pi = (int)((p_[i].x - bmin_.x) / h_);
        int pj = (int)((p_[i].y - bmin_.y) / h_);

        int cindex = pj*nx_ + pi;

        pts_of_cell_[cindex].push_back(i);
    }

    // for verification
    /*for (int i = 0; i < numCells; ++i)
    {
        std::cout << " cell[" << i << "] : ";
        std::vector<int>& nid = pts_of_cell_[i];
        for (int j = 0; j < nid.size(); ++j)
        {
            std::cout << nid[j] << " ";
        }
        std::cout << std::endl;
    }*/
}

void SearchGrid2D::setPoints(int n, float* p)
{
    p_.resize(n);

    memcpy(p_.data()->data(), p, sizeof(mg::Vector2f)*n);
}

void SearchGrid2D::setGrid(float* bmin, float* bmax, float h)
{
    h_ = h;

    nx_ = (int)floor((bmax[0] - bmin[0]) / h) + 1;
    ny_ = (int)floor((bmax[1] - bmin[1]) / h) + 1;

    bmin_ = mg::Vector2f(bmin);

    bmax_.x = bmin[0] + h*nx_;
    bmax_.y = bmin[1] + h*ny_;

    std::cout << "- Search grid : " << nx_ << ", " << ny_ << std::endl;
}

SearchGrid2D::SearchGrid2D(float* bmin, float* bmax, float h)
{
    setGrid(bmin, bmax, h);
}


void SearchGrid3D::find_neighbors(float* query_point, float r, std::vector<int>& neighbor_ind)
{
	neighbor_ind.clear();

	// Find the cell containing the query point
	int ci = (int)((query_point[0] - bmin_.x) / h_);
	int cj = (int)((query_point[1] - bmin_.y) / h_);
	int ck = (int)((query_point[2] - bmin_.z) / h_);
	//std::cout << "(ci, cj, ck)=" << ci << ", " << cj << ", " << ck << std::endl;

	// Find the cell radius
	int b = (int)(r / h_) + 1;
	//std::cout << "b=" << std::endl;

	int k0 = ck - b; if (k0 < 0)  k0 = 0;
	int k1 = ck + b; if (k1 >= nz_)  k1 = nz_ - 1;

	int j0 = cj - b; if (j0 < 0)  j0 = 0;
	int j1 = cj + b; if (j1 >= ny_)  j1 = ny_ - 1;

	int i0 = ci - b; if (i0 < 0)     i0 = 0;
	int i1 = ci + b; if (i1 >= nx_)  i1 = nx_ - 1;

	for (int k = k0; k <= k1; ++k)
	{
		for (int j = j0; j <= j1; ++j)
		{
			for (int i = i0; i <= i1; ++i)
			{
				int cindex = k * nx_ * ny_ + j * nx_ + i;

				neighbor_ind.insert(neighbor_ind.end(), pts_of_cell_[cindex].begin(), pts_of_cell_[cindex].end());
			}
		}
	}

	//return cindex;
}

void SearchGrid3D::build()
{
	const int numCells = nx_ * ny_ * nz_;

	// Initialize
	pts_of_cell_.resize(numCells);
	for (int i = 0; i < numCells; ++i)
	{
		pts_of_cell_[i].clear();
	}

	// Build search structure
	const int numPts = (int)p_.size();
	for (int i = 0; i < numPts; ++i)
	{
		int pi = (int)((p_[i].x - bmin_.x) / h_);
		int pj = (int)((p_[i].y - bmin_.y) / h_);
		int pk = (int)((p_[i].z - bmin_.z) / h_);

		int cindex = pk * nx_ * ny_ + pj * nx_ + pi;

		pts_of_cell_[cindex].push_back(i);
	}

	// for verification
	/*for (int i = 0; i < numCells; ++i)
	{
		std::cout << " cell[" << i << "] : ";
		std::vector<int>& nid = pts_of_cell_[i];
		for (int j = 0; j < nid.size(); ++j)
		{
			std::cout << nid[j] << " ";
		}
		std::cout << std::endl;
	}*/
}

void SearchGrid3D::setPoints(int n, float* p)
{
	p_.resize(n);

	memcpy(p_.data()->data(), p, sizeof(mg::Vector3f)*n);
}

void SearchGrid3D::setGrid(float* bmin, float* bmax, float h)
{
	h_ = h;

	nx_ = (int)floor((bmax[0] - bmin[0]) / h) + 1;
	ny_ = (int)floor((bmax[1] - bmin[1]) / h) + 1;
	nz_ = (int)floor((bmax[2] - bmin[2]) / h) + 1;

	bmin_ = mg::Vector3f(bmin);

	bmax_.x = bmin[0] + h * nx_;
	bmax_.y = bmin[1] + h * ny_;
	bmax_.z = bmin[2] + h * nz_;

	std::cout << "- Search grid : (" << nx_ << ", " << ny_ << ", " << nz_ << ")" <<  std::endl;
}

SearchGrid3D::SearchGrid3D(float* bmin, float* bmax, float h)
{
	setGrid(bmin, bmax, h);
}