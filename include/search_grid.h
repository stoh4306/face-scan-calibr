#ifndef SEARCH_GRID_CLASS_H_
#define SEARCH_GRID_CLASS_H_

//#include <Eigen/Dense>
#include "vec.h"
#include <vector>

//! 2D search grid class 
class SearchGrid2D
{
public:
    //! Default constructor
    SearchGrid2D() {}

    /*! Constructor with a 2d bounding box and a cell size
    * \param bmin   min vertex of the 2d bounding box
    * \param bmax   max vertex of the 2d bounding box
    * \param h      the length of a single cell
    */
    SearchGrid2D(float* bmin, float* bmax, float h);

    /*! Set search grid with a 2d bounding box and a cell size
    * \param bmin   min vertex of the 2d bounding box
    * \param bmax   max vertex of the 2d bounding box
    * \param h      the length of a single cell
    */
    void setGrid(float* bmin, float* bmax, float h);

    /*! Set points to search
    * \param n  number of points
    * \param p  pointer of 2d point array
    */
    void setPoints(int n, float* p);

    //! Build the search structure
    void build();

    /*! Find the neighbor points for a query point and a radius
    * \param query_point    the query point 
    * \param r              search radius
    * \param neighbor_ind   output array of indices of neighbors of the query point
    */
    void find_neighbors(float* query_point, float r, std::vector<int>& neighbor_ind);

public:
    mg::Vector2f bmin_, bmax_;
    float h_;
    int nx_, ny_;

    std::vector< mg::Vector2f > p_;

    std::vector< std::vector<int> > pts_of_cell_;
};


//! 3D search grid class 
class SearchGrid3D
{
public:
	//! Default constructor
	SearchGrid3D() {}

	/*! Constructor with a 3d bounding box and a cell size
	* \param bmin   min vertex of the 3d bounding box
	* \param bmax   max vertex of the 3d bounding box
	* \param h      the length of a single cell
	*/
	SearchGrid3D(float* bmin, float* bmax, float h);

	/*! Set search grid with a 3d bounding box and a cell size
	* \param bmin   min vertex of the 3d bounding box
	* \param bmax   max vertex of the 3d bounding box
	* \param h      the length of a single cell
	*/
	void setGrid(float* bmin, float* bmax, float h);

	/*! Set points to search
	* \param n  number of points
	* \param p  pointer of 3d point array
	*/
	void setPoints(int n, float* p);

	//! Build the search structure
	void build();

	/*! Find the neighbor points for a query point and a radius
	* \param query_point    the query point
	* \param r              search radius
	* \param neighbor_ind   output array of indices of neighbors of the query point
	*/
	void find_neighbors(float* query_point, float r, std::vector<int>& neighbor_ind);

public:
	mg::Vector3f bmin_, bmax_;
	float h_;
	int nx_, ny_, nz_;

	std::vector< mg::Vector3f > p_;

	std::vector< std::vector<int> > pts_of_cell_;
};


#endif