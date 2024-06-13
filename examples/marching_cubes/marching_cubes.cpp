#include "marching_cubes.h"

#include <fstream>
#include <algorithm>
#include <ctime>
#include <string.h>

float MarchingCubes::findNonZeroMinEdgeLength(float* p, unsigned int numTri, unsigned int* t)
{
	float edgeMin = std::numeric_limits<float>::infinity();

	for (unsigned int i = 0; i < numTri; ++i)
	{
		mg::Vector3f v[3];
		v[0] = mg::Vector3f(p+3*t[3*i+0]);
		v[1] = mg::Vector3f(p+3*t[3*i+1]);
		v[2] = mg::Vector3f(p+3*t[3*i+2]);

		//if (i == 0)
		//{
		//	std::cout << t[3 * i + 0] << " : " << v[0] << "\n"
		//		<< t[3 * i + 1] << " : " << v[1] << "\n"
		//		<< t[3 * i + 2] << " : " << v[2] << std::endl;
		//}

		float d2 = v[0].dist2To(v[1]);
		if (d2 > 0.0 && d2 < edgeMin)	edgeMin = d2;

		d2 = v[1].dist2To(v[2]);
		if (d2 > 0.0 && d2 < edgeMin)	edgeMin = d2;

		d2 = v[2].dist2To(v[0]);
		if (d2 > 0.0 && d2 < edgeMin)	edgeMin = d2;
	}
	std::cout << "- edge min = " << sqrt(edgeMin) << std::endl;

	return sqrt(edgeMin);
}

void MarchingCubes::computeBoundingBox(unsigned int numVert, mg::Vector3f* p, float* bmin, float* bmax)
{
	float INF = std::numeric_limits<float>::infinity();
	bmin[0] = bmin[1] = bmin[2] = INF;
	bmax[0] = bmax[1] = bmax[2] = -INF;

	for (unsigned int i = 0; i < numVert; ++i)
	{
		if (p[i].x < bmin[0])	bmin[0] = p[i].x;
		if (p[i].y < bmin[1])	bmin[1] = p[i].y;
		if (p[i].z < bmin[2])	bmin[2] = p[i].z;

		if (p[i].x > bmax[0])	bmax[0] = p[i].x;
		if (p[i].y > bmax[1])	bmax[1] = p[i].y;
		if (p[i].z > bmax[2])	bmax[2] = p[i].z;
	}
	std::cout << "- Bounding box : " << mg::Vector3f(bmin) << " ~ "	<< mg::Vector3f(bmax) << std::endl;
}

int MarchingCubes::findMergeIndex_gridSearch(std::vector<int>& mergeInd, SearchGrid3D& gridSearch, 
	unsigned int numVert, mg::Vector3f* p,
	float r, std::vector< std::vector<int> >& indices, std::vector< std::vector<float> >&dists)
{
	// Find neighbors of vertices
	indices.resize(numVert);
	dists.resize(numVert);

	//std::cout << "r=" << r << std::endl;

	float r2 = r * r;

	for (unsigned int i = 0; i < numVert; ++i)
	{
		gridSearch.find_neighbors(p[i].data(), r, indices[i]);

		std::vector<int> refinedIndices(indices[i].size());

		dists[i].resize(indices[i].size());

		int count = 0;
		for (int j = 0; j < indices[i].size(); ++j)
		{
			float d2 = p[i].dist2To(p[indices[i][j]]);

			if (d2 < r2)
			{
				refinedIndices[count] = indices[i][j];
				dists[i][count] = d2;

				count++;
			}
		}

		refinedIndices.resize(count);
		dists[i].resize(count);

		std::swap(indices[i], refinedIndices);

		//if (i == 0 || i == 1 || i == 2)
		//{
		//	std::cout << i << " : ";
		//	for (int j = 0; j < indices[i].size(); ++j)
		//	{
		//		std::cout << indices[i][j] << "(" << dists[i][j] << "), ";
		//	}
		//	std::cout << std::endl;
		//}
	}

	int currInd = 0;

	// For the first vertex (id=0)
	mergeInd[0] = currInd;
	for (int j = 0; j < indices[0].size(); ++j)
	{
		if (dists[0][j] < r2)
			mergeInd[indices[0][j]] = currInd;
		else
			break;
	}
	currInd++;

	// For the rest vertices(id>0)
	//std::vector<int> subNewInd(max_nn);

	for (int i = 1; i < indices.size(); ++i)
	{
		if (mergeInd[i] != -1) continue;

		mergeInd[i] = currInd;

		//for (int j = 1; j < max_nn; ++j)
		//	subNewInd[j] = mergeInd[indices[i*max_nn + j]];

		for (int j = 0; j < indices[i].size(); ++j)
		{
			if (dists[i][j] < r2)
				mergeInd[indices[i][j]] = currInd;
			else
				break;
		}
		currInd++;
	}
	//std::cout << std::endl;

	std::cout << "- Number of vertices : " << numVert << "->" << currInd << std::endl;

	return currInd;
}

void MarchingCubes::mergeVertices(unsigned int& numVert, float* vertex, unsigned int& numTri, unsigned int* triangle)
{
	float minNonZeroEdge = findNonZeroMinEdgeLength(vertex, numTri, triangle);
	float bmin[3], bmax[3];
	computeBoundingBox(numVert, (mg::Vector3f*)vertex, bmin, bmax);

	const int tempGridSize = 128;
	const float hx = (bmax[0] - bmin[0]) / tempGridSize;
	const float hy = (bmax[1] - bmin[1]) / tempGridSize;
	const float hz = (bmax[2] - bmin[2]) / tempGridSize;

	float h = std::max(std::max(std::max(hx, hy), hz), minNonZeroEdge);// edgeTol;

	SearchGrid3D gridSearch(bmin, bmax, h);

	gridSearch.setPoints(numVert, vertex);
	gridSearch.build();

	std::vector<int> mergeInd(numVert);
	std::fill(mergeInd.begin(), mergeInd.end(), -1);

	int max_nn = 25;
	std::vector< std::vector<int> > indices;
	std::vector< std::vector<float> > dists;
	int numNewVert = findMergeIndex_gridSearch(mergeInd, gridSearch, numVert, (mg::Vector3f*)vertex, 
		0.45f*minNonZeroEdge, indices, dists);

	//auto result2 = std::minmax_element(mergeInd.begin(), mergeInd.end());
	//std::cout << *result2.first << " " << *result2.second << std::endl;
	//std::cout << mergeInd[0] << " " << mergeInd[1] << " " << mergeInd[2] << std::endl;

	//mg::TriMesh_f mergedMesh;
	//mesh.deepCopyTo(mergedMesh);
	std::vector<mg::Vector3f> mergeVertex(numNewVert);
	std::vector<mg::Vector3i> mergeTriangle(numTri);

	//mergedMesh.vertex_.resize(numNewVert);
	mg::Vector3f * p = (mg::Vector3f*)vertex;
	for (unsigned int i = 0; i < numVert; ++i)
	{
		mergeVertex[mergeInd[i]] = p[i];
	}

	int count = 0;
	for (size_t i = 0; i < numTri; ++i)
	{
		int ti, tj, tk;
		ti = mergeInd[triangle[3*i+0]];
		tj = mergeInd[triangle[3*i+1]];
		tk = mergeInd[triangle[3*i+2]];

		if (ti != tj && tj != tk && tk != ti)
		{
			mergeTriangle[count]
				= mg::Vector3i(ti, tj, tk);
			count++;
		}
	}

	mergeTriangle.resize(count);

	// Update vertex and triangle information
	numVert = numNewVert;
	memcpy(vertex, mergeVertex.data(), sizeof(mg::Vector3f)*numVert);

	numTri = count;
	for (unsigned int i = 0; i < numTri; ++i)
	{
		triangle[3 * i + 0] = mergeTriangle[i].x;
		triangle[3 * i + 1] = mergeTriangle[i].y;
		triangle[3 * i + 2] = mergeTriangle[i].z;
	}
}

bool MarchingCubes::exportMeshInObj(const char* filename, unsigned int numVertex, float* vertex, unsigned int numTri, unsigned int* triangle)
{
	float stime = clock();

	std::ofstream outFile(filename);

	if (!outFile.is_open())
	{
		std::cerr << "ERROR! Can't open " << filename << " for writing" << std::endl;
		return false;
	}

	//unsigned int numVertex = (unsigned int)vertex_.size();

	for (unsigned int i = 0; i < numVertex; ++i)
	{
		outFile << "v " << vertex[3*i+0] << " " << vertex[3*i+1] << " " << vertex[3*i+2] << "\n";
	}

	//unsigned int numTri = (unsigned int)triangle_.size();

	if (!genTexture_)
	{
		for (unsigned int i = 0; i < numTri; ++i)
		{
			outFile << "f "
				<< triangle[3 * i + 0] + 1 << " "
				<< triangle[3 * i + 1] + 1 << " "
				<< triangle[3 * i + 2] + 1 << " \n";
		}

		outFile.close();

		float ftime = clock();

		std::cout << "- Mesh exported : " << filename << " (#v=" << numVertex
			<< ", #t=" << numTri << ")" << "..." << (ftime - stime) / CLOCKS_PER_SEC << " sec" << std::endl;
		//std::cout << "- Mesh exported : " << filename << " (#v=" << numVertex
		//    << ",#t=" << numTri << ")" << std::endl;
	}
	else
	{
		outFile << "mtllib test.mtl" << "\n"
			<< "g default" << "\n";

		float Lx = img_dx_ * img_nx_;
		float Ly = img_dy_ * img_ny_;

		for (unsigned int i = 0; i < numVertex; ++i)
		{
			float u, v;
			u = (vertex[3 * i + 0] - img_xmin_) / Lx;
			v = (vertex[3 * i + 1] - img_ymin_) / Ly;
			if (u < 0.0) u = 0.0f; else if (u > 1.0f) u = 1.0f;
			if (v < 0.0) v = 0.0f; else if (v > 1.0f) v = 1.0f;
			outFile << "vt " << u << " " << v << "\n";
		}

		for (unsigned int i = 0; i < numTri; ++i)
		{
			outFile << "f "
				<< triangle[3 * i + 0] + 1 << "/" << triangle[3 * i + 0] + 1 << "/ "
				<< triangle[3 * i + 1] + 1 << "/" << triangle[3 * i + 1] + 1 << "/ "
				<< triangle[3 * i + 2] + 1 << "/" << triangle[3 * i + 2] + 1 << "/ \n";
		}

		outFile.close();

		float ftime = clock();

		std::cout << "- Mesh exported : " << filename << " (#v=" << numVertex
			<< ", #vt=" << numVertex << ", #t=" << numTri << ")" << "..." << (ftime - stime) / CLOCKS_PER_SEC << " sec" << std::endl;

		// mtl
		std::ofstream outFile2("test.mtl");
		
		outFile2 << "newmtl initialShadingGroup" << "\n"
			<< "Ns 154.901961" << "\n"
			<< "Ka 0.000000 0.000000 0.000000" << "\n"
			<< "Kd 1.000000 1.000000 1.00000" << "\n"
			<< "Ks 0.0100000 0.0100000 0.01000" << "\n"
			<< "Ni 1.000000" << "\n"
			<< "d 1.000000" << "\n"
			<< "illum 2" << "\n"
			<< "map_Kd " << "output_texture.png" << "\n";

		outFile2.close();
	}

	return true;
}

bool MarchingCubes::exportMeshInObj(const char* filename)
{
	float stime = clock();

	std::ofstream outFile(filename);

	if (!outFile.is_open())
	{
		std::cerr << "ERROR! Can't open " << filename << " for writing" << std::endl;
		return false;
	}

	unsigned int numVertex = (unsigned int)vertex_.size();
	for (unsigned int i = 0; i < numVertex; ++i)
	{
		outFile << "v " << vertex_[i].x << " " << vertex_[i].y << " " << vertex_[i].z << "\n";
	}

	unsigned int numTri = (unsigned int)triangle_.size();

	for (unsigned int i = 0; i < numTri; ++i)
	{
		outFile << "f " 
			<< triangle_[i].x + 1 << " "
			<< triangle_[i].y + 1 << " "
			<< triangle_[i].z + 1 << " \n";
	}

	outFile.close();

	float ftime = clock();

	std::cout << "- Mesh exported : " << filename << "..." << (ftime - stime) / CLOCKS_PER_SEC << " sec" << std::endl;
	//std::cout << "- Mesh exported : " << filename << " (#v=" << numVertex
	//    << ",#t=" << numTri << ")" << std::endl;

	return true;
}

unsigned char MarchingCubes::findCubeIndex(unsigned int i, unsigned int j, unsigned int k, float isoValue, float* gridValue)
{
	mg::Vector3ui gridSize = volumeSize_ + mg::Vector3ui(1, 1, 1);

	size_t cube_vid[8];
	cube_vid[0] = (size_t)k*gridSize.y*gridSize.x + j * gridSize.x + i;
	cube_vid[1] = cube_vid[0] + 1;
	cube_vid[2] = cube_vid[0] + 1 + gridSize.x;
	cube_vid[3] = cube_vid[0] +     gridSize.x;
	cube_vid[4] = cube_vid[0] + gridSize.x*gridSize.y;
	cube_vid[5] = cube_vid[1] + gridSize.x*gridSize.y;
	cube_vid[6] = cube_vid[2] + gridSize.x*gridSize.y;
	cube_vid[7] = cube_vid[3] + gridSize.x*gridSize.y;

	gridValue[0] = f_[cube_vid[0]];
	gridValue[1] = f_[cube_vid[1]];
	gridValue[2] = f_[cube_vid[2]];
	gridValue[3] = f_[cube_vid[3]];
	gridValue[4] = f_[cube_vid[4]];
	gridValue[5] = f_[cube_vid[5]];
	gridValue[6] = f_[cube_vid[6]];
	gridValue[7] = f_[cube_vid[7]];

	unsigned char cubeindex = 0;
	if (f_[cube_vid[0]] < isoValue) cubeindex |= 1;
	if (f_[cube_vid[1]] < isoValue) cubeindex |= 2;
	if (f_[cube_vid[2]] < isoValue) cubeindex |= 4;
	if (f_[cube_vid[3]] < isoValue) cubeindex |= 8;
	if (f_[cube_vid[4]] < isoValue) cubeindex |= 16;
	if (f_[cube_vid[5]] < isoValue) cubeindex |= 32;
	if (f_[cube_vid[6]] < isoValue) cubeindex |= 64;
	if (f_[cube_vid[7]] < isoValue) cubeindex |= 128;

	return cubeindex;
}

mg::Vector3f MarchingCubes::vertexInterp( float isoValue, mg::Vector3f p1, mg::Vector3f p2, float valp1, float valp2)
{
	float mu;
	mg::Vector3f p;

	if (fabs(isoValue - valp1) < eps) //0.00001)
		return p1;
	if (fabs(isoValue - valp2) < eps) //0.00001)
		return p2;
	if (fabs(valp1 - valp2) < eps) //0.00001)
		return p1;

	mu = (isoValue - valp1) / (valp2 - valp1);
	
	p.x = p1.x + mu * (p2.x - p1.x);
	p.y = p1.y + mu * (p2.y - p1.y);
	p.z = p1.z + mu * (p2.z - p1.z);

	return p;
}

void MarchingCubes::generateSurfaceMesh(float isoValue, float CUT_VAL)
{
	float stime = clock();
	std::cout << "- Extracting mesh...";

	const unsigned int nx = volumeSize_.x;
	const unsigned int ny = volumeSize_.y;
	const unsigned int nz = volumeSize_.z;

	const float dx = voxelSize_.x;
	const float dy = voxelSize_.y;
	const float dz = voxelSize_.z;

	size_t numVoxels = (unsigned int)nx * ny * nz;

	std::vector<mg::Vector3f>	vertlist(12);
	std::vector<float>			gv(8);
	std::vector<mg::Vector3f>   gp(8);

	int triIndex = 0;
	int vertexIndex = 0;
	vertex_.resize(0);
	triangle_.resize(0);

	for (size_t voxId = 0; voxId < numVoxels; ++voxId)
	{
		unsigned int i = voxId % nx;
		unsigned int j = ((voxId - i) / nx) % ny;
		unsigned int k = ((voxId - i) / nx) / ny;

		// Find the cube index for the current voxel=
		unsigned char cubeindex = findCubeIndex(i, j, k, isoValue, gv.data());

		if (fabs(gv[0]) >= CUT_VAL || fabs(gv[1]) >= CUT_VAL || fabs(gv[2]) >= CUT_VAL ||
			fabs(gv[3]) >= CUT_VAL || fabs(gv[4]) >= CUT_VAL || fabs(gv[5]) >= CUT_VAL ||
			fabs(gv[6]) >= CUT_VAL || fabs(gv[7]) >= CUT_VAL 
			|| edgeTable_[cubeindex] == 0) continue;

		// Now, the current is a surface cell
		//std::fill(vertlist.begin(), vertlist.end(), mg::Vector3f(0.0, 0.0, 0.0));

		gp[0] = bmin_ + mg::Vector3f(      i*dx,       j*dy, k*dz);
		gp[1] = bmin_ + mg::Vector3f((i + 1)*dx,       j*dy, k*dz);
		gp[2] = bmin_ + mg::Vector3f((i + 1)*dx, (j + 1)*dy, k*dz);
		gp[3] = bmin_ + mg::Vector3f(      i*dx, (j + 1)*dy, k*dz);
		gp[4] = gp[0] + mg::Vector3f(0.0f, 0.0f, dz);
		gp[5] = gp[1] + mg::Vector3f(0.0f, 0.0f, dz);
		gp[6] = gp[2] + mg::Vector3f(0.0f, 0.0f, dz);
		gp[7] = gp[3] + mg::Vector3f(0.0f, 0.0f, dz);

		// Create vertices on edges
		if (edgeTable_[cubeindex] & 1)
			vertlist[0] = vertexInterp(isoValue, gp[0], gp[1], gv[0], gv[1]);
		if (edgeTable_[cubeindex] & 2)
			vertlist[1] = vertexInterp(isoValue, gp[1], gp[2], gv[1], gv[2]);
		if (edgeTable_[cubeindex] & 4)
			vertlist[2] = vertexInterp(isoValue, gp[2], gp[3], gv[2], gv[3]);
		if (edgeTable_[cubeindex] & 8)
			vertlist[3] = vertexInterp(isoValue, gp[3], gp[0], gv[3], gv[0]);
		if (edgeTable_[cubeindex] & 16)
			vertlist[4] = vertexInterp(isoValue, gp[4], gp[5], gv[4], gv[5]);
		if (edgeTable_[cubeindex] & 32)
			vertlist[5] = vertexInterp(isoValue, gp[5], gp[6], gv[5], gv[6]);
		if (edgeTable_[cubeindex] & 64)
			vertlist[6] = vertexInterp(isoValue, gp[6], gp[7], gv[6], gv[7]);
		if (edgeTable_[cubeindex] & 128)
			vertlist[7] = vertexInterp(isoValue, gp[7], gp[4], gv[7], gv[4]);
		if (edgeTable_[cubeindex] & 256)
			vertlist[8] = vertexInterp(isoValue, gp[0], gp[4], gv[0], gv[4]);
		if (edgeTable_[cubeindex] & 512)
			vertlist[9] = vertexInterp(isoValue, gp[1], gp[5], gv[1], gv[5]);
		if (edgeTable_[cubeindex] & 1024)
			vertlist[10] = vertexInterp(isoValue, gp[2], gp[6], gv[2], gv[6]);
		if (edgeTable_[cubeindex] & 2048)
			vertlist[11] = vertexInterp(isoValue, gp[3], gp[7], gv[3], gv[7]);

		// Generate triangles
		for (int i = 0; triTable_[cubeindex][i] != -1; i += 3) {
			vertex_.resize(vertexIndex + 3);
			vertex_[vertexIndex + 0] = vertlist[triTable_[cubeindex][i]];
			vertex_[vertexIndex + 1] = vertlist[triTable_[cubeindex][i + 1]];
			vertex_[vertexIndex + 2] = vertlist[triTable_[cubeindex][i + 2]];

			triangle_.resize(triIndex + 1);
			triangle_[triIndex].x = vertexIndex;
			triangle_[triIndex].y = vertexIndex + 1;
			triangle_[triIndex].z = vertexIndex + 2;

			vertexIndex += 3;
			triIndex	+= 1;
		}
	}
	float ftime = clock();
	std::cout << "done..." << (ftime - stime) / CLOCKS_PER_SEC << " sec" << std::endl;
}

bool MarchingCubes::loadVolumeData(const char* filename)
{
	std::cout << "- Reading volume data : " << filename << "..." << std::flush;

	std::ifstream inFile(filename, std::ios::binary);
	if (!inFile.is_open())
	{
		std::cout << "Error, can't load the volume data : " << filename << std::endl;
		return false;
	}

	mg::Vector3f bmin, bmax;
	inFile.read((char*)bmin.data(), sizeof(float) * 3);
	inFile.read((char*)bmax.data(), sizeof(float) * 3);

	center_ = 0.5f*(bmin + bmax);
	bmin_ = bmin;

	unsigned int gridNx, gridNy, gridNz;
	inFile.read((char*)&gridNx, sizeof(unsigned int));
	inFile.read((char*)&gridNy, sizeof(unsigned int));
	inFile.read((char*)&gridNz, sizeof(unsigned int));

	volumeSize_ = mg::Vector3ui(gridNx-1, gridNy-1, gridNz-1);

	voxelSize_.x = (bmax.x - bmin.x) / volumeSize_.x;
	voxelSize_.y = (bmax.y - bmin.y) / volumeSize_.y;
	voxelSize_.z = (bmax.z - bmin.z) / volumeSize_.z;

	size_t n = (size_t)gridNx*gridNy*gridNz;

	f_.resize(n);

	inFile.read((char*)f_.data(), sizeof(float)*n);

	inFile.close();

	std::cout << "done" << std::endl;

	std::cout << " . Bounding box : " << bmin << "~" << bmax << std::endl;

	std::cout << " . Resolution = " << volumeSize_ << std::endl;

	// Find min-max
	auto result = std::minmax_element(f_.begin(), f_.end());
	std::cout << " . Data range : " << *(result.first) << "-" << *(result.second) << std::endl;

	return true;
}

MarchingCubes::MarchingCubes() : center_(), bmin_(), volumeSize_(), voxelSize_(), isoValue_(0.0f), genTexture_(false)
{}