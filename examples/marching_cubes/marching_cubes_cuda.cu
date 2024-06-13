#include "marching_cubes.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>

__device__ unsigned char findCubeIndex(uint3 gridSize, float* f, unsigned int i, unsigned int j, unsigned int k, float isoValue, float* gridValue)
{
	unsigned int cube_vid[8];
	cube_vid[0] = (unsigned int)k*gridSize.y*gridSize.x + j * gridSize.x + i;
	cube_vid[1] = cube_vid[0] + 1;
	cube_vid[2] = cube_vid[0] + 1 + gridSize.x;
	cube_vid[3] = cube_vid[0] + gridSize.x;
	cube_vid[4] = cube_vid[0] + gridSize.x*gridSize.y;
	cube_vid[5] = cube_vid[1] + gridSize.x*gridSize.y;
	cube_vid[6] = cube_vid[2] + gridSize.x*gridSize.y;
	cube_vid[7] = cube_vid[3] + gridSize.x*gridSize.y;

	gridValue[0] = f[cube_vid[0]];
	gridValue[1] = f[cube_vid[1]];
	gridValue[2] = f[cube_vid[2]];
	gridValue[3] = f[cube_vid[3]];
	gridValue[4] = f[cube_vid[4]];
	gridValue[5] = f[cube_vid[5]];
	gridValue[6] = f[cube_vid[6]];
	gridValue[7] = f[cube_vid[7]];

	unsigned char cubeindex = 0;
	if (f[cube_vid[0]] < isoValue) cubeindex |= 1;
	if (f[cube_vid[1]] < isoValue) cubeindex |= 2;
	if (f[cube_vid[2]] < isoValue) cubeindex |= 4;
	if (f[cube_vid[3]] < isoValue) cubeindex |= 8;
	if (f[cube_vid[4]] < isoValue) cubeindex |= 16;
	if (f[cube_vid[5]] < isoValue) cubeindex |= 32;
	if (f[cube_vid[6]] < isoValue) cubeindex |= 64;
	if (f[cube_vid[7]] < isoValue) cubeindex |= 128;

	return cubeindex;
}

__device__ float3 vertexInterp(float isoValue, float3 p1, float3 p2, float valp1, float valp2, float eps)
{
	float mu;
	float3 p;

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

__global__ void kernel_extract_mesh(unsigned int nx, unsigned int ny, unsigned int nz, unsigned int numVoxels,
	float bmin_x, float bmin_y, float bmin_z, float dx, float dy, float dz,
	float* f, int* edgeTable, int* triTable,
	unsigned int* vertIndex, unsigned int* triIndex, float* vertex, unsigned int* triangle,
	float isoValue, float eps, float CUT_VAL)
{
	unsigned int tid = threadIdx.x + (unsigned int)blockIdx.x*blockDim.x;
	if (tid < numVoxels)
	{
		float3 vertlist[12];
		float gv[8];
		float3 gp[8];

		unsigned int i = tid % nx;
		unsigned int j = ((tid - i) / nx) % ny;
		unsigned int k = ((tid - i) / nx) / ny;

		uint3 gridSize = make_uint3(nx + 1, ny + 1, nz + 1);

		// Find the cube index for the current voxel
		unsigned char cubeindex = findCubeIndex(gridSize, f, i, j, k, isoValue, gv);

		if (fabs(gv[0]) < CUT_VAL && fabs(gv[1]) < CUT_VAL && fabs(gv[2]) && CUT_VAL &&
			fabs(gv[3]) < CUT_VAL && fabs(gv[4]) < CUT_VAL && fabs(gv[5]) && CUT_VAL &&
			fabs(gv[6]) < CUT_VAL && fabs(gv[7]) < CUT_VAL &&
			edgeTable[cubeindex] > 0)
		{
			// Now, the current is a surface cell
			//std::fill(vertlist.begin(), vertlist.end(), mg::Vector3f(0.0, 0.0, 0.0));
			float3 bmin = make_float3(bmin_x, bmin_y, bmin_z);

			gp[0].x = bmin.x + i * dx;
			gp[0].y = bmin.y + j * dy;
			gp[0].z = bmin.z + k * dz;	// make_float3(i*dx, j*dy, k*dz);

			gp[1].x = bmin.x + (i + 1) * dx;
			gp[1].y = bmin.y + j * dy;
			gp[1].z = bmin.z + k * dz;  //gp[1] = bmin + make_float3((i + 1)*dx, j*dy, k*dz);

			gp[2].x = bmin.x + (i + 1) * dx;
			gp[2].y = bmin.y + (j + 1) * dy;
			gp[2].z = bmin.z + k * dz;	//gp[2] = bmin + make_float3((i + 1)*dx, (j + 1)*dy, k*dz);

			gp[3].x = bmin.x + i * dx;
			gp[3].y = bmin.y + (j + 1) * dy;
			gp[3].z = bmin.z + k * dz;	//gp[3] = bmin + make_float3(i*dx, (j + 1)*dy, k*dz);

			gp[4].x = gp[0].x;
			gp[4].y = gp[0].y;
			gp[4].z = gp[0].z + dz;        //gp[4] = gp[0] + make_float3(0.0f, 0.0f, dz);

			gp[5].x = gp[1].x;
			gp[5].y = gp[1].y;
			gp[5].z = gp[1].z + dz;			//gp[5] = gp[1] + make_float3(0.0f, 0.0f, dz);

			gp[6].x = gp[2].x;
			gp[6].y = gp[2].y;
			gp[6].z = gp[2].z + dz;			//gp[6] = gp[2] + make_float3(0.0f, 0.0f, dz);

			gp[7].x = gp[3].x;
			gp[7].y = gp[3].y;
			gp[7].z = gp[3].z + dz;		 //gp[7] = gp[3] + make_float3(0.0f, 0.0f, dz);

			// Create vertices on edges
			if (edgeTable[cubeindex] & 1)		vertlist[0] = vertexInterp(isoValue, gp[0], gp[1], gv[0], gv[1], eps);
			if (edgeTable[cubeindex] & 2)		vertlist[1] = vertexInterp(isoValue, gp[1], gp[2], gv[1], gv[2], eps);
			if (edgeTable[cubeindex] & 4)		vertlist[2] = vertexInterp(isoValue, gp[2], gp[3], gv[2], gv[3], eps);
			if (edgeTable[cubeindex] & 8)		vertlist[3] = vertexInterp(isoValue, gp[3], gp[0], gv[3], gv[0], eps);
			if (edgeTable[cubeindex] & 16)		vertlist[4] = vertexInterp(isoValue, gp[4], gp[5], gv[4], gv[5], eps);
			if (edgeTable[cubeindex] & 32)		vertlist[5] = vertexInterp(isoValue, gp[5], gp[6], gv[5], gv[6], eps);
			if (edgeTable[cubeindex] & 64)		vertlist[6] = vertexInterp(isoValue, gp[6], gp[7], gv[6], gv[7], eps);
			if (edgeTable[cubeindex] & 128)		vertlist[7] = vertexInterp(isoValue, gp[7], gp[4], gv[7], gv[4], eps);
			if (edgeTable[cubeindex] & 256)		vertlist[8] = vertexInterp(isoValue, gp[0], gp[4], gv[0], gv[4], eps);
			if (edgeTable[cubeindex] & 512)		vertlist[9] = vertexInterp(isoValue, gp[1], gp[5], gv[1], gv[5], eps);
			if (edgeTable[cubeindex] & 1024)	vertlist[10] = vertexInterp(isoValue, gp[2], gp[6], gv[2], gv[6], eps);
			if (edgeTable[cubeindex] & 2048)	vertlist[11] = vertexInterp(isoValue, gp[3], gp[7], gv[3], gv[7], eps);

			// Generate triangles
			for (int ti = 0; triTable[cubeindex * 16 + ti] != -1; ti += 3) {
				//vertex_.resize(vertexIndex + 3);
				int vid = atomicAdd(vertIndex, 3);
				vertex[3 * (vid + 0) + 0] = vertlist[triTable[cubeindex * 16 + ti]].x;
				vertex[3 * (vid + 0) + 1] = vertlist[triTable[cubeindex * 16 + ti]].y;
				vertex[3 * (vid + 0) + 2] = vertlist[triTable[cubeindex * 16 + ti]].z;

				vertex[3 * (vid + 1) + 0] = vertlist[triTable[cubeindex * 16 + ti + 1]].x;
				vertex[3 * (vid + 1) + 1] = vertlist[triTable[cubeindex * 16 + ti + 1]].y;
				vertex[3 * (vid + 1) + 2] = vertlist[triTable[cubeindex * 16 + ti + 1]].z;

				vertex[3 * (vid + 2) + 0] = vertlist[triTable[cubeindex * 16 + ti + 2]].x;
				vertex[3 * (vid + 2) + 1] = vertlist[triTable[cubeindex * 16 + ti + 2]].y;
				vertex[3 * (vid + 2) + 2] = vertlist[triTable[cubeindex * 16 + ti + 2]].z;

				//triangle_.resize(triIndex + 1);
				int tid = atomicAdd(triIndex, 1);
				triangle[3 * tid + 0] = vid;
				triangle[3 * tid + 1] = vid + 2;
				triangle[3 * tid + 2] = vid + 1;

				//if (vid == 0)
				//{
				//	printf("vertex[0] : %f %f %f\n", vertex[vid+0], vertex[vid+1], vertex[vid+2]);
				//	printf("triangle[0] : %d %d %d\n", triangle[0], triangle[1], triangle[2]);
				//}
			}
		}



		//if (f[tid] > -1.0f && f[tid] < 1.0f)
		//{
		//	atomicAdd(sum, f[tid]);
			//if (sum[0] < 10.0f)
			//if(triIndex[0] == 0)
			//{
			//	printf("sum=%f, hurei\n", sum[0]);
				//printf("vertlist[0] : %f %f %f\n", vertlist[0].x, vertlist[0].y, vertlist[0].z);
			//}
		//}


	}
}

void MarchingCubes::generateSurfaceMesh_cuda(float isoValue, float cutValue, const char* filename)
{
	const unsigned int nx = volumeSize_.x;
	const unsigned int ny = volumeSize_.y;
	const unsigned int nz = volumeSize_.z;

	const unsigned int gnx = nx + 1;
	const unsigned int gny = ny + 1;
	const unsigned int gnz = nz + 1;

	const float dx = voxelSize_.x;
	const float dy = voxelSize_.y;
	const float dz = voxelSize_.z;

	unsigned int numVoxels = (unsigned int)nx * ny * nz;

	// Device data allocation
	float* f_d;
	int* edgeTable_d;
	int* triTable_d;
	unsigned int* vertIndex_d;
	unsigned int* triIndex_d;
	float* vertex_d;
	unsigned int* triangle_d;
	float * vertex;
	unsigned int* triangle;

	if (cudaMalloc(&f_d, sizeof(float)*gnx*gny*gnz) != cudaSuccess) { printf("allocation error : f_d\n"); }
	if(cudaMalloc(&edgeTable_d, sizeof(int) * 256) != cudaSuccess) { printf("allocation error : edgeTable_d\n"); }
	if(cudaMalloc(&triTable_d, sizeof(int) * 256 * 16) != cudaSuccess) { printf("allocation error : triTable_d\n"); }
	if(cudaMalloc(&vertIndex_d, sizeof(unsigned int)) != cudaSuccess) { printf("allocation error : vertIndex_d\n"); }
	if(cudaMalloc(&triIndex_d, sizeof(unsigned int)) != cudaSuccess) { printf("allocation error : triIndex_d\n"); }

	// TODO : The following two numbers need to be set properly.
	int maxNumVertex = 1028 * 1028;
	int maxNumTriangles = maxNumVertex;
	vertex_.resize(maxNumVertex);
	triangle_.resize(maxNumTriangles);
	if(cudaMalloc(&vertex_d, sizeof(float) * 3 * maxNumVertex) != cudaSuccess) { printf("allocation error : vertex_d\n"); }
	if(cudaMalloc(&triangle_d, sizeof(unsigned int) * 3 * maxNumTriangles) != cudaSuccess) { printf("allocation error : triangle_d\n"); }
	cudaMallocHost(&vertex, sizeof(float) * 3 * maxNumVertex);
	cudaMallocHost(&triangle, sizeof(unsigned int) * 3 * maxNumTriangles);

	// Copy data from host to device
	if(cudaMemcpy(f_d, f_.data(), sizeof(float)*gnx*gny*gnz, cudaMemcpyHostToDevice) != cudaSuccess) { printf("copy error : f_d\n"); }
	if(cudaMemcpy(edgeTable_d, edgeTable_, sizeof(int) * 256, cudaMemcpyHostToDevice) != cudaSuccess) { printf("copy error : edgeTable_d\n"); }
	if(cudaMemcpy(triTable_d, triTable_, sizeof(int) * 256 * 16, cudaMemcpyHostToDevice) != cudaSuccess) { printf("copy error : triTable_d\n"); }

	unsigned int triIndex = 0;
	unsigned int vertIndex = 0;
	if(cudaMemcpy(vertIndex_d, &vertIndex, sizeof(unsigned int), cudaMemcpyHostToDevice) != cudaSuccess) { printf("copy error : vertIndex_d\n"); }
	if(cudaMemcpy(triIndex_d, &triIndex, sizeof(unsigned int), cudaMemcpyHostToDevice) != cudaSuccess) { printf("allocation error : triIndex_d\n"); }

	// Kernel
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	const int nthreads = 512;
	int nblocks = (numVoxels + nthreads - 1) / nthreads;
	kernel_extract_mesh << < nblocks, nthreads >> > (nx, ny, nz, numVoxels, bmin_.x, bmin_.y, bmin_.z, dx, dy, dz, 
		f_d, edgeTable_d, triTable_d,
		vertIndex_d, triIndex_d, vertex_d, triangle_d, isoValue, eps, cutValue);

	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	float elapsedTime = 0.0f;
	cudaEventElapsedTime(&elapsedTime, start, stop);
	std::cout << "- Kernel computing : " << elapsedTime << " ms" << std::endl;

	// Copy data from device to host
	cudaMemcpy(&vertIndex, vertIndex_d, sizeof(unsigned int), cudaMemcpyDeviceToHost);
	cudaMemcpy(&triIndex, triIndex_d, sizeof(unsigned int), cudaMemcpyDeviceToHost);
	std::cout << "- vertexIndex=" << vertIndex << ", triIndex=" << triIndex << std::endl;

	//vertex_.resize(vertIndex);
	//triangle_.resize(triIndex);
	//cudaMemcpy(vertex_.data(), vertex_d, sizeof(float) * 3 * (vertIndex), cudaMemcpyDeviceToHost);
	//cudaMemcpy(triangle_.data(), triangle_d, sizeof(unsigned int) * 3 * (triIndex), cudaMemcpyDeviceToHost);

	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);
	cudaMemcpy(vertex, vertex_d, sizeof(float) * 3 * (vertIndex), cudaMemcpyDeviceToHost);
	cudaMemcpy(triangle, triangle_d, sizeof(unsigned int) * 3 * (triIndex), cudaMemcpyDeviceToHost);
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&elapsedTime, start, stop);
	std::cout << "- Data copying device to host : " << elapsedTime << " ms" << std::endl;

	//vertex_.resize(vertIndex);
	//triangle_.resize(triIndex);
	//memcpy(vertex_.data(), vertex, sizeof(float) * 3 * vertIndex);
	//memcpy(triangle_.data(), triangle, sizeof(unsigned int) * 3 * triIndex);

	//Merge redundant vertices
	mergeVertices(vertIndex, vertex, triIndex, triangle);

	exportMeshInObj(filename, vertIndex, vertex, triIndex, triangle);

	// Free memory
	cudaFree(f_d);
	cudaFree(edgeTable_d);
	cudaFree(triTable_d);
	cudaFree(vertIndex_d);
	cudaFree(triIndex_d);
	cudaFree(vertex_d);
	cudaFree(triangle_d);
	cudaFreeHost(vertex);
	cudaFreeHost(triangle);
}