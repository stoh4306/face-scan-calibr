#include <vector>
#include <climits>

#define FLANN_USE_CUDA
#include <flann/flann.h>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <cuda.h>
#include <cuda_runtime.h>

float compute_psep_knn_cuda(int N, float* point, int max_nn)
{
    std::vector<float4> tempPoint(N);
    for (int i = 0; i < N; ++i)
    {
        tempPoint[i].x = point[4 * i + 0];
        tempPoint[i].y = point[4 * i + 1];
        tempPoint[i].z = point[4 * i + 2];
    }

    thrust::device_vector<float4> dataPts(N);
    cudaMemcpy((float*)thrust::raw_pointer_cast(&dataPts[0]), tempPoint.data(), sizeof(float4)*N, cudaMemcpyHostToDevice);

    thrust::device_vector<float4> query = dataPts;

    flann::Matrix<float> dataPts_matrix((float*)thrust::raw_pointer_cast(&dataPts[0]), N, 4, 4 * 4);
    flann::Matrix<float> query_matrix((float*)thrust::raw_pointer_cast(&query[0]), N, 4, 4 * 4);

    flann::KDTreeCuda3dIndexParams index_params;
    index_params["input_is_gpu_float4"] = true;

    flann::KDTreeCuda3dIndex< flann::L2<float> > index2(dataPts_matrix, index_params);
    //-----------------------------------------------------------------------------
    // NOTE : The following line should be replaced by the above.
    //        Otherwise, you have an incorrect search result or a unknown crash.
    //-----------------------------------------------------------------------------
    //flann::Index<flann::L2_Simple<float> > index2(data_device_matrix, index_params);

    //std::cout << "\n . Kd-tree built : ";
    //float stime = clock();
    index2.buildIndex();
    //float ftime = clock();
    //std::cout << (ftime - stime) / CLOCKS_PER_SEC << " sec" << std::endl;;

    //std::cout << " . KNN search : ";
    //stime = clock();
    int numNeighbors = max_nn;// 25;
    thrust::device_vector<int> indices_device(N * numNeighbors);
    thrust::device_vector<float> dists_device(N * numNeighbors);
    flann::Matrix<int> indices_device_matrix((int*)thrust::raw_pointer_cast(&indices_device[0]), N, numNeighbors);
    flann::Matrix<float> dists_device_matrix((float*)thrust::raw_pointer_cast(&dists_device[0]), N, numNeighbors);

    flann::SearchParams sp;
    sp.matrices_in_gpu_ram = true;
    index2.knnSearch(query_matrix, indices_device_matrix, dists_device_matrix, numNeighbors, sp);
    //ftime = clock();
    //std::cout << (ftime - stime) / CLOCKS_PER_SEC << " sec...";

    std::vector<float> dists_host(N*numNeighbors);
    cudaMemcpy(dists_host.data(), (float*)thrust::raw_pointer_cast(&dists_device[0]),
        sizeof(float)*N*numNeighbors, cudaMemcpyDeviceToHost);

    double averDist = 0.0;
    for (int i = 0; i < N; ++i)
    {
        averDist += sqrt(dists_host[numNeighbors*i + 1]);
    }
    if (N > 0)
    {
        averDist /= (double)N;
    }

    std::cout << "psep=" << averDist << std::endl;

    return averDist;
}

float test_compute_psep_knn_cuda(int N, float* point, int max_nn)
{
    std::vector<float4> tempPoint(N);
    for (int i = 0; i < N; ++i)
    {
        tempPoint[i].x = point[4 * i + 0];
        tempPoint[i].y = point[4 * i + 1];
        tempPoint[i].z = point[4 * i + 2];
    }

    //thrust::device_vector<float4> dataPts(N);
    //cudaMemcpy((float*)thrust::raw_pointer_cast(&dataPts[0]), tempPoint.data(), sizeof(float4)*N, cudaMemcpyHostToDevice);
    float4* dataPts;
    cudaMalloc(&dataPts, sizeof(float4)*N);
    cudaMemcpy(dataPts, tempPoint.data(), sizeof(float4)*N, cudaMemcpyHostToDevice);

    //thrust::device_vector<float4> query = dataPts;
    float4* query;
    cudaMalloc(&query, sizeof(float4)*N);
    cudaMemcpy(query, dataPts, sizeof(float4)*N, cudaMemcpyDeviceToDevice);

    //flann::Matrix<float> dataPts_matrix((float*)thrust::raw_pointer_cast(&dataPts[0]), N, 4, 4 * 4);
    //flann::Matrix<float> query_matrix((float*)thrust::raw_pointer_cast(&query[0]), N, 4, 4 * 4);
    flann::Matrix<float> dataPts_matrix((float*)dataPts, N, 4, 4 * 4);
    flann::Matrix<float> query_matrix((float*)query, N, 4, 4 * 4);

    flann::KDTreeCuda3dIndexParams index_params;
    index_params["input_is_gpu_float4"] = true;

    flann::KDTreeCuda3dIndex< flann::L2<float> > index2(dataPts_matrix, index_params);
    //-----------------------------------------------------------------------------
    // NOTE : The following line should be replaced by the above.
    //        Otherwise, you have an incorrect search result or a unknown crash.
    //-----------------------------------------------------------------------------
    //flann::Index<flann::L2_Simple<float> > index2(data_device_matrix, index_params);

    //std::cout << "\n . Kd-tree built : ";
    //float stime = clock();
    index2.buildIndex();
    //float ftime = clock();
    //std::cout << (ftime - stime) / CLOCKS_PER_SEC << " sec" << std::endl;;

    //std::cout << " . KNN search : ";
    //stime = clock();
    int numNeighbors = max_nn;// 25;
    //thrust::device_vector<int> indices_device(N * numNeighbors);
    //thrust::device_vector<float> dists_device(N * numNeighbors);
    int* indices_device;
    float* dists_device;
    cudaMalloc(&indices_device, sizeof(int)*N*numNeighbors);
    cudaMalloc(&dists_device, sizeof(float)*N*numNeighbors);

    //flann::Matrix<int> indices_device_matrix((int*)thrust::raw_pointer_cast(&indices_device[0]), N, numNeighbors);
    //flann::Matrix<float> dists_device_matrix((float*)thrust::raw_pointer_cast(&dists_device[0]), N, numNeighbors);
    flann::Matrix<int> indices_device_matrix(indices_device, N, numNeighbors);
    flann::Matrix<float> dists_device_matrix(dists_device, N, numNeighbors);

    flann::SearchParams sp;
    sp.matrices_in_gpu_ram = true;
    index2.knnSearch(query_matrix, indices_device_matrix, dists_device_matrix, numNeighbors, sp);
    //ftime = clock();
    //std::cout << (ftime - stime) / CLOCKS_PER_SEC << " sec...";

    std::vector<float> dists_host(N*numNeighbors);
    cudaMemcpy(dists_host.data(), dists_device,
        sizeof(float)*N*numNeighbors, cudaMemcpyDeviceToHost);

    double averDist = 0.0;
    for (int i = 0; i < N; ++i)
    {
        averDist += sqrt(dists_host[numNeighbors*i + 1]);
    }
    if (N > 0)
    {
        averDist /= (double)N;
    }

    std::cout << "psep=" << averDist << std::endl;

    // Free device memory
    cudaFree(dataPts);
    cudaFree(query);
    cudaFree(indices_device);
    cudaFree(dists_device);

    return averDist;
}

clock_t start_time_;

void start_timer(const std::string& message = "")
{
    if (!message.empty()) {
        printf("%s", message.c_str());
        fflush(stdout);
    }
    start_time_ = clock();
}

double stop_timer()
{
    return double(clock() - start_time_) / CLOCKS_PER_SEC;
}

__device__ float eps = 1.0e-6f;

__device__ void copy_array(float* dest, float* src, size_t N)
{
    for (size_t i = 0; i < N; ++i)
        dest[i] = src[i];
}

__device__ void symSchur2_gpu(int N, float* A, int p, int q, float& c, float& s)
{
    //--------------------------------------
    // ASSUME : A is symmetric and p < q
    //--------------------------------------
    const float A_pq = A[N*p + q];

    if (A_pq > eps || A_pq < -eps)
    {
        const float A_pp = A[N*p + p];
        const float A_qq = A[N*q + q];
        //if (p == 0 && q == 1)
        //{
        //    std::cout << A_pp << " " << A_qq << " " << A_pq << std::endl;
        //}

        float tau = (A_qq - A_pp) / (2.0f*A_pq);

        float t;
        if (tau >= 0.0f)
            t = 1.0f / (tau + sqrt(1.0f + tau*tau));
        else
            t = 1.0f / (tau - sqrt(1.0f + tau*tau));

        c = 1.0f / sqrt(1.0f + t*t);
        s = t*c;
    }
    else
    {
        c = 1.0f;
        s = 0.0f;
    }
}

__device__ void diagonalize(int N, float* B, float* J, float* T)
{
    // B = B * J
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            float temp = 0.0f;

            for (int k = 0; k < N; ++k)
            {
                temp += B[N*i + k] * J[N*k + j];
            }

            T[N*i + j] = temp;
        }
    }

    // B = J^t * B
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            float temp = 0.0f;

            for (int k = 0; k < N; ++k)
            {
                temp += J[N*k + i] * T[N*k + j];
            }

            B[N*i + j] = temp;
        }
    }
}

__device__ void updateRotation(int N, float* V, float* J, float* T)
{
    // B = B * J
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            float temp = 0.0f;

            for (int k = 0; k < N; ++k)
            {
                temp += V[N*i + k] * J[N*k + j];
            }

            T[N*i + j] = temp;
        }
    }

    copy_array(V, T, N*N);
}


__global__ void kernel_compute_normals(int numPts, float* p,
    int numNeighbors, int* indices, float* normals)
{
    int tid = threadIdx.x + blockIdx.x*blockDim.x;

    if (tid < numPts)
    {
        // compute covariance matrix
        float A[9] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };

        float x, y, z, tx, ty, tz;
        x = p[4 * tid + 0];
        y = p[4 * tid + 1];
        z = p[4 * tid + 2];

        for (int j = 0; j < numNeighbors; ++j)
        {
            const int nid = indices[numNeighbors*tid + j];
            tx = p[4 * nid + 0];
            ty = p[4 * nid + 1];
            tz = p[4 * nid + 2];

            A[0] += (x - tx)*(x - tx);
            A[1] += (x - tx)*(y - ty);
            A[2] += (x - tx)*(z - tz);

            A[4] += (y - ty)*(y - ty);
            A[5] += (y - ty)*(z - tz);

            A[8] += (z - tz)*(z - tz);
        }

        A[3] = A[1];
        A[6] = A[2];
        A[7] = A[5];

        float off_diag = 0.0f;

        for (int i = 0; i < 3; ++i)
        {
            for (int j = i + 1; j < 3; j++)
                off_diag += A[3 * i + j] * A[3 * i + j];
        }

        //printf("off(A)=%f\n", off_diag);

        // Compute L2-norm of diagonals
        float diag = 0.0f;
        for (int i = 0; i < 3; ++i)
        {
            diag += A[3 * i + i] * A[3 * i + i];
        }

        float normA = diag + 2.0f*off_diag;

        // Initialize V as the identity matrix
        float I[9] = { 1.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 1.0f };

        float V[9], T[9];

        copy_array(V, I, 9);

        float J[9];

        float tol = eps*normA;

        int it = 0;

        while (off_diag >= tol)
        {
            for (int p = 0; p < 3; ++p)
            {
                for (int q = p + 1; q < 3; ++q)
                {
                    float c, s;
                    symSchur2_gpu(3, A, p, q, c, s);
                    //printf("(%d, %d) : c=%f, s=%f\n", p, q, c, s);

                    copy_array(J, I, 9);

                    J[3 * p + p] = J[3 * q + q] = c;
                    J[3 * p + q] = s;
                    J[3 * q + p] = -s;

                    off_diag -= A[3 * p + q] * A[3 * p + q];

                    diagonalize(3, A, J, T);

                    updateRotation(3, V, J, T);
                }
            }

            it++;
        }

        float lambda[3];
        lambda[0] = A[0];
        lambda[1] = A[4];
        lambda[2] = A[8];

        // Find the smallest eigenvalue
        int sid = 0;
        if (lambda[sid] > lambda[1])    sid = 1;
        if (lambda[sid] > lambda[2])    sid = 2;

        // Set the corresponding eigenvector to the normal
        float tn[3];
        tn[0] = V[3 * 0 + sid];
        tn[1] = V[3 * 1 + sid];
        tn[2] = V[3 * 2 + sid];

        // Check if the normal is going towared to the view point
        float viewPoint[3] = { 0.0f, 0.0f, 0.0f };
        float dot = tn[0] * (tx - viewPoint[0]) + tn[1] * (ty - viewPoint[1]) + tn[2] * (tz - viewPoint[2]);
        if (dot > 0.0f)
        {
            tn[0] = -tn[0];
            tn[1] = -tn[1];
            tn[2] = -tn[2];
        }

        normals[4 * tid + 0] = tn[0];
        normals[4 * tid + 1] = tn[1];
        normals[4 * tid + 2] = tn[2];
        normals[4 * tid + 3] = lambda[sid] / (lambda[0] + lambda[1] + lambda[2]); // 0.0f;

        if (fabs(tn[0]) < eps && fabs(tn[1]) < eps && fabs(tn[2]) < eps)
        {
            printf("p[%d] : zero normals\n", tid);
        }
    }// end if(tid<numPts)
}

void compute_knn_gpu_2(int numPts, float* point, int nknn, float searchRadius,
    float* normals)
{
    thrust::device_vector<float4> dataPts(numPts);
    cudaMemcpy((float*)thrust::raw_pointer_cast(&dataPts[0]), point, sizeof(float4)*numPts, cudaMemcpyHostToDevice);

    thrust::device_vector<float4> query = dataPts;

    flann::Matrix<float> dataPts_matrix((float*)thrust::raw_pointer_cast(&dataPts[0]), numPts, 4, 4 * 4);
    flann::Matrix<float> query_matrix((float*)thrust::raw_pointer_cast(&query[0]), numPts, 4, 4 * 4);

    flann::KDTreeCuda3dIndexParams index_params;
    index_params["input_is_gpu_float4"] = true;

    flann::KDTreeCuda3dIndex< flann::L2<float> > index2(dataPts_matrix, index_params);
    //-----------------------------------------------------------------------------
    // NOTE : The following line should be replaced by the above.
    //        Otherwise, you have an incorrect search result or a unknown crash.
    //-----------------------------------------------------------------------------
    //flann::Index<flann::L2_Simple<float> > index2(data_device_matrix, index_params);

    //start_timer("- Building kd-tree index...");
    index2.buildIndex();
    //printf("done (%g seconds)\n", stop_timer());

    int numNeighbors = nknn;// 25;
    thrust::device_vector<int> indices_device(numPts * numNeighbors);
    thrust::device_vector<float> dists_device(numPts * numNeighbors);
    flann::Matrix<int> indices_device_matrix((int*)thrust::raw_pointer_cast(&indices_device[0]), numPts, numNeighbors);
    flann::Matrix<float> dists_device_matrix((float*)thrust::raw_pointer_cast(&dists_device[0]), numPts, numNeighbors);

    //start_timer("- Searching KNN...");
    //indices.cols = 4;
    //dists.cols = 4;
    flann::SearchParams sp;
    sp.matrices_in_gpu_ram = true;
    index2.knnSearch(query_matrix, indices_device_matrix, dists_device_matrix, numNeighbors, sp);
    //printf("done (%g seconds)\n", stop_timer());
    //std::cout << *(indices_device_matrix[numNeighbors - 1]) << std::endl;

    //----------------------------------------------------------------
    // Compute covariance matrix and normal vector w.r.t the view point
    //----------------------------------------------------------------
    cudaEvent_t start, stop;
    float elapsedTime;

    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start, 0);

    float* d_normals;
    cudaMalloc(&d_normals, sizeof(float4)*numPts);

    const int nthreads = 256;
    const int nblocks = (numPts + nthreads - 1) / nthreads;

    kernel_compute_normals << <nblocks, nthreads >> >(numPts, (float*)thrust::raw_pointer_cast(&dataPts[0]),
        numNeighbors, thrust::raw_pointer_cast(&indices_device[0]), d_normals);

    cudaMemcpy(normals, d_normals, sizeof(float4)*numPts, cudaMemcpyDeviceToHost);

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);
    //printf("kernel comp. time = %f ms\n", elapsedTime);

    cudaFree(d_normals);

    cudaEventDestroy(start);
    cudaEventDestroy(stop);
}

void test_compute_knn_gpu_2(int numPts, float* point, int nknn, float searchRadius,
    float* normals)
{
    //thrust::device_vector<float4> dataPts(numPts);
    //cudaMemcpy((float*)thrust::raw_pointer_cast(&dataPts[0]), point, sizeof(float4)*numPts, cudaMemcpyHostToDevice);
    float4* dataPts;
    cudaMalloc(&dataPts, sizeof(float4)*numPts);
    cudaMemcpy(dataPts, point, sizeof(float4)*numPts, cudaMemcpyHostToDevice);

    //thrust::device_vector<float4> query = dataPts;
    float4* query;
    cudaMalloc(&query, sizeof(float4)*numPts);
    cudaMemcpy(query, dataPts, sizeof(float4)*numPts, cudaMemcpyDeviceToDevice);

    //flann::Matrix<float> dataPts_matrix((float*)thrust::raw_pointer_cast(&dataPts[0]), numPts, 4, 4 * 4);
    //flann::Matrix<float> query_matrix((float*)thrust::raw_pointer_cast(&query[0]), numPts, 4, 4 * 4);
    flann::Matrix<float> dataPts_matrix((float*)dataPts, numPts, 4, 4 * 4);
    flann::Matrix<float> query_matrix((float*)query, numPts, 4, 4 * 4);

    flann::KDTreeCuda3dIndexParams index_params;
    index_params["input_is_gpu_float4"] = true;

    flann::KDTreeCuda3dIndex< flann::L2<float> > index2(dataPts_matrix, index_params);
    //-----------------------------------------------------------------------------
    // NOTE : The following line should be replaced by the above.
    //        Otherwise, you have an incorrect search result or a unknown crash.
    //-----------------------------------------------------------------------------
    //flann::Index<flann::L2_Simple<float> > index2(data_device_matrix, index_params);

    //start_timer("- Building kd-tree index...");
    index2.buildIndex();
    //printf("done (%g seconds)\n", stop_timer());

    int numNeighbors = nknn;// 25;
    //thrust::device_vector<int> indices_device(numPts * numNeighbors);
    //thrust::device_vector<float> dists_device(numPts * numNeighbors);
    int* indices_device;
    float* dists_device;
    cudaMalloc(&indices_device, sizeof(int)*numPts*numNeighbors);
    cudaMalloc(&dists_device, sizeof(float)*numPts*numNeighbors);

    //flann::Matrix<int> indices_device_matrix((int*)thrust::raw_pointer_cast(&indices_device[0]), numPts, numNeighbors);
    //flann::Matrix<float> dists_device_matrix((float*)thrust::raw_pointer_cast(&dists_device[0]), numPts, numNeighbors);
    flann::Matrix<int> indices_device_matrix(indices_device, numPts, numNeighbors);
    flann::Matrix<float> dists_device_matrix(dists_device, numPts, numNeighbors);

    //start_timer("- Searching KNN...");
    //indices.cols = 4;
    //dists.cols = 4;
    flann::SearchParams sp;
    sp.matrices_in_gpu_ram = true;
    index2.knnSearch(query_matrix, indices_device_matrix, dists_device_matrix, numNeighbors, sp);
    //printf("done (%g seconds)\n", stop_timer());
    //std::cout << *(indices_device_matrix[numNeighbors - 1]) << std::endl;

    //----------------------------------------------------------------
    // Compute covariance matrix and normal vector w.r.t the view point
    //----------------------------------------------------------------
    cudaEvent_t start, stop;
    float elapsedTime;

    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start, 0);

    float* d_normals;
    cudaMalloc(&d_normals, sizeof(float4)*numPts);

    const int nthreads = 256;
    const int nblocks = (numPts + nthreads - 1) / nthreads;

    kernel_compute_normals << <nblocks, nthreads >> >(numPts, (float*)thrust::raw_pointer_cast(&dataPts[0]),
        numNeighbors, thrust::raw_pointer_cast(&indices_device[0]), d_normals);

    cudaMemcpy(normals, d_normals, sizeof(float4)*numPts, cudaMemcpyDeviceToHost);

    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime, start, stop);
    //printf("kernel comp. time = %f ms\n", elapsedTime);

    cudaFree(d_normals);

    cudaFree(dataPts);
    cudaFree(query);
    cudaFree(indices_device);
    cudaFree(dists_device);

    cudaEventDestroy(start);
    cudaEventDestroy(stop);
}