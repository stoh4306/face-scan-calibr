#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
//#include <pcl/point_representation.h>
#include <pcl/io/pcd_io.h>

void compute_knn_gpu_2(int numPts, float* point, int nknn, float searchRadius,
        float* normals);
void test_compute_knn_gpu_2(int numPts, float* point, int nknn, float searchRadius,
    float* normals);
float compute_psep_knn_cuda(int N, float*, int max_nn);
float test_compute_psep_knn_cuda(int N, float*, int max_nn);

bool loadPCDFile(const char* input_pcd_file, pcl::PointCloud<pcl::PointXYZRGB>& pointCloud)
{
    // NOTE : We are using PCLPointCloud2 
    //        because we don't know what fields the input data has.
    pcl::PCLPointCloud2 tempCloud;

    // load pcd data
    if (pcl::io::loadPCDFile(input_pcd_file, tempCloud)<0)
    {
        return false;
    }

    pcl::fromPCLPointCloud2(tempCloud, pointCloud);

    return true;
}

bool writePCDFile(const char* filename, pcl::PointCloud<pcl::PointXYZRGBNormal>& pointCloud)
{
    pcl::PCDWriter writer;
    if (writer.writeBinary(std::string(filename), pointCloud) < 0)
    {
        return false;
    }

    return true;
}


void compute_normals_cuda(pcl::PointCloud< pcl::PointXYZRGB > & inCloud,
        pcl::PointCloud< pcl::PointXYZRGBNormal > & outCloud, int nknn) //float radius)
{
    float stime = clock();

    // Extract xyz field
    pcl::PointCloud<pcl::PointXYZ> posCloud;
    pcl::copyPointCloud(inCloud, posCloud);

    std::cout << "#points = " << posCloud.points.size();


    //float radius = 0.01;
    //float psep = compute_psep_knn_cuda(posCloud.size(), (float*)posCloud.points.data(), nknn);
    float psep = test_compute_psep_knn_cuda(posCloud.size(), (float*)posCloud.points.data(), nknn);
    float radius = 2.4f*psep;

    struct Float4 { float x, y, z, w;};
    std::vector<Float4> normals(posCloud.size());
    test_compute_knn_gpu_2(posCloud.size(), (float*)posCloud.points.data(), nknn, radius,
        (float*)normals.data());

    pcl::copyPointCloud(inCloud, outCloud);

    for (int i = 0; i < inCloud.size(); ++i)
    {
        outCloud.points[i].data_n[0] = normals[i].x;
        outCloud.points[i].data_n[1] = normals[i].y;
        outCloud.points[i].data_n[2] = normals[i].z;
        outCloud.points[i].data_n[3] = 0.0f;
        outCloud.points[i].curvature = normals[i].w;
    }

    float ftime = clock();
    std::cout << " -> cuda normal computation(" << (ftime - stime) / CLOCKS_PER_SEC << " sec)" << std::endl;
}


int main(int argc, char** argv)
{
    if (argc <= 3)
    {
        std::cout << "- Usage : " << argv[0] << " input_cloud output_cloud num_neighbors" << "\n"
            << " (ex) compute_normals_cuda.exe p_0000.pcd q_0000.pcd 64" << std::endl;
        return -1;
    }

    const char* inCloudFile = argv[1];
    const char* outCloudFile = argv[2];
    //float radius = atof(argv[3]);
    float nknn = atoi(argv[3]);

    pcl::PointCloud< pcl::PointXYZRGB >         inCloud;
    pcl::PointCloud< pcl::PointXYZRGBNormal >   inCloud2, outCloud;
    
    loadPCDFile(inCloudFile, inCloud);

    compute_normals_cuda(inCloud, outCloud, nknn); // radius);

    writePCDFile(outCloudFile, outCloud);

    return 0;
}