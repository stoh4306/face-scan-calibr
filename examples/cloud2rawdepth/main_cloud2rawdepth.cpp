 #include <fstream>
 #include <iostream>
 #include <vector>
 
#include <Eigen/Dense>

#include "reconst.h"
#include "calibr.h"

bool loadCloudInPly(const char* filename, std::vector<Eigen::Vector3f>& cloud)
{
    std::cout << "- Reading " << filename << "...";

    std::ifstream inFile(filename, std::ios::binary);

    if (!inFile.is_open())
    {
        std::cerr << "Error, can't open PLY file for reading : " << filename << std::endl;
        return false;
    }

    std::stringstream buffer;
    buffer << inFile.rdbuf();

    std::string temp;

    int num_pts = 0;

    bool isBinary = true;
    bool hasPoint = false, hasColor = false, hasNormal = false;

    // Header
    std::vector<std::pair<std::string, std::string>>	vertexProperty;

    do
    {
        if (temp == "format")
        {
            buffer >> temp;
            if (temp == "ascii") isBinary = false;
        }
        else if (temp == "vertex")
        {
            buffer >> num_pts;

            cloud.resize(num_pts);
        }
        else if (temp == "property")
        {
            if (num_pts > 0)
            {
                std::string type, name;
                buffer >> type;
                buffer >> name;

                vertexProperty.push_back(std::make_pair(type, name));

                if (name == "x")
                {
                    hasPoint = true;
                }
                else if (name == "red")
                {
                    hasColor = true;
                }
                else if (name == "nx")
                {
                    hasNormal = true;
                }
            }
        }

        buffer >> temp;

    } while (temp != "end_header");

    std::cout << "num pts=" << num_pts << "...";
    std::cout << "position=" << hasPoint << ", color=" << hasColor << ", normal=" << hasNormal << std::endl;

    if (isBinary)
    {
        std::getline(buffer, temp, '\n');

        unsigned char cdata[3];
        Eigen::Vector3f pn;

        for (int i = 0; i < num_pts; ++i)
        {
            if (hasPoint)
            {
                buffer.read((char*)&cloud[i], sizeof(Eigen::Vector3f));
                //std::cout << i << " : " << cloud[i].transpose() << std::endl;
            }
            if (hasColor)
            {
                buffer.read((char*)cdata, sizeof(unsigned char) * 3);
                //col_[i] = Vector3f(cdata[2], cdata[1], cdata[0]);
            }
            if (hasNormal)
                buffer.read((char*)pn.data(), sizeof(Eigen::Vector3f));

            if (vertexProperty.size() > 9)
            {
                float curv;
                buffer.read((char*)&curv, sizeof(float));
            }
        }
    }
    else // ASCII
    {
        unsigned char cdata[3];
        Eigen::Vector3f pn;

        for (int i = 0; i < num_pts; ++i)
        {
            if (hasPoint)
                buffer >> cloud[i](0) >> cloud[i](1) >> cloud[i](2);
            if (hasColor)
                buffer >> cdata[0] >> cdata[1] >> cdata[2];
            if (hasNormal)
                buffer >> pn(0) >> pn(1) >> pn(2);

            if (vertexProperty.size() > 9)
            {
                float curv;
                buffer >> curv;
            }
        }
    }

    //std::cout << pos_[0] << " " << col_[0] << " " << normal_[0] << std::endl;
    //std::cout << pos_[1] << " " << col_[1] << " " << normal_[1] << std::endl;
    //std::cout << pos_[num_pts - 1] << " " << col_[num_pts - 1] << " " << normal_[num_pts - 1] << std::endl;

    inFile.close();

    return true;
}

int main(int argc, char** argv)
{
    // Load point cloud
    std::vector<Eigen::Vector3f> cloud;
    loadCloudInPly(argv[1], cloud);
    std::cout << "- Point cloud reading done" << std::endl;

    // Load calibration file
    Eigen::Matrix3f K;
    Eigen::Matrix4f CP;
    read_calibrations(argv[2], argv[3], K.data(), CP.data());
    std::cout << "- Calibration data loaded" << std::endl;

    Eigen::Matrix4f T = CP.inverse();

    // Convert to raw depth data
    int n = (int)cloud.size();
    const int width = atoi(argv[4]);
    const int height = atoi(argv[5]);
    std::cout << "- Depth :" << width << " x " << height << std::endl;

    std::vector<float> rawDepth(width*height);
    point2rawdepth(n, cloud.data()->data(), K.data(), T.data(),  width, height, rawDepth.data(), "rawdepth.data");

    return 0;
}