/*==========================================================================
Copyright by:
    Macrograph, Research

    For further information, contact:
        Dr. Seungtaik Oh

        E-mail : stoh@macrograph.co.kr
        Phone  : +82 10 2923 6290
        Address : Seoul, Korea

This copyright notice must be included with all copies of the source codes.
============================================================================*/
#pragma once

#include "vec.h"

#include <vector>
#include <fstream>
#include <sstream>

namespace mg
{
template <typename T>
class TriMesh
{
public:
    void resizeAllArrays(int N)
    {
        vertex_.resize(N);
        vertexNormal_.resize(N);
        texCoor_.resize(N);
        triangle_.resize(N);
        tri_tex_.resize(N);
        tri_normal_.resize(N);
    }

    void deepCopyTo(TriMesh<T>& outMesh)
    {
        outMesh.vertex_ = vertex_;
        outMesh.triangle_ = triangle_;
        outMesh.tri_tex_ = tri_tex_;
        outMesh.tri_normal_ = tri_normal_;
        outMesh.vertexNormal_ = vertexNormal_;
        outMesh.texCoor_ = texCoor_;
    }

    int numVertices()               { return (int)vertex_.size(); }
    int numVertexNormals()          { return (int)vertexNormal_.size(); }
    int numTexCoords()              { return (int)texCoor_.size(); }
    int numTriangles()              { return (int)triangle_.size(); }
    int numTriangleTextures()       { return (int)tri_tex_.size(); }
    int numTriangleNormals()        { return (int)tri_normal_.size(); }

    Vector3<T>* vertexData()        { if (vertex_.size() > 0) return &vertex_[0]; else return NULL; }
    Vector3i*   triangleData()      { if (triangle_.size() > 0) return &triangle_[0]; else return NULL; }
    Vector3<T>* vertexNormalData()  { if (vertexNormal_.size() > 0) return &vertexNormal_[0]; else return NULL; }
    Vector2<T>* texCoorData()       { if (texCoor_.size() > 0) return &texCoor_[0]; else return NULL; }
    Vector3i*   triangleTextureData() { if (tri_tex_.size() > 0) return &tri_tex_[0]; else return NULL; }
    Vector3i*   triangleNormalData() { if (tri_normal_.size() > 0) return &tri_normal_[0]; else return NULL; }

    void computeBoundingBox(mg::Vector3<T>& bmin, mg::Vector3<T>& bmax)
    {
        T INF = std::numeric_limits<T>::infinity();
        bmin = mg::Vector3<T>(INF, INF, INF);
        bmax = -bmin;

        for (int i = 0; i < (int)numVertices(); ++i)
        {
            if (vertex_[i].x < bmin.x)  bmin.x = vertex_[i].x;
            if (vertex_[i].y < bmin.y)  bmin.y = vertex_[i].y;
            if (vertex_[i].z < bmin.z)  bmin.z = vertex_[i].z;

            if (vertex_[i].x > bmax.x)  bmax.x = vertex_[i].x;
            if (vertex_[i].y > bmax.y)  bmax.y = vertex_[i].y;
            if (vertex_[i].z > bmax.z)  bmax.z = vertex_[i].z;
        }
    }

    void computeInertiaMatrix(Vector3<T> center, T* Inertia)
    {
        std::memset(Inertia, 0, sizeof(T) * 9);

        int N = numVertices();
        for (int i = 0; i < N; ++i)
        {
            Vector3<T> p = vertex_[i]-center;

            Inertia[0] += p.y*p.y + p.z*p.z;
            Inertia[1] += -p.x*p.y;
            Inertia[2] += -p.x*p.z;

            Inertia[4] += p.x*p.x + p.z*p.z;
            Inertia[5] += -p.y*p.z;
            
            Inertia[8] += p.x*p.x + p.y*p.y;
        }

        Inertia[3] = Inertia[1];
        Inertia[6] = Inertia[2];
        Inertia[7] = Inertia[5];
    }

    Vector3<T> computeCenterOfMass()
    {
        Vector3<T> center(0.0, 0.0, 0.0);

        int N = numVertices();
        for (int i = 0; i < N; ++i)
        {
            center += vertex_[i];
        }
        
        if (N > 0)
        {
            center *= (T)(1.0 / N);
        }

        return center;
    }

    void computeFaceNormal(std::vector<Vector3<T>>& faceNormal)
    {
        int numTri = (int)triangle_.size();
        faceNormal.resize(numTri);

        for (int t = 0; t < numTri; ++t)
        {
            Vector3<T> v[3];
            v[0] = vertex_[triangle_[t].x];
            v[1] = vertex_[triangle_[t].y];
            v[2] = vertex_[triangle_[t].z];

            computeTriangleNormal(v, faceNormal[t]);
        }
    }

    void computeVertexNormal(std::vector<Vector3<T>>& faceNormal, const Vector3<T>& defaultNormal)
    {
        // Find neighbor triangle
        int numVer = this->numVertices();

        std::vector< std::vector<int> > neiTri(numVer);

        int numTri = this->numTriangles();

        for (int t = 0; t < numTri; ++t)
        {
            neiTri[triangle_[t].x].push_back(t);
            neiTri[triangle_[t].y].push_back(t);
            neiTri[triangle_[t].z].push_back(t);
        }

        // Compute vertex normal
        vertexNormal_.resize(numVer);
        for (int i = 0; i < numVer; ++i)
        {
            int numNeiTri = (int)neiTri[i].size();

            if (numNeiTri > 0)
            {
                vertexNormal_[i] = Vector3<T>();

                std::vector< int >& nTri = neiTri[i];
                for (int j = 0; j < numNeiTri; ++j)
                {
                    vertexNormal_[i] += faceNormal[nTri[j]];
                }

                vertexNormal_[i].normalize();
            }
            else
            {
                vertexNormal_[i] = defaultNormal;
            }
        }
    }

    void getTwoOtherVertices(int vi, int ti, int& v1, int& v2)
    {
        if (triangle_[ti].x == vi)
        {
            v1 = triangle_[ti].y;
            v2 = triangle_[ti].z;
        }
        else if (triangle_[ti].y == vi)
        {
            v1 = triangle_[ti].z;
            v2 = triangle_[ti].x;
        }
        else //if (triangle_[ti].z == vi)
        {
            v1 = triangle_[ti].x;
            v2 = triangle_[ti].y;
        }
    }
	int loadMeshFromUNITY3D(const char* buffer) {
		
		int seek = 0;
		uint32_t pointNum = *((uint32_t*)buffer);
		seek += sizeof(uint32_t);
		uint32_t polygonNum = *((uint32_t*)(buffer+seek));
		seek += sizeof(uint32_t);
		vertex_.resize(pointNum);
		vertexNormal_.resize(pointNum);
		triangle_.resize(polygonNum);
		memcpy(&vertex_[0], (buffer + seek), sizeof(float)*pointNum * 3);
		seek += (sizeof(float)*pointNum * 3);
		memcpy(&vertexNormal_[0], (buffer + seek), sizeof(float)*pointNum * 3);
		seek += (sizeof(float)*pointNum * 3);
		memcpy(&triangle_[0], (buffer + seek), sizeof(int)*polygonNum * 3);
		seek += (sizeof(int)*polygonNum * 3);
		return seek;
	}
	int writeMeshToUNITY3D(char** buffer) {
		uint32_t pointNum = vertex_.size();
		uint32_t polygonNum = triangle_.size();
		char* reBuffer = new char[sizeof(uint32_t)*2 + sizeof(float)*pointNum * 8 + sizeof(int)*polygonNum * 3];
		int seek = 0;
		//uint16_t pointNum = *buffer;
		*((uint32_t*)reBuffer) = pointNum;
		seek += sizeof(uint32_t);
		//uint16_t polygonNum = *(buffer + seek);
		*((uint32_t*)(reBuffer + seek)) = polygonNum;
		seek += sizeof(uint32_t);
		memcpy((reBuffer + seek), &vertex_[0], sizeof(float)*pointNum * 3);
		seek += (sizeof(float)*pointNum * 3);
		memcpy((reBuffer + seek), &vertexNormal_[0], sizeof(float)*pointNum * 3);
		seek += (sizeof(float)*pointNum * 3);
		memcpy((reBuffer + seek), &texCoor_[0], sizeof(float)*pointNum * 2);
		seek += (sizeof(float)*pointNum * 2);

		memcpy((reBuffer + seek), &triangle_[0], sizeof(int)*polygonNum * 3);
		seek += (sizeof(int)*polygonNum * 3);
		*buffer = reBuffer;
		return seek;
	}

    bool loadMeshInOBJ(const char* filename)
    {
        std::ifstream inFile(filename);
        if (!inFile.is_open())
        {
            std::cerr << "Error! Can't open mesh file " << filename << " for reading." << std::endl;
            return false;
        }
        std::cout << "- Loading an OBJ mesh : " << filename;
        // Resize all data arrays
        resizeAllArrays(0);

        std::string temp;
        inFile >> temp;
        while (!inFile.eof())
        {
            if (temp == "v")
            {
                T x, y, z;
                inFile >> x >> y >> z;
                std::getline(inFile, temp, '\n');

                vertex_.push_back(Vector3<T>(x, y, z));
            }
            else if (temp == "vt")
            {
                T x, y;
                inFile >> x >> y;
                std::getline(inFile, temp, '\n');

                texCoor_.push_back(Vector2<T>(x, y));
            }
            else if (temp == "vn")
            {
                T x, y, z;
                inFile >> x >> y >> z;
                std::getline(inFile, temp, '\n');

                vertexNormal_.push_back(Vector3<T>(x, y, z));
            }
            else if (temp == "f")
            {
                std::ostringstream viss[3], niss[3], tiss[3];

                for (int fi = 0; fi < 3; ++fi)
                {
                    inFile >> temp;

                    int count = 0;
                    for (int i = 0; i < (int)temp.size(); ++i)
                    {
                        if (temp[i] == '/')
                        {
                            count++;
                            continue;
                        }

                        if (count == 0)
                        {
                            viss[fi] << temp[i];
                        }
                        else if (count == 1)
                        {
                            tiss[fi] << temp[i];
                        }
                        else if (count == 2)
                        {
                            niss[fi] << temp[i];
                        }
                    }
                }

                if (viss[0].str().size() > 0 && viss[1].str().size() > 0 && viss[2].str().size() > 0)
                {
                    const int vx = atoi(viss[0].str().c_str())-1;
                    const int vy = atoi(viss[1].str().c_str())-1;
                    const int vz = atoi(viss[2].str().c_str())-1;
                    triangle_.push_back(Vector3i(vx, vy, vz));
                }

                if (tiss[0].str().size() > 0 && tiss[1].str().size() > 0 && tiss[2].str().size() > 0)
                {
                    const int tx = atoi(tiss[0].str().c_str())-1;
                    const int ty = atoi(tiss[1].str().c_str())-1;
                    const int tz = atoi(tiss[2].str().c_str())-1;
                    tri_tex_.push_back(Vector3i(tx, ty, tz));
                }

                if (niss[0].str().size() > 0 && niss[1].str().size() > 0 && niss[2].str().size() > 0)
                {
                    const int vx = atoi(viss[0].str().c_str()) - 1;
                    const int vy = atoi(viss[1].str().c_str()) - 1;
                    const int vz = atoi(viss[2].str().c_str()) - 1;

                    const int nx = atoi(niss[0].str().c_str()) - 1;
                    const int ny = atoi(niss[1].str().c_str()) - 1;
                    const int nz = atoi(niss[2].str().c_str()) - 1;

                    tri_normal_.push_back(Vector3i(nx, ny, nz));

                    //if (vx != nx || vy != ny || vz != nz)
                    //{
                    //    std::cerr << "Warning! inconsistent triangle vertex normal!" << std::endl;
                    //}
                }

            }
            inFile >> temp;
            //std::cout << temp << std::endl;
        }

        inFile.close();

        std::cout << "...done\n" 
            << "(#v=" << numVertices() << ", #t=" << numTriangles() << ", #vn=" << numVertexNormals() 
            << ", #vt=" << numTexCoords() << ", #ft=" << numTriangleTextures() << ")" 
            << std::endl;

        return true;
    }

    bool writeMeshInOBJ(const char* filename, bool hasVertexNormal, bool hasTexCoor)
    {
        std::ofstream outFile(filename);

        if (!outFile.is_open())
        {
            std::cerr << "ERROR! Can't open " << filename << " for writing" << std::endl;
            return false;
        }
        int numVertex = (int)vertex_.size();
        for (int i = 0; i < numVertex; ++i)
        {
            outFile << "v " << vertex_[i].x << " " << vertex_[i].y << " " << vertex_[i].z << "\n";
        }

        if (hasVertexNormal)
        {
            int numVNormals = numVertexNormals();
            //std::cout << "numNormals=" << numVNormals << std::endl;

            for (int i = 0; i < numVNormals; ++i)
            {
                outFile << "vn " << vertexNormal_[i].x << " " << vertexNormal_[i].y << " " << vertexNormal_[i].z << "\n";
            }
        }

        if (hasTexCoor)
        {
            int numTex = numTexCoords();
            for (int i = 0; i < numTex; ++i)
            {
                outFile << "vt " << texCoor_[i].x << " " << texCoor_[i].y << "\n";
            }
        }

        int numTri = (int)triangle_.size();
        //std::cout << "numTri=" << numTri << ", " << hasVertexNormal << " " << hasTexCoor << std::endl;

        if (hasVertexNormal && hasTexCoor)
        {
            for (int i = 0; i < numTri; ++i)
            {
                outFile << "f " << triangle_[i].x + 1 << "/" << tri_tex_[i].x + 1 << "/" << tri_normal_[i].x + 1
                    << " " << triangle_[i].y + 1 << "/" << tri_tex_[i].y + 1 << "/" << tri_normal_[i].y + 1
                    << " " << triangle_[i].z + 1 << "/" << tri_tex_[i].z + 1 << "/" << tri_normal_[i].z + 1 << "\n";
            }
        }
        else if (!hasVertexNormal && hasTexCoor)
        {
            for (int i = 0; i < numTri; ++i)
            {
                outFile << "f " << triangle_[i].x + 1 << "/" << tri_tex_[i].x + 1 << "/ "
                    << " " << triangle_[i].y + 1 << "/" << tri_tex_[i].y + 1 << "/ "
                    << " " << triangle_[i].z + 1 << "/" << tri_tex_[i].z + 1 << "/\n";
            }
        }
        else if (hasVertexNormal && !hasTexCoor)
        {
            for (int i = 0; i < numTri; ++i)
            {
                outFile << "f " << triangle_[i].x + 1 << "//" << tri_normal_[i].x + 1 << " "
                    << " " << triangle_[i].y + 1 << "//" << tri_normal_[i].y + 1 << " "
                    << " " << triangle_[i].z + 1 << "//" << tri_normal_[i].z + 1 << " \n";//*/
            }
        }

        else //if (!hasVertexNormal && !hasTexCoor)
        {
            for (int i = 0; i < numTri; ++i)
            {
                outFile << "f " << triangle_[i].x + 1 << " "
                    << " " << triangle_[i].y + 1 << " "
                    << " " << triangle_[i].z + 1 << " \n";
            }
        }

        outFile.close();

        std::cout << "- Mesh exported : " << filename << " (#v=" << numVertex 
            << ",#t=" << numTri << ")" << std::endl;

        return true;
    }

    bool writeMeshInOBJWithMtl(const char* filename, bool hasVertexNormal, bool hasTexCoor, const char* texImgFile)
    {
        std::ofstream outFile(filename);

        if (!outFile.is_open())
        {
            std::cerr << "ERROR! Can't open " << filename << " for writing" << std::endl;
            return false;
        }

        // The name of the shading group used in mtl file
        std::string shadingGroupName = "initialShadingGroup";

        // Find the extension of the filename
        std::string mtlFilePath(filename);
        std::string ext(".obj");
        std::size_t found = mtlFilePath.find(ext);

        // Set the mtl file name
        mtlFilePath.replace(found, ext.length(), std::string(".mtl"));

        found=mtlFilePath.find_last_of("/");
        std::ostringstream mtlFileShortName;
        for (size_t i = found+1; i < mtlFilePath.length(); ++i)
            mtlFileShortName << mtlFilePath[i];
        std::cout << mtlFileShortName.str() << std::endl;

        if (hasTexCoor)
        {
            // Write a mtl file
            std::ofstream mtlOutFile(mtlFilePath.c_str());

            mtlOutFile << "newmtl " << shadingGroupName << "\n"
                << "illum 4\n"
                << "Kd 0.00 0.00 0.00\n"
                << "Ka 0.00 0.00 0.00\n"
                << "Tf 1.00 1.00 1.00\n"
                << "map_Kd " << texImgFile << "\n"//final.jpg
                << "Ni 1.00" << "\n";
            mtlOutFile.close();

            // Add the mtl file info in the OBJ file
            outFile << "mtllib " << mtlFileShortName.str() << "\n";
        }

        int numVertex = (int)vertex_.size();
        for (int i = 0; i < numVertex; ++i)
        {
            outFile << "v " << vertex_[i].x << " " << vertex_[i].y << " " << vertex_[i].z << "\n";
        }

        if (hasVertexNormal)
        {
            int numVNormals = numVertexNormals();
            //std::cout << "numNormals=" << numVNormals << std::endl;

            for (int i = 0; i < numVNormals; ++i)
            {
                outFile << "vn " << vertexNormal_[i].x << " " << vertexNormal_[i].y << " " << vertexNormal_[i].z << "\n";
            }
        }

        if (hasTexCoor)
        {
            int numTex = numTexCoords();
            for (int i = 0; i < numTex; ++i)
            {
                outFile << "vt " << texCoor_[i].x << " " << texCoor_[i].y << "\n";
            }
        }

        int numTri = (int)triangle_.size();
        //std::cout << "numTri=" << numTri << ", " << hasVertexNormal << " " << hasTexCoor << std::endl;

        if (hasVertexNormal && hasTexCoor)
        {
            outFile << "usemtl " << shadingGroupName << "\n";
            for (int i = 0; i < numTri; ++i)
            {
                outFile << "f " << triangle_[i].x + 1 << "/" << tri_tex_[i].x + 1 << "/" << tri_normal_[i].x + 1
                    << " " << triangle_[i].y + 1 << "/" << tri_tex_[i].y + 1 << "/" << tri_normal_[i].y + 1
                    << " " << triangle_[i].z + 1 << "/" << tri_tex_[i].z + 1 << "/" << tri_normal_[i].z + 1 << "\n";
            }
        }
        else if (!hasVertexNormal && hasTexCoor)
        {
            outFile << "usemtl " << shadingGroupName << "\n";

            for (int i = 0; i < numTri; ++i)
            {
                outFile << "f " << triangle_[i].x + 1 << "/" << tri_tex_[i].x + 1 << "/ "
                    << " " << triangle_[i].y + 1 << "/" << tri_tex_[i].y + 1 << "/ "
                    << " " << triangle_[i].z + 1 << "/" << tri_tex_[i].z + 1 << "/\n";
            }
        }
        else if (hasVertexNormal && !hasTexCoor)
        {
            outFile << "usemtl " << shadingGroupName << "\n";
            for (int i = 0; i < numTri; ++i)
            {
                outFile << "f " << triangle_[i].x + 1 << "//" << tri_normal_[i].x + 1 << " "
                    << " " << triangle_[i].y + 1 << "//" << tri_normal_[i].y + 1 << " "
                    << " " << triangle_[i].z + 1 << "//" << tri_normal_[i].z + 1 << " \n";
            }
        }

        else //if (!hasVertexNormal && !hasTexCoor)
        {
            outFile << "usemtl " << shadingGroupName << "\n";

            for (int i = 0; i < numTri; ++i)
            {
                outFile << "f " << triangle_[i].x + 1 << " "
                    << " " << triangle_[i].y + 1 << " "
                    << " " << triangle_[i].z + 1 << " \n";
            }
        }

        outFile.close();

        std::cout << "- Mesh exported : " << filename << " (#v=" << numVertex
            << ",#t=" << numTri << ")" << std::endl;

        return true;
    }

    bool writeMeshInPly(const char* filename, mg::Vector3uc* vertexColor)
    {
        std::ofstream outFile(filename);

        if (!outFile.is_open())
        {
            std::cout << "Error, can't open " << filename << " for writing." << std::endl;
            return false;
        }

        //------------------------------
        // Header
        //------------------------------
        outFile << "ply\n"
            << "format ascii 1.0\n"
            << "comment \n"
            << "comment object : \n";

        outFile << "element vertex " << this->numVertices() << "\n"
            << "property float x\n"
            << "property float y\n"
            << "property float z\n"
            << "property uchar red\n"
            << "property uchar green\n"
            << "property uchar blue\n";

        outFile << "element face " << this->numTriangles() << "\n"
            << "property list uchar int vertex_index\n"
            << "end_header\n";

        //------------------------------
        // Data
        //------------------------------
        for (int i = 0; i < numVertices(); ++i)
        {
            outFile << vertex_[i].x << " " << vertex_[i].y << " " << vertex_[i].z << " "
                << (int)vertexColor[i].z << " " << (int)vertexColor[i].y << " " << (int)vertexColor[i].x << "\n";
        }

        for (int i = 0; i < numTriangles(); ++i)
        {
            outFile << "3 "
                << triangle_[i].x << " "
                << triangle_[i].y << " "
                << triangle_[i].z << "\n";
        }

        outFile.close();

        std::cout << "- Mesh exported : " << filename << " (#v=" << numVertices()
            << ",#t=" << numTriangles() << ")" << std::endl;

        return true;
    }


public:
    std::vector< Vector3<T> > vertex_;
    std::vector< Vector3i >  triangle_;
    std::vector< Vector3i >  tri_tex_;
    std::vector< Vector3i >  tri_normal_;
    std::vector< Vector3<T> > vertexNormal_;
    std::vector< Vector2<T> > texCoor_;
};

template <typename T>
class PolyMesh
{
public:
    void deepCopyTo(PolyMesh<T>& outMesh)
    {
        outMesh.vertex_ = vertex_;
        outMesh.vertex_tex_ = vertex_tex_;
        outMesh.vertex_normal_ = vertex_normal_;

        outMesh.poly_count_ = poly_count_;
        outMesh.poly_index_ = poly_index_;
        outMesh.poly_tex_ = poly_tex_;
        outMesh.poly_normal_ = poly_normal_;
    }

    int numVertices()               { return (int)vertex_.size(); }
    int numVertexNormals()          { return (int)vertex_normal_.size(); }
    int numTexCoords()              { return (int)vertex_tex_.size(); }
    int numPolygons()               { return (int)poly_count_.size(); }

    void get_face_info(std::string face_str,
        std::vector<int>& vid, std::vector<int>& nid, std::vector<int>& tid)
    {
        while (face_str[face_str.size() - 1] == ' ')
            face_str.resize(face_str.size() - 1);

        std::stringstream face_buf;
        face_buf << face_str;

        std::vector<std::string> data;
        std::string temp;

        do
        {
            face_buf >> temp;
            data.push_back(temp);
        } while (!face_buf.eof());

        //std::cout << face_buf.str() << std::endl;

        vid.clear();
        nid.clear();
        tid.clear();


        for (int i = 0; i < (int)data.size(); ++i)
        {
            std::vector<size_t> divider;

            temp = data[i];
            //size_t last_found = -1;

            int count = 0;

            size_t last_found = -1;

            while (temp.length() > 0)
            {
                //std::cout << temp << std::endl;

                size_t found = temp.find('/');
                //std::cout << found << " " << std::string::npos << std::endl;

                if (found != std::string::npos)
                {
                    divider.push_back(found + 1 + last_found);
                    //std::cout << divider.size() << std::endl;

                    temp = temp.substr(found + 1);
                    count++;

                    last_found = found;
                }
                else
                {
                    temp.clear();
                }
            }

            if (divider.size() == 1)
            {
                vid.push_back(atoi(data[i].substr(0, divider[0]).c_str()));
                tid.push_back(atoi(data[i].substr(divider[0] + 1).c_str()));
            }
            else if (divider.size() > 1)
            {
                if (divider[1] - divider[0] > 1)
                {
                    vid.push_back(atoi(data[i].substr(0, divider[0]).c_str()));
                    tid.push_back(atoi(data[i].substr(divider[0] + 1, divider[1] - divider[0] - 1).c_str()));
                    nid.push_back(atoi(data[i].substr(divider[1] + 1).c_str()));
                }
                else
                {
                    vid.push_back(atoi(data[i].substr(0, divider[0]).c_str()));
                    //tid.push_back(atoi(data[i].substr(divider[0] + 1, divider[1] - 1).c_str()));
                    nid.push_back(atoi(data[i].substr(divider[1] + 1).c_str()));
                }
            }

            //for (int j = 0; j < divider.size(); ++j)
            //    std::cout << divider[j] << " ";
            //std::cout << std::endl;

            //while (temp.length() > 0)
            //{
            //    size_t found = temp.find('/');
            //    if (found = std::string::npos)
            //    {
            //        if (count == 0)
            //        {
            //            vid.push_back(atoi(temp.substr(0, found).c_str()));
            //        }
            //        else if (count == 1)
            //        {

            //        }
            //    }
            //    else
            //    {
            //        if (count == 0)
            //        {
            //            vid.push_back(atoi(temp.c_str()));
            //        }

            //        temp = temp.substr(found);
            //    }
            //}// end while(temp.length()>0)
        }// end for(i)
    }

    void extract_data(std::ostringstream& in_buffer)
    {
        // Copy first
        std::stringstream buffer;
        buffer << in_buffer.str();

        //std::cout << buffer.str() << std::endl;

        // Delete arrays
        vertex_.clear();
        vertex_normal_.clear();
        vertex_tex_.clear();

        poly_count_.clear();
        poly_index_.clear();
        poly_tex_.clear();
        poly_normal_.clear();


        std::string temp;
        mg::Vector3<T> v, vn;
        mg::Vector2<T> vt;

        //std::stringstream face_buffer;
        std::stringstream face_buf;
        do
        {
            if (temp == "v")
            {
                buffer >> v.x >> v.y >> v.z;
                std::getline(buffer, temp, '\n');

                vertex_.push_back(v);
            }
            else if (temp == "vn")
            {
                buffer >> vn.x >> vn.y >> vn.z;
                std::getline(buffer, temp, '\n');

                vertex_normal_.push_back(vn);
            }
            else if (temp == "vt")
            {
                buffer >> vt.x >> vt.y;
                std::getline(buffer, temp, '\n');

                vertex_tex_.push_back(vt);
            }
            else if (temp == "f")
            {
                std::getline(buffer, temp, '\n');
                //face_buffer << temp << "\n";

                face_buf << temp;

                std::vector<int> vid, nid, tid;
                get_face_info(temp, vid, nid, tid);

                poly_count_.push_back(vid.size());

                poly_index_.insert(poly_index_.end(), vid.begin(), vid.end());
                poly_tex_.insert(poly_tex_.end(), tid.begin(), tid.end());
                poly_normal_.insert(poly_normal_.end(), nid.begin(), nid.end());

                /*std::cout << "fv : (";
                for (int i = 0; i < vid.size(); ++i)
                {
                std::cout << vid[i] << " ";
                }
                std::cout << ")" << "\n";

                std::cout << "ft : (";
                for (int i = 0; i < tid.size(); ++i)
                {
                std::cout << tid[i] << " ";
                }
                std::cout << ")" << "\n";

                std::cout << "fn : (";
                for (int i = 0; i < nid.size(); ++i)
                {
                std::cout << nid[i] << " ";
                }
                std::cout << ")" << "\n";*/
            }
            buffer >> temp;
            //std::cout << temp << std::endl;
        } while (!buffer.eof());

        /*std::cout
        << "# vertex = " << vertex_.size()
        << ", v[0]=" << vertex_[0] << ", v[" << vertex_.size() - 1 << "]=" << vertex_[vertex_.size() - 1] << "\n"
        << "# vn = " << vertex_normal_.size()
        << ", vn[0]=" << vertex_normal_[0] << ", vn[" << vertex_normal_.size() - 1 << "]=" << vertex_normal_[vertex_normal_.size() - 1] << "\n"
        << "# vt = " << vertex_tex_.size()
        << ", vt[0]=" << vertex_tex_[0] << ", vt[" << vertex_tex_.size() - 1 << "]=" << vertex_tex_[vertex_tex_.size() - 1] << "\n"
        << std::endl;*/

        /*std::cout
        << "poly index = 0 " << "\n"
        << "fv=(" << poly_index_[0] << " " << poly_index_[1] << " " << poly_index_[2] << " " << poly_index_[3] << "\n"
        << "ft=" << poly_tex_[0] << " " << poly_tex_[1] << " " << poly_tex_[2] << " " << poly_tex_[3] << "\n"
        << "fn=" << poly_normal_[0] << " " << poly_normal_[1] << " " << poly_normal_[2] << " " << poly_normal_[3] << "\n"
        << "poly index = " << poly_count_.size() - 1 << "\n"
        << "fv=(" << poly_index_[poly_index_.size() - 4] << " " << poly_index_[poly_index_.size() - 3] << " "
        << poly_index_[poly_index_.size() - 2] << " " << poly_index_[poly_index_.size() - 1] << "\n"
        << "ft=" << poly_tex_[poly_tex_.size() - 4] << " " << poly_tex_[poly_tex_.size() - 3] << " "
        << poly_tex_[poly_tex_.size() - 2] << " " << poly_tex_[poly_tex_.size() - 1] << "\n"
        << "fn=" << poly_normal_[poly_normal_.size() - 4] << " " << poly_normal_[poly_normal_.size() - 3] << " "
        << poly_normal_[poly_normal_.size() - 2] << " " << poly_normal_[poly_normal_.size() - 1] << "\n"
        << std::endl;*/

        //std::cout << face_buffer.str() << std::endl;
    }

    bool loadMeshInOBJ(const char* filename)
    {
        std::ifstream inFile(filename);

        if (!inFile.is_open())
        {
            std::cout << "Error, can't open mesh : " << filename << std::endl;
            return false;
        }

        std::ostringstream buffer;
        buffer << inFile.rdbuf();

        inFile.close();

        // Extract vertex, vertex normal and texture coordinates
        extract_data(buffer);

        return true;
    }

    bool writeMeshInOBJ(const char* filename, bool hasVertexNormal, bool hasTexCoor)
    {
        std::ofstream outFile(filename);

        if (!outFile.is_open())
        {
            std::cout << "Error, can't open " << filename << " for writing" << std::endl;
            return false;
        }

        for (int i = 0; i < (int)vertex_.size(); ++i)
        {
            outFile << "v " << vertex_[i].x << " " << vertex_[i].y << " " << vertex_[i].z << "\n";
        }

        if (hasTexCoor)
        {
            for (int i = 0; i < (int)vertex_tex_.size(); ++i)
            {
                outFile << "vt " << vertex_tex_[i].x << " " << vertex_tex_[i].y << "\n";
            }
        }

        if (hasVertexNormal)
        {
            for (int i = 0; i < (int)vertex_normal_.size(); ++i)
            {
                outFile << "vn " << vertex_normal_[i].x << " " << vertex_normal_[i].y << " " << vertex_normal_[i].z << "\n";
            }
        }

        if (!hasVertexNormal && !hasTexCoor)
        {
            int fid = 0;

            for (int i = 0; i < (int)poly_count_.size(); ++i)
            {
                int nv = poly_count_[i];

                outFile << "f ";
                for (int j = 0; j < nv; ++j)
                    outFile << poly_index_[fid + j] << " ";
                outFile << "\n";

                fid += nv;
            }
        }
        else if (!hasVertexNormal && hasTexCoor)
        {
            int fid = 0;

            for (int i = 0; i < (int)poly_count_.size(); ++i)
            {
                int nv = poly_count_[i];

                outFile << "f ";
                for (int j = 0; j < nv; ++j)
                    outFile << poly_index_[fid + j] << "/" << poly_tex_[fid + j] << " ";
                outFile << "\n";

                fid += nv;
            }
        }
        else if (hasVertexNormal && !hasTexCoor)
        {
            int fid = 0;

            for (int i = 0; i < (int)poly_count_.size(); ++i)
            {
                int nv = poly_count_[i];

                outFile << "f ";
                for (int j = 0; j < nv; ++j)
                    outFile << poly_index_[fid + j] << "//" << poly_normal_[fid + j] << " ";
                outFile << "\n";

                fid += nv;
            }
        }
        else if (hasVertexNormal && hasTexCoor)
        {
            int fid = 0;

            for (int i = 0; i < (int)poly_count_.size(); ++i)
            {
                int nv = poly_count_[i];

                outFile << "f ";
                for (int j = 0; j < nv; ++j)
                    outFile << poly_index_[fid + j] << "/" << poly_tex_[fid + j] << "/" << poly_normal_[fid + j] << " ";
                outFile << "\n";

                fid += nv;
            }
        }

        outFile.close();

        std::cout << "- Mesh exported : " << filename << " (#v=" << vertex_.size()
            << ",#f=" << poly_index_.size() << ", #t=" << vertex_tex_.size() << ")" << std::endl;

        return true;
    }

    /*bool writeMeshInOBJWithMtl(const char* filename, bool hasVertexNormal, bool hasTexCoor, const char* texImgFile)
    {}//*/

    bool writeMeshInPly(const char* filename, mg::Vector3uc* vertexColor)
    {
        std::ofstream outFile(filename);

        if (!outFile.is_open())
        {
            std::cout << "Error, can't open " << filename << " for writing." << std::endl;
            return false;
        }

        //------------------------------
        // Header
        //------------------------------
        outFile << "ply\n"
            << "format ascii 1.0\n"
            << "comment \n"
            << "comment object : \n";

        outFile << "element vertex " << this->numVertices() << "\n"
            << "property float x\n"
            << "property float y\n"
            << "property float z\n"
            << "property uchar red\n"
            << "property uchar green\n"
            << "property uchar blue\n";

        outFile << "element face " << this->numPolygons() << "\n"
            << "property list uchar int vertex_index\n"
            << "end_header\n";

        //------------------------------
        // Data
        //------------------------------
        for (int i = 0; i < numVertices(); ++i)
        {
            outFile << vertex_[i].x << " " << vertex_[i].y << " " << vertex_[i].z << " "
                << (int)vertexColor[i].z << " " << (int)vertexColor[i].y << " " << (int)vertexColor[i].x << "\n";
        }

        int vid = 0;
        for (int i = 0; i < this->numPolygons(); ++i)
        {
            /*outFile << "3 "
                << triangle_[i].x << " "
                << triangle_[i].y << " "
                << triangle_[i].z << "\n";*/

            outFile << poly_count_[i] << " ";
            for (int j = 0; j < poly_count_[i]; ++j)
            {
                outFile << poly_index_[vid + j] << " ";
            }
            outFile << "\n";

            vid += poly_count_[i];
        }

        outFile.close();

        std::cout << "- Mesh exported : " << filename << " (#v=" << numVertices()
            << ",#f=" << poly_index_.size() << ", #t=" << vertex_tex_.size() << ")" << std::endl;

        return true;
    }
public:
    std::vector< mg::Vector3<T> > vertex_;
    std::vector< mg::Vector3<T> > vertex_normal_;
    std::vector< mg::Vector2<T> > vertex_tex_;

    std::vector< int >  poly_count_;
    std::vector< int >  poly_index_;
    std::vector< int >  poly_tex_;
    std::vector< int >  poly_normal_;
};

typedef TriMesh<float> TriMesh_f;
typedef PolyMesh<float> PolyMesh_f;

}