#include <iostream>

#include "marching_cubes.h"

int main(int argc, char** argv)
{
	if (argc <= 5)
	{
		std::cout << "- Usage : " << argv[0] << " arg1  arg2 arg3 arg4 arg5" << "\n"
			<< " . arg1 : tsdf file name " << "\n"
			<< " . arg2 : output mesh file name" << "\n"
			<< " . arg3 : iso-value to extract surface with" << "\n"
			<< " . arg4 : cut-value" << "\n"
			<< " . arg5 : computing option(0: CPU, 1: GPU)" << "\n"
			<< " . arg6 : img_xmin" << "\n"
			<< " . arg7 : img_ymin" << "\n"
			<< " . arg8 : img_dx" << "\n"
			<< " . arg9 : img_dy" << "\n"
			<< " . arg10 : img_nx" << "\n"
			<< " . arg11 : img_ny" << "\n"
			<< std::endl;

		return 0;
	}

	const char* tsdf_file = argv[1];
	const char* outmesh_file = argv[2];
	float isoValue = atof(argv[3]);
	float cutValue = atof(argv[4]);
	int useCUDA = atoi(argv[5]);

	float img_xmin(0.0f), img_ymin(0.0f);
	float img_dx(1.0f), img_dy(1.0f);
	int img_nx(1), img_ny(1);

	if (argc > 6)
	{
		img_xmin = atof(argv[6]);
		img_ymin = atof(argv[7]);
		img_dx   = atof(argv[8]);
		img_dy   = atof(argv[9]);
		img_nx   = atof(argv[10]);
		img_ny   = atof(argv[11]);
	}

	MarchingCubes mc;
	mc.loadVolumeData(tsdf_file);

	mc.img_xmin_ = img_xmin;
	mc.img_ymin_ = img_ymin;
	mc.img_dx_   = img_dx  ;
	mc.img_dy_   = img_dy  ;
	mc.img_nx_   = img_nx  ;
	mc.img_ny_   = img_ny  ;

	if (argc > 6) mc.genTexture_ = true;

	if(useCUDA)
		mc.generateSurfaceMesh_cuda(isoValue, cutValue, outmesh_file);
	else
	{
		mc.generateSurfaceMesh(isoValue, cutValue);
		mc.exportMeshInObj(outmesh_file);
	}

	return 0;
}