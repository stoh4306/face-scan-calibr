#pragma once
#include<opencv2/opencv.hpp>
#include <iostream>
#include <string>

#include<opencv2/opencv.hpp>

typedef char    byte;
typedef int     int32;
typedef short   int16;
typedef unsigned char uByte;
typedef unsigned int uInt32;
typedef unsigned short uInt16;
typedef float Float32;
typedef unsigned long ULONG;

#define RGB_FLAG     (1)
#define ALPHA_FLAG   (2)
#define ZBUFFER_FLAG (4)

#define CHUNK_STACK_SIZE (32)


/* Error code definitions */
#define IFF_NO_ERROR     (0)
#define IFF_OPEN_FAILS   (1)
#define IFF_READ_FAILS   (2)
#define IFF_BAD_TAG      (3)
#define IFF_BAD_COMPRESS (4)
#define IFF_BAD_STACK    (5)
#define IFF_BAD_CHUNK    (6)


/* Define the IFF tags we are looking for in the file. */
const uInt32 IFF_TAG_CIMG = ('C' << 24) | ('I' << 16) | ('M' << 8) | ('G');
const uInt32 IFF_TAG_FOR4 = ('F' << 24) | ('O' << 16) | ('R' << 8) | ('4');
const uInt32 IFF_TAG_TBHD = ('T' << 24) | ('B' << 16) | ('H' << 8) | ('D');
const uInt32 IFF_TAG_TBMP = ('T' << 24) | ('B' << 16) | ('M' << 8) | ('P');
const uInt32 IFF_TAG_RGBA = ('R' << 24) | ('G' << 16) | ('B' << 8) | ('A');
const uInt32 IFF_TAG_CLPZ = ('C' << 24) | ('L' << 16) | ('P' << 8) | ('Z');
const uInt32 IFF_TAG_ESXY = ('E' << 24) | ('S' << 16) | ('X' << 8) | ('Y');
const uInt32 IFF_TAG_ZBUF = ('Z' << 24) | ('B' << 16) | ('U' << 8) | ('F');
const uInt32 IFF_TAG_BLUR = ('B' << 24) | ('L' << 16) | ('U' << 8) | ('R');
const uInt32 IFF_TAG_BLRT = ('B' << 24) | ('L' << 16) | ('R' << 8) | ('T');
const uInt32 IFF_TAG_HIST = ('H' << 24) | ('I' << 16) | ('S' << 8) | ('T');


/* For the stack of chunks */
typedef struct _iff_chunk {
	uInt32 tag;
	uInt32 start;
	uInt32 size;
	uInt32 chunkType;
} iff_chunk;



struct IffIoFileWrapperData
{
	const char *szFilePath;
	FILE *pFile;
};


// Leo 28/Aug/2008: Moved data out of global variables to enable concurrent instances use.
typedef struct _iff_instance {
	iff_chunk chunkStack[CHUNK_STACK_SIZE];
	int chunkDepth;
	/* The current error state. */
	unsigned int iff_error;
	/* The file callbacks */
	//iff_file_io *pFileIO;
	IffIoFileWrapperData *pData;
} iff_instance;
typedef struct _iff_image {
	unsigned width;   //!< Image width.
	unsigned height;  //!< Image height.
	unsigned depth;   //!< Color format, bytes per color pixel.
	unsigned char * rgba;   //!< The color data of size width*height*depth.
	float znear;  //!< The near clipping plane for the z buffer data.
	float zfar;   //!< The far clipping plane for the z buffer data.
	/*!
	Access pixel x,y as zbuffer[width*y + x]. (Starting from the top left.)

	The stored values now are -1/z components in light eye space. When
	reading a z depth value back, one thus need to compute (-1.0 / z)
	to get the real z component in eye space. This format is the same as
	the camera depth format. This format is always used, whatever the light
	type is. In order to be able to compute real 3D distances of the shadow
	samples there's an IFF tag: 'ESXY' 'Eye Space X Y' values are two float
	values, the first one is the width size, the second one is the height size
	matching the map width and map height. When reading the shadow map buffer,
	one can thus convert from the pixel coordinates to the light eye space
	xy coords. If the pixel coordinates are considered to be normalized in
	the [-1.0, 1.0] range, the eye space xy are given by:

	X_light_eye_space = 3D normalized_pixel_x * (width / 2.0)

	Y_light_eye_space = 3D normalized_pixel_y * (heigth / 2.0)

	Once one get the (X,Y,Z) light eye space coordinates, the true 3D
	distance from the light source is:

	true_distance = 3D sqrt(X=B2 + Y=B2 + Z=B2)

	(This was copied from the Alias|Wavefront API Knowledge Base.)
	*/
	float * zbuffer; //!< The z buffer data of size width*height.
	/*!
	Eye x-ratio. This is a depth map specific field used to
	compute the xy eye coordinates from the normalized pixel.
	*/
	float zesx;
	/*!
	Eye y-ratio. This is a depth map specific field used to
	compute the xy eye coordinates from the normalized pixel.
	*/
	float zesy;
	/*!
	This does not work right now! I do not know how to interpret these vectors.
	If you can figure it out, please let me know.
	*/
	float * blurvec; //!< The packed xy motion blur vectors of size 2*width*height.
} iff_image;

int iff_begin_read_chunk(iff_instance *pInstance, iff_chunk *pOutChunk);
int iff_end_read_chunk(iff_instance *pInstance);
//uByte * iff_read_data( iff_instance *pInstance, int size );
std::vector<uByte> iff_read_data(iff_instance *pInstance, int size);
uByte * iff_decompress_rle(iff_instance *pInstance, uInt32 numBytes, uByte * compressedData, uInt32 compressedDataSize, uInt32 * compressedIndex);
std::vector<uByte> iff_decompress_tile_rle(iff_instance *pInstance, uInt16 width, uInt16 height, uInt16 depth, uByte* compressedData, uInt32 compressedDataSize);
uByte * iff_read_uncompressed_tile(iff_instance *pInstance, uInt16 width, uInt16 height, uInt16 depth);

void loadIff(std::string fileName);



int IffOpen(IffIoFileWrapperData *pData);
int IffClose(IffIoFileWrapperData *pData);
size_t IffRead(IffIoFileWrapperData *pRealData, void *pBuffer, size_t elementSize, size_t elementCount);
long IffSeek(IffIoFileWrapperData *pRealData, long offset, int origin);
long IffTell(IffIoFileWrapperData *pRealData);
int iff_loadA(IffIoFileWrapperData* _fileData);
void iff_free(iff_image * image);

int iff_get_short(iff_instance *pInstance, uInt16 *pOut);
int iff_get_long(iff_instance *pInstance, uInt32 *pOut);
int iff_get_float(iff_instance *pInstance, Float32 *pOut);
int iff_begin_read_chunk(iff_instance *pInstance, iff_chunk *pOutChunk);
int iff_end_read_chunk(iff_instance *pInstance);
std::vector<uByte> iff_read_data(iff_instance *pInstance, int size);



