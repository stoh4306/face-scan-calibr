
#include"renew_iff.h"

#ifdef __unix
int fopen_s(FILE** pFile, const char* filename, const char* mode)
{
	*pFile = fopen(filename, mode);
	if(*pFile!=NULL)
		return 0;
	else
		return 1;
}
#endif

inline void ZeroMemory(void* ptr, size_t size) {
	memset(ptr, 0, size);
}

int IffOpen(IffIoFileWrapperData *pData)
{
	int error=0;
	if (!pData->pFile)
	{
		if (error != fopen_s(&pData->pFile, pData->szFilePath, "rb"))
		{
			pData->pFile = NULL;
		}
	}

	return (pData->pFile != NULL);
}

int IffClose(IffIoFileWrapperData *pData)
{

	if (pData->pFile)
	{
		fclose(pData->pFile);
		pData->pFile = NULL;
	}

	return true;
}

size_t IffRead(IffIoFileWrapperData *pRealData, void *pBuffer, size_t elementSize, size_t elementCount)
{

	if (!pRealData->pFile)
	{
		ZeroMemory(pBuffer, elementSize * elementCount);
		return 0;
	}

	size_t res = fread(pBuffer, elementSize, elementCount, pRealData->pFile);

	if (res == 0)
	{
		ZeroMemory(pBuffer, elementSize * elementCount);
	}
	else if (res < elementCount)
	{
		ZeroMemory(reinterpret_cast<byte *>(pBuffer) + res * elementSize, (elementCount - res) * elementSize);
	}

	return res;
}

long IffSeek(IffIoFileWrapperData *pRealData, long offset, int origin)
{
	if (!pRealData->pFile)
	{
		return -1;
	}

	return fseek(pRealData->pFile, offset, origin);
}

long IffTell(IffIoFileWrapperData *pRealData)
{
	if (!pRealData->pFile)
	{
		return -1L;
	}

	return ftell(pRealData->pFile);
}

void loadIff(std::string fileName) {
	IffIoFileWrapperData ioData = { 0 };

	ioData.szFilePath = fileName.c_str();
	ioData.pFile = NULL;

	int result = iff_loadA(&ioData);
}

//iff_image * iff_loadA(iff_file_io *pFileIO, unsigned int *pOutIffError, int *pIsOriginal16Bit, int fWantRGBA, int fWantZBuffer, int fWantBlurVec, int fWantExtraInfo)
int iff_loadA(IffIoFileWrapperData* _fileData)
{
	int fWantRGBA = true;
	int fWantZBuffer = true;
	int fWantBlurVec = false;
	int fWantExtraInfo = false;
	//char* filenamechar = ConvertWCtoC(((IffIoFileWrapperData*)_fileData)->szFilePath);
	std::string fileName = _fileData->szFilePath;
	iff_instance *pInstance;
	iff_chunk chunkInfo;
	iff_image * image;

	// -- Header info.
	uInt32 width, height, depth, npixels;
	uInt32 flags, compress;
	uInt16 tiles;
	uInt16 tbmpfound;

	uInt16 x1, x2, y1, y2, tile_width, tile_height;
	uInt32 tile;
	uInt32 ztile;

	uInt32 iscompressed, i;
	uByte *tileData = nullptr;
	uInt32 remainingDataSize;

	long oldSpot;
	long tellTemp;
	uInt32 fileLength;
	cv::Mat mat;
	cv::Mat mat2;
	cv::Mat colorMat;

	//if (pIsOriginal16Bit) {
	//	*pIsOriginal16Bit = FALSE;
	//}

	//if (!pFileIO) {
	//	if (pOutIffError) {
	//		*pOutIffError = IFF_OPEN_FAILS;
	//	}
	//	return(0);
	//}

	pInstance = (iff_instance*)malloc(sizeof(iff_instance));

	if (!pInstance)
	{
		//if (pOutIffError) {
		//	*pOutIffError = IFF_OPEN_FAILS;
		//}
		return(0);
	}

	memset(pInstance, 0, sizeof(iff_instance));
	pInstance->iff_error = IFF_NO_ERROR;
	pInstance->pData = _fileData;

	// -- Initialize the top of the chunk stack.
	pInstance->chunkDepth = -1;
	image = 0;

	// -- Open the file.    
	if (!IffOpen(_fileData)) {
		if (pInstance->iff_error == IFF_NO_ERROR) {
			pInstance->iff_error = IFF_OPEN_FAILS;
		}
		goto CleanUp;
	}

	// -- File should begin with a FOR4 chunk of type CIMG
	if (!iff_begin_read_chunk(pInstance, &chunkInfo)) {
		if (pInstance->iff_error == IFF_NO_ERROR) {
			pInstance->iff_error = IFF_READ_FAILS;
		}
		goto CleanUp;
	}

	if (chunkInfo.chunkType != IFF_TAG_CIMG) {
		if (pInstance->iff_error == IFF_NO_ERROR) {
			// -- This is not a CIMG, it is not an IFF Image.
			pInstance->iff_error = IFF_BAD_TAG;
		}
		goto CleanUp;
	}

	/*
	* Read the image header
	* OK, we have a FOR4 of type CIMG, look for the following tags
	*		FVER
	*		TBHD	bitmap header, definition of size, etc.
	*		AUTH
	*		DATE
	*/
	while (1) {

		if (!iff_begin_read_chunk(pInstance, &chunkInfo)) {
			if (pInstance->iff_error == IFF_NO_ERROR) {
				pInstance->iff_error = IFF_READ_FAILS;
			}
			goto CleanUp;
		}

		// -- Right now, the only info we need about the image is in TBHD
		// -- so search this level until we find it.
		if (chunkInfo.tag == IFF_TAG_TBHD) {

			uInt16 bytes = 0;

			// -- Header chunk found
			if (!iff_get_long(pInstance, &width)
				|| !iff_get_long(pInstance, &height)
				|| !iff_get_short(pInstance, 0) // -- Don't support 
				|| !iff_get_short(pInstance, 0) // -- Don't support 
				|| !iff_get_long(pInstance, &flags)
				|| !iff_get_short(pInstance, &bytes) // 0 => 8-bit samples, 1 => 16-bit samples. (Our output is 8-bit regardless.)
				|| !iff_get_short(pInstance, &tiles)
				|| !iff_get_long(pInstance, &compress)) {
				if (pInstance->iff_error == IFF_NO_ERROR) {
					pInstance->iff_error = IFF_READ_FAILS;
				}
				goto CleanUp;
			}

			//if (pIsOriginal16Bit && bytes == 1) {
			//	*pIsOriginal16Bit = true;
			//}

#ifdef __IFF_DEBUG_			
			printf("****************************************\n");
			printf("Width: %u\n", width);
			printf("Height: %u\n", height);
			printf("flags: 0x%X\n", flags);
			printf("tiles: %hu\n", tiles);
			printf("compress: %u\n", compress);
			printf("****************************************\n");
#endif

			if (!iff_end_read_chunk(pInstance)) {
				if (pInstance->iff_error == IFF_NO_ERROR) {
					pInstance->iff_error = IFF_READ_FAILS;
				}
				goto CleanUp;
			}

			if (compress > 1) {
				if (pInstance->iff_error == IFF_NO_ERROR) {
					pInstance->iff_error = IFF_BAD_COMPRESS;
				}
				goto CleanUp;
			}

			break;
		}
		else {

#ifdef __IFF_DEBUG_
			// Skipping unused data at FOR4 <size> CIMG depth
			printf("Skipping Chunk: %c%c%c%c\n",
				(((&pInstance->chunkStack[pInstance->chunkDepth].tag)[0] >> 24) & 0xFF),
				(((&pInstance->chunkStack[pInstance->chunkDepth].tag)[0] >> 16) & 0xFF),
				(((&pInstance->chunkStack[pInstance->chunkDepth].tag)[0] >> 8) & 0xFF),
				(((&pInstance->chunkStack[pInstance->chunkDepth].tag)[0] >> 0) & 0xFF));
#endif

			if (!iff_end_read_chunk(pInstance)) {
				if (pInstance->iff_error == IFF_NO_ERROR) {
					pInstance->iff_error = IFF_READ_FAILS;
				}
				goto CleanUp;
			}
		}
	} /* END find TBHD while loop */


	// -- Number of channels.
	depth = 0;

	if (flags & RGB_FLAG) {
		depth += 3;
	}

	if (flags & ALPHA_FLAG) {
		depth += 1;
	}

	if (depth == 0) {
		depth = 1; // If RGB_FLAG isn't set then it's greyscale. (If ALPHA_FLAG is set and RGB_FLAG isn't then I'm not sure what that means. We want to load 1, 3 or 4 channels only.
	}


	npixels = width * height;


	// -- Allocate the image struct.
	image = (iff_image *)malloc(sizeof(iff_image));

	if (!image) {
		if (pInstance->iff_error == IFF_NO_ERROR) {
			pInstance->iff_error = IFF_READ_FAILS;
		}
		goto CleanUp;
	}

	image->width = width;
	image->height = height;
	image->depth = depth;
	image->znear = 0.0;
	image->zfar = 0.0;
	image->zesx = 0.0;
	image->zesy = 0.0;
	image->rgba = 0;
	image->zbuffer = 0;
	image->blurvec = 0;

	// Leo Davidson 27/Aug/2008: If only the header is wanted then we're done.
	//if (!fWantRGBA && !fWantZBuffer && !fWantBlurVec && !fWantExtraInfo)
	//{
	//	goto CleanUp;
	//}

	if (fWantRGBA)
	//if(false)
	{
		image->rgba = (uByte *)malloc(npixels * depth * sizeof(uByte));

		if (!(image->rgba)) {
			if (pInstance->iff_error == IFF_NO_ERROR) {
				pInstance->iff_error = IFF_READ_FAILS;
			}
			goto CleanUp;
		}

		memset(image->rgba, 0, npixels * depth * sizeof(uByte));
	}

	//if (fWantZBuffer && (flags & ZBUFFER_FLAG))
	if((flags & ZBUFFER_FLAG))
	{
		image->zbuffer = (Float32 *)malloc(npixels * sizeof(Float32));

		if (!(image->zbuffer)) {
			if (pInstance->iff_error == IFF_NO_ERROR) {
				pInstance->iff_error = IFF_READ_FAILS;
			}
			goto CleanUp;
		}

		memset(image->zbuffer, 0, npixels * sizeof(Float32));
	}

	// -- Assume the next FOR4 of type TBMP
	tbmpfound = 0;

	//mat = cv::Mat(height, width, CV_32F);
	mat = cv::Mat(height, width, CV_8UC4);
	mat2 = cv::Mat(height, width, CV_8U);
	colorMat = cv::Mat(height, width, CV_8UC3);
	// Read the tiled image data
	while (!tbmpfound) {

		uInt16 tile_area;

		if (!iff_begin_read_chunk(pInstance, &chunkInfo)) {
			if (pInstance->iff_error == IFF_NO_ERROR) {
				pInstance->iff_error = IFF_READ_FAILS;
			}
			goto CleanUp;
		}

		/*
		* OK, we have a FOR4 of type TBMP, (embedded FOR4)
		* look for the following tags
		*		RGBA	color data,	RLE compressed tiles of 32 bbp data
		*		ZBUF	z-buffer data, 32 bit float values
		*		CLPZ	depth map specific, clipping planes, 2 float values
		*		ESXY	depth map specific, eye x-y ratios, 2 float values
		*		HIST
		*		VERS
		*		FOR4 <size>	BLUR (twice embedded FOR4)
		*/
		if (chunkInfo.chunkType == IFF_TAG_TBMP) {
			tbmpfound = 1;

			// Image data found
			tile = 0;
			ztile = 0;


#ifdef __IFF_DEBUG_			
			printf("Reading image tiles\n");
#endif


			if (!(flags & ZBUFFER_FLAG)) {
				ztile = tiles;
			}

			if (depth == 0) {
				tile = tiles;
			}


			// -- Read tiles
			while ((tile < tiles) || (ztile < tiles)) {

				if (!iff_begin_read_chunk(pInstance, &chunkInfo)) {
					if (pInstance->iff_error == IFF_NO_ERROR) {
						pInstance->iff_error = IFF_READ_FAILS;
					}
					goto CleanUp;
				}

				if (!(chunkInfo.tag == IFF_TAG_RGBA) && !(chunkInfo.tag == IFF_TAG_ZBUF)) {
					if (pInstance->iff_error == IFF_NO_ERROR) {
						pInstance->iff_error = IFF_BAD_CHUNK;
					}
					goto CleanUp;
				}

				// Get tile size and location info
				if (!iff_get_short(pInstance, &x1)
					|| !iff_get_short(pInstance, &y1)
					|| !iff_get_short(pInstance, &x2)
					|| !iff_get_short(pInstance, &y2)) {
					if (pInstance->iff_error == IFF_NO_ERROR) {
						pInstance->iff_error = IFF_READ_FAILS;
					}
					goto CleanUp;
				}

				remainingDataSize = chunkInfo.size - 8;

				tile_width = x2 - x1 + 1;
				tile_height = y2 - y1 + 1;
				tile_area = tile_width * tile_height;

#ifdef __IFF_DEBUG_
				printf("Tile x1: %hu  ", x1);
				printf("y1: %hu  ", y1);
				printf("x2: %hu  ", x2);
				printf("y2: %hu\n", y2);
#endif

				// Leo Davidson 27/Aug/2008: Ensure the tile is within the bounds of the image, else we'll write over some other memory.
				if (x1 > x2 || ((unsigned int)(x1 + 1)) >= image->width
					|| y1 > y2 || ((unsigned int)(y1 + 1)) >= image->height) {
					if (pInstance->iff_error == IFF_NO_ERROR) {
						pInstance->iff_error = IFF_BAD_CHUNK;
					}
					goto CleanUp;
				}

				// Leo 27/Aug/2008: This used to assume the data was uncompressed if its size was >= the uncompressed size.
				// I've changed it to require the size to be == the uncompressed size after being given a test image which
				// decoded incorrectly because of compressed tiles that were actually larger than the uncompressed size.
				// Seems rather odd for the encoder to do that, and Photoshop's Maya IFF plugin wouldn't load the image at all,
				// but making this change makes everything work (it seems). Presumably if encoders can generate compressed
				// tiles that are larger than the uncompressed size then they can also generate ones which are exactly equal
				// to it. Is there a way to detect whether a tile is compressed other than by its size?
				if ((int)chunkInfo.size == (tile_width * tile_height * depth + 8)) {
					// -- Compression was not used for this tile.
					iscompressed = 0;
				}
				else {
					iscompressed = 1;
				}

				// -- OK, we found an RGBA chunk, eat it.
				if (chunkInfo.tag == IFF_TAG_RGBA) {

					// Leo Davidson 27/Aug/2008: Increment tile counter if the depth is zero and we skip the chunk, else we'll never stop looping (unless I've missed something?)
					// Leo Davidson 27/Aug/2008: Skip the chunk if we don't want RGBA.
					if (depth == 0 || !fWantRGBA)
					{
						if (!iff_end_read_chunk(pInstance)) {
							if (pInstance->iff_error == IFF_NO_ERROR) {
								pInstance->iff_error = IFF_READ_FAILS;
							}
							goto CleanUp;
						}
						tile++;
						continue;
					}
					std::vector<uByte> tileBuf;
					tileData = 0;

					if (iscompressed) {
						std::vector<uByte> data = iff_read_data(pInstance, remainingDataSize);
						if (data.size() == remainingDataSize) {
							tileBuf = iff_decompress_tile_rle(pInstance, tile_width, tile_height, depth, &data[0], remainingDataSize);
							tileData = &tileBuf[0];
							//free( data );
						}

						//uByte * data = iff_read_data(pInstance, remainingDataSize);
						//if (data) {
						//	tileData = iff_decompress_tile_rle(pInstance, tile_width, tile_height, depth, data, remainingDataSize);
						//	free(data);
						//}
					}
					else {
						tileData = iff_read_uncompressed_tile(pInstance, tile_width, tile_height, depth);
					}

					if (!tileData) {
						if (pInstance->iff_error == IFF_NO_ERROR) {
							pInstance->iff_error = IFF_READ_FAILS;
						}
						goto CleanUp;
					}

					// Dummy block for stupic C variable scope rules.
					{
						/* Dump RGBA data to our data structure */
						uInt16 i;

						// -- base is the tile offset into the packed array
						uInt32 base = image->depth*(image->width*y1 + x1);

						for (i = 0; i < tile_height; i++) {
							memcpy(&image->rgba[base + image->depth*i*width],
								&tileData[depth*i*tile_width],
								tile_width*depth * sizeof(uByte));
						}
					}

					/* End RGBA dump */
					if (!iscompressed) {
						free(tileData);

					}
					tileData = 0;

					if (!iff_end_read_chunk(pInstance)) {
						if (pInstance->iff_error == IFF_NO_ERROR) {
							pInstance->iff_error = IFF_READ_FAILS;
						}
						goto CleanUp;
					}
					tile++;
				} /* END RGBA chunk */

				// -- OK, we found a ZBUF chunk, eat it....hmmm, tasty
				else if (chunkInfo.tag == IFF_TAG_ZBUF) {

					// Leo Davidson 27/Aug/2008: Skip the chunk if we don't want RGBA.
					//if (!fWantZBuffer)
					if(false)
					{
						if (!iff_end_read_chunk(pInstance)) {
							if (pInstance->iff_error == IFF_NO_ERROR) {
								pInstance->iff_error = IFF_READ_FAILS;
							}
							goto CleanUp;
						}
						ztile++;
						continue;
					}
					std::vector<uByte> tileBuf;
					tileData = 0;

					// -- Read in the tile data.
					if (iscompressed) {
						//uByte * data = iff_read_data(pInstance, remainingDataSize);
						//if (data) {
						//	tileData = iff_decompress_tile_rle(pInstance, tile_width, tile_height, 4, data, remainingDataSize);
						//	free(data);
						//}
						std::vector<uByte> data = iff_read_data(pInstance, remainingDataSize);
						if (data.size() == remainingDataSize) {
							tileBuf = iff_decompress_tile_rle(pInstance, tile_width, tile_height, 4, &data[0], remainingDataSize);
							tileData = &tileBuf[0];
							//free( data );
						}
					}
					else {
						tileData = iff_read_uncompressed_tile(pInstance, tile_width, tile_height, 4);
					}

					if (!tileData) {
						if (pInstance->iff_error == IFF_NO_ERROR) {
							pInstance->iff_error = IFF_READ_FAILS;
						}
						goto CleanUp;
					}

					/*
					* Dump DEPTH data into our structure of floats
					*/


					// Dummy block for stupic C variable scope rules.
					{
						int i, j, base;
						uInt32 value, index;
						//uint16_t* matPtr = (uint16_t*)mat.data;
						//uint8_t* mat2Ptr = (uint8_t*)mat2.data;
						base = y1 * width + x1; // Leo 27/Aug/2008: This was being recalculated for every point. Moved it out of the loop.

						for (i = 0; i < tile_height; i++) {
							for (j = 0; j < tile_width; j++) {

								index = 4 * (i*tile_width + j);

								/*
								* If MSB of float is 1, it is negative
								* as it should be, so let's process it
								* if it's not 1, just set the depth to 0
								*
								* NOTE: depth values are stored as
								* -1/z in light eye space in the IFF file.
								*/

								//if( tileData[index+3] & 0x80 ) {

								value = (tileData[index + 3] << 24) +
									(tileData[index + 2] << 16) +
									(tileData[index + 1] << 8) +
									(tileData[index] << 0);
								//if (value != 0)
								//	printf("1");

								image->zbuffer[base + i * width + j] = *((Float32 *)&value);

								auto znear = image->znear;
								auto zfar = image->zfar;

								//if (value != 0)
								//{
								//	int xx = x1 + j;
								//	int yy = y1 + i;
								//	float fval = -(*((Float32 *)&value));
								//	//uint16_t val = (fval * -10000);
								//	//uint16_t val = ((fval - znear) / (zfar - znear)) * 65536;
								//	//uchar uval = ((fval - znear) / (zfar - znear)) * 256;
								//	mat2Ptr[base + i * width + j] = val;
								//	//
								//	//mat.at<cv::uint16_t>(yy, xx) = val;
								//	matPtr[base + i * width + j] = uval;

								//}
								//( -1.0f * (*( (Float32 *)&value ) ));

								//} else {
								//image->zbuffer[base + i*width+j] = 0.0;
								//}

							} /* End DEPTH dump */
						}
					}
					//255,159
					if (!iscompressed) {
						free(tileData);
					}
					tileData = 0;

					if (!iff_end_read_chunk(pInstance)) {
						if (pInstance->iff_error == IFF_NO_ERROR) {
							pInstance->iff_error = IFF_READ_FAILS;
						}
						goto CleanUp;
					}
					ztile++;

				} /* END ZBUF chunk */

			} /* END while TBMP tiles */

		} /* END if TBMP */

		else {

#ifdef __IFF_DEBUG_
			// Skipping unused data IN THE BEGINNING OF THE FILE
			printf("Skipping Chunk in search of TBMP: %c%c%c%c\n",
				(((&pInstance->chunkStack[pInstance->chunkDepth].tag)[0] >> 24) & 0xFF),
				(((&pInstance->chunkStack[pInstance->chunkDepth].tag)[0] >> 16) & 0xFF),
				(((&pInstance->chunkStack[pInstance->chunkDepth].tag)[0] >> 8) & 0xFF),
				(((&pInstance->chunkStack[pInstance->chunkDepth].tag)[0] >> 0) & 0xFF));
#endif

			if (!iff_end_read_chunk(pInstance)) {
				if (pInstance->iff_error == IFF_NO_ERROR) {
					pInstance->iff_error = IFF_READ_FAILS;
				}
				goto CleanUp;
			}
		}
	}

	// Leo Davidson 27/Aug/2008: If nothing else is wanted then we're done.
	//if (!fWantBlurVec && !fWantExtraInfo)
	if(true)
	{
		goto CleanUp;
	}

	if (-1L == (oldSpot = IffTell(_fileData))
		|| -1 == IffSeek(_fileData, 0, SEEK_END)
		|| -1L == (tellTemp = IffTell(_fileData))
		|| -1 == IffSeek(_fileData, oldSpot, SEEK_SET)) {
		if (pInstance->iff_error == IFF_NO_ERROR) {
			pInstance->iff_error = IFF_READ_FAILS;
		}
		goto CleanUp;
	}

	fileLength = tellTemp;

	while (1) {

		if (!iff_begin_read_chunk(pInstance, &chunkInfo)) {
			if (pInstance->iff_error == IFF_NO_ERROR) {
				pInstance->iff_error = IFF_READ_FAILS;
			}
			goto CleanUp;
		}

		if (chunkInfo.tag == IFF_TAG_CLPZ) {

			if (!iff_get_float(pInstance, &image->znear)
				|| !iff_get_float(pInstance, &image->zfar)
				|| !iff_end_read_chunk(pInstance)) {
				if (pInstance->iff_error == IFF_NO_ERROR) {
					pInstance->iff_error = IFF_READ_FAILS;
				}
				goto CleanUp;
			}

#ifdef __IFF_DEBUG_
			printf("Got clipping info: %f %f\n", image->znear, image->zfar);
#endif

		} /* END CLPZ chunk */

		else if (chunkInfo.tag == IFF_TAG_ESXY) {
			if (!iff_get_float(pInstance, &image->zesx)
				|| !iff_get_float(pInstance, &image->zesy)
				|| !iff_end_read_chunk(pInstance)) {
				if (pInstance->iff_error == IFF_NO_ERROR) {
					pInstance->iff_error = IFF_READ_FAILS;
				}
				goto CleanUp;
			}

#ifdef __IFF_DEBUG_
			printf("Got esxy info: %f %f\n", image->zesx, image->zesy);
#endif
		} /* END ESXY chunk */

		else if (chunkInfo.tag == IFF_TAG_FOR4) {

			// Leo 27/Aug/2008: Only dig into the IFF_TAG_BLUR if the caller wants the blur vec.
			if (chunkInfo.chunkType == IFF_TAG_BLUR && fWantBlurVec) {

				// -- FIXME: GET THE BLUR INFO HERE
				if (image->blurvec) {
					free(image->blurvec);
				}

				while (1) {

					if (!iff_begin_read_chunk(pInstance, &chunkInfo)) {
						if (pInstance->iff_error == IFF_NO_ERROR) {
							pInstance->iff_error = IFF_READ_FAILS;
						}
						goto CleanUp;
					}

					if (chunkInfo.tag == IFF_TAG_BLRT) {

						// read in values, uncompressed and in a linear sort
						// of manner, uh huh...
						/*
						printf( "%d\n", iff_get_long( pInstance ) );
						printf( "%d\n", iff_get_long( pInstance ) );
						printf( "%d\n", iff_get_long( pInstance ) );
						printf( "%d\n", iff_get_long( pInstance ) );
						*/

						if (!iff_get_long(pInstance, 0) // -- Don't know what these are
							|| !iff_get_long(pInstance, 0)
							|| !iff_get_long(pInstance, 0)
							|| !iff_get_long(pInstance, 0)) {
							if (pInstance->iff_error == IFF_NO_ERROR) {
								pInstance->iff_error = IFF_READ_FAILS;
							}
							goto CleanUp;
						}

						image->blurvec = (Float32 *)malloc(npixels * 2 * sizeof(Float32));

						if (!(image->blurvec)) {
							if (pInstance->iff_error == IFF_NO_ERROR) {
								pInstance->iff_error = IFF_READ_FAILS;
							}
							goto CleanUp;
						}

						for (i = 0; i < npixels; i++) {
							if (!iff_get_float(pInstance, &image->blurvec[2 * i])
								|| !iff_get_float(pInstance, &image->blurvec[2 * i + 1])) {
								if (pInstance->iff_error == IFF_NO_ERROR) {
									pInstance->iff_error = IFF_READ_FAILS;
								}
								goto CleanUp;
							}
						}

						if (!iff_end_read_chunk(pInstance)) {
							if (pInstance->iff_error == IFF_NO_ERROR) {
								pInstance->iff_error = IFF_READ_FAILS;
							}
							goto CleanUp;
						}

						break;
					}

					else {
#ifdef __IFF_DEBUG_
						printf("Skipping Chunk in search of BLRT: %c%c%c%c\n",
							(((&pInstance->chunkStack[pInstance->chunkDepth].tag)[0] >> 24) & 0xFF),
							(((&pInstance->chunkStack[pInstance->chunkDepth].tag)[0] >> 16) & 0xFF),
							(((&pInstance->chunkStack[pInstance->chunkDepth].tag)[0] >> 8) & 0xFF),
							(((&pInstance->chunkStack[pInstance->chunkDepth].tag)[0] >> 0) & 0xFF));
#endif

						if (!iff_end_read_chunk(pInstance)) {
							if (pInstance->iff_error == IFF_NO_ERROR) {
								pInstance->iff_error = IFF_READ_FAILS;
							}
							goto CleanUp;
						}
					}
				}
#ifdef __IFF_DEBUG_
				printf("Found FOR4 BLUR\n");
#endif
			}

			if (!iff_end_read_chunk(pInstance)) {
				if (pInstance->iff_error == IFF_NO_ERROR) {
					pInstance->iff_error = IFF_READ_FAILS;
				}
				goto CleanUp;
			}
		}

		else {

#ifdef __IFF_DEBUG_
			// Skipping unused data IN THE BEGINNING OF THE FILE
			printf("Skipping Chunk in search of CLPZ: %c%c%c%c\n",
				(((&pInstance->chunkStack[pInstance->chunkDepth].tag)[0] >> 24) & 0xFF),
				(((&pInstance->chunkStack[pInstance->chunkDepth].tag)[0] >> 16) & 0xFF),
				(((&pInstance->chunkStack[pInstance->chunkDepth].tag)[0] >> 8) & 0xFF),
				(((&pInstance->chunkStack[pInstance->chunkDepth].tag)[0] >> 0) & 0xFF));
#endif

			if (!iff_end_read_chunk(pInstance)) {
				if (pInstance->iff_error == IFF_NO_ERROR) {
					pInstance->iff_error = IFF_READ_FAILS;
				}
				goto CleanUp;
			}
		}

		tellTemp = IffTell(_fileData);

		if (-1L == tellTemp) {
			if (pInstance->iff_error == IFF_NO_ERROR) {
				pInstance->iff_error = IFF_READ_FAILS;
			}
			goto CleanUp;
		}

		if ((width*height + tellTemp) > fileLength) {

#ifdef __IFF_DEBUG_
			printf("End of parsable data, time to quit\n");
#endif

			break;
		}

	}

CleanUp:
	//auto val2 = mat.at<cv::uint16_t>(159, 255);
	//printf("val2 :%d", val2);

	//mat.convertTo(mat2, CV_8U, 0.00390625);
	//mat.convertTo(mat2, CV_8U);
	float zmax = -FLT_MAX;
	float zmin = FLT_MAX;
	Float32* matPtr = (Float32*)mat.data;
	uint8_t* mat2Ptr = (uint8_t*)mat2.data;

	cv::Vec3b* colorMatPtr = (cv::Vec3b*)colorMat.data;

	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			auto val = image->zbuffer[i * width + j];
			if (val != 0) {
				if (val > zmax) {
					zmax = val;
				}
				if (val < zmin) {
					zmin = val;
				}
			}
		}
	}

	//image->zbuffer[base + i * width + j] = *((Float32 *)&value);
	//auto znear = image->znear;
	//auto zfar = image->zfar;
	for (int i = 0; i < height; ++i) {
		for (int j = 0; j < width; ++j) {
			auto val = image->zbuffer[i * width + j];
			//uint16_t uval16 = ((val - zmin) / (zmax - zmin)) * -65535;
			if (val >= zmin && val <= zmax) {
				Float32 fval32 = (1.0 / -val);
				uchar uval = 0;
				if (zmax > zmin)
					uval = ((val - zmin) / (zmax - zmin)) * -255 - 1;
				else
					uval = 255;
				mat2Ptr[i * width + j] = uval;
				matPtr[i * width + j] = fval32;
			}

			auto colorVal = &(image->rgba[(i*width + j)*depth]);
			int bit = 1;
			colorMatPtr[i*width + j][2] = *colorVal;
			colorMatPtr[i*width + j][1] = *(colorVal + 1 * bit);
			colorMatPtr[i*width + j][0] = *(colorVal + 2 * bit);
			colorMatPtr[i*width + j][3] = *(colorVal + 3 * bit);
		}
	}
	//if (value != 0)
	//{
	//	int xx = x1 + j;
	//	int yy = y1 + i;
	//	float fval = -(*((Float32 *)&value));
	//	//uint16_t val = (fval * -10000);
	//	//uint16_t val = ((fval - znear) / (zfar - znear)) * 65536;
	//	//uchar uval = ((fval - znear) / (zfar - znear)) * 256;
	//	mat2Ptr[base + i * width + j] = val;
	//	//
	//	//mat.at<cv::uint16_t>(yy, xx) = val;
	//	matPtr[base + i * width + j] = uval;

	//}
	//( -1.0f * (*( (Float32 *)&value ) ));

	//} else {
	//image->zbuffer[base + i*width+j] = 0.0;
	//}

	std::string path = fileName.substr(0, fileName.rfind('\\') + 1);
	std::string pngName = path + std::string("depth_")+fileName.substr(path.size(), fileName.size() - 4 - path.size()) + ".png";



	std::string pngName2 = path + std::string("depth_8bit_")+fileName.substr(path.size(), fileName.size() - 4 - path.size()) + ".png";
	std::string pngName3 = path + std::string("color_")+fileName.substr(path.size(), fileName.size() - 4 - path.size()) + ".png";
	cv::flip(mat, mat, 0);
	cv::imwrite(pngName, mat);
	cv::flip(mat2, mat2, 0);
	cv::imwrite(pngName2, mat2);
	cv::flip(colorMat, colorMat, 0);
	cv::imwrite(pngName3, colorMat);
	free(tileData);

	IffClose(_fileData);

	//if (pOutIffError) {
	//	*pOutIffError = pInstance->iff_error;
	//}

	if (pInstance->iff_error != IFF_NO_ERROR) {
		iff_free(image);
		image = 0;
	}
	else {
		iff_free(image);
		image = 0;
	}

	free(pInstance);
	//delete filenamechar;

	//return(image);
	return 0;
}



void iff_free(iff_image * image)
{
	if (image) {

		if (image->rgba) {
			free(image->rgba);
		}
		if (image->zbuffer) {
			free(image->zbuffer);
		}
		if (image->blurvec) {
			free(image->blurvec);
		}

		free(image);

		image = 0;
	}
}
