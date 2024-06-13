#include"renew_iff.h"
// -- Function prototypes, local to this file.


/*
* -- Basic input functions.
*
*/

int iff_get_short(iff_instance *pInstance, uInt16 *pOut)
{
	uByte buf[2];
	size_t result = 0;

	result = IffRead(pInstance->pData, buf, 2, 1);

	if (result != 1) {
		if (pInstance->iff_error == IFF_NO_ERROR) {
			pInstance->iff_error = IFF_READ_FAILS;
		}
		return(false);
	}

	if (pOut) {
		*pOut = ((buf[0] << 8) + (buf[1]));
	}
	return(true);
}

int iff_get_long(iff_instance *pInstance, uInt32 *pOut)
{
	uByte buffer[4];

	size_t result = IffRead(pInstance->pData, buffer, 4, 1);
	if (result != 1) {
		if (pInstance->iff_error == IFF_NO_ERROR) {
			pInstance->iff_error = IFF_READ_FAILS;
		}
		return(false);
	}

	if (pOut) {
		*pOut = (buffer[0] << 24) + (buffer[1] << 16)
			+ (buffer[2] << 8) + (buffer[3] << 0);
	}
	return(true);
}

int iff_get_float(iff_instance *pInstance, Float32 *pOut)
{
	uByte buffer[4];
	uInt32 value;

	size_t result = IffRead(pInstance->pData, buffer, 4, 1);
	if (result != 1) {
		if (pInstance->iff_error == IFF_NO_ERROR) {
			pInstance->iff_error = IFF_READ_FAILS;
		}
		return(false);
	}

	if (pOut) {
		value = (buffer[3] << 24) + (buffer[2] << 16)
			+ (buffer[1] << 8) + (buffer[0] << 0);

		*pOut = *((Float32 *)&value);
	}
	return(true);
}


/*
* IFF Chunking Routines.
*
*/

int iff_begin_read_chunk(iff_instance *pInstance, iff_chunk *pOutChunk)
{
	long tempTell;

	pInstance->chunkDepth++;
	if ((pInstance->chunkDepth >= CHUNK_STACK_SIZE) || (pInstance->chunkDepth < 0)) {
		if (pInstance->iff_error == IFF_NO_ERROR) {
			pInstance->iff_error = IFF_BAD_STACK;
		}
		memset(pOutChunk, 0, sizeof(iff_chunk));
		return(false);
	}

	tempTell = IffTell(pInstance->pData);

	if (tempTell == -1L) {
		if (pInstance->iff_error == IFF_NO_ERROR) {
			pInstance->iff_error = IFF_BAD_STACK;
		}
		memset(pOutChunk, 0, sizeof(iff_chunk));
		return(false);
	}

	pInstance->chunkStack[pInstance->chunkDepth].start = tempTell;

	if (!iff_get_long(pInstance, &pInstance->chunkStack[pInstance->chunkDepth].tag)
		|| !iff_get_long(pInstance, &pInstance->chunkStack[pInstance->chunkDepth].size)) {
		if (pInstance->iff_error == IFF_NO_ERROR) {
			pInstance->iff_error = IFF_READ_FAILS;
		}
		memset(pOutChunk, 0, sizeof(iff_chunk));
		return(false);
	}

	if (pInstance->chunkStack[pInstance->chunkDepth].tag == IFF_TAG_FOR4) {
		// -- We have a form, so read the form type tag as well. 
		if (!iff_get_long(pInstance, &pInstance->chunkStack[pInstance->chunkDepth].chunkType)) {
			if (pInstance->iff_error == IFF_NO_ERROR) {
				pInstance->iff_error = IFF_READ_FAILS;
			}
			memset(pOutChunk, 0, sizeof(iff_chunk));
			return(false);
		}
	}
	else {
		pInstance->chunkStack[pInstance->chunkDepth].chunkType = 0;
	}
	*pOutChunk = pInstance->chunkStack[pInstance->chunkDepth];
	return(true);
}

int iff_end_read_chunk(iff_instance *pInstance)
{
	uInt32 end;
	int part;

	if ((pInstance->chunkDepth >= CHUNK_STACK_SIZE) || (pInstance->chunkDepth < 0)) {
		if (pInstance->iff_error == IFF_NO_ERROR) {
			pInstance->iff_error = IFF_BAD_STACK;
		}
		return(false);
	}

	end = pInstance->chunkStack[pInstance->chunkDepth].start + pInstance->chunkStack[pInstance->chunkDepth].size + 8;

	if (pInstance->chunkStack[pInstance->chunkDepth].chunkType != 0) {
		end += 4;
	}
	// Add padding 
	part = end % 4;
	if (part != 0) {
		end += 4 - part;
	}

	if (-1 == IffSeek(pInstance->pData, end, SEEK_SET)) {
		if (pInstance->iff_error == IFF_NO_ERROR) {
			pInstance->iff_error = IFF_BAD_STACK;
		}
		return(false);
	}

#ifdef __IFF_DEBUG_
	printf("Closing Chunk: %c%c%c%c\n\n",
		(((&pInstance->chunkStack[pInstance->chunkDepth].tag)[0] >> 24) & 0xFF),
		(((&pInstance->chunkStack[pInstance->chunkDepth].tag)[0] >> 16) & 0xFF),
		(((&pInstance->chunkStack[pInstance->chunkDepth].tag)[0] >> 8) & 0xFF),
		(((&pInstance->chunkStack[pInstance->chunkDepth].tag)[0] >> 0) & 0xFF));
#endif

	pInstance->chunkDepth--;

	return(true);
}

std::vector<uByte> iff_read_data(iff_instance *pInstance, int size)
{
	//std::shared_ptr<uByte> buffer = std::make_shared<std::array<uByte,size> >();//(uByte *)malloc( size * sizeof( uByte ) );
	std::vector<uByte> buffer(size);
	size_t result = 0;

	if (buffer.size() < size) {
		if (pInstance->iff_error == IFF_NO_ERROR) {
			pInstance->iff_error = IFF_READ_FAILS;
		}
		return buffer;;
	}

#ifdef __IFF_DEBUG_
	printf("read_data  size: %d\n", size);
#endif

	result = IffRead(pInstance->pData, &buffer[0], size, 1);

	if (result != 1) {
		if (pInstance->iff_error == IFF_NO_ERROR) {
			pInstance->iff_error = IFF_READ_FAILS;
		}
		//free( buffer );
		return buffer;
	}

	return(buffer);
}



/*
* Compression Routines
*
*/

uByte * iff_decompress_rle(iff_instance *pInstance, uInt32 numBytes,
	uByte * compressedData,
	uInt32 compressedDataSize,
	uInt32 * compressedIndex)
{
	uByte * data = (uByte *)malloc(numBytes * sizeof(uByte));
	uByte nextChar, count;
	int i;
	uInt32 byteCount = 0;

	if (!data) {
		if (pInstance->iff_error == IFF_NO_ERROR) {
			pInstance->iff_error = IFF_READ_FAILS;
		}
		return(0);
	}

	memset(data, 0, numBytes * sizeof(uByte));

#ifdef __IFF_DEBUG_
	printf("Decompressing data %d\n", numBytes);
#endif

	while (byteCount < numBytes) {

		if (*compressedIndex >= compressedDataSize) {
			break;
		}

		nextChar = compressedData[*compressedIndex];
		(*compressedIndex)++;

		count = (nextChar & 0x7f) + 1;
		if ((byteCount + count) > numBytes) break;

		if (nextChar & 0x80) {

			// We have a duplication run

			nextChar = compressedData[*compressedIndex];
			(*compressedIndex)++;

			//assert( ( byteCount + count ) <= numBytes ); 
			for (i = 0; i < count; ++i) {
				data[byteCount] = nextChar;
				byteCount++;
			}
		}
		else {
			// We have a verbatim run 
			for (i = 0; i < count; ++i) {

				data[byteCount] = compressedData[*compressedIndex];
				(*compressedIndex)++;
				byteCount++;
			}
		}
		assert(byteCount <= numBytes);
	}

	return(data);
}

std::vector<uByte> iff_decompress_tile_rle(iff_instance *pInstance,
	uInt16 width, uInt16 height, uInt16 depth,
	uByte * compressedData,
	uInt32 compressedDataSize)
{
	uByte * channels[4];
	std::vector<uByte> data;
	int i, k, row, column;
	uInt32 compressedIndex = 0;

#ifdef __IFF_DEBUG_    
	printf("Decompressing tile [ %hd, ", width);
	printf("%hd, ", height);
	printf("%hd ]\n", depth);
#endif

	// -- Decompress the different channels (RGBA)
	// -- ERROR CHECK, MUST OPERATE ON RGBA !!!
	// Leo Davidson 28/Aug/2008: Allow depths of 3 and 1, and check compressedData too. Not sure why the original code required RGBA as it works fine with RGB and looks like it'd work with grey.
	if (!compressedData || (depth != 4 && depth != 3 && depth != 1)) {
		if (pInstance->iff_error == IFF_NO_ERROR) {
			pInstance->iff_error = IFF_BAD_COMPRESS;
		}
		return(data);
	}

	for (i = (depth - 1); i >= 0; --i) {

		channels[i] = iff_decompress_rle(pInstance, width * height, compressedData, compressedDataSize, &compressedIndex);

		if (!channels[i]) {
			while (--i >= 0) {
				free(channels[i]);
			}
			if (pInstance->iff_error == IFF_NO_ERROR) {
				pInstance->iff_error = IFF_READ_FAILS;
			}
			return(data);
		}
	}

	// -- Pack all of the channels from the decompression into an RGBA array.
	data.resize(width * height * depth); //= //(uByte *)malloc( width * height * depth * sizeof( uByte ) );

	// Leo Davidson 28/Aug/2008: Check malloc result.
	if (data.size() == 0) {
		if (pInstance->iff_error == IFF_NO_ERROR) {
			pInstance->iff_error = IFF_READ_FAILS;
		}
		for (i = 0; i < depth; i++) {
			free(channels[i]);
		}
		return(data);
	}

	for (row = 0; row < height; row++) {
		for (column = 0; column < width; column++) {
			for (k = 0; k < depth; k++) {
				data[depth*(row*width + column) + k] = channels[k][row*width + column];
			}
		}
	}


	for (i = 0; i < depth; i++) {
		free(channels[i]);
	}

	return(data);
}

uByte * iff_read_uncompressed_tile(iff_instance *pInstance,
	uInt16 width, uInt16 height, uInt16 depth)
{
	uByte * data = 0;
	uByte pixel[4];
	int i, j, d, index;
	size_t result;
	data = (uByte *)malloc(width * height * depth * sizeof(uByte));

	// Leo Davidson 28/Aug/2008: Check malloc result.
	if (!data) {
		if (pInstance->iff_error == IFF_NO_ERROR) {
			pInstance->iff_error = IFF_READ_FAILS;
		}
		return(0);
	}


#ifdef __IFF_DEBUG_
	printf("Begin reading uncompressed tile\n", NULL);
#endif

	for (i = 0; i < height; i++)
	{
		index = i * width * depth;
		for (j = 0; j < width; ++j)
		{
			result = IffRead(pInstance->pData, pixel, depth, 1);
			if (result != 1) {
				if (pInstance->iff_error == IFF_NO_ERROR) {
					pInstance->iff_error = IFF_READ_FAILS;
				}
				free(data);
				return(0);
			}
			for (d = (depth - 1); d >= 0; --d)
			{
				data[index] = pixel[d];
				++index;
			}
		}
	}

#ifdef __IFF_DEBUG_
	printf("End reading uncompressed tile\n", NULL);
#endif

	return(data);
}