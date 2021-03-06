/***************************************************************************
 *
 * Authors:    David Strelak (davidstrelak@gmail.com)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#ifndef CUDA_GPU_MOVIE_CORRELATION_KERNELS
#define CUDA_GPU_MOVIE_CORRELATION_KERNELS

#include "reconstruction/movie_alignment_gpu_defines.h"

/**
 * Kernel performing scaling of the 2D FFT images, with possible normalization,
 * filtering and centering
 * @param in input data
 * @param out output data
 * @param noOfImages images to process
 * @param inX input X dim
 * @param inY input Y dim
 * @param outX output X dim (must be less or equal to input)
 * @param outY output Y dim (must be less or equal to input)
 * @param filter to apply. Only if 'applyFilter' is true
 * @param normFactor normalization factor to use (x *= normFactor).
 * Only if 'normalize' is true
 */
template<typename T, typename U, bool applyFilter, bool normalize, bool center>
__global__
void scaleFFT2DKernel(const T* __restrict__ in, T* __restrict__ out,
        int noOfImages, size_t inX, size_t inY, size_t outX, size_t outY,
    const U* __restrict__ filter, U normFactor) {
    // assign pixel to thread
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int idy = blockIdx.y*blockDim.y + threadIdx.y;

    if (idx >= outX || idy >= outY ) return;
    size_t fIndex = idy*outX + idx; // index within single image
    U filterCoef = filter[fIndex];
    U centerCoef = 1-2*((idx+idy)&1); // center FFT, input must have be even
    int yhalf = outY/2;

    size_t origY = (idy <= yhalf) ? idy : (inY - (outY-idy)); // take top N/2+1 and bottom N/2 lines
    for (int n = 0; n < noOfImages; n++) {
        size_t iIndex = n*inX*inY + origY*inX + idx; // index within consecutive images
        size_t oIndex = n*outX*outY + fIndex; // index within consecutive images
        out[oIndex] = in[iIndex];
        if (applyFilter) {
            out[oIndex] *= filterCoef;
        }
        if (normalize) {
            out[oIndex] *= normFactor;
        }
        if (center) {
            out[oIndex] *= centerCoef;
        }
    }
}

/**
 * Function computes correlation between two batches of images
 * @param in1 first batch of images
 * @param in2 second batch of images
 * @param correlations resulting correlations
 * @param xDim of the images
 * @param yDim of the images
 * @param isWithin if true, in1 == in2
 * @param iStart first index in the first batch to start with
 * @param iStop last index in the first batch to process (included)
 * @param jStart first index in the second batch to process
 * @param jStop last index in the second batch to process (included)
 * @param jSize size of the second batch
 */
template<typename T>
__global__
void correlate(const T* __restrict__ in1, const T* __restrict__ in2,
        T* correlations, int xDim, int yDim,
        bool isWithin, int iStart, int iStop, int jStart, int jStop, size_t jSize) {
    // assign pixel to thread
#if TILE > 1
    int id = threadIdx.y * blockDim.x + threadIdx.x;
    int tidX = threadIdx.x % TILE + (id / (blockDim.y * TILE)) * TILE;
    int tidY = (id / TILE) % blockDim.y;
    int idx = blockIdx.x*blockDim.x + tidX;
    int idy = blockIdx.y*blockDim.y + tidY;
#else
    volatile int idx = blockIdx.x*blockDim.x + threadIdx.x;
    volatile int idy = blockIdx.y*blockDim.y + threadIdx.y;
#endif
    int a = 1-2*((idx+idy)&1); // center FFT, input must have be even

    if (idx >= xDim || idy >= yDim ) return;
    size_t pixelIndex = idy*xDim + idx; // index within single image

    bool compute = false;
    int counter = 0;
    for (int i = iStart; i <= iStop; i++) {
        int tmpOffset = i * xDim * yDim;
        T tmp = in1[tmpOffset + pixelIndex];
        for (int j = isWithin ? i + 1 : 0; j < jSize; j++) {
            if (!compute) {
                compute = true;
                j = jStart;
                continue; // skip first iteration
            }
            if (compute) {
                int tmp2Offset = j * xDim * yDim;
                T tmp2 = in2[tmp2Offset + pixelIndex];
                T res;
                res.x = (tmp.x*tmp2.x) + (tmp.y*tmp2.y);
                res.y = (tmp.y*tmp2.x) - (tmp.x*tmp2.y);
                correlations[counter*xDim*yDim + pixelIndex] = res * a;
                counter++;
            }
            if ((iStop == i) && (jStop == j)) {
                return;
            }
        }
    }
}

/**
 * Function to crop (square) center part of an image
 * @param in input images
 * @param out output images
 * @param xDim input X dim
 * @param yDim output Y dim
 * @param noOfImgs no if images to process
 * @param outDim output dimension
 */
template<typename T>
__global__
void cropSquareInCenter(const T* __restrict__ in, T* out, int xDim, int yDim, int noOfImgs,
        int outDim) {
    // assign pixel to thread
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    int idy = blockIdx.y*blockDim.y + threadIdx.y;

    if (idx >= outDim || idy >= outDim ) return;

    int inputImgSize = xDim * yDim;
    int outputImgSize = outDim * outDim;

    int inCenterX = (int)((T) (xDim) / (T)2);
    int inCenterY = (int)((T) (yDim) / (T)2);

    int outCenter = (int)((T) (outDim) / (T)2);

    for (int n = 0; n < noOfImgs; n++) {
        int iX = idx - outCenter + inCenterX;
        int iY = idy - outCenter + inCenterY;
        int inputPixelIdx = (n * inputImgSize) + (iY * xDim) + iX;
        int outputPixelIdx = (n * outputImgSize) + (idy * outDim) + idx;
        out[outputPixelIdx] = in[inputPixelIdx];
    }
}

#endif //* CUDA_GPU_MOVIE_CORRELATION_KERNELS */
