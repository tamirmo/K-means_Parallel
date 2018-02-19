#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "CudaKmeans.h"
#include <time.h>

// Calling device reset
const char* stopCuda() {
	cudaError_t cudaStatus;

	// was at main at first, need to be checked 
	// cudaDeviceReset must be called before exiting
	cudaStatus = cudaDeviceReset();
	if (cudaStatus != cudaSuccess)
		return "cudaDeviceReset failed";

	return NULL;
}

// Setting the cuda device (0)
const char* initCuda() {
	cudaError_t cudaStatus;

	// Choose which GPU to run on (our system has only one)
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) 
		return "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?";

	return NULL;
}

__global__ void increaseTimeKernel(Point *dev_pointArr, double timeInterval, int numOfPoints) {
	// Starting from the thread's id
	// increasing the point index by grid size
	for (int i = blockIdx.x * blockDim.x + threadIdx.x;
		i < numOfPoints; i += blockDim.x * gridDim.x)
	{
		dev_pointArr[i].position.x += timeInterval * dev_pointArr[i].velocity.vx;
		dev_pointArr[i].position.y += timeInterval * dev_pointArr[i].velocity.vy;
	}
}

// Increases each point in the given collection by time with the given velocities
// (for each point x = x + (dt * moment) * vxi ,
// y = y + (dt * moment) * vyi)
const char* increaseTimeCudaStart(Point* points, int numOfPoints, double dt, int moment, Point** gpu_points) {
	double timeInterval = dt * moment;
	Point *dev_points = 0;
	int numOfBlocks;
	cudaError_t cudaStatus;

	numOfBlocks = (numOfPoints + BLOCK_SIZE - 1) / BLOCK_SIZE;
	
	if (numOfBlocks > MAX_BLOCKS)
		numOfBlocks = MAX_BLOCKS;

	if (*gpu_points == NULL) {
		// Allocate GPU buffers for the points array
		cudaStatus = cudaMalloc((void**)&dev_points, numOfPoints * sizeof(Point));
		if (cudaStatus != cudaSuccess)
			return "cudaMalloc failed!";
		
		*gpu_points = dev_points;
	}
	else
		dev_points = *gpu_points;

	// Copy data array from host memory to GPU buffers
	cudaStatus = cudaMemcpy(dev_points, points, numOfPoints * sizeof(Point), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		cudaFree(dev_points);
		return "cudaMemcpy failed!";
	}

	// Launch a kernel increasing time for one part of the points on the GPU with
	increaseTimeKernel << <numOfBlocks, BLOCK_SIZE >> >(dev_points, timeInterval, numOfPoints);

	// Check for any errors launching the kernel
	cudaStatus = cudaGetLastError();
	if (cudaStatus != cudaSuccess) {
		cudaFree(dev_points);
		return "countKernel launch failed";
	}

	return NULL;
}

const char* increaseTimeCudaEnd(Point* dev_points, Point* pointsArr, int numOfPoints) {
	cudaError_t cudaStatus;

	// cudaDeviceSynchronize waits for the kernel to finish, and returns
	// any errors encountered during the launch.
	cudaStatus = cudaDeviceSynchronize();
	if (cudaStatus != cudaSuccess) {
		cudaFree(dev_points);
		return "cudaDeviceSynchronize returned error code after launching countKernel!";
	}

	// Copy histogram result vector from GPU buffer to host memory
	cudaStatus = cudaMemcpy(pointsArr, dev_points, numOfPoints * sizeof(Point), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		cudaFree(dev_points);
		return "cudaMemcpy failed!";
	}

	// Free all GPU memory
	return NULL;
}

const char* freeCuda(Point* gpu_points) {
	if (cudaFree(gpu_points) != cudaSuccess)
		return "cudaFree error";
	return NULL;
}