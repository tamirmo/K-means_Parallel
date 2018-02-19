#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "CudaKmeans.h"
#include <time.h>

void calculateClustersDiameterCuda(Point* points, int numOfPoints, Cluster* clusters, int numOfClusters) {
#pragma omp parallel for
	for (int clusterIndex = 0; clusterIndex < numOfClusters; clusterIndex++) {
		double maxDistance = 0, currDistance;

		// Going over all points and calculating distance with each other point in the cluster
		// to get the maximum distance
		for (int i = 0; i < numOfPoints - 1; i++)
			if (points[i].cluster->id == clusterIndex)
				for (int j = i + 1; j < numOfPoints; j++)
					// Calculating distance for points in the same cluster
					if (points[j].cluster->id == clusterIndex) {
						currDistance = getPointsDistance(&(points[i].position), &(points[j].position));
						if (currDistance > maxDistance)
							maxDistance = currDistance;
					}

		clusters[clusterIndex].diameter = maxDistance;
	}
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
		//cudaFree(&dev_data, &dev_threadsCounterArray, &dev_histogram);
		return "cudaDeviceSynchronize returned error code after launching countKernel!";
	}

	// Copy histogram result vector from GPU buffer to host memory
	cudaStatus = cudaMemcpy(pointsArr, dev_points, numOfPoints * sizeof(Point), cudaMemcpyDeviceToHost);
	if (cudaStatus != cudaSuccess) {
		//cudaFree(&dev_data, &dev_threadsCounterArray, &dev_histogram);
		return "cudaMemcpy failed!";
	}

	// Free all GPU memory
	//cudaFree(&dev_data, &dev_threadsCounterArray, &dev_histogram);
	return NULL;
}

const char* stopCuda() {
	cudaError_t cudaStatus;

	// was at main at first, need to be checked 
	// cudaDeviceReset must be called before exiting
	cudaStatus = cudaDeviceReset();
	if (cudaStatus != cudaSuccess) {
		//cudaFree(&dev_data, &dev_threadsCounterArray, &dev_histogram);
		return "cudaDeviceReset failed";
	}

	return NULL;
}

const char* initCuda() {
	cudaError_t cudaStatus;

	// Choose which GPU to run on (our system has only one)
	cudaStatus = cudaSetDevice(0);
	if (cudaStatus != cudaSuccess) {
		//cudaFree(&dev_data, &dev_threadsCounterArray, &dev_histogram);
		return "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?";
	}

	return NULL;
}

const char* freeCuda(Point* gpu_points) {
	if (cudaFree(gpu_points) != cudaSuccess)
		return "cudaFree error";
	return NULL;
}

void calculateClustersCentersCuda(Point* points, int numOfPoints, Cluster* clusters, int numOfClusters) {
#pragma omp parallel for
	for (int clusterIndex = 0; clusterIndex < numOfClusters; clusterIndex++) {
		double sumY = 0, sumX = 0;

		// When a cluster has no points, we keep it's center intact
		if (clusters[clusterIndex].numOfPoints != 0) {

			for (int pointIndex = 0; pointIndex < numOfPoints; pointIndex++) {
				// Calculating sum of all point in the cluster
				if (points[pointIndex].cluster->id == clusterIndex) {
					sumX += points[pointIndex].position.x;
					sumY += points[pointIndex].position.y;
				}
			}

			// Each cluster's center is the average of points positions
			clusters[clusterIndex].center.x = sumX / clusters[clusterIndex].numOfPoints;
			clusters[clusterIndex].center.y = sumY / clusters[clusterIndex].numOfPoints;
		}
	}
}

// Assigning all points to clusters
// Returns: TRUE if the point's cluster has changed, FALSE if not
Boolean assignClustersToPointsCuda(Point* points, int numOfPoints, Cluster* clusters, int numOfClusters) {
	int pointIndex;
	// Indicating if at least one point has changed cluster
	Boolean pointChanged = FALSE;

#pragma omp parallel for
	for (pointIndex = 0; pointIndex < numOfPoints; pointIndex++)
		// Assigning the curr point and checking if changed cluster
		if (assignClusterToPoint(&(points[pointIndex]), clusters, numOfClusters))
			pointChanged = TRUE;

	return pointChanged;
}