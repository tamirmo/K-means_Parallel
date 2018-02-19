#pragma once
#include "IOHandler.h"
#include "Cluster.h"
#include "Point.h"
#include <math.h>

#define BLOCK_SIZE 1024
// This number is the maximum block size for cuda as mentioned:
// http://docs.nvidia.com/cuda/cuda-c-programming-guide/index.html
#define MAX_BLOCKS pow(2, 31) - 1
// The portion of points increing via cuda (other in OMP) 
#define CUDA_INC_TIME_PORTION 0.9

// Setting the cuda device (0)
const char* initCuda();

// Calling device reset
const char* stopCuda();

const char*  freeCuda(Point* gpu_points);

// Increases each point in the given collection by time with the given velocities
// (for each point x = x + (dt * moment) * vxi ,
// y = y + (dt * moment) * vyi
// This function starts the kernel without calling sync
const char* increaseTimeCudaStart(Point* points, int numOfPoints, double dt, int moment, Point** gpu_points);

//  This function calls sync and returns the increased points to pointsArr
const char* increaseTimeCudaEnd(Point* dev_points, Point* pointsArr, int numOfPoints);