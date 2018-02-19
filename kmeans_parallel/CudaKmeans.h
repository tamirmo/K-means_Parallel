#pragma once
#include "IOHandler.h"
#include "Cluster.h"
#include "Point.h"
#include <math.h>

#define BLOCK_SIZE 256
#define MAX_BLOCKS pow(2, 31) - 1

const char* initCuda();

const char* stopCuda();

const char*  freeCuda(Point* gpu_points);

void calculateClustersDiameterCuda(Point* points, int numOfPoints, Cluster* clusters, int numOfClusters);

// Increases each point in the given collection by time with the given velocities
// (for each point x = x + (dt * moment) * vxi ,
// y = y + (dt * moment) * vyi)
const char* increaseTimeCudaStart(Point* points, int numOfPoints, double dt, int moment, Point** gpu_points);

const char* increaseTimeCudaEnd(Point* dev_points, Point* pointsArr, int numOfPoints);

void calculateClustersCentersCuda(Point* points, int numOfPoints, Cluster* clusters, int numOfClusters);

// Assigning all points to clusters
// Returns: TRUE if the point's cluster has changed, FALSE if not
Boolean assignClustersToPointsCuda(Point* points, int numOfPoints, Cluster* clusters, int numOfClusters);