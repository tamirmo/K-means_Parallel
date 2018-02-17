#pragma once
#include "Point.h"
#include "IOHandler.h"

struct Cluster {
	int numOfPoints;
	Position center;
	double diameter;
	int id;
};

Boolean assignClusterToPoint(Point* point, Cluster* clusters, int clustersCount);

// Assigning K clusters that their centers are the first K point
void initClusters(InputParams* inputParams, Point* points, Cluster** clusters);