#pragma once
#include "IOHandler.h"
#include "Cluster.h"
#include "Point.h"


void calculateClustersDiameterOmp(Point* points, int numOfPoints, Cluster* clusters, int numOfClusters);

// Increases each point in the given collection by time with the given velocities
// (for each point x = x + (dt * moment) * vxi ,
// y = y + (dt * moment) * vyi)
void increaseTimeOmp(Point* points, int numOfPoints, double dt, int moment);

void calculateClustersCentersOmp(Point* points, int numOfPoints, Cluster* clusters, int numOfClusters);

// Assigning all points to clusters
// Returns: TRUE if the point's cluster has changed, FALSE if not
Boolean assignClustersToPointsOmp(Point* points, int numOfPoints, Cluster* clusters, int numOfClusters);