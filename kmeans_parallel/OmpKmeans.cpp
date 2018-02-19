#include "OmpKmeans.h"
#include <stdlib.h>
#include <omp.h>
#include "Point.h"

void calculateClustersDiameterOmp(Point* points, int numOfPoints, Cluster* clusters, int numOfClusters) {
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

void increaseTimeOmp(Point* points, int numOfPoints, double dt, int moment) {
	double timeInterval = dt * moment;

#pragma omp parallel for
	for (int i = 0; i < numOfPoints; i++) {
		points[i].position.x += timeInterval * points[i].velocity.vx;
		points[i].position.y += timeInterval * points[i].velocity.vy;
	}
}

void calculateClustersCentersOmp(Point* points, int numOfPoints, Cluster* clusters, int numOfClusters) {
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
Boolean assignClustersToPointsOmp(Point* points, int numOfPoints, Cluster* clusters, int numOfClusters) {
	// Indicating if at least one point has changed cluster
	Boolean pointChanged = FALSE;

	// A counter array for points each thread assigns
	int* clustersPointsCounter = (int*)calloc(omp_get_max_threads() * numOfClusters, sizeof(int));
	
	if (clustersPointsCounter == NULL) {
		printf("\nAllocation error assignClustersToPointsOmp\n");
		exit(0);
	}

#pragma omp parallel for
	for (int pointIndex = 0; pointIndex < numOfPoints; pointIndex++) {
		// Assigning the curr point and checking if changed cluster
		if (assignClusterToPoint(&(points[pointIndex]), clusters, numOfClusters))
			pointChanged = TRUE;

		int pointClusterId = points[pointIndex].cluster->id;
		int threadId = omp_get_thread_num();

		// Increasing the counter for the chosen cluster of the current thread 
		clustersPointsCounter[threadId * numOfClusters + pointClusterId] += 1;
	}

	// Sum of point for each cluster
#pragma omp parallel for
	for (int clusterId = 0; clusterId < numOfClusters; clusterId++) {
		clusters[clusterId].numOfPoints = 0;
		for (int threadId = 0; threadId < omp_get_num_threads(); threadId++) {
			clusters[clusterId].numOfPoints += clustersPointsCounter[threadId * numOfClusters + clusterId];
		}
	}

	return pointChanged;
}