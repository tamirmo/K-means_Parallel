#include "Cluster.h"
#include <stdlib.h>

// Calculating the distance between each cluster and the given point,
// finding the minumim and assigning the point to the suitable cluster
// Returns: TRUE if the point's cluster has changed, FALSE if not
Boolean assignClusterToPoint(Point* point, Cluster* clusters, int clustersCount) {
	int clusterIndex;
	double minDistance = INITIAL_DISTANCE, currDistance;
	Cluster* newCluster = NULL;

	for (clusterIndex = 0; clusterIndex < clustersCount; clusterIndex++) {
		currDistance = getPointsDistance(&(point->position), &(clusters[clusterIndex].center));

		if (currDistance < minDistance || minDistance == INITIAL_DISTANCE) {
			minDistance = currDistance;
			newCluster = &(clusters[clusterIndex]);
		}
	}

	// If this is the fist assignment
	if (point->cluster == NULL ||
		// Or the point has changed cluster
		point->cluster->id != newCluster->id) {

		// If this is not the first assignment
		//if (point->cluster != NULL)
			// Updating the points count in the cluster
			//point->cluster->numOfPoints--;

		//newCluster->numOfPoints++;
		point->cluster = newCluster;

		return TRUE;
	}

	return FALSE;
}

// Assigning K clusters that their centers are the first K point
void initClusters(InputParams* inputParams, Point* points, Cluster** clusters) {
	int i;

	if (*clusters == NULL)
		// Assigning an array of K clusters
		*clusters = (Cluster*)malloc(sizeof(Cluster) * inputParams->K);

	if (*clusters != NULL) {
		for (i = 0; i < inputParams->K; i++) {
			(*clusters)[i].center.x = points[i].position.x;
			(*clusters)[i].center.y = points[i].position.y;
			(*clusters)[i].id = i;
			(*clusters)[i].numOfPoints = 0;
		}

		// Clusters has changed, clearing the clusters of the points:
		for (i = 0; i < inputParams->N; i++)
			points[i].cluster = NULL;
	}
	else{
		printf("\nAllocation error initClusters\n");
		exit(1);
	}
}