#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "IOFilesHandler.h"

#define INITIAL_DISTANCE -1
#define INVALID_TIME -1

// Increases each point in the given collection by time with the given velocities
// (for each point x = x + (dt * moment) * vxi ,
// y = y + (dt * moment) * vyi)
void increaseTime(Point* points, Velocity* velocities, int numOfPoints, double dt, int moment) {
	double timeInterval = dt * moment;
	for (int i = 0; i < numOfPoints; i++) {
		points[i].position.x += timeInterval * velocities[i].vx;
		points[i].position.y += timeInterval * velocities[i].vy;
	}
}

void freeRecourses(Input *fileInput, Cluster* clusters, int numOfClusters) {
	// TODO: Check if all alocations were freed

	free(fileInput->points);
	free(fileInput->velocities);

	free(fileInput);
	free(clusters);
}

// Assigning K clusters that their centers are the first K point
void initClusters(Input *fileInput, Cluster** clusters) {
	int i;

	if(*clusters == NULL)
		// Assigning an array of K clusters
		*clusters = (Cluster*)malloc(sizeof(Cluster) * fileInput->K);
	
	if (*clusters != NULL) {
		for (i = 0; i < fileInput->K; i++) {
			(*clusters)[i].center.x = fileInput->points[i].position.x;
			(*clusters)[i].center.y = fileInput->points[i].position.y;
			(*clusters)[i].id = i;
			(*clusters)[i].numOfPoints = 0;
		}

		// Clusters has changed, clearing the clusters of the points:
		for (i = 0; i < fileInput->N; i++)
			fileInput->points[i].cluster = NULL;
	}
	else
	{
		// TODO: Alloc failed
	}
}

// Calculates distance between two points
double getPointsDistance(Position* p1, Position* p2)
{
	double x = p2->x - p1->x;
	double y = p2->y - p1->y;
	return sqrt(x*x + y*y);
}

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
		if (point->cluster != NULL)
			// Updating the points count in the cluster
			point->cluster->numOfPoints--;
		
		newCluster->numOfPoints++;
		point->cluster = newCluster;

		return TRUE;
	}

	return FALSE;
}

void calculateClustersCenters(Input *fileInput, Cluster* clusters) {
	int clusterIndex, pointIndex;
	double sumX, sumY;

	for (clusterIndex = 0; clusterIndex < fileInput->K; clusterIndex++) {
		// When a cluster has no points, we keep it's center intact
		if (clusters[clusterIndex].numOfPoints != 0 ) {
			sumX = sumY = 0;

			for (pointIndex = 0; pointIndex < fileInput->N; pointIndex++) {
				// Calculating sum of all point in the cluster
				if (fileInput->points[pointIndex].cluster->id == clusterIndex) {
					sumX += fileInput->points[pointIndex].position.x;
					sumY += fileInput->points[pointIndex].position.y;
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
Boolean assignClustersToPoints(Input *fileInput, Cluster* clusters) {
	int pointIndex;
	// Indicating if at least one point has changed cluster
	Boolean pointChanged = FALSE;

	for (pointIndex = 0; pointIndex < fileInput->N; pointIndex++) 
		// Assigning the curr point and checking if changed cluster
		if (assignClusterToPoint(&(fileInput->points[pointIndex]), clusters, fileInput->K))
			pointChanged = TRUE;

	return pointChanged;
}

void calculateClustersDiameter(Input *fileInput, Cluster* clusters) {
	int clusterIndex, i, j;
	double maxDistance, currDistance;

	for (clusterIndex = 0; clusterIndex < fileInput->K; clusterIndex++) {
		maxDistance = 0;

		// Going over all points and calculating distance with each other point in the cluster
		// to get the maximum distance
		for (i = 0; i < fileInput->N - 1; i++)
			if (fileInput->points[i].cluster->id == clusterIndex)
				for (j = i+1; j < fileInput->N; j++)
					// Calculating distance for points in the same cluster
					if(fileInput->points[j].cluster->id == clusterIndex) {
						currDistance = getPointsDistance(&(fileInput->points[i].position), &(fileInput->points[j].position));
						if (currDistance > maxDistance)
							maxDistance = currDistance;
					}

		clusters[clusterIndex].diameter = maxDistance;
	}
}

// The clusters quality is an average of diameters of the cluster divided by distance to other clusters
double getClustersQuality(Input *fileInput, Cluster* clusters) {
	double q = 0;
	int i, j;

	calculateClustersDiameter(fileInput, clusters);

	for (i = 0; i < fileInput->K; i++) 
		for (j = 0; j < fileInput->K; j++) 
			if (j != i) 
				q += clusters[i].diameter / getPointsDistance(&(clusters[i].center), &(clusters[j].center));

	// Devide q by the number of addends
	q /= (fileInput->K-1) * (fileInput->K);

	return q;
}

// Applying kmeans algo for the given input
// Returns: The quality of the clusters calculated
double kmeans(Input *fileInput, Cluster* clusters) {
	int iteration;
	// Indicating if at least one point changed cluster
	// (initializing to true for the first round)
	Boolean pointChanged = TRUE;

	// Assigning clusters until no point has changed
	// or the limit iterations count has reached
	for (iteration = 0; iteration < fileInput->LIMIT && pointChanged; iteration++) {
		pointChanged = assignClustersToPoints(fileInput, clusters);
		
		// Updating the new centers
		calculateClustersCenters(fileInput, clusters);
	}

	return getClustersQuality(fileInput, clusters);
}

void main() {
	Input *fileInput;
	Cluster* clusters = NULL;
	Output output;
	double q;
	// TODO: MPI (pid)
	int timeInterval = 0;
	int numOfProcesses = 1;
	double n;
	double currTime = INVALID_TIME;
	Boolean isFinished = FALSE;

	// Reading from file
	fileInput = readInputFile(INPUT_FILE_NAME);

	for (n = timeInterval; n < (fileInput->T / fileInput->dT) && !isFinished; n += numOfProcesses) {
		initClusters(fileInput, &clusters);

		q = kmeans(fileInput, clusters);

		// Checking termination condition:
		if (q < fileInput->QM)
			isFinished = TRUE;
			// TODO: Send result to master process
		else
			// Have not reached the desired quality
			increaseTime(fileInput->points, fileInput->velocities, fileInput->N, fileInput->dT, numOfProcesses);
			// TODO: Send null as result for master process
		
		currTime = n * numOfProcesses * fileInput->dT;
	}

	output.q = q;
	output.clusters = clusters;
	output.K = fileInput->K;
	output.t = currTime;

	// Writing the result
	writeOutputFile(&output, OUTPUT_FILE_NAME);

	print(clusters, fileInput->K);

	// TODO: Call free
	freeRecourses(fileInput, clusters, fileInput->K);

	printf("\nFinished!");

	getchar();
}