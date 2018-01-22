#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define FILE_NAME "input.txt"
#define INITIAL_DISTANCE -1

typedef enum Boolean { FALSE, TRUE } Boolean;

struct Position {
	double x;
	double y;
};

struct Cluster {
	int numOfPoints;
	Position center;
	double diameter;
	int id;
};

struct Velocity {
	double vx;
	double vy;
};

struct Point {
	Position position;
	Cluster* cluster;
	Velocity velocity;
};

struct Input {
	// Number of points
	int N;
	// Number of clusters to find
	int K;
	//  The maximum number of iterations for K - MEAN algorithm
	int LIMIT;
	// Quality measure to stop
	double QM;
	//  Defines the end of time interval[0, T]
	int T;
	// Defines moments t = n*dT, n = { 0, 1, 2, … , T / dT } 
	// for which calculate the clusters and the quality
	double dT;
	Point* points;
	Velocity* velocities;
};

// Reads the input file (holding all points and parameters for k-means)
Input* readInputFile(const char *fileName){
	int index;
	double coordX, coordY, vx, vy;
	FILE *file;
	Input* input = NULL;
	// Indicating if points array allocation was successful
	Boolean allocSuccessful = FALSE;

	fopen_s(&file, fileName, "r");
	if (file != NULL){
		input = (Input*)malloc(sizeof(Input));

		if (input != NULL) {
			// First line in the file has many input parameters
			fscanf_s(file, "%d %d %d %lf %d %lf", 
				&(input->N), 
				&(input->K), 
				&(input->T), 
				&(input->dT), 
				&(input->LIMIT), 
				&(input->QM));

			input->velocities = (Velocity*)calloc(input->N, sizeof(Velocity));
			input->points = (Point*)calloc(input->N, sizeof(Point));

			if (input->velocities != NULL &&
				input->points != NULL) {
				allocSuccessful = TRUE;
				for (index = 0; index < input->N; index++) {
					// Reading each point & velocity and placing it in the array:

					fscanf_s(file, "%lf %lf %lf %lf", &coordX, &coordY, &vx, &vy);
					input->points[index].position.x = coordX;
					input->points[index].position.y = coordY;

					input->velocities[index].vx = vx;
					input->velocities[index].vy = vy;
				}
			}
			fclose(file);

			if (!allocSuccessful){
				// TODO: 
			}
		}
		// There is no one process per Point in file
		else{
			// TODO: 
		}
	}
	// Could not read points file
	else{
		// TODO: 
	}

	return input;
}

void print(Input* input) {
	printf("N = %d K = %d T =  %d \ndT = %lf LIMIT = %d QM = %lf\n",
		(input->N),
		(input->K),
		(input->T),
		(input->dT),
		(input->LIMIT),
		(input->QM));
	
	printf("Points:\n");

	for (int i = 0; i < input->N; i++) {
		printf("%lf %lf\t%lf %lf\n", 
			input->points[i].position.x, 
			input->points[i].position.y, 
			input->velocities[i].vx, 
			input->velocities[i].vy);
	}
}

void print(Cluster* clusters, int clusterCount) {

	printf("Clusters:\n");

	for (int i = 0; i < clusterCount; i++) {
		printf("Id: %d\t(%lf, %lf) num of points = %d\n",
			clusters[i].id,
			clusters[i].center.x,
			clusters[i].center.y,
			clusters[i].numOfPoints);
	}
}

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

void freeRecourses(Input *fileInput, Cluster** clusters, int numOfClusters) {
	int i;
	
	// TODO: Check if all alocations were freed

	for (i = 0; i < fileInput->N; i++) {
		free(fileInput->points);
		free(fileInput->velocities);
	}

	for (i = 0; i < numOfClusters; i++) {
		free(clusters[i]);
	}

	free(fileInput);
	free(clusters);
}

// Assigning K clusters that their centers are the first K point
Cluster* initClusters(Input *fileInput) {
	int i;

	// Assigning an array of K clusters
	Cluster* clusters = (Cluster*)malloc(sizeof(Cluster) * fileInput->K);
	
	if (clusters != NULL) {
		for (i = 0; i < fileInput->K; i++) {
			clusters[i].center.x = fileInput->points[i].position.x;
			clusters[i].center.y = fileInput->points[i].position.y;
			clusters[i].id = i;
			clusters[i].numOfPoints = 0;
		}
	}

	return clusters;
}

// Calculates distance between two points
double getPointsDistance(Position* p1, Position* p2)
{
	return sqrt(pow(p2->x - p1->x, 2) + pow(p2->y - p1->y, 2));
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
	double maxDistance = 0, currDistance;

	for (clusterIndex = 0; clusterIndex < fileInput->K; clusterIndex++) {
		
		// Going over all points and calculating distance with each other point in the cluster
		// to get the maximum distance
		for (i = 0; i < fileInput->N; i++)
			for (j = 0; j < fileInput->N; j++)
				// Calculating distance for points in the same cluster
				if (fileInput->points[i].cluster->id == clusterIndex &&
					fileInput->points[j].cluster->id == clusterIndex) {
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

	for (i = 0; i < fileInput->K; i++) {
		for (j = 0; j < fileInput->K; j++) {
			if (j != i) {
				q += clusters[i].diameter / getPointsDistance(&(clusters[i].center), &(clusters[j].center));
			}
		}
	}

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

		print(clusters, fileInput->K);

		printf("\n\n***** iteration ***** \n\n");
	}

	printf("\nIteration: %d\n\n", iteration);

	return getClustersQuality(fileInput, clusters);
}

void main() {
	Input *fileInput;
	Cluster* clusters;
	double q;

	// Reading from file and testing increase time
	fileInput = readInputFile(FILE_NAME);
	print(fileInput);
	//increaseTime(fileInput->points, fileInput->velocities, fileInput->N, fileInput->dT, 1);
	//printf("\nAfter some time:\n\n");
	print(fileInput);


	clusters = initClusters(fileInput);

	q = kmeans(fileInput, clusters);

	printf("\nq=%lf:\n\n", q);

	print(clusters, fileInput->K);

	increaseTime(fileInput->points, fileInput->velocities, fileInput->N, fileInput->dT, 1);
	printf("\nAfter some time:\n\n");
	print(fileInput);

	q = kmeans(fileInput, clusters);

	printf("\nq=%lf:\n\n", q);

	print(clusters, fileInput->K);

	// TODO: Call free
	//freeRecourses(fileInput, &clusters, );

	getchar();
}