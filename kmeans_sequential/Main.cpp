#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <mpi.h>
#include "IOHandler.h"
#include "Master.h"

#define INITIAL_DISTANCE -1
#define INVALID_TIME -1

// Increases each point in the given collection by time with the given velocities
// (for each point x = x + (dt * moment) * vxi ,
// y = y + (dt * moment) * vyi)
void increaseTime(Point* points, int numOfPoints, double dt, int moment) {
	double timeInterval = dt * moment;
	for (int i = 0; i < numOfPoints; i++) {
		points[i].position.x += timeInterval * points[i].velocity.vx;
		points[i].position.y += timeInterval * points[i].velocity.vy;
	}
}

void freeRecourses(InputParams* inputParams, 
					Point* points, 
					Cluster* clusters, 
					Output* outputs, 
					int myId) {
	// TODO: Check if all alocations were freed

	free(points);
	if (isMasterRank(myId))
	{
		printf("\nfree1");
		fflush(stdout);
	}
	free(inputParams);
	if (isMasterRank(myId))
	{
		printf("\nfree2");
		fflush(stdout);
	}
	free(clusters);
	if (isMasterRank(myId))
	{
		printf("\nfree3");
		printf("\nk = %d", outputs[0].K);
		printf("\nq = %lf", outputs[0].q);
		printf("\nt = %lf", outputs[0].t);
		fflush(stdout);
	}
	free(outputs);
	if (isMasterRank(myId))
	{
		printf("\nfree4");
		fflush(stdout);
	}
}

// Assigning K clusters that their centers are the first K point
void initClusters(InputParams* inputParams, Point* points, Cluster** clusters) {
	int i;

	if(*clusters == NULL)
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

void calculateClustersCenters(InputParams* inputParams, Point* points, Cluster* clusters) {
	int clusterIndex, pointIndex;
	double sumX, sumY;

	for (clusterIndex = 0; clusterIndex < inputParams->K; clusterIndex++) {
		// When a cluster has no points, we keep it's center intact
		if (clusters[clusterIndex].numOfPoints != 0 ) {
			sumX = sumY = 0;

			for (pointIndex = 0; pointIndex < inputParams->N; pointIndex++) {
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
Boolean assignClustersToPoints(InputParams* inputParams, Point* points, Cluster* clusters) {
	int pointIndex;
	// Indicating if at least one point has changed cluster
	Boolean pointChanged = FALSE;

	for (pointIndex = 0; pointIndex < inputParams->N; pointIndex++)
		// Assigning the curr point and checking if changed cluster
		if (assignClusterToPoint(&(points[pointIndex]), clusters, inputParams->K))
			pointChanged = TRUE;

	return pointChanged;
}

void calculateClustersDiameter(InputParams* inputParams, Point* points, Cluster* clusters) {
	int clusterIndex, i, j;
	double maxDistance, currDistance;

	for (clusterIndex = 0; clusterIndex < inputParams->K; clusterIndex++) {
		maxDistance = 0;

		// Going over all points and calculating distance with each other point in the cluster
		// to get the maximum distance
		for (i = 0; i < inputParams->N - 1; i++)
			if (points[i].cluster->id == clusterIndex)
				for (j = i+1; j < inputParams->N; j++)
					// Calculating distance for points in the same cluster
					if(points[j].cluster->id == clusterIndex) {
						currDistance = getPointsDistance(&(points[i].position), &(points[j].position));
						if (currDistance > maxDistance)
							maxDistance = currDistance;
					}

		clusters[clusterIndex].diameter = maxDistance;
	}
}

// The clusters quality is an average of diameters of the cluster divided by distance to other clusters
double getClustersQuality(InputParams* inputParams, Point* points, Cluster* clusters) {
	double q = 0;
	int i, j;

	calculateClustersDiameter(inputParams, points, clusters);

	for (i = 0; i < inputParams->K; i++)
		for (j = 0; j < inputParams->K; j++)
			if (j != i) 
				q += clusters[i].diameter / getPointsDistance(&(clusters[i].center), &(clusters[j].center));

	// Devide q by the number of addends
	q /= (inputParams->K-1) * (inputParams->K);

	return q;
}

// Applying kmeans algo for the given input
// Returns: The quality of the clusters calculated
double kmeans(InputParams* inputParams, Point* points, Cluster** clusters) {
	int iteration;
	// Indicating if at least one point changed cluster
	// (initializing to true for the first round)
	Boolean pointChanged = TRUE;

	initClusters(inputParams, points, clusters);

	// Assigning clusters until no point has changed
	// or the limit iterations count has reached
	for (iteration = 0; iteration < inputParams->LIMIT && pointChanged; iteration++) {
		pointChanged = assignClustersToPoints(inputParams, points, *clusters);
		
		// Updating the new centers
		calculateClustersCenters(inputParams, points, *clusters);
	}

	return getClustersQuality(inputParams, points, *clusters);
}

void broadcastInput(MPITypes *types, InputParams* inputParams, Point** points, int myId) {
	// TODO: Broadcast all slaves the array
	MPI_Bcast(inputParams, 1, types->MPI_Input_Param, MASTER_RANK, MPI_COMM_WORLD);
	printf("\nCheckpoint 1.2 ");
	fflush(stdout);
	
	// Master has the points already, slaves need to allocate memory
	if (!isMasterRank(myId)) {
		(*points) = (Point*)calloc((inputParams)->N, sizeof(Point));

		if (*points == NULL) {
			// TODO: Alloc failed
			printf("\nAlloc failed points inside broadcastInput");
			fflush(stdout);
		}
	}

	MPI_Bcast((*points), inputParams->N, types->MPI_Point, MASTER_RANK, MPI_COMM_WORLD);
	printf("\nCheckpoint 1.3 ");
	fflush(stdout);
}

// TODO: Think of a name for this function
void parallel_kmeans(int myId, int numOfProcesses) {
	InputParams* inputParams = (InputParams*)malloc(sizeof(InputParams));
	Point* points = NULL;
	Cluster* myClusters = NULL;
	Output result, *slaveOutputs = NULL;
	int n;
	Boolean isFinished = FALSE;
	double q, currTime = INVALID_TIME;
	MPITypes types;
	
	printf("\nCheckpoint 1 : id = %d", myId);
	fflush(stdout);

	if (inputParams == NULL) {
		// TODO: Alloc failed
		printf("\nCheckpoint 1 : id = %d", myId);
		fflush(stdout);
	}

	if (isMasterRank(myId)) {
		slaveOutputs = (Output*)malloc(sizeof(Output) * numOfProcesses);
		if (slaveOutputs == NULL) {
			// TODO: Alloc failed
		}
		printf("\nCheckpoint MASTER 1: id = %d", myId);
		fflush(stdout);
		// Reading from file
		readInputFile(INPUT_FILE_NAME, &inputParams, &points);
		if (inputParams == NULL) {
			printf("\nFILE ERROR MASTER 1: id = %d", myId);
			fflush(stdout);
		}
	}
	printf("\nCheckpoint 1.0 : id = %d", myId);
	createMpiTypes(&types);
	printf("\nCheckpoint 1.1 : id = %d", myId);
	fflush(stdout);
	broadcastInput(&types, inputParams, &points, myId);
	createOutputType(&(types.MPI_Output), &(types.MPI_Cluster), inputParams->K);

	printf("\nCheckpoint 2: id = %d", myId);
	fflush(stdout);

	if (!isMasterRank(myId)) {
		print(inputParams, points);
	}

	for (n = myId; n < (inputParams->T / inputParams->dT) && !isFinished; n += numOfProcesses) {
		q = kmeans(inputParams, points, &myClusters);
		currTime = n * numOfProcesses * inputParams->dT;

		printf("\nCheckpoint 3: id = %d", myId);
		fflush(stdout);

		// Checking termination condition:
		if (q < inputParams->QM) {
			printf("\nFound QM: id = %d", myId);
			fflush(stdout);

			// Preparing the result and stopping loop:

			isFinished = TRUE;
			result.t = currTime;
			result.K = inputParams->K;
			result.q = q;
			result.clusters = myClusters;
		}

		// Waiting till all hosts finish. No point in moving to next time, 
		// other proc migth finished
		MPI_Barrier(MPI_COMM_WORLD);
		
		printf("\nCheckpoint 4: id = %d", myId);
		fflush(stdout);

		MPI_Gather(
			&result, 
			1, types.MPI_Output,
			slaveOutputs,
			1, types.MPI_Output,
			MASTER_RANK,
			MPI_COMM_WORLD);

		printf("\nCheckpoint 5: id = %d", myId);
		fflush(stdout);

		if (isMasterRank(myId))
			/*masterCheckResults(inputParams,
				&result, &isFinished,
				numOfProcesses,
				slaveOutputs);*/

		printf("\nCheckpoint 6: id = %d", myId);
		fflush(stdout);

		MPI_Bcast(&isFinished, 1, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);
			
		printf("\nCheckpoint 7: id = %d isFinished = %d", myId, isFinished);
		fflush(stdout);

		if (!isFinished)
			// Have not reached the desired quality
			increaseTime(points, inputParams->N, inputParams->dT, numOfProcesses);

		printf("\nCheckpoint 8: id = %d", myId);
		fflush(stdout);
	}

	if (isMasterRank(myId))
		masterPrintResults(inputParams, &result);
	
	printf("\nCheckpoint 9: id = %d", myId);
	fflush(stdout);

	// TODO: Call free
	freeRecourses(inputParams,
		points,
		myClusters,
		slaveOutputs,
		myId);
}

int main(int argc, char *argv[]) {
	int numOfProcesses, myId;

	printf("\nbefore init!");
	fflush(stdout);

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myId);
	MPI_Comm_size(MPI_COMM_WORLD, &numOfProcesses);

	printf("\nid = %d Started!", myId);
	fflush(stdout);

	parallel_kmeans(myId, numOfProcesses);

	printf("\nid = %d Finished!", myId);
	fflush(stdout);

	MPI_Finalize();

	return 0;
}