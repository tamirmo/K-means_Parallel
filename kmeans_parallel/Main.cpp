#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
#include <time.h>
#include "OmpKmeans.h"
#include "CudaKmeans.h"
#include "IOHandler.h"
#include "Master.h"

#define INVALID_TIME -1

// Increases each point in the given collection by time with the given velocities
// (for each point x = x + (dt * moment) * vxi ,
// y = y + (dt * moment) * vyi)
void increaseTime(Point* points, Point** gpu_points, int numOfPoints, double dt, int moment) {
	int cudaNumOfPoints = numOfPoints * CUDA_INC_TIME_PORTION, ompNumOfPoints = numOfPoints - cudaNumOfPoints;
	const char* cudaErr = NULL;

	// Starting cuda kernel
	cudaErr = increaseTimeCudaStart(points, cudaNumOfPoints, dt, moment, gpu_points);

	if (cudaErr != NULL){
		printf(cudaErr);
		exit(1);
	}

	// Increasing the left point in omp 
	increaseTimeOmp(points + cudaNumOfPoints, ompNumOfPoints, dt, moment);

	// Wating for cuda to finish the kernel
	cudaErr = increaseTimeCudaEnd(*gpu_points, points, cudaNumOfPoints);

	if (cudaErr != NULL) {
		printf(cudaErr);
		exit(1);
	}
}

void freeRecourses(InputParams* inputParams, 
					Point* points, 
					Cluster* clusters, 
					Output* outputs, 
					Point* pointsCuda) {
	// TODO: Check if all alocations were freed

	free(points);
	free(inputParams);
	free(clusters);
	free(outputs);
	freeCuda(pointsCuda);
}

void calculateClustersCenters(InputParams* inputParams, Point* points, Cluster* clusters) {
	calculateClustersCentersOmp(points, inputParams->N, clusters, inputParams->K);
}

// Assigning all points to clusters
// Returns: TRUE if the point's cluster has changed, FALSE if not
Boolean assignClustersToPoints(InputParams* inputParams, Point* points, Cluster* clusters) {
	Boolean pointChangedOmp = FALSE;
	
	pointChangedOmp = assignClustersToPointsOmp(points, inputParams->N, clusters, inputParams->K);

	return pointChangedOmp;
}

void calculateClustersDiameter(InputParams* inputParams, Point* points, Cluster* clusters) {
	calculateClustersDiameterOmp(points, inputParams->N, clusters, inputParams->K);
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
	// Broadcast all slaves the input read
	MPI_Bcast(inputParams, 1, types->MPI_Input_Param, MASTER_RANK, MPI_COMM_WORLD);
	
	// Master has the points already, slaves need to allocate memory
	if (!isMasterRank(myId)) {
		(*points) = (Point*)calloc((inputParams)->N, sizeof(Point));

		if (*points == NULL) {
			printf("\nAlloc failed points inside broadcastInput");
			fflush(stdout);
			exit(1);
		}
	}

	MPI_Bcast((*points), inputParams->N, types->MPI_Point, MASTER_RANK, MPI_COMM_WORLD);
}

// Getting the finished proccess index and sending the clusters if needed
void handleFinishedProcessIndex(int myId, int finishedProccess, Boolean *isFinished, Cluster** clusters, int clustersCount, MPITypes *types) {
	MPI_Status status;

	// If the result comes from one of slaves,
	// master needs to receive it
	if (!isMasterRank(myId)) {
		printf("\nSlave got id %d", finishedProccess);
		fflush(stdout);
	}

	// If we have the required result
	if (finishedProccess != NO_SLAVE_FINISHED) {
		// If the result comes from one of slaves,
		// master needs to receive it
		if (isMasterRank(myId)) {
			if (!isMasterRank(finishedProccess)) {
				printf("\nMaster receiving from %d", finishedProccess);
				fflush(stdout);
				MPI_Recv(*clusters, clustersCount, types->MPI_Cluster, finishedProccess, 0, MPI_COMM_WORLD, &status);
			}
			else {
				printf("\nMaster result taken from %d", finishedProccess);
				fflush(stdout);
			}
		}
		// Slaves code:
		else {
			// If master wants our result
			if (finishedProccess == myId) {
				printf("\nSlave sending result %d", finishedProccess);
				fflush(stdout);
				MPI_Send(*clusters, clustersCount, types->MPI_Cluster, MASTER_RANK, 0, MPI_COMM_WORLD);
			}
		}
		*isFinished = TRUE;
	}
}

// Handling proccesses left wihtout a time to check at the end
// (we need them to call barrier for the last time before exiting)
void handleLeftProccesses(int myId, 
	int numOfProcesses, 
	Boolean isFinished, 
	InputParams* inputParams, 
	Output *result, 
	int finishedProccess, 
	MPITypes *types) {
	// Calculating the number is processes left t the last round
	// checking a time
	int lastProcessesLeft = (int)(inputParams->T / inputParams->dT) % numOfProcesses;

	// When the prcesses finished it's times
	if (!isFinished && myId >= lastProcessesLeft) {

		// Waiting for all others to finish:
		printf("\nCheckpoint 8.1: id = %d, isFinished = %d, lastProcessesLeft = %d ", myId, isFinished, lastProcessesLeft);
		fflush(stdout);
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Gather(
			result,
			1, types->MPI_Output,
			NULL,
			1, types->MPI_Output,
			MASTER_RANK,
			MPI_COMM_WORLD);
		MPI_Bcast(&finishedProccess, 1, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);
	}
}

// Checking if the current result (q) is good enough (< QM)
// and if so, filling result fields
Boolean isFoundResult(double q, Output *result, InputParams* inputParams, double currTime){
	Boolean isFinished = FALSE;

	// Checking termination condition:
	if (q < inputParams->QM) {

		// Preparing the result and stopping loop:

		isFinished = TRUE;
		(*result).t = currTime;
		(*result).K = inputParams->K;
		(*result).q = q;
	}

	return isFinished;
}

// Communicating with master and checking if someone has finished
// (also sending the result if this proccess has finished)
// Returning the id of the finished proccess
int exchangeProccessesResults(Boolean *isFinished, Output *result, Cluster* clusters,
								Output *slaveOutputs, InputParams* inputParams, 
								int myId, int numOfProcesses,
								MPITypes types) {
	int finishedProccess;
	//printf("\nCheckpoint 4: time = %lf id = %d", MPI_Wtime(), myId);
	//fflush(stdout);

	MPI_Gather(
		result,
		1, types.MPI_Output,
		slaveOutputs,
		1, types.MPI_Output,
		MASTER_RANK,
		MPI_COMM_WORLD);

	//printf("\nCheckpoint 5: id = %d", myId);
	//fflush(stdout);

	if (isMasterRank(myId))
		finishedProccess = masterCheckResults(inputParams,
			result, isFinished,
			numOfProcesses,
			slaveOutputs);

	//printf("\nCheckpoint 6: id = %d", myId);
	//fflush(stdout);

	// Master sends all the index of the finished process
	MPI_Bcast(&finishedProccess, 1, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);

	// Checking if we are finished and if we need to send the clusters
	handleFinishedProcessIndex(myId, finishedProccess, isFinished, &clusters, inputParams->K, &types);

	return finishedProccess;
}

void printResultsAndFree(int myId, Output* result, Cluster* clusters, InputParams* inputParams, Point* points, Point* gpu_points, Output* slaveOutputs) {
	if (isMasterRank(myId))
		masterPrintResults(result, clusters);

	printf("\nCheckpoint 9: id = %d", myId);
	fflush(stdout);

	freeRecourses(inputParams, points, clusters, slaveOutputs, gpu_points);

	stopCuda();

}

// TODO: Think of a name for this function
void parallel_kmeans(int myId, int numOfProcesses) {
	InputParams* inputParams = (InputParams*)malloc(sizeof(InputParams));
	Point* points = NULL, *gpu_points = NULL;;
	Cluster* clusters = NULL;
	Output result = {0,0,0}, *slaveOutputs = NULL;
	int n, finishedProccess;
	Boolean isFinished = FALSE;
	double q, currTime = INVALID_TIME;
	MPITypes types;

	initCuda();

	if (inputParams == NULL) {
		printf("\nAllocation failed parallel_kmeans - inputParams");
		exit(1);
	}

	if (isMasterRank(myId)) {
		slaveOutputs = (Output*)calloc(numOfProcesses, sizeof(Output));
		if (slaveOutputs == NULL) {
			printf("\nAllocation failed parallel_kmeans - slaveOutputs");
			exit(1);
		}
		printf("\nCheckpoint MASTER 1: id = %d", myId);
		fflush(stdout);
		// Reading from file
		readInputFile(INPUT_FILE_NAME, &inputParams, &points);
		if (inputParams == NULL) {
			printf("\n*****FILE ERROR MASTER*****\n");
			fflush(stdout);
			exit(1);
		}
	}

	createMpiTypes(&types);
	broadcastInput(&types, inputParams, &points, myId);

	// Each process starts from t = rank
	increaseTime(points, &gpu_points, inputParams->N, inputParams->dT, myId);

	for (n = myId; n < (inputParams->T / inputParams->dT) && !isFinished; n += numOfProcesses) {
		q = kmeans(inputParams, points, &clusters);
		currTime = n * inputParams->dT;

		printf("\nCheckpoint 3: id = %d n = %d", myId, n);
		fflush(stdout);

		// Checking if the current quality is good enough
		isFinished = isFoundResult(q, &result, inputParams, currTime);
		
		finishedProccess = exchangeProccessesResults(&isFinished, &result, clusters,
									slaveOutputs, inputParams, 
									myId, numOfProcesses, types);

		if (!isFinished)
			// Have not reached the desired quality
			increaseTime(points, &gpu_points, inputParams->N, inputParams->dT, numOfProcesses);

		printf("\nCheckpoint 8: id = %d", myId);
		fflush(stdout);
	}

	// There migth be some proccesses left wihtout a time to check at the end,
	// we need them to call Bcast for the last time before exiting
	handleLeftProccesses(myId, numOfProcesses, 
		isFinished, 
		inputParams, 
		&result, 
		finishedProccess, 
		&types);

	printResultsAndFree(myId, &result, clusters, inputParams, points, gpu_points, slaveOutputs);
}

int main(int argc, char *argv[]) {
	int numOfProcesses, myId;
	clock_t start, end;
	double cpu_time_used;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myId);
	MPI_Comm_size(MPI_COMM_WORLD, &numOfProcesses);

	printf("\nid = %d Started! numOfPr = %d", myId, numOfProcesses);
	fflush(stdout);

	start = clock();

	parallel_kmeans(myId, numOfProcesses);

	end = clock();
	cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;

	printf("\nid = %d Finished! in %lf", myId, cpu_time_used);
	fflush(stdout);

	MPI_Finalize();

	return 0;
}