#pragma once

#include <mpi.h>
#include "Point.h"

struct Cluster;

#define INPUT_FILE_NAME "input.txt"
#define OUTPUT_FILE_NAME "output.txt"
#define INVALID_OUTPUT_K 0

typedef enum Boolean { FALSE, TRUE } Boolean;

struct InputParams {
	// Number of points
	int N;
	// Number of clusters to find
	int K;
	//  The maximum number of iterations for K - MEAN algorithm
	int LIMIT;
	// Quality measure to stop
	double QM;
	//  Defines the end of time interval[0, T]
	double T;
	// Defines moments t = n*dT, n = { 0, 1, 2, … , T / dT } 
	// for which calculate the clusters and the quality
	double dT;
};

struct Output {
	// Number of clusters
	int K;
	// Defines the time the algo has stopped at
	double t;
	// The quality of the clusters
	double q;
};

struct MPITypes {
	MPI_Datatype MPI_Velocity;
	MPI_Datatype MPI_Input_Param;
	MPI_Datatype MPI_Cluster;
	MPI_Datatype MPI_Position;
	MPI_Datatype MPI_Output;
	MPI_Datatype MPI_Point;
};

// Reads the input file (holding all points and parameters for k-means)
void readInputFile(const char *fileName, InputParams **inputParams, Point **points);
void writeOutputFile(Output *output, Cluster* clusters, const char* fileName);
// Indicating if the given output struct has valid result in it
Boolean isOutputValid(Output *output);

void print(InputParams* inputParams, Point* points);
void print(Cluster* clusters, int clusterCount);

// Creates MPI_Types:

void createPositionType(MPI_Datatype *MPI_Position);
void createClusterType(MPI_Datatype *MPI_Cluster, MPI_Datatype *MPI_Position);
void createInputParamsType(MPI_Datatype *MPI_Input_Param);
void createVelocityType(MPI_Datatype *MPI_Velocity);
void createPointType(MPI_Datatype *MPI_Point, MPI_Datatype *MPI_Velocity, MPI_Datatype *MPI_Position);
void createMpiTypes(MPITypes *MPITypes);
void createOutputType(MPI_Datatype *MPI_Output);