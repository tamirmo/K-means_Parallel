#pragma once

#define INPUT_FILE_NAME "input.txt"
#define OUTPUT_FILE_NAME "output.txt"

// TODO: Think of a better place for boolean
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

struct Output {
	// Number of clusters
	int K;
	// Defines the time the algo has stopped at
	double t;
	// The quality of the clusters
	double q;

	// The clusters of the solution
	Cluster* clusters;
};

// Reads the input file (holding all points and parameters for k-means)
Input* readInputFile(const char *fileName);
void writeOutputFile(Output *output, const char* fileName);

void print(Cluster* clusters, int clusterCount);
void print(Input* input);