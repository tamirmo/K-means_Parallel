#include <stdio.h>
#include <stdlib.h>

#include "IOFilesHandler.h"

// Reads the input file (holding all points and parameters for k-means)
Input* readInputFile(const char *fileName) {
	int index;
	double coordX, coordY, vx, vy;
	FILE *file;
	Input* input = NULL;
	// Indicating if points array allocation was successful
	Boolean allocSuccessful = FALSE;

	fopen_s(&file, fileName, "r");
	if (file != NULL) {
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

			if (!allocSuccessful) {
				// TODO: 
			}
		}
		// There is no one process per Point in file
		else {
			// TODO: 
		}
	}
	// Could not read points file
	else {
		// TODO: 
	}

	return input;
}

void writeOutputFile(Output *output, const char* fileName) {
	int index;
	FILE *file;

	fopen_s(&file, fileName, "w+");
	if (file != NULL) {
		// First line in the file has details of the found clusters
		fprintf_s(file, "First occurrence at t = %lf with q = %lf\n",
			(output->t),
			(output->q));

		fprintf_s(file, "Centers of the clusters:\n");

		for (index = 0; index < output->K; index++)
			fprintf_s(file, "%lf %lf\n",
				output->clusters[index].center.x,
				output->clusters[index].center.y);

		fclose(file);
	}
	// Could not open file
	else {
		// TODO: 
	}
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
		printf("Id: %d\t(%lf, %lf) num of points = %d diam = %lf\n",
			clusters[i].id,
			clusters[i].center.x,
			clusters[i].center.y,
			clusters[i].numOfPoints,
			clusters[i].diameter);
	}
}