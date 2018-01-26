#include <stdio.h>
#include <stdlib.h>

#include "IOHandler.h"

// Reads the input file (holding all points and parameters for k-means)
void readInputFile(const char *fileName, InputParams **inputParams, Point **points) {
	int index;
	double coordX, coordY, vx, vy;
	FILE *file;
	// Indicating if points array allocation was successful
	Boolean allocSuccessful = FALSE;

	fopen_s(&file, fileName, "r");
	if (file != NULL) {
		*inputParams = (InputParams*)malloc(sizeof(InputParams));
		if (inputParams != NULL) {
			// First line in the file has many input parameters
			fscanf_s(file, "%d %d %d %lf %d %lf",
				&((*inputParams)->N),
				&((*inputParams)->K),
				&((*inputParams)->T),
				&((*inputParams)->dT),
				&((*inputParams)->LIMIT),
				&((*inputParams)->QM));
			
			(*points) = (Point*)calloc((*inputParams)->N, sizeof(Point));
			if ((*points) != NULL) {
				allocSuccessful = TRUE;
				
				for (index = 0; index < (*inputParams)->N; index++) {
					// Reading each point & velocity and placing it in the array:
					fscanf_s(file, "%lf %lf %lf %lf", &coordX, &coordY, &vx, &vy);
					(*points)[index].position.x = coordX;
					(*points)[index].position.y = coordY;
					(*points)[index].velocity.vx = vx;
					(*points)[index].velocity.vy = vy;
				}
			}
			fclose(file);

			if (!allocSuccessful) {
				// TODO: 
				printf("\nAlloc velocities, points in readInputFile failed");
				fflush(stdout);
			}
		}
		// There is no one process per Point in file
		else {
			// TODO: 
			printf("\nAlloc inputParams,points in readInputFile failed");
			fflush(stdout);
		}
	}
	// Could not read points file
	else {
		// TODO: 
		// TODO: 
		printf("\nRead file failed");
		fflush(stdout);
	}
}

void writeOutputFile(Output *output, Cluster* clusters, const char* fileName) {
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
				clusters[index].center.x,
				clusters[index].center.y);

		fclose(file);
	}
	// Could not open file
	else {
		// TODO: 
	}
}

// Indicating if the given output struct has valid result in it
Boolean isOutputValid(Output *output) {
	Boolean isOutputValid = FALSE;

	if (output->K != INVALID_OUTPUT_K) {
		isOutputValid = TRUE;
	}

	return isOutputValid;
}

// Creates a MPI_Type for Position 
void createPositionType(MPI_Datatype *MPI_Position) {
	Position point;
	MPI_Datatype type[2] = { MPI_DOUBLE, MPI_DOUBLE };
	int blocklen[2] = { 1, 1 };
	MPI_Aint disp[2];

	disp[0] = (char *)&point.x - (char *)&point;
	disp[1] = (char *)&point.y - (char *)&point;
	MPI_Type_create_struct(2, blocklen, disp, type, MPI_Position);
	MPI_Type_commit(MPI_Position);
}

// Creates a MPI_Type for Cluster 
void createClusterType(MPI_Datatype *MPI_Cluster, MPI_Datatype *MPI_Position) {
	Cluster cluster;
	MPI_Datatype type[4] = { MPI_INT, *MPI_Position, MPI_DOUBLE, MPI_INT };
	int blocklen[4] = { 1, 1, 1, 1 };
	MPI_Aint disp[4];

	disp[0] = (char *)&cluster.numOfPoints - (char *)&cluster;
	disp[1] = (char *)&cluster.center - (char *)&cluster;
	disp[2] = (char *)&cluster.diameter - (char *)&cluster;
	disp[3] = (char *)&cluster.id - (char *)&cluster;
	MPI_Type_create_struct(4, blocklen, disp, type, MPI_Cluster);
	MPI_Type_commit(MPI_Cluster);
}

// Creates a MPI_Type for InputParams 
void createInputParamsType(MPI_Datatype *MPI_Input_Param) {
	InputParams inputParam;
	MPI_Datatype type[6] = { MPI_INT, MPI_INT, MPI_INT, MPI_DOUBLE, MPI_INT, MPI_DOUBLE };
	int blocklen[6] = { 1, 1, 1, 1, 1, 1 };
	MPI_Aint disp[6];

	disp[0] = (char *)&inputParam.N - (char *)&inputParam;
	disp[1] = (char *)&inputParam.K - (char *)&inputParam;
	disp[2] = (char *)&inputParam.LIMIT - (char *)&inputParam;
	disp[3] = (char *)&inputParam.QM - (char *)&inputParam;
	disp[4] = (char *)&inputParam.T - (char *)&inputParam;
	disp[5] = (char *)&inputParam.dT - (char *)&inputParam;
	MPI_Type_create_struct(6, blocklen, disp, type, MPI_Input_Param);
	MPI_Type_commit(MPI_Input_Param);
}

// Creates a MPI_Type for Velocity 
void createVelocityType(MPI_Datatype *MPI_Velocity) {
	Velocity velocity;
	MPI_Datatype type[2] = { MPI_DOUBLE, MPI_DOUBLE };
	int blocklen[2] = { 1, 1 };
	MPI_Aint disp[2];

	disp[0] = (char *)&velocity.vx - (char *)&velocity;
	disp[1] = (char *)&velocity.vy - (char *)&velocity;
	MPI_Type_create_struct(2, blocklen, disp, type, MPI_Velocity);
	MPI_Type_commit(MPI_Velocity);
}

void createPointType(MPI_Datatype *MPI_Point, MPI_Datatype *MPI_Velocity, MPI_Datatype *MPI_Position) {
	Point point;
	MPI_Datatype type[3] = { *MPI_Position, MPI_INT, *MPI_Velocity };
	int blocklen[3] = { 1, 1, 1 };
	MPI_Aint disp[3];

	disp[0] = (char *)&(point.position) - (char *)&point;
	disp[1] = (char *)&(point.cluster) - (char *)&point;
	disp[2] = (char *)&(point.velocity) - (char *)&point;
	MPI_Type_create_struct(3, blocklen, disp, type, MPI_Point);
	MPI_Type_commit(MPI_Point);
}

void createOutputType(MPI_Datatype *MPI_Output) {
	Output output;
	MPI_Datatype type[3] = { MPI_INT, MPI_DOUBLE, MPI_DOUBLE };
	int blocklen[3] = { 1, 1, 1};
	MPI_Aint disp[4];

	disp[0] = (char *)&(output.K) - (char *)&output;
	disp[1] = (char *)&(output.t) - (char *)&output;
	disp[2] = (char *)&(output.q) - (char *)&output;
	MPI_Type_create_struct(3, blocklen, disp, type, MPI_Output);
	MPI_Type_commit(MPI_Output);
}

void createMpiTypes(MPITypes *MPITypes) {
	createInputParamsType(&(MPITypes->MPI_Input_Param));
	createPositionType(&MPITypes->MPI_Position);
	createVelocityType(&(MPITypes->MPI_Velocity));
	createClusterType(&(MPITypes->MPI_Cluster), &(MPITypes->MPI_Position));
	createPointType(&(MPITypes->MPI_Point), &(MPITypes->MPI_Velocity), &(MPITypes->MPI_Position));
	createOutputType(&(MPITypes->MPI_Output));
}

void print(InputParams* inputParams, Point* points) {
	printf("N = %d K = %d T =  %d \ndT = %lf LIMIT = %d QM = %lf\n",
		(inputParams->N),
		(inputParams->K),
		(inputParams->T),
		(inputParams->dT),
		(inputParams->LIMIT),
		(inputParams->QM));

	printf("Points:\n");

	for (int i = 0; i < inputParams->N; i++) {
		printf("%lf %lf\t%lf %lf\n",
			points[i].position.x,
			points[i].position.y,
			points[i].velocity.vx,
			points[i].velocity.vy);
	}
	fflush(stdout);
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
	fflush(stdout);
}