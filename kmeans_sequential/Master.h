#pragma once
#include "IOHandler.h"

#define MASTER_RANK 0
#define NO_SLAVE_FINISHED -1

void masterPrintResults(Output *result, Cluster* clusters);
// Cheking the output array with results from slaves.
// Returning the index of the first finished slave or NO_SLAVE_FINISHED
int masterCheckResults(InputParams* inputParams,
	Output *result,
	Boolean *isFinished,
	int numOfProcesses,
	Output *slaveOutputs);
Boolean isMasterRank(int rank);