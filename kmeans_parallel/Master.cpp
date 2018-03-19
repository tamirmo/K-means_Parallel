#include "Master.h"

Boolean isMasterRank(int rank) {
	Boolean isMasterRank = FALSE;
	if (rank == MASTER_RANK)
		isMasterRank = TRUE;
	return isMasterRank;
}

// Proccess the results goten from slaves at each iteration and returns 
// the index of the first result that is not null or -1 if none finished
int getFirstFinishedIndex(int numOfProcesses, Output* slaveClusters) {
	int finishedSlaveIndex = NO_SLAVE_FINISHED;
	int currProcess;

	for (currProcess = 0; currProcess < numOfProcesses; currProcess++) {
		if (isOutputValid(&(slaveClusters[currProcess]))) {
			finishedSlaveIndex = currProcess;
			// No need to continue, we want the first one
			break;
		}
	}

	return finishedSlaveIndex;
}

void masterPrintResults(Output *result, Cluster* clusters) {
	// Writing the result
	writeOutputFile(result, clusters, OUTPUT_FILE_NAME);
	print(clusters, result->K);
}

// Cheking the output array with results from slaves.
// Returning the index of the first finished slave or NO_SLAVE_FINISHED
int masterCheckResults(InputParams* inputParams,
	Output *result,
	Boolean *isFinished,
	int numOfProcesses,
	Output *slaveOutputs) {
	// Assuming no one finished
	int finishedProcessIndex = NO_SLAVE_FINISHED;

	// If master has finished, no need to process more,
	// his result will be taken
	if (!(*isFinished)) {
		// Checking if there is a process with result
		finishedProcessIndex = getFirstFinishedIndex(numOfProcesses, slaveOutputs);
		if (finishedProcessIndex != NO_SLAVE_FINISHED) {
			*isFinished = TRUE;
			*result = slaveOutputs[finishedProcessIndex];
		}
	}
	else
		finishedProcessIndex = MASTER_RANK;

	return finishedProcessIndex;
}