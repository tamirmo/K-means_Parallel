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

void masterPrintResults(InputParams* inputParams, Output *result) {
	// TODO: Print to file "result"
	// Writing the result
	writeOutputFile(result, OUTPUT_FILE_NAME);
	print(result->clusters, inputParams->K);
}

void masterCheckResults(InputParams* inputParams,
	Output *result,
	Boolean *isFinished,
	int numOfProcesses,
	Output *slaveOutputs) {
	// If master has finished, no need to process more,
	// his result will be taken
	if (!isFinished) {
		// Checking if there is a process with result
		int finishedProcessIndex = getFirstFinishedIndex(numOfProcesses, slaveOutputs);
		if (finishedProcessIndex != NO_SLAVE_FINISHED) {
			*isFinished = TRUE;
			*result = slaveOutputs[finishedProcessIndex];
		}
	}
}