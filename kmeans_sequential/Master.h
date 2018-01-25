#pragma once
#include "IOHandler.h"

#define MASTER_RANK 0
#define NO_SLAVE_FINISHED 0

void masterPrintResults(InputParams* inputParams, Output *result);
void masterCheckResults(InputParams* inputParams,
	Output *result,
	Boolean *isFinished,
	int numOfProcesses,
	Output *slaveOutputs);
Boolean isMasterRank(int rank);