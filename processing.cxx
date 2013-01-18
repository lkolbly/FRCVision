#include <cv.h>
#include <stdio.h>
#include <pthread.h>
#include <windows.h>
#include "processing.hxx"
#include "networking.hxx"

void *processingMain(void *arg)
{
	threadData_t *td = (threadData_t*)arg;
	printf("Starting processing thread.\n");
	
	int last_processed = 0;

	// Wait until the image file has been updated.
	while (1) {
		pthread_mutex_lock(&td->mutex);
		//printf("Processing: %i/%i\n", last_processed, td->var);
		if (last_processed < td->var) {
			printf("Processing file.\n");
			Sleep(1000);
			last_processed = ++td->var;
		}
		pthread_mutex_unlock(&td->mutex);
		//Sleep(500);
	}

	return 0;
}
