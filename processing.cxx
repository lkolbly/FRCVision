#include <cv.h>
#include <stdio.h>
#include <pthread.h>
#include <windows.h>
#include "processing.hxx"
#include "semaphores.hxx"

void *processingMain(void *arg)
{
	threadData_t *td = (threadData_t*)arg;
	printf("Starting processing thread.\n");

	// Wait until the image file has been updated.
	while (1) {
		pthread_mutex_lock(&td->image_file_lock);
		if (!td->has_processed_image) {
			MoveFile("out.jpg", "queue.jpg"); // Copy out the file
			td->has_processed_image = 1;
			pthread_mutex_unlock(&td->image_file_lock);

			printf("Processing file.\n");
			Sleep(1000);
		} else {
			pthread_mutex_unlock(&td->image_file_lock);
		}
		//Sleep(500);
	}

	return 0;
}
