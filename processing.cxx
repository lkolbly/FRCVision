#include <stdio.h>
#include <pthread.h>
#include <windows.h>
#include "processing.hxx"
#include "semaphores.hxx"

processedImagery_t processFile(const char *in_fname)
{
	processedImagery_t v;
	v.img_data = cv::imread(in_fname);
	return v;
}

void *processingMain(void *arg)
{
	threadData_t *td = (threadData_t*)arg;
	printf("Starting processing thread.\n");

	// TODO: Multiple images.
	while (1) {
		pthread_mutex_lock(&td->image_file_lock);
		if (td->processing_result.uid != td->collection_cfg.uid) {
			//printf("%i\n", td->processing_result.uid);
			MoveFile("new.jpg", "queue.jpg"); // Copy out the file
			//td->has_processed_image = 1;
			pthread_mutex_unlock(&td->image_file_lock);

			printf("Processing file.\n");
			processedImagery_t processed_imagery = processFile("queue.jpg");
			//Sleep(2000);

			pthread_mutex_lock(&td->processed_data_lock);
			processed_imagery.uid = td->collection_cfg.uid;
			memcpy(&td->processing_result, &processed_imagery, sizeof(processedImagery_t));
			pthread_mutex_unlock(&td->processed_data_lock);
		} else {
			pthread_mutex_unlock(&td->image_file_lock);
		}
		//Sleep(500);
	}

	return 0;
}
