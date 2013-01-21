#ifndef SEMAPHORES_HXX
#define SEMAPHORES_HXX

#include <pthread.h>

typedef struct processedImagery_t {
	// Some sort of image? Maybe?
} processedImagery_t;

typedef struct threadData_t {
	pthread_mutex_t mutex;
	int var;
	
	// Information about the image downloading thread
	pthread_mutex_t image_file_lock; // Used by the processing thread
	char *image_filename;
	int has_processed_image;
	
	// Information about the processing thread
	pthread_mutex_t processed_data_lock;
	processedImagery_t *processed_imagery;
	
} threadData_t;

#endif
