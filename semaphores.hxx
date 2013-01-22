#ifndef SEMAPHORES_HXX
#define SEMAPHORES_HXX

#include <pthread.h>
#include <opencv2/opencv.hpp>

class processedImagery_t
{
public:
	// Some sort of image? Maybe?
	char *filename; // Filename of the final image
	int uid;
	cv::Mat img_data;
};

typedef struct imageryCollectionConfig_t {
	// Probably camera IP addresses or something here...
	char *image_filename; // The filename of the local file cache
} imageryCollection_t;

typedef struct threadData_t {
	pthread_mutex_t mutex;
	int var;
	
	// Information about the image downloading ("collection") thread
	pthread_mutex_t image_file_lock; // Used by the processing thread
	imageryCollectionConfig_t collection_cfg;
	int has_processed_image;
	
	// Information about the processing thread
	pthread_mutex_t processed_data_lock;
	processedImagery_t processing_result;
	
} threadData_t;

#endif
