#include <stdio.h>
#include <pthread.h>
#include <windows.h>
#include "processing.hxx"
#include "semaphores.hxx"

using namespace std;
using namespace cv;

unsigned char *processedImagery_t::render_contours(unsigned int &len_out)
{
	char *data = (char*)malloc(1024);
	data[0] = 0;
	data[1] = 1;
	int offset=4;

	short ncontours = contours.size();
	memcpy(data+2, &ncontours, 2);
	
	for (int i=0; i<contours.size(); i++) {
		//printf("%i: %i has %i points\n", i, hierarchy[i][0], contours[i].size());
		data = (char*)realloc(data, offset+5+4*contours[i].size());
		short npoints = contours[i].size();
		memcpy(data+offset, &npoints, 2);
		offset+=2;
		for (int j=0; j<contours[i].size(); j++) {
			//printf(" - %i/%i: %i %i\n", i,j, contours[i][j].x, contours[i][j].y);
			short x = contours[i][j].x;
			short y = contours[i][j].y;
			memcpy(data+offset, &x, 2);
			memcpy(data+offset, &y, 2);
		}
	}
	
	return (unsigned char *)data;
}

processedImagery_t processFile(const char *in_fname)
{
	processedImagery_t v;
	v.img_data = imread(in_fname);
	
	// Make it grayscale
	Mat gray;
	cvtColor(v.img_data, gray, CV_BGR2GRAY);
	blur(gray, gray, Size(3,3));
	
	// Find contours
	Mat canny_output;
	vector<vector<Point> > contours;
	vector<Vec4i> hierarchy;
	int thresh = 100;
	Canny(gray, canny_output, thresh, thresh*2, 3);
	findContours(canny_output, contours, hierarchy, CV_RETR_TREE, CV_CHAIN_APPROX_SIMPLE, Point(0,0));
	
	v.contours = contours;
	/*
	for (int i=0; i<contours.size(); i++) {
		printf("%i: %i has %i points\n", i, hierarchy[i][0], contours[i].size());
		for (int j=0; j<contours[i].size(); j++) {
			printf(" - %i/%i: %i %i\n", i,j, contours[i][j].x, contours[i][j].y);
		}
	}
	*/
	
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
