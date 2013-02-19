#include "rectangle_finder.hxx"
#include <stdio.h>
#include <pthread.h>
#include <windows.h>
#include <algorithm>
#include "processing.hxx"
#include "semaphores.hxx"

using namespace std;
using namespace cv;

int processing_Debug_Counter = 0;

#define USE_DUMMY_TARGET_DATA 0

unsigned char *processedImagery_t::render_contours(unsigned int &len_out)
{
	char *data = new char[1024];
	data[0] = 0;
	data[1] = 1;
	int offset=4;

	short ncontours = contours.size();
	memcpy(data+2, &ncontours, 2);
	
	for (size_t i=0; i<contours.size(); i++) {
		//printf("%i: %i has %i points\n", i, hierarchy[i][0], contours[i].size());
		data = (char*)realloc(data, offset+5+4*contours[i].size());
		short npoints = contours[i].size();
		memcpy(data+offset, &npoints, 2);
		offset+=2;
		for (size_t j=0; j<contours[i].size(); j++) {
			//printf(" - %i/%i: %i %i\n", i,j, contours[i][j].x, contours[i][j].y);
			short x = contours[i][j].x;
			short y = contours[i][j].y;
			memcpy(data+offset, &x, 2);
			memcpy(data+offset, &y, 2);
		}
	}
	
	return (unsigned char *)data;
}

bool sortTargets(Target a, Target b)
{
	Vec2i pa((a.px_left+a.px_right)/2,(a.px_top+a.px_bottom)/2);
	Vec2i pb((b.px_left+b.px_right)/2,(b.px_top+b.px_bottom)/2);
	if (pa[0] > pb[0]) {
		return false;
	}
	return true;
}

processedImagery_t *processFile(const char *in_fname)
{
	//printf("Processing image '%s'\n", in_fname);

	//processedImagery_t *v = (processedImagery_t*)malloc(sizeof(processedImagery_t));
	processedImagery_t *v = new processedImagery_t;
	FILE *f = fopen(in_fname, "rb");
	if (!f) {
		printf("Could not open %s\n", in_fname);
		delete v;
		return NULL;
	}
	fclose(f);
	v->img_data = imread(in_fname);
	if (v->img_data.data == NULL) {
		delete v;
		return NULL;
	}
	
#if USE_DUMMY_TARGET_DATA
	for (int i=0; i<4; i++) {
		Target t;
		t.centroid_distance = 10.0+i;
		t.azimuth = 20.0+i;
		t.elevation = 30.0+i;
		t.px_left = 10+i;
		t.px_top = 20+i;
		t.px_right = 30+i;
		t.px_bottom = 40+i;
		v->targets.push_back(t);
	}
#else
	vector<Rectangle3d> rectangles = findRectanglesInImage(v->img_data);
	for (size_t i=0; i<rectangles.size(); i++) {
		Rectangle3d r = rectangles[i];
		Target t;
		t.centroid_distance = r.centroid_dist();
		t.azimuth = r.azimuth();
		t.elevation = r.elevation();
		t.offcenter_angle = r.offcenter_angle();
		Vec4i bounds = r.get_image_bounds();
		t.px_left = bounds[0];
		t.px_top = bounds[1];
		t.px_right = bounds[2];
		t.px_bottom = bounds[3];
		v->targets.push_back(t);
	}
	//std::sort(v.targets.begin(), v.targets.end(), sortTargets);
	//printf("I have %i targets!\n", v->targets.size());
#endif

	return v;
}

void *processingMain(void *arg)
{
	threadData_t *td = (threadData_t*)arg;
	printf("Starting processing thread.\n");

	// TODO: Multiple images.
	while (1) {
		pthread_mutex_lock(&td->image_file_lock);
		if (!td->processing_result) {
			td->processing_result = new processedImagery_t;
			td->processing_result->uid = td->collection_cfg.uid;
		}
		if (td->processing_result->uid != td->collection_cfg.uid) {
			CopyFile("storage-tmp.jpg", "processing-tmp.jpg", false); // Copy out the file
			pthread_mutex_unlock(&td->image_file_lock);

			//printf("Processing file.\n");
			//Sleep(50);
			processedImagery_t *processed_imagery = processFile("processing-tmp.jpg");
			if (!processed_imagery) {
				continue;
			}
			//Sleep(50);

			pthread_mutex_lock(&td->processed_data_lock);
			processed_imagery->uid = td->collection_cfg.uid;
			td->processing_result->uid = td->collection_cfg.uid;
			//memcpy(&td->processing_result, &processed_imagery, sizeof(processedImagery_t));
			if (td->processing_result) {
				delete td->processing_result;
			}
			td->processing_result = processed_imagery;
			pthread_mutex_unlock(&td->processed_data_lock);
			//printf("It's bloody well unlocked.\n");
		} else {
			pthread_mutex_unlock(&td->image_file_lock);
		}
		//Sleep(500);
	}

	return 0;
}
