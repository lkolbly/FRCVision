#ifndef PROCESSING_HXX
#define PROCESSING_HXX

#include <Vector>
#include <opencv2/opencv.hpp>

using namespace std;

#define OBJ_HOLLOW_RECT 0

typedef struct trackingObject_t {
	int type;
	double width;
	double height;
} trackingList_t;

typedef struct pContours_t {
	std::vector<std::vector<cv::Point> > contour_data;
} pContours_t;

class processedImagery_t
{
public:
	// Some sort of image? Maybe?
	char *filename; // Filename of the final image
	int uid;
	cv::Mat img_data;
	
	vector<vector<cv::Point> > contours;
	
	unsigned char *render_contours(unsigned int &len);
	
	vector<trackingObject_t> tracking_targets;
};

void *processingMain(void *arg);

#endif
