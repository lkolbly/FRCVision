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

class Target
{
public:
	double centroid_distance;
	double azimuth, elevation;
	int px_left, px_top, px_right, px_bottom;
};

#if 0
class processedImagery_t
{
public:
#endif
typedef struct {
	// Some sort of image? Maybe?
	char *filename; // Filename of the final image
	int uid;
	cv::Mat img_data;
	
	vector<vector<cv::Point> > contours;
	
	vector<Target> targets;
	
	unsigned char *render_contours(unsigned int &len);
	
	vector<trackingObject_t> tracking_targets;
} processedImagery_t;
#if 0
};
#endif

void *processingMain(void *arg);

#endif
