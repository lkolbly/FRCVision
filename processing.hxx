#ifndef PROCESSING_HXX
#define PROCESSING_HXX

#include <Vector>
#include <opencv2/opencv.hpp>

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
	
	std::vector<std::vector<cv::Point> > contours;
	
	unsigned char *render_contours(unsigned int &len);
};

void *processingMain(void *arg);

#endif
