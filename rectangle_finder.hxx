#ifndef RECTANGLE2_FINDER_HXX
#define RECTANGLE2_FINDER_HXX

#include <opencv2/opencv.hpp>
#include <vector>

using namespace cv;

class Polygon {
public:
  std::vector<Vec4i> m_edges;
  cv::Vec2i m_endpoints[2];

  // Finds the min. distance between our endpoints and the line
  double dist_to_line(cv::Vec4i line);

  // Adds the line to us, attaching the nearest neighbors.
  void add_line(cv::Vec4i line);
  void fromCorners(std::vector<cv::Vec2i> corners);

  int contains(cv::Vec4i line);

  // Find all the corners of this polygon
  std::vector<cv::Vec2i> find_corners(void);
  
  // Find the bound rectangle of this polygon
  cv::Vec4i get_bounds(void);

  // Find the RMS difference from another polygon
  double difference(Polygon *other);

  int should_add_line(cv::Vec4i line);

  std::vector<cv::Vec4i> get_edges(void);
};

class Rectangle3d {
private:
  Polygon m_p;
  cv::Vec2f m_corners[4];
  double m_dist[4];
  cv::Vec2f m_coefs[4];
  
  cv::Vec2f m_camera_FOV;

public:
  Rectangle3d(Polygon p);
  void solve(double w, double h, double fovx, double fovy,
	     double known_w, double known_h);
  cv::Vec2f coef_from_point(double x, double y, double w, double h,
			double fovx, double fovy);
  double find_squareness(cv::Vec2f *coefs, double *dist);

  std::vector<cv::Vec3f> get_points(double *dist);
  std::vector<cv::Vec3f> get_points(void);
  cv::Vec3f get_centroid(void);

  // Returns the distance to the geometric mean
  double distance(void);
  double centroid_dist(void);
  cv::Vec3f normal(void);

  double aspect_ratio(void);
  double azimuth(void);
  double elevation(void);
  double offcenter_angle(void);
  cv::Vec4i get_image_bounds(void);
};

std::vector<Rectangle3d> findRectanglesInImage(Mat img);

#endif
