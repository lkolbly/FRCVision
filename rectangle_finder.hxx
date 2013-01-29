#ifndef RECTANGLE_FINDER_HXX
#define RECTANGLE_FINDER_HXX

#include <opencv/cv.h>
#include <vector>

using namespace std;
using namespace cv;

class Polygon {
public:
  vector<Vec4i> m_edges;
  Vec2i m_endpoints[2];

  // Finds the min. distance between our endpoints and the line
  double dist_to_line(Vec4i line);

  // Adds the line to us, attaching the nearest neighbors.
  void add_line(Vec4i line);

  int contains(Vec4i line);

  // Find all the corners of this polygon
  vector<Vec2i> find_corners(void);

  // Find the RMS difference from another polygon
  double difference(Polygon other);

  int should_add_line(Vec4i line);

  vector<Vec4i> get_edges(void);
};

class Rectangle3d {
private:
  Polygon m_p;
  Vec2f m_corners[4];
  double m_dist[4];
  Vec2f m_coefs[4];

public:
  Rectangle3d(Polygon p);
  void solve(double w, double h, double fovx, double fovy,
	     double known_w, double known_h);
  Vec2f coef_from_point(double x, double y, double w, double h,
			double fovx, double fovy);
  double find_squareness(Vec2f *coefs, double *dist);

  vector<Vec3f> get_points(double *dist);
  vector<Vec3f> get_points(void);

  // Returns the distance to the geometric mean
  double distance(void);
  double centroid_dist(void);

  double aspect_ratio(void);
};

vector<Rectangle3d> findRectanglesInImage(Mat img);

#endif
