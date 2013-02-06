/* This is a standalone program. Pass an image name as a first parameter
of the program.  Switch between standard and probabilistic Hough transform
by changing "#if 1" to "#if 0" and back */
// all angles are stored in degrees
#include <opencv/cv.h>
#include <opencv2/highgui/highgui.hpp>
#include <math.h>
#include <iostream>
#include <unistd.h>
#include <stdio.h>
#include <vector>
#include <algorithm>
#include "rectangle_finder.hxx"

#define DEG2RAD(x) ((x)/57.2957795)
#define RAD2DEG(x) ((x)*57.2957795)

using namespace cv;

double min(double a, double b)
{
  return (a<b)?a:b;
}

double pnt_dist(double x1, double y1, double x2, double y2)
{
  double d = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
  //printf("%f,%f,%f,%f => %f\n", x1,y1,x2,y2, d);
  return d;
}

double vec_dist(Vec2f v1, Vec2f v2)
{
	return pnt_dist(v1[0], v1[1], v2[0], v2[1]);
}

#if 0
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
  double difference(Polygon *other);

  int should_add_line(Vec4i line);

  vector<Vec4i> get_edges(void);
};
#endif

vector<Vec4i> Polygon::get_edges(void) {
	return m_edges;
}

void Polygon::fromCorners(std::vector<cv::Vec2i> corners)
{
	for (int i=0; i<corners.size(); i++) {
		m_edges.push_back(Vec4i(corners[i][0], corners[i][1], corners[(i+1)%corners.size()][0], corners[(i+1)%corners.size()][1]));
	}
}

double Polygon::difference(Polygon *other)
{
  double sum = 0.0;
  for (size_t i=0; i<m_edges.size(); i++) {
    // Find the closest point on the other polygon
    double min_diff = 10000.0;
    for (size_t j=0; j<other->get_edges().size(); j++) {
      double diff = pnt_dist(m_edges[i][0], m_edges[i][1], other->get_edges()[j][0], other->get_edges()[j][1]);
      if (diff < min_diff) {
		min_diff = diff;
      }
    }

    sum += min_diff*min_diff;
  }

  return sqrt(1.0 / m_edges.size() * sum);
}

double Polygon::dist_to_line(Vec4i line)
{
  double d1 = pnt_dist(line[0], line[1], m_endpoints[0][0], m_endpoints[0][1]);
  double d2 = pnt_dist(line[0], line[1], m_endpoints[1][0], m_endpoints[1][1]);
  double d3 = pnt_dist(line[2], line[3], m_endpoints[0][0], m_endpoints[0][1]);
  double d4 = pnt_dist(line[2], line[3], m_endpoints[1][0], m_endpoints[1][1]);
  return min(min(d1, d2), min(d3, d4));
}

int min_index(double d1, double d2, double d3, double d4)
{
  if (d1 < d2 && d1 < d3 && d1 < d4) {
    return 1;
  } else if (d2 < d3 && d2 < d4) {
    return 2;
  } else if (d3 < d4) {
    return 3;
  }
  return 4;
}

void Polygon::add_line(Vec4i line)
{
  double d1 = pnt_dist(line[0], line[1], m_endpoints[0][0], m_endpoints[0][1]);
  double d2 = pnt_dist(line[0], line[1], m_endpoints[1][0], m_endpoints[1][1]);
  double d3 = pnt_dist(line[2], line[3], m_endpoints[0][0], m_endpoints[0][1]);
  double d4 = pnt_dist(line[2], line[3], m_endpoints[1][0], m_endpoints[1][1]);
  int m = min_index(d1,d2,d3,d4);//min(min(d1, d2), min(d3, d4));
  if (m == 1) {
    Vec4i l2;
    l2[2] = line[0];
    l2[3] = line[1];
    l2[0] = line[2];
    l2[1] = line[3];
    m_edges.insert(m_edges.begin(), l2);
  } else if (m == 2) {
    m_edges.push_back(line);
  } else if (m == 3) {
    m_edges.insert(m_edges.begin(), line);
  } else {
    Vec4i l2;
    l2[2] = line[0];
    l2[3] = line[1];
    l2[0] = line[2];
    l2[1] = line[3];
    m_edges.push_back(l2);
  }
  m_endpoints[0][0] = m_edges[0][0];
  m_endpoints[0][1] = m_edges[0][1];
  m_endpoints[1][0] = m_edges[m_edges.size()-1][2];
  m_endpoints[1][1] = m_edges[m_edges.size()-1][3];
}

int Polygon::contains(Vec4i line)
{
  for (size_t i=0; i<m_edges.size(); i++) {
    Vec4i l = m_edges[i];
    //printf("%i %i %i %i\n", l[0], line[0], l[2], line[2]);
    if (l[0] == line[0] && l[2] == line[2]) {
      return 1;
    }
    if (l[0] == line[2] && l[2] == line[0]) {
      return 1;
    }
  }
  return 0;
}

vector<Vec2i> Polygon::find_corners(void)
{
  vector<Vec2i> corners;
  for (size_t i=0; i<m_edges.size(); i++) {
    Vec2i c;
    c[0] = m_edges[i][0];
    c[1] = m_edges[i][1];
    corners.push_back(c);
  }
  return corners;
}

double hypot(double x,double y)
{
    double t;
    x = abs(x);
    y = abs(y);
    t = min(x,y);
    x = max(x,y);
    t = t/x;
    return x*sqrt(1+t*t);
}

double dist_pnt_to_line(Vec2i p, Vec4i l)
{
  // Thank you, Wikipedia!
  double normal_len = hypot(l[2]-l[0], l[3]-l[1]);
  double d = fabs((p[0]-l[0])*(l[3]-l[1]) - (p[1]-l[1])*(l[2]-l[0])) / normal_len;
  return d;
}

// This function determines if said line could belong to us based on two things:
// - If the distance from one of our endpoints to the line is < 5, then
//   it could be part of us.
// - If they are within 10 degrees of being parallel, and the distance is <20,
//   and the angle off of the line to the new endpoint is <30 degrees.
int Polygon::should_add_line(Vec4i line)
{
  double d = dist_to_line(line);
  if (d < 20) {
    return 1;
  }

  return 0;
}

Vec4i Polygon::get_bounds(void)
{
	int left=1000, right=0, top=0, bottom=1000;
	for (int i=0; i<m_edges.size(); i++) {
		printf("Edge: %i %i %i %i\n", m_edges[i][0], m_edges[i][1], m_edges[i][2], m_edges[i][3]);
		if (m_edges[i][0] < left) left = m_edges[i][0];
		if (m_edges[i][0] > right) right = m_edges[i][0];

		if (m_edges[i][2] < left) left = m_edges[i][2];
		if (m_edges[i][2] > right) right = m_edges[i][2];

		if (m_edges[i][1] < bottom) bottom = m_edges[i][1];
		if (m_edges[i][1] > top) top = m_edges[i][1];

		if (m_edges[i][3] < bottom) bottom = m_edges[i][3];
		if (m_edges[i][3] > top) top = m_edges[i][3];
	}
	return Vec4i(left, right, bottom, top);
}

vector<Polygon> detectRectangles(vector<Vec4i> lines)
{
  // We iterate through every line. For the beginning and ending points,
  // we find nearby lines, and generate shapes that match.
  // After discarding shapes that don't have 4 sides, we are left with
  // only quadrangles.
  vector<Polygon> rectangles;
  //printf("We're detecting rectangles in %lu lines.\n", lines.size());
  for (size_t i=0; i<lines.size(); i++) {
    Polygon p;
    p.add_line(lines[i]);
    for (size_t j=0; j<lines.size(); j++) {
      if (!p.contains(lines[j])) {
	//double d = p.dist_to_line(lines[j]);
	//printf("%i (%i,%i,%i,%i) has distance %f\n", j, lines[j][0], lines[j][1], lines[j][2], lines[j][3], d);
	//if (d < 5) {
	if (p.should_add_line(lines[j])) {
	  p.add_line(lines[j]);
	  j = 0;
	}

#if 0 // We don't actually do this check. But we should.
	// Check to see if the endpoint is nearby one of our endpoint lines.
	// If so, then we should shorten said endpoint line.
	double d=dist_pnt_to_line(Vec2i(lines[j][0], lines[j][1]), m_edges[0]);
	if (d < 5) {
	  // Add it!
	  m_edges[0][0] = lines[j][0];
	  m_edges[0][1] = lines[j][1];
	}
#endif
      }
    }
    //printf("p has %lu edges.\n", p.m_edges.size());
    if (p.m_edges.size() == 4) {
      // Find the corners
      vector<Vec2i> corners = p.find_corners();
      //printf(" - ");
      for (size_t j=0; j<corners.size(); j++) {
	//printf("%i,%i ", corners[j][0], corners[j][1]);
      }
      //printf("\n");

      rectangles.push_back(p);
    }
  }

  // Now collapse all of the rectangles
  vector<Polygon> unique_rectangles;
  for (size_t i=0; i<rectangles.size(); i++) {
    double min_diff = 10000.0;
    for (int j=i-1; j>=0; j--) {
	double diff = rectangles[i].difference(&rectangles[j]);
	if (diff < min_diff) {
	  min_diff = diff;
	}
	//printf("%lu %i %f\n", i, j, diff);
    }

    // If we're the first one like this one that we've seen...
   //printf("%lu %f\n", i, min_diff);
    if (min_diff > 2.0) {
      unique_rectangles.push_back(rectangles[i]);
    }
  }

  // Print out the 'unique' rectangles
 //printf("Printing out unique rectangles...\n");
  for (size_t i=0; i<unique_rectangles.size(); i++) {
    vector<Vec2i> corners = unique_rectangles[i].find_corners();
   //printf(" - ");
    for (size_t j=0; j<corners.size(); j++) {
     //printf("%i,%i ", corners[j][0], corners[j][1]);
    }
   //printf("\n");
  }

  return unique_rectangles;
}

#if 0
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
#endif

Rectangle3d::Rectangle3d(Polygon p)
{
  m_p = p;
  for (int i=0; i<4; i++) {
    m_corners[i][0] = m_p.m_edges[i][0];
    m_corners[i][1] = m_p.m_edges[i][1];
  }
}

Vec2f Rectangle3d::coef_from_point(double x, double y, double w, double h,
				   double fovx, double fovy)
{
  Vec2f v;
  v[0] = tan(DEG2RAD(x/w*fovx - fovx/2.0));
  v[1] = tan(DEG2RAD(y/h*fovy - fovy/2.0));
  //printf("%f,%f => %f,%f (in %f,%f from %f)\n", x,y, v[0], v[1], w,h, DEG2RAD(x/w*fovx - fovx/2.0));
  return v;
}

Vec3f my_normalize(Vec3f v)
{
  double len = sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
  return v;
}

// Find the angle at b by a and c.
double angle3d(Vec3f a, Vec3f b, Vec3f c)
{
  //Vec3f v1 = //my_normalize(a - b);
  //Vec3f v2 = //my_normalize(c - b);
  Vec3f v1;
  Vec3f v2;
  v1[0] = a[0] - b[0];
  v1[1] = a[1] - b[1];
  v1[2] = a[2] - b[2];
  v2[0] = c[0] - b[0];
  v2[1] = c[1] - b[1];
  v2[2] = c[2] - b[2];
  v1 = my_normalize(v1);
  v2 = my_normalize(v2);
  double dp = v1.dot(v2);
  //printf("%f,%f,%f %f,%f,%f => %f\n", v1[0],v1[1],v1[2], v2[0],v2[1],v2[2],dp);
  return RAD2DEG(acos(dp));
}

double dist3d(Vec3f a, Vec3f b)
{
  return sqrt((b[0]-a[0])*(b[0]-a[0]) + (b[1]-a[1])*(b[1]-a[1]));
}

vector<Vec3f> Rectangle3d::get_points(double *dist)
{
  vector<Vec3f> pnts;
  for (int i=0; i<4; i++) {
    Vec3f p;
    pnts.push_back(p);
    if (dist[i] < 0.0) dist[i] = 0.0;
    pnts[i][0] = m_coefs[i][0] * dist[i];
    pnts[i][1] = m_coefs[i][1] * dist[i];
    pnts[i][2] = dist[i];
  }
  return pnts;
}

vector<Vec3f> Rectangle3d::get_points(void)
{
  return get_points(m_dist);
}

double Rectangle3d::find_squareness(Vec2f *coefs, double *dist)
{
  // Find the current trial points
#if 0
  Vec3f pnts[4];
  for (int i=0; i<4; i++) {
    if (dist[i] < 0.0) dist[i] = 0.0;
    pnts[i][0] = coefs[i][0] * dist[i];
    pnts[i][1] = coefs[i][1] * dist[i];
    pnts[i][2] = dist[i];
  }
#endif
  vector<Vec3f> pnts = get_points(dist);

  // Find the current angles
  double a[4];
  a[0] = angle3d(pnts[3], pnts[0], pnts[1]);
  a[1] = angle3d(pnts[0], pnts[1], pnts[2]);
  a[2] = angle3d(pnts[1], pnts[2], pnts[3]);
  a[3] = angle3d(pnts[2], pnts[3], pnts[0]);

  // Figure out the MS inaccuracy...
  double sum = 0.0;
  for (int i=0; i<4; i++) {
    sum += (a[i]-90.0)*(a[i]-90.0);
    //printf("%f => %f\n", dist[i], a[i]);
  }
  for (int i=0; i<4; i++) {
    //printf("%f,%f,%f\n", pnts[i][0], pnts[i][1], pnts[i][2]);
  }
  sum = sum / 4.0;
  //printf("Mean Square inaccuracy = %f\n", sum);

  return sum;
}

void Rectangle3d::solve(double w, double h, double fovx, double fovy,
			double known_w, double known_h)
{
  m_camera_FOV[0] = fovx;
  m_camera_FOV[1] = fovy;
  //std::cout << "Solving a rectangle...";
  double k = 0.01; // This is used in the iterative solver.

  // Find the four lines that extend to the corners from the camera
  Vec2f coefs[4];
  for (int i=0; i<4; i++) {
    m_coefs[i] = coef_from_point((double)m_corners[i][0], (double)m_corners[i][1],
			       w, h, fovx, fovy);
  }

  // Solve to figure out the distance along said lines to make the angles 90deg
  // Yes, it's an iterative solver. Sorry, Ms. Harrelson.
  double dist[4];
  dist[0] = 1.0; // Fix one of the distances for now.
  dist[1] = dist[2] = dist[3] = 1.0;
  double last_inaccuracy = 0.0;
  while (1) {
    double sum = find_squareness(coefs, dist);
    printf("Inaccuracy: %f (k=%f)\n", sum, k);
    if (last_inaccuracy-sum < 0.001) {
      break;
    }
	if (k < 0.000001) {
		break;
	}
    //usleep(1000000);

    // Figet with the distances, and go with the optimal one.
    int did_change = 0;
    for (int i=0; i<4; i++) {
      double orig = dist[i];
      dist[i] = orig + k;
      double incr = find_squareness(coefs, dist);
      dist[i] = orig - k;
      double decr = find_squareness(coefs, dist);
      dist[i] = orig;
      if (sum-incr>sum-decr && incr<sum) {
	dist[i] += k;
	sum = incr;
	did_change = 1;
      } else if (decr < sum) {
	dist[i] -= k;
	sum = decr;
	did_change = 1;
      }
    }
    if (!did_change) {
      k = k/2;
    }

    // Eventually we break, when we're close enough...
    if (sum < 0.01) {
      break;
    }

    // Adjust the distances based on the current angles
#if 0
    double overall_diff = 0.0;
    if (a[0] > 90.0) {
      overall_diff = k*(a[0]-90.0)/2.0;
    } else {
      overall_diff = -k*(90.0-a[0])/2.0;
    }
    for (int i=1; i<4; i++) { // Remember: point 0 is fixed.
     //printf("%f => %f => %f\n", dist[i], a[i], k*(a[i]-90.0)+overall_diff);
      if (a[i] > 90.0) {
	dist[i] -= k*(a[i]-90.0);
      } else {
	dist[i] += k*(90.0-a[i]);
      }
      dist[i] += overall_diff;
    }
   //printf("\n");
#endif
  }

  // Copy out the distances.
  for (int i=0; i<4; i++)
    m_dist[i] = dist[i];

  // Now we can scale everyone to match the width to the known width
  // We'll assume the 'width' is the 'longest' edge.
  while (1) {
    // Find the longest edge
    double est_w = 0.0;
    vector<Vec3f> pnts = get_points();
    for (int i=0; i<4; i++) {
      double d = dist3d(pnts[i], pnts[(i+1)%3]);
      if (d > est_w) {
	est_w = d;
      }
    }

    // Scale all of the distances by the needed correction factor
    double factor = known_w / est_w;
    for (int i=0; i<4; i++) {
      m_dist[i] *= factor;
    }

    break;
  }

  // Some debugging info
 //printf("Adjusted distances to %f,%f,%f,%f\n", m_dist[0], m_dist[1], m_dist[2], m_dist[3]);
    double est_w = 0.0;
    double est_h = 0.0;
    vector<Vec3f> pnts = get_points();
    for (int i=0; i<4; i++) {
      double d = dist3d(pnts[i], pnts[(i+1)%3]);
      if (d > est_w) {
	est_w = d;
	est_h = dist3d(pnts[(i+1)%3], pnts[(i+2)%3]);
      }
    }
   //printf("Final estimated size: %fx%f, ar=%f\n", est_w, est_h, est_w/est_h);
}

double Rectangle3d::aspect_ratio(void)
{
  double est_w = 0.0;
  double est_h = 0.0;
  vector<Vec3f> pnts = get_points();
  for (int i=0; i<4; i++) {
    double d = dist3d(pnts[i], pnts[(i+1)%3]);
    if (d > est_w) {
      est_w = d;
      est_h = dist3d(pnts[(i+1)%3], pnts[(i+2)%3]);
    }
  }
  return est_w / est_h;
}

Vec3f Rectangle3d::get_centroid(void)
{
  vector<Vec3f> pnts = get_points();
  Vec3f centroid(0,0,0);
  for (int i=0; i<4; i++) {
    centroid += pnts[i];
  }
  centroid[0] = centroid[0] / 4.0;
  centroid[1] = centroid[1] / 4.0;
  centroid[2] = centroid[2] / 4.0;
  return centroid;
}

double Rectangle3d::centroid_dist(void)
{
  Vec3f centroid = get_centroid();
  return sqrt(centroid[2]*centroid[2]+centroid[1]*centroid[1]+centroid[0]*centroid[0]);
}

double Rectangle3d::azimuth(void)
{
	Vec3f p = get_centroid();
	return m_camera_FOV[0] * p[0] / (p[2] * sin(DEG2RAD(m_camera_FOV[0])));
}

double Rectangle3d::elevation(void)
{
	Vec3f p = get_centroid();
	return m_camera_FOV[1] * p[1] / (p[2] * sin(DEG2RAD(m_camera_FOV[1])));
}

Vec4i Rectangle3d::get_image_bounds(void)
{
	return m_p.get_bounds();
}

int rectangle_Finder_Debug_Counter = 0;

void outputPicture(const char *name, Mat img)
{
	char buf[1024];
	snprintf(buf, 1024, "debug/%s_%04i.jpg", name, rectangle_Finder_Debug_Counter);
	imwrite(buf, img);
}

double getReferenceAngle(Vec2f a, Vec2f b) {
	return RAD2DEG(atan2(a[1]-b[1], a[0]-b[0]));
}

class HullVertexComparator
{
public:
	Vec2f m_centroid;
	bool operator() (Vec2f a, Vec2f b) {
		double a1 = RAD2DEG(atan2(a[1]-m_centroid[1], a[0]-m_centroid[0]));
		double a2 = RAD2DEG(atan2(b[1]-m_centroid[1], b[0]-m_centroid[0]));
		printf("Centroid=<%f,%f> a=<%f,%f>=>%f b=<%f,%f>=>%f\n", m_centroid[0], m_centroid[1], a[0], a[1], a1, b[0], b[1], a2);
		return a1>a2;
	};
};

Vec2f normalize(Vec2f v)
{
	double l = sqrt(v[0]*v[0] + v[1]*v[1]);
	Vec2f v2(v[0]/l, v[1]/l);
	return v2;
}

Mat myHaughTransform(vector<Vec2f> hull, const char *img_name)
{
	printf("Haugh transforming a hull of size %i\n", hull.size());

	Vec2f centroid(0,0);
	for (int j=0; j<hull.size(); j++) {
		centroid[0] += hull[j][0];
		centroid[1] += hull[j][1];
	}
	centroid[0] = centroid[0] / hull.size();
	centroid[1] = centroid[1] / hull.size();
	HullVertexComparator h;
	h.m_centroid = centroid;
	std::sort(hull.begin(), hull.end(), h);
	
	for (int i=0; i<hull.size(); i++) {
		printf("Angle: %f\n", atan2(hull[i][1], hull[i][0]));
	}
	
	double max_x = 0.0;
	double max_y = 0.0;
	for (size_t i=0; i<hull.size(); i++) {
		if (hull[i][0] > max_x) max_x = hull[i][0];
		if (hull[i][1] > max_y) max_y = hull[i][1];
	}

	printf("Max: %f,%f\n", max_x, max_y);
	Mat img2((int)max_y+10, (int)max_x+10, CV_8UC1);
	//circle(img2, Point(centroid[0],centroid[1]), 10, 64000);
	hull.push_back(hull[0]); // For convinience...
	for (size_t i=0; i<hull.size()-1; i++) {
		line(img2, Point((int)hull[i][0],(int)hull[i][1]), Point((int)hull[i+1][0],(int)hull[i+1][1]), 255, 1);
		//line(img2, Point(0,0), Point((int)hull[i+1][0],(int)hull[i+1][1]), 64000, 1);
		printf("%f %f\n", hull[i][0], hull[i][1]);
	}

	outputPicture(img_name, img2);

	return img2;

	// First off, a few config variables
	int num_angle_buckets = 360;
	int num_distance_buckets = 640;
	
	// Where we're putting the (initial) result
	Mat img(num_angle_buckets, num_distance_buckets, CV_16UC1);
	for (int i=0; i<num_angle_buckets; i++) {
		for (int j=0; j<num_distance_buckets; j++) {
			img.at<unsigned short>(i,j) = 0;
		}
	}
	
	//printf("%f => %f\n", 9.0, sqrt(9.0));
	//exit(0);
	
	// Let's get cracking!
	double max_distance = 0.0;
	for (size_t i=0; i<hull.size(); i++) {
		double d = std::sqrt(hull[i][0]*hull[i][0] + hull[i][1]*hull[i][1]);
		if (d > max_distance) {
			max_distance = d;
		}
		//printf("%f %f\n", d, max_distance);
	}
	//printf("max_distance=%f\n", max_distance);
	unsigned int fullest_bucket = 0;
	unsigned int last_fullest_bucket = 0;
	for (size_t i=0; i<hull.size(); i++) {
		for (int j=0; j<num_angle_buckets; j++) {
			// Find the shortest vector to the line
			double angle = DEG2RAD(360.0 / (double)num_angle_buckets * (double)j);
			double nx = cos(angle);
			double ny = sin(angle);
			double dp = hull[i][0]*nx + hull[i][1]*ny;
			double x = hull[i][0] - dp*nx;
			double y = hull[i][1] - dp*ny;
			double dist = sqrt(x*x + y*y);
			int d_bucket = (int)(dist / max_distance * num_distance_buckets);
			img.at<unsigned short>(j, d_bucket) += 1;
			if (img.at<unsigned short>(j, d_bucket) > fullest_bucket) {
				fullest_bucket = img.at<unsigned short>(j, d_bucket);
			}
			if (last_fullest_bucket != fullest_bucket)
				printf("New max: %i\n", fullest_bucket);
			last_fullest_bucket = fullest_bucket;


#if 0
			double angle = DEG2RAD(360.0 / (double)num_angle_buckets * (double)j);
			double d_to_pnt = sqrt(hull[i][0]*hull[i][0] + hull[i][1]*hull[i][1]);
			double inside_angle = atan2(hull[i][1], hull[i][0]) - angle + DEG2RAD(90);
			double dist = fabs(d_to_pnt*cos(inside_angle));
			int d_bucket = (int)(dist / max_distance * num_distance_buckets);
			img.at<unsigned short>(j, d_bucket) += 1;
			if (img.at<unsigned short>(j, d_bucket) > fullest_bucket) {
				fullest_bucket = img.at<unsigned short>(j, d_bucket);
			}
			printf("%i,%i => %u\n", j, d_bucket, img.at<unsigned short>(j, d_bucket));
			if (last_fullest_bucket != fullest_bucket)
				printf("New max: %i\n", fullest_bucket);
			last_fullest_bucket = fullest_bucket;
			if (d_bucket < 0) {
				printf("angle/dist for pnt %i (d2=%f, <%f,%f>), angle %i: a=%f/d=%f (from cos(%f)) is %i (%i from %f)\n", i, d_to_pnt, hull[i][0], hull[i][1], j, RAD2DEG(angle), dist, RAD2DEG(inside_angle), img.at<unsigned short>(j, dist / max_distance * num_distance_buckets), (int)(dist / max_distance * num_distance_buckets), max_distance);
			}
#endif
		}
	}

	// For debugging purposes, let's rescale the image
	printf("Rescaling result with max=%u\n", fullest_bucket);
	for (int i=0; i<num_angle_buckets; i++) {
		for (int j=0; j<num_distance_buckets; j++) {
			if (img.at<unsigned short>(i, j) > 1) {
				//printf("Ratio: max=%i, %f/%f\n", max, ((double)img.at<unsigned short>(i, j)), 64500.0 / (double)max);
				img.at<unsigned short>(i, j) = (unsigned short)(((double)img.at<unsigned short>(i, j)) / (double)fullest_bucket * 255);
			}
		}
	}

	// Process it a little bit, to make it a little prettier
	GaussianBlur(img, img, Size(15, 15), 2);

	outputPicture(img_name, img);
}

vector<Rectangle3d> findRectanglesInImage(Mat src)
{
	
	rectangle_Finder_Debug_Counter++;
    Mat dst, color_dst;
	
	cvtColor(src, src, CV_BGR2GRAY);

    resize(src, src, Size(640,480));

	GaussianBlur(src, src, Size(5, 5), 0, 0);
	//outputPicture("blur", src);

	threshold( src, src, 127, 255, 0 );
	//outputPicture("threshold", src);
	
	// Let's mess with "goodfeaturestotrack"
	Mat src2;
	src.copyTo(src2);
	cvtColor(src2, src2, CV_GRAY2BGR);
    cv::FeatureDetector* blobDetect = new cv::GoodFeaturesToTrackDetector(1000,0.1,1.0,3,false,0.04);
	std::vector<cv::KeyPoint> kps;
	vector<Vec2f> kps2;
    blobDetect->detect(src, kps);
	for (int i=0; i<kps.size(); i++) {
        //cv::KeyPoint blob = kps[i];
        //rectangle(src2, cv::Point(blob.pt.x-2, blob.pt.y-2), cv::Point(blob.pt.x+blob.size+2, blob.pt.y+blob.size+2), cv::Scalar(0,255,255), 1, 8);
		Vec2f k;
		k[0] = kps[i].pt.x;
		k[1] = kps[i].pt.y;
		kps2.push_back(k);
    }
	
	// Canny it!
	Mat canny_output(src.size(), src.type());
	Canny(src, canny_output, 1, 2);
	outputPicture("canny", canny_output);
	canny_output.convertTo(canny_output, CV_16UC1);
	
	// Take the canny, make up a bunch of points on it, and then call it a point cloud.
	vector<Vec2f> canny_based_point_cloud;
	for (int i=0; i<canny_output.size().width; i++) {
		for (int j=0; j<canny_output.size().height; j++) {
			//printf("%i,%i => %u\n", i, j, canny_output.at<unsigned short>(j,i));
			if (canny_output.at<unsigned short>(j,i) > 1) {
				canny_based_point_cloud.push_back(Vec2f(i,j));
			}
		}
	}
	kps2 = canny_based_point_cloud;
	printf("I found %u points from the canny.\n", kps2.size());

	for (int i=0; i<kps2.size(); i++) {
		rectangle(src2, cv::Point(kps2[i][0]-2, kps2[i][1]-2), cv::Point(kps2[i][0]+2, kps2[i][1]+2), cv::Scalar(0,255,255), 1, 8);
	}
	outputPicture("gftt", src2);
	
	// Test of cornerHarris
	Mat harris_output(src2.size(), CV_32FC1);
	Mat harris2;
	cvtColor(src2, harris2, CV_BGR2GRAY);
	preCornerDetect(harris2, harris_output, 3);
	outputPicture("harris", harris_output);

	Mat dilated_corners;
	dilate(harris_output, dilated_corners, Mat(), Point(0,0));
	Mat corner_mask = harris_output == dilated_corners;
	outputPicture("corner", corner_mask);
	
	// Get rid of all blobs that are far away
	vector<vector<Vec2f> > point_clouds;
	for (int i=0; i<kps2.size(); i++) {
		bool did_add = false;
		// We're deciding which point cloud to add kps2[i] to
		for (int j=0; j<point_clouds.size(); j++) {
			// Find the nearest one in the point cloud
			double min_dist = 10000.0;
			for (int k=0; k<point_clouds[j].size(); k++) {
				Vec2f p1 = point_clouds[j][k];
				Vec2f p2 = kps2[i];
				Vec2f diff;
				diff[0] = p2[0]-p1[0];
				diff[1] = p2[1]-p1[1];
				double dist = sqrt(diff[0]*diff[0] + diff[1]*diff[1]);
				if (dist < min_dist) min_dist = dist;
			}
			if (min_dist < 40) {
				// Add it
				point_clouds[j].push_back(kps2[i]);
				did_add = true;
				break;
			}
		}

		if (!did_add) {
			// Create a new point cloud, containing this point
			vector<Vec2f> cloud;
			cloud.push_back(kps2[i]);
			point_clouds.push_back(cloud);
		}
	}

	// Merge nearby point clouds
	printf("Number of point clouds: %i\n", point_clouds.size());
	vector<vector<Vec2f> > new_point_clouds;
	for (int i=0; i<point_clouds.size(); i++) {
		vector<Vec2f> cloud;
		bool did_merge_down = false;
		for (int k=0; k<point_clouds[i].size(); k++) {
			cloud.push_back(point_clouds[i][k]);
		}
		for (int j=0; j<point_clouds.size(); j++) {
			//printf("Hello? %i %i %i\n", i, j, point_clouds.size());
			if (i == j) {
				continue;
			}
			double min_dist = 10000.0;
			for (int k=0; k<point_clouds[i].size(); k++) {
				for (int m=0; m<point_clouds[j].size(); m++) {
					Vec2f p1 = point_clouds[j][m];
					Vec2f p2 = point_clouds[i][k];
					Vec2f diff;
					diff[0] = p2[0]-p1[0];
					diff[1] = p2[1]-p1[1];
					double dist = sqrt(diff[0]*diff[0] + diff[1]*diff[1]);
					if (dist < min_dist) {
						min_dist = dist;
					}
				}
			}

			if (min_dist < 20) {
				// Tack it onto cloud
				for (int k=0; k<point_clouds[j].size(); k++) {
					cloud.push_back(point_clouds[j][k]);
				}
				
				if (j < i) {
					did_merge_down = true;
				}
			}
			
#if 0
			if (min_dist < 40) {
				// Merge them.
				if (j > i) {
				printf("Merging clouds %i and %i dist=%f.\n", i, j, min_dist);
				vector<Vec2f> cloud;
				for (int k=0; k<point_clouds[i].size(); k++) {
					cloud.push_back(point_clouds[i][k]);
				}
				for (int k=0; k<point_clouds[j].size(); k++) {
					cloud.push_back(point_clouds[j][k]);
				}
				new_point_clouds.push_back(cloud);
				}
			} else {
				// Don't merge them (i.e. do nothing).
				vector<Vec2f> cloud;
				for (int k=0; k<point_clouds[i].size(); k++) {
					cloud.push_back(point_clouds[i][k]);
				}
				new_point_clouds.push_back(cloud);
			}
#endif
		}
		if (!did_merge_down) {
			new_point_clouds.push_back(cloud);
		}
	}

	// Convex hull this sucker
	point_clouds = new_point_clouds;
	vector<vector<Vec2f> > hulls;
	printf("There are %i (%i) point clouds.\n", point_clouds.size(), new_point_clouds.size());
	for (int i=0; i<point_clouds.size(); i++) {
		printf("Point cloud %i has %i points.\n", i, point_clouds[i].size());
		vector<Vec2f> hull;
		if (point_clouds[i].size() > 0) {
			convexHull(point_clouds[i], hull);
			for (int i=0; i<hull.size()-1; i++) {
				line(src2, Point(hull[i]), Point(hull[i+1]), Scalar(0,0,255), 1, 8);
			}
			hulls.push_back(hull);
		}
	}

	// Find the four corners for each hull
	for (int i=0; i<hulls.size(); i++) {
		Vec2f centroid(0,0);
		for (int j=0; j<hulls[i].size(); j++) {
			centroid[0] += hulls[i][j][0];
			centroid[1] += hulls[i][j][1];
		}
		centroid[0] = centroid[0] / hulls[i].size();
		centroid[1] = centroid[1] / hulls[i].size();

		// Now, order the points so that we can apply the other formula
		for (int j=0; j<hulls[i].size(); j++) {
		}
		HullVertexComparator h;
		h.m_centroid = centroid;
		std::sort(hulls[i].begin(), hulls[i].end(), h);
		for (int j=0; j<hulls[i].size(); j++) {
			double a = getReferenceAngle(hulls[i][j], centroid);
			printf("Centroid: <%f,%f> Hull %i, pnt %i: <%f,%f> => %f\n", centroid[0], centroid[1], i, j, hulls[i][j][0], hulls[i][j][1], a);
		}
		circle(src2, Point(centroid[0], centroid[1]), 5, Scalar(0,255,255), 2);

		// Find the centroid
		// But first, find the area.
#if 1
		double A = 0.0;
#if 0
		for (int j=0; j<hulls[i].size()-1; j++) {
		}
#endif

		// Now, find the centroid
		//Vec2f centroid(0,0);
		centroid[0] = centroid[1] = 0;
		for (int j=0; j<hulls[i].size()-1; j++) {
			double a = hulls[i][j][0] * hulls[i][j+1][1] - hulls[i][j+1][0] * hulls[i][j][1];
			centroid[0] += (hulls[i][j][0] + hulls[i][j+1][0]) * a;
			centroid[1] += (hulls[i][j][1] + hulls[i][j+1][1]) * a;
			A += a;
			printf("centroid=<%f,%f> a=%f A=%f\n", centroid[0], centroid[1], a, A);
		}

		// Do the last vertex
		double a = hulls[i][hulls[i].size()-1][0] * hulls[i][0][1] - hulls[i][0][0] * hulls[i][hulls[i].size()-1][1];
		centroid[0] += (hulls[i][hulls[i].size()-1][0] + hulls[i][0][0]) * a;
		centroid[1] += (hulls[i][hulls[i].size()-1][1] + hulls[i][0][1]) * a;
		A += a;

		A = 0.5 * A;
		centroid[0] *= 1.0/(6.0 * A);
		centroid[1] *= 1.0/(6.0 * A);
		printf("We have a centroid, and I'm not telling what it is. (%f,%f A=%f)\n", centroid[0], centroid[1], A);
#endif

		circle(src2, Point(centroid[0], centroid[1]), 10, Scalar(0,255,255), 2);
	}

	// Find the rectangle edges for each hull...
	vector<Rectangle3d> new_rects_3d;
	for (int i=0; i<hulls.size(); i++) {
		if (hulls[i].size() < 4) {
			continue;
		}

		// Loop through and find the four angles closest to 90deg
		// For convinience, add #0 to the end

		double closest_dps[4] = {90.0, 90.0, 90.0, 90.0};
		int closest_indices[4] = {-1,-1,-1,-1};
		double SKIP_DIST = 20.0;
		int SKIP_CNT = 3;
		for (int j=0; j<hulls[i].size(); j++) {
			int index0 = (j-SKIP_CNT+hulls[i].size())%hulls[i].size();//(j-SKIP_CNT + hulls[i].size()) % hulls[i].size();
			int index1 = j;
			int index2 = (j+SKIP_CNT)%hulls[i].size();//(j+SKIP_CNT) % hulls[i].size();
			Vec2f v1 = normalize(hulls[i][index1] - hulls[i][index2]);
			Vec2f v2 = normalize(hulls[i][index0] - hulls[i][index1]);
			double dp = fabs(v1[0]*v2[0] + v1[1]*v2[1]);
			for (int k=0; k<4; k++) {
				if (dp < closest_dps[k]) {
					// Check to make sure that we're not too close to another one
					bool is_too_close = false;
					for (int m=0; m<4; m++) {
						if (m != k) {
							// Check the distance
							if (closest_indices[m] > -1) {
								double dist = vec_dist(hulls[i][index1],hulls[i][closest_indices[m]]);
								printf("Distance between %i and %i is %f\n", index1, closest_indices[m], dist);
								if (dist < SKIP_DIST) {
									is_too_close = true;
									break;
								}
							}
						}
					}
					if (!is_too_close) {
						closest_dps[k] = dp;
						closest_indices[k] = j;
						break;
					}
				}
			}
		}

		// Render these four points
		for (int j=0; j<4; j++) {
			printf("Closest DPS %i: %f\n", j, closest_dps[j]);
			circle(src2, Point(hulls[i][closest_indices[j]][0], hulls[i][closest_indices[j]][1]), 10, Scalar(255,0,255), 2);
		}
		
		// Pass them into the Rectangle3D framework
		Polygon poly;
		vector<Vec2i> corners;
		for (int j=0; j<4; j++) {
			corners.push_back(Vec2i(hulls[i][closest_indices[j]][0],hulls[i][closest_indices[j]][1]));
		}
		poly.fromCorners(corners);
		Rectangle3d r(poly);
		r.solve((double)src.size().width, (double)src.size().height, 60.0,60.0, 7.0,2.75);

		Vec4i b = r.get_image_bounds();
		printf("%i %i %i %i\n", b[1],b[0],b[3],b[2]);
		rectangle(src2, Point(b[0], b[2]), Point(b[1], b[3]), Scalar(0,255,0), 1, 8);
		char render_text[2048];
		snprintf(render_text, 2048, "%i %i %i %i", b[0], b[1], b[2], b[3]);
		putText(src2, render_text, Point(b[0], b[2]), FONT_HERSHEY_SIMPLEX, 1.0, Scalar(255,255,255));
		new_rects_3d.push_back(r);
		
		// Output some information about this point cloud
		char bufname[2048];
		snprintf(bufname, 2048, "pc_%02i", i);
		Mat pc_img;
		for (int j=0; j<hulls[i].size(); j++) {
			circle(src2, Point(hulls[i][j][0], hulls[i][j][1]), 5, Scalar(255,0,255), 2);
		}
		outputPicture(bufname, pc_img);
	}
	outputPicture("newalg", src2);
	return new_rects_3d;
	//printf("EXITING... Press ctrl-C\n"); while(1);

	// Let's do that again, but using a haugh transform
	vector<Mat> contour_images;
	for (int i=0; i<hulls.size(); i++) {
		char img_name[2048];
		snprintf(img_name, 2048, "haugh_%02i", i);
		contour_images.push_back(myHaughTransform(hulls[i], img_name));
	}

	outputPicture("cvxhull", src2);

	// Alright, back to the "canonical" algorithm
    Canny( src, dst, 50, 200, 3, true );
    cvtColor( dst, color_dst, CV_GRAY2BGR );

	//outputPicture("canny", color_dst);

    vector<Vec4i> lines;
	// Rho (distance), theta, threshold, minlen, maxgap
    //HoughLinesP( dst, lines, 1, CV_PI/3600, 10, 10, 15 );
    HoughLinesP( contour_images[0], lines, 1, CV_PI/3600, 10, 10, 15 );
	Scalar color1(0,0,255);
	Scalar color2(255,0,0);
	Scalar colors[2] = {color1, color2};
    for( size_t i = 0; i < lines.size(); i++ )
    {
        line( color_dst, Point(lines[i][0], lines[i][1]),
            Point(lines[i][2], lines[i][3]), colors[i%2], 3, 8 );
    }

	//outputPicture("lines", color_dst);

    // Detect rectangles
    printf("We found %lu lines.\n", lines.size());
    vector<Polygon> rectangles_2d = detectRectangles(lines);
    vector<Rectangle3d> rectangles_3d;
    for (size_t i=0; i<rectangles_2d.size(); i++) {
      Rectangle3d r(rectangles_2d[i]);
      r.solve((double)src.size().width, (double)src.size().height, 60.0,60.0, 7.0,2.75);
      rectangles_3d.push_back(r);
	  
	  // Add the edges of the rects to color_dst
	  Vec4i b = r.get_image_bounds();
	  printf("%i %i %i %i\n", b[1],b[0],b[2],b[3]);
	  rectangle(color_dst, Point(b[1], b[0]), Point(b[2], b[3]), Scalar(0,255,0), 1, 8);
    }
	
	outputPicture("lines", color_dst);

    // Print information on all of them with roughly the right aspect ratio
    for (size_t i=0; i<rectangles_3d.size(); i++) {
      double ar = rectangles_3d[i].aspect_ratio();
      if (ar > 2.0 && ar < 3.0) {
		printf("Distance to rect %lu is %f, ar=%f\n", i, rectangles_3d[i].centroid_dist(), ar);
      }
    }
	
	printf("Exiting at the end of the rectangle processor.\n");
	exit(0);
	
	return rectangles_3d;
}

#if 0
int main(int argc, char** argv)
{
    std::cout << "Hello?\n";
 //printf("Welcome to the thing...\n");
    Mat src, dst, color_dst;
    if( argc != 2 || !(src=imread(argv[1], 0)).data)
        return -1;

    resize(src, src, Size(640,480));
    Canny( src, dst, 50, 200, 3 );
    cvtColor( dst, color_dst, CV_GRAY2BGR );

   //printf("Hello?\n");

#if 0
    vector<Vec2f> lines;
    HoughLines( dst, lines, 1, CV_PI/180, 100 );

    for( size_t i = 0; i < lines.size(); i++ )
    {
        float rho = lines[i][0];
        float theta = lines[i][1];
        double a = cos(theta), b = sin(theta);
        double x0 = a*rho, y0 = b*rho;
        Point pt1(cvRound(x0 + 1000*(-b)),
                  cvRound(y0 + 1000*(a)));
        Point pt2(cvRound(x0 - 1000*(-b)),
                  cvRound(y0 - 1000*(a)));
        line( color_dst, pt1, pt2, Scalar(0,0,255), 3, 8 );
    }
#else
    vector<Vec4i> lines;
    HoughLinesP( dst, lines, 1, CV_PI/180, 50, 10, 20 );
    for( size_t i = 0; i < lines.size(); i++ )
    {
        line( color_dst, Point(lines[i][0], lines[i][1]),
            Point(lines[i][2], lines[i][3]), Scalar(0,0,255), 3, 8 );
    }

    // Detect rectangles
   //printf("We found %lu lines.\n", lines.size());
    vector<Polygon> rectangles_2d = detectRectangles(lines);
    vector<Rectangle3d> rectangles_3d;
    for (size_t i=0; i<rectangles_2d.size(); i++) {
      Rectangle3d r(rectangles_2d[i]);
      r.solve((double)src.size().width, (double)src.size().height, 60.0,60.0, 7.0,2.75);
      rectangles_3d.push_back(r);
    }

    // Print information on all of them with roughly the right aspect ratio
    for (size_t i=0; i<rectangles_3d.size(); i++) {
      double ar = rectangles_3d[i].aspect_ratio();
      if (ar > 2.0 && ar < 3.0) {
	printf("Distance to rect %lu is %f, ar=%f\n", i, rectangles_3d[i].centroid_dist(), ar);
      }
    }
#endif
#if 1
    namedWindow( "Source", 1 );
    imshow( "Source", src );

    namedWindow( "Detected Lines", 1 );
    imshow( "Detected Lines", color_dst );

    waitKey(0);
#endif
    return 0;
}
#endif
