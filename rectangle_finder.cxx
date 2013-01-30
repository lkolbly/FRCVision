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
  if (d < 5) {
    return 1;
  }

  return 0;
}

Vec4i Polygon::get_bounds(void)
{
	int left=1000, right=0, top=0, bottom=1000;
	for (int i=0; i<m_edges.size(); i++) {
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
	printf("%lu %i %f\n", i, j, diff);
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
    //printf("Inaccuracy: %f (k=%f)\n", sum, k);
    if (last_inaccuracy-sum < 0.001) {
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

vector<Rectangle3d> findRectanglesInImage(Mat src)
{
    Mat dst, color_dst;

    resize(src, src, Size(640,480));
    Canny( src, dst, 50, 200, 3 );
    cvtColor( dst, color_dst, CV_GRAY2BGR );
	
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
