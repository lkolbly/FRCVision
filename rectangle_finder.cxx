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

#define POLYGON_DEBUG 0
#define ENABLE_DEBUG_IMAGES 0

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

Vec2f my_normalize(Vec2f v)
{
	//printf("Normalizing? Maybe?\n");
	double len = sqrt(v[0]*v[0] + v[1]*v[1]);
	return Vec2f(v[0]/len, v[1]/len);
}

vector<Vec4i> Polygon::get_edges(void) {
	return m_edges;
}

void Polygon::fromCorners(std::vector<cv::Vec2i> corners)
{
	for (size_t i=0; i<corners.size(); i++) {
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
  int m = min_index(d1,d2,d3,d4);
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
	for (size_t i=0; i<m_edges.size(); i++) {
#if POLYGON_DEBUG
		printf("Edge: %i %i %i %i\n", m_edges[i][0], m_edges[i][1], m_edges[i][2], m_edges[i][3]);
#endif
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
				if (p.should_add_line(lines[j])) {
					p.add_line(lines[j]);
					j = 0;
				}
			}
		}
		if (p.m_edges.size() == 4) {
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
		}

		// If we're the first one like this one that we've seen...
		if (min_diff > 2.0) {
			unique_rectangles.push_back(rectangles[i]);
		}
	}
	return unique_rectangles;
}

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
  //v[0] = tan(DEG2RAD(x/w*fovx));
  //v[1] = tan(DEG2RAD(y/h*fovy));
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
  }
  sum = sum / 4.0;

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
	//Vec2f coefs[4];
	for (int i=0; i<4; i++) {
		m_coefs[i] = coef_from_point((double)m_corners[i][0], (double)m_corners[i][1], w, h, fovx, fovy);
		//printf("%i: %f,%f => coefs=%f,%f\n", i, m_corners[i][0], m_corners[i][1], m_coefs[i][0], m_coefs[i][1]);
	}

	// Solve to figure out the distance along said lines to make the angles 90deg
	// Yes, it's an iterative solver. Sorry, Ms. Harrelson.
	double dist[4];
	dist[0] = 1.0; // Fix one of the distances for now.
	dist[1] = dist[2] = dist[3] = 1.0;
	double last_inaccuracy = 0.0;
	while (1) {
		double sum = find_squareness(m_coefs, dist);
		//printf("Inaccuracy: %f (k=%f)\n", sum, k);
		/*if (last_inaccuracy-sum < 0.001) {
			break;
		}*/
		if (k < 0.000001) {
			break;
		}

		// Figet with the distances, and go with the optimal one.
		int did_change = 0;
		for (int i=0; i<4; i++) {
			double orig = dist[i];
			dist[i] = orig + k;
			double incr = find_squareness(m_coefs, dist);
			dist[i] = orig - k;
			double decr = find_squareness(m_coefs, dist);
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
	}

	// Copy out the distances.
	for (int i=0; i<4; i++) {
		m_dist[i] = dist[i];
	}

	// Now we can scale everyone to match the width to the known width
	// We'll assume the 'width' is the 'longest' edge.
	// Find the longest edge
	scale_distances:;
	double est_w = 0.0;
	vector<Vec3f> pnts = get_points();
	for (int i=0; i<4; i++) {
		double d = dist3d(pnts[i], pnts[(i+1)%4]);
		//printf("Known=%f, est=%f, d=%f at pnt=%i\n", known_w, est_w, d, i);
		//printf("%f=>%f,%f,%f %f=>%f,%f,%f\n", m_dist[i], pnts[i][0],pnts[i][1],pnts[i][2], m_dist[(i+1)%4], pnts[(i+1)%4][0],pnts[(i+1)%4][1],pnts[(i+1)%4][2]);
		if (d > est_w) {
			est_w = d;
		}
	}

	// Scale all of the distances by the needed correction factor
	double factor = known_w / est_w;
	for (int i=0; i<4; i++) {
		m_dist[i] *= factor;
	}
	if (est_w/known_w > 1.5 || known_w/est_w > 1.5) {
		goto scale_distances;
	}

	// Some debugging info
	//printf("Adjusted distances to %f,%f,%f,%f\n", m_dist[0], m_dist[1], m_dist[2], m_dist[3]);
#if 0
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
#endif
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

Vec3f calc_normal(Vec3f p1, Vec3f p2, Vec3f p3)
{
	Vec3f a(p1[0]-p2[0], p1[1]-p2[1], p1[2]-p2[2]);// = p1 - p2;
	Vec3f b(p3[0]-p2[0], p3[1]-p2[1], p3[2]-p2[2]);// = p3 - p2;
	Vec3f n;
	n[0] = a[1]*b[2] + a[2]*b[1];
	n[1] = a[2]*b[0] + a[0]*b[2];
	n[2] = a[1]*b[0] + a[0]*b[1];
	return n;
}

// Compute the (average) three-space normal for the quadrangle
Vec3f Rectangle3d::normal(void)
{
	vector<Vec3f> pnts = get_points();
	Vec3f normals[4];
	for (int i=0; i<4; i++) {
	normals[i] = calc_normal(pnts[(i+0)%4], pnts[(i+1)%4], pnts[(i+2)%4]);
	}
	Vec3f n(0,0,0);
	for (int i=0; i<4; i++) {
		n[0] += normals[i][0];
		n[1] += normals[i][1];
		n[2] += normals[i][2];
	}
	n[0] /= 4.0;
	n[1] /= 4.0;
	n[2] /= 4.0;
	//printf("Computing a normal...\n");
	return n;
}

double Rectangle3d::offcenter_angle(void)
{
	Vec3f centroid = get_centroid();
	Vec2f c(centroid[0], centroid[1]);
	Vec2f tr = -c;
	Vec3f n3 = normal();
	Vec2f n(n3[0],n3[1]);
	//printf("Hello?\n");
	n = my_normalize(n);
	tr = my_normalize(tr);
	double dp = tr[0]*n[0] + tr[1]*n[1];
	//printf("dp=%f\n", dp);
	return 90.0 - RAD2DEG(acos(dp)) / 2.0;
}

Vec4i Rectangle3d::get_image_bounds(void)
{
	return m_p.get_bounds();
}

int rectangle_Finder_Debug_Counter = 0;

void outputPicture(const char *name, Mat img)
{
#if ENABLE_DEBUG_IMAGES
	char buf[1024];
	snprintf(buf, 1024, "debug/%s_%04i.jpg", name, rectangle_Finder_Debug_Counter);
	imwrite(buf, img);
#endif
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
	for (size_t j=0; j<hull.size(); j++) {
		centroid[0] += hull[j][0];
		centroid[1] += hull[j][1];
	}
	centroid[0] = centroid[0] / hull.size();
	centroid[1] = centroid[1] / hull.size();
	HullVertexComparator h;
	h.m_centroid = centroid;
	std::sort(hull.begin(), hull.end(), h);
	
	for (size_t i=0; i<hull.size(); i++) {
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
	vector<Vec2f> kps2;
#if 0
    cv::FeatureDetector* blobDetect = new cv::GoodFeaturesToTrackDetector(1000,0.1,1.0,3,false,0.04);
	std::vector<cv::KeyPoint> kps;
    blobDetect->detect(src, kps);
	for (size_t i=0; i<kps.size(); i++) {
        //cv::KeyPoint blob = kps[i];
        //rectangle(src2, cv::Point(blob.pt.x-2, blob.pt.y-2), cv::Point(blob.pt.x+blob.size+2, blob.pt.y+blob.size+2), cv::Scalar(0,255,255), 1, 8);
		Vec2f k;
		k[0] = kps[i].pt.x;
		k[1] = kps[i].pt.y;
		kps2.push_back(k);
    }
#endif
	
	// Canny it!
	Mat canny_output(src.size(), src.type());
	Canny(src, canny_output, 1, 2);
#if ENABLE_DEBUG_IMAGES
	outputPicture("canny", canny_output);
#endif
	canny_output.convertTo(canny_output, CV_16UC1);
	
	// Take the canny, make up a bunch of points on it, and then call it a point cloud.
	vector<Vec2f> canny_based_point_cloud;
	for (int i=0; i<canny_output.size().width; i++) {
		for (int j=0; j<canny_output.size().height; j++) {
			if (canny_output.at<unsigned short>(j,i) > 1) {
				canny_based_point_cloud.push_back(Vec2f(i,j));
			}
		}
	}
	kps2 = canny_based_point_cloud;

#if ENABLE_DEBUG_IMAGES
	printf("I found %u points from the canny.\n", kps2.size());
	for (size_t i=0; i<kps2.size(); i++) {
		rectangle(src2, cv::Point(kps2[i][0]-2, kps2[i][1]-2), cv::Point(kps2[i][0]+2, kps2[i][1]+2), cv::Scalar(0,255,255), 1, 8);
	}
	outputPicture("gftt", src2);
#endif

#if 0
	// Test of cornerHarris
	Mat harris_output(src2.size(), CV_32FC1);
	Mat harris2;
	cvtColor(src2, harris2, CV_BGR2GRAY);
	preCornerDetect(harris2, harris_output, 3);
	outputPicture("harris", harris_output);
#endif
	
#if 0
	Mat dilated_corners;
	dilate(harris_output, dilated_corners, Mat(), Point(0,0));
	Mat corner_mask = harris_output == dilated_corners;
	outputPicture("corner", corner_mask);
#endif
	
	// Get rid of all blobs that are far away
	int POINT_CLOUD_MIN_DIST = 5; // Distance which makes two point clouds the same
	vector<vector<Vec2f> > point_clouds;
	for (size_t i=0; i<kps2.size(); i++) {
		bool did_add = false;
		// We're deciding which point cloud to add kps2[i] to
		for (size_t j=0; j<point_clouds.size(); j++) {
			// Find the nearest one in the point cloud
			double min_dist = 10000.0;
			for (size_t k=0; k<point_clouds[j].size(); k++) {
				Vec2f p1 = point_clouds[j][k];
				Vec2f p2 = kps2[i];
				Vec2f diff;
				diff[0] = p2[0]-p1[0];
				diff[1] = p2[1]-p1[1];
				double dist = sqrt(diff[0]*diff[0] + diff[1]*diff[1]);
				if (dist < min_dist) min_dist = dist;
				if (min_dist < POINT_CLOUD_MIN_DIST) break;
			}
			if (min_dist < POINT_CLOUD_MIN_DIST) {
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
	//printf("Number of point clouds: %i\n", point_clouds.size());
	vector<vector<Vec2f> > new_point_clouds;
	for (size_t i=0; i<point_clouds.size(); i++) {
		vector<Vec2f> cloud;
		bool did_merge_down = false;
		for (size_t k=0; k<point_clouds[i].size(); k++) {
			cloud.push_back(point_clouds[i][k]);
		}
		for (size_t j=0; j<point_clouds.size(); j++) {
			//printf("Hello? %i %i %i\n", i, j, point_clouds.size());
			if (i == j) {
				continue;
			}
			double min_dist = 10000.0;
			for (size_t k=0; k<point_clouds[i].size(); k++) {
				for (size_t m=0; m<point_clouds[j].size(); m++) {
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

			if (min_dist < POINT_CLOUD_MIN_DIST) {
				// Tack it onto cloud
				for (size_t k=0; k<point_clouds[j].size(); k++) {
					cloud.push_back(point_clouds[j][k]);
				}
				
				if (j < i) {
					did_merge_down = true;
				}
			}
		}
		if (!did_merge_down) {
			new_point_clouds.push_back(cloud);
		}
	}

	// Convex hull this sucker
	point_clouds = new_point_clouds;
	vector<vector<Vec2f> > hulls;
	//printf("There are %i (%i) point clouds.\n", point_clouds.size(), new_point_clouds.size());
	for (size_t i=0; i<point_clouds.size(); i++) {
		//printf("Point cloud %i has %i points.\n", i, point_clouds[i].size());
		vector<Vec2f> hull;
		if (point_clouds[i].size() > 0) {
			convexHull(point_clouds[i], hull);
			for (size_t i=0; i<hull.size()-1; i++) {
				line(src2, Point(hull[i]), Point(hull[i+1]), Scalar(0,0,255), 1, 8);
			}
			hulls.push_back(hull);
		}
	}

#if 0
	// Find the four corners for each hull
	for (size_t i=0; i<hulls.size(); i++) {
		Vec2f centroid(0,0);
		for (size_t j=0; j<hulls[i].size(); j++) {
			centroid[0] += hulls[i][j][0];
			centroid[1] += hulls[i][j][1];
		}
		centroid[0] = centroid[0] / hulls[i].size();
		centroid[1] = centroid[1] / hulls[i].size();

		// Now, order the points so that we can apply the other formula
		HullVertexComparator h;
		h.m_centroid = centroid;
		std::sort(hulls[i].begin(), hulls[i].end(), h);
		for (size_t j=0; j<hulls[i].size(); j++) {
			double a = getReferenceAngle(hulls[i][j], centroid);
			//printf("Centroid: <%f,%f> Hull %i, pnt %i: <%f,%f> => %f\n", centroid[0], centroid[1], i, j, hulls[i][j][0], hulls[i][j][1], a);
		}
		circle(src2, Point(centroid[0], centroid[1]), 5, Scalar(0,255,255), 2);

		// Find the centroid
		// But first, find the area.
		double A = 0.0;

		// Now, find the centroid
		//Vec2f centroid(0,0);
		centroid[0] = centroid[1] = 0;
		for (size_t j=0; j<hulls[i].size()-1; j++) {
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

		circle(src2, Point(centroid[0], centroid[1]), 10, Scalar(0,255,255), 2);
	}
#endif

	// Find the rectangle edges for each hull...
	vector<Rectangle3d> new_rects_3d;
	for (size_t i=0; i<hulls.size(); i++) {
		if (hulls[i].size() < 4) {
			continue;
		}

		// Loop through and find the four angles closest to 90deg
		// For convinience, add #0 to the end

		double closest_dps[4] = {90.0, 90.0, 90.0, 90.0};
		int closest_indices[4] = {-1,-1,-1,-1};
		double SKIP_DIST = 20.0;
		int SKIP_CNT = 1;
		for (size_t j=0; j<hulls[i].size(); j++) {
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
								//printf("Distance between %i and %i is %f\n", index1, closest_indices[m], dist);
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
			//printf("Closest DPS %i: %f\n", j, closest_dps[j]);
			circle(src2, Point(hulls[i][closest_indices[j]][0], hulls[i][closest_indices[j]][1]), 10, Scalar(255,0,255), 2);

			char render_text[2048];
			snprintf(render_text, 2048, "%.2f", closest_dps[j]);
			putText(src2, render_text, Point(hulls[i][closest_indices[j]][0], hulls[i][closest_indices[j]][1]), FONT_HERSHEY_SIMPLEX, 0.5, Scalar(255,255,0));
		}
		
		// Check to make sure that we actually have four points
		bool is_valid_rect = true;
		for (int j=0; j<4; j++) {
			if (closest_dps[j] > 3.0) {
				// Oh no!
				is_valid_rect = false;
				break;
			}
		}
		if (!is_valid_rect) continue; // Onto the next point cloud!
		
		// Pass them into the Rectangle3D framework
		Polygon poly;
		vector<Vec2i> corners;
		for (size_t j=0; j<4; j++) {
			corners.push_back(Vec2i(hulls[i][closest_indices[j]][0],hulls[i][closest_indices[j]][1]));
		}
		poly.fromCorners(corners);
		Rectangle3d r(poly);
		//r.solve((double)src.size().width, (double)src.size().height, 75.75*0.47,48.72*0.47, 53,11.75);
		r.solve((double)src.size().width, (double)src.size().height, 45.75,28.72, 53,11.75);

		Vec4i b = r.get_image_bounds();

		// Check to make sure that "0" isn't in the bounds.
		// It's naive, but it'll also remove many false-positives from the test cases.
		for (size_t j=0; j<4; j++) {
			if (b[j] <= 0) {
				is_valid_rect = false;
			}
		}
		if (!is_valid_rect) continue;

		//printf("%i %i %i %i\n", b[1],b[0],b[3],b[2]);
		rectangle(src2, Point(b[0], b[2]), Point(b[1], b[3]), Scalar(0,255,0), 1, 8);
		char render_text[2048];
		snprintf(render_text, 2048, "%i %i %i %i", b[0], b[1], b[2], b[3]);
		putText(src2, render_text, Point(b[1], b[3]), FONT_HERSHEY_SIMPLEX, 1.0, Scalar(255,255,255));
		new_rects_3d.push_back(r);
		
		// Output some information about this point cloud
		char bufname[2048];
		snprintf(bufname, 2048, "pc_%02i", i);
		Mat pc_img;
		for (size_t j=0; j<hulls[i].size(); j++) {
			circle(src2, Point(hulls[i][j][0], hulls[i][j][1]), 5, Scalar(255,0,255), 2);
		}
		outputPicture(bufname, pc_img);
	}
	outputPicture("newalg", src2);
	return new_rects_3d;
}
