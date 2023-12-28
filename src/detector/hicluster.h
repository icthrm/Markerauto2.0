#include<iostream>
#include<math.h>
#include<opencv2/core/core.hpp>
#include<opencv2/imgproc/imgproc.hpp>
#include<opencv2/highgui/highgui.hpp>
// #include <mlpack/core.hpp>
// #include <mlpack/methods/mean_shift.hpp>
using namespace std;
using namespace cv;
// using namespace mlpack;

cv::Mat compute_distance(cv::Mat loc);
cv::Mat hicluster(cv::Mat img,double r0);
