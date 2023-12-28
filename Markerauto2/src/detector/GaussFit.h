#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
using namespace std;
using namespace cv;

cv::Mat compute_center_Gauss(cv::Mat pic);
float compute_gauss_error(cv::Mat pic,cv::Mat dir);
