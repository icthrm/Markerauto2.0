#include<opencv2/core/core.hpp>
#include<opencv2/imgproc/imgproc.hpp>
#include<opencv2/highgui/highgui.hpp>
#include<opencv2/imgproc/imgproc.hpp>
#include<opencv2/highgui/highgui.hpp>
#include<math.h>
using namespace cv;

#define _USE_MATH_DEFINES

cv::Mat normalize1(cv::Mat img);
double roundness(cv::Mat label_img);
void find_local_peak(cv::Mat img, int m,int n,int m_w,int n_w,int &out_m,int &out_n);
bool have255(cv::Mat img);
float sum_pix_float(cv::Mat img);
