#include<opencv2/core/core.hpp>
#include<opencv2/imgproc/imgproc.hpp>
#include<opencv2/highgui/highgui.hpp>
#include<iostream>
using namespace std;
using namespace cv;

cv::Mat expandImage(cv::Mat img,int row,int col);
double corr_p_p(cv::Mat p1,cv::Mat p2);
