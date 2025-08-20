// #include"detector.h"
#include<iostream>
#include<vector>
#include<opencv2/core/core.hpp>
#include<opencv2/imgproc/imgproc.hpp>
#include<opencv2/highgui/highgui.hpp>
#include<assert.h>
using namespace std;
using namespace cv;

float GetMidValue(cv::Mat input);
cv::Mat insertmat_dim2to3(cv::Mat m3,cv::Mat m2,int d);
void hardvalsmall(cv::Mat w,double t);
cv::Mat hardval2(cv::Mat w,double t);
cv::Mat ToOne(cv::Mat img);
void removeSmall(cv::Mat img,double threshold);
cv::Mat avgImg(cv::Mat img,int k);
cv::Mat Img_in(cv::Mat img);
cv::Mat Img_in2(cv::Mat img);
void drawPoint(cv::Mat img,cv::Mat res);
void sortMat(Mat &stats, int colId);
void drawimg(Mat img,string name,int time_wait);
