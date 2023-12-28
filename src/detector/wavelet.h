#include<iostream>
#include<assert.h>
#include <time.h>
#include<opencv2/core/core.hpp>
#include<opencv2/imgproc/imgproc.hpp>
#include<opencv2/highgui/highgui.hpp>
#include"fundmental.h"
#include"BaseFun.h"

using namespace std;
using namespace cv;

cv::Mat fullZero2(cv::Mat k);
cv::Mat atrous_cv(cv::Mat ai,cv::Mat k);
cv::Mat wavelet(cv::Mat a,cv::Mat A);
cv::Mat waveletprocess2(cv::Mat img,int J);
