/*
  This file contains definitions for functions to compute transforms from
  image feature correspondences

  Copyright (C) 2006-2010  Rob Hess <hess@eecs.oregonstate.edu>

  @version 1.1.2-20100521
*/

#include "xform.h"
#include "util/exception.h"
#include "dataf/keypoint.h"
// #include <cxcore.h>
#include<opencv2/opencv.hpp>
#include <stdlib.h>
#include <ctime>

/************************* Local Function Prototypes *************************/

cv::Point2f persp_xform_pt(cv::Point2f pt, cv::Mat& T)
{
    cv::Mat XY, UV;
    double xy[3] = { pt.x, pt.y, 1.0 }, uv[3] = { 0 };
    cv::Point2f rslt;

//     XY = cv::Mat(3, 1, CV_64FC1, xy, CV_AUTO_STEP);
//     UV = cv::Mat(3, 1, CV_64FC1, uv, CV_AUTO_STEP);
    XY = cv::Mat(3, 1, CV_64FC1, xy);
    UV = cv::Mat(3, 1, CV_64FC1, uv);
//     cvInitMatHeader(&XY, 3, 1, CV_64FC1, xy, CV_AUTOSTEP);
//     cvInitMatHeader(&UV, 3, 1, CV_64FC1, uv, CV_AUTOSTEP);
   // std::cout<<T.size().height;
    UV = T * XY;
// //     cvMatMul(T, &XY, &UV);
    rslt = cv::Point2f(uv[0] / uv[2], uv[1] / uv[2]);

    return rslt;
}
