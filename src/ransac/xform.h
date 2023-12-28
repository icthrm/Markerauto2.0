#ifndef XFORM_H__
#define XFORM_H__

// #include <cxcore.h>
#include <opencv2/core/core.hpp>
#include<opencv2/core/core_c.h>
// #include <opencv2/highgui/highgui.hpp>
#include "dataf/dataf.h"
#include <vector>
#include <dataf/keypoint.h>

/* RANSAC error tolerance in pixels */
#define RANSAC_ERR_TOL 3

/** pessimistic estimate of fraction of inlers for RANSAC */
#define RANSAC_INLIER_FRAC_EST 0.25

/** estimate of the probability that a correspondence supports a bad model */
#define RANSAC_PROB_BAD_SUPP 0.10


extern cv::Point2f persp_xform_pt(cv::Point2f pt, cv::Mat& T);


#endif
