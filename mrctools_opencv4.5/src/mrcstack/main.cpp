#include "opts.h"
#include <iostream>
#include <fstream>
#include "dataf/dataf.h"
#include "dataf/calibration.h"
#include "mrcimg/mrc2img.h"
#include <string>
#include<unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <mrcimg/img_util.h>

using namespace std;

int main(int argc, char **argv)
{
	struct options opts;
	if(GetOpts(argc, argv, &opts) < 0){
		return 0;
	}
	
	std::cout<<"MRC transformation."<<std::endl;
	
	util::MrcStack mrcs, mrcnew;
    
	if(!mrcs.Open(opts.input.c_str())){
		return 0;
	}
	
	mrcs.CopyToNewStack(mrcnew);
	
	cv::Mat mid = mrcs.GetStackImage(mrcs.Size()/2);
// 	IplImage* mid = mrcs.GetIplImage(mrcs.Size()/2);
	double avg, dev, min, max;    
    util::GetImgAvgSdvCutoffMinMax(mid, &avg, &dev, &min, &max);
// 	cvReleaseImage(&mid);
	mrcnew.SetHeader(util::MrcStack::MODE_FLOAT, min, avg, max);
	mrcnew.WriteHeaderToFile(opts.output.c_str());
	
	ifstream xf(opts.xform.c_str());
	EX_TRACE("Transformation parameters: %s\n", opts.xform.c_str())
    
	float a11, a12, a21, a22, shiftx, shifty;
	float cx = mrcs.Width()*.5, cy = mrcs.Height()*.5;
	int m  = opts.mode;
	
	EX_TRACE("%-4s  %-12s\t%-12s\t%-12s\t%-12s\t\n", "NO.", "AVG", "DEV", "MIN", "MAX")
	
	for(int i = 0; i < mrcs.Size(); i++){
		cv::Mat img = mrcs.GetStackImage(i);
        cv::Mat newimg = cv::Mat(img.size(), img.type());
// 		IplImage* newimg = cvCreateImage(cvGetSize(img), img->depth, img->nChannels);
		double avg, dev, min, max;
        util::GetImgAvgSdvCutoffMinMax(img, &avg, &dev, &min, &max);

		xf>>a11>>a12>>a21>>a22>>shiftx>>shifty;
		EX_TRACE("%-4d  %-12.3f\t%-12.3f\t%-12.3f\t%-12.3f\t\n", i, avg, dev, min, max)
//
// 		CvMat* xfmatrix = cvCreateMat(2, 3, CV_32F);
		cv::Mat xfmatrix = cv::Mat(2, 3, CV_32F);
		cv::Mat invmatrix = cv::Mat(2, 3, CV_32F);
		xfmatrix.at<float>(0,0)=a11; xfmatrix.at<float>(0,1)=a12; xfmatrix.at<float>(0,2)=cx+shiftx-a11*cx-a12*cy;
		xfmatrix.at<float>(1,0)=a21; xfmatrix.at<float>(1,1)=a22; xfmatrix.at<float>(1,2)=cy+shifty-a21*cx-a22*cy;
// 		CvMat* invmatrix = cvCreateMat(2, 3, CV_32F);
// 		cvmSet(xfmatrix,0,0,a11); cvmSet(xfmatrix,0,1,a12); cvmSet(xfmatrix,0,2,cx+shiftx-a11*cx-a12*cy);
// 		cvmSet(xfmatrix,1,0,a21); cvmSet(xfmatrix,1,1,a22); cvmSet(xfmatrix,1,2,cy+shifty-a21*cx-a22*cy);
// 		int avgInt = static_cast<int> (avg);
		cv::Scalar fill_color(255, 255, 255);
// 		std::cout<<"........................."<<avgScalar.val[1]<<std::endl;
		if(m==2){
			cv::warpAffine(img, newimg, xfmatrix, img.size(), cv::INTER_LANCZOS4, cv::BORDER_CONSTANT, cv::Scalar::all(avg));
		}
		else{
			cv::warpAffine(img, newimg, xfmatrix, img.size(), cv::INTER_CUBIC, cv::BORDER_CONSTANT, cv::Scalar::all(avg));
		}
		mrcnew.AppendStackImageToFile(&newimg);
// 		cvReleaseImage(&newimg);
// 		cvReleaseMat(&xfmatrix);
// 		cvReleaseImage(&img);
	}
//
       mrcs.Close();	
//     txt.close();
}

