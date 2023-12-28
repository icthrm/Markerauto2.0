#ifndef IMG_UTIL_H__
#define IMG_UTIL_H__

#include "dataf/keypoint.h"
#include "util/exception.h"
#include <vector>
#include <fstream>
#include "opencv2/core.hpp"
#include "opencv2/imgproc.hpp"
#include "opencv2/highgui.hpp"

namespace util{

// static void DrawFeature(cv::Mat& img, const util::feature& feat);
// 
// inline void DrawX(cv::Mat& img, int x, int y)
// {
//     cv::Scalar color = CV_RGB(255, 255, 255);
//     int r = 3;
//     cv::line(img, cv::Point(x - r, y - r), cv::Point(x + r, y + r), color, 3, 8, 0);
//     cv::line(img, cv::Point(x + r, y - r), cv::Point(x - r, y + r), color, 3, 8, 0);
// }
// 
// inline void DrawLine(cv::Mat& img, const util::_point& pt1, const util::_point& pt2)
// {
//     cv::Scalar color = CV_RGB(255, 255, 255);
//     cv::line(img, cv::Point(pt1.x, pt1.y), cv::Point(pt2.x, pt2.y), color);
// }
// 
// inline void DrawFeatures(cv::Mat& img, const std::vector<util::feature>& feats)
// {
// //     cv::Scalar color = CV_RGB(255, 255, 255);
// 
//     for(int i = 0; i < feats.size(); i++) {
//         DrawFeature(img, feats[i]);
//     }
// }
// 
// static void DrawFeature(cv::Mat& img, const util::feature& feat)
// {
//     cv::Scalar color = CV_RGB(255, 255, 255);
// 
//     int len, hlen, blen, start_x, start_y, end_x, end_y, h1_x, h1_y, h2_x, h2_y;
//     double scl, ori;
//     double scale = 5.0;
//     double hscale = 0.75;
//     cv::Point start, end, h1, h2;
// 
//     /* compute points for an arrow scaled and rotated by feat's scl and ori */
//     start_x = cvRound(feat.kpt.x);
//     start_y = cvRound(feat.kpt.y);
//     scl = feat.kpt.sigma;
//     ori = feat.kpt.orient;
//     len = cvRound(scl * scale);
//     hlen = cvRound(scl * hscale);
//     blen = len - hlen;
//     end_x = cvRound(len *  cos(ori)) + start_x;
//     end_y = cvRound(len * -sin(ori)) + start_y;
//     h1_x = cvRound(blen *  cos(ori + CV_PI / 18.0)) + start_x;
//     h1_y = cvRound(blen * -sin(ori + CV_PI / 18.0)) + start_y;
//     h2_x = cvRound(blen *  cos(ori - CV_PI / 18.0)) + start_x;
//     h2_y = cvRound(blen * -sin(ori - CV_PI / 18.0)) + start_y;
//     start = cv::Point(start_x, start_y);
//     end = cv::Point(end_x, end_y);
//     h1 = cv::Point(h1_x, h1_y);
//     h2 = cv::Point(h2_x, h2_y);
// 
//     cv::line(img, start, end, color, 1, 8, 0);
//     cv::line(img, end, h1, color, 1, 8, 0);
//     cv::line(img, end, h2, color, 1, 8, 0);
// }
// 
inline void SaveImage(const cv::Mat& img, const char* filename) {
    EX_TRACE("Save Image %s\n", filename)
    cv::Mat copy(img.size(), CV_8UC1);
	
	for(int y = 0; y < img.size().height; y++){
		float* src = (float*)(img.ptr()+y*img.step);
		char* start = (char*)(copy.ptr()+y*copy.step);
		for(int x = 0; x < img.size().width; x++){
			*start++ = (*src++)*255;
		}
	}

    cv::imwrite(filename, copy);
}
// 
// // typedef struct _IplImage
// // {
// //     int  nSize;             /* sizeof(IplImage) */
// //     int  ID;                /* version (=0)*/
// //     int  nChannels;         /* Most of OpenCV functions support 1,2,3 or 4 channels */
// //     int  alphaChannel;      /* Ignored by OpenCV */
// //     int  depth;             /* Pixel depth in bits: IPL_DEPTH_8U, IPL_DEPTH_8S, IPL_DEPTH_16S,
// //                                IPL_DEPTH_32S, IPL_DEPTH_32F and IPL_DEPTH_64F are supported.  */
// //     char colorModel[4];     /* Ignored by OpenCV */
// //     char channelSeq[4];     /* ditto */
// //     int  dataOrder;         /* 0 - interleaved color channels, 1 - separate color channels.
// //                                cvCreateImage can only create interleaved images */
// //     int  origin;            /* 0 - top-left origin,
// //                                1 - bottom-left origin (Windows bitmaps style).  */
// //     int  align;             /* Alignment of image rows (4 or 8).
// //                                OpenCV ignores it and uses widthStep instead.    */
// //     int  width;             /* Image width in pixels.                           */
// //     int  height;            /* Image height in pixels.                          */
// //     struct _IplROI *roi;    /* Image ROI. If NULL, the whole image is selected. */
// //     struct _IplImage *maskROI;      /* Must be NULL. */
// //     void  *imageId;                 /* "           " */
// //     struct _IplTileInfo *tileInfo;  /* "           " */
// //     int  imageSize;         /* Image data size in bytes
// //                                (==image->height*image->widthStep
// //                                in case of interleaved data)*/
// //     char *imageData;        /* Pointer to aligned image data.         */
// //     int  widthStep;         /* Size of aligned image row in bytes.    */
// //     int  BorderMode[4];     /* Ignored by OpenCV.                     */
// //     int  BorderConst[4];    /* Ditto.                                 */
// //     char *imageDataOrigin;  /* Pointer to very origin of image data
// //                                (not necessarily aligned) -
// //                                needed for correct deallocation */
// // }
// // IplImage;
// 
// // inline void SeriesSaveToFile(const cv::Mat& img, const char* filename){
// // 	std::ofstream ofs;
// // 	ofs.open(filename, std::ofstream::binary);
// // 	ofs.write(reinterpret_cast<const char*>(&(img.size().width)), sizeof(int));
// // 	ofs.write(reinterpret_cast<const char*>(&(img->height)), sizeof(int));
// // 	ofs.write(reinterpret_cast<const char*>(&(img->depth)), sizeof(int));
// // 	ofs.write(reinterpret_cast<const char*>(&(img->nChannels)), sizeof(int));
// // 	
// // 	ofs.write(reinterpret_cast<const char*>(&(img->nSize)), sizeof(int));
// // 	ofs.write(reinterpret_cast<const char*>(&(img->ID)), sizeof(int));
// // 	ofs.write(reinterpret_cast<const char*>(&(img->alphaChannel)), sizeof(int));
// // 	ofs.write(img->colorModel, sizeof(char)*4);
// // 	ofs.write(img->channelSeq, sizeof(char)*4);
// // 	ofs.write(reinterpret_cast<const char*>(&(img->dataOrder)), sizeof(int));
// // 	ofs.write(reinterpret_cast<const char*>(&(img->origin)), sizeof(int));
// // 	ofs.write(reinterpret_cast<const char*>(&(img->align)), sizeof(int));
// // 	//roi   	NULL
// // 	//maskROI   NULL
// // 	//imageId	NULL
// // 	//tileInfo	NULL
// // 	ofs.write(reinterpret_cast<const char*>(&(img->imageSize)), sizeof(int));
// // 	ofs.write(img->imageData, sizeof(char)*img->imageSize);
// // 	ofs.write(reinterpret_cast<const char*>(&(img.step)), sizeof(int));
// // 	ofs.write(reinterpret_cast<const char*>(&(img->BorderMode)), sizeof(int)*4);
// // 	ofs.write(reinterpret_cast<const char*>(&(img->BorderConst)), sizeof(int)*4);
// // 	//imageDataOrigin   point to the same place of imageData
// // 	ofs.close();
// // }
// // 
// // inline void SeriesReadFromFile(IplImage** img, const char* filename){
// // 	std::ifstream ifs;
// // 	ifs.open(filename, std::ofstream::binary);
// // 	int tmp[4];//width, height, depth, nChannels;
// // 	ifs.read(reinterpret_cast<char*>(&tmp), sizeof(int)*4);
// // 	
// // 	*img = cvCreateImage(cvSize(tmp[0], tmp[1]), tmp[2], tmp[3]);
// // 	ifs.read(reinterpret_cast<char*>(&((*img)->nSize)), sizeof(int));
// // 	ifs.read(reinterpret_cast<char*>(&((*img)->ID)), sizeof(int));
// // 	ifs.read(reinterpret_cast<char*>(&((*img)->alphaChannel)), sizeof(int));
// // 	ifs.read((*img)->colorModel, sizeof(char)*4);
// // 	ifs.read((*img)->channelSeq, sizeof(char)*4);
// // 	ifs.read(reinterpret_cast<char*>(&((*img)->dataOrder)), sizeof(int));
// // 	ifs.read(reinterpret_cast<char*>(&((*img)->origin)), sizeof(int));
// // 	ifs.read(reinterpret_cast<char*>(&((*img)->align)), sizeof(int));
// // 	//roi   	NULL
// // 	//maskROI   NULL
// // 	//imageId	NULL
// // 	//tileInfo	NULL
// // 	ifs.read(reinterpret_cast<char*>(&((*img)->imageSize)), sizeof(int));
// // 	ifs.read((*img)->imageData, sizeof(char)*(*img)->imageSize);
// // 	ifs.read(reinterpret_cast<char*>(&((*img)->widthStep)), sizeof(int));
// // 	ifs.read(reinterpret_cast<char*>(&((*img)->BorderMode)), sizeof(int)*4);
// // 	ifs.read(reinterpret_cast<char*>(&((*img)->BorderConst)), sizeof(int)*4);
// // }
// 
// inline cv::Mat GetSubImage(cv::Mat& image, cv::Rect roi)
// {
//     IplImage *result;
// 	
// 	if(roi.x < 0 || roi.x+roi.width >= image->width || roi.y < 0 || roi.y+roi.height >= image->height){
// 		return NULL;
// 	}
// 	
//     cvSetImageROI(image,roi);
// 	CvRect roi2 = cvGetImageROI(image);
// 	if(roi2.width != roi.width || roi2.height != roi.height){
// 		cvResetImageROI(image);
// 		return NULL;
// 	}
//     result = cvCreateImage( cvSize(roi.width, roi.height), image->depth, image->nChannels );
//     cvCopy(image, result);
//     cvResetImageROI(image);
//     return result;
// }
// 
// inline void ScaleImage(IplImage*& image, float scale)
// {
// 	if(scale == 1){
// 		return;
// 	}
// 	
// 	IplImage* cpy = cvCreateImage(cvSize(image->width*scale, image->height*scale), image->depth, image->nChannels);
// 	cvResize(image, cpy, CV_INTER_CUBIC);
// 	cvReleaseImage(&image);
// 	image = cpy;
// }
// 
// inline void Reversal(cv::Mat& img)
// {
//     double minimum, maximum;
//     cvMinMaxLoc(img, &minimum, &maximum);
//     double p = maximum + minimum;
//     for(int y = 0; y < img->height; y++){
// 		float* ptr = (float*)(img->imageData+y*img.step);
// 		for(int x = 0; x < img.size().width; x++){
// 			*ptr =p-*ptr;
// 			ptr++;
// 		}
// 	}
// }
// 
// // 	cvErode( img, img, NULL, 5);
// //      cvDilate(img, img, NULL, 5);
// 
// inline void MedianSmooth(IplImage *img)
// {
//     IplImage* cpy = cvCreateImage(cvGetSize(img), img->depth, img->nChannels);
//     cvSmooth(img, cpy, CV_MEDIAN, 3);
//     cvCopy(cpy, img);
//     cvReleaseImage(&cpy);
// }
// 
// #define HISTOGRAM_BIN			256//	1024
// inline void HistogramStretch(IplImage *img, int binsize = HISTOGRAM_BIN)
// {
//     cv::Mat mat(img, 0);
//     double minimum, maximum;
//     cvMinMaxLoc(img, &minimum, &maximum);
//     minimum = minimum<0?minimum:0;
//     maximum = maximum>1?maximum:1;
//     float range[] = {minimum, maximum};
//     const float* ranges = {range};
//     cv::Mat hist;
//     const int channel = 0;
// 	
//     cv::calcHist(&mat, 1, &channel, cv::Mat(), hist, 1, &binsize, &ranges);	/*warning: if has some problem, look here */
//     cv::normalize(hist, hist, 1, 0, CV_L1);
//     float* vlu = (float*)hist.data;
// 
//     for(int i = 1; i < binsize; i++) {
//         vlu[i] += vlu[i-1];
//     }
// 
//     float delta = (maximum-minimum)/(binsize);
// 
//     float res = delta/2/(maximum-minimum);
//     for(int x = 0; x < img->height; x++) {
//         for(int y = 0; y < img.size().width; y++) {
//             int loc = (int)((CV_IMAGE_ELEM(img, float, x, y)-minimum)/delta);
//             CV_IMAGE_ELEM(img, float, x, y) = (maximum-minimum)*(vlu[loc]+delta)+minimum;
//         }
//     }
// }
// #undef HISTOGRAM_BIN

inline void ConvertTo1(cv::Mat& img, bool cutoff = false)
{
    double minimum, maximum;
    cv::Scalar mean, sdv;
    cv::meanStdDev(img, mean, sdv);
    cv::minMaxLoc(img, &minimum, &maximum);

#define CUTOFF (3)

    maximum = (maximum > mean.val[0]+CUTOFF*sdv.val[0])? mean.val[0]+CUTOFF*sdv.val[0] : maximum;
    minimum = (minimum > mean.val[0]-CUTOFF*sdv.val[0])? minimum : mean.val[0]-CUTOFF*sdv.val[0];
	
	float span;
	(maximum-minimum)==0?span = 1:(span = maximum-minimum);

	for(int y = 0; y < img.size().height; y++){
		float* ptr = (float*)(img.ptr()+y*img.step);
		for(int x = 0; x < img.size().width; x++){
			if(cutoff){
				if(fabs(*ptr-mean.val[0]) <= CUTOFF*sdv.val[0]){
					*ptr = (*ptr-minimum)/span;
				}
				else {
					if(*ptr > mean.val[0]) {
						*ptr = 1;
					}
					else {
						*ptr = 0;
					}
				}
			}
			else{
				*ptr = (*ptr-minimum)/span;
			}
			ptr++;
		}
	}
	
#undef CUTOFF
}

inline void ConvertTo1(cv::Mat& img, float minimum, float maximum)
{
    float span = maximum-minimum;

    for(int y = 0; y < img.size().height; y++){
        float* ptr = (float*)(img.ptr()+y*img.step);
        for(int x = 0; x < img.size().width; x++){
            *ptr = (*ptr-minimum)/span;
            ptr++;
        }
    }
}

inline void ConvertTo1(const cv::Mat& img, float average, float deviation, float minimum, float maximum, bool cutoff = false)
{
#define CUTOFF (3)//(4)//(5)
    
    float span;
    (maximum-minimum)==0?span = 1:(span = maximum-minimum);
    
    for(int y = 0; y < img.size().height; y++){
        float* ptr = (float*)(img.ptr()+y*img.step);
        for(int x = 0; x < img.size().width; x++){
            if(cutoff){
                if(fabs(*ptr-average) <= CUTOFF*deviation) {
                    *ptr = (*ptr-minimum)/span;
                }
                else {
                    if(*ptr > average) {
                        *ptr = 1;
                    }
                    else {
                        *ptr = 0;
                    }
                }
            }
            else{
                *ptr = (*ptr-minimum)/span;
            }
            ptr++;
        }
    }
}

inline void GetImgAvgSdvCutoffMinMax(const cv::Mat& img, double* average, double* deviation, double* minimum, double* maximum, double cutoff = 3)
{
	cv::Scalar mean, sdv;

    cv::meanStdDev(img, mean, sdv);
    cv::minMaxLoc(img, minimum, maximum);

//     cvAvgSdv(img, &mean, &sdv);
//     cvMinMaxLoc(img, minimum, maximum);

	*maximum = (*maximum > mean.val[0]+cutoff*sdv.val[0])? mean.val[0]+cutoff*sdv.val[0] : *maximum;
    *minimum = (*minimum > mean.val[0]-cutoff*sdv.val[0])? *minimum : mean.val[0]-cutoff*sdv.val[0];

	*average = mean.val[0];
	*deviation = sdv.val[0];
}

// inline void CopyToUChar256(const cv::Mat& img, IplImage** cpy, int w, int b, double average, double deviation, double minimum, double maximum, bool cutoff)
// {
// #define CUTOFF (3)//(4)//(5)
//
// 	*cpy = cvCreateImage(cvSize(img.size().width, img->height), IPL_DEPTH_8U, 1);
//
// 	for(int y = 0; y < img->height; y++){
// 		float* ptr = (float*)(img->imageData+y*img.step);
// 		uchar* dst = (uchar*)((*cpy)->imageData+y*(*cpy)->widthStep);
// 		for(int x = 0; x < img.size().width; x++){
// 			if(cutoff){
// 				if(fabs(*ptr-average) <= CUTOFF*deviation) {
// 					*dst = (uchar)((*ptr-minimum)*(w-b)/((maximum-minimum)==0?1:(maximum-minimum))+b);
// 					*dst = *dst > w ? w : *dst;
// 				}
// 				else {
// 					if(*ptr > average) {
// 						*dst = w;
// 					}
// 					else {
// 						*dst = b;
// 					}
// 				}
// 			}
// 			else{
// 				*dst = (uchar)((*ptr-minimum)*(w-b)/((maximum-minimum)==0?1:(maximum-minimum))+b);
// 				*dst = *dst > w ? w : *dst;
// 			}
// 			ptr++;
// 			dst++;
// 		}
// 	}
// }
// 
// inline void CopyToUChar256(const cv::Mat& img, IplImage** cpy, int w, int b, bool cutoff = false)
// {
//     double minimum, maximum, average, deviation;
// 	
// 	GetImgAvgSdvCutoffMinMax(img, &average, &deviation, &minimum, &maximum, 3);
//     
// 	CopyToUChar256(img, cpy, w, b, average, deviation, minimum, maximum, cutoff);
// 
// #undef CUTOFF
// }


// inline void CopyToUChar256(const cv::Mat& img, IplImage** cpy, int w, int b, bool cutoff = false)
// {
//     *cpy = cvCreateImage(cvSize(img.size().width, img->height), IPL_DEPTH_8U, 1);
//     double minimum, maximum;
//     cv::Scalar mean, sdv;
//     cvAvgSdv(img, &mean, &sdv);
//     cvMinMaxLoc(img, &minimum, &maximum);
// 
// #define CUTOFF (3)//(4)//(5)
// 
//     maximum = (maximum > mean.val[0]+CUTOFF*sdv.val[0])? mean.val[0]+CUTOFF*sdv.val[0] : maximum;
//     minimum = (minimum > mean.val[0]-CUTOFF*sdv.val[0])? minimum : mean.val[0]-CUTOFF*sdv.val[0];
//     
//     for(int y = 0; y < img->height; y++){
// 		float* ptr = (float*)(img->imageData+y*img.step);
// 		uchar* dst = (uchar*)((*cpy)->imageData+y*(*cpy)->widthStep);
// 		for(int x = 0; x < img.size().width; x++){
// 			if(cutoff){
// 				if(fabs(*ptr-mean.val[0]) <= CUTOFF*sdv.val[0]) {
// 					*dst = (uchar)((*ptr-minimum)*(w-b)/((maximum-minimum)==0?1:(maximum-minimum))+b);
// 					*dst = *dst > w ? w : *dst;
// 				}
// 				else {
// 					if(*ptr > mean.val[0]) {
// 						*dst = w;
// 					}
// 					else {
// 						*dst = b;
// 					}
// 				}
// 			}
// 			else{
// 				*dst = (uchar)((*ptr-minimum)*(w-b)/((maximum-minimum)==0?1:(maximum-minimum))+b);
// 				*dst = *dst > w ? w : *dst;
// 			}
// 			ptr++;
// 			dst++;
// 		}
// 	}
// 
// #undef CUTOFF
// }

// inline void LaplaceSharpen(cv::Mat& img)
// {
//     CvMat* kernel;
//     kernel = cvCreateMat(3,3,CV_32F);
//     cvmSet(kernel, 0, 0, 0);
//     cvmSet(kernel, 0, 1, -1);
//     cvmSet(kernel, 0, 2, 0);
//     cvmSet(kernel, 1, 0, -1);
//     cvmSet(kernel, 1, 1, 5);
//     cvmSet(kernel, 1, 2, -1);
//     cvmSet(kernel, 2, 0, 0);
//     cvmSet(kernel, 2, 1, -1);
//     cvmSet(kernel, 2, 2, 0);
// 
//     cvFilter2D(img, img, kernel);
// }

}

#endif
