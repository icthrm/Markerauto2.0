#include "detector.h"
// #include "fundmental.h"
#include "mrcimg/img_util.h"
#include "ransac/xform.h"
#include <util/exception.h>
#include <iostream>
#include <stack>
#include <algorithm>
#include <cstdio>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>
#include <malloc.h>

// #include "BaseFun.h"
// #include "wavelet.h"
#include "kd_tree.h"
// #include "GaussFit.h"
// #include "option.h"

Detector Detector::detector;

//CV_IMAGE_ELEM
static inline float pixval32f( cv::Mat& img, int r, int c )
{
    return ( (float*)(img.data + img.step*r))[c];
}

static void DrawFiducialMarkerPositions(cv::Mat& img, float diameter, const std::vector<util::point2d>& fids)
{
	std::cout << "ok" << std::endl;
	cv::Scalar color = CV_RGB(0, 255, 255);
	int r = int(diameter+.5), t = int(.1*r+.5);
	if(t <= 0){
		t = 1;
	}
	for(int k = 0; k < fids.size(); k++){
		int x = int(fids[k].x+.5), y = int(fids[k].y+.5);
		cv::circle(img, cv::Point(x, y), r, color, t);
	}
}


Detector::Detector():fid_diameter(DEFAULT_FIDD)
{

}

Detector::~Detector()
{
	if(!fid_avgtmplt.empty()){
        fid_avgtmplt.release();
// 		cvReleaseImage(&fid_avgtmplt);
	}
	if(!fid_tmplt.empty()){
        fid_tmplt.release();
// 		cvReleaseImage(&fid_tmplt);
	}
}

void Detector::GenerateDiameterSpace(std::vector<float>& diavec, float dia_begin, float dia_end, int bin_size)
{
	diavec.clear();
	float step = (dia_end-dia_begin)/(bin_size-1);
	for(int i = 0; i < bin_size; i++){
		diavec.push_back(dia_begin+i*step);
	}
}


void Detector::CalculateCorralation(const cv::Mat& src, const cv::Mat& tmplt, cv::Mat* dst)
{
    *dst=cv::Mat::zeros(cv::Size(src.size().width-tmplt.size().width+1, src.size().height-tmplt.size().height+1), CV_32FC1);
    //*dst = cvCreateImage(cvSize(src.size().width-tmplt.size().width+1, src.size().height-tmplt.size().height+1), IPL_DEPTH_32F, 1);

    cv::matchTemplate(src, tmplt, *dst, cv::TM_CCORR_NORMED);//CV_TM_CCOEFF_NORMED);//
}

void Detector::FindMaxPeak(cv::Mat& corr, int seed_x, int seed_y, int idiameter, int* peak_x, int* peak_y)
{
	*peak_x = seed_x, *peak_y = seed_y;
    float max_value = corr.at<float>(seed_y, seed_x);
	//float max_value = CV_IMAGE_ELEM(corr, float, seed_y, seed_x);
	for(int x = seed_x; x < seed_x+idiameter && x < corr.size().width; x++){
		for(int y = seed_y; y < seed_y+idiameter && y < corr.size().height; y++){
			if(corr.at<float>(y, x) > max_value){
				//max_value = CV_IMAGE_ELEM(corr, float, y, x);
                max_value = corr.at<float>(y, x);
				*peak_x = x;
				*peak_y = y;
			}
		}
	}
}


float Detector::GetAverageOfMarkerPatch(cv::Mat& img, int seed_x, int seed_y, int diameter)
{
	int iradius = int(diameter*.5+.5);
	int iradius2 = iradius*iradius*.81;
	if(iradius2 == 0){
		iradius2 = iradius;
	}

	int centre_x = seed_x+iradius, centre_y = seed_y+iradius;
	float value = 0;
	int count = 0;
	for(int x = seed_x; x < seed_x+diameter && x < img.size().width; x++){
		for(int y = seed_y; y < seed_y+diameter && y < img.size().height; y++){
			if((x-centre_x)*(x-centre_x)+(y-centre_y)*(y-centre_y) < iradius2){
// 				value += CV_IMAGE_ELEM(img, float, y, x);
                value += img.at<float>(y, x);
				count++;
			}
		}
	}

	return value/count;
}


void Detector::FindFiducialMarkerPositions(cv::Mat& img, const cv::Mat& tmplt, float diameter, std::vector< util::point2d >& fids, float* positive_ratio, bool limit_contrast, bool forsubtomo)
{
#define DIA_MIN_THRE		0.36
#define DIA_MAX_THRE		1.414

	fid_diameter = diameter;

	cv::Mat corr;

	CalculateCorralation(img, tmplt, &corr);


	util::ConvertTo1(corr, false);//true);
// 	cvNormalize(corr, corr, 1, 0, CV_MINMAX);
// 	util::ConvertTo1(corr, true);
// 	util::SaveImage(corr, "corr.pgm");

	cv::Scalar corr_avg, corr_std;
	//cvAvgSdv(corr, &corr_avg, &corr_std);
    cv::meanStdDev(corr, corr_avg, corr_std);

	cv::Scalar pixel_avg, pixel_std;
	//cvAvgSdv(img, &pixel_avg, &pixel_std);
    cv::meanStdDev(img, pixel_avg, pixel_std);


#define CUTOFF_THRE       0//2.5

	if(sampling < 2500){
		corr_threshold = 0;//corr_avg.val[0]+CUTOFF_THRE*corr_std.val[0];
		pixel_threshold = 0;//pixel_avg.val[0]+CUTOFF_THRE*pixel_std.val[0];
	}
	else{
		corr_threshold = corr_avg.val[0]+CUTOFF_THRE*corr_std.val[0];
		pixel_threshold = pixel_avg.val[0]+CUTOFF_THRE*pixel_std.val[0];
	}

#undef 	CUTOFF_THRE
// 	corr_real_threshold = 0;
// 	pixel_real_threshold = 0;

	std::vector<ImgRef> refvec;
	int idia = int(diameter), iradius = int(diameter*.5);
	int sqrt2idia = int(idia*1.414);

	for(int x = 0; x < corr.size().width; x += sqrt2idia){
		for(int y = 0; y < corr.size().height; y += sqrt2idia){
			int peak_x, peak_y;
			FindMaxPeak(corr, x, y, sqrt2idia, &peak_x, &peak_y);
//             std::cout<<"x:"<<x<<"y:"<<y<<std::endl;
//             std::cout<<"peak_x:"<<peak_x<<"peak_y:"<<peak_y<<std::endl;
			//if(corr.at<float>(peak_y, peak_x) > corr_threshold && CV_IMAGE_ELEM(img, float, peak_y+iradius, peak_x+iradius) > pixel_threshold){
            if(corr.at<float>(peak_y, peak_x) > corr_threshold && img.at<float>(peak_y+iradius, peak_x+iradius) > pixel_threshold){
				refvec.push_back(ImgRef(corr, peak_x, peak_y));
			}
		}
	}

	std::vector<util::point2d> raw_fids;
	std::vector<std::pair<float, float> > scores;
	GetRawPositions(img, refvec, diameter, raw_fids, scores);

// 	if(test){
// 		for(int i = 0; i < scores.size(); i++){
// 			std::cout<<scores[i].first<<" "<<scores[i].second<<std::endl;
// 		}
//
// 		for(int i = 0; i < scores.size(); i++){
// 			std::cout<<scores[i].first*scores[i].second<<std::endl;
// 		}
//
// 		std::cout<<raw_fids.size()<<std::endl;
// 	}

	std::vector<util::point2d> new_fids;

 	RefineFiducialMarkersByGaussianDistribution(raw_fids, scores, new_fids);

	std::cout<<"Stage 1:"<<new_fids.size()<<std::endl;

// 	DrawFiducialMarkerPositions(img, new_fids);
// 	util::SaveImage(img, "img.pgm");
//
 	raw_fids.clear();
 	refvec.clear();
//
    cv::Mat img_cpy = img.clone();
 	//IplImage* img_cpy = cvCloneImage(img);
//
 	std::vector<util::point2d> candidates;
	for(int i = 0; i < new_fids.size(); i++){
		util::RECT region;
		RegionGrow(img_cpy, diameter, new_fids[i], region);
		float w = region.right-region.left;
		float h = region.bottom-region.top;
		if((w < DIA_MAX_THRE*diameter && h < DIA_MAX_THRE*diameter) && (!limit_contrast || (w >= DIA_MIN_THRE*diameter && h >= DIA_MIN_THRE*diameter))){  // && w >= DIA_MIN_THRE*diameter && h >= DIA_MIN_THRE*diameter){
			candidates.push_back(new_fids[i]);
		}
	}

	for(int i = 0; i < candidates.size(); i++){
		util::point2d fid;
		if(GetCenter(img, candidates[i], diameter, fid)){
			fids.push_back(fid);
		}
	}
// 	DrawFiducialMarkerPositions(img, fids);
// 	util::SaveImage(img, "img.pgm");
	if(forsubtomo){
		int magnif = 1;

		if(diameter <= 16){
			magnif = 16;
		}
		else if(diameter <= 32){
			magnif = 8;
		}
		else if(diameter <= 128){
			magnif = 4;
		}
		else if(diameter <= 150){
			magnif = 4;
		}
		else{
			magnif = 2;
		}
// 		std::cout<<"!"<<std::endl;
		RefinePositionsForSubTomo(img, tmplt, fids, diameter, magnif);		//WARNING
	}
	std::cout<<"Stage 2:"<<fids.size()<<std::endl;

	*positive_ratio = (float)fids.size()/new_fids.size();

	//cvReleaseImage(&img_cpy);

#undef DIA_MIN_THRE
#undef 	DIA_THRE
}


void Detector::RefineFiducialMarkersByGaussianDistribution(const std::vector< util::point2d >& raw_fids,
														   const std::vector< std::pair< float, float > >& scores, std::vector< util::point2d >& new_fids)
{
	new_fids.clear();
	std::vector<float> sscores;
	for(int i = 0; i < scores.size(); i++){
		sscores.push_back(scores[i].first*scores[i].second);
	}
	float avg = 0, stdev = 0;
	for(int i = 0; i < sscores.size(); i++){
		avg += sscores[i];
		stdev += sscores[i]*sscores[i];
	}
	avg /= sscores.size();
	stdev = sqrt(stdev/sscores.size()-avg*avg);

// #define CUTOFF_THRE  2.5//2
// 	float thre = avg+CUTOFF_THRE*stdev;
// #undef 	CUTOFF_THRE
	float thre = avg+(sampling>2500?2.5:2)*stdev;

	for(int i = 0; i < raw_fids.size(); i++){
		if(sscores[i] > thre){
			new_fids.push_back(raw_fids[i]);
		}
	}
}


void Detector::GetRawPositions(cv::Mat& img, std::vector<ImgRef>& refvec, float diameter,
							   std::vector<util::point2d>& raw_fids, std::vector<std::pair<float, float> >& scores)
{
	std::sort(refvec.begin(), refvec.end(), std::greater<ImgRef>());

	int dim = 2;
	int nk = 5;
	float radius = diameter/2;
	int nPts = refvec.size();
    ANNpoint qryPt = annAllocPt(dim);
    ANNpointArray dataPts = annAllocPts(nPts, dim);
    ANNidxArray nNIdx = new ANNidx[nk]; // near neighbor indices
    ANNdistArray dists = new ANNdist[nk]; // near neighbor distances
    ANNkd_tree* kdTree; // search structure

    for(int i = 0; i < refvec.size(); i++){
		dataPts[i][0] = refvec[i].x;
		dataPts[i][1] = refvec[i].y;
	}

	kdTree = new ANNkd_tree(dataPts, nPts, dim);
	float dists_thre = diameter;//*1.414;

	for(int i = 0; i < refvec.size(); i++){
		if(refvec[i].Value() < 0){
			continue;
		}

		qryPt[0] = refvec[i].x;
		qryPt[1] = refvec[i].y;
		kdTree->annkSearch(qryPt, nk, nNIdx, dists, 0.5);
		for(int j = 1; j < nk; j++){
			if(dists[j] < dists_thre){
				refvec[nNIdx[j]].Value() = -1;
 			}
		}

		raw_fids.push_back(util::point2d(refvec[i].x+radius, refvec[i].y+radius));
		std::pair<float, float> score;
		score.first = refvec[i].Value();

// 		int peak_x, peak_y; float half_radius = radius*.5;
// 		FindMaxPeak(img, refvec[i].x+half_radius, refvec[i].y+half_radius, half_radius, &peak_x, &peak_y);
// 		score.second = CV_IMAGE_ELEM(img, float, peak_y, peak_x);
		score.second = GetAverageOfMarkerPatch(img, refvec[i].x, refvec[i].y, int(diameter+.5));
		scores.push_back(score);
	}

	delete [] nNIdx;
    delete [] dists;
    delete kdTree;
    annDeallocPts(dataPts);
    annDeallocPt(qryPt);
}

bool Detector::RegionGrow(cv::Mat& img, float diameter, const util::point2d& seed, util::RECT& region)
{
	int idia = int(diameter*1.2533+.5);		//1.2533 = sqrt(pi/2)

	cv::Rect roi;
	roi.x = int(seed.x - idia*.5+.5);
	roi.y = int(seed.y - idia*.5+.5);
	roi.width = idia;
	roi.height = idia;
	cv::Mat sub = util::GetSubImage(img, roi);

	if(sub.empty()){
		return false;
	}

	cv::Scalar pixel_avg, pixel_std;
	cv::meanStdDev(sub, pixel_avg, pixel_std);

// 	std::cout<<pixel_avg.val[0]<<" "<<pixel_std.val[0]<<" "<<CV_IMAGE_ELEM(img, float, int(seed.y), int(seed.x))<<std::endl;


	int seed_y = seed.y;
	int seed_x = seed.x;
// 	float max_value = ( (float*)(img->imageData + img->widthStep*seed_y))[seed_x];//seed_y, seed_x);

    double d = pixel_avg.val[0]+.5*pixel_std.val[0];//max_value*RG_THRE;
    //std::cout<<"d:"<<d<<std::endl;
    std::stack<util::point2d> seedd;
    seedd.push(seed);

    region.left = seed.x;
    region.right = seed.x;
    region.top = seed.y;
    region.bottom = seed.y;

    int width = img.size().width;
    int height = img.size().height;
    while(!seedd.empty()){
		util::point2d point = seedd.top();
        seedd.pop();
		int intx = int(point.x);
		int inty = int(point.y);

		if(!(intx > 0 && intx < width-1 && inty > 0 && inty < height-1)){
			break;
		}

		((float*)(img.data + inty*img.step))[intx] = 0;
		float value = ((float*)(img.data + inty*img.step))[intx-1];		//(x-1, y)

		if(value > d){
			if(point.x-1 < region.left){
				region.left = point.x-1;
			}
			util::point2d tmp(point.x-1, point.y);
			seedd.push(tmp);
		}

		value = ((float*)(img.data + inty*img.step))[intx+1];			//(x+1, y)
		if(value > d){
			if(point.x+1 > region.right){
				region.right = point.x+1;
			}
			util::point2d tmp(point.x+1, point.y);
			seedd.push(tmp);
		}

		value = ((float*)(img.data + (inty-1)*img.step))[intx];		//(x, y-1)
		if(value > d){
			if(point.y-1 < region.top){
				region.top = point.y-1;
			}
			util::point2d tmp(point.x, point.y-1);
			seedd.push(tmp);
		}

		value = ((float*)(img.data + (inty+1)*img.step))[intx];		//(x, y+1)
		if(value > d){
			if(point.y+1 > region.bottom){
				region.bottom = point.y+1;
			}
			util::point2d tmp(point.x, point.y+1);
			seedd.push(tmp);
		}
    }

    if(region.left < 0 || region.right >= width || region.top < 0 || region.bottom >= height){
        return false;
	}

    return true;

#undef BG_FID_THRE
#undef RG_THRE
}

bool Detector::GetCenter(cv::Mat& img, const util::point2d& seed, float diameter, util::point2d& centre)
{
	int idia = int(diameter+.5);		//1.2533 = sqrt(pi/2)

	cv::Rect roi;
	roi.x = int(seed.x - idia*.5+.5);
	roi.y = int(seed.y - idia*.5+.5);
	roi.width = idia;
	roi.height = idia;
	cv::Mat sub = util::GetSubImage(img, roi);

	if(sub.empty()){
		return false;
	}
	std::vector<float> values;
	for(int i = 0; i < sub.size().width; i++){
		for(int j = 0; j < sub.size().height; j++){
			float value = ((float*)(sub.data + j*sub.step))[i];
			values.push_back(value);
		}
	}
	std::sort(values.begin(), values.end(), std::greater<float>());
    float d = values[values.size()*.785*.81];		//0.785 = pi/4

	float sx = 0, sy = 0, sum = 0;

	for(int i = 0; i < sub.size().width; i++){
		for(int j = 0; j < sub.size().height; j++){
			float value = ((float*)(sub.data + j*sub.step))[i];
            if(value > d){
                sx += i*value;
                sy +=j*value;
                sum += value;
            }
		}
	}


	if(sum != 0){
		sx /= sum;
		sy /= sum;
		centre.x = roi.x+sx;
		centre.y = roi.y+sy;
		return true;
	}
	return false;
}


void Detector::CreateTemplate(cv::Mat* tmplt, float diameter)
{
#define TMP_ROUND(x)		int(x+0.5)
#define TEMPLATE_RATIO		16
	float radius = diameter*.5;
	int dia = TMP_ROUND(radius*TEMPLATE_RATIO);
    *tmplt = cv::Mat::zeros(cv::Size(TMP_ROUND(diameter+1), TMP_ROUND(diameter+1)), CV_32FC1);
    cv::Mat tmp(cv::Size(2*dia+1, 2*dia+1), CV_32FC1, cv::Scalar::all(-.5f));
    cv::circle(tmp, cv::Point(dia, dia), dia, CV_RGB(.5f, .5f, .5f), -1, 8);//CV_RGB(1, 1, 1), -1, 8);//
	cv::resize(tmp, *tmplt, cv::Size(TMP_ROUND(diameter+1), TMP_ROUND(diameter+1)), 0, 0, cv::INTER_CUBIC);

}



void Detector::GaussianSmoothBasedOnSampling(cv::Mat& img, int sampling)
{
	if(sampling > 2500){

	}
	else if(sampling > 1250){
		//cvSmooth(img, img, CV_GAUSSIAN, 3);
        cv::GaussianBlur(img, img, cv::Size(3, 3), 0, 0);
	}
	else{
		cv::GaussianBlur(img, img, cv::Size(5, 5), 0, 0);
        //cvSmooth(img, img, CV_GAUSSIAN, 5);
	}
}

void Detector::CreateAverageTemplate(cv::Mat& img, const std::vector< util::point2d >& fids, float diameter, cv::Mat* tmplt, bool zoom_out, int magnif)
{
	int TEMPLATE_BOUND = 4;
	float radius = diameter*.5;
	int dia = TMP_ROUND(radius*magnif);
// 	std::cout << img << std::endl;
// 	*tmplt = cvCreateImage(cvSize(TMP_ROUND(diameter+11), TMP_ROUND(diameter+11)), IPL_DEPTH_32F, 1);

    cv::Mat tmp(cv::Size(2*dia+1, 2*dia+1), CV_32FC1, cv::Scalar::all(-0.0f));
	//IplImage* tmp = cvCreateImage(cvSize(2*dia+1, 2*dia+1), IPL_DEPTH_32F, 1);
	//cvFillImage(tmp, 0.0f);

    for(int i = 0; i < fids.size(); i++){
		cv::Rect rect;
		rect.x = TMP_ROUND(fids[i].x-radius-TEMPLATE_BOUND);
		rect.y = TMP_ROUND(fids[i].y-radius-TEMPLATE_BOUND);
		rect.width = TMP_ROUND(diameter+TEMPLATE_BOUND*2);
		rect.height = TMP_ROUND(diameter+TEMPLATE_BOUND*2);

		*tmplt = util::GetSubImage(img, rect);
		if(tmplt->empty()){
			continue;
		}

		cv::Mat tmp_patch=cv::Mat::zeros(cv::Size(rect.width*magnif, rect.height*magnif), CV_32FC1);

// 		IplImage* tmp_patch = cvCreateImage(cvSize(rect.width*magnif, rect.height*magnif), IPL_DEPTH_32F, 1);
        //cvResize(*tmplt, tmp_patch, CV_INTER_CUBIC);
 		cv::resize(*tmplt, tmp_patch, tmp_patch.size(), 0, 0, cv::INTER_CUBIC);

		int fid_xr = int(magnif*fids[i].x)-rect.x*magnif;
		int fid_yr = int(magnif*fids[i].y)-rect.y*magnif;
		rect.x = fid_xr-dia;
		rect.y = fid_yr-dia;
		rect.width = 2*dia+1;
		rect.height = 2*dia+1;
		*tmplt = util::GetSubImage(tmp_patch, rect);
		cv::add(*tmplt, tmp, tmp);
    }
// std::cout << fids.size() << std::endl;

//     util::Reversal(tmp);
// 	util::ConvertTo1(tmp, false);
// 	util::SaveImage(tmp, "tmp.pgm");
    tmp.convertTo(tmp, -1, 1.f/(int)fids.size());
	if(!zoom_out){
		*tmplt = cv::Mat::zeros(cv::Size(TMP_ROUND(diameter+1), TMP_ROUND(diameter+1)), CV_32FC1);
		cv::resize(tmp, *tmplt, cv::Size(TMP_ROUND(diameter+1), TMP_ROUND(diameter+1)), 0, 0, cv::INTER_CUBIC);
	}
	else{
		*tmplt = tmp;
	}
}


void Detector::RefinePositionsForSubTomo(cv::Mat& img, const cv::Mat& tmplt, std::vector< util::point2d >& fids, float diameter, int magnif)
{
// 	diameter = 5.184*2;
	float radius = diameter*.5;
	int dia = TMP_ROUND(radius*magnif);

	int TEMPLATE_BOUND = diameter*.1;
	if(TEMPLATE_BOUND<=3){
		TEMPLATE_BOUND = 3;
	}
// std::cout << img << std::endl;
	cv::Mat cpytmplt = cv::Mat(cv::Size(2*dia+1, 2*dia+1), tmplt.depth(), tmplt.channels());
	cv::resize(tmplt, cpytmplt, cv::Size(2*dia+1, 2*dia+1), 0, 0, cv::INTER_CUBIC);
	while(true){
		cv::Mat avgtmplt;
		CreateAverageTemplate(img, fids, diameter, &avgtmplt, true, magnif);
// 				std::cout << diameter<<"       "<<magnif << std::endl;
		cv::add(cpytmplt, avgtmplt, avgtmplt);
		std::vector< util::point2d > finfids;
		float res = 0;
// std::cout << avgtmplt.size() << std::endl;
		for(int i = 0; i < fids.size(); i++){
// 	  std::cout << fids[i].x<<"  "<< fids[i].y<< std::endl;
			cv::Rect rect;
			rect.x = TMP_ROUND(fids[i].x-radius-TEMPLATE_BOUND);
			rect.y = TMP_ROUND(fids[i].y-radius-TEMPLATE_BOUND);
			rect.width = TMP_ROUND(diameter+TEMPLATE_BOUND*2);
			rect.height = TMP_ROUND(diameter+TEMPLATE_BOUND*2);
// std::cout << rect.x<<"  "<< rect.y<<"  "<< rect.width<<"  "<< rect.height<<"  " << std::endl;
			cv::Mat fidtmplt = util::GetSubImage(img, rect);
			if(fidtmplt.empty()){
				continue;
			}
// 			std::cout << fidtmplt << std::endl;
			cv::Mat tmp_patch = cv::Mat::zeros(cv::Size(rect.width*magnif, rect.height*magnif), CV_32FC1);
			cv::resize(fidtmplt, tmp_patch, cv::Size(rect.width*magnif, rect.height*magnif), 0, 0, cv::INTER_CUBIC);
// 			cvReleaseImage(&fidtmplt);

			cv::Mat corr;
			CalculateCorralation(tmp_patch, avgtmplt, &corr);

// 			std::cout << corr.size() << std::endl;

			int peak_x, peak_y;
			FindMaxPeak(corr, 0, 0, corr.size().width, &peak_x, &peak_y);
			peak_x += dia; peak_y += dia;
// std::cout <<"peak_x:"<< peak_x<<"peak_y:"<<peak_y << std::endl;
// std::cout<<(float)peak_x/magnif<<std::endl;
			util::point2d pt;
			pt.x = rect.x+(float)peak_x/magnif;
			pt.y = rect.y+(float)peak_y/magnif;

			res += sqrt((fids[i].x-pt.x)*(fids[i].x-pt.x)+ (fids[i].y-pt.y)*(fids[i].y-pt.y));

			finfids.push_back(pt);

// 			cvReleaseImage(&corr);
// 			cvReleaseImage(&tmp_patch);
		}

// 		cvReleaseImage(&avgtmplt);
		fids = finfids;

		res /= fids.size();

// 		std::cout << "*********" << res << std::endl;

		if(res < 0.173){
			break;
		}
	}

// 	cvReleaseImage(&cpytmplt);

#undef TMP_ROUND
#undef TEMPLATE_RATIO
}


float Detector::InitEstimateDiameter(const cv::Mat& img)
{

	cv::Size size(img.size().width, img.size().height);
    cv::Mat cpy = cv::Mat::zeros(size, CV_32FC1);
    //cv::Mat cpy = cvCreateImage(cvSize(img->width, img->height), IPL_DEPTH_32F, 1);
    cv::resize(img, cpy, size, 0, 0, cv::INTER_CUBIC);
	//cvResize(img, cpy, CV_INTER_CUBIC);
    cv::GaussianBlur(cpy, cpy, cv::Size(9, 9), 0, 0);
	//cvSmooth(cpy, cpy, CV_GAUSSIAN, 9);

	std::vector<float> whv;

	int count = 0;
    //return 5;

	while(count < 6){
		int max_i, max_j;
		float max_value = -9999;
		for(int i = 0; i < cpy.size().width; i++){
			for(int j = 0; j < cpy.size().height; j++){
				float cpy_ij = ((float*)(cpy.ptr() + j*cpy.step))[i];
				if(cpy_ij > max_value){
					max_value = cpy_ij;
					max_i = i; max_j = j;
				}
			}
		}

		util::RECT region;
		BlindRegionGrow(cpy, max_i, max_j, region);

		float width = region.right-region.left;
		float height = region.bottom-region.top;

		if(width < 0.7*height || 0.7*width > height){
			continue;
		}

		if(region.left < img.size().width*.025 || region.right > img.size().width*.975 || region.top < img.size().width*.025 || region.bottom > img.size().width*.975){
			continue;
		}


		whv.push_back(width*.5+height*.5);

		std::cout<<region.right-region.left<<"\t"<<region.bottom-region.top<<std::endl;

		count++;
	}

	//cvReleaseImage(&cpy);

	std::sort(whv.begin(), whv.end());

 	return (whv[2]+whv[3])*.5;
}

void Detector::BlindRegionGrow(cv::Mat& img, int seed_x, int seed_y, util::RECT& region)
{
	float max_value = ( (float*)(img.ptr() + img.step*seed_y))[seed_x];//seed_y, seed_x);

    double d = (max_value*.65+.5*.35);//max_value*RG_THRE;
    std::stack<util::point2d> seedd;
	util::point2d seed;
	seed.x = seed_x; seed.y = seed_y;

    seedd.push(seed);


    region.left = seed.x;
    region.right = seed.x;
    region.top = seed.y;
    region.bottom = seed.y;

    int width = img.size().width;
    int height = img.size().height;
    while(!seedd.empty()){
		util::point2d point = seedd.top();
        seedd.pop();
		int intx = int(point.x);
		int inty = int(point.y);

		if(!(intx >= 0 && intx <= width-1 && inty >= 0 && inty <= height-1)){
			continue;
		}

		((float*)(img.ptr() + inty*img.step))[intx] = 0;

		float value = 0;
		if(intx-1 >= 0){
			value= ((float*)(img.ptr() + inty*img.step))[intx-1];		//(x-1, y)
		}

		if(value > d){
			if(point.x-1 < region.left){
				region.left = point.x-1;
			}
			util::point2d tmp(point.x-1, point.y);
			seedd.push(tmp);
		}

		if(intx+1 <= width-1){
			value = ((float*)(img.ptr() + inty*img.step))[intx+1];			//(x+1, y)
		}

		if(value > d){
			if(point.x+1 > region.right){
				region.right = point.x+1;
			}
			util::point2d tmp(point.x+1, point.y);
			seedd.push(tmp);
		}

		if(inty-1 >= 0){
			value = ((float*)(img.ptr() + (inty-1)*img.step))[intx];		//(x, y-1)
		}

		if(value > d){
			if(point.y-1 < region.top){
				region.top = point.y-1;
			}
			util::point2d tmp(point.x, point.y-1);
			seedd.push(tmp);
		}

		if(inty+1 <= height-1){
			value = ((float*)(img.ptr() + (inty+1)*img.step))[intx];		//(x, y+1)
		}

		if(value > d){
			if(point.y+1 > region.bottom){
				region.bottom = point.y+1;
			}
			util::point2d tmp(point.x, point.y+1);
			seedd.push(tmp);
		}
    }
}

bool Detector::DetermineDiameter(const cv::Mat& img, float diameter, float* refined_diameter, cv::Mat* avgtmplt)
{
#define OUTPUT_THRESHOLD	0.5

	EX_TRACE("Input diameter value = %.2f\n", diameter)

	float dbdia_begin;
	float dbdia_end;
	int bin_size;

	if(diameter < 0){
		EX_TRACE("Program will estimate the initial value of diameter for the user!\n")
		diameter = InitEstimateDiameter(img);
		diameter *= .9;
		EX_TRACE("Initial estimation of the diameter = %.2f\n", diameter)

		if(diameter < 30){
			dbdia_begin = diameter*1;
			dbdia_end = diameter*2;
			bin_size = 11;
		}
		else{
			dbdia_begin = diameter*1;
			dbdia_end = diameter*1.8;
			bin_size = int(diameter*.4)+1;
		}
	}
	else if(diameter < 20){
		dbdia_begin = diameter*.5;
		dbdia_end = diameter*1.5;
		bin_size = 11;
	}
	else if(diameter >= 20 && diameter < 30){					//diameter >= 30
		dbdia_begin = diameter-10;
		dbdia_end = diameter+10;
		bin_size = 11;
	}
	else if(diameter >= 30 && diameter < 40){					//diameter >= 30
		dbdia_begin = diameter-15;
		dbdia_end = diameter+15;
		bin_size = 16;
	}
	else if(diameter >= 40 && diameter < 60){
		dbdia_begin = diameter-18;
		dbdia_end = diameter+18;
		bin_size = 19;
	}
	else if(diameter >= 60 && diameter < 80){
		dbdia_begin = diameter-26;
		dbdia_end = diameter+26;
		bin_size = 27;
	}
	else if(diameter >= 80){
		dbdia_begin = diameter-32;
		dbdia_end = diameter+32;
		bin_size = 33;
	}

	EX_TIME_BEGIN("Refining the fiducial marker diameter...")

	while(true){
		std::vector<float> diavec;
		std::vector<float> ccvec;
		std::vector<cv::Mat> tmpltvec;
		GenerateDiameterSpace(diavec, dbdia_begin, dbdia_end, bin_size);
//         for (auto it = diavec.begin(); it != diavec.end(); it++)
//             std::cout << *it << " ";

 		bin_size = 11;
// 		diavec[0] = 8;
		for(int i = 0; i < diavec.size(); i++){
			std::vector<util::point2d> fids;
			cv::Mat manutmplt;
			CreateTemplate(&manutmplt, diavec[i]);

            cv::Size size(img.size().width, img.size().height);
            cv::Mat cpy=cv::Mat::zeros(size, CV_32FC1);
			//IplImage* cpy = cvCreateImage(cvSize(img.size().width, img.size().height), IPL_DEPTH_32F, 1);
		// 	cvSmooth(img, cpy, CV_GAUSSIAN, 5);
			cv::resize(img, cpy, size, 0, 0, cv::INTER_CUBIC);

			int csampling = int(cpy.size().width*cpy.size().height/(diavec[i]*diavec[i])*.5);		//in case of electorn noise in large scale


			GaussianSmoothBasedOnSampling(cpy, csampling);
			sampling = csampling;
			float positive_value;
			FindFiducialMarkerPositions(cpy, manutmplt, diavec[i], fids, &positive_value, true, false);


			cv::Mat tmplt, ccv;

			int magnif;
			if(diavec[i] <= 16){
				magnif = 16;
			}
			else if(diavec[i] <= 32){
				magnif = 8;
			}
			else if(diavec[i] <= 128){
				magnif = 4;
			}
			else if(diavec[i] <= 150){
				magnif = 4;
			}
			else{
				magnif = 2;
			}

			CreateAverageTemplate(cpy, fids, diavec[i], &tmplt, false);

  			CalculateCorralation(manutmplt, tmplt, &ccv);
			float ncc = pixval32f(ccv, 0, 0);
// 			float ncc = GetCorrelationScore(cpy, fids, tmplt);
			if(tmplt.size().width > 16){
				ccvec.push_back(ncc);//*positive_value);//(positive_value > 0.7 ? 1 : sqrt(positive_value)));
			}
			else{
				ccvec.push_back(ncc*positive_value);
			}

// 			util::SaveImage(tmplt, "tmplt.pgm");
			std::cout<<"diameter(value "<<diavec[i]<<"): cc score = "<<ncc<<", acceptive ratio = "<<positive_value<<std::endl;

			tmpltvec.push_back(tmplt);
		}

		int max_idx = 0;
		float max_cc = ccvec[0];
		for(int i = 1; i < ccvec.size(); i++){
			if(ccvec[i] > max_cc){
				max_idx = i;
				max_cc = ccvec[i];
			}
		}

		if(max_idx > 0 && max_idx < diavec.size()-1){
			dbdia_begin = diavec[max_idx-1];
			dbdia_end = diavec[max_idx+1];
		}
		else if(max_idx == 0){
			dbdia_begin = diavec[0];
			dbdia_end = diavec[1];
		}
		else{
			dbdia_begin = diavec[diavec.size()-2];
			dbdia_end = diavec[diavec.size()-1];
		}

		if(dbdia_end-dbdia_begin < OUTPUT_THRESHOLD){
			*avgtmplt = tmpltvec[max_idx];
// 			for(int i = 0; i < max_idx-1; i++){
// 				//cvReleaseImage(&(tmpltvec[i]));
// 			}
// 			for(int i = max_idx+1; i < tmpltvec.size(); i++){
// 				//cvReleaseImage(&(tmpltvec[i]));
// 			}
			break;
		}

// 		for(int i = 0; i < tmpltvec.size(); i++){
// 			//cvReleaseImage(&(tmpltvec[i]));
// 		}
 	}

	*refined_diameter = (dbdia_begin+dbdia_end)*.5;

	std::cout<<"Totally sampling number will be "<<int(img.size().width*img.size().height/ (*refined_diameter* *refined_diameter)*.5)<<std::endl;
 	EX_TIME_END("Refined fiducial marker diameter: value = %.2f", *refined_diameter)

#undef OUTPUT_THRESHOLD
}

void Detector::SetFidDiameter(float __fid_diameter)
{
	fid_diameter = __fid_diameter;
}

void Detector::SetFidRadius(float fid__radius)
{
	fid_radius = fid__radius;
}

void Detector::SetWaveimg(cv::Mat &wave)
{
	wave_img = wave.clone();
}

void Detector::SetSampling(cv::Mat& img, float __fid_diameter)
{
	sampling = int(img.size().width*img.size().height/(__fid_diameter*__fid_diameter)*.5);
}

void Detector::SetFidAverageTemplate(const cv::Mat& __fid_tmplt)
{
    //fid_tmplt.copyTo(__fid_tmplt);
    fid_tmplt=__fid_tmplt.clone();
	//fid_tmplt = cvCloneImage(__fid_tmplt);
}

void Detector::Process(cv::Mat& img, std::vector< util::point2d >& fids, bool use_avgtmplt)
{
	EX_TIME_BEGIN("\nCompute Positions of Fiducial Markers (X Y)")

	cv::Mat tmplt;
	if(!use_avgtmplt){
		CreateTemplate(&tmplt, fid_diameter);
	}
	else{
        tmplt=fid_avgtmplt.clone();
	}

	GaussianSmoothBasedOnSampling(img, sampling);

	float p_v;
	FindFiducialMarkerPositions(img, tmplt, fid_diameter, fids, &p_v, false, true);

    EX_TIME_END("Totally found %ld fiducial markers", fids.size())
}

float Detector::FindLocalPeak(cv::Mat& img, const cv::Mat& tmplt, const util::point2d& seed, float roi_radio, util::point2d* loc)
{
	cv::Rect rect;
	rect.x = seed.x-roi_radio;
	rect.y = seed.y-roi_radio;
	rect.width = roi_radio*2;
	rect.height = roi_radio*2;

	cv::Mat reg = util::GetSubImage(img, rect);
	if(reg.empty()){
		return -1;
	}
	cv::Mat corr;
	CalculateCorralation(reg, tmplt, &corr);
// 	std::cout<<corr<<std::endl;
// 	std::cout<<reg<<std::endl;
// 	std::cout<<tmplt<<std::endl;
// 	util::ConvertTo1(corr, false);//true);

	int peak_x, peak_y;
	FindMaxPeak(corr, 0, 0, corr.size().width, &peak_x, &peak_y);
// 	float max_value = CV_IMAGE_ELEM(corr, float, peak_y, peak_x);
    float max_value = corr.at<float>(peak_y, peak_y);

	float radius = tmplt.size().width*.5;
	loc->x = rect.x+peak_x+radius;
	loc->y = rect.y+peak_y+radius;

	util::point2d rloc;

	bool success = GetCenter(img, *loc, (float)tmplt.size().width, rloc);
	*loc = rloc;

	if(success){
		return max_value;
	}
	else{
		return -1;
	}
}

void Detector::InitLocalDetector(float diameter, bool use_avgtmplt)
{
	if(use_avgtmplt){
		util::SeriesReadFromFile(&detector.fid_avgtmplt, "avgtmplt");

        detector.fid_tmplt = detector.fid_avgtmplt.clone();
		//detector.fid_tmplt = cvCloneImage(detector.fid_avgtmplt);
		util::ConvertTo1(detector.fid_tmplt, false);

		detector.SetFidDiameter(detector.fid_avgtmplt.size().width);				//replace the fid_diameter in detector by average template
	}
	else{
		detector.SetFidDiameter(diameter);
		detector.CreateTemplate(&(detector.fid_avgtmplt), diameter);
	}
}

float Detector::LocalMarkerDetect(cv::Mat& img, const util::point2d& seed, util::point2d* loc)
{
#define ROI_RADIUS_RATIO		.6f
	return detector.FindLocalPeak(img, detector.fid_tmplt, seed, detector.fid_diameter*ROI_RADIUS_RATIO, loc);
}


void Detector::PredictMissFiducial(const util::point2d& seed, HSetVector& hsets, int z_idx1, int z_idx2, std::vector< util::point2d >& slocs)
{
	slocs.clear();
	h_set* txset = hsets.GetHSetWithIdx(z_idx1, z_idx2);
	if(txset){
// 		std::cout<<txset->h[0]<<std::endl;
		for(int m = 0; m < txset->size(); m++){
			cv::Mat txm = txset->h[m];
			if(txm.empty()){
				return;
			}
			cv::Point2f tmp, p;
//             CvPoint2D32f tmp, p;
			p.x = seed.x; p.y = seed.y;
			util::_point tmpp;
			tmp = persp_xform_pt(p, txm);
			tmpp.x = tmp.x; tmpp.y = tmp.y;
			slocs.push_back(tmpp);
		}
	}
	else{
		txset = hsets.GetHSetWithIdx(z_idx2, z_idx1);

		if(!txset){
			return;
		}

		for(int m = 0; m < txset->size(); m++){
			cv::Mat txm = txset->h[m];
			cv::Mat txinv = cv::Mat::zeros(3, 3, CV_64FC1);
			cv::invert(txm, txinv);
			cv::Point2f tmp, p;
			p.x = seed.x; p.y = seed.y;
			util::_point tmpp;
			tmp = persp_xform_pt(p, txinv);
			tmpp.x = tmp.x; tmpp.y = tmp.y;
			slocs.push_back(tmpp);
		}
	}
}

void Detector::LocalDetectorMain(util::MrcStack& mrcs, util::TrackSpace& trackspace, HSetVector& hsets, float diameter, util::FiducialStack* fidsk, util::ImgMatchVector* imvector, const std::vector<int>& m_index)
{
	fidsk->ReSize(m_index.size());
    fidsk->SetWxH(mrcs.Width(), mrcs.Height());

	int zidx = trackspace.Size()/2;
	Detector::InitLocalDetector(diameter);

	std::vector<LocalSeed>* lsearchinfos = new std::vector<LocalSeed>[trackspace.Size()];

	for(util::TrackSpace::Iterator itr = trackspace.Z_Iterator(zidx); !itr.IsNULL(); itr++){
		util::TrackSpace::TrackNode node_itr = util::TrackSpace::TrackNode(itr);

		for(; !node_itr.IsNULL(); node_itr++){
			util::TrackSpace::TrackNode tmp = node_itr;
			tmp++;

			if(tmp.IsNULL()){
				if(node_itr.Z()+1 == trackspace.Size()){
					continue;
				}
				LocalSeed lsed;
				lsed.z_bp.x = node_itr.X(); lsed.z_bp.y = node_itr.Y();
				lsed.z_begin = node_itr.Z(); lsed.z_next = -1;
				lsed.is_final = false;
				lsed.z_img_idx = node_itr.Z()+1;
				lsearchinfos[lsed.z_img_idx].push_back(lsed);
			}
			else if(tmp.Z() != node_itr.Z()+1){
				LocalSeed lsed;
				lsed.z_bp.x = node_itr.X(); lsed.z_bp.y = node_itr.Y();
				lsed.z_np.x = tmp.X(); lsed.z_np.y = tmp.Y();
				lsed.z_begin = node_itr.Z(); lsed.z_next = tmp.Z();
				if(lsed.z_begin+2 != lsed.z_next){
					lsed.is_final = false;
				}
				else{
					lsed.is_final = true;
				}
				lsed.z_img_idx = node_itr.Z()+1;

				lsearchinfos[lsed.z_img_idx].push_back(lsed);
			}
		}

		node_itr = util::TrackSpace::TrackNode(itr);

		for(; !node_itr.IsNULL(); node_itr--){
			util::TrackSpace::TrackNode tmp = node_itr;
			tmp--;

			if(tmp.IsNULL()){
				if(node_itr.Z() == 0){
					continue;
				}
				LocalSeed lsed;
				lsed.z_bp.x = node_itr.X(); lsed.z_bp.y = node_itr.Y();
				lsed.z_begin = node_itr.Z(); lsed.z_next = -1;
				lsed.is_final = false;
				lsed.z_img_idx = node_itr.Z()-1;
				lsearchinfos[lsed.z_img_idx].push_back(lsed);
			}

			else if(tmp.Z() != node_itr.Z()-1){
				LocalSeed lsed;
				lsed.z_bp.x = node_itr.X(); lsed.z_bp.y = node_itr.Y();
				lsed.z_np.x = tmp.X(); lsed.z_np.y = tmp.Y();
				lsed.z_begin = node_itr.Z(); lsed.z_next = tmp.Z();
				if(lsed.z_begin-2 != lsed.z_next){
					lsed.is_final = false;
				}
				else{
					lsed.is_final = true;
				}
				lsed.z_img_idx = node_itr.Z()-1;

				lsearchinfos[lsed.z_img_idx].push_back(lsed);
			}
		}
	}

	for(int i = 0/*0 zidx+1*/; i < trackspace.Size(); i++){
		int index=m_index[i];
		cv::Mat slice = mrcs.GetStackImage(index);
// 		cv::Mat slice = mrcs.GetStackImage(i);
		util::Reversal(slice);
		util::ConvertTo1(slice);
		util::MedianSmooth(slice);

		for(int j = 0; j < lsearchinfos[i].size(); j++){
			LocalSeed& lsed = lsearchinfos[i][j];

			std::vector<util::point2d> slocs;

			PredictMissFiducial(lsed.z_bp, hsets, lsed.z_begin, lsed.z_img_idx, slocs);

			h_set* txset = hsets.GetHSetWithIdx(lsed.z_begin, lsed.z_img_idx);

			util::point2d nxloc;
			float max_score = -999;
			for(int m = 0; m < slocs.size(); m++){
				util::point2d p;
				float score = LocalMarkerDetect(slice, slocs[m], &p);
				if(score > max_score){
					max_score = score;
					nxloc = p;
				}
			}
			if(max_score > 0.85){								//WARNING threshold
				fidsk->V(i).push_back(nxloc);
			}
			else{
				continue;
			}

            util::img_match* imatch;
			bool noex;

			if(lsed.is_final){
				noex = true;
				imatch = imvector->GetMatchSetWithIdx(lsed.z_begin, lsed.z_img_idx, noex);
				if(!imatch){
					imvector->MallocNewMatch();
					imatch = &((*imvector)[imvector->Size()-1]);
					imatch->idx1 = lsed.z_begin; imatch->idx2 = lsed.z_img_idx;
				}
				if(noex){
					imatch->pairs.push_back(std::make_pair(lsed.z_bp, nxloc));
				}
				else{
					imatch->pairs.push_back(std::make_pair(nxloc, lsed.z_bp));
				}

				noex = true;
				imatch = imvector->GetMatchSetWithIdx(lsed.z_img_idx, lsed.z_next, noex);
				if(!imatch){
					imvector->MallocNewMatch();
					imatch = &((*imvector)[imvector->Size()-1]);
					imatch->idx1 = lsed.z_img_idx; imatch->idx2 = lsed.z_next;
				}
				if(noex){
					imatch->pairs.push_back(std::make_pair(nxloc, lsed.z_np));
				}
				else{
					imatch->pairs.push_back(std::make_pair(lsed.z_np, nxloc));
				}
			}
			else if(lsed.z_img_idx == trackspace.Size()-1 || lsed.z_img_idx == 0){
				noex = true;
				imatch = imvector->GetMatchSetWithIdx(lsed.z_begin, lsed.z_img_idx, noex);
				if(!imatch){
					imvector->MallocNewMatch();
					imatch = &((*imvector)[imvector->Size()-1]);
					imatch->idx1 = lsed.z_begin; imatch->idx2 = lsed.z_img_idx;
				}
				if(noex){
					imatch->pairs.push_back(std::make_pair(lsed.z_bp, nxloc));
				}
				else{
					imatch->pairs.push_back(std::make_pair(nxloc, lsed.z_bp));
				}

				imatch->pairs.push_back(std::make_pair(lsed.z_bp, nxloc));
			}
			else{
				float sign = 1;

				if(lsed.z_begin < zidx){
					sign = -1;
				}

				for(int step = 1; step <= 2 && lsed.z_img_idx+step < trackspace.Size() && lsed.z_img_idx-step >= 0; step++){
					PredictMissFiducial(nxloc, hsets, lsed.z_img_idx, lsed.z_img_idx+step*sign, slocs);
					std::vector<util::point2d> candidates;
					for(int m = 0; m < slocs.size(); m++){
						for(util::TrackSpace::Iterator itr = trackspace.Z_Iterator(lsed.z_img_idx+step*sign); !itr.IsNULL(); itr++){
							float xdelt = slocs[m].x-itr.X(), ydelt = slocs[m].y-itr.Y();
							if(xdelt*xdelt+ydelt*ydelt < 1.5*1.5){
								candidates.push_back(util::point2d(itr.X(), itr.Y()));
								break;
							}
						}
					}

					if(candidates.size() == 1){
						noex = true;
						imatch = imvector->GetMatchSetWithIdx(lsed.z_begin, lsed.z_img_idx, noex);
						if(!imatch){
							imvector->MallocNewMatch();
							imatch = &((*imvector)[imvector->Size()-1]);
							imatch->idx1 = lsed.z_begin; imatch->idx2 = lsed.z_img_idx;
						}
						if(noex){
							imatch->pairs.push_back(std::make_pair(lsed.z_bp, nxloc));
						}
						else{
							imatch->pairs.push_back(std::make_pair(nxloc, lsed.z_bp));
						}

						noex = true;
						imatch = imvector->GetMatchSetWithIdx(lsed.z_img_idx, lsed.z_img_idx+step*sign, noex);
						if(!imatch){
							imvector->MallocNewMatch();
							imatch = &((*imvector)[imvector->Size()-1]);
							imatch->idx1 = lsed.z_img_idx; imatch->idx2 = lsed.z_img_idx+step*sign;
						}
						if(noex){
							imatch->pairs.push_back(std::make_pair(nxloc, candidates[0]));
						}
						else{
							imatch->pairs.push_back(std::make_pair(candidates[0], nxloc));
						}
					}
				}
			}
		}

// 		cvReleaseImage(&slice);
	}
        delete []lsearchinfos;
}


cv::Mat Detector::atrous_cv(cv::Mat ai, cv::Mat k)
{
	/**************
         * @atrous变换处理函数
	     * @输入一个图像ai与数组k，返回处理结果，所用卷积核为数组与自身转置的克罗内克积 ********/

	cv::Mat ker,out;
	ker=cv::Mat::zeros(k.cols,k.cols, CV_32FC1);
	for(int i=0;i<k.cols;i++){
		ker.row(i)=k.at<float>(0,i)*k;
	}
	ai.convertTo(ai,CV_64FC1,1,0);
	ker.convertTo(ker,CV_64FC1,1,0);

	cv::filter2D(ai, out, -1, ker, cv::Point(-1, -1), 0, cv::BORDER_REFLECT);
	out.convertTo(out,CV_32FC1,1,0);

	return out;


}

cv::Mat Detector::wavelet(cv::Mat a, cv::Mat A)
{

	    /********
	          * @小波图像计算函数
	          * @a，A为尺度为i，i+1的近似图像，返回尺度i+1经过阈值化的图像
	          *******/


	 cv::Mat w,temp;
	 float temp1;
	 double t;

	 w=a-A;
	 temp1=util::GetMidValue(w);
	 temp=temp1 * cv::Mat::ones(a.rows, a.cols, CV_32FC1)-w;

//      std::cout<<"tmp:"<<temp.at<float>(0,0)<<std::endl;
// 	 std::cout<<temp.size()<<std::endl;

     for(int i=0;i<temp.rows;i++)
	 {
		   for(int j=0;j<temp.cols;j++)
		   {
			     if(temp.at<float>(i, j)<0) temp.at<float>(i, j) = -temp.at<float>(i, j);
		}
	}

// 	 temp=cv::abs(temp);
// std::cout<<"tmp:"<<temp.at<float>(0,0)<<std::endl;

	 temp1=util::GetMidValue(temp);
	 t=temp1*3/0.67;
	//cout<<"t is "<<t<<endl;
	//if(outpro == 1) cv::imwrite("beforehaidvalsmall.jpg",w*255);
	 util::hardvalsmall(w,t);
	//if(outpro == 1) cv::imwrite("afterhaidvalsmall.jpg",w*255);

	 return w;


}








void Detector::waveletprocess2(cv::Mat img, cv::Mat& wave_img, int J)
{

	/**************
     * @小波变换处理函数
     * @输入一个图像，与小波变换深度，返回小波相关图像
     * sigma为第一幅图的方差,同时传回
     ********/

     assert(J>1);//小波变换深度小于一时终止程序

     cv::Mat inimg,a,k,temp,temp1,temp2,temp3,temp4,aaa,kk;
     double mean1;
     int m, n;

     //atrous核初始化
     k=(cv::Mat_<float>(1,5)<<1.0 / 16, 1.0 / 4, 3.0 / 8, 1.0 / 4, 1.0 / 16);


     //Ai初始化，用于记录所有近似图像
     //int size_Ai[3]={J,img.rows,img.cols};
     //Mat Ai(3,size_Ai,CV_32FC1,cv::Scalar(0));

     //W初始化，用以记录所有小波图像
     //int size_Wi[3]={J-1,img.rows,img.cols};
     //Mat Wi(3,size_Wi,CV_32FC1,cv::Scalar(0));
//           std::cout<<"inimg"<<img<<std::endl;
     img.convertTo(img, CV_8UC1, 255, 0);//转为uchar
     inimg=util::Img_in(img);//图像取反
// 	 util::Reversal(img);
//      	std::cout<<"inimg"<<img<<std::endl;

	 img.convertTo(img, CV_8UC1, 255, 0);                      //转为uchar
     inimg.convertTo(inimg, CV_8UC1, 1, 0);
     a=inimg.clone();
     a.convertTo(a,CV_32FC1,1,0);
     m=a.rows;
     n=a.cols;
//      //Ai=insertmat_dim2to3(Ai,a,0);//记录取反后的原始图像
//
     //计算近似图像组
     temp=a.clone();
     for(int i=1;i<J;i++){
          if(i == J-1) temp1=temp.clone();

		  temp2=atrous_cv(temp, k);
		  temp=temp2;
		  kk=cv::Mat::zeros(1, 2*k.cols-1, CV_32FC1);
		  for(int j=0;j<k.cols;j++)
		  {
				kk.at<float>(0, 2*j)=k.at<float>(0, j);
		  }
//           kk=fullZero2(k,k.cols);
          k=kk;
     }

     wave_img = wavelet(temp1, temp);
//      wave_img = temp2



}


void Detector::removeSmall(cv::Mat img, double threshold)
{
	double ae;
    std::vector<std::vector<cv::Point>> contours;
    std::vector<cv::Vec4i> hierarchy;

    img.convertTo(img,CV_8UC1, 1, 0);
    img = img > 200;                                        //大于200的置为255
    //寻找轮廓
    findContours(img, contours, hierarchy, cv::RETR_TREE, cv::CHAIN_APPROX_SIMPLE, cv::Point(0,0));
    for(int i=0;i<contours.size();i++){
        ae = contourArea(contours[i]);
        if(ae<threshold){
            cv::drawContours(img, contours, i, cv::Scalar(0), -1);
        }
    }
}

double Detector::roundness(cv::Mat label_img)
{
	//计算圆度
    cv::Mat img;
    double a,b;
    std::vector<std::vector<cv::Point>> contours;
    std::vector<cv::Vec4i> hierarchy;

    label_img.convertTo(img,CV_8UC1,1,0);
    //寻找轮廓
    findContours(img, contours, hierarchy, cv::RETR_TREE, cv::CHAIN_APPROX_SIMPLE, cv::Point(0,0));

    a=cv::contourArea(contours[0])*4*M_PI;
    b = pow(cv::arcLength(contours[0],true), 2);
    if (b == 0) return 0;
    return a/b;
}

void Detector::template_make(cv::Mat img, int scale, cv::Mat& template1, cv::Mat& wavelet_img, int& d_end)
{
	// 模板点制备函数
    int j, m, n, k, margin, J, img_m, img_n, num_labels, max_index_hist, temp_index, r_end, num_template;
    int *fid_index;
    float threshold_pixel, threshold_shape, threshold_remove;
    double round, d_mean, minValue_d, maxValue_d;
    cv::Mat img1, img_template, img_2value, labels, stats, centroids, center_int, sub_img, sub_wave, mean, dev, d, hist_d, arr, add_template, temp, wave_ori, wave_img;
    cv::Point minIdx_d, maxIdx_d;

    margin = 2;
    J = scale + 1;
    img_m = img.rows;
    img_n = img.cols;

	threshold_pixel = 0.5;
    threshold_shape = 0.75;
    if (img_m < 2000 || img_n < 2000) threshold_remove = 4;
    else threshold_remove = 8;

	int tmp = img.rows;
	waveletprocess2(img, wave_img, J);
// 	std::cout << wave_img << std::endl;
/*
for(int i=0;i<wave_img.rows;i++)
	 {
		   for(int j=0;j<wave_img.cols;j++)
		   {
			     if(wave_img.at<float>(i, j)<0) std::cout << wave_img.at<float>(i, j) << std::endl;
		}
	}*/


	wave_ori = util::hardval2(wave_img, 2);
	img_2value = wave_ori.clone();
	removeSmall(img_2value, threshold_remove);

	num_labels = cv::connectedComponentsWithStats(img_2value, labels, stats, centroids); // 得到每个连通分支的数据
// 	std::cout<<centroids<<std::endl;
    // stats为左上角坐标，水平/竖直长度，像素个数,int。
    // centroids:double,labels:int
    centroids.convertTo(center_int, CV_32SC1, 1, 0);
    j = 0;
    fid_index = new int[num_labels];
//     std::cout << "img" << img << std::endl;
	for (int i = 1; i < num_labels; i++)
    {
//         std::cout<<stats.at<int>(i,1)<<" "<<stats.at<int>(i,1)+stats.at<int>(i,3)+1<<" "<<stats.at<int>(i,0)<<" "<<stats.at<int>(i,0)+stats.at<int>(i,2)+1<<std::endl;//测试数组越界用

        m = stats.at<int>(i, 1) + stats.at<int>(i, 3) + 1;
        if (m > img.rows) m = tmp;
        sub_img = img(cv::Range(stats.at<int>(i, 1), m), cv::Range(stats.at<int>(i, 0), stats.at<int>(i, 0) + stats.at<int>(i, 2)));
        sub_wave = img_2value(cv::Range(stats.at<int>(i, 1), m), cv::Range(stats.at<int>(i, 0), stats.at<int>(i, 0) + stats.at<int>(i, 2)));

        cv::meanStdDev(sub_img, mean, dev);
//         std::cout << sub_img.at<float>(center_int.at<int>(i, 1) - stats.at<int>(i, 1), center_int.at<int>(i, 0) - stats.at<int>(i, 0)) << "    " << mean.at<double>(0, 0) << "    " << dev.at<double>(0, 0)<< std::endl;
        // remove1: remove the sub_cubic by pixel，判断sub_img中点，均值mean，方差dev是否会合适
        if (sub_img.at<float>(center_int.at<int>(i, 1) - stats.at<int>(i, 1), center_int.at<int>(i, 0) - stats.at<int>(i, 0)) > threshold_pixel || mean.at<double>(0, 0) > threshold_pixel || dev.at<double>(0, 0) < 0.1)
		{
// 	  std::cout << "remove1"<< std::endl;
			continue;
		}

        // remove2: remove the roundness too small,第一个判断是把十分不合理的从小波图去除，第二个判断更强，表示可以用作模板的patches的要求更高
        round = roundness(sub_wave);
        if (round < threshold_shape || threshold_shape == 0.85)
		{
// 	  std::cout << "remove2"<< std::endl;
			continue;
		}

        fid_index[j] = i;
        j++;
    }
    if (j == 0)
    {
        std::cout << "The wavelet detailed coefficients do not get proper information" << std::endl
             << "Please select the scale again." << std::endl;
    }
    d = cv::Mat::zeros(j, 1, CV_32FC1);
    d_mean = 0;

	for (int i = 0; i < j; i++)
    {
        d.at<float>(i, 0) = stats.at<int>(fid_index[i], 2);
        if (stats.at<int>(fid_index[i], 2) < stats.at<int>(fid_index[i], 3))
            d.at<float>(i, 0) = stats.at<int>(fid_index[i], 3);
        d_mean += d.at<float>(i, 0);
    }
    d_mean = (d_mean + 0.0) / (j + 0.0);
    removeSmall(img_2value, d_mean * 0.8);
    cv::minMaxLoc(d, &minValue_d, &maxValue_d, &minIdx_d, &maxIdx_d); // 最值函数
    float sranges[] = {(float)(minValue_d), (float)(maxValue_d + 0.001)};
    const float *ranges[] = {sranges};
    int histSize[] = {10};
    int channels[] = {0};

    // d的频率直方图,arr与hist均为float型
    cv::calcHist(&d, 1, {0}, cv::Mat(), hist_d, 1, histSize, ranges, true, false);

    arr = cv::Mat::zeros(11, 1, CV_32FC1);
    arr.at<float>(0, 0) = (float)minValue_d;
    arr.at<float>(10, 0) = (float)maxValue_d;
    for (int i = 1; i < 10; i++)
    {
        arr.at<float>(i, 0) = (float)(minValue_d + i * (maxValue_d - minValue_d) / 10);
    }

    max_index_hist = 0;
    for (int i = 1; i < 10; i++)
    {
        if (hist_d.at<float>(i, 0) > hist_d.at<float>(max_index_hist, 0))
            max_index_hist = i;
    }
    temp_index = (int)((arr.at<float>(max_index_hist, 0) + arr.at<float>(max_index_hist + 1, 0)) / 2);

    if (temp_index % 2 == 1)
        d_end = temp_index + 2 * margin;
    else
        d_end = temp_index + 1 + 2 * margin;
    r_end = int(((d_end + 0.0) / 2) + 0.5);
// 	std::cout << "temp:"<< d_end  << std::endl;
    add_template = cv::Mat::zeros(2 * r_end + 1, 2 * r_end + 1, CV_32FC1);
// std::cout << "add_template:"<< j << std::endl;
    for (int i = 0; i < j; i++)
    {
        // 判断是否是边缘点
        k = fid_index[i];
        if (center_int.at<int>(k, 1) - r_end < 0 || center_int.at<int>(k, 0) - r_end < 0 || center_int.at<int>(k, 1) + r_end + 1 > img_m || center_int.at<int>(k, 0) + r_end + 1 > img_n)
		{
// 	  std::cout << "j" << std::endl;
			continue;
		}
        temp = img(cv::Range(center_int.at<int>(k, 1) - r_end, center_int.at<int>(k, 1) + r_end + 1), cv::Range(center_int.at<int>(k, 0) - r_end, center_int.at<int>(k, 0) + r_end + 1));
//         std::cout << "temp:"<< temp << std::endl;
        cv::add(temp, add_template, add_template);
        num_template += 1;
    }
//     std::cout << "temp:"<< center_int.at<int>(k, 0) + r_end + 1 << std::endl;
//     std::cout << "add_template:"<< temp << std::endl;
    add_template.convertTo(add_template, CV_32FC1, 1.0 / num_template, 0);

//     delete (fid_index);

    img_2value.copyTo(wavelet_img);
    add_template.copyTo(template1);

}

void Detector::template_average(cv::Mat& template1)
{
	// This function used to average template to make the center point more center
    int height, width;
    cv::Mat rotate_matrix, rotated_image1, rotated_image2, rotated_image3, template_new;
    cv::Point2f center;

    // dividing height and width by 2 to get the center of the image
    height = template1.rows;
    width = template1.cols;

    // get the center coordinates of the image to create the 2D rotation matrix
    center.x = (width - 1) / 2;
    center.y = (height - 1) / 2;

    // using cv2.getRotationMatrix2D() to get the rotation matrix
    rotate_matrix = getRotationMatrix2D(center, 90, 1);

    // rotate the image using cv2.warpAffine
    warpAffine(template1, rotated_image1, rotate_matrix, cv::Size(width, height));
    warpAffine(rotated_image1, rotated_image2, rotate_matrix, cv::Size(width, height));
    warpAffine(rotated_image2, rotated_image3, rotate_matrix, cv::Size(width, height));

    cv::add(template1, rotated_image1, template_new);
    cv::add(template_new, rotated_image2, template_new);
    cv::add(template_new, rotated_image3, template_new);
    template_new.convertTo(template_new, CV_32FC1, 1.0 / 4, 0);

    template1 = template_new.clone();
//     return template_new;

}

cv::Mat Detector::refine_fid_by_gaussian_distribution_markerauto_no_wave(cv::Mat candidate)
{
	 // 利用ncc*pixel的分布筛选candidate
    // :param candidate:坐标集合，[x,y1,ncc,pixel,none,index]
    // :return:new_fid每一维数据分别是(x,y1,ncc,avg_pixel,ncc*avg_pixel,index)
    int i, num, fid_index;
    float avg, stdev, thre;
    cv::Mat new_score, meantemp, stdtemp, new_fid;

    // 构造ncc*pixel
    num = candidate.rows;

    new_score = cv::Mat::zeros(num, 1, CV_32FC1);
    for (i = 0; i < num; i++)
        new_score.at<float>(i, 0) = candidate.at<float>(i, 2) * candidate.at<float>(i, 3);

    avg = 0.0;
    stdev = 0.0;
    for (i = 0; i < num; i++)
    {
        avg += new_score.at<float>(i, 0);
        stdev += new_score.at<float>(i, 0) * new_score.at<float>(i, 0);
    }
    avg /= num;
    stdev = sqrt(stdev / num - avg * avg);
    thre = avg + 3 * stdev;

    new_fid = cv::Mat::zeros(num, 6, CV_32FC1);
    fid_index = 0;
    for (i = 0; i < num; i++)
    {
        if (new_score.at<float>(i, 0) > thre)
        {
            new_fid.at<float>(fid_index, 0) = candidate.at<float>(i, 0);
            new_fid.at<float>(fid_index, 1) = candidate.at<float>(i, 1);
            new_fid.at<float>(fid_index, 2) = candidate.at<float>(i, 2);
            new_fid.at<float>(fid_index, 3) = candidate.at<float>(i, 3);
            new_fid.at<float>(fid_index, 4) = new_score.at<float>(i, 0);
            new_fid.at<float>(fid_index, 5) = fid_index;
            fid_index++;
        }
    }
    new_fid = new_fid(cv::Range(0, fid_index), cv::Range(0, 6));

    return new_fid;
}

cv::Mat Detector::refine_fid_by_gaussian_distribution_markerauto_wave(cv::Mat candidate)
{
	// 利用ncc*pixel的分布筛选candidate
    // :param candidate:坐标集合，[x,y1,ncc,pixel,none,index]
    // :return:new_fid每一维数据分别是(x,y1,ncc,avg_pixel,ncc*avg_pixel,index)
    int i, num, fid_index;
    float avg, stdev, thre;
    cv::Mat new_score, meantemp, stdtemp, new_fid;

    // 构造ncc*pixel
    num = candidate.rows;
    new_score = cv::Mat::zeros(num, 1, CV_32FC1);
    for (i = 0; i < num; i++)
        new_score.at<float>(i, 0) = candidate.at<float>(i, 2) * candidate.at<float>(i, 3);

    avg = 0.0;
    stdev = 0.0;
    for (i = 0; i < num; i++)
    {
        avg += new_score.at<float>(i, 0);
        stdev += new_score.at<float>(i, 0) * new_score.at<float>(i, 0);
    }
    avg /= num;
    stdev = sqrt(stdev / num - avg * avg);
    thre = avg - 0.5 * stdev;

    new_fid = cv::Mat::zeros(num, 6, CV_32FC1);
    fid_index = 0;
    for (i = 0; i < num; i++)
    {
        if (new_score.at<float>(i, 0) > thre)
        {
            new_fid.at<float>(fid_index, 0) = candidate.at<float>(i, 0);
            new_fid.at<float>(fid_index, 1) = candidate.at<float>(i, 1);
            new_fid.at<float>(fid_index, 2) = candidate.at<float>(i, 2);
            new_fid.at<float>(fid_index, 3) = candidate.at<float>(i, 3);
            new_fid.at<float>(fid_index, 4) = new_score.at<float>(i, 0);
            new_fid.at<float>(fid_index, 5) = fid_index;
            fid_index++;
        }
    }
    new_fid = new_fid(cv::Range(0, fid_index), cv::Range(0, 6));

    return new_fid;
}



void Detector::markerauto_work_flow(cv::Mat& fidstack, cv::Mat img_ori, cv::Mat template_ori, int radius_int, int diameter_int, cv::Mat wave_img)
{
	int idiameter, no_index, wave_index, corr_m, corr_n, index_remove,candidate_index_location,dist_thr,new_fid_index;
    double img_mean, img_std_dev, corr_mean, corr_std_dev, img_threshold, corr_threshold;
    cv::Point xy;
    cv::Mat img, img_draw, fid, template1, corr, meantemp, stdtemp, img_draw_temp, candidate_no_wave, candidate_wave, remove_point, fid_no_wave, fid_wave, new_fid_temp, new_fid_no_wave, new_fid,L;
    Node_kd node;
    clock_t start_time, end_time;
    cv::cvtColor(img_ori, img_draw, cv::COLOR_GRAY2RGB);
    // pre-processing
    // 归一化与取反操作
    img = util::normalize1(img_ori);
    img = util::Img_in2(img);
// 	std::cout<<img<<std::endl;
    cv::Mat template_o = template_ori.clone();
    template1 = util::ToOne(template_o);
    template1 = util::Img_in2(template1);
//     std::cout << template_ori << std::endl;
// std::cout << "peak_n:" << template_ori << std::endl;
    // 模板匹配并归一化
    cv::matchTemplate(img, template1, corr, 3);
    corr = util::ToOne(corr);
// 	std::cout<<corr<<std::endl;
    // cv::imwrite("corr.jpg",corr);
    cv::meanStdDev(img, meantemp, stdtemp);
    img_mean = meantemp.at<double>(0, 0);
    img_std_dev = stdtemp.at<double>(0, 0);
    cv::meanStdDev(corr, meantemp, stdtemp);
    corr_mean = meantemp.at<double>(0, 0);
    corr_std_dev = stdtemp.at<double>(0, 0);
    img_threshold = img_mean;

#define CUTOFF_THRE       2//2.5

    corr_threshold = corr_mean + CUTOFF_THRE * corr_std_dev;
// std::cout << corr_mean << std::endl;
//     int idiameter1 = int(2 * radius_int + 1);
	idiameter = int(2 * radius_int + 1);
// 	idiameter=detector.fid_diameter;
//     std::cout<<"r is "<<radius_int<<std::endl;
// std::cout << "idiameter1:" << idiameter1<<      "idiameter:" << idiameter << std::endl;
    no_index = 0;
    wave_index = 0;
    index_remove = 0;
    img_draw_temp = img_draw.clone();
    corr_m = corr.rows;
    corr_n = corr.cols;
    // cout<<corr_m<<endl<<corr_n<<endl<<idiameter<<endl<<(int)((corr_m/idiameter)*(corr_n/idiameter)+10)<<endl;
    candidate_wave = cv::Mat::zeros((int)(((corr_m + 0.0) / idiameter) * ((corr_n + 0.0) / idiameter) + 10), 6, CV_32FC1);
// 	std::cout << corr_m<<"  "<<corr_n<<"  "<<idiameter << std::endl;
// 	std::cout << corr_threshold<<"  "<<img_threshold << std::endl;
    candidate_no_wave = cv::Mat::zeros((int)(((corr_m + 0.0) / idiameter) * ((corr_n + 0.0) / idiameter) + 10), 6, CV_32FC1);
    remove_point = cv::Mat::zeros((int)(((corr_m + 0.0) / idiameter) * ((corr_n + 0.0) / idiameter) + 10), 2, CV_32FC1);
	int tmp=0,tmp2=0;
    for (int i = 0; i < corr_m; i = i + idiameter)
    {
// 		std::cout<<"ok"<<std::endl;
		tmp++;
        for (int j = 0; j < corr_n; j = j + idiameter)
        {
			int peak_m, peak_n;
// 			FindMaxPeak(corr,  i,  j,  idiameter, &peak_m, &peak_n);
            util::find_local_peak(corr, i, j, idiameter, idiameter, peak_m, peak_n);
// 			std::cout<<peak_m<<"      "<<peak_n<<std::endl;
            // 筛选一，ncc达标，峰值达标，是小波图像中的点
//             if (corr.at<float>(peak_m, peak_n) > corr_threshold && get_ave_pixel(img, peak_n, peak_m, radius_int) > img_threshold)
//             {
// 	std::cout << "peak_n:" << corr.at<float>(peak_m, peak_n) << std::endl;
			if (corr.at<float>(peak_m, peak_n) > corr_threshold && util::get_ave_pixel(img, peak_n, peak_m, radius_int) > img_threshold)
            {
				tmp2++;
// 								            std::cout << wave_img.size()<<"     "<<peak_m + 2 * radius_int + 1<<"            "
// 											<<peak_n + 2 * radius_int + 1<< std::endl;
                if (util::have255(wave_img(cv::Range(peak_m, peak_m + 2 * radius_int + 1), cv::Range(peak_n, peak_n + 2 * radius_int + 1))))
                {
                    candidate_wave.at<float>(wave_index, 0) = peak_n;
                    candidate_wave.at<float>(wave_index, 1) = peak_m;
                    candidate_wave.at<float>(wave_index, 2) = corr.at<float>(peak_m, peak_n);
                    candidate_wave.at<float>(wave_index, 3) = util::get_ave_pixel(img, peak_n, peak_m, radius_int);
                    candidate_wave.at<float>(wave_index, 4) = 1;
                    candidate_wave.at<float>(wave_index, 5) = wave_index;
                    xy.x = peak_n + radius_int;
                    xy.y = peak_m + radius_int;
                    cv::circle(img_draw_temp, xy, radius_int, cv::Scalar(0, 255, 0));
                    wave_index++;
                }
                else
                {
                    candidate_no_wave.at<float>(no_index, 0) = peak_n;
                    candidate_no_wave.at<float>(no_index, 1) = peak_m;
                    candidate_no_wave.at<float>(no_index, 2) = corr.at<float>(peak_m, peak_n);
                    candidate_no_wave.at<float>(no_index, 3) = util::get_ave_pixel(img, peak_n, peak_m, radius_int);
                    candidate_no_wave.at<float>(no_index, 4) = 1;
                    candidate_no_wave.at<float>(no_index, 5) = no_index;
                    xy.x = peak_n + radius_int;
                    xy.y = peak_m + radius_int;
                    cv::circle(img_draw_temp, xy, radius_int, cv::Scalar(255, 0, 0));
                    no_index++;
                }
            }
            else
            {
// 				std::cout << "peak_n:" << peak_n<< "peak_m:" << peak_m << std::endl;
                remove_point.at<float>(index_remove, 0) = (float)peak_n;
                remove_point.at<float>(index_remove, 1) = (float)peak_m;
                index_remove++;
            }
        }
    }
//     std::cout << tmp2<< std::endl;
//     std::cout << "candidate_wave" << candidate_wave.size()<< "candidate_no_wave" << candidate_no_wave.size() << std::endl;
    candidate_wave = candidate_wave(cv::Range(0, wave_index), cv::Range(0, 6));
    candidate_no_wave = candidate_no_wave(cv::Range(0, no_index), cv::Range(0, 6));
// 	std::cout <<remove_point.size()<<"    "<<index_remove<<   std::endl;
    remove_point = remove_point(cv::Range(0, index_remove), cv::Range(0, 2));
// 	std::cout <<wave_index<<"  "<<no_index<<"    "<<index_remove<<   std::endl;
    fid_no_wave = candidate_no_wave;
    fid_wave = candidate_wave;
// 	std::cout << "candidate_wave" << candidate_wave.size()<< "candidate_no_wave" << candidate_no_wave.size() << std::endl;
//     	std::cout << "fid_no_wave:" << fid_no_wave.rows<< "fid_wave:" << fid_wave.rows << std::endl;
//     std::cout << "new_fid_temp1:" << fid_wave.rows << std::endl;
    new_fid_no_wave = refine_fid_by_gaussian_distribution_markerauto_no_wave(fid_no_wave);
	new_fid_temp = cv::Mat::zeros(new_fid_no_wave.rows + fid_wave.rows, 6, CV_32FC1);
// 	std::cout << "new_fid_no_wave:" << new_fid_no_wave.rows<< "fid_wave:" << fid_wave.rows << "new_fid_temp:" << new_fid_temp.rows << std::endl;
	for (int i = 0; i < new_fid_no_wave.rows; i++){
        for(int j=0;j<6;j++){
            new_fid_temp.at<float>(i,j) = new_fid_no_wave.at<float>(i,j);
        }
    }
//     std::cout << "new_fid_temp1:" << new_fid_temp.rows << std::endl;
    for (int i = 0; i < fid_wave.rows; i++){
        for(int j=0;j<6;j++){
            new_fid_temp.at<float>(i + new_fid_no_wave.rows,j) = fid_wave.at<float>(i,j);
        }
    }
// 	std::cout << "new_fid_temp2:" << new_fid_temp.rows << std::endl;
	new_fid = refine_fid_by_gaussian_distribution_markerauto_wave(new_fid_temp);
	fid = new_fid.clone();
// 	std::cout << "new_fid" << fid.rows << std::endl;

	//Remove repeated kd tree
    candidate_index_location = 5;// fid 的第几个分量为index
    img_draw_temp=img_draw.clone();
    dist_thr = diameter_int;
    new_fid = cv::Mat::zeros(fid.rows,6,CV_32FC1);
    new_fid_index = 0;
    node.construct(2, fid.clone(), &node, 0);
    for(int i = 0; i < fid.rows; i++){
        if(fid.at<float>(i, 4) < 0) continue;
        L = cv::Mat::zeros(0, 0, CV_32FC1);
        L = node.search(&node, fid.row(i), L, 5);
        for(int j = 0; j < L.rows; j++){
            if(node.distance(fid.row(i), L.row(j))<dist_thr)
            fid.at<float>((int)L.at<float>(j, candidate_index_location), 4) = -1;
        }
        new_fid.row(new_fid_index) = fid.row(i).clone() + 0.0;
        new_fid_index++;
        xy.x = (int)fid.at<float>(i,0) + radius_int;
        xy.y = (int)fid.at<float>(i,1) + radius_int;
        cv::circle(img_draw_temp, xy, radius_int, cv::Scalar(0, 0, 255));
        node.clear_flag(&node);
    }
    new_fid=new_fid(cv::Range(0, new_fid_index), cv::Range(0, new_fid.cols));
    fid = new_fid;

    fidstack = new_fid.clone();

// 	std::cout << "fid" << fid << std::endl;
//     return fid;


}


cv::Mat Detector::Gauss_A(cv::Mat loc)
{
    int i;
    cv::Mat A(4,4,CV_32FC1,cv::Scalar(0));
    for(i=0;i<loc.rows;i++){
        A.at<float>(0,0)+=(loc.at<float>(i,0)*loc.at<float>(i,0)+loc.at<float>(i,1)*loc.at<float>(i,1))*(loc.at<float>(i,0)*loc.at<float>(i,0)+loc.at<float>(i,1)*loc.at<float>(i,1));

        A.at<float>(0,1)+=loc.at<float>(i,0)*(loc.at<float>(i,0)*loc.at<float>(i,0)+loc.at<float>(i,1)*loc.at<float>(i,1));

        A.at<float>(0,2)+=loc.at<float>(i,1)*(loc.at<float>(i,0)*loc.at<float>(i,0)+loc.at<float>(i,1)*loc.at<float>(i,1));

        A.at<float>(0,3)+=(loc.at<float>(i,0)*loc.at<float>(i,0)+loc.at<float>(i,1)*loc.at<float>(i,1));


        A.at<float>(1,1)+=loc.at<float>(i,0)*loc.at<float>(i,0);

        A.at<float>(1,2)+=loc.at<float>(i,0)*loc.at<float>(i,1);

        A.at<float>(1,3)+=loc.at<float>(i,0);

        A.at<float>(2,2)+=loc.at<float>(i,1)*loc.at<float>(i,1);

        A.at<float>(2,3)+=loc.at<float>(i,1);

        A.at<float>(3,3)++;
    }

    A.at<float>(1,0)=A.at<float>(0,1);
    A.at<float>(2,0)=A.at<float>(0,2);
    A.at<float>(2,1)=A.at<float>(1,2);
    A.at<float>(3,0)=A.at<float>(0,3);
    A.at<float>(3,1)=A.at<float>(1,3);
    A.at<float>(3,2)=A.at<float>(2,3);

    return A;

}


cv::Mat Detector::Gauss_b(cv::Mat info)
{
    int i;
    cv::Mat b(4,1,CV_32FC1,cv::Scalar(0));

    for(i=0;i<info.rows;i++){
        b.at<float>(0,0)+=(info.at<float>(i,0)*info.at<float>(i,0)+info.at<float>(i,1)*info.at<float>(i,1))*log(info.at<float>(i,2));

        b.at<float>(0,1)+=info.at<float>(i,0)*log(info.at<float>(i,2));

        b.at<float>(0,2)+=info.at<float>(i,1)*log(info.at<float>(i,2));

        b.at<float>(0,3)+=log(info.at<float>(i,2));
    }
    return b;
}



cv::Mat Detector::compute_center_Gauss(cv::Mat pic)
{
	int m, n, i, j;
    float sigma, x0, y0, para_A;
    cv::Mat paraDir,loc,info,X,Y,A,b,x;

    m=pic.rows;
    n=pic.cols;

    X.create(m*n,1,CV_32FC1);
    Y.create(m*n,1,CV_32FC1);
    loc.create(m*n,2,CV_8UC1);
    info.create(m*n,3,CV_32FC1);
    paraDir.create(4,1,CV_32FC1);

    for(i=0;i<m;i++){
        for(j=0;j<n;j++){
            info.at<float>(i*n+j,0)=(i+0.0)/(m+0.0);
            info.at<float>(i*n+j,1)=(j+0.0)/(n+0.0);
            info.at<float>(i*n+j,2)=pic.at<float>(i,j);
        }
    }
    A = Gauss_A(info);
    b = Gauss_b(info);
    cv::solve(A,b,x,cv::DECOMP_LU);//因为只有4*4，所以用LU

    sigma=sqrt(0.5/abs(x.at<float>(0,0))*m*n);
    x0=-0.5*x.at<float>(0,1)/x.at<float>(0,0)*m;
    y0=-0.5*x.at<float>(0,2)/x.at<float>(0,0)*n;
    para_A=exp(x.at<float>(0,3)-(x.at<float>(0,2)*x.at<float>(0,2)+x.at<float>(0,1)*x.at<float>(0,1))/(4.0*x.at<float>(0,1)));

    paraDir.at<float>(0,0)=sigma;
    paraDir.at<float>(0,1)=x0;
    paraDir.at<float>(0,2)=y0;
    paraDir.at<float>(0,3)=para_A;

    return paraDir;
}

float Detector::compute_gauss_error(cv::Mat pic, cv::Mat dir)
{
	// 判断gaussian拟合结果的好坏
    // :param picture: 被拟合子图
    // :param x0: 拟合参数x0
    // :param y0: 拟合参数y0
    // :param sigma: 拟合参数sigma
    // :param para_A: 拟合参数A
    // :return:均方根误差
    int i,j;
    cv::Mat f;
    float rmse = 0;
    f = util::gaussian(dir.at<float>(0,3), dir.at<float>(0,1), dir.at<float>(0,2), dir.at<float>(0,0), pic.rows, pic.cols);
//     std::cout << "f" << f << std::endl;
    f=f-pic;
    for(i=0;i<pic.rows;i++){
        for(j=0;j<pic.cols;j++){
                rmse+=f.at<float>(i,j);
        }
    }
//     std::cout << "rmse" << rmse << std::endl;
    rmse = sqrt(abs(rmse)/(pic.rows*pic.cols+0.0));
    return rmse;
}


void Detector::location_fid(std::vector<util::point2d>& fids, cv::Mat img, cv::Mat location, float width, cv::Mat& refined_xy, cv::Mat& score_xy)
{
	/*****************
    标记点定位阶段
    :param img: 定位的原图片
    :param location: 识别阶段的坐标点,都为左上角
    :param width: 进行圆心refine的子图的大小
    :return:refined_xy, 定位后fids的坐标
    *****************/
    cv::Mat img_inv,sub_img,gauss_out,dir;
    img_inv = img.clone();
    location.convertTo(location,CV_32FC1,1,0);
//     img.convertTo(img_inv,CV_32FC1,-1,255);
	util::Reversal(img_inv);

    refined_xy=cv::Mat::zeros(location.rows,2,CV_32FC1);
    score_xy=cv::Mat::zeros(location.rows,1,CV_32FC1);
    // fid元素解释： fid=(x,y1,ncc,avg_pixel,ncc*avg_pixel,index)
    for(int i=0; i<location.rows; i++){
		util::point2d point;
        //cout<<location.at<int>(i,1)+width<<" "<<location.at<int>(i,0)+width<<endl;
        sub_img=img_inv(cv::Range(location.at<float>(i,1), location.at<float>(i,1)+width),cv::Range(location.at<float>(i,0),location.at<float>(i,0)+width));
        sub_img=util::change_aTob(sub_img, 0, 2);
//
        gauss_out=compute_center_Gauss(sub_img);
//         std::cout << "gauss_out" << gauss_out << std::endl;
// std::cout << "gauss_out" << gauss_out.size() << std::endl;
//         //cout<<gauss_out<<endl;//sigma,x0,y0,para_A
        refined_xy.at<float>(i,0) = gauss_out.at<float>(0,1) + location.at<float>(i,0);
        refined_xy.at<float>(i,1)=gauss_out.at<float>(0,2)+location.at<float>(i,1);

        point.x = gauss_out.at<float>(0,1) + location.at<float>(i,0);
        point.y = gauss_out.at<float>(0,2)+location.at<float>(i,1);
        fids.push_back(point);
//         score_xy.at<float>(i,0)=compute_gauss_error(sub_img, gauss_out);
    }

}

void Detector::Wavelet(cv::Mat img, int scale, cv::Mat& wavelet_img)
{
	int img_m, img_n;
	img_m = img.rows;
    img_n = img.cols;

    cv::Mat wave_img, wave_ori;

    int threshold_remove;
    if (img_m < 2000 || img_n < 2000) {
		threshold_remove = 4;
	}
	else threshold_remove = 8;

	int J = scale + 1;
	waveletprocess2(img, wave_img, J);
    wave_ori = util::hardval2(wave_img, 2);
	cv::Mat img_2value = wave_ori.clone();
	removeSmall(img_2value, threshold_remove);


	cv::Mat labels, stats, centroids, center_int;
	int num_labels;
	num_labels = cv::connectedComponentsWithStats(img_2value, labels, stats, centroids); // 得到每个连通分支的数据
    centroids.convertTo(center_int, CV_32SC1, 1, 0);
    int j = 0;
    int *fid_index = new int[num_labels];

    int tmp = img.rows, m;
    cv::Mat sub_img, sub_wave;
    cv::Mat mean, dev;
    float threshold_pixel = 0.5;
    float threshold_shape = 0.75;
    double round;

    for (int i = 1; i < num_labels; i++)
    {

        m = stats.at<int>(i, 1) + stats.at<int>(i, 3) + 1;
        if (m > img.rows) m = tmp;
        sub_img = img(cv::Range(stats.at<int>(i, 1), m), cv::Range(stats.at<int>(i, 0), stats.at<int>(i, 0) + stats.at<int>(i, 2)));
        sub_wave = img_2value(cv::Range(stats.at<int>(i, 1), m), cv::Range(stats.at<int>(i, 0), stats.at<int>(i, 0) + stats.at<int>(i, 2)));

        cv::meanStdDev(sub_img, mean, dev);
        if (sub_img.at<float>(center_int.at<int>(i, 1) - stats.at<int>(i, 1), center_int.at<int>(i, 0) - stats.at<int>(i, 0)) > threshold_pixel || mean.at<double>(0, 0) > threshold_pixel || dev.at<double>(0, 0) < 0.1)
		{
			continue;
		}

        round = roundness(sub_wave);
        if (round < threshold_shape || threshold_shape == 0.85)
		{
			continue;
		}

        fid_index[j] = i;
        j++;
    }

	cv::Mat d = cv::Mat::zeros(j, 1, CV_32FC1);
    double d_mean = 0;

	for (int i = 0; i < j; i++)
    {
        d.at<float>(i, 0) = stats.at<int>(fid_index[i], 2);
        if (stats.at<int>(fid_index[i], 2) < stats.at<int>(fid_index[i], 3))
            d.at<float>(i, 0) = stats.at<int>(fid_index[i], 3);
        d_mean += d.at<float>(i, 0);
    }
    d_mean = (d_mean + 0.0) / (j + 0.0);
    removeSmall(img_2value, d_mean * 0.8);

    wavelet_img = img_2value.clone();

}

void Detector::RefinePositions(cv::Mat& img, const cv::Mat& tmplt, std::vector<util::point2d>& fids, float diameter)
{
	int magnif = 1;

	if(diameter <= 16){
		magnif = 16;
	}
	else if(diameter <= 32){
		magnif = 8;
	}
	else if(diameter <= 128){
		magnif = 4;
	}
	else if(diameter <= 150){
		magnif = 4;
	}
	else{
		magnif = 2;
	}
	cv::Mat tmp = tmplt.clone();
	util::Reversal(tmp);
// std::cout << tmplt << std::endl;
	RefinePositionsForSubTomo(img, tmp, fids, diameter, magnif);
}

void Detector::Process_new(cv::Mat& slice, cv::Mat& img1, cv::Mat& img, std::vector<util::point2d>& fids, bool use_avgtmplt, float ratio)
{
	EX_TIME_BEGIN("\nCompute Positions of Fiducial Markers (X Y)")
	cv::Mat tmplt;
	cv::Mat wave_img, fids_tmp;
	float diameter = fid_diameter;
// 	float diameter = 35.91;
	float radius = fid_radius;
	if(!use_avgtmplt){
		CreateTemplate(&tmplt, fid_diameter);
	}
	else{
        tmplt=fid_avgtmplt.clone();
	}
// 	std::cout<<img1<<std::endl;
// 	std::cout << diameter<< "   "<<radius <<   std::endl;
// 	cv::Mat img_no;
	cv::Mat img_no = img1.clone();
// 	std::cout << diameter<< "   "<<radius <<   std::endl;
// 	std::cout<<img_no<<std::endl;
	detector.Wavelet(img1, detector.scale, wave_img);
// 	std::cout<<tmplt<<std::endl;
// 		std::cout<<img<<std::endl;
	detector.markerauto_work_flow(fids_tmp, img, tmplt, radius, diameter, wave_img);
// 	std::cout<<fids_tmp.size()<<std::endl;
// std::cout << ori_fid.size() << std::endl;
	cv::Mat ori_fid = fids_tmp * ratio;
// 	float temp_ori_m=11.34;
	int temp_ori_m=fid_oritmplt.rows;
	cv::Mat score_error, xy_refine;
// 	std::cout << ori_fid.size() << std::endl;
// 	std::cout<<"ori_fid:"<<ori_fid<<std::endl;
// 		std::cout<<"temp_ori_m:"<<temp_ori_m<<std::endl;
// 		std::cout<<"xy_refine:"<<xy_refine<<std::endl;
	detector.location_fid(fids, slice, ori_fid, temp_ori_m, xy_refine, score_error);


// 	std::vector<util::point2d> candidate;
// 	for (int i = 0;i < ori_fid.rows;i++)
// 	{
// 		util::point2d point;
// 		point.x = ori_fid.at<float>(i,0)+radius;
//         point.y = ori_fid.at<float>(i,1)+radius;
//         candidate.push_back(point);
// 	}
//
// 	util::Reversal(slice);
// 	util::ConvertTo1(slice, false);
// 	util::MedianSmooth(slice);
//
// 	diameter=36;
//
// 	for(int i = 0; i < candidate.size(); i++){
// // 		std::cout << candidate[i].x<<"  "<< candidate[i].y<< std::endl;
// 		util::point2d fid;
// 		if(GetCenter(slice, candidate[i], diameter, fid)){
// 			fids.push_back(fid);
// 		}
// 	}
//
// // 	RefinePositions(slice, fid_oritmplt, fids, diameter);
// 	RefinePositions(slice, tmplt, fids, diameter);

// 	std::cout << fids_tmp << std::endl;


	EX_TIME_END("Totally found %ld fiducial markers", fids.size())
}


void Detector::DetectorMain(util::MrcStack& mrcr, util::FiducialStack* fidsk, float& diameter, const std::vector< int >& m_index, float ratio)
{
    Detector detector;
    fidsk->ReSize(mrcr.Size());
//     fidsk->ReSize(m_index.size());
 	fidsk->SetRatio(ratio);
    fidsk->SetWxH(mrcr.Width()*ratio, mrcr.Height()*ratio);
 	diameter *= ratio;

 	EX_TRACE("MRC Image rescale ratio: %.2f\n", ratio)
    cv::Mat mid_slice = mrcr.GetStackImage(mrcr.Size()/2);
 	util::ScaleImage(mid_slice, ratio);
 	util::Reversal(mid_slice);
 	util::ConvertTo1(mid_slice, false);
 	util::MedianSmooth(mid_slice);
//
 	float refined_diameter;
 	cv::Mat avgtmplt;
//     std::cout<<"ok"<<std::endl;
 	detector.DetermineDiameter(mid_slice, diameter, &refined_diameter, &avgtmplt);
 	detector.SetFidDiameter(refined_diameter);
 	detector.SetSampling(mid_slice, refined_diameter);
 	detector.SetFidAverageTemplate(avgtmplt);
 	diameter = refined_diameter;

 	detector.fid_avgtmplt = avgtmplt.clone();
 	util::ConvertTo1(detector.fid_avgtmplt, false);

 	util::SeriesSaveToFile(avgtmplt, "avgtmplt");

 	util::Reversal(avgtmplt);
 	util::ConvertTo1(avgtmplt, false);
 	util::SaveImage(avgtmplt, "avgtmplt(reversal).pgm");

// 	for(int i = 0; i < m_index.size(); i++){
    for(int i = 0; i < mrcr.Size(); i++){
//         int j=m_index[i];
        EX_TIME_BEGIN("%sProcessing MRC[%d]", _DASH, i)

//         cv::Mat slice = mrcr.GetStackImage(j);
        cv::Mat slice = mrcr.GetStackImage(i);
		util::ScaleImage(slice, ratio);
		util::Reversal(slice);
        util::ConvertTo1(slice, false);
        util::MedianSmooth(slice);
//         util::HistogramStretch(slice);

// 	std::stringstream ss;
// 	ss<<"histo"<<i<<".pgm";
// 	util::SaveImage(slice, ss.str().c_str());

        std::vector<util::point2d>& fids = fidsk->V(i);
        detector.Process(slice, fids, true);

        EX_TIME_END("Processing MRC[%d]", i)
    }

}

void Detector::Test(util::MrcStack& mrcr, const util::FiducialStack& fidsk, float diameter, const char* folder)
{
    EX_TIME_BEGIN("%s\nDetector Testing", _DASH)

	if(diameter < 0){
		diameter = detector.fid_diameter;
	}

    for(int i = 0; i < mrcr.Size(); i++){
        cv::Mat slice = mrcr.GetStackImage(i);

        cv::Mat init = cv::Mat::zeros(cv::Size(slice.size().width*fidsk.Ratio(), slice.size().height*fidsk.Ratio()), CV_32FC1);
        cv::resize(slice, init, cv::Size(slice.size().width*fidsk.Ratio(), slice.size().height*fidsk.Ratio()), 0, 0, cv::INTER_CUBIC);
        slice = init;
        util::ConvertTo1(slice, true);

        const std::vector<util::point2d>& fids = fidsk.V(i);

		diameter = diameter;//*fidsk.Ratio();			!WARNING

        cv::Mat tmp;
        cv::cvtColor(slice,tmp,cv::COLOR_GRAY2BGR);
        DrawFiducialMarkerPositions(slice, diameter, fids);

        if(access(folder,0) == -1){		//create file folder
            mkdir(folder,0777);
        }
        std::ostringstream oss;
        oss <<"/home/xzh/文档/markerauto/mk_detection/build2/bin/fids/"<<"("<<i<<").pgm";
//         oss <<folder<<"/"<<mrcr.Name()<<"("<<i<<").pgm";
        try {
            util::SaveImage(slice, oss.str().c_str());
        } catch(ex::Exception& e){
            EX_TRACE("%s\n", e.Msg())
        }


//         std::ofstream fidfile;
// 		fidfile.open("/home/xzh/文档/markerauto/mk_detection/build2/bin/fids" + "/"+"fiducial_"+io.inputfilename+".txt");
	//     for(i=0;i<fid.rows;i++){
	//         fidfile<<"["<<fid.at<int>(i,0)<<" "<<fid.at<int>(i,1)<<"]"<<endl;
	//         xy.x=fid.at<int>(i,0);
	//         xy.y=fid.at<int>(i,1);
	//         circle(ori_img_draw, xy, int(radius_int*mul_para),cv::Scalar(0, 0, 255));
	//     }




    }
    EX_TIME_END("Detector Testing")
}


void Detector::TransMain(util::MrcStack& mrcs, HSetVector& hsets, util::FiducialStack* fids_high, util::FiducialStack* fidsk, util::FiducialStack* addfisk, util::FiducialStack* pro_fidsk, util::FiducialStack* allfids, float diameter, const std::vector<int>& m_index)
{
    int num=fidsk->Size();
    pro_fidsk->ReSize(num);
    addfisk->ReSize(num);
    allfids->ReSize(num);
    Detector::InitLocalDetector(diameter);
    for(int i=0;i<fidsk->Size();i++)
    {
        std::vector<util::point2d>& fids2 = pro_fidsk->V(i);
//         std::cout<<"time"<<i<<std::endl;
        std::vector<util::point2d>& fids = fidsk->V(i);
        cv::Mat slice = mrcs.GetStackImage(m_index[i]);
		util::Reversal(slice);
		util::ConvertTo1(slice);
		util::MedianSmooth(slice);
        std::vector<util::point2d>& orifids = fids_high->V(i);
//         std::cout<<"size:"<<fids.size()<<std::endl;
        for(int j=0;j<fids.size();j++)
        {
            std::vector<util::point2d> slocs;
//             std::cout<<"x:"<<fids[j].x<<" y:"<<fids[j].y<<std::endl;
            PredictMissFiducial(fids[j], hsets, i, i, slocs);
//             std::cout<<"x:"<<slocs[0].x<<" y:"<<slocs[0].y<<std::endl;
            util::point2d nxloc;
            bool exist=false;
			float max_score = -999;
			for(int m = 0; m < slocs.size(); m++){
				util::point2d p;
				float score = LocalMarkerDetect(slice, slocs[m], &p);
				if(score > max_score){
					max_score = score;
					nxloc = p;
				}
// 				else{
//                     std::cout<<"***************************"<<std::endl;
//                     std::cout<<"j:"<<j<<"m:"<<m<<std::endl;
//                     std::cout<<p.x<<" "<<p.y<<std::endl;
//                     std::cout<<"score:"<<score<<std::endl;
//                     std::cout<<"***************************"<<std::endl;
//
//                 }
			}
// 			std::cout<<max_score<<std::endl;
			int K=10;
			for(int k=0;k<orifids.size();k++){
                float xdelt = orifids[k].x-nxloc.x, ydelt = orifids[k].y-nxloc.y;
// 				std::cout<<xdelt*xdelt+ydelt*ydelt<<std::endl;
                if(xdelt*xdelt+ydelt*ydelt < K*K+K*K){
                    exist=true;
					break;
                }
//                 else{
// 					std::cout<<xdelt*xdelt+ydelt*ydelt<<std::endl;
// 				}
            }
			fids2.push_back(nxloc);
			if(max_score > 0.85 && !exist){								//WARNING threshold
// 				std::cout<<"存入点："<<nxloc.x<<" "<<nxloc.y<<std::endl;
                addfisk->V(i).push_back(nxloc);
// 				fids2.push_back(nxloc);
			}
// 			else if(max_score<0.85){
// //                 std::cout<<"***************************"<<std::endl;
// //                 std::cout<<nxloc.x<<" "<<nxloc.y<<std::endl;
// //                 std::cout<<"score:"<<max_score<<std::endl;
// //                 std::cout<<"***************************"<<std::endl;
// 				continue;
// 			}
// 			else if(exist)
//             {
// //                 std::cout<<"该点存在："<<nxloc.x<<" "<<nxloc.y<<std::endl;
//             }
//             fids2.push_back(nxloc);
        }
//         std::cout<<"x:"<<std::endl;
//         for(int n=0;n<fids2.size();n++)
//         {
//             std::cout<<fids2[n].x<<","<<std::endl;
//         }
//         std::cout<<"y:"<<std::endl;
//         for(int n=0;n<fids2.size();n++)
//         {
//             std::cout<<fids2[n].y<<","<<std::endl;
//         }
    }
    for(int i=0;i<fids_high->Size();i++)
    {
        std::vector<util::point2d>& fids_all = allfids->V(i);
        std::vector<util::point2d>& fids_add = addfisk->V(i);
        std::vector<util::point2d>& fids_H = fids_high->V(i);
        for(int j=0;j<fids_H.size();j++)
        {
            fids_all.push_back(fids_H[j]);
        }
        for(int j=0;j<fids_add.size();j++)
        {
// 			std::cout<<fids_add[j].x<<"     "<<fids_add[j].y<<std::endl;
            fids_all.push_back(fids_add[j]);
        }
    }
}

void Detector::TransMain2(util::ImgMatchVector& imvector, HSetVector& hsets, util::FiducialStack* fidsk, util::FiducialStack* pro_fidsk)
{
    int num=fidsk->Size();
    pro_fidsk->ReSize(num);
//     std::vector<util::point2d> fids_new=new std::vector<util::point2d>[size];

    for(int i=0;i<hsets.Size();i++)
    {
        std::cout<<"*************"<<i<<std::endl;
        std::vector<util::point2d> fids;
        for(int j=0;j<imvector[i].size();j++)
        {
            util::point2d fid_new;
            fid_new.x=(imvector[i].pairs)[j].first.x;
            fid_new.y=(imvector[i].pairs)[j].first.y;
//             fids_new.push_back(fid_new);
            std::vector<util::point2d> slocs;
            PredictMissFiducial(fid_new, hsets, i, i, slocs);
            fids.push_back(slocs[0]);
//             std::cout<<slocs[0].x<<","<<std::endl;
//             std::cout<<slocs[0].y<<","<<std::endl;
//             std::cout<<"x:"<<slocs[0].x<<" y:"<<slocs[0].y<<std::endl;
        }
        std::cout<<"x:"<<std::endl;
        for(int n=0;n<fids.size();n++)
        {
            std::cout<<fids[n].x<<","<<std::endl;
        }
        std::cout<<"y:"<<std::endl;
        for(int n=0;n<fids.size();n++)
        {
            std::cout<<fids[n].y<<","<<std::endl;
        }

//         std::vector<util::point2d>& fids2 = pro_fidsk->V(i);
//         std::cout<<"time"<<i<<std::endl;
//         std::vector<util::point2d>& fids = fidsk->V(i);
//         for(int j=0;j<fids.size();j++)
//         {
//             std::vector<util::point2d> slocs;
// //             std::cout<<"x:"<<fids[j].x<<" y:"<<fids[j].y<<std::endl;
//             PredictMissFiducial(fids[j], hsets, i, i, slocs);
// //             std::cout<<"x:"<<slocs[0].x<<" y:"<<slocs[0].y<<std::endl;
//             fids2.push_back(slocs[0]);
//         }
    }

}



// void Detector::DetectorMain2(util::MrcStack& mrcr, util::FiducialStack* fidsk, float& diameter, const std::vector< int >& m_index, float ratio)
// {
// 	int scale, size, min_shape, resize_index, mul_para, diameter_int, m, n, radius_int, temp_ori_m;
// 	cv::Mat ori_img, img_resized, img, img1, template1, wave_img, ori_template, fid,ori_fid,score_error,ori_img_draw;
//     cv::Point xy;
//     int dense;
//     Detector detector;
//     fidsk->ReSize(mrcr.Size());
//  	fidsk->SetRatio(ratio);
//     fidsk->SetWxH(mrcr.Width()*ratio, mrcr.Height()*ratio);
//  	diameter *= ratio;
//
//
// 	scale=2;    //deepth
// 	dense = 1;                                                 // number pioints
//  	EX_TRACE("MRC Image rescale ratio: %.2f\n", ratio)
//
//  	cv::Mat mid_slice = mrcr.GetStackImage(mrcr.Size()/2);
// //  	ori_img = mid_slice.clone();
//  	if (dense) size = 2500;
// 	else size = 2000;
// // 	cv::Mat mid_slice = mrcr.GetStackImage(mrcr.Size()/2);
// // 	std::cout<<"ok"<<std::endl;
//  	util::ScaleImage(mid_slice, ratio);
//  	util::Reversal(mid_slice);
//  	util::ConvertTo1(mid_slice, false);
//  	util::MedianSmooth(mid_slice);
// //
//  	float refined_diameter;
//  	cv::Mat avgtmplt;
//  	detector.DetermineDiameter(mid_slice, diameter, &refined_diameter, &avgtmplt);
//  	detector.SetFidDiameter(refined_diameter);
//  	detector.SetSampling(mid_slice, refined_diameter);
//  	detector.SetFidAverageTemplate(avgtmplt);
//  	diameter = refined_diameter;
//
//  	detector.fid_avgtmplt = avgtmplt.clone();
//  	util::ConvertTo1(detector.fid_avgtmplt, false);
//
//  	util::SeriesSaveToFile(avgtmplt, "avgtmplt");
//
//  	util::Reversal(avgtmplt);
//  	util::ConvertTo1(avgtmplt, false);
//  	util::SaveImage(avgtmplt, "avgtmplt(reversal).pgm");
//
// 	radius_int = (int)((refined_diameter + 0.0) / 2 + 0.5);
// // 	std::cout<<radius_int<<std::endl;
// // 	detector.SetFidDiameter(diameter_int);
// 	detector.SetFidRadius(radius_int);
//
// // 	std::cout << template1 <<   std::endl;
//
//
//  	for(int i = 0; i < mrcr.Size(); i++){
//
// 		EX_TIME_BEGIN("%sProcessing MRC[%d]", _DASH, i)
// //
// // // 		dense = 1;
// 		std::vector<util::point2d>& fids = fidsk->V(i);
// 		cv::Mat mid_slice = mrcr.GetStackImage(i);
//
//
// 		ori_img = mid_slice.clone();
//
// 		if (dense) size = 2500;
// 		else size = 2000;
// 		min_shape = ori_img.rows;
// 		if (ori_img.cols < ori_img.rows) {
// 			min_shape = ori_img.cols;
// 		}
// 		resize_index = 0;
// 		if (min_shape > size)
// 		{
// 			for (int i = 2; i < 7; i = i + 2)
// 			{
// 				if (int(min_shape / i) > size)
// 					continue;
// 				if (!(int(min_shape / i) > size))
// 				{
// 					ratio = i;
// 					resize_index = 1;
// 					break;
// 				}
// 			}
// 			cv::resize(ori_img, img_resized, cv::Size(ori_img.cols / ratio, ori_img.rows / ratio), 0, 0, cv::INTER_AREA);
// 		}
// 		else
// 		{
// 			img_resized = ori_img;
// 			ratio = 1;
// 		}
//
// 		cv::GaussianBlur(img_resized, img, cv::Size(5, 5), 0);
// 		img1 = img.clone();       //keep img
// 		util::ToOne(img1);
// // 		std::cout << "template1" <<   std::endl;
// 		detector.Process_new(mid_slice, img1, fids, true);
//
//
//
// 		EX_TIME_END("Processing MRC[%d]", i)
//
// 	}
//
//
// }
// void Detector::DetectorMain2(util::MrcStack& mrcr, util::FiducialStack* fidsk, float& diameter, const std::vector< int >& m_index, float ratio)
void Detector::DetectorMain2(util::MrcStack& mrcr, util::FiducialStack* fidsk, float& diameter, const std::vector<int>& m_index, float ratio)
{
// 	int scale, size, min_shape, resize_index, mul_para, diameter_int, m, n, radius_int, temp_ori_m;
// 	cv::Mat ori_img, img_resized, img, img1, template1, wave_img, ori_template, fid,ori_fid,score_error,ori_img_draw;
// 	cv::Point xy;
// 	int dense;
// 	Detector detector;
// // 	fidsk->ReSize(1);
// 	fidsk->ReSize(mrcr.Size());
// 	fidsk->SetRatio(ratio);
// 	fidsk->SetWxH(mrcr.Width()*ratio, mrcr.Height()*ratio);
// 	diameter *= ratio;
//
// 	dense=1;
// 	scale=2;


// 	float diameter_int, radius_int;
	int scale, size, min_shape, resize_index, mul_para, diameter_int, m, n, radius_int, temp_ori_m;
	cv::Mat ori_img, img_resized, img, img1, template1, wave_img, ori_template, fid,ori_fid,score_error,ori_img_draw;
    cv::Point xy;
    int dense;
    Detector detector;
    fidsk->ReSize(mrcr.Size());
 	fidsk->SetRatio(ratio);
    fidsk->SetWxH(mrcr.Width()*ratio, mrcr.Height()*ratio);
 	diameter *= ratio;


	scale=2;    //deepth
	dense = 0;                                                 // number pioints
 	EX_TRACE("MRC Image rescale ratio: %.2f\n", ratio)



	bool template_all=true;

	if(template_all)     //一个模板
	{
		cv::Mat mid_slice = mrcr.GetStackImage(mrcr.Size()/2);
		util::ConvertTo1(mid_slice, true);
// 		util::ToOne(mid_slice);
		cv::Mat copy(mid_slice.size(), CV_8UC1);
//
		for(int y = 0; y < mid_slice.size().height; y++){
			float* src = (float*)(mid_slice.ptr()+y*mid_slice.step);
			char* start = (char*)(copy.ptr()+y*copy.step);
			for(int x = 0; x < mid_slice.size().width; x++){
				*start++ = (*src++)*255;
			}
		}
// 		std::cout<<"copy"<<std::endl;

		ori_img = copy.clone();
// 		util::ConvertTo1(mid_slice, false);
// 		util::ToOne(ori_img);

// 		std::cout<<"mid_slice"<<mid_slice<<std::endl;
// 		std::cout<<"ori_img"<<ori_img<<std::endl;

		if (dense) size = 2500;
		else size = 2000;
		min_shape = mid_slice.rows;
	// 	    std::cout<<"min_shape"<<min_shape<<std::endl;
		if (ori_img.cols < ori_img.rows) {
			min_shape = ori_img.cols;
		}
		resize_index = 0;
		if (min_shape > size)
		{
			for (int i = 2; i < 7; i = i + 2)
			{
				if (int(min_shape / i) > size)
					continue;
				if (!(int(min_shape / i) > size))
				{
					ratio = i;
					resize_index = 1;
					break;
				}
			}
				cv::resize(ori_img, img_resized, cv::Size(mid_slice.cols / ratio, mid_slice.rows / ratio), 0, 0, cv::INTER_AREA);
		}
		else
		{
			img_resized = ori_img;
			ratio = 1;
		}
// 		util::ConvertTo1(mid_slice, false);
		img_resized.convertTo(img_resized, CV_32FC1);
// 		std::cout<<img_resized.size()<<std::endl;
		cv::GaussianBlur(img_resized, img, cv::Size(5, 5), 0);
// 		std::cout<<img.at<float>(0,0)<<std::endl;
		img1 = img.clone();       //keep img
		util::ToOne(img1);
// 		util::ConvertTo1(img1, false);
// 		std::cout<<img1<<std::endl;
// std::cout<<template1<<std::endl;


		detector.template_make(img1, detector.scale, template1, wave_img, diameter_int);
// 		std::cout<<template1<<std::endl;
		if (template1.rows < 1)
		{
			std::cout << "There is a problem when making template of "<< std::endl;
		}
		detector.template_average(template1);
		cv::Mat tmp = template1.clone();
		std::cout<<diameter_int<<std::endl;
// 		util::Reversal(tmp);
		util::ConvertTo1(tmp, false);
		util::SaveImage(tmp, "avgtmplt(reversal).pgm");
		m=template1.rows;
		n=template1.cols;
		cv::resize(template1, ori_template, cv::Size(ratio * m, ratio * n), 0, 0, cv::INTER_AREA);
		util::SeriesSaveToFile(template1, "avgtmplt");
// 		util::SaveImage(template1, "avgtmplt(reversal).pgm");
	// 	diameter_int=18;
		radius_int = ((diameter_int + 0.0) / 2 + 0.5);
		detector.SetFidDiameter(diameter_int);
		detector.SetFidRadius(radius_int);
		detector.fid_oritmplt = ori_template.clone();
		detector.SetFidAverageTemplate(template1);
		detector.fid_avgtmplt = template1.clone();
		detector.SetWaveimg(wave_img);

	}


// int tmp=0;

	for(int i = 0; i < mrcr.Size(); i++){
// 	for(int i = 0; i < 1; i++){

		EX_TIME_BEGIN("%sProcessing MRC[%d]", _DASH, i)
//
// // 		dense = 1;
		std::vector<util::point2d>& fids = fidsk->V(i);
		cv::Mat mid_slice = mrcr.GetStackImage(i);
		cv::Mat ori_img;
// std::cout<<mid_slice<<std::endl;

// 		ori_img = mid_slice.clone();


		util::ConvertTo1(mid_slice, true);
// 		util::ToOne(mid_slice);
		cv::Mat copy(mid_slice.size(), CV_8UC1);
//
		for(int y = 0; y < mid_slice.size().height; y++){
			float* src = (float*)(mid_slice.ptr()+y*mid_slice.step);
			char* start = (char*)(copy.ptr()+y*copy.step);
			for(int x = 0; x < mid_slice.size().width; x++){
				*start++ = (*src++)*255;
			}
		}
// 		std::cout<<copy.at<char>(0,0)<<std::endl;

		ori_img = copy.clone();







		if (dense) size = 2500;
		else size = 2000;
		min_shape = ori_img.rows;
		if (ori_img.cols < ori_img.rows) {
			min_shape = ori_img.cols;
		}
		resize_index = 0;
		if (min_shape > size)
		{
			for (int i = 2; i < 7; i = i + 2)
			{
				if (int(min_shape / i) > size)
					continue;
				if (!(int(min_shape / i) > size))
				{
					ratio = i;
					resize_index = 1;
					break;
				}
			}
			cv::resize(ori_img, img_resized, cv::Size(ori_img.cols / ratio, ori_img.rows / ratio), 0, 0, cv::INTER_AREA);
		}
		else
		{
			img_resized = ori_img;
			ratio = 1;
		}

		img_resized.convertTo(img_resized, CV_32FC1);
// 		std::cout<<img_resized<<std::endl;
		cv::GaussianBlur(img_resized, img, cv::Size(5, 5), 0);
		img1 = img.clone();       //keep img
		util::ToOne(img1);

// 		std::cout << img1 << std::endl;



		if(!template_all)
		{
			detector.template_make(img1, scale, template1, wave_img, diameter_int);
// 			std::cout << template1 << std::endl;
	// 		std::cout << "diameter_int:" << diameter_int << std::endl;
			detector.template_average(template1);
	//
			m=template1.rows;
			n=template1.cols;
			cv::resize(template1, ori_template, cv::Size(ratio * m, ratio * n), 0, 0, cv::INTER_AREA);
			util::SeriesSaveToFile(template1, "template1");
			util::SaveImage(template1, "template1(reversal).pgm");
			detector.fid_oritmplt = ori_template.clone();
// 			detector.SetFidAverageTemplate(template1);
// 			detector.fid_avgtmplt = template1.clone();
// 			radius_int = (int)((diameter_int + 0.0) / 2 + 0.5);


			detector.fid_avgtmplt = template1.clone();
			detector.SetFidAverageTemplate(template1);
			radius_int = (int)((diameter_int + 0.0) / 2 + 0.5);
			detector.SetFidDiameter(diameter_int);
			detector.SetFidRadius(radius_int);
			detector.SetWaveimg(wave_img);


// 		std::cout<<detector.wave_img<<std::endl;

// 		detector.SetFidDiameter(diameter_int);
// 		detector.markerauto_work_flow(fid, img, template1, radius_int, diameter_int, wave_img);


// 		ori_fid = fid * ratio;
// 		cv::Mat xy_refine;
// 		temp_ori_m=ori_template.rows;
// // 		location_fid(ori_img, ori_fid, temp_ori_m,fid,score_error);
// // 		std::cout<<"mid_slice:"<<mid_slice<<std::endl;
// 		std::cout<<"ori_fid:"<<ori_fid<<std::endl;
// 		std::cout<<"temp_ori_m:"<<temp_ori_m<<std::endl;
// 		std::cout<<"xy_refine:"<<xy_refine<<std::endl;
// 		detector.location_fid(fids, mid_slice, ori_fid, temp_ori_m, xy_refine, score_error);


		}




		detector.Process_new(mid_slice, img1, img, fids, true, ratio);


		EX_TIME_END("Processing MRC[%d]", i)

	}
}














