#ifndef DETECTOR_H__
#define DETECTOR_H__

#include <opencv2/opencv.hpp>
// #include "opencv2/core/core_c.h"
// #include "opencv2/imgproc/imgproc_c.h"
#include<opencv2/core/core.hpp>
#include<opencv2/imgproc/imgproc.hpp>
#include<opencv2/highgui/highgui.hpp>
#include "dataf/dataf.h"
#include "mrcimg/mrc2img.h"
#include "modelmatch/match_core.h"
#include <ANN/ANN.h>


/** double image size before pyramid construction? */
#define IMG_DBL				 1
#define DEFAULT_FIDD		 10//30//10
#define CORR_THRESHOLD		 0.5f

class Detector{
private:
class ImgRef{
public:
	float value;
	int x;
	int y;
public:
	ImgRef(){}
	ImgRef(cv::Mat& __img, int __x, int __y):value(((float*)(__img.data + __img.step*__y))[__x]), x(__x), y(__y){}
	float& Value(){return value;}
	bool operator < (const ImgRef& ref) const{
		return value < ref.value;
	}
	bool operator > (const ImgRef& ref) const{
		return value > ref.value;
	}
	bool operator <= (const ImgRef& ref) const{
		return value <= ref.value;
	}
	bool operator >= (const ImgRef& ref) const{
		return value >= ref.value;
	}
	const ImgRef& operator = (const ImgRef& ref){
		value = ref.value;
		x = ref.x;
		y = ref.y;
		return *this;
	}
};

struct LocalSeed{
public:
	int z_begin, z_next;			//match direct is z_begin to z_next;
	int z_img_idx;
	bool is_final;
	util::point2d z_bp;
	util::point2d z_np;
};

private:
    float fid_diameter;				//refer to the original mrc
	float fid_radius;
	float corr_threshold;
	float pixel_threshold;
	int sampling;
	int scale=2;
    cv::Mat fid_tmplt;
    cv::Mat fid_avgtmplt;			//saved tmplt is the same to scale_ratio
    cv::Mat fid_oritmplt;			//saved tmplt is the same to scale_ratio
    cv::Mat wave_img;
    std::vector<util::point2d> fids;

    static Detector detector;

private:
    void FindMaxPeak(cv::Mat& corr, int seed_x, int seed_y, int idiameter, int* peak_x, int* peak_y);
    float GetAverageOfMarkerPatch(cv::Mat& img, int seed_x, int seed_y, int diameter);
    void GetRawPositions(cv::Mat& img, std::vector<ImgRef>& refvec, float diameter, std::vector<util::point2d>& raw_fids, std::vector<std::pair<float, float> >& scores);
    void RefineFiducialMarkersByGaussianDistribution(const std::vector<util::point2d>& raw_fids,
													 const std::vector<std::pair<float, float> >& scores, std::vector<util::point2d>& new_fids);
    void CreateTemplate(cv::Mat* tmplt , float diameter);
    void GaussianSmoothBasedOnSampling(cv::Mat& img, int sampling);
    void CalculateCorralation(const cv::Mat& src, const cv::Mat& tmplt, cv::Mat* dst);
    bool RegionGrow(cv::Mat& img, float diameter, const util::point2d& seed, util::RECT& region);			//change the img
    void BlindRegionGrow(cv::Mat& img, int seed_x, int seed_y, util::RECT& region);		//change the img
    bool GetCenter(cv::Mat& img, const util::point2d& seed, float diameter, util::point2d& centre);
    void FindFiducialMarkerPositions(cv::Mat& img, const cv::Mat& tmplt, float diameter, std::vector<util::point2d>& fids, float* positive_ratio, bool limit_contrast = false, bool forsubtomo = true);//true);
    void GenerateDiameterSpace(std::vector<float>& diavec, float dia_begin, float dia_end, int bin_size = 11);
    void CreateAverageTemplate(cv::Mat& img, const std::vector<util::point2d>& fids, float diameter, cv::Mat* tmplt, bool zoom_out = false, int magnif = 16);
    void RefinePositionsForSubTomo(cv::Mat& img, const cv::Mat& tmplt, std::vector<util::point2d>& fids, float diameter, int magnif = 16);
	float compute_gauss_error(cv::Mat pic,cv::Mat dir);


	void waveletprocess2(cv::Mat img, cv::Mat& wave_img, int J);
    cv::Mat atrous_cv(cv::Mat ai, cv::Mat k);
	cv::Mat wavelet(cv::Mat a,cv::Mat A);
    void removeSmall(cv::Mat img,double threshold);
	double roundness(cv::Mat label_img);
	cv::Mat refine_fid_by_gaussian_distribution_markerauto_no_wave(cv::Mat candidate);
	cv::Mat refine_fid_by_gaussian_distribution_markerauto_wave(cv::Mat candidate);
	cv::Mat compute_center_Gauss(cv::Mat pic);
	cv::Mat Gauss_A(cv::Mat loc);
	cv::Mat Gauss_b(cv::Mat info);


private:
	Detector();
	~Detector();

private:
    void SetFidDiameter(float __fid_diameter);
	void SetFidRadius(float __radius);
    void SetSampling(cv::Mat& img, float __fid_diameter);
	void SetWaveimg(cv::Mat& Wave);
    void SetFidAverageTemplate(const cv::Mat& __fid_tmplt);
    float InitEstimateDiameter(const cv::Mat& img);
    bool DetermineDiameter(const cv::Mat& img, float diameter, float* refined_diameter, cv::Mat* avgtmplt);
    void Process(cv::Mat& img, std::vector<util::point2d>& fids, bool use_avgtmplt = true);

	void template_make(cv::Mat img, int scale, cv::Mat &template1, cv::Mat &wavelet_img, int &d_end);
	void template_average(cv::Mat& template1);
	void markerauto_work_flow(cv::Mat& fidstack, cv::Mat img_ori, cv::Mat template_ori, int radius_int, int diameter_int, cv::Mat wave_img);
	void location_fid(std::vector<util::point2d>& fids, cv::Mat img, cv::Mat location, float width, cv::Mat &refined_xy, cv::Mat &score_xy);
	void Process_new(cv::Mat& slice, cv::Mat& copy, cv::Mat& img, std::vector<util::point2d>& fids, bool use_avgtmplt = true, float ratio = 1);
	void Wavelet(cv::Mat img, int scale, cv::Mat &wavelet_img);
	void RefinePositions(cv::Mat& img, const cv::Mat& tmplt, std::vector<util::point2d>& fids, float diameter);//true);

    float FindLocalPeak(cv::Mat& img, const cv::Mat& tmplt, const util::point2d& seed, float roi_radio, util::point2d* loc);

private:
    static void InitLocalDetector(float diameter, bool use_avgtmplt = true);
	static float LocalMarkerDetect(cv::Mat& img, const util::point2d& seed, util::point2d* loc);
    static void PredictMissFiducial(const util::point2d& seed, HSetVector& hsets, int z_idx1, int z_idx2, std::vector<util::point2d>& slocs);


public:
    static void DetectorMain(util::MrcStack& mrcs, util::FiducialStack* fidstk, float& diameter, const std::vector<int>& m_index, float ratio = 1);
// 	static void DetectorMain(util::MrcStack& mrcs, util::FiducialStack* fidstk, float& diameter, const std::vector<int>& m_index, float ratio = 1);
	static void DetectorMain2(util::MrcStack& mrcs, util::FiducialStack* fidstk, float& diameter, const std::vector<int>& m_index, float ratio = 1);
    static void LocalDetectorMain(util::MrcStack& mrcs, util::TrackSpace& trackspace, HSetVector& hsets, float diameter, util::FiducialStack* fidsk, util::ImgMatchVector* imvector, const std::vector<int>& m_index);

	static void TransMain(util::MrcStack& mrcs, HSetVector& hsets, util::FiducialStack* fids_high, util::FiducialStack* fidsk, util::FiducialStack* addfisk, util::FiducialStack* pro_fidsk,util::FiducialStack* allfids, float diameter, const std::vector<int>& m_index);
    static void TransMain2(util::ImgMatchVector& imvector, HSetVector& hsets, util::FiducialStack* fidsk, util::FiducialStack* pro_fidsk);

    static void Test(util::MrcStack& mrcs, const util::FiducialStack& fidstk, float diameter, const char* folder= "fids");
};


#endif
