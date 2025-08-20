#ifndef BUNDLE_CORE_H__
#define BUNDLE_CORE_H__

#include "dataf/dataf.h"
#include "params.h"
#include "sfm_driver.h"
#include "util/exception.h"
#include "geometry_data.h"
#include <vector>
#include "dataf/camera.h"



/** The implement of bundle adjustment.*/
class PPBundleApp{
private:
    bundle::PSFMDriver driver;
    bundle::GeometryData data;		///< the used data
    int width, height;
    
    const static int m_min_max_matches = 16;           /* Minimum number of matches needed to register an image */
    
private:
    mx::pproj_params* result;
    std::vector<bundle::PointData> m_point_data;   		///> Information about 3D points in the scene
    
private:
    PPBundleApp();
    ~PPBundleApp();

private:
    void InitGeometryData(util::TrackSpace& trackspace);
    
    void InitBoundary(int _width, int _height);
    
    void InitCameraParam(float angle, mx::pproj_params* camera) const;
	
	void InitCameraParam(double alpha, double beta, double gamma, double t0, double t1, mx::pproj_params* camera);
    
    void InitMotParam(float angle, double* params) const;
    
    /** @brief Triangulate n points, mots: s, alpha, beta, gamma, t0, t1 */
	v3_t Triangulate(int num_points, const v2_t* pv, const double* mots, double *error_out);
	
	v3_t Triangulate(v2_t p, v2_t q, const mx::pproj_params& c1, const mx::pproj_params& c2, double* proj_error);
    /** @brief Triangulate a subtrack */
    v3_t TriangulateNViews(const bundle::ImageKeyVector& views, int* added_order, mx::pproj_params* cameras, double& error);
    
    void Process(std::vector<mx::pproj_params>* cameras_tmp,float percent = 0 ,bool first = true);		//percent indicates the length threshold of tracks		0~1
    void PrintCameras(const char* filename = "cameras.params") const;
    void PrintPoints(const char* filename = "points.mot") const;
    void PrintTrackBin(const char* filename = "track.hist") const;
    void PrintRepErrorBin(mx::pproj_params* camparams, int ncams, const int added_order[],
				 const v3_t* points, const char* foldername = "verbose") const;
    void DumpOutCameras(std::vector<mx::pproj_params>* cameras) const;
    void DumpOutPoints(std::vector<v3_t>* points) const;
//     void create(const std::vector< float >& angles, std::vector<mx::pproj_params>* cameras, std::vector<v3_t>& points) ;
    
public:
    static void BundleMain(util::TrackSpace& tspace_, int width_, int height_, std::vector<mx::pproj_params>* cameras, std::vector<mx::pproj_params>* cameras_tmp, std::vector<v3_t>* points, bool first=true);
    
    static void CreateMain(const std::vector< float >& angles, std::vector<mx::pproj_params>* cameras, std::vector<v3_t>& points, 
                                   util::FiducialStack* addfids);
	
	static void PrintCamerasAsIMOD(const std::vector<mx::pproj_params>& cameras, const float rotation, const float scale, const char* xffilename, 
								   const char* xanglefilename, const char* yanglefilename, const char* exfilename);

public:
    static void Test();
};

#endif
