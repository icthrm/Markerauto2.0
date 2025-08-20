#ifndef SFM_DRIVER_H__
#define SFM_DRIVER_H__

#include "img_projs.h"
#include "dataf/quat_vec.h"
#include "params.h"
#include "dataf/camera.h"
#include "geometry_data.h"
#include "matrix/vector.h"
#include "sba.h"
#include "lm.h"

#define FULLQUATSZ	4

namespace bundle {

class PSFMDriver {
private:
    //defination of parameters used in sba API
    static const int cnp = 6; 					///s, alpha, beta, gamma, t0, t1
    static const int pnp = 3; 				///< euclidean 3D points
    static const int mnp = 2; 				///< image points are 2D
    double opts[SBA_OPTSSZ], info[SBA_INFOSZ+1];
    struct pglobs glob;
    static const int verbose = 0;
    static constexpr double min_proj_error_threshold = 1.2;//1.2;//1.5;
    static constexpr double max_proj_error_threshold = 3.2;//2.5;//4.5;

private:
    void Initialize();

    int Run(double* motstruct, double* imgpts, char* vmask,
            const int n3Dpts, const int ncams, const int nconcam, const int ncon3Dpts);
    
        
    int Run_ceres(double* motstruct, double* imgpts, char* vmask,
            const int n3Dpts, const int ncams, const int nconcam, const int ncon3Dpts, std::vector<std::vector<double>>& W_init, bool use_MAND);
    
    int Run_l1(double* motstruct, double* imgpts, char* vmask,
            const int n3Dpts, const int ncams, const int nconcam, const int ncon3Dpts);

    void sba_L1_smoothing_grad(double *p, double *imgpts, 
            char *vmask, int ncams, int n3Dpts,double *newp);

public:
    PSFMDriver();
    ~PSFMDriver();

    void SetGlobMask(bool s, bool alpha, bool beta, bool gamma, bool t0, bool t1, bool p);
//     double Run(int* added_order, mx::pproj_params* pparams, const int start_camera, const int ncams, const int nconcam,
//                v3_t* init_pts, const int n3Dpts, const int ncon3Dpts, std::vector<ImageKeyVector>& pt_views, GeometryData& data, bool remove_outliers = true);
    
    double Run(int* added_order, mx::pproj_params* pparams, const int start_camera, const int ncams, const int nconcam,
               v3_t* init_pts, const int n3Dpts, const int ncon3Dpts, std::vector<ImageKeyVector>& pt_views, std::vector<std::vector<double>>& W_init, GeometryData& data, bool remove_outliers = true, bool use_MAND=true);

    void RefineCameraParameters(v3_t* points, v2_t* projs, int num_points, mx::pproj_params* pparams);

    double compute_norm1(double *a,int len);

    void Toone_recover_new(double *motstr,double *imgpts,int ncams,int n3Dpts,int n2Dpts,double mean_data,double std_data,double *gamma_old,double *motstr_recover,double *imgpts_recover);

    void double_copy(double *changed_arr,double *notchanged_arr,int len);

    void data_initial(double* motstruct, double* imgpts, char* vmask, double* vec,
            int ncams,int n3Dpts,double ratio_noise,double ratio_outlier,double ratio_noise_cams,double range_img, double& mean, double& std);

    double Armijo_f(double *p, double mu, double* grad, double sigma,
            double rho, double *imgpts, char * vmask, double *weight, int ncams, int n3Dpts);

    void update_p(double *p, double mu, double *grad, double sigma, double rho, double *imgpts,
            char *vmask, double *weight, int ncams, int n3Dpts);

    void compute_gradient(double *p, double *imgpts, double mu, char *vmask, double *weight,
            int ncams, int n3Dpts, double *grad);

    void compute_uv(double *camera, double *point, double *uv);

    void compute_right_R(double *camera,double *point,double *right_R_u,double *right_R_v);
    // void compute_dfdR(double*p,double *imgpts,char *vmask,double*weight,double mu,int ncams,int n3Dpts,int k,double *dfdR);
    
    void compute_right_X(double *camera,double *point,double *right_X_u,double *right_X_v);
    
    void compute_dfdX(double*p,double *imgpts,char *vmask,double *weight,double mu,int ncams,int n3Dpts,int k,double *dfdX);

    void compute_dfdR(double*p,double *imgpts,char *vmask,double*weight,double mu,int ncams,int n3Dpts,int k,double *dfdR);

    double f(double *p, double mu, double *imgpts, char *vmask, double *weight, int ncams, int n3Dpts);

    void f(double *p, double mu, double *imgpts, char *vmask, double *weight, int ncams, int n3Dpts, std::string filename);

    double f_L2(double *p, double *imgpts, char *vmask, double *weight, int ncams, int n3Dpts);

    void write_txt(int *a,int len,std::string filename);
    void write_txt(double *a,int len,std::string filename);
    void write_txt(char *a,int len,std::string filename);
    // void write_txt_app(T *a,int len,std::string filename);
    void Toone_new(double *one_avg,double *one_std,double *mot,double *imgpts,char *vmask,int ncams,int n3Dpts);


    static void Project(const double param[6], const double M[3], double n[2]);
    static v2_t Project(const mx::pproj_params& pparam, v3_t pt);

    int Cnp() const {
        return cnp;
    }
    int Mnp() const {
        return mnp;
    }
    int Pnp() const {
        return pnp;
    }
    const double* Info() const {
        return info;
    }
};

class PosePointParametersBlock
{
public:
    PosePointParametersBlock(){}
    void create(int pose_num, int point_num)
    {
        poseNum = pose_num;
        pointNum = point_num;
        values = new double[pose_num * 6 + point_num * 3];
    }
    PosePointParametersBlock(int pose_num, int point_num): poseNum(pose_num), pointNum(point_num)
    {
        values = new double[pose_num * 6 + point_num * 3];
    }
    ~PosePointParametersBlock() { delete[] values; }

    double* pose(int idx) {  return values + idx * 6; }

    double* point(int idx) { return values + poseNum * 6 + idx * 3; }


    int poseNum;
    int pointNum;
    double *values;    //存储全部的信息

};

}

#include "sfm_driver2.h"

#endif
