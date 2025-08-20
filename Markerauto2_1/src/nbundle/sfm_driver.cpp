#include "sfm_driver.h"
#include "img_projs.h"
#include <cmath>
// #include <random> 
#include <iostream>
#include "util/exception.h"
#include "util/matrix.h"
#include "matrix/matrix.h"
#include "util/qsort.h"
#include <cfloat>
#include <matrix/vector.h>
#include "geometry_data.h"
#include <ext/algorithm>
#include <numeric>
#include <vector>
#include <fstream> // file
#include <lm.h>
#include "ceres/ceres.h"
// #include <ext/hash_map>
#include <random>  // Add this at the top of your file

#define MAXITER 		800

// using ceres::AutoDiffCostFunction;
// using ceres::CostFunction;
// using ceres::Problem;
// using ceres::Solve;
// using ceres::Solver;

// A templated cost functor that implements the residual r = 10 -
// x. The method operator() is templated so that we can then use an
// automatic differentiation wrapper around it to generate its
// derivatives.

using namespace mx;
using namespace std;

bundle::PSFMDriver::PSFMDriver()
{
    Initialize();
}
bundle::PSFMDriver::~PSFMDriver() {}

void bundle::PSFMDriver::Initialize()
{
    /* call sparse LM routine */
    opts[0]=SBA_INIT_MU;
    opts[1]=SBA_STOP_THRESH;
    opts[2]=SBA_STOP_THRESH;
    opts[3]=SBA_STOP_THRESH;
    //opts[3]=0.05*numprojs; // uncomment to force termination if the average reprojection error drops below 0.05
    opts[4]=0.0;
    //opts[4]=1E-12; // uncomment to force termination if the relative reduction in the RMS reprojection error drops below 1E-05
}

#ifndef MIN_POINTS
#define MIN_POINTS		6
#endif

void bundle::PSFMDriver::Project(const double param[6], const double M[3], double n[2])
{
    calcImgParallelProj(param, M, n);
}

void bundle::PSFMDriver::SetGlobMask(bool s, bool alpha, bool beta, bool gamma, bool t0, bool t1, bool p)
{
	glob.cpmask[0] = s;
	glob.cpmask[1] = alpha;
	glob.cpmask[2] = beta;
	glob.cpmask[3] = gamma;
	glob.cpmask[4] = t0;
	glob.cpmask[5] = t1;
	glob.pmask = p;
}

template<typename T>
void write_txt_app(T *a,int len,std::string filename){
    int i;
    std::ofstream out;
    out.open(filename,std::ios_base:: app);
    for(i=0;i<len;i++) out<<a[i]<<std::endl;
    out.close();
}

void bundle::PSFMDriver::write_txt(double *a,int len,std::string filename){
    int i;
    std::ofstream out;
    out.open(filename);
    for(i=0;i<len;i++) out<<a[i]<<std::endl;
    out.close();
}

void bundle::PSFMDriver::write_txt(char *a,int len,std::string filename){
    int i;
    std::ofstream out;
    out.open(filename);
    for(i=0;i<len;i++) out<<a[i]+0<<std::endl;
    out.close();
}

void bundle::PSFMDriver::write_txt(int *a,int len,std::string filename){
    int i;
    std::ofstream out;
    out.open(filename);
    for(i=0;i<len;i++) out<<a[i]<<std::endl;
    out.close();
}

void bundle::PSFMDriver::double_copy(double *changed_arr,double *notchanged_arr,int len){
    for(int i=0;i<len;i++) changed_arr[i] = notchanged_arr[i];
}

void bundle::PSFMDriver::Toone_new(double *one_avg,double *one_std,double *mot,double *imgpts,char *vmask,int ncams,int n3Dpts){
    //gui yi hua 
    int i,j,num_img,num_keys;
    double mean_temp,std_temp,sin_gamma_i,cos_gamma_i;
    std::vector<float> points2d;

    j=0;
    for(i=0;i<n3Dpts*ncams;i++){
        if(!((int)vmask[i] == 0)){
            points2d.push_back(imgpts[j]);
            j++;
            points2d.push_back(imgpts[j]);
            j++;
        }
    }
    for(i=0;i<ncams;i++){
        points2d.push_back(mot[i*6+4]);
        points2d.push_back(mot[i*6+5]);
        //points2d.push_back(camparams[i].s);
    }
    mean_temp = 0;
    for(i=0;i<points2d.size();i++){
        mean_temp+=points2d[i];
    }
    mean_temp=mean_temp/points2d.size();
    std_temp = 0;
    for(i=0;i<points2d.size();i++){
        std_temp+=pow((points2d[i]-mean_temp),2);
    }
    std_temp=std_temp/points2d.size();
    std_temp=sqrt(std_temp);//compute mean and std

    j=0;
    for(i=0;i<n3Dpts*ncams;i++){
        if(!((int)vmask[i] == 0)){
            imgpts[j]=(imgpts[j]-mean_temp)/std_temp;
            j++;
            imgpts[j]=(imgpts[j]-mean_temp)/std_temp;
            j++;
        }
    }
    for(i=0;i<ncams;i++){
        sin_gamma_i = sin(mot[6*i+3]);
        cos_gamma_i = cos(mot[6*i+3]);
        mot[6*i+4]=(mot[6*i+4] + mean_temp*(cos_gamma_i+sin_gamma_i))/std_temp;
        mot[6*i+5]=(mot[6*i+5] + mean_temp*(-sin_gamma_i+cos_gamma_i))/std_temp;
    }
    for(i=6*ncams;i<6*ncams+3*n3Dpts;i++){
        mot[i]=mot[i]/std_temp;
    }
    *one_avg=mean_temp;
    *one_std=std_temp;
}

static bool pairCompare(const std::pair<double, std::pair<int, int> >& firstElem, const std::pair<double, std::pair<int, int> >& secondElem){
  return firstElem.first > secondElem.first;
}

/** the first index of @c pt_views is the index of @c added_order*/
double bundle::PSFMDriver::Run(int* added_order, mx::pproj_params* pparams, const int start_camera, const int ncams, const int nconcam,
                              v3_t* init_pts, const int n3Dpts, const int ncon3Dpts, std::vector<ImageKeyVector>& pt_views, std::vector<std::vector<double>>& W_init, GeometryData& data, bool remove_outliers, bool use_MAND)
{
    int total_removed_points = 0;
    int num_outliers = 0;

    double dist_total = 0.0;
    int num_dists = 0;
	double global_error;

    int *remap = new int [n3Dpts];		//store the index of 3D points in nz_pts, if no, store -1
    double* motstruct = new double[ncams*cnp + n3Dpts*pnp];
    double* nz_pts = motstruct+ncams*cnp;

    char* vmask = new char[n3Dpts*ncams];

    int num_projections = 0;
    for(int i = 0; i < n3Dpts; i++) {
        num_projections += (int)pt_views[i].size();
    }

    double* projections = new double[mnp*num_projections];
	int totalkeys[ncams];
	int outacc[ncams];

    int num_3d_pts;
    

    do {
        if((num_3d_pts = n3Dpts - total_removed_points) < MIN_POINTS) {
            EX_PRINT("# Too few points remaining, exiting!\n")

            dist_total = DBL_MAX;
            break;
        }

        int arr_idx = 0;
        int nz_count = 0;
		memset(totalkeys, 0, sizeof(int)*ncams);

        /* Set up the vmask and projections */
        memset(vmask, 0, sizeof(char)*num_3d_pts*ncams);

        int fixed_pt_num = ncon3Dpts;

        for(int i = 0; i < n3Dpts; i++) {
            int num_views =(int)pt_views[i].size();

            if(num_views > 0) {
                for(int j = 0; j < num_views; j++) {
                    int c = pt_views[i][j].first;
                    int v = added_order[c];
                    int k = pt_views[i][j].second;
                    vmask[nz_count * ncams + c] = 1;
					totalkeys[c]++;										//patch
                    projections[2 * arr_idx + 0] = data.GetKey(v,k).m_x;
                    projections[2 * arr_idx + 1] = data.GetKey(v,k).m_y;

                    arr_idx++;
                }

                remap[i] = nz_count;
                memcpy(&nz_pts[nz_count*3], init_pts[i].p, sizeof(double)*3);
                nz_count++;
            }
            else {
                if(i < ncon3Dpts) {
                    fixed_pt_num--;
                }
                remap[i] = -1;
            }
        }

        double* motparams = motstruct;
        for(int i = 0; i < ncams; i++) {
            mx::MotCopyFormPProjParams(motparams, pparams[i+start_camera]);
            motparams += cnp;
        }

        dist_total = 0.0;
        num_dists = 0;

        std::ofstream fffout;
        fffout.open("cams.txt");
        for(int i=0;i<ncams*cnp + n3Dpts*pnp;i++){
            fffout<<motstruct[i]<<std::endl;
        }
        fffout.close();
        fffout.open("pts2d.txt");
        for(int i=0;i<mnp*num_projections;i++){
            fffout<<projections[i]<<std::endl;
        }
        fffout.close();
        fffout.open("vmask.txt");
        for(int i=0;i<n3Dpts*ncams;i++){
            fffout<<(int)vmask[i]<<std::endl;
        }
        fffout.close();

//         std::string res_file = "/home/xzh/文档/markerauto/mk_all/build/bin/projections.txt";
//         std::ofstream outputfile_res(res_file);
//
//         for(int i=0;i<ncams*cnp + n3Dpts*pnp;i++)
//         {
//             outputfile_res<<motstruct[i]<<std::endl;
//     //         if(i%20==0) outputfile_res<<"\n";
//         }
//
//         outputfile_res.close();
        use_MAND = false;
        EX_BEGIN_CLOCK()
        if(1){
            Run_ceres(motstruct, projections, vmask, num_3d_pts, ncams, nconcam, fixed_pt_num, W_init, use_MAND);
        }
        else Run_l1(motstruct, projections, vmask, num_3d_pts, ncams, nconcam, fixed_pt_num);
        std::cout<<"Run_l1 END"<<std::endl;
//         Run(motstruct, projections, vmask, num_3d_pts, ncams, nconcam, fixed_pt_num);
        // Run_l1;
        EX_END_CLOCK()
        EX_TRACE("# SFM using %d 3D pts(%d fixed), %d frames(%d fixed) and %d image projections(%g p/p), error %g [initial %g](elapse: %ld)\n",
                 num_3d_pts, ncon3Dpts, ncams, nconcam, arr_idx, ((double)arr_idx)/num_3d_pts, sqrt(info[1]/arr_idx)>1?sqrt(info[1]/arr_idx):info[1]/arr_idx, 
				 sqrt(info[0]/arr_idx)>1?sqrt(info[0]/arr_idx):info[0]/arr_idx, EX_ELAPSE());
		
		global_error = sqrt(info[1]/arr_idx);

        motparams = motstruct;
        for(int i = 0; i < ncams; i++){
            mx::PProjPCopyFormMotParams(&pparams[i+start_camera], motparams);
            motparams += cnp;
        }
        double abs_re2=0;
        std::vector<std::pair<int, int> > outliers;
		std::vector<std::pair<double, std::pair<int, int> > > outidx;

        if(!remove_outliers){
            goto end;
        }

        for(int i = 0; i < ncams; i++){
            double params[6];
			mx::MotCopyFormPProjParams(params, pparams[i+start_camera]);

            int num_keys = data.GetNumKeys(added_order[i]);

            int num_pts_proj = 0;
            for(int j = 0; j < num_keys; j++) {
                if(data.GetKey(added_order[i], j).m_extra >= 0) {
                    num_pts_proj++;
                }
            }

            double *dists = new double[num_pts_proj];
            int pt_count = 0;

            std::vector<Keypoint>::iterator iter;

            for(iter = data.m_image_data[added_order[i]].m_keys.begin(); iter != data.m_image_data[added_order[i]].m_keys.end(); iter++) {
                const Keypoint &key = *iter;

                if(key.m_extra >= 0) {
                    int pt_idx = key.m_extra;
                    double X[3], pr[2];
                    memcpy(X, &nz_pts[remap[pt_idx]*3], sizeof(double)*3);


                    Project(params, X, pr);
                    abs_re2=pr[0]-key.m_x+pr[1]-key.m_y+abs_re2;
                    double dx = pr[0]-key.m_x;
                    double dy = pr[1]-key.m_y;

                    double dist = sqrt(dx * dx + dy * dy);
                    dist_total += dist;
                    num_dists++;

                    dists[pt_count] = dist;

                    pt_count++;
                }
            }

            /* Estimate the median of the distances */
            double med = kth_element_copy(num_pts_proj, int(0.5/*0.7/*0.8/* 0.9 */* num_pts_proj), dists);

#define NUM_STDDEV 2//2.0//3.0//6.0
            double thresh = 1.2 * NUM_STDDEV * med;/* k * stddev */
            thresh = CLAMP(thresh, min_proj_error_threshold, max_proj_error_threshold);
			if(global_error > max_proj_error_threshold){		//considering that the global error; large noise 
				thresh = global_error;
			}

            /* Compute the average reprojection error for this camera */
            double sum = 0.0;
            for(int j = 0; j < num_pts_proj; j++){
                sum += dists[j];
            }

            double avg = sum/num_pts_proj;
            EX_PRINT("# Camera %d[%d] (%d pts), mean error: %0.3f [median %0.3f(0.7 quantile %0.3f), error thresh %0.3f]\n",
					 i, added_order[i], num_pts_proj, avg, med, kth_element_copy(num_pts_proj, int(0.7 * num_pts_proj), dists), thresh);
                    // i, added_order[i], num_pts_proj, avg, kth_element_copy(num_pts_proj, int(0.5 * num_pts_proj), dists), med, thresh);

            pt_count = 0;
			outidx.clear();
            for(int j = 0; j < num_keys; j++) {
                int pt_idx = data.GetKey(added_order[i],j).m_extra;

                if(pt_idx < 0)
                    continue;
				
                if(dists[pt_count] > thresh){ //|| dists[pt_count] > max_proj_error_threshold) {
                    /* Remove this point from consideration */
//                     outliers.push_back(std::pair<int, int>(pt_idx, i));		
					outidx.push_back(std::make_pair(dists[pt_count], std::pair<int, int>(pt_idx, i)));		//? sort of outliers??? WARNING
                }
                pt_count++;
            }
			
			std::sort(outidx.begin(), outidx.end(), pairCompare);
			
			for(int i = 0; i < outidx.size(); i++){
				outliers.push_back(outidx[i].second);
			}
            
            delete [] dists;
        }


        /* Remove outlying points */
        num_outliers = 0;
		memset(outacc, 0, sizeof(int)*ncams);
        for(int i = 0; i <(int)outliers.size(); i++) {
            int idx = outliers[i].first;

            if(idx < ncon3Dpts) {
                continue;
            }

            if(!pt_views[idx].size()) {
                continue;
            }

            for(ImageKeyVector::iterator itr = pt_views[idx].begin(); itr != pt_views[idx].end(); ) {
                int v = (*itr).first;
                int k = (*itr).second;
                if(v == outliers[i].second && totalkeys[added_order[v]]-outacc[added_order[v]] > 6){
                    if(data.GetKey(added_order[v], k).m_extra != idx) {
                        EX_ERROR("Error!  Entry for(%d,%d) should be %d, but is %d\n",
                                 added_order[v], k, idx, data.GetKey(added_order[v], k).m_extra);
                    }
                    data.GetKey(added_order[v], k).m_extra = -2;
                    pt_views[idx].erase(itr);
                    num_outliers++;
					outacc[added_order[v]]++;
                    break;
                }
                else {
                    itr++;
                }
            }
            
            for(int i = 0; i < ncams; i++){
				totalkeys[i] -= outacc[i];
            }

            if(pt_views[idx].size() < 2) {
                for(ImageKeyVector::iterator itr = pt_views[idx].begin(); itr != pt_views[idx].end(); itr++) {
                    int v = (*itr).first;
                    int k = (*itr).second;

                    data.GetKey(added_order[v], k).m_extra = -2;
                }

                pt_views[idx].clear();
                total_removed_points++;
            }
        }

        outliers.clear();

        EX_PRINT("# Removing %d outliers\n", num_outliers);

end:
        for(int i = 0; i < n3Dpts; i++) {
            if(remap[i] != -1) {
                memcpy(init_pts[i].p, &nz_pts[remap[i]*3], sizeof(double)*3);
            }
        }

    } while(num_outliers > 0);

    delete [] vmask;
    delete [] projections;

    delete [] remap;
    delete [] motstruct;

//     std::cout<<"num_dists:"<<num_dists<<"   dist_total"<<dist_total<<std::endl;

    return dist_total/num_dists;
}

v2_t bundle::PSFMDriver::Project(const mx::pproj_params& pparam, v3_t pt)
{
	v2_t p2;
    Project(pparam.mot, pt.p, p2.p);

    return p2;
}

void bundle::PSFMDriver::RefineCameraParameters(v3_t* points, v2_t* projs, int num_points, mx::pproj_params* pparams)			//
{
    double* motstruct = new double[cnp + num_points*pnp];
    double* nz_pts = motstruct+cnp;

    char* vmask = new char[num_points];
    double* projections = new double[mnp*num_points];

    memset(vmask, 0, sizeof(char)*num_points);

    double* proj = projections;
    double* z_3dpt = nz_pts;
    for(int i = 0; i < num_points; i++) {
        memcpy(proj, projs[i].p, sizeof(double)*2);
        proj += 2;
        memcpy(z_3dpt, points[i].p, sizeof(double)*3);
        z_3dpt +=3;
        vmask[i] = 1;
    }

    mx::MotCopyFormPProjParams(motstruct, *pparams);

    EX_BEGIN_CLOCK()
    Run(motstruct, projections, vmask, num_points, 1, 0, num_points);
    EX_END_CLOCK()
    EX_TRACE("# CameraRefine using %d pts, %d frames(%d fixed) and %d image projections, error %g [initial %g](elapse: %ld)\n",
             num_points, 1, 0, num_points, info[1]/num_points, info[0]/num_points, EX_ELAPSE());

    delete [] motstruct;
    delete [] vmask;
    delete [] projections;

	mx::PProjPCopyFormMotParams(pparams, motstruct);

    return;
}

int bundle::PSFMDriver::Run(double* motstruct, double* imgpts, char* vmask,
                           const int n3Dpts, const int ncams, const int nconcam, const int ncon3Dpts)
{
    int nvars = ncams*cnp+n3Dpts*pnp;
    int n = sba_motstr_levmar_x(n3Dpts, ncon3Dpts, ncams, nconcam, vmask, motstruct,
                                cnp, pnp, imgpts, NULL, mnp, img_ParallelProj_x,
                                img_ParallelProj_jac_x, &glob, MAXITER, verbose, opts, info);


    return n;
}

template<int PoseBlockSize>
class ReprojectionError: public ceres::SizedCostFunction<2, PoseBlockSize, 3>
{
public:
    ReprojectionError(double observation_x,double observation_y):
            _observation_x(observation_x),
            _observation_y(observation_y){}

    virtual bool Evaluate(double const* const* parameters,
                          double* residuals,
                          double** jacobians) const;


private:
    double _observation_x;
    double _observation_y;

};

template<>
bool ReprojectionError<6>::Evaluate(const double * const *parameters, double *residuals, double **jacobians) const
{
    util::point2d p;
    mx::pproj_params cam;
    cam.s=parameters[0][0];
    cam.alpha=parameters[0][1];
    cam.beta=parameters[0][2];
    cam.gamma=parameters[0][3];
    cam.t0=parameters[0][4];
    cam.t1=parameters[0][5];

    double t0 = cam.t0;
    double t1 = cam.t1;

    Eigen::Map<const Eigen::Vector3d> point(parameters[1]);

    double X, Y, Z;
    X=parameters[1][0];
    Y=parameters[1][1];
    Z=parameters[1][2];
    
    double point3d[3], pr[2];
    point3d[0]=X;point3d[1]=Y;point3d[2]=Z;
    bundle::PSFMDriver::Project(cam.mot, point3d, pr);

//     p=cam.Project(point);

// //        std::cout<<"X: "<<X<<", Y:"<<Y<<", Z:"<<Z<<std::endl;
//    std::cout<<"x: "<<_observation_x<<", y:"<<_observation_y<<std::endl;
//    std::cout<<"p_x: "<<pr[0]<<", p_y:"<<pr[1]<<std::endl;

    residuals[0] = pr[0]  - _observation_x;
    residuals[1] = pr[1]  - _observation_y;

//    std::cout<<"obs_x:"<<residuals[0]<<",obs_y:"<<residuals[1]<<std::endl;

    double cos_alpha = cos(cam.alpha);
    double sin_alpha = sin(cam.alpha);
    double cos_beta = cos(cam.beta);
    double sin_beta = sin(cam.beta);
    double cos_gamma = cos(cam.gamma);
    double sin_gamma = sin(cam.gamma);
    double s_s = 1/cam.s;

    if(jacobians !=NULL)
    {
        if(jacobians[0] != NULL)
        {
//            Eigen::Map<Eigen::Matrix<double, 2, 6, Eigen::RowMajor> > J_se3(jacobians[0]);
            Eigen::Map<Eigen::Matrix<double, 2, 6, Eigen::RowMajor>> J(jacobians[0]);
            J(0,0)=(-cos_gamma*(cos_beta*X+sin_alpha*sin_beta*Y-cos_alpha*sin_beta*Z)+sin_gamma*(cos_alpha*Y+sin_alpha*Z))*s_s*s_s;;
            J(0,1) = (cos_gamma*(cos_alpha*sin_beta*Y+sin_alpha*sin_beta*Z)-sin_gamma*(-sin_alpha*Y+cos_alpha*Z))*s_s;
            J(0,2) = cos_gamma*(-sin_beta*X+sin_alpha*cos_beta*Y-cos_alpha*cos_beta*Z)*s_s;
            J(0,3) = -sin_gamma*((cos_beta*X+sin_alpha*sin_beta*Y-cos_alpha*sin_beta*Z)*s_s-t0)-cos_gamma*((cos_alpha*Y+sin_alpha*Z)*s_s-t1);
            J(0,4) = -cos_gamma;
            J(0,5) = sin_gamma;
            J(1,0) = (-sin_gamma*(cos_beta*X+sin_alpha*sin_beta*Y-cos_alpha*sin_beta*Z)-cos_gamma*(cos_alpha*Y+sin_alpha*Z))*s_s*s_s;
            J(1,1) = (sin_gamma*(cos_alpha*sin_beta*Y+sin_alpha*sin_beta*Z)+cos_gamma*(-sin_alpha*Y+cos_alpha*Z))*s_s;
            J(1,2) = sin_gamma*(-sin_beta*X+sin_alpha*cos_beta*Y-cos_alpha*cos_beta*Z)*s_s;
            J(1,3) = cos_gamma*((cos_beta*X+sin_alpha*sin_beta*Y-cos_alpha*sin_beta*Z)*s_s-t0)-sin_gamma*((cos_alpha*Y+sin_alpha*Z)*s_s-t1);
            J(1,4) = -sin_gamma;
            J(1,5) = -cos_gamma;


        }


        if(jacobians[1] != NULL)
        {
            jacobians[1][0] = cos_gamma*cos_beta*s_s;
            jacobians[1][1] = (cos_gamma*sin_alpha*sin_beta-sin_gamma*cos_alpha)*s_s;
            jacobians[1][2] = (-cos_gamma*cos_alpha*sin_beta-sin_gamma*sin_alpha)*s_s;
            jacobians[1][3] = sin_gamma*cos_beta*s_s;
            jacobians[1][4] = (sin_gamma*sin_alpha*sin_beta+cos_gamma*cos_alpha)*s_s;
            jacobians[1][5] = (-sin_gamma*cos_alpha*sin_beta+cos_gamma*sin_alpha)*s_s;
        }

    }


    return true;

}

class MADN_loss : public ceres::LossFunction {
public:
    explicit MADN_loss(double delta) : delta_(delta) {}

    virtual void Evaluate(double s, double rho[3]) const {
//        std::cout<<"delta_:"<<delta_<<std::endl;
        rho[0]=s*delta_;
        rho[1] = delta_;       // first derivative of loss
        rho[2] = 0.0;       // second derivative of loss
    }

private:
    const double delta_;
};

double func_weight(double dij)
{
    if(dij < 0) return 1;
    else if(dij >1) return 0;
    else return (1 - dij *dij) * (1 - dij *dij);
}

double get_median_tmp(std::vector<double> mjlist)
{
    double median = 0.0;
    std::vector<double> tmp;
    for(int i=0;i<mjlist.size();i++)
    {
        if(mjlist[i]!=0) tmp.push_back(mjlist[i]);
    }
    std::sort(tmp.begin(), tmp.end());
//    for(int i = 0;i<mjlist.size();i++)
//    {
//        std::cout << mjlist[i]<<" ";
//    }
//    std::cout << std::endl;
    int n=tmp.size();
    if(tmp.size()%2!=0) median = tmp[tmp.size()*0.5];
    else {
        median = (tmp[n / 2] + tmp[n / 2 - 1]) / 2;
    }
    return median;
//    std::cout<<"median:" << median;
//    std::cout << std::endl;
}

double get_median(std::vector<double> mjlist)
{
    double median = 0.0;
    int n=mjlist.size();
    std::sort(mjlist.begin(), mjlist.end());
//    for(int i = 0;i<mjlist.size();i++)
//    {
//        std::cout << mjlist[i]<<" ";
//    }
//    std::cout << std::endl;
    if(mjlist.size()%2!=0) median = mjlist[mjlist.size()*0.5];
    else {
        median = (mjlist[n / 2] + mjlist[n / 2 - 1]) / 2;
    }
    return median;

}

std::vector<std::vector<double>> MAND(double *motstr, double *imgpts, char *vmask, int ncams, int n3Dpts){

    const int cmp=6;
    mx::pproj_params* camparams = new mx::pproj_params[ncams];
    double *Point3D=new double[n3Dpts*3];

    std::vector<double> ei_list;  // 储存每个点对应的ei(残差平方和列表的平均值)
//    double *mp_list = new double[n3Dpts*ncams]; // 储存每个点对应链的mp(残差值中位数)
    std::vector<double> mp_list; // 储存每个点对应链的mp(残差值中位数)
    std::vector<double> ma_list;  // 储存所有点的ei/mp的中位数
    double *MADN_list = new double[n3Dpts*ncams];  // 储存每个点对应的MADN
    double *di_list = new double[n3Dpts*ncams];  // 储存每个点对应的di
//    double *e_squares_list = new double[n3Dpts*ncams]; // 储存平方和列表,一行表示一个点对应链的平方和
    std::vector<std::vector<double> > eimp_list;
    std::vector<std::vector<double> > e_squares_list;

    for(int i=0;i<n3Dpts;i++)
    {
        Point3D[3*i]=motstr[ncams*6+3*i];
        Point3D[3*i+1]=motstr[ncams*6+3*i+1];
        Point3D[3*i+2]=motstr[ncams*6+3*i+2];
    }

    for(int i=0;i<ncams;i++)
    {
        mx::pproj_params cam_tmp;
        cam_tmp.s=motstr[i*6];
        cam_tmp.alpha=motstr[i*6+1];
        cam_tmp.beta=motstr[i*6+2];
        cam_tmp.gamma=motstr[i*6+3];
        cam_tmp.t0=motstr[i*6+4];
        cam_tmp.t1=motstr[i*6+5];
        camparams[i]=cam_tmp;
    }

    std::vector<double> f, u;
    int gap = 0;
    double e = 0;
    for(int i=0;i<n3Dpts;i++)
    {
        util::point2d p;
        Eigen::Vector3d point_tmp=Eigen::Vector3d(Point3D[3*i],Point3D[3*i+1],Point3D[3*i+2]);
        double X[3], pr[2];
        X[0]=point_tmp[0]; X[1]=point_tmp[1]; X[2]=point_tmp[2];
//        std::cout<<"Point"<<i<<std::endl;
        std::vector<double> residual_list_ipoints;
        for(int j=0;j<ncams;j++)
        {
            double x_o,y_o;
            mx::pproj_params cam_tmp = camparams[j];
            bundle::PSFMDriver::Project(cam_tmp.mot, X, pr);
//             p = cam_tmp.Project(point_tmp);
//            std::cout<<"x:"<<p.x<<" y:"<<p.y<<std::endl;
            if(vmask[i*ncams+j]){
                x_o=imgpts[i*(ncams*2)+2*j+0-gap*2];
                y_o=imgpts[i*(ncams*2)+2*j+1-gap*2];
                u.push_back(x_o);
                u.push_back(y_o);

                f.push_back(pr[0]);
                f.push_back(pr[1]);

                e=(pr[0]-x_o)*(pr[0]-x_o)+(pr[1]-y_o)*(pr[1]-y_o);

                residual_list_ipoints.push_back(e);
//                e_squares_list[i]
            }
            else{
                f.push_back(0.0);
                f.push_back(0.0);
                u.push_back(0.0);
                u.push_back(0.0);
                gap++;

                residual_list_ipoints.push_back(0);
            }
        }
        e_squares_list.push_back(residual_list_ipoints);
    }

    for(int i=0;i<e_squares_list[0].size();i++){
        double median=0;
        std::vector<double> e_tmp;
        for(int j=0;j<e_squares_list.size();j++){
            e_tmp.push_back(e_squares_list[j][i]);
        }
//        std::cout<<"size:"<<i<<" "<<e_tmp.size()<<std::endl;
//        std::cout<<"e_squares_list"<<i<<std::endl;
        median=get_median(e_tmp);
////        std::cout<<std::endl;
//        std::cout<<"median:"<<median<<std::endl;
        mp_list.push_back(median);
    }

    for(int i=0;i<n3Dpts;i++){
        double ma = 0.0;
        double median = 0.0;
        std::vector<double> eimp_i;
        for(int j=0;j<ncams;j++){
            if(mp_list[j])
            {
                ma=e_squares_list[i][j]/mp_list[j];
                eimp_i.push_back(ma);
            } else eimp_i.push_back(0);
        }
        median=get_median(eimp_i);
        ma_list.push_back(median);
        eimp_list.push_back(eimp_i);
    }

    for(int i=0;i<n3Dpts;i++){
        for(int j=0;j<ncams;j++){
            ei_list.push_back(e_squares_list[i][j]);
        }
    }
    
    std::vector<double> E_MADN_list;
    std::vector<std::vector<double> > ej_list;
    std::cout<<"ei_list:"<<e_squares_list.size()<<std::endl;
    for(int i=0;i<ncams;i++)
    {
        std::vector<double> ej;
        for(int j=0;j<n3Dpts;j++)
        {
            ej.push_back(e_squares_list[j][i]);
        }
        ej_list.push_back(ej);
    }

//    std::cout<<"ej_list:"<<ej_list.size()<<std::endl;

    for(int i=0;i<ej_list.size();i++)
    {
        double med = get_median_tmp(ej_list[i]);
        std::vector<double> abs_x_med_j;
        for(int j=0;j<ej_list[i].size();j++)
        {
            if(ej_list[i][j]!=0) abs_x_med_j.push_back(abs(ej_list[i][j]-med));
        }
        double MADN_j = get_median(abs_x_med_j) / 0.6745;
        E_MADN_list.push_back(MADN_j);
    }


    std::vector<std::vector<double>> W_MADN_j;
std::vector<double> W_MADN_row;
double W_MADN_tmp[ncams][n3Dpts];
    for(int i=0;i<E_MADN_list.size();i++)
    {
        std::vector<double> W_i;
        double dij, wij;
        for(int j=0;j<ej_list[i].size();j++)
        {
            if(ej_list[i][j]!=0)
            {
                dij = (ej_list[i][j]) / (4.685 * E_MADN_list[i]);
                wij = func_weight(dij);
                W_MADN_tmp[i][j]=wij;
            }
            else {
                W_MADN_tmp[i][j]=0;
            }
            dij = (ej_list[i][j]) / (4.685 * E_MADN_list[i]);
            wij = func_weight(dij);
//             W_MADN_row.push_back(wij);
//             W_MADN_tmp[i][j]=wij;
//             if(wij<0) std::cout<<"wij:"<<wij<<"*******************"<<std::endl;
            W_i.push_back(wij);
        }
        W_MADN_j.push_back(W_i);
    }
    
    for(int i=0;i<n3Dpts;i++)
    {
        for(int j=0;j<ncams;j++)
        {
//             std::cout<<W_MADN_tmp[j][i]<<" ";
        }
//         std::cout<<std::endl;
    }
    std::vector<std::vector<double>> W_MADN;
    std::vector<double > W_tmp;

    for(int i=0;i<n3Dpts;i++)
    {
        std::vector<double > W_tmp;
        for(int j=0;j<ncams;j++)
        {
            W_tmp.push_back(W_MADN_tmp[j][i]);
//             if(W_MADN_j[j][i]<0) std::cout<<"*******************"<<W_MADN_j[j][i]<<std::endl;
        }
        for(int j=0;j<W_tmp.size();j++)
        {
//             std::cout<<W_tmp[j]<<" ";
        }
//         std::cout<<std::endl;
        W_MADN.push_back(W_tmp);
    }
    int num_outliers1=0;

    double med_1 = get_median(ei_list);
    std::vector<double> abs_x_med;
    for(int i=0;i<ei_list.size();i++){
        abs_x_med.push_back(abs(ei_list[i]-med_1));
    }

    double MADN = get_median(abs_x_med) / 0.6745;



    std::vector<std::vector<double>> W;

    for(int i=0;i<n3Dpts;i++)
    {
        std::vector<double> W_i;
        double dij, wij;
        for(int j=0;j<ncams;j++)
        {
            dij = (eimp_list[i][j] - ma_list[i]) / (4.685 * MADN);
            wij = func_weight(dij);
            W_i.push_back(wij);
        }
        W.push_back(W_i);
    }
int num_outliers=0;
// std::cout<<"W:"<<std::endl;
   for(int i=0;i<n3Dpts;i++)
   {
//        std::cout<<"Point"<<i<<std::endl;
       for(int j=0;j<ncams;j++){
//            std::cout<<W_MADN[i][j]<<" ";
           if(W_MADN[i][j]>0.1) num_outliers++;
           
       }
//        std::cout<<std::endl;
   }
//    std::cout<<"num_outliers:"<<num_outliers<<std::endl;
//     return W;
    return W_MADN;

}


int bundle::PSFMDriver::Run_ceres(double* motstruct, double* imgpts, char* vmask, const int n3Dpts, const int ncams, const int nconcam, const int ncon3Dpts, std::vector<std::vector<double> >& W_init, bool use_MAND)
{
    int arr_num=0, arr_vmask=0;
    for(int i=0;i<n3Dpts*ncams;i++)
    {
        if(vmask[i]) arr_vmask++;
    }
    std::cout<<"arr_vmask:"<<arr_vmask<<std::endl;
    int n=0;
    mx::pproj_params* camparams = new mx::pproj_params[ncams];
    double *Point3D=new double[n3Dpts*3];

    for(int i=0;i<n3Dpts;i++)
    {
        Point3D[3*i]=motstruct[ncams*6+3*i];
        Point3D[3*i+1]=motstruct[ncams*6+3*i+1];
        Point3D[3*i+2]=motstruct[ncams*6+3*i+2];
    }

    PosePointParametersBlock states;
//    PosePointParametersBlock init_states;
    states.create(ncams, n3Dpts);
    for(int i=0;i<ncams;i++)
    {
        double *cameras_true (states.pose(i));
//        double *cameras=A;
        cameras_true[0]=motstruct[i*6];
        cameras_true[1]=motstruct[i*6+1];
        cameras_true[2]=motstruct[i*6+2];
        cameras_true[3]=motstruct[i*6+3];
        cameras_true[4]=motstruct[i*6+4];
        cameras_true[5]=motstruct[i*6+5];

    }

//    for(int i=0;i<ncams;i++)
//    {
//        std::cout<<states.pose(i)[0]<<" "<<states.pose(i)[1]<<std::endl;
//    }

    for(int i=0;i<n3Dpts;i++)
    {
        Eigen::Map<Eigen::Vector3d> true_pt(states.point(i));  //类似于引用地址
//        Eigen::Map<Eigen::Vector3d> init_pt(init_states.point(i));
        true_pt = Eigen::Vector3d(motstruct[ncams*6+3*i],
                                  motstruct[ncams*6+3*i+1],
                                  motstruct[ncams*6+3*i+2]);
//        init_pt = Eigen::Vector3d(motstr[ncams*6+3*i],
//                                  motstr[ncams*6+3*i+1],
//                                  motstr[ncams*6+3*i+2]);

    }


    for(int i=0;i<ncams;i++)
    {
        mx::pproj_params cam_tmp;
        cam_tmp.s=motstruct[i*6];
        cam_tmp.alpha=motstruct[i*6+1];
        cam_tmp.beta=motstruct[i*6+2];
        cam_tmp.gamma=motstruct[i*6+3];
        cam_tmp.t0=motstruct[i*6+4];
        cam_tmp.t1=motstruct[i*6+5];
        camparams[i]=cam_tmp;
//        Camera.cam.push_back(cam_tmp);
    }

    std::vector<double> f, u;
    int gap = 0;
    std::vector<util::point2d> img_2d;
    for(int i=0;i<n3Dpts;i++)
    {
        double X[3], pr[2];
        util::point2d p;
        Eigen::Vector3d point_tmp=Eigen::Vector3d(Point3D[3*i],Point3D[3*i+1],Point3D[3*i+2]);
        X[0]=Point3D[3*i]; X[1]=Point3D[3*i+1]; X[2]=Point3D[3*i+2];
//        std::cout<<"Point"<<i<<std::endl;
        for(int j=0;j<ncams;j++)
        {
            double x_o,y_o;
            mx::pproj_params cam_tmp = camparams[j];
            Project(cam_tmp.mot, X, pr);
//             p = cam_tmp.Project(point_tmp);
//            std::cout<<"x:"<<p.x<<" y:"<<p.y<<std::endl;
            if(vmask[i*ncams+j]){
                x_o=imgpts[i*(ncams*2)+2*j+0-gap*2];
                y_o=imgpts[i*(ncams*2)+2*j+1-gap*2];

                util::point2d point_tmp;
                point_tmp.x=x_o;
                point_tmp.y=y_o;
                img_2d.push_back(point_tmp);

//                std::cout<<"x_o:"<<x_o<<" y_o:"<<y_o<<std::endl;
                u.push_back(x_o);
                u.push_back(y_o);

                f.push_back(pr[0]);
                f.push_back(pr[1]);
            }
            else{
                f.push_back(0.0);
                f.push_back(0.0);
                u.push_back(0.0);
                u.push_back(0.0);
                gap++;
            }
        }
    }


    double dij=0.0;
    std::vector<std::vector<double>> W;
//     bool use_MAND=true;
    if(use_MAND){
        W_init=MAND(motstruct, imgpts, vmask, ncams, n3Dpts);
    }
    else{
        for(int i=0;i<n3Dpts;i++)
        {
            std::vector<double> W_i;
            for(int j=0;j<ncams;j++)
            {
               W_i.push_back(1);
            }
            W_init.push_back(W_i);
        }
    }
    
    ceres::Problem problem;
    gap=0;
    double x_o,y_o;
    int sum_tmp=0;
    for(int i=0; i<n3Dpts; i++){
        for(int j=0;j<ncams; j++){
            dij=W_init[i][j];
            if(vmask[i*ncams+j]){
                sum_tmp++;

                x_o=imgpts[i*(ncams*2)+2*j+0-gap*2];
                y_o=imgpts[i*(ncams*2)+2*j+1-gap*2];

//                std::cout<<"x_o"<<x_o<<" y_o"<<y_o<<std::endl;

//                std::cout<<"Point:"<<i<<" "<<j<<std::endl;
//                std::cout<<"x:"<<point3d[0]<<" y:"<<point3d[1]<<" z:"<<point3d[2]<<std::endl;
                ceres::CostFunction* cost_function;
                cost_function=new ReprojectionError<6>(x_o, y_o);
//                ceres::Problem::AddResidualBlock()
//std::cout<<"dij:"<<dij<<std::endl;
// if(dij!=1) dij=1;
                ceres::LossFunction* lossFunc = new MADN_loss(dij);
//                ceres::LossFunction* lossFunc = new ceres::HuberLoss(0.5);
//                problem.AddResidualBlock(cost_function, NULL, states.pose(j), states.point(i));
                problem.AddResidualBlock(cost_function,lossFunc, states.pose(j), states.point(i));
//                problem.AddResidualBlock(cost_function,new ceres::CauchyLoss(1), states.pose(j), states.point(i));
            }
            else{
                gap++;
            }


        }
    }
    std::cout<<"sum_tmp:"<<sum_tmp<<std::endl;
    double init_residual=0;
    for (int i = 0; i < u.size(); i=i+1) {
        double delt;
        delt=sqrt((u[i]-f[i])*(u[i]-f[i]));
//        delt=abs(u[i]-f[i]);
        init_residual=init_residual+delt;
    }

    std::cout<<"init_ResErr:"<<init_residual/(arr_vmask)<<std::endl;
    
    ceres::Solver::Options options;
    //配置增量方程的解法
//    options.linear_solver_type = ceres::SPARSE_SCHUR;
    options.trust_region_strategy_type = ceres::LEVENBERG_MARQUARDT;
    options.minimizer_progress_to_stdout = true;
    options.max_num_iterations = 1000;

    //第三步，创建Summary对象用于输出迭代结果
    ceres::Solver::Summary summary;

//    ceres::Solver(options,&problem,&summary);
    ceres::Solve(options, &problem, &summary);

    std::cout << summary.BriefReport() << "\n";
    
    
    gap=0;
    f.clear();
double abs_re=0;
    std::vector<std::vector<int>> W_mask;
    std::vector<double> residual_e;

    std::vector<std::vector<int > > W_ord;
    for(int i=0;i<n3Dpts;i++)
    {
        std::vector<int > W_ord_tmp;
        for(int j=0;j<ncams;j++)
        {
            W_ord_tmp.push_back(0);
        }
        W_ord.push_back(W_ord_tmp);
    }

    for(int i=0;i<n3Dpts;i++)
    {
        std::vector<int> mask_tmp;
        util::point2d p;
        Eigen::Vector3d point_tmp=Eigen::Vector3d(states.point(i)[0],states.point(i)[1],states.point(i)[2]);
//        std::cout<<"Point"<<i<<"    "<<point_tmp[0]<<"    "<<point_tmp[1]<<"    "<<point_tmp[2]<<std::endl;
        double X[3], pr[2];
        X[0]=point_tmp[0]; X[1]=point_tmp[1]; X[2]=point_tmp[2];
        for(int j=0;j<ncams;j++)
        {
            double x_o,y_o;
            mx::pproj_params cam_tmp ;
            cam_tmp.s= states.pose(j)[0];
            cam_tmp.alpha= states.pose(j)[1];
            cam_tmp.beta= states.pose(j)[2];
            cam_tmp.gamma= states.pose(j)[3];
            cam_tmp.t0= states.pose(j)[4];
            cam_tmp.t1= states.pose(j)[5];
//             Project(cam_tmp.mot, X, pr);
//             p = cam_tmp.Project(point_tmp);
//            std::cout<<"x:"<<pr[0]<<" y:"<<pr[1]<<std::endl;
            if(vmask[i*ncams+j]){
                if(W_init[i][j]>0.9) W_ord[i][j]=1;
                else W_ord[i][j]=2;
                if(W_init[i][j]>0)
                {
                    Project(cam_tmp.mot, X, pr);
//                     std::cout<<"x:"<<pr[0]<<" y:"<<pr[1]<<std::endl;
                    mask_tmp.push_back(1);
                    arr_num++;
                    x_o=imgpts[i*(ncams*2)+2*j+0-gap*2];
                    y_o=imgpts[i*(ncams*2)+2*j+1-gap*2];
                                        abs_re=pr[0]-x_o+pr[1]-y_o+abs_re;
//                std::cout<<"x_o:"<<x_o<<" y_o:"<<y_o<<std::endl;
                    f.push_back(pr[0]-x_o);
                    f.push_back(pr[1]-y_o);

                    double dx, dy;
                    dx=pr[0]-x_o;
                    dy=pr[1]-y_o;
// std::cout<<"dx:"<<dx<<" dy:"<<dy<<std::endl;
                    residual_e.push_back(sqrt(dx*dx+dy*dy));
                }
            }
            else{
                mask_tmp.push_back(0);
                    f.push_back(0.0);
                    f.push_back(0.0);
                    gap++;
            }
        }
        W_mask.push_back(mask_tmp);
    }
//     std::cout<<"***********************abs_re"<<abs_re<<std::endl;
// int oooo=0;
//     for(int i=0;i<W_mask.size();i++)
//     {
//         for(int j=0;j<W_mask[i].size();j++)
//         {
//             std::cout<<W_mask[i][j]<<"   ";
//             if(W_mask[i][j]>0) oooo++;
//         }
//         std::cout<<std::endl;
//         std::cout<<oooo<<std::endl;
//     }

    double residual=0;
//     std::cout<<f.size()<<std::endl;
    for (int i = 0; i < f.size(); i=i+1) {
        double delt;
        delt=sqrt(f[i]*f[i]);
        residual+=delt;
    }
    double avg_res=0, std_sum=0, std=0;

    double residual_tmp=0;
    for (int i = 0; i < residual_e.size(); i=i+1) {
        double delt;
        delt=residual_e[i];
        residual_tmp=residual_tmp+delt;
    }

    avg_res=residual_tmp/arr_num;
    for (int i = 0; i < residual_e.size(); i=i+1) {
        double delt=0;
        delt=residual_e[i]-avg_res;
        std_sum=delt*delt+std_sum;
    }

    std = sqrt(std_sum/arr_num);

    std::cout<<std::endl;
//     std::cout<<"ResErr:"<<residual/(arr_num)<<std::endl;
    std::cout<<"ResErr_new:"<<residual_tmp/(arr_num)<<"   std:"<<std<<std::endl;
    std::cout<<"residual:"<<residual<<"        arr_num:"<<arr_num<<std::endl;
    
    
//     std::cout<<"***in***"<<std::endl;
    for(int j=0;j<ncams;j++) {
       motstruct[j*6] = states.pose(j)[0];
       motstruct[j*6+1] = states.pose(j)[1];
       motstruct[j*6+2] = states.pose(j)[2];
       motstruct[j*6+3] = states.pose(j)[3];
       motstruct[j*6+4] = states.pose(j)[4];
       motstruct[j*6+5] = states.pose(j)[5];

//        std::cout<<states.pose(j)[0]<<"      "<<states.pose(j)[1]<<"      "<<states.pose(j)[2]<<"      "<<states.pose(j)[3]<<"      "<<states.pose(j)[4]
//                     <<"      "<<states.pose(j)[5]<<std::endl;

   }
   for(int j=0;j<n3Dpts;j++) {
       motstruct[ncams*6+j*3] = states.point(j)[0];
       motstruct[ncams*6+j*3+1] = states.point(j)[1];
       motstruct[ncams*6+j*3+2] = states.point(j)[2];
//               std::cout<<states.point(j)[0]<<"      "<<states.point(j)[1]<<"      "<<states.point(j)[2]<<std::endl;
   }

    
    
    return n;
}


// void bundle::PSFMDriver::data_initial(int ncams,int n3Dpts,double ratio_noise,double ratio_outlier,double ratio_noise_cams,double range_img){
void bundle::PSFMDriver::data_initial(double* motstruct, double* imgpts_input, char* vmask_input, double* vec, int ncams,int n3Dpts,double ratio_noise,double ratio_outlier,double ratio_noise_cams,double range_img, double& mean,double& std){
    int i,j,n2dpts,ratio_miss_pts,num_outer,*index_outer;
    double *cams,*pts3d,*imgpts,*cams_noise,*real_2dpoints,*pts2d_noise,*pts2d_noise_outlier,*p,*pts3d_simulate;
    char *vmask;
    std::ifstream fin;
    std::random_device rd;
    std::default_random_engine rng {rd()};

    cams = new double[ncams*6];
    pts3d = new double[n3Dpts*3];
    vmask = new char[n3Dpts*ncams];


    for(int i=0;i<ncams*6;i++){
        cams[i] = motstruct[i];
        // std::cout<<cams[i]<<std::endl;
    }

    for(int i = 0; i < n3Dpts * 3; i++) {
        pts3d[i] = motstruct[ncams * 6 + i];
        // std::cout << pts3d[i] << std::endl;
    }
    n2dpts =0;
    for(int i=0;i<n3Dpts*ncams;i++){
        vmask[i] = vmask_input[i];
        // vmask[i] -= '0';
        if(!((int)vmask[i] == 0)){ n2dpts++;}
    }

    n2dpts *=2;
    pts2d_noise_outlier = new double[n2dpts];

    for (int i = 0; i < n2dpts; i++) {
        pts2d_noise_outlier[i] = imgpts_input[i];
    }

    double *weight_one = new double[ncams*n3Dpts] ;
    for(i=0;i<ncams*n3Dpts;i++) weight_one[i] = 1.0;



    //将模拟数据g归一化并存储到文件夹中
    p = new double[ncams*6+n3Dpts*3];
    double one_avg_random,one_std_random;
    pts3d_simulate = new double[n3Dpts*3];
    for (i = 0; i < n3Dpts*3; i++){
        pts3d_simulate[i] = pts3d[i];
    }
    for(i=0;i<ncams*6;i++) p[i] = cams[i];
    for(i=ncams*6;i<ncams*6+n3Dpts*3;i++) p[i] = pts3d_simulate[i-ncams*6];
    write_txt(p,ncams*6+n3Dpts*3,"simulate_data_noOne/mot_random.txt");
    write_txt(pts2d_noise_outlier,n2dpts,"simulate_data_noOne/pts2d_random.txt");
    write_txt(vmask,ncams*n3Dpts,"simulate_data_noOne/vmask_random.txt");


    Toone_new(&one_avg_random,&one_std_random,p,pts2d_noise_outlier,vmask,ncams,n3Dpts);
    double *one_avg = new double[1];
    double *one_std = new double[1];
    double *vec_gamma = new double[ncams];
    for(i=0;i<ncams;i++) vec_gamma[i] = p[i*6+3];
    one_avg[0] = one_avg_random;
    one_std[0] = one_std_random;
    write_txt_app<double>(one_avg,1,"simulate_data_noOne/one_avg.txt");
    write_txt_app<double>(one_std,1,"simulate_data_noOne/one_std.txt");

    mean = one_avg[0];
    std = one_std[0];


    write_txt(vec_gamma,ncams,"simulate_data_noOne/vec_gamma.txt");
    std::memcpy(vec, vec_gamma, (ncams) * sizeof(double));

    int *ncams_vec,*n3Dpts_vec;
    ncams_vec = new int[1];
    n3Dpts_vec = new int[1];
    ncams_vec[0] = ncams;
    n3Dpts_vec[0] = n3Dpts;
    write_txt_app<int>(ncams_vec,1,"output_temp/ncams_vec.txt");
    write_txt_app<int>(n3Dpts_vec,1,"output_temp/n3dpts_vec.txt");

    double *error_begin = new double[1];
    error_begin[0] = f(p,0,pts2d_noise_outlier,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2);
    std::cout<<"n2Dpts is "<<n2dpts<<std::endl;
    std::cout<<"initial error is "<<f(p,0,pts2d_noise_outlier,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2)<<std::endl;
    /*for(i=0;i<ncams*6+n3Dpts*3;i++) std::cout<<p[i]<<" ";
    std::cout<<std::endl;
    for(i=0;i<n2dpts;i++) std::cout<<pts2d_noise_outlier[i]<<" ";
    std::cout<<std::endl;*/
    write_txt_app<double>(error_begin,1,"output_temp/error_begin.txt");
    write_txt(p,ncams*6+n3Dpts*3,"simulate_data/mot_random.txt");
    write_txt(pts2d_noise_outlier,n2dpts,"simulate_data/pts2d_random.txt");
    write_txt(vmask,ncams*n3Dpts,"simulate_data/vmask_random.txt");



    double *mot = new double[ncams*6+n3Dpts*3];
    double *mot_temp = new double[ncams*6+n3Dpts*3];
    fin.open("simulate_data/mot_random.txt", std::ios::in);
        for (i = 0; i < ncams*6+n3Dpts*3; i++) fin >> mot[i];
        fin.close();
    double_copy(mot_temp,mot,ncams*6+n3Dpts*3);
    std::cout<<"initial error is "<<f(mot_temp,0,pts2d_noise_outlier,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2)<<std::endl;




    // delete[] cams;
    // delete[] pts3d;
    // delete[] vmask;
    // delete[] pts2d_noise_outlier;
    // delete[] p;
    // delete[] pts3d_simulate;
    // delete[] weight_one;
    // delete[] ncams_vec;
    // delete[] n3Dpts_vec;
    // delete[] vec_gamma;
}

// void bundle::PSFMDriver::data_initial(double* motstruct, double* imgpts_input, char* vmask_input, int ncams,int n3Dpts,double ratio_noise,double ratio_outlier,double ratio_noise_cams,double range_img){
//     int n2dpts,ratio_miss_pts,num_outer,*index_outer;
//     double *cams,*pts3d,*p,*pts3d_simulate;
//     char *vmask;
//     std::ifstream fin;
//     // std::random_device rd;
//     // std::default_random_engine rng {rd()};

//     cams = new double[ncams*6];
//     pts3d = new double[n3Dpts*3];
//     vmask = new char[n3Dpts*ncams];

//     for(int i=0;i<ncams*6;i++){
//         cams[i] = motstruct[i];
//         // std::cout<<"*****************"<<cams[i] <<std::endl;
//     }

//     for(int i=ncams*6;i<ncams*6+n3Dpts*3;i++){
//         pts3d[i] = motstruct[i];
//         // std::cout<<"*****************"<<pts3d[i] <<std::endl;
//     }
//     n2dpts =0;
//     for(int i=0;i<n3Dpts*ncams;i++){
//         vmask[i] = vmask_input[i];
//         // vmask[i] -= '0';
//         if(!((int)vmask[i] == 0)){ n2dpts++;}
//         // std::cout<<"*****************"<<(int)vmask[i] <<std::endl;
//     }

   
//     double *weight_one = new double[ncams*n3Dpts] ;
//     for(int i=0;i<ncams*n3Dpts;i++) weight_one[i] = 1.0;
//     // n2dpts =0;
//     // fin.open("changedata_output/cams.txt", std::ios::in);
//     // for (i = 0; i < ncams*6; i++) fin >> cams[i];
//     // fin.close();
//     // fin.open("changedata_output/points3d.txt", std::ios::in);
//     // for (int i = 0; i < n3Dpts*3; i++) fin >> pts3d[i];
//     // fin.close();
//     // fin.open("changedata_output/vmask.txt", std::ios::in);
//     // for (int i = 0; i < n3Dpts*ncams; i++) {
//     //     fin >> vmask[i];
//     //     vmask[i] -= '0';
//     //     if(!((int)vmask[i] == 0)){ n2dpts++;}
//     // }
//     // fin.close();
//     n2dpts *=2;
//     double* pts2d_noise_outlier = new double[n2dpts];
//     for (int i = 0; i < n2dpts; i++) {
//         pts2d_noise_outlier[i] = imgpts_input[i];
//     }
//     // fin.open("changedata_output/points2d.txt", std::ios::in);
//     // for (int i = 0; i < n2dpts; i++) {
//     //     fin >> pts2d_noise_outlier[i];
//     // }
//     // fin.close();
//     // std::cout<<"*****************"<<std::endl;

    

//     //将模拟数据g归一化并存储到文件夹中
//     p = new double[ncams*6+n3Dpts*3];
//     double one_avg_random,one_std_random;
//     pts3d_simulate = new double[n3Dpts*3];
//     for (int i = 0; i < n3Dpts*3; i++){
//         pts3d_simulate[i] = pts3d[i]; 
//     }
//     for(int i=0;i<ncams*6;i++) p[i] = cams[i];
//     for(int i=ncams*6;i<ncams*6+n3Dpts*3;i++) p[i] = pts3d_simulate[i-ncams*6];
//     write_txt(p,ncams*6+n3Dpts*3,"simulate_data_noOne/mot_random.txt");
//     write_txt(pts2d_noise_outlier,n2dpts,"simulate_data_noOne/pts2d_random.txt");
//     write_txt(vmask,ncams*n3Dpts,"simulate_data_noOne/vmask_random.txt");


//     Toone_new(&one_avg_random,&one_std_random,p,pts2d_noise_outlier,vmask,ncams,n3Dpts);
//     double *one_avg = new double[1];
//     double *one_std = new double[1];
//     double *vec_gamma = new double[ncams];
//     for(int i=0;i<ncams;i++) vec_gamma[i] = p[i*6+3]; 
//     one_avg[0] = one_avg_random;
//     one_std[0] = one_std_random;
//     write_txt_app<double>(one_avg,1,"simulate_data_noOne/one_avg.txt");
//     write_txt_app<double>(one_std,1,"simulate_data_noOne/one_std.txt");
//     write_txt(p,ncams*6+n3Dpts*3,"simulate_data/mot_random.txt");
//     write_txt(pts2d_noise_outlier,n2dpts,"simulate_data/pts2d_random.txt");
//     write_txt(vmask,ncams*n3Dpts,"simulate_data/vmask_random.txt");
//     write_txt(vec_gamma,ncams,"simulate_data_noOne/vec_gamma.txt");

//     int *ncams_vec,*n3Dpts_vec;
//     ncams_vec = new int[1];
//     n3Dpts_vec = new int[1];
//     ncams_vec[0] = ncams;
//     n3Dpts_vec[0] = n3Dpts;
//     write_txt_app<int>(ncams_vec,1,"output_temp/ncams_vec.txt");
//     write_txt_app<int>(n3Dpts_vec,1,"output_temp/n3dpts_vec.txt");

//     double *error_begin = new double[1];
//     error_begin[0] = f(p,0,pts2d_noise_outlier,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2);
//     std::cout<<"initial error is "<<f(p,0,pts2d_noise_outlier,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2)<<std::endl;
//     write_txt_app<double>(error_begin,1,"output_temp/error_begin.txt");
    


//     delete[] cams;
//     delete[] pts3d;
//     delete[] vmask;
//     delete[] pts2d_noise_outlier;
//     delete[] p;
//     delete[] pts3d_simulate;
//     delete[] weight_one;
//     delete[] ncams_vec;
//     delete[] n3Dpts_vec;
//     delete[] vec_gamma;
// }



double bundle::PSFMDriver::compute_norm1(double *a,int len){
    int i;
    double output;
    output=0;
    for(i=0;i<len;i++){
        if(a[i]>0)output+=a[i];
        if(a[i]<0)output=output-a[i];
        //output+=abs(a[i]);
    }
    return output;
}

void bundle::PSFMDriver::compute_gradient(double *p,double *imgpts, double mu,char *vmask,double *weight, int ncams, int n3Dpts, double *grad) {
    /*********
    计算当前光滑近似函数的梯度
    Input：
    p为自变量，先是相机参数（s,omiga,beta,gamma,t0,t1），后是维点坐标（x，y，z）
    imgpts为点链，排列为第i个点在第j个相机下的坐标
    mu为光滑近似参数，vmask为掩码
    ncams，n3Dpts为相机个数与三维点个数
    Output：
    grad：当前近似函数的梯度
    ***********/
    int i,j,ncol;
    double *cameras,*points,*dfdR,*dfdX,du,dv;

    ncol = ncams*6+n3Dpts*3;
    memset(grad,0,ncol*sizeof(double));

    /*cameras = new double[ncams*6];
    points = new double[n3Dpts*3];
    memcpy(cameras, p,sizeof(double) * (ncams *6));
    memcpy(points, p + (ncams *6), sizeof(double) * (n3Dpts*3));*/

    dfdR = new double[6];
    memset(dfdR,0,6*sizeof(double));
    dfdX = new double[3];
    memset(dfdR,0,3*sizeof(double));

    for(i=0;i<ncams;i++){

        compute_dfdR(p,imgpts,vmask,weight,mu,ncams,n3Dpts,i,dfdR);
        memcpy(grad+i*6,dfdR,6*sizeof(double));
        /*****************
        //for(j=0;j<6;j++) grad[i*6+j]=dfdR[j];
        //for(j=0;j<6;j++) std::cout<<dfdR[j]<<" ";
        ******************/
    }

    for(i=0;i<n3Dpts;i++){
        compute_dfdX(p,imgpts,vmask,weight,mu,ncams,n3Dpts,i,dfdX);
        memcpy(grad+ncams*6+i*3,dfdX,3*sizeof(double));

    }
 
    //delete(cameras);
    //delete(points);
    delete(dfdR);
    delete(dfdX);
}

double bundle::PSFMDriver::f(double *p,double mu,double *imgpts, char *vmask,double *weight,int ncams,int n3Dpts){
    //计算光滑近似函数的值
    int i, j, index_uv = 0;
    double f_out = 0, u_fid, v_fid;

    std::vector<double> camera(6);
    std::vector<double> point3D(3);
    std::vector<double> uv(2);

    std::ofstream out;
    // out.open("pts2d_residual.txt");

    for (i = 0; i < n3Dpts; i++) {
        for (j = 0; j < ncams; j++) {
            if ((int)vmask[i * ncams + j] == 1) {
                std::copy(p + j * 6, p + j * 6 + 6, camera.begin());
                std::copy(p + 6 * ncams + i * 3, p + 6 * ncams + i * 3 + 3, point3D.begin());
                compute_uv(camera.data(), point3D.data(), uv.data());
                
                u_fid = imgpts[index_uv++];
                v_fid = imgpts[index_uv++];
                
                out << uv[0] - u_fid << std::endl;
                out << uv[1] - v_fid << std::endl;
                
                f_out += weight[i * ncams + j] * sqrt(pow(uv[0] - u_fid, 2) + mu * mu);
                f_out += weight[i * ncams + j] * sqrt(pow(uv[1] - v_fid, 2) + mu * mu);
            }
        }
    }
    // out.close();

    // return 0;

    return f_out;
}
void bundle::PSFMDriver::f(double *p,double mu,double *imgpts, char *vmask,double *weight,int ncams,int n3Dpts, std::string filename){
    //计算光滑近似函数的值
    int i,j,k,index_uv;
    double *camera,*point3D,*uv,f_out,u_fid,v_fid;

    camera = new double[6];
    point3D = new double[3];
    uv = new double[2];

    f_out =0;
    index_uv = 0;
    std::ofstream out;
    // out.open(filename);

    for(i=0;i<n3Dpts;i++){
        for(j=0;j<ncams;j++){
            if((int)vmask[i*ncams+j] == 1){
                memcpy(camera,p+j*6,6*sizeof(double));
                memcpy(point3D,p+6*ncams+i*3,3*sizeof(double));
                compute_uv(camera,point3D,uv);
                //std::cout<<" "<<uv[0]<<" "<<uv[1]<<std::endl;
                u_fid = imgpts[index_uv];
                index_uv++;
                v_fid = imgpts[index_uv];
                index_uv++;
                out<<uv[0]-u_fid<<std::endl;out<<uv[1]-v_fid<<std::endl;
                //f_out +=abs(uv[0]-u_fid)+abs(uv[1]-v_fid);
                f_out +=weight[i*ncams+j]*sqrt(pow(uv[0]-u_fid,2)+mu*mu);
                f_out +=weight[i*ncams+j]*sqrt(pow(uv[1]-v_fid,2)+mu*mu);
                //std::cout<<"f"<<f_out<<"u"<<uv[0]<<"v"<<uv[1]<<" ";
                //std::cout<<"f"<<f_out<<"f+"<<sqrt(pow(uv[0]-u_fid,2))<<" ";
            }
        }
    }
    // // out.close();
    // delete(camera);
    // delete(point3D);
    // delete(uv);

    //return f_out;
}

void bundle::PSFMDriver::compute_right_R(double *camera,double *point,double *right_R_u,double *right_R_v){
    //被compute_dfdR调用，因为公式太长了
    double s,alpha,beta,gamma,t0,t1,x,y,z,sin_alpha,cos_alpha,sin_beta,cos_beta,sin_gamma,cos_gamma,temp1,temp2;

    s = camera[0];
    alpha = camera[1];
    beta = camera[2];
    gamma = camera[3];
    t0 = camera[4];
    t1 = camera[5];

    x=point[0];
    y=point[1];
    z=point[2];

    sin_alpha=sin(alpha);
    cos_alpha=cos(alpha);
    sin_beta=sin(beta);
    cos_beta=cos(beta);
    sin_gamma=sin(gamma);
    cos_gamma=cos(gamma);

    temp1=(cos_beta*x+sin_alpha*sin_beta*y-cos_alpha*sin_beta*z)/(s*s);
    temp2=(cos_alpha*y+sin_alpha*z)/(s*s);
    right_R_u[0]=-cos_gamma*temp1+sin_gamma*temp2;
    right_R_v[0]=-sin_gamma*temp1-cos_gamma*temp2;

    temp1=(cos_alpha*sin_beta*y+sin_alpha*sin_beta*z)/s;
    temp2=(-sin_alpha*y+cos_alpha*z)/s;
    right_R_u[1]=cos_gamma*temp1-sin_gamma*temp2;
    right_R_v[1]=sin_gamma*temp1+cos_gamma*temp2;

    temp1=(-sin_beta*x+sin_alpha*cos_beta*y-cos_alpha*cos_beta*z)/s;
    right_R_u[2]=cos_gamma*temp1;
    right_R_v[2]=sin_gamma*temp1;

    temp1=(cos_beta*x+sin_alpha*sin_beta*y-cos_alpha*sin_beta*z)/s-t0;
    temp2=(cos_alpha*y+sin_alpha*z)/s-t1;
    right_R_u[3]= - sin_gamma*temp1-cos_gamma*temp2;
    right_R_v[3]= cos_gamma*temp1-sin_gamma*temp2;

    right_R_u[4]= - cos_gamma;
    right_R_u[5]= sin_gamma;

    right_R_v[4]= - sin_gamma;
    right_R_v[5]= - cos_gamma;
}

void bundle::PSFMDriver::compute_dfdR(double*p,double *imgpts,char *vmask,double*weight,double mu,int ncams,int n3Dpts,int k,double *dfdR){
    //计算近似函数对第k个相机参数的偏导，dfdR为输出，长度为6
    int i,j,index_temp,index_u,index_v,index_p,num_mask0;
    double temp_u,temp_v,*camera,*point3D,*uv,*right_R_u,*right_R_v;

    camera=new double[6];//当前相机参数
    point3D=new double[3];
    uv=new double[2];
    right_R_u=new double[6];
    right_R_v=new double[6];

    for(i=0;i<6;i++){
        dfdR[i]=0.0;
        camera[i] = p[k*6+i];
    }
    for(i=0;i<n3Dpts;i++){
        if((int)vmask[i*ncams+k] == 1){
            //取出第i个点的坐标
            point3D[0]=p[ncams*6+3*i];
            point3D[1]=p[ncams*6+3*i+1];
            point3D[2]=p[ncams*6+3*i+2];

            //计算当前投影
            compute_uv(camera,point3D,uv);

            //取出第i个点在第k个相机投影的坐标
            num_mask0 = 0;
            for(j=0;j<i*ncams+k;j++){
                num_mask0+=(int)vmask[j];
            }
            num_mask0 = i*ncams+k - num_mask0;
            index_u=2*i*ncams+2*k-2*num_mask0;
            index_v=2*i*ncams+2*k+1-2*num_mask0;

            //计算左项
            temp_u = weight[i*ncams+k]*(uv[0]-imgpts[index_u])/sqrt(pow((uv[0]-imgpts[index_u]),2)+mu*mu);
            temp_v = weight[i*ncams+k]*(uv[1]-imgpts[index_v])/sqrt(pow((uv[1]-imgpts[index_v]),2)+mu*mu);

            //计算右项
            compute_right_R(camera,point3D,right_R_u,right_R_v);

            //相乘后相加
            dfdR[0]+=temp_u*right_R_u[0]+temp_v*right_R_v[0];
            dfdR[1]+=temp_u*right_R_u[1]+temp_v*right_R_v[1];
            dfdR[2]+=temp_u*right_R_u[2]+temp_v*right_R_v[2];
            dfdR[3]+=temp_u*right_R_u[3]+temp_v*right_R_v[3];
            dfdR[4]+=temp_u*right_R_u[4]+temp_v*right_R_v[4];
            dfdR[5]+=temp_u*right_R_u[5]+temp_v*right_R_v[5];

            //std::cout<<temp_u<<" "<<right_R_u[0]<<" ";
        }
    }

    delete(camera);
    delete(point3D);
    delete(uv);
    delete(right_R_u);
    delete(right_R_v);
}

void bundle::PSFMDriver::compute_right_X(double *camera,double *point,double *right_X_u,double *right_X_v){
    //被compute_dfdX调用，因为公式太长了
    double s,alpha,beta,gamma,t0,t1,x,y,z,sin_alpha,cos_alpha,sin_beta,cos_beta,sin_gamma,cos_gamma,temp1,temp2;

    s = camera[0];
    alpha = camera[1];
    beta = camera[2];
    gamma = camera[3];
    t0 = camera[4];
    t1 = camera[5];

    x=point[0];
    y=point[1];
    z=point[2];

    sin_alpha=sin(alpha);
    cos_alpha=cos(alpha);
    sin_beta=sin(beta);
    cos_beta=cos(beta);
    sin_gamma=sin(gamma);
    cos_gamma=cos(gamma);

    right_X_u[0]=cos_gamma*cos_beta/s;
    right_X_u[1]=cos_gamma*sin_alpha*sin_beta/s-sin_gamma*cos_alpha/s;
    right_X_u[2]=-cos_gamma*cos_alpha*sin_beta/s-sin_gamma*sin_alpha/s;

    right_X_v[0]=sin_gamma*cos_beta/s;
    right_X_v[1]=sin_gamma*sin_alpha*sin_beta/s+cos_gamma*cos_alpha/s;
    right_X_v[2]=-sin_gamma*cos_alpha*sin_beta/s+cos_gamma*sin_alpha/s;
}

void bundle::PSFMDriver::compute_dfdX(double*p,double *imgpts,char *vmask,double *weight,double mu,int ncams,int n3Dpts,int k,double *dfdX){
    //计算近似函数对第k个三维点的偏导，dfdX为输出，长度为3
    int i,j,index_temp,index_u,index_v,index_p,num_mask0;
    double temp_u,temp_v,*camera,*point3D,*uv,*right_X_u,*right_X_v;

    camera=new double[6];//当前相机参数
    point3D=new double[3];
    uv=new double[2];
    right_X_u=new double[3];
    right_X_v=new double[3];

    for(i=0;i<3;i++){
        dfdX[i]=0;
        point3D[i]=p[6*ncams+3*k+i];
    }

    for(i=0;i<ncams;i++){

        if((int)vmask[k*ncams+i] == 1){

            //取出第i个相机的参数
            for(j=0;j<6;j++)camera[j]=p[i*6+j];

            //计算当前投影
            compute_uv(camera,point3D,uv);

            //取出第k个点在第i个相机投影的坐标
            num_mask0 = 0;
            for(j=0;j<k*ncams+i;j++){
                num_mask0+=(int)vmask[j];
            }
            num_mask0 = k*ncams+i - num_mask0;
            index_u=2*k*ncams+2*i-2*num_mask0;
            index_v=2*k*ncams+2*i+1-2*num_mask0;


            //计算左项
            temp_u = weight[k*ncams+i]* (uv[0]-imgpts[index_u])/sqrt(pow((uv[0]-imgpts[index_u]),2)+mu*mu);
            temp_v = weight[k*ncams+i]* (uv[1]-imgpts[index_v])/sqrt(pow((uv[1]-imgpts[index_v]),2)+mu*mu);

            //计算右项
            compute_right_X(camera,point3D,right_X_u,right_X_v);

            dfdX[0]+=temp_u*right_X_u[0]+temp_v*right_X_v[0];
            dfdX[1]+=temp_u*right_X_u[1]+temp_v*right_X_v[1];
            dfdX[2]+=temp_u*right_X_u[2]+temp_v*right_X_v[2];
        }
    }

    delete(camera);
    delete(point3D);
    delete(uv);
    delete(right_X_u);
    delete(right_X_v);
}

void bundle::PSFMDriver::compute_uv(double *camera,double *point,double *uv){
    //计算投影u
    double u,v,s,alpha,beta,gamma,t0,t1,x,y,z,temp1,temp2;

    s = camera[0];
    alpha = camera[1];
    beta = camera[2];
    gamma = camera[3];
    t0 = camera[4];
    t1 = camera[5];

    x=point[0];
    y=point[1];
    z=point[2];

    temp1=(cos(beta)*x+sin(alpha)*sin(beta)*y-cos(alpha)*sin(beta)*z)/s-t0;
    temp2=(cos(alpha)*y+sin(alpha)*z)/s-t1;

    u=cos(gamma)*temp1-sin(gamma)*temp2;
    v=sin(gamma)*temp1+cos(gamma)*temp2;

    uv[0]=u;
    uv[1]=v;
}

double bundle::PSFMDriver::f_L2(double *p,double *imgpts, char *vmask,double *weight,int ncams,int n3Dpts){
    //计算L2残差
    int i,j,k,index_uv;
    double *camera,*point3D,*uv,f_out,u_fid,v_fid;

    camera = new double[6];
    point3D = new double[3];
    uv = new double[2];

    f_out =0;
    index_uv = 0;
    for(i=0;i<n3Dpts;i++){
        for(j=0;j<ncams;j++){
            if((int)vmask[i*ncams+j] == 1){
                memcpy(camera,p+j*6,6*sizeof(double));
                memcpy(point3D,p+6*ncams+i*3,3*sizeof(double));
                compute_uv(camera,point3D,uv);
                u_fid = imgpts[index_uv];
                index_uv++;
                v_fid = imgpts[index_uv];
                index_uv++;

                //f_out +=abs(uv[0]-u_fid)+abs(uv[1]-v_fid);
                f_out +=sqrt(weight[i*ncams+j]*pow(uv[0]-u_fid,2) + weight[i*ncams+j]*pow(uv[1]-v_fid,2));
                // f_out +=weight[i*ncams+j]*pow(uv[1]-v_fid,2);
            }
        }
    }

    delete(camera);
    delete(point3D);
    delete(uv);

    return f_out;
}

double bundle::PSFMDriver::Armijo_f(double *p,double mu,double* grad,double sigma,double rho,double *imgpts,char * vmask,double *weight,int ncams,int n3Dpts){
    /*************
    回退法寻找符合Armijo准则的步长
    ************/
    int i,ncols;
    bool is_continue;
    double alpha,*p_new,*p_new_temp,dot,right,f_right;

    ncols=6*ncams+3*n3Dpts;
    p_new = new double[ncols];

    alpha = 1.0;
    dot = 0.0;
    for(i=0;i<ncols;i++){
        p_new[i]=p[i]-alpha*grad[i];
        dot += (-grad[i]*grad[i]);
    }

    std::clock_t start, end;
    start = clock();//debug;

    f_right=f(p,mu,imgpts,vmask,weight,ncams,n3Dpts);
    right =f_right  + sigma*alpha*dot;
    is_continue = (f(p_new,mu,imgpts,vmask,weight,ncams,n3Dpts) > right);
    
    while( is_continue ){
        
        end=clock();
        //std::cout<<(double(end-start)/CLOCKS_PER_SEC>60)<<std::endl;
        //if(double(end-start)/CLOCKS_PER_SEC>10) break;//debug


        alpha *= rho;
        for(i=0;i<ncols;i++){
            p_new[i]=p[i]-alpha*grad[i];
        }
        right =f_right  + sigma*alpha*dot;
        is_continue = (f(p_new,mu,imgpts,vmask,weight,ncams,n3Dpts) > right);
    }

    delete(p_new);

    return alpha;
}

void bundle::PSFMDriver::update_p(double *p,double mu,double *grad,double sigma,double rho,double *imgpts,char *vmask,double *weight,int ncams,int n3Dpts){
    /*******************
    梯度法迭代，方向为负梯度，步长由Armijo准则得出
    ********************/
    int i,ncols;
    double alpha;

    ncols=6*ncams+3*n3Dpts;
    alpha = Armijo_f(p,mu,grad,sigma,rho,imgpts,vmask,weight,ncams,n3Dpts);

    for(i=0;i<ncols;i++){
        p[i]=p[i]-alpha*grad[i];
    }

}



// void bundle::PSFMDriver::Toone_recover_new(double *motstr,double *imgpts,int ncams,int n3Dpts,int n2Dpts,double mean_data,double std_data,double *gamma_old,double *motstr_recover,double *imgpts_recover){
//     int i;
//     double sin_gamma_i,cos_gamma_i;
//     for(i=0;i<ncams;i++){
//         sin_gamma_i = sin(gamma_old[i]);
//         cos_gamma_i = cos(gamma_old[i]);
//         motstr_recover[i*6]=motstr[i*6];
//         motstr_recover[i*6+1]=motstr[i*6+1];
//         motstr_recover[i*6+2]=motstr[i*6+2];
//         motstr_recover[i*6+3]=motstr[i*6+3];
//         motstr_recover[i*6+4]=motstr[i*6+4]*std_data-mean_data*(cos_gamma_i+sin_gamma_i);
//         motstr_recover[i*6+5]=motstr[i*6+5]*std_data-mean_data*(-sin_gamma_i+cos_gamma_i);
//     }
//     for(i=6*ncams;i<ncams*6+n3Dpts*3;i++){
//         motstr_recover[i]=motstr[i]*std_data;
//     }
//     for(i=0;i<n2Dpts;i++){
//         imgpts_recover[i] = imgpts[i]*std_data+mean_data;
//     }

// }


void bundle::PSFMDriver::Toone_recover_new(double *motstr,double *imgpts,int ncams,int n3Dpts,int n2Dpts,double mean_data,double std_data,double *gamma_old,double *motstr_recover,double *imgpts_recover){
    int i;
    double sin_gamma_i,cos_gamma_i,sin_new,cos_new;
    for(i=0;i<ncams;i++){
        sin_gamma_i = sin(gamma_old[i]);
        cos_gamma_i = cos(gamma_old[i]);
        sin_new = sin(motstr[i*6+3]);
        cos_new = cos(motstr[i*6+3]);
        motstr_recover[i*6]=motstr[i*6];
        motstr_recover[i*6+1]=motstr[i*6+1];
        motstr_recover[i*6+2]=motstr[i*6+2];
        motstr_recover[i*6+3]=motstr[i*6+3];
        motstr_recover[i*6+4]=motstr[i*6+4]*std_data-mean_data*(cos_new+sin_new);
        motstr_recover[i*6+5]=motstr[i*6+5]*std_data-mean_data*(-sin_new+cos_new);
    }
    for(i=6*ncams;i<ncams*6+n3Dpts*3;i++){
        motstr_recover[i]=motstr[i]*std_data;
    }
    for(i=0;i<n2Dpts;i++){
        imgpts_recover[i] = imgpts[i]*std_data+mean_data;
    }

}

// void bundle::PSFMDriver::sba_L1_smoothing_grad(double *p, double *imgpts, char *vmask, int ncams, int n3Dpts,double *newp) {
void bundle::PSFMDriver::sba_L1_smoothing_grad(double *p, double *imgpts, char *vmask, int ncams, int n3Dpts,double *newp) {
    int i,k,maxstep,ncols;
    double sigma,rho,mu,epsilon,omiga,norm1_grad,*grad,error_L1,*error_vec,error_L2,f_now,*weight;
    std::vector<double> error_everystep_L1,error_everystep_L2,norm1_grad_everystep;
    std::ofstream fffout,ffout1;
    std::clock_t start, end;


    sigma=0.05;
    rho=0.5;
    mu=1;
    epsilon=0.001;
    omiga=1000;
    maxstep=1000;
    //手动设置的参数

    k=0;
    ncols = ncams*6 + n3Dpts*3;
    grad = new double[ncols];
    error_vec = new double[ncams*n3Dpts];
    weight = new double[ncams*n3Dpts];
    memcpy(newp,p,sizeof(double)*ncols);
    for(i=0;i<ncams*n3Dpts;i++) weight[i] = 1.0;

    //初始化
    compute_gradient(newp,imgpts,mu,vmask,weight,ncams,n3Dpts,grad);
    norm1_grad=compute_norm1(grad,ncols);

    error_L1=f(newp,0,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2;
    error_L2=f_L2(newp,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2;
    f_now=f(newp,mu,imgpts,vmask,weight,ncams,n3Dpts);

    std::cout<<"/****************"<<"intial"<<"/****************/"<<std::endl<<"mu is "<<mu<<",f is "<<f_now<<",error_L1 is "<<error_L1<<",error_L2 is "<<error_L2<<",norm_grad is "<<norm1_grad<<std::endl;

    start = clock();//debug;
    fffout.open("out/output.txt");
    ffout1.open("out/output_x_everystep.txt");
    
    while(norm1_grad>epsilon){
        
        //*************out every step*****************//
        end=clock();
        if(double(end-start)/CLOCKS_PER_SEC>600) break;

        
        error_L1=f(newp,0,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2;
        error_L2=f_L2(newp,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2;
        f_now=f(newp,mu,imgpts,vmask,weight,ncams,n3Dpts);
        error_everystep_L1.push_back(error_L1);
        error_everystep_L2.push_back(error_L2);
        norm1_grad_everystep.push_back(norm1_grad);
        
        fffout<<"/****************/"<<k<<"/****************/"<<std::endl<<"mu is "<<mu<<",f is "<<f_now<<",error_L1 is "<<error_L1<<",error_L2 is "<<error_L2<<",norm_grad is "<<norm1_grad<<std::endl;
        ffout1<<"/****************/"<<k<<"/****************/"<<std::endl;
        for(i=0;i<ncams*6+n3Dpts*3;i++)ffout1<<newp[i]<<std::endl;
        //std::cout<<"/****************/"<<k<<"/****************/"<<std::endl<<"mu is "<<mu<<",f is "<<f_now<<",error_L1 is "<<error_L1<<",error_L2 is "<<error_L2<<",norm_grad is "<<norm1_grad<<std::endl;
        
        //*************out every step*****************//

        update_p(newp,mu,grad,sigma,rho,imgpts,vmask,weight,ncams,n3Dpts);
        compute_gradient(newp,imgpts,mu,vmask,weight,ncams,n3Dpts,grad);
        norm1_grad=compute_norm1(grad,ncols);
        if(norm1_grad<omiga*mu) mu*=sigma;
        k++;
        
        if(k>maxstep) break;
    }
 
    fffout.close();
    ffout1.close();
    std::cout<<"/****************final****************/"<<std::endl<<"mu is "<<mu<<",f is "<<f(newp,mu,imgpts,vmask,weight,ncams,n3Dpts)<<",error_L1 is "<<f(newp,0,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2<<",error_L2 is "<<f_L2(newp,imgpts,vmask,weight,ncams,n3Dpts)/(ncams*n3Dpts)/2<<",norm_grad is "<<norm1_grad<<std::endl;
    

    fffout.open("out/error_L2_step.txt");
    for(i=0;i<error_everystep_L2.size();i++) fffout<<error_everystep_L2[i]<<std::endl;
    fffout.close();

    // delete(grad);
    // delete(error_vec);
    // delete(weight);
}


int bundle::PSFMDriver::Run_l1(double* motstruct, double* imgpts, char* vmask1, const int n3Dpts1, const int ncams, const int nconcam, const int ncon3Dpts)
{
    std::cout<<"L1"<<std::endl;
    // int n3Dpts = 39;
    int n3Dpts = n3Dpts1;
    double *p_recover,*pts2d_recover,*vec_gamma;
    p_recover = new double[ncams*6+n3Dpts*3];
    vec_gamma = new double[ncams];
    std::cout<<"n3Dpts"<<n3Dpts<<std::endl;
    std::ifstream fin;
    double ratio_noise,ratio_outlier,ratio_noise_cams,range_img,mean_one,std_one;
    data_initial(motstruct, imgpts, vmask1, vec_gamma, ncams, n3Dpts, 0, 0, 0, 1024, mean_one, std_one);
    // double *time_smoothing,*error_smoothing;
    double time_smoothing,error_smoothing;
    std::clock_t start, end;
    int n2Dpts = 0;
    double *weight_one = new double[n3Dpts*ncams];
    char *vmask = new char[n3Dpts*ncams];
    double *mot = new double[ncams*6+n3Dpts*3];

    fin.open("simulate_data/mot_random.txt", std::ios::in);
    for (int i = 0; i < ncams*6+n3Dpts*3; i++) fin >> mot[i];
    fin.close();
    fin.open("simulate_data/vmask_random.txt", std::ios::in);
    for (int i = 0; i < n3Dpts*ncams; i++) {
        fin >> vmask[i];
        vmask[i] -= '0';
        weight_one[i] = 1.0;
        if((int)vmask[i] == 1) n2Dpts++;
    }
    fin.close();
    n2Dpts *=2;

    std::cout<<"n2Dpts:"<<n2Dpts<<std::endl;
    double *pts2d = new double[n2Dpts];
    fin.open("simulate_data/pts2d_random.txt", std::ios::in);
    for (int i = 0; i < n2Dpts; i++) fin >> pts2d[i];
    fin.close();
    double *mot_temp = new double[ncams*6+n3Dpts*3];
    double *mot_new = new double[ncams*6+n3Dpts*3];

    double_copy(mot_temp,mot,ncams*6+n3Dpts*3);
    std::cout<<"initial error is "<<f(mot_temp,0,pts2d,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2)<<std::endl;
    start = clock();
    sba_L1_smoothing_grad(mot_temp,pts2d,vmask,ncams,n3Dpts,mot_new);
    end = clock();

//     // std::cout<<mean_one<<"    "<<std_one<<std::endl;
    pts2d_recover = new double[n2Dpts];

    time_smoothing = double(end-start)/CLOCKS_PER_SEC;
    error_smoothing = f(mot_new,0,pts2d,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2);
    std::cout<<"Final L1 error of smoothing is "<<error_smoothing<<", running time is "<<time_smoothing<<"s"<<std::endl;
    // write_txt_app<double>(time_smoothing,1,"output_temp/time_smoothing.txt");
    // write_txt_app<double>(error_smoothing,1,"output_temp/error_smoothing.txt");
    // double *p_recover,*pts2d_recover,*vec_gamma;
    // p_recover = new double[ncams*6+n3Dpts*3];
    // pts2d_recover = new double[n2Dpts];
    // vec_gamma = new double[ncams];
    // std::cout<<"2222222"<<std::endl;
    // mean_one = 1;
    // fin.open("simulate_data_noOne/one_avg.txt", std::ios::in);
    // std::cout<<"???"<<std::endl;
    // if (!fin.is_open()) {
    //     std::cerr << "Failed to open one_avg.txt!" << std::endl;
    //     exit(1);
    // }
    // std::cout<<"!!!!!"<<std::endl;
    // // for (int i = 0; i < 1; i++) fin >> mean_one;
    // std::cout<<"*************"<<std::endl;
    // // fin.close();
    // std::cout<<"333333333"<<std::endl;
    // fin.open("simulate_data_noOne/one_std.txt", std::ios::in);
    // // for (int i = 0; i < 1; i++) fin >> std_one;
    // // fin.close();
    // std::cout<<"4444444444"<<std::endl;
    // fin.open("simulate_data_noOne/vec_gamma.txt", std::ios::in);
    // for (int i = 0; i < ncams; i++) fin >> vec_gamma[i];
    // fin.close();
    // std::cout<<"555555"<<std::endl;
    // std::cout<<"avg is "<<mean_one<<", std is "<<std_one<<std::endl;

    Toone_recover_new(mot_new,pts2d,ncams,n3Dpts,n2Dpts,mean_one,std_one,vec_gamma,p_recover,pts2d_recover);
    error_smoothing = f(p_recover,0,pts2d_recover,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts*2);
    std::cout<<"Final recovery L1 error is "<<error_smoothing<<std::endl;
    double out_error_l2;
    out_error_l2=f_L2(p_recover,pts2d_recover,vmask,weight_one,ncams,n3Dpts)/(ncams*n3Dpts)/2;
    std::cout<<"Final recovery L2 error is "<<out_error_l2<<std::endl;
//     // write_txt_app<double>(error_smoothing,1,"error_Smoothing.txt");

    std::memcpy(motstruct, p_recover, (ncams * 6 + n3Dpts * 3) * sizeof(double));
    std::memcpy(imgpts, pts2d_recover, n2Dpts * sizeof(double));

    

    delete[] p_recover;
    delete[] pts2d_recover;
    delete[] vec_gamma;
    delete[] mot;
    delete[] pts2d;
    delete[] mot_temp;
    delete[] weight_one;
    delete[] mot_new;
    // delete[] time_smoothing;
    // delete[] error_smoothing;
    delete[] vmask;

    // for(int i=0;i<ncams*6;i++){
    //     std::cout<<mot_new[i]<<std::endl;
    // }
    
}


