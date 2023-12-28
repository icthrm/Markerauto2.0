#include "sfm_driver.h"
#include "img_projs.h"
#include <cmath>
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
#include <lm.h>
#include "ceres/ceres.h"
// #include <ext/hash_map>

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

        EX_BEGIN_CLOCK()
//         Run(motstruct, projections, vmask, num_3d_pts, ncams, nconcam, fixed_pt_num);
        Run_ceres(motstruct, projections, vmask, num_3d_pts, ncams, nconcam, fixed_pt_num, W_init, use_MAND);
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

    std::cout<<"num_dists:"<<num_dists<<"   dist_total"<<dist_total<<std::endl;

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

    std::string res_file = "/home/xzh/桌面/projections.txt";
        std::ofstream outputfile_res(res_file);

        for(int i=0;i<img_2d.size();i++)
        {
            outputfile_res<<img_2d[i].x<<"   "<<img_2d[i].y<<std::endl;
    //         if(i%20==0) outputfile_res<<"\n";
        }

        outputfile_res.close();

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


   	std::string res_file1 = "/home/xzh/桌面/残差/hiv_10643/all/W.txt";
    std::ofstream outputfile_res1(res_file1);

    for(int i=0;i<W_ord.size();i++)
    {
        for(int j=0;j<W_ord[i].size();j++)
        {
            outputfile_res1<<W_ord[i][j]<<"  "<<std::endl;
        }
//         outputfile_res<<"\n";
    }

    outputfile_res1.close();
    
    
    return n;
}


