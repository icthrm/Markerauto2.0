#include "img_projs.h"
#include <cmath>
#include <iostream>
#include "dataf/quat_vec.h"
#include "sba.h"

// using namespace bundle;

void calcImgParallelProj(const double param[6], const double M[3], double n[2])
{
	double s, alpha, beta, gamma, t0, t1;
	double X, Y, Z;
	s = param[0]; alpha = param[1]; beta = param[2]; gamma = param[3]; t0 = param[4]; t1 = param[5];
	X = M[0]; Y = M[1]; Z = M[2];
	
	double cos_alpha = cos(alpha);
	double sin_alpha = sin(alpha);
	double cos_beta = cos(beta);
	double sin_beta = sin(beta);
	double cos_gamma = cos(gamma);
	double sin_gamma = sin(gamma);
	
	double tm1, tm2;
	tm1 = (cos_beta*X+sin_alpha*sin_beta*Y-cos_alpha*sin_beta*Z)/s-t0;
	tm2 = (cos_alpha*Y+sin_alpha*Z)/s-t1;
//     std::cout<<"gamma:"<<beta<<std::endl;
//         if(std::abs(beta)>0.7)
//     {
//     std::cout<<"s:"<<param[0]<<"  alpha:"<<param[1]<<"  beta:"<<param[2]
//             <<"  gamma:"<<param[3]<<"  t0:"<<param[4]<<"  t1:"<<param[5]<<std::endl;
//         std::cout<<"X:"<<X<<"Y:"<<Y<<"Z:"<<Z<<std::endl;
//         std::cout<<"cos_alpha:"<<cos_alpha<<"sin_alpha:"<<sin_alpha<<std::endl;
//         std::cout<<"cos_beta:"<<cos_beta<<"sin_beta:"<<sin_beta<<std::endl;
//         std::cout<<"cos_gamma:"<<cos_gamma<<"sin_gamma:"<<sin_gamma<<std::endl;
//         std::cout<<"tm1:"<<tm1<<"tm2:"<<tm2<<std::endl;
//     }
//     std::cout<<"tm1:"<<tm1<<"tm2:"<<tm2<<std::endl;
	n[0] = cos_gamma*tm1-sin_gamma*tm2;
	n[1] = sin_gamma*tm1+cos_gamma*tm2;
}


void calccalcImgParallelProjJac(double param[6], double M[3], double jacmP[2][6], double jacmS[2][3])
{
	double s, alpha, beta, gamma, t0, t1;
	double X, Y, Z;
	s = param[0]; alpha = param[1]; beta = param[2]; gamma = param[3]; t0 = param[4]; t1 = param[5];
	X = M[0]; Y = M[1]; Z = M[2];
	
	double cos_alpha = cos(alpha);
	double sin_alpha = sin(alpha);
	double cos_beta = cos(beta);
	double sin_beta = sin(beta);
	double cos_gamma = cos(gamma);
	double sin_gamma = sin(gamma);
	double s_s = 1/s;
	jacmP[0][0] = (-cos_gamma*(cos_beta*X+sin_alpha*sin_beta*Y-cos_alpha*sin_beta*Z)+sin_gamma*(cos_alpha*Y+sin_alpha*Z))*s_s*s_s;
	jacmP[0][1] = (cos_gamma*(cos_alpha*sin_beta*Y+sin_alpha*sin_beta*Z)-sin_gamma*(-sin_alpha*Y+cos_alpha*Z))*s_s;
	jacmP[0][2] = cos_gamma*(-sin_beta*X+sin_alpha*cos_beta*Y-cos_alpha*cos_beta*Z)*s_s;
	jacmP[0][3] = -sin_gamma*((cos_beta*X+sin_alpha*sin_beta*Y-cos_alpha*sin_beta*Z)*s_s-t0)-cos_gamma*((cos_alpha*Y+sin_alpha*Z)*s_s-t1);
	jacmP[0][4] = -cos_gamma;
	jacmP[0][5] = sin_gamma;
	jacmP[1][0] = (-sin_gamma*(cos_beta*X+sin_alpha*sin_beta*Y-cos_alpha*sin_beta*Z)-cos_gamma*(cos_alpha*Y+sin_alpha*Z))*s_s*s_s;
	jacmP[1][1] = (sin_gamma*(cos_alpha*sin_beta*Y+sin_alpha*sin_beta*Z)+cos_gamma*(-sin_alpha*Y+cos_alpha*Z))*s_s;
	jacmP[1][2] = sin_gamma*(-sin_beta*X+sin_alpha*cos_beta*Y-cos_alpha*cos_beta*Z)*s_s;
	jacmP[1][3] = cos_gamma*((cos_beta*X+sin_alpha*sin_beta*Y-cos_alpha*sin_beta*Z)*s_s-t0)-sin_gamma*((cos_alpha*Y+sin_alpha*Z)*s_s-t1);
	jacmP[1][4] = -sin_gamma;
	jacmP[1][5] = -cos_gamma;
	
	jacmS[0][0] = cos_gamma*cos_beta*s_s;
	jacmS[0][1] = (cos_gamma*sin_alpha*sin_beta-sin_gamma*cos_alpha)*s_s;
	jacmS[0][2] = (-cos_gamma*cos_alpha*sin_beta-sin_gamma*sin_alpha)*s_s;
	jacmS[1][0] = sin_gamma*cos_beta*s_s;
	jacmS[1][1] = (sin_gamma*sin_alpha*sin_beta+cos_gamma*cos_alpha)*s_s;
	jacmS[1][2] = (-sin_gamma*cos_alpha*sin_beta+cos_gamma*sin_alpha)*s_s;
}


/** FULL BUNDLE ADJUSTMENT, I.E. SIMULTANEOUS ESTIMATION OF CAMERA AND STRUCTURE PARAMETERS */

/**Given a parameter vector p made up of the 3D coordinates of n points and the parameters of m cameras, compute in
 * hx the prediction of the measurements, i.e. the projections of 3D points in the m images. The measurements
 * are returned in the order (hx_11^T, .. hx_1m^T, ..., hx_n1^T, .. hx_nm^T)^T, where hx_ij is the predicted
 * projection of the i-th point on the j-th camera.
 * Notice that depending on idxij, some of the hx_ij might be missing
 *
 */
void img_ParallelProj_x(double* p, sba_crsm* idxij, int* rcidxs, int* rcsubs, double* hx, void* adata)
{
    int cnp = 6, pnp = 3, mnp = 2;
    double *pa, *pb, *ppt, *pmeas, *pcalib;
    //int n;
    int m, nnz;

    //n=idxij->nr;
    m=idxij->nc;
    pa = p;					//the point to params of cameras
    pb = p+m*cnp;			//the point to params of 3D points

    for(register int j = 0; j < m; ++j){
        /* j-th camera parameters */
        pcalib = pa+j*cnp;

        nnz = sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */

        for(register int i = 0; i < nnz; ++i){
            ppt = pb+rcsubs[i]*pnp;
            pmeas = hx+idxij->val[rcidxs[i]]*mnp; // set pmeas to point to hx_ij

            calcImgParallelProj(pcalib, ppt, pmeas); // evaluate Q in pmeas
        }
    }
}

/**Given a parameter vector p made up of the 3D coordinates of n points and the parameters of m cameras, compute in
 * jac the jacobian of the predicted measurements, i.e. the jacobian of the projections of 3D points in the m images.
 * The jacobian is returned in the order (A_11, ..., A_1m, ..., A_n1, ..., A_nm, B_11, ..., B_1m, ..., B_n1, ..., B_nm),
 * where A_ij=dx_ij/db_j and B_ij=dx_ij/db_i (see HZ).
 * Notice that depending on idxij, some of the A_ij, B_ij might be missing
 *
 */
void img_ParallelProj_jac_x(double* p, sba_crsm* idxij, int* rcidxs, int* rcsubs, double* jac, void* adata)
{
	int cnp = 6, pnp = 3, mnp = 2;
    double *pa, *pb, *ppt, *pcalib, *pA, *pB;
    //int n;
    int m, nnz, Asz, Bsz, ABsz;

    //n=idxij->nr;
    m=idxij->nc;
    pa=p;
    pb=p+m*cnp;
    Asz=mnp*cnp;
    Bsz=mnp*pnp;
    ABsz=Asz+Bsz;
	
	struct pglobs *gl;
  
	gl=(struct pglobs *)adata;

    for(register int j = 0; j < m; ++j){
        /* j-th camera parameters */
        pcalib=pa+j*cnp;

        nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */

        for(register int i = 0; i < nnz; ++i){
            ppt=pb + rcsubs[i]*pnp;
            pA=jac + idxij->val[rcidxs[i]]*ABsz; // set pA to point to A_ij
            pB=pA  + Asz; // set pB to point to B_ij

            calccalcImgParallelProjJac(pcalib, ppt, (double(*)[6])pA, (double(*)[3])pB); // evaluate dQ/da, dQ/db in pA, pB
			
			for(int ii = 0; ii < mnp; ++ii, pA+=cnp){
				for(int jj = 0; jj < cnp; ++jj){
					if(gl->cpmask[jj]){
						pA[jj] = 0;
					}
				}
			}
			
			for(int ii = 0; ii < mnp; ++ii, pB+=pnp){
				for(int jj = 0; jj < pnp; ++jj){
					if(gl->pmask){
						pB[jj] = 0;
					}
				}
			}
        }
    }
}


