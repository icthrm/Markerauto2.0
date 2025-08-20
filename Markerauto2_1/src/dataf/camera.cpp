#include "camera.h"
#include "quat_vec.h"
#include <cstring>
#include "matrix/matrix.h"
#include"nbundle/params.h"

void mx::CameraInitFormKRT(camera_params* camparams, const double K[9], const double R[9], const double T[3])
{
    memcpy(camparams->R.D(), R, sizeof(double)*9);
    camparams->calib[0] = K[0]; 			// fu
    camparams->calib[1] = K[2]; 			// u0
    camparams->calib[2] = K[5]; 			// v0
    camparams->calib[3] = K[8]; 			// ar
    camparams->calib[4] = 0; 				// s
    double quatRT[7];
    RotMatrix2Quat(R, quatRT);
    memcpy(quatRT+4, T, sizeof(double)*3);
    Quat2Vec(quatRT, 7, camparams->mot, 6);
    camparams->R = mx::Matrix<3,3,double>(R);
    camparams->T = mx::Matrix<3,1,double>(T);
}

void mx::CameraCopyFormMotParams(camera_params* camparams, const double* motparams, int cnp)
{
//     if(camparams->fixedflag){
// 	return;
//     }
    const double* params = motparams;
    assert(cnp == 11 || cnp ==6);
    if(cnp == 11){
	memcpy(camparams->calib, params, sizeof(double)*5);
	params += 5;
    }
    memcpy(camparams->mot, params, sizeof(double)*6);
    double quatR[4];
    _MK_QUAT_FRM_VEC(quatR, params);
    double R_[9];
    Quat2RotMatrix(quatR, R_);
    camparams->R = mx::Matrix<3,3,double>(R_);
    params += 3;
    camparams->T.D(0) = params[0];
    camparams->T.D(1) = params[1];
    camparams->T.D(2) = params[2];
    
//     camparams->fixedflag = true;
}

void mx::MotCopyFormCameraParams(double* motparams, const camera_params& camparams, int cnp)
{
    double* params = motparams;
    assert(cnp == 11 || cnp ==6);
//     assert(camparams.fixedflag);
    if(cnp == 11){
	memcpy(params, camparams.calib, sizeof(double)*5);
	params += 5;
    }
    memcpy(params, camparams.mot, sizeof(double)*6);
}

void mx::MxKFromCameraParams(const camera_params& camparams, double K[9])
{
    K[0] = camparams.calib[0]; K[1] = 0; K[2] = camparams.calib[1];
    K[3] = 0; K[4] = camparams.calib[0]/camparams.calib[3]; K[5] = camparams.calib[2];
    K[6] = 0; K[7] = 0; K[8] = 1;
}

void mx::MxPFromCameraParams(const camera_params& camparams, double P[12])
{
    double K[9] = {camparams.calib[0], 0, camparams.calib[1],
		    0, camparams.calib[0]/camparams.calib[3], camparams.calib[2],
		    0, 0, 1};
    double Rt[12] = {camparams.R.D(0), camparams.R.D(1), camparams.R.D(2), camparams.T.D(0),
		    camparams.R.D(3), camparams.R.D(4), camparams.R.D(5), camparams.T.D(1),
		    camparams.R.D(6), camparams.R.D(7), camparams.R.D(8), camparams.T.D(2)};
//     matrix_print(3, 3, K);
//     matrix_print(3, 4, R);
    matrix_product(3, 3, 3, 4, K, Rt, P);
}

void mx::MxCFromCameraParams(const camera_params& camparams, double C[4])
{
    double R_inv[9];
    matrix_invert(3, camparams.R.D(), R_inv);
    matrix_product(3, 3, 3, 1, R_inv, camparams.T.D(), C);
    C[0] = -C[0]; C[1] = -C[1]; C[2] = -C[2];
    C[3] = 1;
}

void mx::MxPinvFromCameraParams(const camera_params& camparams, double P_inv[12])
{
    double P[12];
    MxPFromCameraParams(camparams, P);
    double P_t[12];
    matrix_transpose(3, 4, P, P_t);
    double PxP_t[9];
    matrix_product(3, 4, 4, 3, P, P_t, PxP_t);
    matrix_invert_inplace(3, PxP_t);
    matrix_product(4, 3, 3, 3, P_t, PxP_t, P_inv);
}

void mx::PrintMxCofCamera(const camera_params& camparams, std::ostream& o)
{
    double C[4];
    MxCFromCameraParams(camparams, C);
    o<<C[0]<<" "<<C[1]<<" "<<C[2]<<" "<<C[3]<<std::endl;
}

void mx::MotCopyFormPProjParams(double* motparams, const mx::pproj_params& pparams)
{
	memcpy(motparams, pparams.mot, sizeof(double)*6);
}

void mx::PProjPCopyFormMotParams(mx::pproj_params* pparams, const double* motparams)
{
	memcpy(pparams->mot, motparams, sizeof(double)*6);
}

void mx::InitMot(float angle, double* params)
{
    bundle::PInitating::InitSR123T12(params, 1.0/bundle::PInitating::f, 0, angle, 0, 0, 0);
}


void mx::InitCamera(float angle, mx::pproj_params* camera)
{
    double mot[6];
    InitMot(angle, mot);
	mx::PProjPCopyFormMotParams(camera, mot);
}

// void mx::MotCopyFormPoints(double* point, const int& points)
// {
// 	memcpy(motparams, pparams.mot, sizeof(double)*6);
// }

void mx::MotCopyFormPoints(double* X, const double points[3])
{
    for(int i=0;i<3;i++)
    {
        X[0]=points[0];
        X[1]=points[1];
        X[2]=points[2];
    }
}



