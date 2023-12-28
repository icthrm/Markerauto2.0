#ifndef CAMERA_H__
#define CAMERA_H__

#include "util/matrix.h"

namespace mx{

struct camera_params{
    double calib[5];
    double mot[6];
    mx::Matrix<3,3,double> R;
    mx::Matrix<3,1,double> T;
    bool fixedflag;
    bool valid;
};

struct pproj_params{
	union{
		double mot[6];
		struct{
			double s;
			double alpha;
			double beta;
			double gamma;
			double t0;
			double t1;
		};
	};
	bool fixedflag;
    bool valid;
};

void PProjPCopyFormMotParams(pproj_params* pparams, const double* motparams);
void MotCopyFormPProjParams(double* motparams , const pproj_params& pparams);
void MxKFromCameraParams(const camera_params& camparams, double K[9]);
void CameraInitFormKRT(camera_params* camparams, const double K[9], const double R[9], const double T[3]);
void CameraCopyFormMotParams(camera_params* camparams, const double* motparams, int cnp = 11);
void MotCopyFormCameraParams(double* motparams , const camera_params& camparams, int cnp = 11);
void MxPFromCameraParams(const camera_params& camparams, double P[12]);
void MxPinvFromCameraParams(const camera_params& camparams, double P_inv[12]);
void MxCFromCameraParams(const camera_params& camparams, double C[4]);
void PrintMxCofCamera(const camera_params& camparams, std::ostream& o);

}

#endif