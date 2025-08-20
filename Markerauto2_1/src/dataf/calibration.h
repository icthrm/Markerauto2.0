#ifndef CALIBRATION_H__
#define CALIBRATION_H__

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include "matrix/matrix.h"

namespace util{
    
/** base on image; first shift, then rotation, scale or stretch */
struct calibrated{
    double x_shift;
    double y_shift;
    double global_rot;				//in image rotation of x and y axis;make the x axis in horizon
    double x_stretch_rot;			//after the global rotation, stretch x as x_stretch_rot
    double y_stretch_rot;			//after the global rotation, stretch y as y_stretch_rot
    double tilt_vector[3];			//a direct vector of tilt-plane & yOz plane
    double tilt_angle;
    double scale;
    bool valid;
};

inline void PrintCalibrationAsImod(const std::vector<calibrated>& calib_v, const float rotation, 
				   const float scale, const char* xffilename, const char* anglefilename, const char* exfilename){
    std::ofstream out(xffilename); 
    std::ofstream agout(anglefilename);
    std::ofstream exceptf(exfilename); 
    for(int i = 0; i < calib_v.size(); i++){
	double cos_x_rot, sin_x_rot;
	cos_x_rot = cos(calib_v[i].global_rot);
	sin_x_rot = sin(calib_v[i].global_rot);
	double cos_y_rot, sin_y_rot;
	cos_y_rot = cos(calib_v[i].global_rot);
	sin_y_rot = sin(calib_v[i].global_rot);
	double A1[4] = {cos_x_rot, -sin_y_rot, sin_x_rot, cos_y_rot};
	matrix_invert_inplace(2, A1);
	
	double cos_x_stretch_rot, sin_x_stretch_rot;
	cos_x_stretch_rot = cos(calib_v[i].x_stretch_rot);
	sin_x_stretch_rot = sin(calib_v[i].x_stretch_rot);
	double cos_y_stretch_rot, sin_y_stretch_rot;
	cos_y_stretch_rot = cos(calib_v[i].y_stretch_rot);
	sin_y_stretch_rot = sin(calib_v[i].y_stretch_rot);
	double A2[4] = {cos_x_stretch_rot, -sin_y_stretch_rot, sin_x_stretch_rot, cos_y_stretch_rot};
	matrix_invert_inplace(2, A2);
	
	double A[4];	//A = A2*A1
	matrix_product(2, 2, 2, 2, A2, A1, A);
	
	double dx = A[0]*calib_v[i].x_shift+A[1]*calib_v[i].y_shift;
	double dy = A[2]*calib_v[i].x_shift+A[3]*calib_v[i].y_shift;

	if(calib_v[i].valid){
	    out<<(A[0]*cos(rotation)+A[1]*sin(rotation))/calib_v[i].scale
               <<"\t"<<(-A[0]*sin(rotation)+A[1]*cos(rotation))/calib_v[i].scale
               <<"\t"<<(A[2]*cos(rotation)+A[3]*sin(rotation))/calib_v[i].scale
               <<"\t"<<(-A[2]*sin(rotation)+A[3]*cos(rotation))/calib_v[i].scale
               <<"\t"<<dx*scale/calib_v[i].scale<<"\t"<<dy*scale/calib_v[i].scale<<std::endl;
	}
	else{
	    out<<1<<"\t"<<0<<"\t"<<0<<"\t"<<1<<"\t"<<0<<"\t"<<0<<std::endl;
	    exceptf<<i<<" "<<std::endl;
	}
	agout<<calib_v[i].tilt_angle/M_PI*180<<std::endl;
    }
    out.close();
    agout.close();
}

}

#endif