#ifndef PARAMS_H__
#define PARAMS_H__

#include "dataf/dataf.h"
#include "sba.h"
#include <vector>
#include <list>
#include <sys/stat.h>
#include <iostream>
#include "matrix/matrix.h"
#include "matrix/vector.h"
#include "dataf/quat_vec.h"
#include "util/matrix.h"
#include "dataf/camera.h"
#include <algorithm>

#define CNP		6+5
#define PNP		3
#define MNP		2
#define FULLQUATSZ	4

namespace bundle {

class PInitating {			//Parallel projection initating
public:
    static constexpr double f = 1;//500000;//200000;		//focal length
    static constexpr double z = 0;			//rate = z/f
    //f = 100000~1000000, z = 10~100
public:
	static void InitSR123T12(double mot[6], double s, double alpha, double beta, double gamma, double t0, double t1)		//alpha, beta, gamma  ~pi
	{
		mot[0] = s;
		mot[1] = alpha/180*M_PI;
		mot[2] = beta/180*M_PI;
		mot[3] = gamma/180*M_PI;
		mot[4] = t0;
		mot[5] = t1;
	}

    static void InverseProjection(const util::_point ori, double P_inv[12], double C[4], double pt3D[3])
    {
        double m[3];
        m[0] = ori.x;
        m[1] = ori.y;
        m[2] = 1;
        double P_invxm[4];
        matrix_product(4, 3, 3, 1, P_inv, m, P_invxm);
        double u = (C[2]-z*C[3])/(z*P_invxm[3]-P_invxm[2]);
        double k = u*P_invxm[3]+C[3];
        pt3D[0] = (u*P_invxm[0]+C[0])/k;
        pt3D[1] = (u*P_invxm[1]+C[1])/k;
        pt3D[2] = z;
// 	std::cout<<"("<<ori.x<<","<<ori.y<<")"<<std::endl;
// 	matrix_print(1, 3, pt3D);
// 	std::cout<<"-------------------------------------------"<<std::endl;
    }

    static void InverseProjectionToPlane(const double m[3], const double P_inv[12], const double C[4], const double abcd[4], double v3[4])
    {
        double P_invxm[4];
        matrix_product(4, 3, 3, 1, P_inv, m, P_invxm);
        double r1, r2, u;
        matrix_product(1, 4, 4, 1, P_invxm, abcd, &r1);
        matrix_product(1, 4, 4, 1, C, abcd, &r2);
        u = -r2/r1;
        v3[0] = P_invxm[0]*u+C[0];
        v3[1] = P_invxm[1]*u+C[1];
        v3[2] = P_invxm[2]*u+C[2];
        v3[3] = P_invxm[3]*u+C[3];
    }

    static util::_point ContraMapping(util::_point ori, int width, int height, float angle)
    {
        int half_width = width/2;
        if(ori.x < half_width) {
            return util::_point(half_width-(half_width-ori.x)/cos(angle/180*M_PI), ori.y);
        }
        else if(ori.x > half_width) {
            return util::_point(half_width+(ori.x-half_width)/cos(angle/180*M_PI), ori.y);
        }
        else {
            return ori;
        }
    }
};


}

#endif
