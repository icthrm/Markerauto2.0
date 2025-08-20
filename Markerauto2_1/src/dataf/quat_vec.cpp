#include "quat_vec.h"
#include <cmath>
#include <iostream>
#include <string.h>

void Quat2Vec(const double* inp, int nin, double* outp, int nout)
{
    double mag, sg;
    register int i;

    /* intrinsics & distortion */
    if(nin>7){ 			// are they present?
        for(i=0; i<nin-7; ++i){
            outp[i]=inp[i];
	}
    }
    else{
        i=0;
    }

    /* rotation */
    /* normalize and ensure that the quaternion's scalar component is non-negative;
     * if not, negate the quaternion since two quaternions q and -q represent the
     * same rotation
     */
    mag=sqrt(inp[i]*inp[i] + inp[i+1]*inp[i+1] + inp[i+2]*inp[i+2] + inp[i+3]*inp[i+3]);
    sg=(inp[i]>=0.0)? 1.0 : -1.0;
    mag=sg/mag;
    outp[i]  =inp[i+1]*mag;
    outp[i+1]=inp[i+2]*mag;
    outp[i+2]=inp[i+3]*mag;
    i+=3;

    /* translation*/
    for( ; i<nout; ++i){
        outp[i]=inp[i+1];
    }
}

void Vec2Quat(const double* inp, int nin, double* outp, int nout)
{
    const double *v;
    double q[FULLQUATSZ];
    register int i;

    /* intrinsics & distortion */
    if(nin>7-1){		 // are they present?
        for(i=0; i<nin-(7-1); ++i){
            outp[i]=inp[i];
	}
    }
    else{
        i=0;
    }

    /* rotation */
    /* recover the quaternion from the vector */
    v=inp+i;
    _MK_QUAT_FRM_VEC(q, v);
    outp[i]  =q[0];
    outp[i+1]=q[1];
    outp[i+2]=q[2];
    outp[i+3]=q[3];
    i+=FULLQUATSZ;

    /* translation */
    for ( ; i<nout; ++i){
        outp[i]=inp[i-1];
    }
}

void QuatMultiply(const double q1[FULLQUATSZ], const double q2[FULLQUATSZ], double p[FULLQUATSZ])
{
double t1, t2, t3, t4, t5, t6, t7, t8, t9;
//double t10, t11, t12;

  t1=(q1[0]+q1[1])*(q2[0]+q2[1]);
  t2=(q1[3]-q1[2])*(q2[2]-q2[3]);
  t3=(q1[1]-q1[0])*(q2[2]+q2[3]);
  t4=(q1[2]+q1[3])*(q2[1]-q2[0]);
  t5=(q1[1]+q1[3])*(q2[1]+q2[2]);
  t6=(q1[1]-q1[3])*(q2[1]-q2[2]);
  t7=(q1[0]+q1[2])*(q2[0]-q2[3]);
  t8=(q1[0]-q1[2])*(q2[0]+q2[3]);

#if 0
  t9 =t5+t6;
  t10=t7+t8;
  t11=t5-t6;
  t12=t7-t8;

  p[0]= t2 + 0.5*(-t9+t10);
  p[1]= t1 - 0.5*(t9+t10);
  p[2]=-t3 + 0.5*(t11+t12);
  p[3]=-t4 + 0.5*(t11-t12);
#endif

  /* following fragment it equivalent to the one above */
  t9=0.5*(t5-t6+t7+t8);
  p[0]= t2 + t9-t5;
  p[1]= t1 - t9-t6;
  p[2]=-t3 + t9-t8;
  p[3]=-t4 + t9-t7;
}

void RotMatrix2Quat(const double R[9], double quat[FULLQUATSZ])
{
    double tmp[FULLQUATSZ];
    int maxpos = -1;
    double mag = -1.0;
    // find the maximum of the 4 quantities
    tmp[0] = 1.0+R[0]+R[4]+R[8];
    tmp[1] = 1.0+R[0]-R[4]-R[8];
    tmp[2] = 1.0-R[0]+R[4]-R[8];
    tmp[3] = 1.0-R[0]-R[4]+R[8];
    
    for(int i = 0; i < 4; i++){
        if (tmp[i] > mag) {
            mag = tmp[i];
            maxpos = i;
        }
    }

    if(maxpos == 0){
        quat[0] = sqrt(tmp[0])*0.5;
        quat[1] = (R[7]-R[5])/(4.0*quat[0]);
        quat[2] = (R[2]-R[6])/(4.0*quat[0]);
        quat[3] = (R[3]-R[1])/(4.0*quat[0]);
    }
    else if(maxpos == 1){
        quat[1] = sqrt(tmp[1])*0.5;
        quat[0] = (R[7]-R[5])/(4.0*quat[1]);
        quat[2] = (R[3]+R[1])/(4.0*quat[1]);
        quat[3] = (R[2]+R[6])/(4.0*quat[1]);
    }
    else if(maxpos == 2){
        quat[2] = sqrt(tmp[2])*0.5;
        quat[0] = (R[2]-R[6])/(4.0*quat[2]);
        quat[1] = (R[3]+R[1])/(4.0*quat[2]);
        quat[3] = (R[7]+R[5])/(4.0*quat[2]);
    }
    else if(maxpos == 3){
        quat[3] = sqrt(tmp[3])*0.5;
        quat[0] = (R[3]-R[1])/(4.0*quat[3]);
        quat[1] = (R[2]+R[6])/(4.0*quat[3]);
        quat[2] = (R[7]+R[5])/(4.0*quat[3]);
    }

    // enforce unit length
    mag=quat[0]*quat[0]+quat[1]*quat[1]+quat[2]*quat[2]+quat[3]*quat[3];
    mag = 1.0/sqrt(mag);
    quat[0] *= mag;
    quat[1] *= mag;
    quat[2] *= mag;
    quat[3] *= mag;
}

void Quat2RotMatrix(const double quat[FULLQUATSZ], double R[9])
{
    double mag;
    // ensure unit length
    double tmp[FULLQUATSZ];
    memcpy(tmp, quat, sizeof(double)*FULLQUATSZ);
    mag=tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]+tmp[3]*tmp[3];
    if(mag!=1.0){
        mag=1.0/sqrt(mag);
        tmp[0]*=mag;
        tmp[1]*=mag;
        tmp[2]*=mag;
        tmp[3]*=mag;
    }

    R[0]=tmp[0]*tmp[0]+tmp[1]*tmp[1]-tmp[2]*tmp[2]-tmp[3]*tmp[3];
    R[1]=2*(tmp[1]*tmp[2]-tmp[0]*tmp[3]);
    R[2]=2*(tmp[1]*tmp[3]+tmp[0]*tmp[2]);

    R[3]=2*(tmp[1]*tmp[2]+tmp[0]*tmp[3]);
    R[4]=tmp[0]*tmp[0]+tmp[2]*tmp[2]-tmp[1]*tmp[1]-tmp[3]*tmp[3];
    R[5]=2*(tmp[2]*tmp[3]-tmp[0]*tmp[1]);

    R[6]=2*(tmp[1]*tmp[3]-tmp[0]*tmp[2]);
    R[7]=2*(tmp[2]*tmp[3]+tmp[0]*tmp[1]);
    R[8]=tmp[0]*tmp[0]+tmp[3]*tmp[3]-tmp[1]*tmp[1]-tmp[2]*tmp[2];
}