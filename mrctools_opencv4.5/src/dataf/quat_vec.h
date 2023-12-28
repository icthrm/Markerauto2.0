#ifndef QUAT_VEC_H__
#define QUAT_VEC_H__

#define FULLQUATSZ	4

/** unit quaternion from vector part */
#define _MK_QUAT_FRM_VEC(q, v){                                     \
  (q)[1]=(v)[0]; (q)[2]=(v)[1]; (q)[3]=(v)[2];                      \
  (q)[0]=sqrt(1.0 - (q)[1]*(q)[1] - (q)[2]*(q)[2]- (q)[3]*(q)[3]);  \
}

/**convert a vector of camera parameters so that rotation is represented by
 * the vector part of the input quaternion. The function converts the
 * input quaternion into a unit one with a non-negative scalar part. Remaining
 * parameters are left unchanged.
 *
 * Input parameter layout: intrinsics (5, optional), distortion (5, optional), rot. quaternion (4), translation (3)
 * Output parameter layout: intrinsics (5, optional), distortion (5, optional), rot. quaternion vector part (3), translation (3)
 */
void Quat2Vec(const double *inp, int nin, double *outp, int nout);

/**convert a vector of camera parameters so that rotation is represented by
 * a full unit quaternion instead of its input 3-vector part. Remaining
 * parameters are left unchanged.
 *
 * Input parameter layout: intrinsics (5, optional), distortion (5, optional), rot. quaternion vector part (3), translation (3)
 * Output parameter layout: intrinsics (5, optional), distortion (5, optional), rot. quaternion (4), translation (3)
 */
void Vec2Quat(const double *inp, int nin, double *outp, int nout);

/**fast multiplication of the two quaternions in q1 and q2 into p
 * this is the second of the two schemes derived in pg. 8 of
 * T. D. Howell, J.-C. Lafon, The complexity of the quaternion product, TR 75-245, Cornell Univ., June 1975.
 *
 * total additions increase from 12 to 27 (28), but total multiplications decrease from 16 to 9 (12)
 */
void QuatMultiply(const double q1[FULLQUATSZ], const double q2[FULLQUATSZ], double p[FULLQUATSZ]);

/** compute the quaternion corresponding to a rotation matrix; see A8 in Horn's paper */
void RotMatrix2Quat(const double R[9], double quat[FULLQUATSZ]);

/** compute the rotation matrix corresponding to a quaternion; see Horn's paper */
void Quat2RotMatrix(const double quat[FULLQUATSZ], double R[9]);


#endif