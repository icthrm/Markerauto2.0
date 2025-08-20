#ifndef IMG_PROJS_H__
#define IMG_PROJS_H__

#include "sba.h"
#include <iostream>


/** pointers to additional data, used for computed image projections and their jacobians */
struct pglobs{
	bool cpmask[6];		//mask of camera params: s, alpha, beta, gama, t0, t1;		if mask = 1; the jac will be masked;
	bool pmask;			//mask of points
};	


void calcImgParallelProj(const double param[6], const double M[3], double n[2]);
void calccalcImgParallelProjJac(const double param[6], const double M[3], double jacmKRT[2][11], double jacmS[2][3]);


void img_ParallelProj_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata);
void img_ParallelProj_jac_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata);

#endif
