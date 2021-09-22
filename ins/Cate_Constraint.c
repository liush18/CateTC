/*------------------------------------------------------------------------------
* Cate_Constraint.c : Cate Constraint
*
*          Copyright (C) 2019-2019 by Liu Shuai, All rights reserved.
*
* version :
* history : 2019/09/08  new
*------------------------------------------------------------------------------*/
#include "Cate.h"
#include "State.h"

#define ZVU_R	0.1

/* zero velocity update */
static int zvu_res(const int nv, const int nz, cate_t *cate, double *Cb2e, double *v, double *H, double *R)
{
	int i, j, nx = cate->nx, nv_n = nv;
	double vb[3], OV[9], CV[9], CV_t[9], *R_;

	matmul("TN", 3, 1, 3, 1.0, Cb2e, cate->solv.ve, 0.0, vb);
	Mat_Skew(cate->solv.ve, OV);			/* skew symmetric matrix of ve */
	matmul("TN", 3, 3, 3, 1.0, Cb2e, OV, 0.0, CV); /* from ecef to body */
	Mat_Trans(CV, 3, 3, CV_t);				/* transposition */

	//printf("\n\nCb2e:\n");
	//matprint(Cb2e, 3, 3, 10, 7);
	//printf("\n\nCV_t:\n");
	//matprint(CV_t, 3, 3, 10, 7);

	for (i = 0; i < 3; i++) {
		H[IIV + i + nx * nv_n] = Cb2e[i];
		H[IIA + i + nx * nv_n] = CV_t[i];
	}
	v[nv_n++] = -vb[0];

	if (nz == 3) {
		for (i = 0; i < 3; i++) {
			H[IIV + i + nx * nv_n] = Cb2e[i + 3];
			H[IIA + i + nx * nv_n] = CV_t[i + 3];
		}
		v[nv_n++] = -vb[1];
	}

	for (i = 0; i < 3; i++) {
		H[IIV + i + nx * nv_n] = Cb2e[i + 6];
		H[IIA + i + nx * nv_n] = CV_t[i + 6];
	}
	v[nv_n++] = -vb[2];

	R_ = zeros(nv_n, nv_n);
	Mat_Sub(R_, nv_n, R, nv, 0, 0, 1.0);
	matcpy(R, R_, nv_n, nv_n);

	for (i = nv; i < nv_n; i++) {
		for (j = nv; j < nv_n; j++) {
			R[j + i * nv_n] = i == j ? SQR(ZVU_R) : 0.0;
		}
	}

	free(R_);
	return nv_n;
}

/* 0:static, 1:kinematic */
static int detec_static(const cateopt_t *opt, const double *ve)
{

}

/* free ins zero velocity update for xb and zb */
extern int Constraint(cate_t *cate, double *Cb2e, double *v, double *H, double *R, const int nv)
{
	int nx = cate->nx, nv_n, nz = 2;

	if (norm(cate->solv.ve, 3) < cate->opt.maxvel) nz = 3;

	nv_n = zvu_res(nv, nz, cate, Cb2e, v, H, R);

	return nv_n;
}