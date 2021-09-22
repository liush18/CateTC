/*------------------------------------------------------------------------------
* Cate_NavFins.c : free ins
*
*          Copyright (C) 2019-2019 by Liu Shuai, All rights reserved.
*
* version :
* history : 2019/04/18  new
*------------------------------------------------------------------------------*/
#include "Cate.h"
#include "State.h"

static int revs = 0;		/* 0:forward, 1:backward */
static int nimu = 0;
static int isolvf = 0;
static int isolvb = 0;
static solvd_t *solvf;
static solvd_t *solvb;

static void ins_sol(imud_t *data, double *Cb2e, const double tor_i, solvd_t *solv)
{
	int i;
	double w;
	double *wie2e, *wie2b, *wib2b, *web2b, *fib2e, *ve, *ge, *pos;
	double *old_Cb2e, *W, *W2, *A, *ave_Cb2e, *Wie2e;					/* W = Web2b , W2 = (Web2b)^2 */

	/* attitude update */
	old_Cb2e = mat(3, 3); W = mat(3, 3); W2 = mat(3, 3); A = eye(3);
	wie2e = zeros(3, 1); wie2b = mat(3, 1); wib2b = mat(3, 1); web2b = mat(3, 1);

	matcpy(old_Cb2e, Cb2e, 3, 3);
	wie2e[2] = OMGE * tor_i;
	matmul("TN", 3, 1, 3, 1.0, old_Cb2e, wie2e, 0.0, wie2b);	/* wie2b = Cb2e'*wie2e */
	for (i = 0; i < 3; i++) wib2b[i] = data->gyro[i] * tor_i;
	Mat_Add(3, 1, 1.0, wib2b, -1.0, wie2b, 0.0, web2b);		/* web2b = wib2b - wie2b */

	w = norm(web2b, 3);
	Mat_Skew(web2b, W);
	matmul("NN", 3, 3, 3, 1.0, W, W, 0.0, W2);
	Mat_Add(3, 3, sin(w) / w, W, (1 - cos(w)) / SQR(w), W2, 1.0, A); /* A = I + sinw/w*W + (1-cosw)/w/w*W*W */
	matmul("NN", 3, 3, 3, 1.0, old_Cb2e, A, 0.0, Cb2e);

	/* velocity update */
	fib2e = mat(3, 1); ave_Cb2e = mat(3, 3); ve = mat(3, 1); ge = mat(3, 1); Wie2e = mat(3, 3);
	Mat_Add(3, 3, 0.5, Cb2e, 0.5, old_Cb2e, 0.0, ave_Cb2e);
	matmul("NN", 3, 1, 3, 1.0, ave_Cb2e, data->acce, 0.0, fib2e); /* f_ib_e = (Cb2e+old_Cb2e)*0.5*f_ib_b */

	Grav_ECEF(solv->re, ge);
	Mat_Skew(wie2e, Wie2e);
	matmul("NN", 3, 1, 3, -2.0, Wie2e, solv->ve, 0.0, ve);
	Mat_Add(3, 1, tor_i, fib2e, tor_i, ge, 1.0, ve);
	Mat_Add(3, 1, 1.0, solv->ve, 1.0, ve, 0.0, ve);			/* ve = old_ve + (f_ib_e + ge - 2*Omega_ie*old_ve)*tor_i */

	/* position update */
	pos = mat(3, 1);
	Mat_Add(3, 1, 0.5*tor_i, solv->ve, 0.5*tor_i, ve, 1.0, solv->re);	/* re = old_re + (old_ve+ve)*tor_i*0.5 （5.38）*/

	/* get attitude angle */
	ecef2pos(solv->re, pos);
	Cb2e_to_Att(Cb2e, pos, solv->att);
	matcpy(solv->ve, ve, 3, 1);

	free(wie2e); free(wie2b); free(wib2b); free(web2b); free(fib2e); free(ve); free(ge); free(pos);
	free(old_Cb2e); free(W); free(W2); free(A); free(ave_Cb2e); free(Wie2e);

}



static int const_fins(cate_t *cate, double *Cb2e, double *acce, double *gyro)
{
	int nx = cate->nx, nv = 0;
	double *v, *H, *R;

	v = zeros(6, 1); H = zeros(nx, 6); R = zeros(6, 6);

	UdState_INS(cate, cate->P, Cb2e, acce, gyro, nx);
	nv = Constraint(cate, Cb2e, v, H, R, nv);

	Kalman(cate->xx, cate->P, v, H, R, nx, nv);

	Udsol(NULL, NULL, 0, cate, Cb2e);

	free(v); free(H); free(R);
	return 1;
}

/*
 * rev ( 0:forward,			1:backward	)
 * flag( 0:no constraint,	1:constraint)
*/
extern void FinsPos(cate_t *cate, imud_t *data, double *Cb2e, const double tor_i, const int flag)
{

	//if (cate->opt.coupled_type > CSOL_PINS && LEVERCORR)
	//	LeverG2I(cate->solv.re, cate->solv.ve, Cb2e, data->gyro, cate->opt.lever);
	ins_sol(data, Cb2e, tor_i, &cate->solv);
	//if (cate->opt.coupled_type > CSOL_PINS && LEVERCORR)
	//	LeverI2G(cate->solv.re, cate->solv.ve, Cb2e, data->gyro, cate->opt.lever);

	cate->solv.stat = CSOL_PINS;

	if (flag && cate->opt.zvu) {
		const_fins(cate, Cb2e, data->acce, data->gyro);
	}
}

#if 0
/*
 * mode (0:forward/backward, 1:combined ) 
*/
static void pins_proc(cate_t *cate, double *Cb2e, imu_t *imusol, int mode)
{
	imud_t imu = { 0 };
	double tHz = 1.0 / cate->opt.Hz;

	while (InputIMU(imusol, &imu, revs, &nimu, cate->solv.stat, 1)) {

		if (cate->opt.corr_imu) { Corr_IMUd(&imu, &cate->solv, cate->opt.corr_imu_k); }

		if (!revs) {	/* forward */
			cate->solv.time = imu.time;
			FinsPos(cate, &imu, Cb2e, tHz, revs, 1);
		}
		else {	/* backward */
			if (nimu > 0) {
				cate->solv.time = (imusol->data + nimu)->time;
			}
			FinsPos(cate, &imu, Cb2e, tHz, revs, 0);
		}

		cate->solv.stat = CSOL_PINS;

		if (mode == 0) {
			cOutSols(&cate->solv);
		}
		else if (!revs) {	/* combined-forward */
			solvf[isolvf++] = cate->solv;
		}
		else {	/* combined-backward */
			solvb[isolvb++] = cate->solv;
		}
	}
}

static void comb_pins(cate_t *cate)
{
	solvd_t solvs = { {0} };
	int i, j, k;
	double tHz2 = 1.0 / cate->opt.Hz / 2, tt;

	for (i = 0, j = isolvb - 1; i < isolvf&& j >= 0; i++, j--) {

		/* 前向滤波历元在后向滤波历元之前 */
		if ((tt = timediff(solvf[i].time, solvb[j].time)) < -tHz2) {
			solvs = solvf[i];
			j++;/* 只使用了前向滤波的i历元，没使用后向滤波的j历元，而下一次for会j--，
				相当于跳过没使用的j历元了，为了不跳过，这里j++，下面的i--相同道理 */
		}
		else if (tt > tHz2) { /* 后向滤波历元在前向滤波历元之前 */
			solvs = solvb[j];
			i--;
		}
		else if (solvf[i].stat < solvb[j].stat) { /* 两个滤波时间间隔在允许范围内时 */
			solvs = solvf[i];
		}
		else if (solvf[i].stat > solvb[j].stat) {
			solvs = solvb[j];
		}
		else {
			solvs = solvf[i];	/* 能运行到这里，说明两个滤波方式除了结果不一致，其他一样，随便一种滤波的方式初始化solvs  */
			solvs.time = timeadd(solvs.time, -tt / 2.0);

			/* 直接平均的方式平滑 */
			for (k = 0; k < 3; k++) {
				solvs.re[k] = 0.5*(solvf[i].re[k] + solvb[j].re[k]);
				solvs.ve[k] = 0.5*(solvf[i].ve[k] + solvb[j].ve[k]);
				solvs.att[k] = 0.5*(solvf[i].att[k] + solvb[j].att[k]);
			}
			cate->solv = solvs;
		}

		cOutSols(&solvs);
	}
}

/* 
 * 沌口前向pos=[-2252179.2756,5018496.1117,3217680.8927],
 *			vel=[0.0,0.0,0.0],
 *			att=[-0.0476477,0.0205796,-0.1250026];
 * 沌口后向(纯惯导模式)pos=[-2253336.1727,5017517.1391,3215729.9247],
 *			vel=[-2.2935,4.4881,-16.3500],
 *			att=[-0.0367802,0.0094927,2.8980073]
 */
extern int PinsExec(cate_t *cate, double *Cb2e, imu_t *imusol)
{
	revs = nimu = 0;

	if (cate->popt.soltype == 0) {		/* forward */
		pins_proc(cate, Cb2e, imusol, 0);
		//printf("\n");matprint(cate->solv.re, 3, 1, 15, 4);
		//printf("\n"); matprint(cate->solv.ve, 3, 1, 15, 4);
		//double att[3], pos[3];
		//ecef2pos(cate->solv.re, pos);
		//Cb2e_to_Att(Cb2e, pos, att);
		//printf("\n"); matprint(att, 3, 1, 12, 7);
	}
	else if (cate->popt.soltype == 1) {	/* backward */
		revs = 1; nimu = imusol->n - 1;
		pins_proc(cate, Cb2e, imusol, 0);
	}
	else {	/* combined */
		solvf = (solvd_t *)malloc(sizeof(solvd_t)*imusol->n);
		solvb = (solvd_t *)malloc(sizeof(solvd_t)*imusol->n);

		/* 前后向平滑，考虑不用前向的结果作为后向的初始化，而是要重新后向对准初始化 */
		if (solvf&&solvb) {
			isolvf = isolvb = 0;
			pins_proc(cate, Cb2e, imusol, 1);	/* forward */
			revs = 1; nimu = imusol->n - 1;
			pins_proc(cate, Cb2e, imusol, 1);	/* backward */
		}
		else {
			follow(1, "PinsExec memory allocation error\n");
		}
		comb_pins(cate);
		free(solvf);
		free(solvb);
	}

	return 1;
}

#endif