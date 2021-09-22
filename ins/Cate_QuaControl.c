/* ------------------------------------------------------------------------------
 * Cate_QuaControl.c : Quality Control
 *
 *			Copyright (C) 2019-2019 by Liu Shuai, All rights reserved.
 *
 * version : 
 * history : 2019/09/22  new
 *			 2019/10/15  add ObsMeas_Detec
 * ------------------------------------------------------------------------------ */

#include "Cate.h"

/* gf measurements */
static double Meas_GF(const obsd_t *obs, const nav_t *nav)
{
	const double *lam = nav->lam[obs->sat - 1];

	if (lam[0] == 0.0 || lam[1] == 0.0 || obs->L[0] == 0.0 || obs->L[1] == 0.0) return 0.0;

	return lam[0] * obs->L[0] - lam[1] * obs->L[1];
}

/* cycle slip detection by gf */
static void CycSlip_Detec_GF(cate_t *cate, rtk_t *rtk, const obsd_t *obs, const int n, const nav_t *nav)
{
	int i, j, sat;
	double g0, g1, dE, el, thres;

	dE = fabs(cate->tt / cate->opt.si);		/* delta epoch */

	for (i = 0; i < n&&i < MAXOBS; i++) {
		sat = obs[i].sat;

		el = rtk->ssat[sat - 1].azel[1] * R2D;/* degree */

		/* cycle slip threshold according sample interval and satellite elevation angle */
		if (el > 15.0) thres = cate->copt.gft;
		else thres = (-1.0 / 15.0*el + 2.0)*cate->copt.gft;

		if ((g0 = Meas_GF(obs + i, nav)) == 0.0) continue;	/* current epoch gf measurements */
		g1 = cate->csat[sat - 1].gf;

		if (g1 != 0.0 && fabs(g0 - g1) > MIN(thres*dE, 1.5)) {
			for (j = 0; j < cate->popt.nf; j++)
				rtk->ssat[sat - 1].slip[j] = 1;
		}

		cate->csat[sat-1].gf = g0;
	}

}

/* mw measurements */
static double Meas_MW(const obsd_t *obs, const nav_t *nav)
{
	int i = 0, j = 1;
	const double *lam = nav->lam[obs->sat - 1];
	double P1, P2, P1_C1, P2_C2, lam1, lam2, res;

	if (obs->L[i] == 0.0) return 0.0;
	if (obs->L[j] == 0.0) return 0.0;
	if (obs->P[i] == 0.0) return 0.0;
	if (obs->P[j] == 0.0) return 0.0;
	if (lam[i] * lam[j] == 0.0) return 0.0;

	P1 = obs->P[i];
	P2 = obs->P[j];
	P1_C1 = nav->cbias[obs->sat - 1][1];
	P2_C2 = nav->cbias[obs->sat - 1][2];

	if (obs->code[0] == CODE_L1C) P1 += P1_C1; /* C1->P1 */
	if (obs->code[1] == CODE_L2C) P2 += P2_C2; /* C2->P2 */

	lam1 = lam[i];
	lam2 = lam[j];
	res = (obs->L[i] - obs->L[j]) - (lam2 - lam1) / (lam1 + lam2)*(P1 / lam1 + P2 / lam2);

	return res;
}

/* cycle slip detection by mw */
static void CycSlip_Detec_MW(cate_t *cate, rtk_t *rtk, const obsd_t *obs, const int n, const nav_t *nav)
{
	int i, j, sat;
	double w0, w1, el, thres;

	for (i = 0; i < n&&i < MAXOBS; i++) {
		sat = obs[i].sat;

		if (timediff(obs->time, cate->csat[sat - 1].pepoch) > cate->opt.si) {
			cate->csat[sat - 1].mw = 0.0;
			cate->csat[sat - 1].imw = 0;
		}

		el = rtk->ssat[sat - 1].azel[1] * R2D;

		/* cycle slip threshold according sample interval and satellite elevation angle */
		if (el > 20.0) thres = cate->copt.mwt;
		else thres = (-0.1*el + 3.0)*cate->copt.mwt;

		if ((w0 = Meas_MW(obs + i, nav)) == 0.0) continue;	/* current epoch gf measurements */
		w1 = cate->csat[sat - 1].mw;

		if (w1 != 0.0 && fabs(w0 - w1) > MIN(thres, 6.0)) {
			for (j = 0; j < cate->popt.nf; j++)
				rtk->ssat[sat - 1].slip[j] = 1;
		}

		if ((j=cate->csat[sat-1].imw) > 0) {
			cate->csat[sat - 1].mw = (w1*j + w0) / (j + 1);
			cate->csat[sat - 1].imw++;
		}
		else {
			cate->csat[sat - 1].mw = w0;
			cate->csat[sat - 1].imw++;
		}

		cate->csat[sat - 1].pepoch = obs[i].time;
	}

}

/* calculate cycle slip threshold according gnss sample interval */
extern void CycSlip_Thres(const double samp, double *gft, double *mwt)
{
	follow(2, "CycSlip_Thres...\n");

	if (samp > 0.0) {

		if (samp <= 1.0) {
			*gft = 0.05;
			*mwt = 2.5;
		}
		else if (samp <= 20) {
			*gft = (0.1) / (20.0)*samp + 0.05;
			*mwt = (2.5) / (20.0)*samp + 2.5;
		}
		else if (samp <= 60) {
			*gft = 0.15;
			*mwt = 5.0;
		}
		else if (samp <= 100) {
			*gft = 0.25;
			*mwt = 7.5;
		}
		else {
			*gft = 0.35;
			*mwt = 7.5;
		}
	}
}

/* cycle slip detection by gf and mw */
extern void CycSlip_Detec(cate_t *cate, rtk_t *rtk, const obsd_t *obs, const int n, const nav_t *nav )
{
	int i, j;

	for (i = 0; i < SATNUM; i++) {
		for (j = 0; j < cate->popt.nf; j++)
			rtk->ssat[i].slip[j] = 0;
	}
		
	CycSlip_Detec_GF(cate, rtk, obs, n, nav);
	CycSlip_Detec_MW(cate, rtk, obs, n, nav);
}

/* obs measurements detection
 * nobs		: number of input obs
 * ref : gamp
*/
extern int ObsMeasSPP_Detec(const prcopt_t *opt, obsd_t *obs, const int nobs)
{
	double tP = 0.0;
	int i, j, n, sat, sys;

	for (i = n = 0; i < nobs&&i < MAXOBS; i++) {
		sat = obs[i].sat;
		sys = satsys(sat, NULL);

		/* exclude satellites */
		if (!(sys&opt->navsys)) continue;
		if (opt->exsats[sat - 1] == 1)	 continue;

		for (j = 0; j < NFREQ; j++) {
			tP += SQR(obs[i].P[j]);
		}
		if (tP < 1.0) continue;

		obs[n++] = obs[i];
	}

	return n; /* number of valid obs */
}

/* obs measurements detection
 * nobs		: number of input obs
 * ref : gamp
*/
extern int ObsMeasPPP_Detec(obsd_t *obs, const int nobs)
{
	int i, n, sat;

	for (i = n = 0; i < nobs&&i < MAXOBS; i++) {
		sat = obs[i].sat;

		if (obs[i].L[0] * obs[i].L[1] == 0.0) continue;
		if (fabs(obs[i].P[0] - obs[i].P[1]) >= 200.0)continue;

		obs[n++] = obs[i];
	}

	return n; /* number of valid obs */
}

/* doppler meas detection */
extern int ObsMeasDOP_Detec(const obsd_t *obs, const int nobs)
{
	int i, n, sat;

	for (i = n = 0; i < nobs&&i < MAXOBS; i++) {
		sat = obs[i].sat;

		if (obs[i].D[0] == 0.0) continue;
		n++;
	}

	return (n == nobs) ? 1 : 0;
}


extern void LeverI2G(double *re, double *ve, const double *Cb2e, const double *gyro, const double *lever)
{
	double cl[3], Wib2b[9], CW[9], cwl[3];
	double wie2e[3] = { 0,0,OMGE }, Wie2e[9], WC[9], wcl[3];

	/* rGe = rIe + Cb2e lb */
	matmul("NN", 3, 1, 3, 1.0, Cb2e, lever, 0.0, cl);
	Mat_Add(3, 1, 1.0, re, 1.0, cl, 0.0, re);

	/* vGe = vIe + Cb2e Wib2b lb - Wie2e Cb2e lb */
	Mat_Skew(gyro, Wib2b);
	matmul("NN", 3, 3, 3, 1.0, Cb2e, Wib2b, 0.0, CW);
	matmul("NN", 3, 1, 3, 1.0, CW, lever, 0.0, cwl);

	Mat_Skew(wie2e, Wie2e);
	matmul("NN", 3, 3, 3, 1.0, Wie2e, Cb2e, 0.0, WC);
	matmul("NN", 3, 1, 3, 1.0, WC, lever, 0.0, wcl);

	Mat_Add(3, 1, 1.0, cwl, -1.0, wcl, 1.0, ve);
}

extern void LeverG2I(double *re, double *ve, const double *Cb2e, const double *gyro, const double *lever)
{
	double cl[3], Wib2b[9], CW[9], cwl[3];
	double wie2e[3] = { 0,0,OMGE }, Wie2e[9], WC[9], wcl[3];

	/* rIe = rGe - Cb2e lb */
	matmul("NN", 3, 1, 3, 1.0, Cb2e, lever, 0.0, cl);
	Mat_Add(3, 1, 1.0, re, -1.0, cl, 0.0, re);

	/* vIe = vGe - Cb2e Wib2b lb + Wie2e Cb2e lb */
	Mat_Skew(gyro, Wib2b);
	matmul("NN", 3, 3, 3, 1.0, Cb2e, Wib2b, 0.0, CW);
	matmul("NN", 3, 1, 3, 1.0, CW, lever, 0.0, cwl);

	Mat_Skew(wie2e, Wie2e);
	matmul("NN", 3, 3, 3, 1.0, Wie2e, Cb2e, 0.0, WC);
	matmul("NN", 3, 1, 3, 1.0, WC, lever, 0.0, wcl);

	Mat_Add(3, 1, -1.0, cwl, 1.0, wcl, 1.0, ve);
}