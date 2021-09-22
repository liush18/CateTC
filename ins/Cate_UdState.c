/*------------------------------------------------------------------------------
 * Cate_UdState.c : Update State
 *
 *			Copyright (C) 2019-2019 by Liu Shuai, All rights reserved.
 *
 * version :
 * history : 2019/10/31  new
 *------------------------------------------------------------------------------*/
#include "Cate.h"
#include "State.h"

#define ROUND(x)		(int)floor((x)+0.5)
#define GAP_RESION		120					/* default gap to reset ionos parameters (ep) */

extern int lscpnx(const cateopt_t *opt)
{
	return CI_N(opt);
}

extern int ttcpnx(const cateopt_t *opt, const prcopt_t *popt)
{
	return CI_N(opt) + CG_N(popt);
}

static void udins_px(cate_t *cate, const cateopt_t *opt, const solvd_t *solv)
{
	int i;

	if (norm(cate->xx, 3) == 0.0) {
		for (i = 0; i < 3; i++) {
			Init_Px(cate, solv->re[i], CVAR_POS, IIP + i);
			Init_Px(cate, solv->ve[i], CVAR_VEL, IIV + i);
			Init_Px(cate, solv->att[i], CVAR_ATT, IIA + i);
			Init_Px(cate, opt->ba, CVAR_BA, IIBA + i);
			Init_Px(cate, opt->bg, CVAR_BG, IIBG + i);

			if (opt->est_scale) {
				Init_Px(cate, opt->sa, CVAR_SA, IISA(opt) + i);
				Init_Px(cate, opt->sg, CVAR_SG, IISG(opt) + i);
			}
		}
	}
	else {
		for (i = 0; i < 3; i++) {
			Init_Px(cate, solv->re[i], CVAR_POS, IIP + i);
			Init_Px(cate, solv->ve[i], CVAR_VEL, IIV + i);
			Init_Px(cate, solv->att[i], CVAR_ATT, IIA + i);
		}
	}
}

static void udins_fg(const cateopt_t *opt, const double *Cb2e, const double *acce, const double *gyro,
	const int nx, double *F, double *G)
{
	int i;
	double *I, *IG, alpha_ie[3] = { 0 }, Alpha_ie[9];
	double CF[9], temp_v[3];
	double Fa[9], Fg[9], CFa[9], CFg[9];
	double Rho[9] = { 0 };

	for (i = 0; i < nx*nx; i++) {
		G[i] = 0;
	}

	I = eye(3); IG = eye(3);

	/* matrix of earth rotation angular velocity */
	alpha_ie[2] = OMGE;
	Mat_Skew(alpha_ie, Alpha_ie);

	/* CF = (Cb2e * fib2b)× */
	matmul("NN", 3, 1, 3, 1.0, Cb2e, acce, 0.0, temp_v);
	Mat_Skew(temp_v, CF);

	/* CFa = Cb2e * Fa, CFg = Cb2e * Fg */
	Vec_to_DiagMat(acce, Fa, 3);
	Vec_to_DiagMat(gyro, Fg, 3);
	matmul("NN", 3, 3, 3, 1.0, Cb2e, Fa, 0.0, CFa);
	matmul("NN", 3, 3, 3, 1.0, Cb2e, Fg, 0.0, CFg);

	/* correlation time */
	if (opt->markov)	Rho[0] = Rho[4] = Rho[8] = 1 / opt->corretime_ag;

	/* system state matrix F */
	/* col.2 vel */
	Mat_Sub(F, nx, I, 3, IIP, IIV, 1.0);				/* I */
	Mat_Sub(F, nx, Alpha_ie, 3, IIV, IIV, -2.0);		/* -2 \Omega_{ie}^e */
	/* col.3 att */
	Mat_Sub(F, nx, CF, 3, IIV, IIA, -1.0);			/* -(C_b^e * f_{ib}^b)\times */
	Mat_Sub(F, nx, Alpha_ie, 3, IIA, IIA, -1.0);		/* -\Omega_{ie}^e */
	/* col.4 ba */
	Mat_Sub(F, nx, Cb2e, 3, IIV, IIBA, 1.0);			/* C_b^e */
	Mat_Sub(F, nx, Rho, 3, IIBA, IIBA, -1.0);				/* -Rho */
	/* col.5 bg */
	Mat_Sub(F, nx, Cb2e, 3, IIA, IIBG, 1.0);			/* C_b^e */
	Mat_Sub(F, nx, Rho, 3, IIBG, IIBG, -1.0);				/* -Rho */

	/* system noise matrix G */
	/* col.1 wa */
	Mat_Sub(G, nx, Cb2e, 3, IIV, IIV, 1.0);
	/* col.2 wg */
	Mat_Sub(G, nx, Cb2e, 3, IIA, IIA, 1.0);
	/* col.3 wba */
	Mat_Sub(G, nx, IG, 3, IIBA, IIBA, 1.0);
	/* col.4 wbg */
	Mat_Sub(G, nx, IG, 3, IIBG, IIBG, 1.0);

	if (opt->est_scale) {
		/* F */
	/* col.6 sa */
		Mat_Sub(F, nx, CFa, 3, IIV, IISA(opt), 1.0);				/* C_b^e F_a */
		Mat_Sub(F, nx, Rho, 3, IISA(opt), IISA(opt), -1.0);			/* -Rho */
		/* col.7 sg */
		Mat_Sub(F, nx, CFg, 3, IIA, IISG(opt), 1.0);				/* C_b^e F_g */
		Mat_Sub(F, nx, Rho, 3, IISG(opt), IISG(opt), -1.0);			/* -Rho */

		/* G */
		/* col.5 wsa */
		Mat_Sub(G, nx, IG, 3, IISA(opt), IISA(opt), 1.0);
		/* col.6 wsg */
		Mat_Sub(G, nx, IG, 3, IISG(opt), IISG(opt), 1.0);
	}

	free(I); free(IG);
}

static void udins_q(const cateopt_t *opt, const int nx, const double t, double *Q)
{
	int i;
	double tt = fabs(t);
	cprn_t prn = opt->prn;

	for (i = 0; i < 3; i++) {
		Init_SqMat(nx, Q, SQR(prn.wa)*tt, IIV + i);
		Init_SqMat(nx, Q, SQR(prn.wg)*tt, IIA + i);
		Init_SqMat(nx, Q, SQR(prn.wba)*tt, IIBA + i);
		Init_SqMat(nx, Q, SQR(prn.wbg)*tt, IIBG + i);

		if (opt->est_scale) {
			Init_SqMat(nx, Q, SQR(prn.wsa)*tt, IISA(opt) + i);
			Init_SqMat(nx, Q, SQR(prn.wsg)*tt, IISG(opt) + i);
		}
	}
}

/* update p prediction */
static void udins_pp(double *P, const int nx, const double tt, const double *F, double *G, double *Q)
{
	double *I, *Ft, *Phi;
	double *T1, *T2;

	I = eye(nx); Ft = zeros(nx, nx); Phi = zeros(nx, nx);
	Mat_Num(F, tt, nx, nx, Ft);
	Mat_Add(nx, nx, 1.0, I, 1.0, Ft, 0.0, Phi);
	
	T1 = zeros(nx, nx); T2 = zeros(nx, nx);
	matmul("NN", nx, nx, nx, 1.0, Phi, P, 0.0, T1);		/* Phi*P */
	matmul("NT", nx, nx, nx, 1.0, T1, Phi, 0.0, P);		/* Phi*P*Phi' */
	matmul("NN", nx, nx, nx, 1.0, G, Q, 0.0, T2);		/* G*Q */
	matmul("NT", nx, nx, nx, 1.0, T2, G, 1.0, P);		/* Phi*P*Phi'+G*Q*G' */
	Mat_Sym(nx, P);

	free(I); free(Ft); free(Phi);
	free(T1); free(T2);
}

extern void UdState_INS(cate_t *cate, double *P, double *Cb2e, double *acce, double *gyro, const int nx)
{
	int i, j;
	double *F, *G, *Q;

	F = zeros(nx, nx); G = zeros(nx, nx); Q = zeros(nx, nx);
	udins_px(cate, &cate->opt, &cate->solv);
	udins_fg(&cate->opt, Cb2e, acce, gyro, nx, F, G);
	udins_q(&cate->opt, nx, cate->tt, Q);
	
	for (i = 0; i < nx; i++) {
		for (j = 0; j < nx; j++) {
			P[j + i * nx] = cate->P[j + i * cate->nx];
		}
	}
	udins_pp(P, nx, cate->tt, F, G, Q);
	
	free(F); free(G); free(Q);
}

static void udgnss_clk(cate_t *cate, rtk_t *rtk)
{
	int i, j, nx=GNC+GNR;
	double dtr, *F, *P, *T;
	double *I, *Ft, *Phi;

	for (i = 0; i < GNC; i++) {
		if (cate->popt.sateph == EPHOPT_PREC) {
			dtr = rtk->sol.dtr[0];
		}
		else {
			dtr = i == 0 ? rtk->sol.dtr[0] : rtk->sol.dtr[0] + rtk->sol.dtr[i];
		}
		Init_Px(cate, dtr*CLIGHT, CVAR_CLK, GIC(i, &cate->opt));
	}
	Init_Px(cate, rtk->sol.drif, CVAR_DRIF, GIR(&cate->opt));

	F = zeros(nx, nx); P = zeros(nx, nx); 
	I = eye(nx); Ft = zeros(nx, nx); Phi = zeros(nx, nx);
	T = zeros(nx, nx);
	for (i = 0; i < GNC; i++) {
		F[i + GNC * nx] = 1.0;
	}
	for (i = 0; i < nx; i++) {
		for (j = 0; j < nx; j++) {
			P[j + i * nx] = cate->P[GIC(j, &cate->opt) + GIC(i, &cate->opt) * cate->nx];
		}
	}
	
	Mat_Num(F, cate->tt, nx, nx, Ft);
	Mat_Add(nx, nx, 1.0, I, 1.0, Ft, 0.0, Phi);
	matmul("NN", nx, nx, nx, 1.0, Phi, P, 0.0, T);
	matmul("NT", nx, nx, nx, 1.0, T, Phi, 0.0, P);
	for (i = 0; i < nx; i++) {
		for (j = 0; j < nx; j++) {
			cate->P[GIC(j, &cate->opt) + GIC(i, &cate->opt) * cate->nx] = P[j + i * nx];
		}
	}
	free(F); free(P); free(T); free(I); free(Ft); free(Phi);
}

static void udgnss_trop(cate_t *cate)
{
	double pos[3], azel[] = { 0.0,PI / 2.0 }, ztd, var;
	int i = GIT(&cate->opt), j;

	if (cate->xx[i] == 0.0) {
		ecef2pos(cate->solv.re, pos);
		ztd = sbstropcorr(cate->solv.time, pos, azel, &var);
		Init_Px(cate, ztd, var, i);

		if (cate->popt.tropopt >= TROPOPT_ESTG) {
			for (j = i + 1; j < i + 3; j++) Init_Px(cate, 1E-6, CVAR_GRA, j);
		}
	}
	else {
		cate->P[i + i * cate->nx] += SQR(cate->popt.prn[2])*fabs(cate->tt);

		if (cate->popt.tropopt >= TROPOPT_ESTG) {
			for (j = i + 1; j < i + 3; j++) {
				cate->P[j + j * cate->nx] += SQR(cate->popt.prn[2] * 0.1)*fabs(cate->tt);
			}
		}
	}
}

static void udgnss_iono(cate_t *cate, rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
	const double *lam;
	double ion, sinel, pos[3], *azel;
	char *p;
	int i, j, k, gap_resion = GAP_RESION;

	if ((p = strstr(rtk->opt.pppopt, "-GAP_RESION="))) {
		sscanf(p, "-GAP_RESION=%d", &gap_resion);
	}
	for (i = 0; i < SATNUM; i++) {
		j = GII(i + 1, &rtk->opt, &cate->opt);
		if (cate->xx[j] != 0.0 && (int)rtk->ssat[i].outc[0] > gap_resion) {
			cate->xx[j] = 0.0;
		}
	}
	for (i = 0; i < n; i++) {
		j = GII(obs[i].sat, &rtk->opt, &cate->opt);
		if (cate->xx[j] == 0.0) {
			k = satsys(obs[i].sat, NULL) == SYS_GAL ? 2 : 1;
			lam = nav->lam[obs[i].sat - 1];
			if (obs[i].P[0] == 0.0 || obs[i].P[k] == 0.0 || lam[0] == 0.0 || lam[k] == 0.0) {
				continue;
			}
			ion = (obs[i].P[0] - obs[i].P[k]) / (1.0 - SQR(lam[k] / lam[0]));
			ecef2pos(rtk->sol.rr, pos);
			azel = rtk->ssat[obs[i].sat - 1].azel;
			ion /= ionmapf(pos, azel);
			Init_Px(cate, ion, CVAR_IONO, j);
		}
		else {
			sinel = sin(MAX(rtk->ssat[obs[i].sat - 1].azel[1], 5.0*D2R));
			cate->P[j + j * cate->nx] += SQR(rtk->opt.prn[1] / sinel)*fabs(cate->tt);
		}
	}
}

static void udgnss_amb(cate_t *cate, rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
	const double *lam;
	double L[NFREQ], P[NFREQ], Lc, Pc, bias[MAXOBS], offset = 0.0, pos[3] = { 0 };
	double ion, dantr[NFREQ] = { 0 }, dants[NFREQ] = { 0 };
	int i, j, k, l, f, sat, slip[MAXOBS] = { 0 }, clk_jump = 0;

	/* handle day-boundary clock jump */
	if (rtk->opt.posopt[5]) {
		clk_jump = ROUND(time2gpst(obs[0].time, NULL) * 10) % 864000 == 0;
	}
	CycSlip_Detec(cate, rtk, obs, n, nav);
	ecef2pos(cate->solv.re, pos);

	for (f = 0; f < GNF(&rtk->opt); f++) {

		/* reset phase-bias if expire obs outage counter */
		for (i = 0; i < SATNUM; i++) {
			if (++rtk->ssat[i].outc[f] > (unsigned int)rtk->opt.maxout || clk_jump) {
				Init_Px(cate, 0.0, 0.0, GIB(i+1, f, &rtk->opt, &cate->opt));
			}
		}
		for (i = k = 0; i < n&&i < MAXOBS; i++) {
			sat = obs[i].sat;
			j = GIB(sat, f, &rtk->opt, &cate->opt);
			Cate_corr_meas(obs + i, nav, rtk->ssat[sat - 1].azel, &rtk->opt, dantr, dants,
				0.0, L, P, &Lc, &Pc);

			bias[i] = 0.0;

			if (rtk->opt.ionoopt == IONOOPT_IFLC) {
				bias[i] = Lc - Pc;
				slip[i] = rtk->ssat[sat - 1].slip[0] || rtk->ssat[sat - 1].slip[1];
			}
			else if (L[f] != 0.0&&P[f] != 0.0) {
				slip[i] = rtk->ssat[sat - 1].slip[f];
				l = satsys(sat, NULL) == SYS_GAL ? 2 : 1;
				lam = nav->lam[sat - 1];
				if (obs[i].P[0] == 0.0 || obs[i].P[l] == 0.0 ||
					lam[0] == 0.0 || lam[l] == 0.0 || lam[f] == 0.0) continue;
				ion = (obs[i].P[0] - obs[i].P[l]) / (1.0 - SQR(lam[l] / lam[0]));
				bias[i] = L[f] - P[f] + 2.0*ion*SQR(lam[f] / lam[0]);
			}
			if (cate->xx[j] == 0.0 || slip[i] || bias[i] == 0.0) continue;

			offset += bias[i] - cate->xx[j];
			k++;
		}
		/* correct phase-code jump to ensure phase-code coherency */
		if (k >= 2 && fabs(offset / k) > 0.0005*CLIGHT) {
			for (i = 0; i < SATNUM; i++) {
				j = GIB(sat, f, &rtk->opt, &cate->opt);
				if (cate->xx[j] != 0.0){ cate->xx[j] += offset / k; }
			}
		}
		for (i = 0; i < n&&i < MAXOBS; i++) {
			sat = obs[i].sat;
			j = GIB(sat, f, &rtk->opt, &cate->opt);

			cate->P[j + j * cate->nx] += SQR(rtk->opt.prn[0])*fabs(rtk->tt);
			/* 1、相位偏差等于0，不重新初始化，没有相位偏差怎么初始化呢？ */
			/* 2、相位偏差不等于0时，
				(1) x!=0(即上一历元估计过),slip==0,没周跳，不初始化；
				(2) x!=0(即上一历元估计过),slip==1,有周跳，初始化(有周跳相当于原来的偏差不准，所以重新初始化)；
				(3) x==0(新加入的估计卫星),slip==0,没周跳，初始化(新加入的估计卫星能不初始化吗)；
				(3) x==0(新加入的估计卫星),slip==1,有周跳，初始化 */
			/* 综上，只要bias存在，有周跳或新加入卫星，一定初始化 */
			if (bias[i] == 0.0 || (cate->xx[j] != 0.0 && !slip[i])) continue;

			Init_Px(cate, bias[i], CVAR_BIAS, j);
		}
	}
}

extern void UdState_TTCP(cate_t *cate, rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav,
	double *Cb2e, double *acce, double *gyro)
{
	int nx = CI_N(&cate->opt);
	double *P = zeros(nx, nx);
	
	udgnss_clk(cate, rtk);

	if (rtk->opt.tropopt == TROPOPT_EST || rtk->opt.tropopt == TROPOPT_ESTG) {
		udgnss_trop(cate);
	}
	if (rtk->opt.ionoopt == IONOOPT_EST) {
		udgnss_iono(cate, rtk, obs, n, nav);
	}
	udgnss_amb(cate, rtk, obs, n, nav);
	
	UdState_INS(cate, P, Cb2e, acce, gyro, nx);
	Mat_Sub(cate->P, cate->nx, P, nx, 0, 0, 1.0);

	free(P);
}

extern void Udsol(rtk_t *rtk, const obsd_t *obs, const int n, cate_t *cate, double *Cb2e)
{
	int i, j;
	double  datt[3];

	for (i = 0; i < 3; i++) {
		datt[i] = cate->xx[IIA + i] - cate->solv.att[i];
	}
	Out_Matrix(cate->P, cate->nx, cate->nx);
	/* correction */
	for (i = 0; i < 3; i++) {
		cate->solv.re[i] = cate->xx[IIP + i];
		cate->solv.qr[i] = (float)cate->P[IIP + i + (IIP + i) * cate->nx];
		cate->solv.ve[i] = cate->xx[IIV + i];
		cate->solv.qv[i] = (float)cate->P[IIV + i + (IIV + i) * cate->nx];
		cate->solv.att[i] = cate->xx[IIA + i];
		cate->solv.qa[i] = (float)cate->P[IIA + i + (IIA + i) * cate->nx];

		cate->solv.ba[i] = cate->xx[IIBA + i];
		cate->solv.bg[i] = cate->xx[IIBG + i];
		if (cate->opt.est_scale) {
			cate->solv.sa[i] = cate->xx[IISA(&cate->opt) + i];
			cate->solv.sg[i] = cate->xx[IISG(&cate->opt) + i];
		}
	}
	cate->solv.qr[3] = (float)cate->P[1];
	cate->solv.qr[4] = (float)cate->P[2 + cate->nx];
	cate->solv.qr[5] = (float)cate->P[2];

	cate->solv.qv[3] = (float)cate->P[IIV + 1 + IIV * cate->nx];
	cate->solv.qv[4] = (float)cate->P[IIV + 2 + (IIV + 1)*cate->nx];
	cate->solv.qv[5] = (float)cate->P[IIV + 2 + IIV * cate->nx];

	cate->solv.qa[3] = (float)cate->P[IIA + 1 + IIA * cate->nx];
	cate->solv.qa[4] = (float)cate->P[IIA + 2 + (IIA + 1)*cate->nx];
	cate->solv.qa[5] = (float)cate->P[IIA + 2 + IIA * cate->nx];

	if (cate->opt.corr_dcm) {
		double pos[3];
		ecef2pos(cate->solv.re, pos);
		Att_to_Cb2e(cate->solv.att, pos, Cb2e);
	}

#if 0
	if (opt->corr_dcm == CORR_DCM) {
		double *I, dAtt[9], temp1[9], temp2[9];
		I = eye(3);
		Mat_Skew(datt, dAtt);
		Mat_Add(3, 3, 1.0, I, -1.0, dAtt, 0.0, temp1);
		matmul("NN", 3, 3, 3, 1.0, temp1, Cb2e, 0.0, temp2);
		Mat_Norm(temp2, 3, Cb2e);
		free(I);
	}
#endif

	if (cate->solv.stat < CSOL_TTCP) {
		return;
	}

	for (i = 0; i < n&&i < MAXOBS; i++) {
		if (!cate->csat[obs[i].sat - 1].vs) continue;
		cate->csat[obs[i].sat - 1].lock++;
		cate->csat[obs[i].sat - 1].outc = 0;
	}

	cate->solv.ns = 0;
	for (i = 0; i < n&&i < MAXOBS; i++) {
		for (j = 0; j < rtk->opt.nf; j++) {
			if (!rtk->ssat[obs[i].sat - 1].vsat[j]) continue;
			rtk->ssat[obs[i].sat - 1].lock[j]++;
			rtk->ssat[obs[i].sat - 1].outc[j] = 0;
			if (j == 0) cate->solv.ns++;
		}
	}

	for (i = 0; i < NSYS; i++) { cate->solv.off[i] = cate->xx[GIC(i, &cate->opt)]; }

	cate->solv.drif = cate->xx[GIR(&cate->opt)];

	if (rtk->opt.tropopt == TROPOPT_EST || rtk->opt.tropopt == TROPOPT_ESTG) {
		cate->solv.trop = cate->xx[GIT(&cate->opt)];
	}
	if (rtk->opt.ionoopt == IONOOPT_EST) {
		for (i = 0; i < SATNUM; i++) {
			cate->solv.iono[i] = cate->xx[GII(i + 1, &cate->popt, &cate->opt)];
		}
	}
	for (i = 0; i < SATNUM; i++) {
		for (j = 0; j < NFREQ; j++) {
			cate->solv.amb[j][i] = cate->xx[GIB(i + 1, j, &cate->popt, &cate->opt)];
		}
	}
}