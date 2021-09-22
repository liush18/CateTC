/*------------------------------------------------------------------------------
 * Cate_CoupNav.c : Coupled Navigation in ECEF
 *
 *			Copyright (C) 2019-2019 by Liu Shuai, All rights reserved.
 *
 * version :
 * history : 2019/04/22  new
 *------------------------------------------------------------------------------*/

#include "Cate.h"
#include "State.h"

#define ROUND(x)		(int)floor((x)+0.5)
#define GAP_RESION		120					/* default gap to reset ionos parameters (ep) */
#define MAX_ITER		8					/* max number of iterations */
#define THRES_REJECT	4.0					/* reject threshold of posfit-res (sigma) */

static void calcu_doppler(const obsd_t *obs, const nav_t *nav, const double *azel, 
							const prcopt_t *opt, double *D, double *Dc)
{
	const double *lam = nav->lam[obs->sat - 1];
	double C1, C2;
	int i, sys;

	for (i = 0; i < NFREQ; i++) {
		D[i] = 0.0;
		if (lam[i] == 0.0 || obs->D[i] == 0.0) continue;
		if (testsnr(0, 0, azel[1], obs->SNR[i] * 0.25, &opt->snrmask)) continue;

		D[i] = obs->D[i] * lam[i];
	}
	if (D[0] != 0.0) {
		*Dc = D[0];
	}
	else {
		*Dc = 0.0;
	}

	*Dc = 0.0;
	sys = satsys(obs->sat, NULL);
	i = (sys&(SYS_GAL | SYS_SBS)) ? 2 : 1; /* L1/L2 or L1/L5 */
	if (lam[0] == 0.0 || lam[i] == 0.0) return;

	C1 = SQR(lam[i]) / (SQR(lam[i]) - SQR(lam[0]));
	C2 = -SQR(lam[0]) / (SQR(lam[i]) - SQR(lam[0]));

#if 0
	if (D[0] != 0.0&&D[i] != 0.0) {
		*Dc = C1 * D[0] + C2 * D[i];
	}
	else if (D[0] != 0.0) {
		*Dc = D[0];
	}
	else if (D[i] != 0.0) {
		*Dc = D[i];
	}
	else {
		*Dc = 0.0;
	}
#endif
}

/* phase and code residuals --------------------------------------------------
* post
* obs
* n			当前历元的观测个数（观测卫星数）
* rs		(6,n)，每一列为观测卫星的坐标及速度
* dts		(2,n)，每一列为钟差和钟漂(卫星)
* var_rs	(1,n)，每颗卫星位置和钟差误差的平方，用于卫星可用性检测及设置R阵
* svh
*/
static int ttcp_res(int post, const obsd_t *obs, int n, const double *rs,
					const double *dts, const double *var_rs, const int *svh,
					const double *dr, int *exc, const nav_t *nav, 
					const double *x, rtk_t *rtk, double *v, double *H, double *R,
					double *azel, cate_t *cate)
{
	const double *lam;
	prcopt_t *popt = &rtk->opt; const cateopt_t *opt = &cate->opt;
	double y, r, cdtr, bias, C, rr[3], pos[3], e[3], dtdx[3], L[NFREQ], P[NFREQ], D[NFREQ], Lc, Pc, Dc;
	double var[MAXOBS * 2], dtrp = 0.0, dion = 0.0, vart = 0.0, vari = 0.0, dcb;
	double dantr[NFREQ] = { 0 }, dants[NFREQ] = { 0 };
	double ve[MAXOBS * 2 * NFREQ] = { 0 }, vmax = 0.0;
	char str[32];
	int ne = 0, obsi[MAXOBS * 2 * NFREQ] = { 0 }, frqi[MAXOBS * 2 * NFREQ], maxobs, maxfrq, rej;
	int i, j, k, sat, sys, nv = 0, nx = cate->nx, stat = 1;
	int sv = 3;
	double vs[3], rate;

	time2str(obs[0].time, str, 2);

	for (i = 0; i < MAXSAT; i++) {
		for (j = 0; j < popt->nf; j++) { rtk->ssat[i].vsat[j] = 0; }
		if (!post) {
			for (j = 0; j < NFREQ; j++) {
				cate->csat[i].resc_prio[j] = 0.0;
				cate->csat[i].resp_prio[j] = 0.0;
				cate->csat[i].resd_prio[j] = 0.0;
			}
		}
		else {
			for (j = 0; j < NFREQ; j++) {
				cate->csat[i].resc_post[j] = 0.0;
				cate->csat[i].resp_post[j] = 0.0;
				cate->csat[i].resd_post[j] = 0.0;
			}
		}
	}
	for (i = 0; i < 3; i++) rr[i] = x[i] + dr[i];
	ecef2pos(rr, pos);

	for (i = 0; i < n&&i < MAXOBS; i++) {
		sat = obs[i].sat;
		lam = nav->lam[sat - 1];
		if (lam[j / 2] == 0.0 || lam[0] == 0.0) continue;

		if ((r = geodist(rs + i * 6, rr, e)) <= 0.0 ||
			satazel(pos, e, azel + i * 2) < popt->elmin) {
			exc[i] = 1;
			continue;
		}
		if (!(sys = satsys(sat, NULL)) || !rtk->ssat[sat - 1].vs ||
			satexclude(obs[i].sat, var_rs[i], svh[i], popt) || exc[i]) {
			exc[i] = 1;
			continue;
		}
		/* tropospheric and ionospheric model */
		if (!Cate_model_trop(obs[i].time, pos, azel + i * 2, popt, x, dtdx, nav, &dtrp, &vart, cate) ||
			!Cate_model_iono(obs[i].time, pos, azel + i * 2, popt, sat, x, nav, &dion, &vari, cate)) {
			continue;
		}
		/* satellite and receiver antenna model */
		if (popt->posopt[0]) Cate_satantpcv(rs + i * 6, rr, nav->pcvs + sat - 1, dants);
		antmodel(popt->pcvr, popt->antdel[0], azel + i * 2, popt->posopt[1], dantr);

		/* phase windup model */
		if (!Cate_model_phw(rtk->sol.time, sat, nav->pcvs[sat - 1].type,
			popt->posopt[2] ? 2 : 0, rs + i * 6, rr, &rtk->ssat[sat - 1].phw)) {
			continue;
		}
		/* corrected phase and code measurements */
		Cate_corr_meas(obs + i, nav, azel + i * 2, &rtk->opt, dantr, dants,
			rtk->ssat[sat - 1].phw, L, P, &Lc, &Pc);

		calcu_doppler(obs + i, nav, azel + i * 2, popt, D, &Dc);

		/* stack phase and code residuals {L1,P1,D1,L2,P2,D1,...} */
		for (j = 0; j < sv*GNF(popt); j++) {

			dcb = bias = 0.0;

			if (popt->ionoopt == IONOOPT_IFLC) {
				y = (j%sv == 0) ? Lc : ((j%sv == 1) ? Pc : Dc);
				if (y == 0.0) { continue; }
			}
			else {
				y = (j%sv == 0) ? L[j/sv] : ((j%sv == 1) ? P[j/sv] : D[j/sv]);
				if (y == 0.0) { continue; }
			}
			C = SQR(lam[j / sv] / lam[0])*ionmapf(pos, azel + i * 2)*(j % sv == 0 ? -1.0 : 1.0);

			if (j%sv == 2) {
				for (k = 0; k < nx; k++) { H[k + nx*nv] = (k > 2 && k < 6) ? -e[k - 3] : 0.0; }
				/* satellite velocity relative to receiver in ecef */
				for (k = 0; k < 3; k++) vs[k] = rs[k + 3 + i * 6] - x[3 + k];
				/* range rate with earth rotation correction */
				rate = dot(vs, e, 3) + OMGE / CLIGHT * (rs[4 + i * 6] * rr[0] + rs[1 + i * 6] * x[3] -
					rs[3 + i * 6] * rr[1] - rs[i * 6] * x[4]);
				H[GIR(opt) + nx*nv] = 1.0;
			}
			else {
				for (k = 0; k < nx; k++) { H[k + nx * nv] = (k < 3) ? -e[k] : 0.0; }
			}

			/* receiver clock */
			if (j%sv != 2) {
				k = sys == SYS_GLO ? 1 : (sys == SYS_GAL ? 2 : (sys == SYS_CMP ? 3 : 0));
				cdtr = x[GIC(k, opt)];
				H[GIC(k, opt) + nx * nv] = 1.0;
			}
			
			if (j%sv!=2 &&(popt->tropopt == TROPOPT_EST || popt->tropopt == TROPOPT_ESTG)) {
				for (k = 0; k < (popt->tropopt >= TROPOPT_ESTG ? 3 : 1); k++) {
					H[GIT(opt) + k + nx * nv] = dtdx[k];
				}
			}

			if (j%sv != 2 && popt->ionoopt == IONOOPT_EST) {
				if (cate->xx[GII(sat, popt, opt)] == 0.0) continue;
				H[GII(sat, popt, opt) + nx * nv] = C;
			}

			if (j % sv == 0) {	/* phase bias */
				if ((bias = x[GIB(sat, j/sv, popt, opt)]) == 0.0) continue;
				H[GIB(sat, j/sv, popt, opt) + nx * nv] = 1.0;
			}
			
			/* residual */
			if (j%sv == 2) {
				v[nv] = -y - (rate + x[GIR(opt)] - CLIGHT * dts[1 + i * 2]);
				if (v[nv] > 0.1 || v[nv] < -0.1)
					continue;
			}
			else {
				v[nv] = y - (r + cdtr - CLIGHT * dts[i * 2] + dtrp + C * dion + dcb + bias);
			}
			
			if (!post) {
				if (j % sv == 0) { cate->csat[sat - 1].resc_prio[j/sv] = v[nv]; }
				else if (j%sv == 1) { cate->csat[sat - 1].resp_prio[j/sv] = v[nv]; }
				else { cate->csat[sat - 1].resd_prio[j/sv] = v[nv]; }
			}
			else {
				if (j % sv == 0) { cate->csat[sat - 1].resc_post[j / sv] = v[nv]; }
				else if (j%sv == 1){ cate->csat[sat - 1].resp_post[j / sv] = v[nv]; }
				else{ cate->csat[sat - 1].resd_post[j / sv] = v[nv]; }
			}

			if (j%sv == 2) {
				var[nv] = SQR(0.1);
			}
			else {
				var[nv] = Cate_varerr(obs[i].sat, sys, azel[1 + i * 2], j / sv, j % sv, popt) +
					vart + SQR(C)*vari + var_rs[i];
			}

			/* reject satellite by pre-fit residuals */
			if (!post && popt->maxinno > 0.0 && fabs(v[nv]) > popt->maxinno) {
				follow(2, "outlier (%d) rejected %s sat=%2d %s%d res=%9.4f el=%4.1f\n",
					post, str, sat, j % sv ? (j%sv==1? "P":"D") : "L", j / sv + 1, v[nv], azel[1 + i * 2] * R2D);
				exc[i] = 1; rtk->ssat[sat - 1].rejc[j / sv]++;
				continue;
			}
			/* record large post-fit residuals */
			if (post&&fabs(v[nv]) > sqrt(var[nv])*THRES_REJECT) {
				obsi[ne] = i; frqi[ne] = j; ve[ne] = v[nv]; ne++;
			}
			if (j % sv == 0) rtk->ssat[sat - 1].vsat[j / sv] = 1;
			nv++;
		}
	}
	/* reject satellite with large and max post-fit residual */
	if (post&&ne > 0) {
		vmax = ve[0]; maxobs = obsi[0]; maxfrq = frqi[0]; rej = 0;
		for (j = 1; j < ne; j++) {
			if (fabs(vmax) >= fabs(ve[j])) continue;
			vmax = ve[j]; maxobs = obsi[j]; maxfrq = frqi[j]; rej = j;
		}
		sat = obs[maxobs].sat;
		follow(2, "outlier (%d) rejected %s sat=%2d %s%d res=%9.4f el=%4.1f\n",
			post,str,sat,maxfrq%sv?(maxfrq%sv==1?"P":"D"):"L",maxfrq/sv+1,vmax,azel[1+maxobs*2]*R2D);
		exc[maxobs] = 1; rtk->ssat[sat - 1].rejc[maxfrq/sv]++; stat = 0;
		ve[rej] = 0;
	}

	for (i = 0; i < nv; i++) for (j = 0; j < nv; j++) {
		R[i + j * nv] = i == j ? var[i] : 0.0;
	}

	return post ? stat : nv;
}

extern void TC_Pos(cate_t *cate, double *Cb2e, double *acce, double *gyro,
	rtk_t *rtk, obsd_t *obs, int n, const nav_t *nav)
{
	const cateopt_t *opt = &cate->opt;
	const prcopt_t *popt = &rtk->opt;
	double *rs, *dts, *var, *v, *H, *R, *azel, *xp, *Pp, dr[3] = { 0 };
	int i, nv, svh[MAXOBS], exc[MAXOBS] = { 0 };

	rs = mat(6, n); dts = mat(2, n); var = mat(1, n); azel = zeros(2, n);

	/* temporal update of ekf states */
	UdState_TTCP(cate, rtk, obs, n, nav, Cb2e, acce, gyro);

	/* satellite positions and clocks */
	satposs(obs[0].time, obs, n, nav, rtk->opt.sateph, rs, dts, var, svh);

	/* exclude measurements of eclipsing satellite (block IIA) */
	if (popt->posopt[3]) {
		Cate_testeclipse(obs, n, nav, rs);
	}
	/* earth tides correction */
	if (popt->tidecorr) {
		tidedisp(gpst2utc(obs[0].time), cate->solv.re, popt->tidecorr == 1 ? 1 : 7, &nav->erp,
			popt->odisp[0], dr);
	}
	
	nv = SATNUM * 3 * popt->nf;
	xp = mat(cate->nx, 1); Pp = mat(cate->nx, cate->nx); 
	v = mat(nv, 1); H = mat(cate->nx, nv); R = mat(nv, nv);
	
	for (i = 0; i < MAX_ITER; i++) {
		if (i) printf("\n%d\n", i);

		matcpy(xp, cate->xx, cate->nx, 1);
		matcpy(Pp, cate->P, cate->nx, cate->nx);

		if (!(nv = ttcp_res(0, obs, n, rs, dts, var, svh, dr, exc, nav, xp, rtk, v, H, R, azel, cate))) {
			break;
		}
		if (!Kalman(xp, Pp, v, H, R, cate->nx, nv)) {
			break;
		}
		if (ttcp_res(i+1, obs, n, rs, dts, var, svh, dr, exc, nav, xp, rtk, v, H, R, azel, cate)) {
			matcpy(cate->xx, xp, cate->nx, 1);
			matcpy(cate->P, Pp, cate->nx, cate->nx);
			cate->solv.stat = obs->D[0] == 0 ? CSOL_TTCP : CSOL_TTCD;
			break;
		}
	}

	Udsol(rtk, obs, n, cate, Cb2e);

	free(xp); free(Pp);
	free(H); free(v); free(R);
	free(rs); free(dts); free(var); free(azel);
}