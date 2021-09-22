/*------------------------------------------------------------------------------
* Cate_Solutions.c : Cate Solutions Output
*
*          Copyright (C) 2019-2019 by Liu Shuai, All rights reserved.
*
* version :
* history : 2019/04/22  new
*------------------------------------------------------------------------------*/

#include "Cate.h"
#include "State.h"

static char stime[]	= "$TIME,week,tow,yyyy/mm/dd hh:mm:ss.ssss,stat,valid sat";
static char spos []	= "$POS,posx,posy,posz";
static char svel []	= "$VEL,velx,vely,velz";
static char satt []	= "$ATT,pitch,roll,yaw(deg)";
static char sbias[]	= "$BIAS,accex,accey,accez(mg),gyrox,gyroy,gyroz(deg/hr)";
static char sscal[]	= "$SCALE,accex,accey,accez,gyrox,gyroy,gyroz";
static char sclk [] = "$CLK,bias1,bias2,bias3,bias4(m),drif(m/s)";
static char strop[] = "$TROP,ztd";
static char siono[] = "$IONO,sat,iono";
static char ssat [] = "$SAT,sat,frq,az,el,resp_prio,resc_prio,resd_prio,resp_post,resc_post,resd_post,vsat,snr,slip,lock,outc,slipc,rejc,amb";

static FILE *fp_sol = NULL;
static FILE *fp_stat = NULL;
static FILE *fp_flw = NULL;		/* file pointer of follow */

extern int showproc(char *format, ...)
{
	va_list arg;
	va_start(arg, format); vfprintf(stderr, format, arg); va_end(arg);
	fprintf(stderr, "\n");
	return 0;
}

static void followopen(const char *path)
{
	if (!(fp_flw = fopen(path, "w"))) {
		fp_flw = stderr;
	}
}

static void followclose(void)
{
	if (fp_flw && fp_flw != stderr) {
		fclose(fp_flw);
	}
	fp_flw = NULL;
}

extern void follow(int level, const char *format, ...)
{
	va_list ap;
	if (level <= 1) {
		va_start(ap, format); vfprintf(stderr, format, ap); va_end(ap);
	}
	if (!fp_flw) return;
	fprintf(fp_flw, "%d ", level);
	va_start(ap, format); vfprintf(fp_flw, format, ap); va_end(ap);
	fflush(fp_flw);
}

/* open outfile path (0:error, 1:ok) */
extern int File_Open(const cateopt_t *opt, const char *path)
{
	char buff[1024];

	if (*path && !(fp_sol = fopen(path, "a"))) {
		follow(1, "ERROR : File_Open-open file error, %s\n", path);
		return 0;
	}
	if (opt->outsolstat) {
		strcpy(buff, path);
		strcat(buff, ".stat");
		if (!(fp_stat = fopen(buff, "w"))) {
			follow(1, "WARNING : File_Open-open file error, %s\n", buff);
		}
		else {
			fprintf(fp_stat, "%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n\n", 
				stime, spos, svel, satt, sbias, sscal, sclk, strop, siono, ssat);
		}
	}

	if (opt->outdebug) {
		strcpy(buff, path);
		strcat(buff, ".flw");
		followopen(buff);
	}

	follow(2, "Open files...\n");

	return 1;
}

extern void File_Close(const cateopt_t *opt)
{
	follow(2, "Close files...\n");

	if (fp_sol) {
		fclose(fp_sol);
		fp_sol = NULL;
	}
	if (fp_stat) {
		fclose(fp_stat);
		fp_stat = NULL;
	}
	followclose();
}

extern int Out_Header(const char *outfile)
{
	FILE *fp;
	char *sep = " ";

	follow(2, "Out_Header...\n");

	if (!(fp = fopen(outfile, "w")))
		return 0;
	fprintf(fp, "%s\n", " % position(m in ecef), velocity(m/s in ecef), attitude(deg)");
	fprintf(fp, "%15s%15s%15s%15s%11s%11s%11s%13s%13s%13s%4s%4s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s%9s\n", 
		"%     GPST    ", "x","y", "z", "vx", "vy", "vz", "pitch", "roll", "yaw", "Q", "ns",
		"sdx", "sdy", "sdz", "sdzy", "sdyz", "sdzx", 
		"sdvx", "sdvy", "sdvz", "sdvzy", "sdvyz", "sdvzx", 
		"sdp", "sdr", "sdy", "sdpr", "sdry", "sdyp");
	fclose(fp);

	return 1;
}

static int out_solstat(cate_t *cate, rtk_t *rtk, char *buff)
{
	ssat_t *ssat;
	char str[1024], *p = buff, id[32];
	double tow;
	int i, week;

	if (cate->solv.stat < CSOL_LSCS) return 0;

	time2str(cate->solv.time, str, 3);
	tow = time2gpst(cate->solv.time, &week);

	p += sprintf(p, "$TIME,%d,%.3f,%s,%d,%2d\n", week, tow, str, cate->solv.stat, cate->solv.ns);
	p += sprintf(p, "$POS,%.4f,%.4f,%.4f\n", cate->solv.re[0], cate->solv.re[1], cate->solv.re[2]);
	p += sprintf(p, "$VEL,%.4f,%.4f,%.4f\n", cate->solv.ve[0], cate->solv.ve[1], cate->solv.ve[2]);
	p += sprintf(p, "$ATT,%.6f,%.6f,%.6f\n", cate->solv.att[0]*R2D, cate->solv.att[1]*R2D, cate->solv.att[2]*R2D);
	p += sprintf(p, "$BIAS,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
		cate->solv.ba[0] / MG2AC, cate->solv.ba[1] / MG2AC, cate->solv.ba[2] / MG2AC,
		cate->solv.bg[0] / DH2RS, cate->solv.bg[1] / DH2RS, cate->solv.bg[2] / DH2RS);
	if (cate->opt.est_scale) {
		p += sprintf(p, "$SCALE,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
			cate->solv.sa[0], cate->solv.sa[1], cate->solv.sa[2],
			cate->solv.sg[0], cate->solv.sg[1], cate->solv.sg[2]);
	}

	if (cate->solv.stat < CSOL_TTCP) return (int)(p - buff);

	p += sprintf(p, "$CLK,%.4f,%.4f,%.4f,%.4f,%.4f\n", 
		cate->solv.off[0], cate->solv.off[1], cate->solv.off[2], cate->solv.off[3], cate->solv.drif);

	if (rtk->opt.tropopt == TROPOPT_EST || rtk->opt.tropopt == TROPOPT_ESTG) {
		p += sprintf(p, "$TROP,%.4f\n", cate->solv.trop);
	}
	if (rtk->opt.ionoopt == IONOOPT_EST) {
		for (i = 0; i < SATNUM; i++) {
			ssat = rtk->ssat + i;
			if (!ssat->vs) continue;
			satno2id(i + 1, id);
			p += sprintf(p, "$IONO,%s,%.4f\n", id, cate->solv.iono[i]);
		}
	}

	return (int)(p - buff);
}

extern void cOutSolStat(cate_t *cate, rtk_t *rtk)
{
	ssat_t *ssat;
	csat_t *csat;
	char buff[MAXSOLMSG + 1], id[32];
	int i, j, n, nfreq, nf=GNF(&cate->popt);

	if (!fp_stat || !cate->opt.outsolstat) return;
	if (!(n = out_solstat(cate, rtk, buff))) {
		return;
	}
	buff[n] = '\0';
	fprintf(fp_stat, "%s", buff);
	if (cate->solv.stat < CSOL_TTCP) {
		fprintf(fp_stat, "\n");
		return;
	}

	nfreq = cate->opt.coupled_type >= PMODE_TTCP ? nf : 1;

	for (i = 0; i < SATNUM; i++) {
		ssat = rtk->ssat + i;
		csat = cate->csat + i;
		if (!ssat->vs) continue;
		satno2id(i + 1, id);
		for (j = 0; j < nfreq; j++) {
			fprintf(fp_stat, "$SAT,%s,f=%d,az=%.2f,el=%.2f,ppr=%.4f,cpr=%.4f,dpr=%.4f,ppt=%.4f,cpt=%.4f,dpt=%.4f,v=%d,snr=%.1f,slip=%d,lock=%d,outc=%d,slipc=%d,rejc=%d,amb=%.4f\n",
				id, j + 1, ssat->azel[0] * R2D, ssat->azel[1] * R2D,
				csat->resp_prio[j], csat->resc_prio[j], csat->resd_prio[j],
				csat->resp_post[j], csat->resc_post[j], csat->resd_post[j],
				ssat->vsat[j], ssat->snr[j] * 0.25, ssat->slip[j], ssat->lock[j],
				ssat->outc[j], ssat->slipc[j], ssat->rejc[j], cate->solv.amb[j][i]);
		}
	}
	fprintf(fp_stat, "\n");
}

extern void cOutSols(solvd_t *solv)
{
	int i, week;
	double tow;
	char buff[MAXSOLMSG + 1], *sep = " ";
	char *p = buff;

	tow = time2gpst(solv->time, &week);

	p += sprintf(p, "%4d%s%10.3f%s", week, sep, tow, sep);
	for (i = 0; i < 3; i++) {
		p += sprintf(p, "%14.4f%s", solv->re[i], sep);
	}
	for (i = 0; i < 3; i++) {
		p += sprintf(p, "%10.5f%s", solv->ve[i], sep);
	}
	for (i = 0; i < 3; i++) {
		p += sprintf(p, "%12.6f%s", solv->att[i] * R2D, sep);
	}
	p += sprintf(p, "%3d%s%3d%s", solv->stat, sep, solv->ns, sep);
	/* position */
	for (i = 0; i < 3; i++) {
		p += sprintf(p, "%8.5f%s", SQRT(solv->qr[i]), sep);
	}
	for (i = 0; i < 3; i++) {
		p += sprintf(p, "%8.5f%s", SQRC(solv->qr[i+3]), sep);
	}
	/* velocity */
	for (i = 0; i < 3; i++) {
		p += sprintf(p, "%8.5f%s", SQRT(solv->qv[i]), sep);
	}
	for (i = 0; i < 3; i++) {
		p += sprintf(p, "%8.5f%s", SQRC(solv->qv[i + 3]), sep);
	}
	/* attitude */
	for (i = 0; i < 3; i++) {
		p += sprintf(p, "%8.5f%s", SQRT(solv->qa[i]), sep);
	}
	for (i = 0; i < 3; i++) {
		p += sprintf(p, "%8.5f%s", SQRC(solv->qa[i + 3]), sep);
	}

	buff[(int)(p - buff)] = '\0';

	fprintf(fp_sol, "%s\n", buff);
}