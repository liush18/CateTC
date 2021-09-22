/*------------------------------------------------------------------------------
* Cate.h : Header
*
*          Copyright (C) 2019-2019 by Liu Shuai, All rights reserved.
*
* version :
* history : 2019/04/17  new
*------------------------------------------------------------------------------*/
#include "rtklib.h"

#ifndef CATE_H
#define CATE_H

#define DH2RS		(D2R/3600.0)		/* deg/hr to rad/s */
#define DSH2RSS		(D2R/60.0)			/* deg/sqrt(h) to rad/sqrt(s) */
#define MG2AC		9.80665E-3			/* mg to m/s^2 */
#define UG2AC		9.80665E-6			/* ug to m/s^2 */
#define MGA2AC		1E-5				/* mgal to m/s^2 */
#define UGA2AC		1E-8				/* ugal to m/s^2 */

#define E2_WGS84	(2*FE_WGS84-FE_WGS84*FE_WGS84)	/* square of eccentricity of the earth */
#define FM			3986005E8						/* geocentric gravitational constant */
#define J2			0.00108263						/* coefficient of series expansion */
#define J4			-2.37091222e-6	
#define J6			6.08347e-9

#define SQR(x)		((x)*(x))
#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define SQRC(x)		((x)<0.0?-sqrt(-x):sqrt(x))		/* sqrt of covariance */
#define MAX(x,y)    ((x)>(y)?(x):(y))
#define MIN(x,y)    ((x)<(y)?(x):(y))

/* alignment or get dcm from attitude */
#define ALIGNATT		0		/* attitude */
#define ALIGNCOR		1		/* coarse */

#define LEVERCORR		1		/* lever correction */

/* satellite numbers */
#define SATNUM			32			/* satllite numbers */

/* coupled types */
#define PMODE_GNSS		0
#define PMODE_PINS		1		/* pure ins */
#define PMODE_LSCP		2		/* loosely coupled */
#define PMODE_TTCP		3		/* tightly coupled */

/* solution type */
#define CSOL_NONE	0		/* coupled solution status: no solution */
#define CSOL_PINS	1		/* coupled solution status: pure ins */
#define CSOL_LSCS	2		/* coupled solution status: loosely coupled of spp and ins */
#define CSOL_LSCP	3		/* coupled solution status: loosely coupled of ppp and ins */
#define CSOL_TTCP	4		/* coupled solution status: tightly coupled of ppp and ins */
#define CSOL_TTCD	5		/* coupled solution status: tightly coupled of ppp and ins with doppler */

#define	MAXEXFILE2	10
#define DOPPLER		1
/* type definitions---------------------------------------*/

typedef struct {
	int week;
	double second;
}gpst_t;

typedef struct {	/* coupled prn */
	int flag;
	double wa;
	double wg;
	double wba;
	double wbg;
	double wsa;
	double wsg;
}cprn_t;

typedef struct {	/* input processing options */
/* gnss */
	gtime_t ts;
	gtime_t te;
	double ti;		/* gnss processing interval */
	double tu;		/* gnss processing unit time */
	double si;		/* gnss sample interval */

/* ins */
	double stini;	/* start time of initial alignment */
	double edini;	/* end time of initial alignment */
	double stsol;	/* start time of inertial navigation solution */
	double edsol;	/* end time of inertial navigation solution */
	int gweek;		/* gps week of imu data starting from */
	double Hz;		/* sampling frequency */

	double fre[3];	/* initial coordinate in ecef for forward filter */
	double fve[3];	/* initial velocity for forward filter */
	double fatt[3];	/* initial attitude for forward filter */

	double bre[3];	/* initial coordinate in ecef for backward filter */
	double bve[3];	/* initial velocity for backward filter */
	double batt[3];	/* initial attitude for backward filter */

	double ba, bg;	/* gyroscope and accelerometer bias */
	double sa, sg;	/* gyroscope and accelerometer scale factors */
	double lever[3];/* lever arm in body frame (RFU) */
	double corretime_ag;/* correlation time of ba bg */
	int align_type;		/* 0:att,1:coarse */
	//int soltype;        /* solution type (0:forward,1:backward,2:combined) */

	double maxvel;		/* max velocity under stationary state */

	int coupled_type;	/* loosely coupled or tightly coupled */
	int markov;			/* 0:off,1:on */
	int est_scale;		/* 0:off,1:on */
	int corr_imu;		/* 0:off,1:on */
	int corr_dcm;		/* 0:off,1:on */
	int zvu;		/* 0:off,1:on */
	int outsolstat;		/* 0:off,1:on */
	int outdebug;		/* 0:off,1:on */
	double corr_imu_k;	/* bias partical feedback coefficient */

	
	
	cprn_t prn;
}cateopt_t;

typedef struct {		/* calculated processing options */
	double gft;			/* gf threshold according gnss sample interval */
	double mwt;			/* mw threshold according gnss sample interval */
}calcuopt_t;

typedef struct {		/* input files */
	char ginfile[MAXEXFILE2][MAXSTRPATH];	/* 0:O, 1:N, 2:clk, 3:sp3 */
	char goutfile[MAXSTRPATH];				/* out file */

	char inimu[MAXSTRPATH];					/* path of imu data input */
	char inalign[MAXSTRPATH];				/* path of imu data for initial alignment */
	char insolv[MAXSTRPATH];				/* path of imu data for inertial solving */
	char solution[MAXSTRPATH];				/* path of inertial solution */
}filepath_t;

typedef struct {	/* imu measurement data record */
	gtime_t time;	/* imu sampling time (ins) */
	double gyro[3];	/* gyroscope measurements */
	double acce[3]; /* accelerometer measurements */
}imud_t;

typedef struct {	/* imu measurement data */
	int n, nmax;	/* number of imu data/allocated */
	imud_t *data;	/* imu measurement data records */
}imu_t;

typedef struct {
	gtime_t time;
	unsigned char stat;	/* CSOL_??? */
	unsigned char ns;	/* number of valid satellites */

	double re[3];
	double ve[3];
	double att[3];
	double ba[3];	/* bias */
	double bg[3];
	double sa[3];	/* scale factor */
	double sg[3];
	double off[6];	/* m */
	double drif;
	double trop;
	double iono[SATNUM];		/* vertical ionospheric delay L1(m) */
	double amb[NFREQ][SATNUM];

	float qr[6];	/* position variance/covariance (m^2) */
					/* {c_xx,c_yy,c_zz,c_xy,c_yz,c_zx} */
	float qv[6];	/* velocity variance/covariance (m^2/s^2) */
	float qa[6];	/* attitude angle variance/covariance (rad^2) */

}solvd_t;

typedef struct {				/* coupled satellite status */
	gtime_t pepoch;				/* previous epoch time */
	unsigned char vs;			/* 1:valid satellite flag single, CycSlip_Detec中重置，Coupled_IF_Res中更新 */

	int imw;					/* 当前卫星构成平均mw观测值的mw个数 */
	double gf;					/* 周跳探测前为前一历元的gf，周跳探测后为当前历元的gf */
	double mw;					/* 周跳探测前为之前所有历元的mw平均值，周跳探测后为当前历元与之前所有历元的mw平均值 */

	double resc_prio[NFREQ];		/* priori residuals of carrier-phase (m) */
	double resp_prio[NFREQ];		/* priori residuals of pseudorange (m) */
	double resd_prio[NFREQ];		/* priori residuals of doppler (m/s) */
	double resc_post[NFREQ];		/* psot residuals of carrier-phase (m) */
	double resp_post[NFREQ];		/* post residuals of pseudorange (m) */
	double resd_post[NFREQ];		/* post residuals of doppler (m/s)*/

	int lock;					/* lock counter of satellite */
	unsigned int outc;			/* obs outage counter of satellite */

}csat_t;

typedef struct {
	int nx;				/* state dimension, meas dimension */
	double tt;			/* time difference between current and previous (s) */
	double *xx, *P;	/* states  */
	csat_t csat[MAXSAT];
	cateopt_t opt;
	solvd_t solv;
	prcopt_t popt;
	solopt_t sopt;
	calcuopt_t copt;
}cate_t;


/* global variables ----------------------------------------------------------*/
extern const cateopt_t cateopt_default;
extern opt_t Coupled_SysOpts[];

/* functions ----------------------------------------------------------*/
/* Cate_IMU.c --------------------------------------------------------------------------*/
extern int Read_IMU(const char *infile, imu_t *imu, const cateopt_t *iopt);
extern int Divide_IMU(imu_t *imu, const cateopt_t *iopt, imu_t *initalign, imu_t *inersolv);
extern void Corr_IMUd(imud_t *data, solvd_t *solv, const double k);
extern int Proc_INS(const cateopt_t *opt, const filepath_t *flp, double *Cb2e);
extern int InputIMU(imu_t *imus, imud_t *imu, int revs, int *n, int solq, int flag);

/* Cate_Common_Fun.c -------------------------------------------------------------- */
extern void Vec_to_DiagMat(const double *a, double *A, int n);
extern void Mat_Add(int n, int m, const double a, const double *A, const double b, const double *B, const double c, double *C);
extern void	Mat_Trans(const double *A, int n, int m, double *B);
extern void	Mat_Norm(const double *A, int n, double *B);
extern void	Mat_Num(const double *A, const double B, int n, int m, double *C);
extern void	Mat_Skew(const double *vec, double *mat);
extern void Mat_Sub(double *A, const int n, const double *B, const int m, const int i, const int j, const double b);
extern void Mat_Sym(const int n, double *A);
extern void iMatfprint(const int A[], int n, int m, FILE *fp);
extern void iMatprint(const int A[], int n, int m);
extern void Out_Matrix(double *M, int n, int m);

/* attitude -------------------------------------------------------------------*/
extern void Pos_to_Cn2e(const double *pos, double *Cn2e);
extern void Cb2n_to_Att(const double *Cb2n, double *att);
extern void Att_to_Cb2n(const double *att, double *Cb2n);
extern void Cb2e_to_Att(const double *Cb2e, const double *pos, double *att);
extern void Att_to_Cb2e(const double *att, const double *pos, double *Cb2e);

/* common functions ----------------------------------------------------------*/
extern double randn();
extern void	Grav_ECEF(const double *r, double *ge);
extern void Init_Px(cate_t *cate, double xi, double var, int i);
extern void Init_SqMat(const int n, double *M, double var, int i);

extern int Kalman(double *xx, double *P, double *v, double *H, double *R, const int nx, const int nv);

extern void followopen(const char *path);
extern void followclose(void);
extern void follow(int level, const char *format, ...);
extern int showproc(char *format, ...);

/* Cate_Alignment.c --------------------------------------------------------------*/
extern int Init_Align(imu_t *imu, const double *re, double Hz, double *r);

/* Cate_NavPINS.c ----------------------------------------------------------------*/
extern void FinsPos(cate_t *cate, imud_t *data, double *Cb2e, const double tor_i, const int flag);
extern int PinsExec(cate_t *cate, double *Cb2e, imu_t *imusol);

/* Cate_Options.c -------------------------------------------------------------------- */
extern int Get_SysOpts(const char *path, opt_t *opts, prcopt_t *popt, solopt_t *sopt, filopt_t *fopt, 
	cateopt_t *opt, filepath_t *file);
extern void Reset_SysOpts();

/* Cate_CoupNav.c --------------------------------------------------------------------- */
extern int lscpnx(const cateopt_t *opt);
extern int lscp_pos(cate_t *cate, char *gnssol, double *Cb2e, imu_t *imusol);
extern int ttcpnx(const cateopt_t *opt, const prcopt_t *popt);
extern void TC_Pos(cate_t *cate, double *Cb2e, double *acce, double *gyro,
	rtk_t *rtk, obsd_t *obs, int n, const nav_t *nav);

/* Constraint */
extern int Constraint(cate_t *cate, double *Cb2e, double *v, double *H, double *R, const int nv);

/* Out_Out.c */
extern int File_Open(const cateopt_t *iopt, const char *is);
extern void File_Close(const cateopt_t *iopt);
extern int Out_Header(const char *outfile);
extern void cOutSols(solvd_t *solv);
extern void cOutSolStat(cate_t *cate, rtk_t *rtk);

/* Cate_QuaControl.c ------------------------------------------------------------------- */
extern void CycSlip_Thres(const double samp, double *gft, double *mwt);
extern void CycSlip_Detec(cate_t *ekf, rtk_t *rtk, const obsd_t *obs, const int n, const nav_t *nav);
extern int ObsMeasSPP_Detec(const prcopt_t *opt, obsd_t *obs, const int nobs);
extern int ObsMeasPPP_Detec(obsd_t *obs, const int nobs);
extern int ObsMeasDOP_Detec(const obsd_t *obs, const int nobs);
extern void Corr_Close(cate_t *cate, const cateopt_t *opt, solvd_t *solv, double *x, double *Cb2e, const int vsn);
extern void LeverG2I(double *re, double *ve, const double *Cb2e, const double *gyro, const double *lever);
extern void LeverI2G(double *re, double *ve, const double *Cb2e, const double *gyro, const double *lever);


/* Cate_UdSate.c -----------------------------------------------------------------------------*/
extern void UdState_INS(cate_t *cate, double *P, double *Cb2e, double *acce, double *gyro, const int nx);
extern void UdState_TTCP(cate_t *cate, rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav,
	double *Cb2e, double *acce, double *gyro);
extern void Udsol(rtk_t *rtk, const obsd_t *obs, const int n, cate_t *cate, double *Cb2e);

/* postpos.c -----------------------------------------------------------------------------*/
extern int cPostPos(prcopt_t *popt, solopt_t *sopt, const filopt_t *fopt,
					const cateopt_t *opt, filepath_t *flp, char **infile, 
					int *index, int n);


/* ppp.c -----------------------------------------------------------------------------*/
extern void Cate_testeclipse(const obsd_t *obs, int n, const nav_t *nav, double *rs);
extern int Cate_model_trop(gtime_t time, const double *pos, const double *azel,
							const prcopt_t *opt, const double *x, double *dtdx,
							const nav_t *nav, double *dtrp, double *var, cate_t *cate);
extern int Cate_model_iono(gtime_t time, const double *pos, const double *azel,
							const prcopt_t *opt, int sat, const double *x,
							const nav_t *nav, double *dion, double *var, cate_t *cate);
extern void Cate_satantpcv(const double *rs, const double *rr, const pcv_t *pcv, double *dant);
extern int Cate_model_phw(gtime_t time, int sat, const char *type, int opt,
							const double *rs, const double *rr, double *phw);
extern void Cate_corr_meas(const obsd_t *obs, const nav_t *nav, const double *azel,const prcopt_t *opt, 
							const double *dantr,const double *dants, double phw, double *L, double *P,
							double *Lc, double *Pc);

extern double Cate_varerr(int sat, int sys, double el, int freq, int type,
	const prcopt_t *opt);



#endif