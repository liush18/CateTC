/*------------------------------------------------------------------------------
 * Sate.h :
 *
 *			Copyright (C) 2019-2019 by Liu Shuai, All rights reserved.
 *
 * version :
 * history : 2019/11/01  new
 *------------------------------------------------------------------------------*/

#ifndef STATE_H
#define STATE_H

#define SQR(x)      ((x)*(x))

 /* (coupled) uncertain for P */
#define CVAR_POS		SQR(60.0)				/* init variance receiver position (m^2) */
#define CVAR_VEL		SQR(10.0)				/* init variance of receiver vel ((m/s)^2) */
#define CVAR_ATT		SQR(10.0*D2R)			/* deg to rad */
#define CVAR_BA			SQR(10.0*MG2AC)			/* mg to m/s^2 */
#define CVAR_BG			SQR(1.0*AS2R)			/* deg/hr to rad/s */
#define CVAR_SA			SQR(1E-4)
#define CVAR_SG			SQR(1E-4)
#define CVAR_CLK		SQR(60.0)				/* init variance receiver clock (m^2) */
#define CVAR_DRIF		SQR(10.0)				/* init variance receiver clock drift ((s/s)^2) */
#define CVAR_ZTD		SQR( 0.6)				/* init variance ztd (m^2) */
#define CVAR_GRA		SQR(0.01)				/* init variance gradient (m^2) */
#define CVAR_IONO		SQR(60.0)				/* init variance iono-delay (m^2) */
#define CVAR_DCB		SQR(30.0)				/* init variance dcb (m^2) */
#define CVAR_BIAS		SQR(60.0)				/* init variance phase-bias (m^2) */
#define cVAR_GLO_IFB	SQR( 0.6)				/* variance of glonass ifb */


/* (coupled) system noise default prn(process-noise std) for Q */
#define PRN_WA		(100*UG2AC)
#define PRN_WG		(0.1*DSH2RSS)
#define PRN_WBA		(1E-3)
#define PRN_WBG		(1E-3)
#define PRN_WSA		(1E-6)
#define PRN_WSG		(1E-6)

#define MEAS_POS	0.1		/* position measurement noise */
#define MEAS_VEL	0.05	/* velocity measurement noise */

#define MEAS_CP		0.009	/* carrier-phase noise */
#define MEAS_PR		0.9		/* pseudorange noise */
#define MEAS_DOP	0.09	/* doppler-freq noise */

/* (coupled) state numbers */
#define INP			3		/* delta pos */
#define INV			3		/* delta vel */
#define INA			3		/* delta attitude */
#define INBA		3		/* acce bias */
#define INBG		3		/* gyro bias */
#define INSA		3		/* acce scale */
#define INSG		3		/* gyro scale */
#define GNF(popt)   ((popt)->ionoopt==IONOOPT_IFLC?1:(popt)->nf) /* rtklib popt */
#define GNC			(NSYS)	/* clock off */
#define GNR			1		/* clock drift */
#define GNT(popt)	((popt)->tropopt<TROPOPT_EST?0:((popt)->tropopt==TROPOPT_EST?1:3))		/* rtklib opt, delta troposphere */
#define GNI(popt)	((popt)->ionoopt==IONOOPT_EST?SATNUM:0)	/* ionosphere, rtklib popt */
#define GNB(popt)	(GNF(popt)*SATNUM)/* bias , rtklib opt*/

#define CI_N(opt)	((opt)->est_scale?(INP+INV+INA+INBA+INBG+INSA+INSG):(INP+INV+INA+INBA+INBG)) /* cate opt */
#define CG_N(popt)	(GNC+GNR+GNT(popt)+GNI(popt)+GNB(popt))		/* rtklib opt */

/* (coupled) state index from ... */
#define IIP			(0)								/*  0  1  2 */
#define IIV			(INP)							/*  3  4  5 */
#define IIA			(INP+INV)						/*  6  7  8 */
#define IIBA		(INP+INV+INA)					/*  9 10 11 */
#define IIBG		(INP+INV+INA+INBA)				/* 12 13 14 */
#define IISA(opt)	((opt)->est_scale?(INP+INV+INA+INBA+INBG):(INP+INV+INA+INBA))			/* cate opt, 15 16 17 */
#define IISG(opt)	((opt)->est_scale?(INP+INV+INA+INBA+INBG+INSA):(INP+INV+INA+INBA+INBG))	/* cate opt, 18 19 20 */
#define GIC(s,opt)	(CI_N(opt)+s)					/* cate opt, s from 0, 21 22 or 15 16 */
#define GIR(opt)	(CI_N(opt)+GNC)					/* cate opt, 23 or 17 */
#define GIT(opt)	(CI_N(opt)+GNC+GNR)				/* cate opt, 24 or 18 */
#define GII(s,popt,opt)	(CI_N(opt)+GNC+GNR+GNT(popt)+(s)-1)	/* cate opt, from sat 1 */
#define GIB(s,f,popt,opt)(CI_N(opt)+GNC+GNR+GNT(popt)+GNI(popt)+SATNUM*(f)+(s)-1)	/* cate opt, s from sat 1, 25..56 or 19..50 */

#endif