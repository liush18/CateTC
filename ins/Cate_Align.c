/*------------------------------------------------------------------------------
 * Cate_Align.c : Alignment
 *
 *			Copyright (C) 2019-2019 by Liu Shuai, All rights reserved.
 *
 * version :
 * history : 2019/04/22  new
 *------------------------------------------------------------------------------*/
#include "Cate.h"

#define ECEF 1

/* get inital alignment matrix in ecef coordinate -------------------------*/
static void Get_m_ECEF(const double *re, double Hz, double *me)
{
	double dt = 1.0 / Hz;
	double ge[3], wie2e[3] = { 0 }, gwe[3], gwge[3];
	
	Grav_ECEF(re, ge);
	wie2e[2] = OMGE;

	cross3(ge, wie2e, gwe);
	cross3(gwe, ge, gwge);

	/* me = [ge		ge*wie2e	ge*wie2e*ge] */
	for (int i = 0; i < 9; i++)
	{
		if (i < 3) me[i] = ge[i];
		else if (i < 6) me[i] = gwe[i - 3];
		else me[i] = gwge[i - 6];
	}

}

/* get inital alignment matrix in ENU coordinate -------------------------*/
static void Get_m_ENU(imu_t imu, double *mn)
{

}

/* get inital alignment matrix in body coordinate -------------------------*/
static void Get_m_Body(imu_t *imu, double *mb)
{
	int i, j;
	double fb[3] = { 0 }, wie2b[3] = { 0 }, fwb[3], fwfb[3];
	double temp[3];
	for (i = 0; i < imu->n; i++)
	{
		for (j = 0; j < 3; j++)
		{
			fb[j] = fb[j] + imu->data[i].acce[j];
			wie2b[j] = wie2b[j] + imu->data[i].gyro[j];
		}
	}
	Mat_Num(fb, 1.0 / imu->n, 3, 1, temp);
	for (i = 0; i < 3; i++) fb[i] = temp[i];
	Mat_Num(wie2b, 1.0 / imu->n, 3, 1, temp);
	for (i = 0; i < 3; i++) wie2b[i] = temp[i];

	cross3(fb, wie2b, fwb);
	cross3(fwb, fb, fwfb);

	/* mb = [-fb	-fb*wie2b	fb*wie2b*fb] */
	for (i = 0; i < 9; i++)
	{
		if (i < 3) mb[i] = -fb[i];
		else if (i < 6) mb[i] = -fwb[i - 3];
		else mb[i] = fwfb[i - 6];
	}

}


/* coarse alignment with analysis formula in ecef cooedinate system -------------------
*				me = rb2e*mb
----------------------------------------------------------------------------------------*/
static void Align_ECEF(imu_t *imu, const double *re, double Hz, double *Cb2e)
{
	double me[9], mb[9], temp[9];
	Get_m_ECEF(re, Hz, me);
	Get_m_Body(imu, mb);
	matinv(mb, 3);
	matmul("NN", 3, 3, 3, 1.0, me, mb, 0.0, temp);
	Mat_Norm(temp, 3, Cb2e);
	//matcpy(Cb2e, temp, 3, 3);
}

static void Align_ENU()
{

}

extern int Init_Align(imu_t *imu, const double *re, double Hz, double *DCM)
{
	follow(2, "Init_Align...\n");

#if ECEF
	Align_ECEF(imu, re, Hz, DCM);
#else
	Align_ENU();
#endif

	return 1;
}