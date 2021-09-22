/*------------------------------------------------------------------------------
* Cate_IMU.c : Read or output IMU 
*
*          Copyright (C) 2019-2019 by Liu Shuai, All rights reserved.
*
* version :
* history : 2019/04/17  new
*------------------------------------------------------------------------------*/
#include "Cate.h"

#define		MAXIMUBUFF		1024
#define		MAXIMUCOLM		1024		/* max number of imu data column */
#define		NINCIMU			262144		/* inclimental number of imu data */

/* add imu data to imu --------------------------------------------------------*/
static int Add_IMUd(imu_t *imu, imud_t *data)
{
	imud_t *imu_data;

	if (imu->nmax <= imu->n)
	{
		if (imu->nmax <= 0) imu->nmax = NINCIMU; else imu->nmax *= 2;
		if (!(imu_data = (imud_t *)realloc(imu->data, sizeof(imud_t) * imu->nmax)))
		{
			free(imu->data); imu->data = NULL; imu->n = imu->nmax = 0;
			follow(1, "Error:Add_IMUd-Reallocate memory!\n");
			return 0;
		}
		imu->data = imu_data;
	}

	imu->data[imu->n++] = *data;
	return 1;
}

/* read imu data --------------------------------------------------------*/
/* read imu file --------------------------------------------------------
* open imu file and read
* args   :	const char		*infile		I		imu file path
*			imu_t			*imu		IO		imu measurement data
			const iopt_t	*iopt		I		sins processing options
* return : states (1:ok, 0:error)
*------------------------------------------------------------------------*/
extern int Read_IMU(const char *infile, imu_t *imu, const cateopt_t *opt)
{
	FILE *fp;
	imud_t *data;
	char *p, buff[MAXIMUBUFF];
	double dat[MAXIMUCOLM] = { 0 };
	int i = 0, j = 0, stat = 1;
	int ntime = 0, ngyro = 1, nacce = 4;

	follow(2, "Read_imu...\n");

	if (!(fp = fopen(infile, "r"))) {
		follow(1, "Error:Read_IMU-Open IMU file=%s\n", infile);
		return 0;
	}

	if (!(data = (imud_t *)malloc(sizeof(imud_t)))) {
		follow(1, "Error:Read_IMU-Allocate memory, file=%s\n", infile);
		return 0;
	}

	while (fgets(buff, MAXIMUBUFF, fp) && stat > 0)
	{
		i = 0;
		p = strtok(buff, " ");

		dat[0] = atof(p);
		while (1)
		{
			i++;
			p = strtok(NULL, " ");
			if (p == NULL) break;
			dat[i] = atof(p);
		}
		//data->gpst.week = opt->gweek;
		//data->gpst.second = dat[ntime];
		data->time = gpst2time(opt->gweek, dat[ntime]);

		for (j = 0; j < 3; j++)
		{
			data->gyro[j] = dat[ngyro+j];
			data->acce[j] = dat[nacce+j];
		}

		stat = Add_IMUd(imu, data);
	}

	free(data); fclose(fp);
	return stat;
}

/* divide the data into static parts for initial alignment ----------------------
				   and dynamic parts for inertial navigation --------------------
* args   :	imu_t			*imu			IO		imu measurement data
*			const inpprcopt_t	*iopt			I		sins processing options
*			imu_t			*initalign		IO		data for initial alignment
*			imu_t			*inersolv		IO		data for inertial navigation
* return :	states (1:ok, 0:error)
*-------------------------------------------------------------------------------*/
extern int Divide_IMU(imu_t *imu, const cateopt_t *iopt, imu_t *initalign, imu_t *inersolv)
{
	imud_t *data;
	double tow;
	if (!(data = (imud_t *)malloc(sizeof(imu->data)))) { 
		printf("Error:Divide_IMU-Allocate memory\n");
		return 0; 
	}

	for (int i = 0; i < imu->n; i++)
	{
		tow = time2gpst(imu->data[i].time, NULL);
		if (tow >= (iopt->stini - 0.001) && tow <= (iopt->edini + 0.001))
		{
			if (!Add_IMUd(initalign, imu->data + i)) {
				printf("Error:Divide_IMU-Add %d imu!\n", i);
				return 0;
			}
		}
		if (tow >= (iopt->stsol - 0.001) && tow <= (iopt->edsol + 0.001))
		{
			if (!Add_IMUd(inersolv, imu->data + i)) {
				printf("Error:Divide_IMU-Add %d imu!\n", i);
				return 0;
			}
		}	
	}

	if (!(initalign->data = (imud_t *)realloc(initalign->data, sizeof(imud_t)*initalign->n)))
	{
		printf("Error:Divide_IMU-Reallocate memory of alignment imu!\n");
		free(initalign->data); initalign->data = NULL; initalign->n = initalign->nmax = 0;
		return 0;
	}

	if (!(inersolv->data = (imud_t *)realloc(inersolv->data, sizeof(imud_t)*inersolv->n)))
	{
		printf("Error:Divide_IMU-Reallocate memory of solving imu!\n");
		free(inersolv->data); inersolv->data = NULL; inersolv->n = inersolv->nmax = 0;
		return 0;
	}

	free(data);

	return 1;
}

/* correct imu -------------------------------------------------
* args   :	imud_t			*data		I		output file path
*			const inspar_t	*par		I		sins parameter setting
*-------------------------------------------------------------------------------*/
extern void Corr_IMUd(imud_t *data, solvd_t *solv, const double k)
{
	double fb[3], wb[3], IMa[9], IMg[9];
	int i, j;

	for (i = 0; i < 3; i++) {
		for (j = 0; j < 3; j++) {
			IMa[j + i * 3] = i == j ? (1 + solv->sa[i]) : 0.0;
			IMg[j + i * 3] = i == j ? (1 + solv->sg[i]) : 0.0;
		}
	}
	matinv(IMa, 3);
	matinv(IMg, 3);
	/* 部分校正 */
	Mat_Add(3, 1, 1.0, data->acce, -k, solv->ba, 0.0, fb);
	Mat_Add(3, 1, 1.0, data->gyro, -k, solv->bg, 0.0, wb);

	matmul("NN", 3, 1, 3, 1.0, IMa, fb, 0.0, data->acce);
	matmul("NN", 3, 1, 3, 1.0, IMg, wb, 0.0, data->gyro);

}

static char *solq2str(int solq)
{
	switch (solq) {
	case 0:return "no solution"; break;
	case 1:return "pure inertial navigation"; break;
	case 2:return "loosely coupled of spp and ins"; break;
	case 3:return "loosely coupled of ppp and ins"; break;
	case 4:return "tightly coupled of ppp and ins"; break;
	case 5:return "tightly coupled of ppp and ins with doppler"; break;
	default:return "error"; break;
	}
}

/* 
 * revs : (0:forward,1:backwoard)
 * flag : output processing or not 
 * 
 */
extern int InputIMU(imu_t *imus, imud_t *imu, int revs, int *n, int solq, int flag)
{
	if (*n >= imus->n || *n < 0) {
		return 0;
	}
	
	if (flag && 0 <= *n && *n < imus->n && (*n % 100 == 0)) {
		showmsg("coupled processing:%s, Q=%d\r", time_str((imus->data + *n)->time, 0), solq);
	}

	if (!revs) {
		*imu = *(imus->data + *n);
		*n = *n + 1;
	}
	else {
		*imu = *(imus->data + *n);
		*n = *n - 1;
	}

	return 1;
}