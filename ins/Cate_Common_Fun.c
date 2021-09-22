/*------------------------------------------------------------------------------
* Cate_Common_Fun.c : Cate Common Functions
*
*          Copyright (C) 2019-2019 by Liu Shuai, All rights reserved.
*
* version :
* history : 2019/04/17  new
*------------------------------------------------------------------------------*/
#include "Cate.h"

const cateopt_t cateopt_default = { /* defaults processing options */
	{0,0.0},{0,0.0},0.0,0.0,0.0,

	0.0,0.0,0.0,0.0,0,100,

	0.0,0.0,0.0,		/* initial coordinate in ecef */
	0.0,0.0,0.0,		/* initial velocity */
	0.0,0.0,0.0,		/* initial attitude */
	0.0,0.0,0.0,		/* initial coordinate in ecef for backward filter  */
	0.0,0.0,0.0,		/* initial velocity for backward filter  */
	0.0,0.0,0.0,		/* initial attitude for backward filter  */
	0.0,0.0,0.0,0.0,	/* ba, bg, sa, sg */
	0.0,0.0,0.0,		/* lever arm */
	3600,				/* corretime_ag */
	ALIGNCOR,//0,			/* alignment type, ins soltype */

	0.2,				/* max velocity under stationary state */

	PMODE_LSCP, 1, 0, 	/* coupled_type, markov, est_scale */
	0, 0, 0,			/* corr_imu, corr_dcm, constraint */
	0, 0,				/* outsolstat, outdebug */
	1.0,				/* bias partical feedback coefficient */

	{1,100.0,0.1,1E-3,1E-3,1E-6,1E-6}				/* cprn_t */
};

/* functions ----------------------------------------------------------------*/

/* K(nx,nv)=P*H'*(H*P*H'+R)^-1  */
static int KF_K(const int nx, const int nv, const double *Pp, const double *H, const double *R, double *K)
{
	double *T3, *T4;
	T3 = zeros(nx, nv); T4 = zeros(nv, nv);

	matcpy(T4, R, nv, nv);
	matmul("NT", nx, nv, nx, 1.0, Pp, H, 0.0, T3);		/* P*H' */
	matmul("NN", nv, nv, nx, 1.0, H, T3, 1.0, T4);		/* H*P*H'+R */
	if (matinv(T4, nv) < 0) return 0;
	matmul("NN", nx, nv, nv, 1.0, T3, T4, 0.0, K);

	free(T3); free(T4);

	return 1;
}

static void KF_P(const int nx, const int nv, const double *K, const double *H, const double *R, const double *Pp, double *P)
{
	double *T1, *T2, *T3;
	T1 = eye(nx); T2 = zeros(nx, nx); T3 = zeros(nx, nv);

	matmul("NN", nx, nx, nv, -1.0, K, H, 1.0, T1);	/* T1=I-KH */
	matmul("NN", nx, nx, nx, 1.0, T1, Pp, 0.0, T2);	/* T2=(I-KH)Pp */

#if 1	/* Joseph 形式 */
	matmul("NT", nx, nx, nx, 1.0, T2, T1, 0.0, P);	/* P=(I-KH)Pp(I-KH)' */
	matmul("NN", nx, nv, nv, 1.0, K, R, 0.0, T3);	/* T3=KR */
	matmul("NT", nx, nx, nv, 1.0, T3, K, 1.0, P);	/* P=(I-KH)Pp(I-KH)'+KRK' */
#endif

	free(T1); free(T2); free(T3);
}

static int Kalman_(double *x, double *P, const double *v, const double *H, const double *R, const int nx, const int nv)
{
	double *Xp, *K, *Pp;
	Xp = zeros(nx, 1); K = zeros(nx, nv); Pp = zeros(nx, nx);

	matcpy(Pp, P, nx, nx);
	/* K(nx,nv)=P*H'*(H*P*H'+R)^-1  */
	if (!KF_K(nx, nv, Pp, H, R, K)) { return 0; }

	/* x=x+K*v */
	matmul("NN", nx, 1, nv, 1.0, K, v, 1.0, Xp);
	matcpy(x, Xp, nx, 1);

	/* p=(I-KH)P or Joseph形式 */
	KF_P(nx, nv, K, H, R, Pp, P);
	free(Xp); free(K); free(Pp);
	return 1;
}

extern int Kalman(double *xx, double *P, double *v, double *H, double *R, const int nx, const int nv)
{
	int i, j, k, *ix;
	double *x_, *P_, *H_, *HH;

	ix = imat(nx, 1);
	for (i = k = 0; i < nx; i++) {
		if (xx[i] != 0.0 && P[i + i * nx] > 0.0) {
			ix[k++] = i;
		}
	}

	x_ = zeros(k, 1);
	P_ = mat(k, k); H_ = mat(k, nv); HH = mat(nv, k);
	for (i = 0; i < k; i++) {
		for (j = 0; j < k; j++) {
			P_[i + j * k] = P[ix[i] + ix[j] * nx];
		}
		for (j = 0; j < nv; j++) { H_[i + j * k] = H[ix[i] + j * nx]; }
	}

	Mat_Trans(H_, k, nv, HH);

	if (!Kalman_(x_, P_, v, HH, R, k, nv)) {
		return 0;
	}

	for (i = 0; i < k; i++) {
		xx[ix[i]] = xx[ix[i]] + x_[i];
		for (j = 0; j < k; j++) {
			P[ix[i] + ix[j] * nx] = P_[i + j * k];
		}
	}

	free(ix); free(x_); free(P_); free(H_); free(HH);
	return 1;
}

/* vector to diagonal matrix --------------------------------------------
* args   : double *a	    I  vector A (n x 1)
*		   double *A		I  matrix B (n x n)
*		   int n			I  dimension of vector a
* return : none
*-----------------------------------------------------------------------------*/
extern void Vec_to_DiagMat(const double *a, double *A, int n)
{
	int i, j;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			A[j + i * n] = i == j ? a[i] : 0.0;
		}
	}
}

/* matrix addition --------------------------------------------------------
* matrix plus matrix
* args   : double *A,*B     I  matrix A (n x m), B (n x m)
*			int    n,m	    I  size of matrix A,B
*          double *C        IO matrix C (n x m)
* return : none
*-----------------------------------------------------------------------------*/
/* C=a*A+b*B+c*C */
extern void Mat_Add(int n, int m, const double a, const double *A, const double b, const double *B,
	const double c, double *C)
{
	int s = n * m;
	while (--s >= 0)
	{
		C[s] = a*A[s] + b*B[s] + c*C[s];
	}
}

/* transpose matrix --------------------------------------------------------
* matrix transposed B = A'
* args   :	double *A	    I  matrix A (n x m)
*			int    n,m	    I  size of matrix A,B    
*			double *B       IO matrix transposed B (m x n)
* return :	none
*-----------------------------------------------------------------------------*/
extern void Mat_Trans(const double *A, int n, int m, double *B)
{
	int s = n * m;
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < n; j++)
		{
			B[i + m * j] = A[n * i + j];
		}
	}
}

/* substitude a small square matrix into another large matrix(square or not) ----------------------
* substitude matrix
* args   :	double *A	    IO large matrix A (n x *)
*			int    n	    I  row of square matrix A
*			double *B       I  small matrix B (m x m)
*			int    m		I  size of square matrix B
*			int		i		I  from row i of A
*			int		j		I  from column j of B
*			int		b		I  coefficient of matrix B
* return :	none
*-----------------------------------------------------------------------------*/
extern void Mat_Sub(double *A, const int n, const double *B, const int m, const int i, const int j, const double b)
{
	int i1, j1, i2, j2;
	for (i1 = 0, i2 = i; i1 < m; i1++, i2++) {
		for (j1 = 0, j2 = j; j1 < m; j1++, j2++) {
			A[i2 + j2 * n] = b * B[i1 + j1 * m];
		}
	}
}

/* normalize matrix ---------------------------------------------------------
* normalize matrix
* args   :	double *A       I   matrix A (n x n)
			int    n		I	size of matrix A
*			double *B       O   normlized matrix (n x n) || b || = 1
* return :	none
*-----------------------------------------------------------------------------*/
extern void Mat_Norm(const double *A, int n, double *B)
{

	double *tA = zeros(n, n), *iA = zeros(n, n), *pA = zeros(n, n), *nA = zeros(n, n);
	matcpy(nA, A, n, n);
	for (int i = 0; i < 5; i++)
	{
		Mat_Trans(nA, n, n, tA); 
		matinv(tA, n);
		Mat_Add(n, n, 1.0, nA, 1.0, tA, 0.0, pA);
		Mat_Num(pA, 0.5, n, n, nA); 
	}
	matcpy(B, nA, n, n);
	free(tA); free(iA); free(pA); free(nA);

}

/* number and matrix --------------------------------------------------------
* multiply matrix by number
* args   :	double	*a		I   vector or matrix a (n x m)
*			double	b		I	number
*			int		n,m     I   number of rows and columns of vector or matrix
*			double	*C      O   product of number and matrix (n x m) 
* return :	none
*-----------------------------------------------------------------------------*/
extern void Mat_Num(const double *A, const double b, int n, int m, double *C)
{
	int s = n * m;
	while (--s >= 0)
	{
		C[s] = A[s] * b;
	}
}

/* Skew symmetric matrix of vector */
extern void Mat_Skew(const double *vec, double *mat)
{
	mat[0] = 0;			mat[3] = -vec[2];	mat[6] = vec[1];
	mat[1] = vec[2];	mat[4] = 0;			mat[7] = -vec[0];
	mat[2] = -vec[1];	mat[5] = vec[0];	mat[8] = 0;
}

/* matrix symmetry */
extern void Mat_Sym(const int n, double *A)
{
	double *B = zeros(n, n);
	Mat_Trans(A, n, n, B);
	Mat_Add(n, n, 0.5, A, 0.5, B, 0.0, A);
	free(B);
}

extern void iMatfprint(const int A[], int n, int m, FILE *fp)
{
	int i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < m; j++)
			fprintf(fp, " %d", A[i + j * n]);
		fprintf(fp, "\n");
	}
}

extern void iMatprint(const int A[], int n, int m)
{
	iMatfprint(A, n, m, stdout);
}

/* attitude -------------------------------------------------------------------*/

/* get Cn2e from postion(B,L) */
extern void Pos_to_Cn2e(const double *pos, double *Cn2e)
{
	double sinb = sin(pos[0]), cosb = cos(pos[0]);
	double sinl = sin(pos[1]), cosl = cos(pos[1]);
	Cn2e[0] = -sinl;	Cn2e[3] = -sinb * cosl;		Cn2e[6] = cosb * cosl;
	Cn2e[1] = cosl;		Cn2e[4] = -sinb * sinl;		Cn2e[7] = cosb * sinl;
	Cn2e[2] = 0;		Cn2e[5] = cosb;				Cn2e[8] = sinb;
}

/* attitude to matrix ------------------------------------*/
extern void Att_to_Cb2n(const double *att, double *Cb2n)
{
	double si = sin(att[0]);
	double sj = sin(att[1]);
	double sk = sin(att[2]);

	double ci = cos(att[0]);
	double cj = cos(att[1]);
	double ck = cos(att[2]);

	Cb2n[0] = cj * ck - si * sj*sk; Cb2n[3] = -ci * sk; Cb2n[6] = sj * ck + si * cj*sk;
	Cb2n[1] = cj * sk + si * sj*ck; Cb2n[4] = ci * ck;  Cb2n[7] = sj * sk - si * cj*ck;
	Cb2n[2] = -ci * sj;				Cb2n[5] = si;		Cb2n[8] = ci * cj;

}

/* convert matrix to attitude , att = [pitch, roll, yaw] */
extern void Cb2n_to_Att(const double *Cb2n, double *att)
{
	att[0] = asin(Cb2n[5]);
	if (fabs(Cb2n[5]) <= 0.999999)
	{
		att[1] = -atan2(Cb2n[2], Cb2n[8]);
		att[2] = -atan2(Cb2n[3], Cb2n[4]);
	}
	else
	{
		att[1] = atan2(Cb2n[6], Cb2n[0]);
		att[2] = 0;
	}
}

/* get attitude angle from Cb2e, , att = [pitch, roll, yaw] */
extern void Cb2e_to_Att(const double *Cb2e, const double *pos, double *att)
{
	double *Cn2e, *Cb2n;
	Cn2e = zeros(3, 3); Cb2n = zeros(3, 3);

	Pos_to_Cn2e(pos, Cn2e);
	matmul("TN", 3, 3, 3, 1.0, Cn2e, Cb2e, 0.0, Cb2n);		/* Cb2n=Ce2n*Cb2e */
	Mat_Norm(Cb2n, 3, Cb2n);
	Cb2n_to_Att(Cb2n, att);

	free(Cn2e); free(Cb2n);
	Cn2e = Cb2n = NULL;
}

/* get Cb2e from attitude angle */
extern void Att_to_Cb2e(const double *att, const double *pos, double *Cb2e)
{
	double *Cn2e, *Cb2n;

	Cn2e = zeros(3, 3); Cb2n = zeros(3, 3);
	Pos_to_Cn2e(pos, Cn2e);
	Att_to_Cb2n(att, Cb2n);
	matmul("NN", 3, 3, 3, 1.0, Cn2e, Cb2n, 0.0, Cb2e);

	free(Cn2e); free(Cb2n);
	Cn2e = Cb2n = NULL;
}

/* random numbers --------------------------------------------------------------
* generate random numbers from -1 to 1 
* return :	random number
*-------------------------------------------------------------------------------*/
extern double randn()
{
	double r;
	r = rand() % 10000;
	r = (r / 10000 - 0.5) / 0.5;
	return r;
}

/* gravity in ecef ----------------------------------------------------------
* gravity calculation in ecef coordinate
* args   :	double		*r			I		ecef position pointer
*			double		*ge			O		ecef gravity
* return :	none
*-----------------------------------------------------------------------------*/
extern void Grav_ECEF(const double *r, double *ge)
{
	double rho;
	double r1, r2, r4, r6;
	double t, t2, t4, t6;
	double a1, a2, a3, a4, a5;
	double b1, b2, b3;
	double c1, c2, c3, c4;
	double d1, d2, d3, d4;

	rho = SQRT(SQR(r[0]) + SQR(r[1]) + SQR(r[2]));

	r1 = RE_WGS84 / rho; 
	r2 = SQR(r1); 
	r4 = SQR(r2);
	r6 = r2 * r4;

	t = r[2] / rho;
	t2 = SQR(t);
	t4 = SQR(t2);
	t6 = t2 * t4;

	a1 = -FM / SQR(rho);
	a2 = 1 + 3.0 / 2.0*J2*r2 - 15.0 / 8.0*J4*r4 + 35.0 / 16.0*J6*r6;
	a3 = -9.0 / 2.0*J2*r2 + 75.0 / 4.0*J4*r4 - 735.0 / 16.0*J6*r6;
	a4 = -175.0 / 8.0*J4*r4 + 2205.0 / 16.0*J6*r6;
	a5 = -1617.0 / 16.0*J6*r6;

	b1 = 3.0*J2*r2 - 15.0 / 2.0*J4*r4 + 105.0 / 8.0*J6*r6;
	b2 = 35.0 / 2.0*J4*r4 - 945.0 / 12.0*J6*r6;
	b3 = 693.0 / 8.0*J6*r6;

	c1 = a2; 
	c2 = a3 - b1; 
	c3 = a4 - b2; 
	c4 = a5 - b3;

	d1 = a2 + b1;
	d2 = c2 + b2;
	d3 = c3 + b3;
	d4 = c4;

	ge[0] = a1 / rho * (c1 + c2 * t2 + c3 * t4 + c4 * t6)*r[0] + SQR(OMGE)*r[0];
	ge[1] = a1 / rho * (c1 + c2 * t2 + c3 * t4 + c4 * t6)*r[1] + SQR(OMGE)*r[1];
	ge[2] = a1 / rho * (d1 + d2 * t2 + d3 * t4 + d4 * t6)*r[2];
}

/* initialize state and covariance -------------------------------------------*/
extern void Init_Px(cate_t *cate, double xi, double var, int i)
{
	int j;
	cate->xx[i] = xi;	/* xx 用作记录，不做状态进行估计 */
	for (j = 0; j < cate->nx; j++) {
		cate->P[i + j * cate->nx] = cate->P[j + i * cate->nx] = i == j ? var : 0.0;
	}
}

/* initialize square matrix, such as Q -------------------------------------------*/
extern void Init_SqMat(const int n, double *M, double var, int i)
{
	int j;
	for (j = 0; j < n; j++) {
		M[i + j * n] = M[j + i * n] = i == j ? var : 0.0;
	}
}

static int om = 0;
extern void Out_Matrix(double *M, int n, int m)
{

	//if (om == 0) {
	//	char *path = ".\\data\\ins-output\\M.txt";
	//	remove(path);
	//	om++;
	//}

	FILE *fpx1;
	if (fpx1 = fopen(".\\zhuankou\\cpnav\\M.txt", "w")) {
		matfprint(M, n, m, 50, 40, fpx1); fprintf(fpx1, "\n");
		fclose(fpx1);
	}
	

}
