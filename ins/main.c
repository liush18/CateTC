#if 1

#include "Cate.h"

/* (0:ok,0>:error,1:aborted) */
extern int Cate_Conf(const char *path)
{
	int i, n, index[MAXEXFILE] = { 0 }, stat = 0;
	char *infile[MAXEXFILE];
	prcopt_t popt;
	solopt_t sopt;
	filopt_t filopt = { "" };
	cateopt_t opt;
	filepath_t flp = { 0 };

	showproc("Cate_Conf...");

	Reset_SysOpts();
	if (!Get_SysOpts(path, Coupled_SysOpts, &popt, &sopt, &filopt, &opt, &flp)) {
		getchar();
		return 0;
	}

	for (i = 0, n = 0; i < MAXEXFILE2; i++) {
		if (*flp.ginfile[i]) {
			infile[i] = flp.ginfile[i];
			n++;
			index[i] = i;
		}
	}

	stat = cPostPos(&popt, &sopt, &filopt, &opt, &flp, infile, index, n);
	
	return stat;
}

int main()
{
	const char path[MAXSTRPATH] = ".\\zhuankou\\zhuankou.conf";

	Cate_Conf(path);

	showproc("\nPress any key to exit!");
	getchar();
	return 0;
}

#endif