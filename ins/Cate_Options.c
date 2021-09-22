/*------------------------------------------------------------------------------
* Cate_Options.c : load coupled options
*
*          Copyright (C) 2019-2019 by Liu Shuai, All rights reserved.
*
* version : 
* history : 2019/09/09  new
*------------------------------------------------------------------------------*/

#include "Cate.h"
#include "State.h"

static cateopt_t iopt_;
static filepath_t flp_;
static char gts_[1024];
static char gte_[1024];
static double ba_, bg_;
static double uncatt_, uncba_, uncbg_;
static double prnwa_, prnwg_;

static prcopt_t prcopt_;
static solopt_t solopt_;
static filopt_t filopt_;
static int antpostype_[2];
static double elmask_, elmaskar_, elmaskhold_;
static double antpos_[2][3];
static char exsats_[1024];
static char snrmask_[NFREQ][1024];

/* system options table ------------------------------------------------------*/
#define SEPOPT	"0:space,1:comma"
#define ALIOPT	"0:att,1:coarse"		/* alignment opt */
#define CPTOPT	"0:gnss,1:fins,2:lscp,3:ttcp"	/* coupled navigation type */

#define SWTOPT  "0:off,1:on"
#define MODOPT  "0:single,1:dgps,2:kinematic,3:static,4:movingbase,5:fixed,6:ppp-kine,7:ppp-static,8:ppp-fixed"
#define FRQOPT  "1:l1,2:l1+l2,3:l1+l2+l5,4:l1+l5"
#define TYPOPT  "0:forward,1:backward,2:combined"
#define IONOPT  "0:off,1:brdc,2:sbas,3:dual-freq,4:est-stec,5:ionex-tec,6:qzs-brdc,7:qzs-lex,8:stec"
#define TRPOPT  "0:off,1:saas,2:sbas,3:est-ztd,4:est-ztdgrad,5:ztd"
#define EPHOPT  "0:brdc,1:precise,2:brdc+sbas,3:brdc+ssrapc,4:brdc+ssrcom"
#define NAVOPT  "1:gps+2:sbas+4:glo+8:gal+16:qzs+32:comp"
#define GAROPT  "0:off,1:on,2:autocal"
#define SOLOPT  "0:llh,1:xyz,2:enu,3:nmea"
#define TSYOPT  "0:gpst,1:utc,2:jst"
#define TFTOPT  "0:tow,1:hms"
#define DFTOPT  "0:deg,1:dms"
#define HGTOPT  "0:ellipsoidal,1:geodetic"
#define GEOOPT  "0:internal,1:egm96,2:egm08_2.5,3:egm08_1,4:gsi2000"
#define STAOPT  "0:all,1:single"
#define STSOPT  "0:off,1:state,2:residual"
#define ARMOPT  "0:off,1:continuous,2:instantaneous,3:fix-and-hold"
#define POSOPT  "0:llh,1:xyz,2:single,3:posfile,4:rinexhead,5:rtcm,6:raw"
#define TIDEOPT "0:off,1:on,2:otl"
#define PHWOPT  "0:off,1:on,2:precise"


extern opt_t Coupled_SysOpts[] = {
	/* gnss options (0:int,1:double,2:string,3:enum) */
	{"pos0-start",			2,	(void *)&gts_,					""		},
	{"pos0-end",			2,	(void *)&gte_,					""		},
	{"pos0-procint",		1,	(void *)&iopt_.ti,				""		},
	{"pos0-procunit",		1,	(void *)&iopt_.tu,				""		},
	{"pos0-sampint",		1,	(void *)&iopt_.si,				""		},

	{"infile1-obs",			2,	(void *)&flp_.ginfile[0],		""		},
	{"infile1-nav",			2,	(void *)&flp_.ginfile[1],		""		},
	{"infile1-clk",			2,	(void *)&flp_.ginfile[2],		""		},
	{"infile1-sp3",			2,	(void *)&flp_.ginfile[3],		""		},
	{"outfile1-sol",		2,	(void *)&flp_.goutfile,			""		},

	/* ins options */
	{"time-align_start",	1,	(void *)&iopt_.stini,			""		},
	{"time-align_end",		1,	(void *)&iopt_.edini,			""		},
	{"time-solv_start",		1,	(void *)&iopt_.stsol,			""		},
	{"time-solv_end",		1,	(void *)&iopt_.edsol,			""		},
	{"time-start_week",		0,	(void *)&iopt_.gweek,			""		},
	{"ins-frequency",		1,	(void *)&iopt_.Hz,				""		},

	{"ins-init_fpos1",		1,	(void *)&iopt_.fre[0],			""		},
	{"ins-init_fpos2",		1,	(void *)&iopt_.fre[1],			""		},
	{"ins-init_fpos3",		1,	(void *)&iopt_.fre[2],			""		},
	{"ins-init_fvel1",		1,	(void *)&iopt_.fve[0],			""		},
	{"ins-init_fvel2",		1,	(void *)&iopt_.fve[1],			""		},
	{"ins-init_fvel3",		1,	(void *)&iopt_.fve[2],			""		},
	{"ins-init_fattp",		1,	(void *)&iopt_.fatt[0],			""		},
	{"ins-init_fattr",		1,	(void *)&iopt_.fatt[1],			""		},
	{"ins-init_fatty",		1,	(void *)&iopt_.fatt[2],			""		},

	{"ins-init_bpos1",		1,	(void *)&iopt_.bre[0],			""		},
	{"ins-init_bpos2",		1,	(void *)&iopt_.bre[1],			""		},
	{"ins-init_bpos3",		1,	(void *)&iopt_.bre[2],			""		},
	{"ins-init_bvel1",		1,	(void *)&iopt_.bve[0],			""		},
	{"ins-init_bvel2",		1,	(void *)&iopt_.bve[1],			""		},
	{"ins-init_bvel3",		1,	(void *)&iopt_.bve[2],			""		},
	{"ins-init_battp",		1,	(void *)&iopt_.batt[0],			""		},
	{"ins-init_battr",		1,	(void *)&iopt_.batt[1],			""		},
	{"ins-init_batty",		1,	(void *)&iopt_.batt[2],			""		},

	{"ins-bias_acce",		1,	(void *)&ba_,					""		},
	{"ins-bias_gyro",		1,	(void *)&bg_,					""		},
	{"ins-scal_acce",		1,	(void *)&iopt_.sa,				""		},
	{"ins-scal_gyro",		1,	(void *)&iopt_.sg,				""		},

	{"ins-leverRight",		1,	(void *)&iopt_.lever[0],		""		},
	{"ins-leverFront",		1,	(void *)&iopt_.lever[1],		""		},
	{"ins-leverUp",			1,	(void *)&iopt_.lever[2],		""		},

	{"ins-corretime",		1,	(void *)&iopt_.corretime_ag,	""		},
	{"ins-align_type",		3,	(void *)&iopt_.align_type,		ALIOPT	},
	//{"ins-soltype",			3,	(void *)&iopt_.soltype,			TYPOPT	},

	{"ins2-maxvel",			1,	(void *)&iopt_.maxvel,			""		},

	{"coupled-type",		3,	(void *)&iopt_.coupled_type,	CPTOPT	},
	{"coupled-markov",		3,	(void *)&iopt_.markov,			SWTOPT	},
	{"coupled-est_scale",	3,	(void *)&iopt_.est_scale,		SWTOPT	},
	{"coupled-corr_imu",	3,	(void *)&iopt_.corr_imu,		SWTOPT	},
	{"coupled-corr_dcm",	3,	(void *)&iopt_.corr_dcm,		SWTOPT	},
	{"coupled-zvu",			3,	(void *)&iopt_.zvu,				SWTOPT	},
	{"coupled-outsolstat",	3,	(void *)&iopt_.outsolstat,		SWTOPT	},
	{"coupled-outdebug",	3,	(void *)&iopt_.outdebug,		SWTOPT	},
	{"coupled-coor_imu_k",	1,	(void *)&iopt_.corr_imu_k,		""		},

	{"prn-flag",			3,	(void *)&iopt_.prn.flag,		SWTOPT	},
	{"prn-wa",				1,	(void *)&prnwa_,				""		},
	{"prn-wg",				1,	(void *)&prnwg_,				""		},
	{"prn-wba",				1,	(void *)&iopt_.prn.wba,			""		},
	{"prn-wbg",				1,	(void *)&iopt_.prn.wbg,			""		},
	{"prn-wsa",				1,	(void *)&iopt_.prn.wsa,			""		},
	{"prn-wsg",				1,	(void *)&iopt_.prn.wsg,			""		},

	{"infile2-imu_all",		2,	(void *)&flp_.inimu,			""		},
	{"infile2-imu_align",	2,	(void *)&flp_.inalign,			""		},
	{"infile2-imu_solv",	2,	(void *)&flp_.insolv,			""		},
	{"outfile2-solution",	2,	(void *)&flp_.solution,			""		},

	/* rtklib */
	{"pos1-posmode",    3,  (void *)&prcopt_.mode,       MODOPT },
	{"pos1-frequency",  3,  (void *)&prcopt_.nf,         FRQOPT },
	{"pos1-soltype",    3,  (void *)&prcopt_.soltype,    TYPOPT },
	{"pos1-elmask",     1,  (void *)&elmask_,            "deg"  },
	{"pos1-snrmask_r",  3,  (void *)&prcopt_.snrmask.ena[0],SWTOPT},
	{"pos1-snrmask_b",  3,  (void *)&prcopt_.snrmask.ena[1],SWTOPT},
	{"pos1-snrmask_L1", 2,  (void *)snrmask_[0],         ""     },
	{"pos1-snrmask_L2", 2,  (void *)snrmask_[1],         ""     },
	{"pos1-snrmask_L5", 2,  (void *)snrmask_[2],         ""     },
	{"pos1-dynamics",   3,  (void *)&prcopt_.dynamics,   SWTOPT },
	{"pos1-tidecorr",   3,  (void *)&prcopt_.tidecorr,   TIDEOPT},
	{"pos1-ionoopt",    3,  (void *)&prcopt_.ionoopt,    IONOPT },
	{"pos1-tropopt",    3,  (void *)&prcopt_.tropopt,    TRPOPT },
	{"pos1-sateph",     3,  (void *)&prcopt_.sateph,     EPHOPT },
	{"pos1-posopt1",    3,  (void *)&prcopt_.posopt[0],  SWTOPT },
	{"pos1-posopt2",    3,  (void *)&prcopt_.posopt[1],  SWTOPT },
	{"pos1-posopt3",    3,  (void *)&prcopt_.posopt[2],  PHWOPT },
	{"pos1-posopt4",    3,  (void *)&prcopt_.posopt[3],  SWTOPT },
	{"pos1-posopt5",    3,  (void *)&prcopt_.posopt[4],  SWTOPT }, /* raim_fde */
	{"pos1-posopt6",    3,  (void *)&prcopt_.posopt[5],  SWTOPT },
	{"pos1-exclsats",   2,  (void *)exsats_,             "prn ..."},
	{"pos1-navsys",     0,  (void *)&prcopt_.navsys,     NAVOPT },

	{"pos2-armode",     3,  (void *)&prcopt_.modear,     ARMOPT },
	{"pos2-gloarmode",  3,  (void *)&prcopt_.glomodear,  GAROPT },
	{"pos2-bdsarmode",  3,  (void *)&prcopt_.bdsmodear,  SWTOPT },
	{"pos2-arthres",    1,  (void *)&prcopt_.thresar[0], ""     },
	{"pos2-arthres1",   1,  (void *)&prcopt_.thresar[1], ""     },
	{"pos2-arthres2",   1,  (void *)&prcopt_.thresar[2], ""     },
	{"pos2-arthres3",   1,  (void *)&prcopt_.thresar[3], ""     },
	{"pos2-arthres4",   1,  (void *)&prcopt_.thresar[4], ""     },
	{"pos2-arlockcnt",  0,  (void *)&prcopt_.minlock,    ""     },
	{"pos2-arelmask",   1,  (void *)&elmaskar_,          "deg"  },
	{"pos2-arminfix",   0,  (void *)&prcopt_.minfix,     ""     },
	{"pos2-armaxiter",  0,  (void *)&prcopt_.armaxiter,  ""     },
	{"pos2-elmaskhold", 1,  (void *)&elmaskhold_,        "deg"  },
	{"pos2-aroutcnt",   0,  (void *)&prcopt_.maxout,     ""     },
	{"pos2-maxage",     1,  (void *)&prcopt_.maxtdiff,   "s"    },
	{"pos2-syncsol",    3,  (void *)&prcopt_.syncsol,    SWTOPT },
	{"pos2-slipthres",  1,  (void *)&prcopt_.thresslip,  "m"    },
	{"pos2-rejionno",   1,  (void *)&prcopt_.maxinno,    "m"    },
	{"pos2-rejgdop",    1,  (void *)&prcopt_.maxgdop,    ""     },
	{"pos2-niter",      0,  (void *)&prcopt_.niter,      ""     },
	{"pos2-baselen",    1,  (void *)&prcopt_.baseline[0],"m"    },
	{"pos2-basesig",    1,  (void *)&prcopt_.baseline[1],"m"    },

	{"out-solformat",   3,  (void *)&solopt_.posf,       SOLOPT },
	{"out-outhead",     3,  (void *)&solopt_.outhead,    SWTOPT },
	{"out-outopt",      3,  (void *)&solopt_.outopt,     SWTOPT },
	{"out-outvel",      3,  (void *)&solopt_.outvel,     SWTOPT },
	{"out-timesys",     3,  (void *)&solopt_.times,      TSYOPT },
	{"out-timeform",    3,  (void *)&solopt_.timef,      TFTOPT },
	{"out-timendec",    0,  (void *)&solopt_.timeu,      ""     },
	{"out-degform",     3,  (void *)&solopt_.degf,       DFTOPT },
	{"out-fieldsep",    2,  (void *)solopt_.sep,        ""     },
	{"out-outsingle",   3,  (void *)&prcopt_.outsingle,  SWTOPT },
	{"out-maxsolstd",   1,  (void *)&solopt_.maxsolstd,  "m"    },
	{"out-height",      3,  (void *)&solopt_.height,     HGTOPT },
	{"out-geoid",       3,  (void *)&solopt_.geoid,      GEOOPT },
	{"out-solstatic",   3,  (void *)&solopt_.solstatic,  STAOPT },
	{"out-nmeaintv1",   1,  (void *)&solopt_.nmeaintv[0],"s"    },
	{"out-nmeaintv2",   1,  (void *)&solopt_.nmeaintv[1],"s"    },
	{"out-outstat",     3,  (void *)&solopt_.sstat,      STSOPT },

	{"stats-eratio1",   1,  (void *)&prcopt_.eratio[0],  ""     },
	{"stats-eratio2",   1,  (void *)&prcopt_.eratio[1],  ""     },
	{"stats-errphase",  1,  (void *)&prcopt_.err[1],     "m"    },
	{"stats-errphaseel",1,  (void *)&prcopt_.err[2],     "m"    },
	{"stats-errphasebl",1,  (void *)&prcopt_.err[3],     "m/10km"},
	{"stats-errdoppler",1,  (void *)&prcopt_.err[4],     "Hz"   },
	{"stats-stdbias",   1,  (void *)&prcopt_.std[0],     "m"    },
	{"stats-stdiono",   1,  (void *)&prcopt_.std[1],     "m"    },
	{"stats-stdtrop",   1,  (void *)&prcopt_.std[2],     "m"    },
	{"stats-prnaccelh", 1,  (void *)&prcopt_.prn[3],     "m/s^2"},
	{"stats-prnaccelv", 1,  (void *)&prcopt_.prn[4],     "m/s^2"},
	{"stats-prnbias",   1,  (void *)&prcopt_.prn[0],     "m"    },
	{"stats-prniono",   1,  (void *)&prcopt_.prn[1],     "m"    },
	{"stats-prntrop",   1,  (void *)&prcopt_.prn[2],     "m"    },
	{"stats-prnpos",    1,  (void *)&prcopt_.prn[5],     "m"    },
	{"stats-clkstab",   1,  (void *)&prcopt_.sclkstab,   "s/s"  },

	{"ant1-postype",    3,  (void *)&antpostype_[0],     POSOPT },
	{"ant1-pos1",       1,  (void *)&antpos_[0][0],      "deg|m"},
	{"ant1-pos2",       1,  (void *)&antpos_[0][1],      "deg|m"},
	{"ant1-pos3",       1,  (void *)&antpos_[0][2],      "m|m"  },
	{"ant1-anttype",    2,  (void *)prcopt_.anttype[0],  ""     },
	{"ant1-antdele",    1,  (void *)&prcopt_.antdel[0][0],"m"   },
	{"ant1-antdeln",    1,  (void *)&prcopt_.antdel[0][1],"m"   },
	{"ant1-antdelu",    1,  (void *)&prcopt_.antdel[0][2],"m"   },

	{"ant2-postype",    3,  (void *)&antpostype_[1],     POSOPT },
	{"ant2-pos1",       1,  (void *)&antpos_[1][0],      "deg|m"},
	{"ant2-pos2",       1,  (void *)&antpos_[1][1],      "deg|m"},
	{"ant2-pos3",       1,  (void *)&antpos_[1][2],      "m|m"  },
	{"ant2-anttype",    2,  (void *)prcopt_.anttype[1],  ""     },
	{"ant2-antdele",    1,  (void *)&prcopt_.antdel[1][0],"m"   },
	{"ant2-antdeln",    1,  (void *)&prcopt_.antdel[1][1],"m"   },
	{"ant2-antdelu",    1,  (void *)&prcopt_.antdel[1][2],"m"   },
	{"ant2-maxaveep",   0,  (void *)&prcopt_.maxaveep    ,""    },
	{"ant2-initrst",    3,  (void *)&prcopt_.initrst,    SWTOPT },

	{"misc-timeinterp", 3,  (void *)&prcopt_.intpref,    SWTOPT },
	{"misc-sbasatsel",  0,  (void *)&prcopt_.sbassatsel, "0:all"},
	{"misc-rnxopt1",    2,  (void *)prcopt_.rnxopt[0],   ""     },
	{"misc-rnxopt2",    2,  (void *)prcopt_.rnxopt[1],   ""     },
	{"misc-pppopt",     2,  (void *)prcopt_.pppopt,      ""     },

	{"file-satantfile", 2,  (void *)&filopt_.satantp,    ""     },
	{"file-rcvantfile", 2,  (void *)&filopt_.rcvantp,    ""     },
	{"file-staposfile", 2,  (void *)&filopt_.stapos,     ""     },
	{"file-geoidfile",  2,  (void *)&filopt_.geoid,      ""     },
	{"file-ionofile",   2,  (void *)&filopt_.iono,       ""     },
	{"file-dcbfile",    2,  (void *)&filopt_.dcb,        ""     },
	{"file-eopfile",    2,  (void *)&filopt_.eop,        ""     },
	{"file-blqfile",    2,  (void *)&filopt_.blq,        ""     },
	{"file-tempdir",    2,  (void *)&filopt_.tempdir,    ""     },
	{"file-geexefile",  2,  (void *)&filopt_.geexe,      ""     },
	{"file-solstatfile",2,  (void *)&filopt_.solstat,    ""     },
	{"file-tracefile",  2,  (void *)&filopt_.trace,      ""     },

	{"",0,NULL,""}
};

/* string to enum ------------------------------------------------------------*/
static int Str_to_Enum(const char *str, const char *comment, int *val)
{
	const char *p;
	char s[32];

	for (p = comment;; p++) {
		if (!(p = strstr(p, str))) break;
		if (*(p - 1) != ':') continue;
		for (p -= 2; '0' <= *p&&*p <= '9'; p--);
		return sscanf(p + 1, "%d", val) == 1;
	}
	sprintf(s, "%30.30s:", str);
	if ((p = strstr(comment, s))) { /* number */
		return sscanf(p, "%d", val) == 1;
	}
	return 0;
}

/* discard space characters at tail ------------------------------------------*/
static void Chop(char *str)
{
	char *p;
	if ((p = strchr(str, '#'))) *p = '\0'; /* comment */
	for (p = str + strlen(str) - 1; p >= str && !isgraph((int)*p); p--) *p = '\0';
}

static opt_t *Search_SysOpts(const char *name, const opt_t *opts)
{
	int i;

	for (i = 0; *opts[i].name; i++) {
		if (strstr(opts[i].name, name)) return (opt_t *)(opts + i);
	}
	return NULL;
}

static int Str_to_SysOpts(opt_t *opt, const char *str)
{
	switch (opt->format) {
	case 0: *(int    *)opt->var = atoi(str); break;
	case 1: *(double *)opt->var = atof(str); break;
	case 2: strcpy((char *)opt->var, str);  break;
	case 3: return Str_to_Enum(str, opt->comment, (int *)opt->var);
	default: return 0;
	}
	return 1;
}

/* load coupled options -------------------------------------*/
static int Load_SysOpts(const char *file, opt_t *opts)
{
	FILE *fp;
	opt_t *opt;
	char buff[2048], *p;
	int n = 0;

	showproc("Load_SysOpts...");

	if (!(fp = fopen(file, "r"))) {
		follow(1, "ERROR : Load_Coupled_SysOpts: options file open error (%s)\n", file);
		return 0;
	}
	while (fgets(buff, sizeof(buff), fp)) {
		n++;
		Chop(buff);

		if (buff[0] == '\0') continue;

		if (!(p = strstr(buff, "="))) {
			follow(1, "WARNING : Invalid option %s (%s:%d)\n", buff, file, n);
			continue;
		}
		*p++ = '\0';
		Chop(buff);
		if (!(opt = Search_SysOpts(buff, opts))) continue;

		if (!str2opt(opt, p)) {
			follow(1, "WARNING : Invalid option value %s (%s:%d)\n", buff, file, n);
			continue;
		}
	}
	fclose(fp);

	return 1;
}

/* system options buffer to options ------------------------------------------*/
static void Buff_to_SysOpts()
{
	double epoch[6] = { 0.0 };
	double pos[3], *rr;
	char buff[1024], *p, *id;
	int i=0, j, sat, *ps;

	sscanf(gts_, "%d %lf/%lf/%lf %lf:%lf:%lf", &i, epoch + 0, epoch + 1, epoch + 2, 
		epoch + 3, epoch + 4, epoch + 5);
	if (i) { iopt_.ts = epoch2time(epoch); }

	sscanf(gte_, "%d %lf/%lf/%lf %lf:%lf:%lf", &i, epoch + 0, epoch + 1, epoch + 2, 
		epoch + 3, epoch + 4, epoch + 5);
	if (i) { iopt_.te = epoch2time(epoch); }

	iopt_.ba = ba_ * MG2AC;
	iopt_.bg = bg_ * DH2RS;

	if (iopt_.prn.flag == 1) {
		iopt_.prn.wa = prnwa_ * UG2AC;		/* ug/sqrt(hz) to m/s^2/sqrt(hz) */
		iopt_.prn.wg = prnwg_ * DSH2RSS;	/* deg/sqrt(h) to rad/sqrt(s) */
	}
	else {
		iopt_.prn.wa = PRN_WA;iopt_.prn.wg = PRN_WG;
		iopt_.prn.wba = PRN_WBA;iopt_.prn.wbg = PRN_WBG;
		iopt_.prn.wsa = PRN_WSA;iopt_.prn.wsg = PRN_WSG;
	}

	/* rtklib options */
	prcopt_.elmin = elmask_ * D2R;
	prcopt_.elmaskar = elmaskar_ * D2R;
	prcopt_.elmaskhold = elmaskhold_ * D2R;

	for (i = 0; i < 2; i++) {
		ps = i == 0 ? &prcopt_.rovpos : &prcopt_.refpos;
		rr = i == 0 ? prcopt_.ru : prcopt_.rb;

		if (antpostype_[i] == 0) { /* lat/lon/hgt */
			*ps = 0;
			pos[0] = antpos_[i][0] * D2R;
			pos[1] = antpos_[i][1] * D2R;
			pos[2] = antpos_[i][2];
			pos2ecef(pos, rr);
		}
		else if (antpostype_[i] == 1) { /* xyz-ecef */
			*ps = 0;
			rr[0] = antpos_[i][0];
			rr[1] = antpos_[i][1];
			rr[2] = antpos_[i][2];
		}
		else *ps = antpostype_[i] - 1;
	}
	/* excluded satellites */
	for (i = 0; i < MAXSAT; i++) prcopt_.exsats[i] = 0;
	if (exsats_[0] != '\0') {
		strcpy(buff, exsats_);
		for (p = strtok(buff, " "); p; p = strtok(NULL, " ")) {
			if (*p == '+') id = p + 1; else id = p;
			if (!(sat = satid2no(id))) continue;
			prcopt_.exsats[sat - 1] = *p == '+' ? 2 : 1;
		}
	}
	/* snrmask */
	for (i = 0; i < NFREQ; i++) {
		for (j = 0; j < 9; j++) prcopt_.snrmask.mask[i][j] = 0.0;
		strcpy(buff, snrmask_[i]);
		for (p = strtok(buff, ","), j = 0; p&&j < 9; p = strtok(NULL, ",")) {
			prcopt_.snrmask.mask[i][j++] = atof(p);
		}
	}
	/* number of frequency (4:L1+L5) */
	if (prcopt_.nf == 4) {
		prcopt_.nf = 3;
		prcopt_.freqopt = 1;
	}
}

/* return 0:err, 1:ok */
extern int Get_SysOpts(const char *path, opt_t *opts, prcopt_t *popt, solopt_t *sopt, filopt_t *fopt, 
	cateopt_t *iopt, filepath_t *file)
{
	double epoch[6] = { 0.0 };
	int flag = 0;

	showproc("Get_Coupled_SysOpts : file = %s", path);

	if (!Load_SysOpts(path, Coupled_SysOpts)) return 0;

	Buff_to_SysOpts();

	*popt = prcopt_;
	*sopt = solopt_;
	*fopt = filopt_;
	*iopt = iopt_;
	*file = flp_;

	return 1;
}

extern void Reset_SysOpts()
{
	int i, j;

	showproc("Reset_SysOpts...");

	iopt_ = cateopt_default;
	for (i = 0; i < MAXEXFILE2; i++) flp_.ginfile[i][0] = '\0';
	flp_.goutfile[0] = '\0';
	flp_.inimu[0] = '\0';
	flp_.inalign[0] = '\0';
	flp_.insolv[0] = '\0';
	flp_.solution[0] = '\0';

	prcopt_ = prcopt_default;
	solopt_ = solopt_default;
	filopt_.satantp[0] = '\0';
	filopt_.rcvantp[0] = '\0';
	filopt_.stapos[0] = '\0';
	filopt_.geoid[0] = '\0';
	filopt_.dcb[0] = '\0';
	filopt_.blq[0] = '\0';
	filopt_.solstat[0] = '\0';
	filopt_.trace[0] = '\0';
	for (i = 0; i < 2; i++) antpostype_[i] = 0;
	elmask_ = 15.0;
	elmaskar_ = 0.0;
	elmaskhold_ = 0.0;
	for (i = 0; i < 2; i++) for (j = 0; j < 3; j++) {
		antpos_[i][j] = 0.0;
	}
	exsats_[0] = '\0';
}