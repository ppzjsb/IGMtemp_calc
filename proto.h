
/* Function prototypes */

/* main.c */
void InitDumpMemory();
void write_data();

/* utils.c */
double dmax(double x,double y);
double dmin(double x,double y);

/* cooling.c */
void InitCool();
void ReadIonizeParams(char *fname);
void IonizeParamsTable();
void MakeCoolingTable();
void InitCoolMemory();
#ifdef DARK_PHOTON
double DoCooling(double u,double rho,double rho_old,double dt,double *nHp_guess,double *nHep_guess,double *nHepp_guess);
#else
double DoCooling(double u,double rho,double dt,double *nHp_guess,double *nHep_guess,double *nHepp_guess);
#endif
double solve_rates(double u, double rho, double dt, double *nHp_guess, double *nHep_guess, double *nHepp_guess);

/* UVB.c */
void IonizationRates(double logxe);
void InitUVB();
void InitUVBMemory();
void MakeUVBTable();
void ReadSecondaryTable();
void InitSecondaryMemory();

/* gausslegendre.c */
double plgndr(int l, int m, double x);
void gausslegendre(void);
void absc_and_weight(void);
void InitGLMemory(void);

/* elec_interp.c */
void MakeSecondaryTable();
