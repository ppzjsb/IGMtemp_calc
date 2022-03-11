
/* Function prototypes for the entire code */

/* main.c */
void InitDumpMemory();
void write_data();

/* utils.c */
void endrun(int ierr);
double dmax(double x,double y);
double dmin(double x,double y);

/* cooling.c */
void InitCool(double temp, double logxe);
void ReadIonizeParams(char *fname);
void IonizeParamsTable();
void MakeCoolingTable();
void InitCoolMemory();
double DoCooling(double u,double rho,double dt,double *nHp_guess,double *nHep_guess,double *nHepp_guess);
double solve_rates(double u, double rho, double dt, double *nHp_guess, double *nHep_guess, double *nHepp_guess);

/* UVB.c */
void IonizationRates(double temp, double logxe);
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
