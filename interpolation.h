REAL interpolation(PARA_DATA *para, REAL *d0, REAL x_1, REAL y_1, REAL z_1,
                   int p, int q, int r);

REAL interpolation_bilinear(REAL x_1, REAL y_1, REAL z_1,
                            REAL d000, REAL d010, REAL d100, REAL d110,
                            REAL d001, REAL d011, REAL d101, REAL d111);
REAL interpolation_linear(REAL x_1, REAL y_1, REAL d000, REAL d010, REAL d100, REAL d110);

REAL advection_upwind(PARA_DATA *para, REAL **var,REAL *d0,int i, int j, int k);

void interpolation_coef(PARA_DATA *para,REAL **var, REAL *coef,REAL x1, REAL y1, REAL z1,int p, int q, int r);
 REAL interpolation_temp(PARA_DATA *para, REAL **var,REAL *d0, REAL x_1, REAL y_1, REAL z_1,
                   int i0, int j0, int k0, int p, int q, int r);

