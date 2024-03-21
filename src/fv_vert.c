/*
 * 
 * Vertex-based Finite volume solver for conservation law.
 *
 * By Wei Leng, Lsec, 2012.
 *
 * Input:
 *   1) nodes coordinates
 *   2) triangle to nodes
 *
 * Output:
 *   error
 *
 * TODO:
 *   1) no boundary condition is implemented
 *
 * 
 *  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <netcdf.h>
#include "phg/netcdf-utils.h"
#include "mpi.h"
#include "ins.h"

static int dbg = 0;


#define DUAL_GRID 0			/* 0: mass,
					 * 1: circum */


void phgInfo(int verbose_level, const char *fmt, ...);
void phgError(int code, const char *fmt, ...);


/* ------------------------------------------------------------
 * 
 *  Finite volume solver Interface
 *
 * ------------------------------------------------------------ */

typedef void (*DOF_USER_FUNC_T)(double x, double y, double z, 
				double t, double *values);

void fv_solver_init(void **fv_data,
		    const char *mesh_file,
		    double (*verts)[4],
		    const DOF_USER_FUNC_T func_f);
void fv_update(void *fv_data,
	       const double *H,
	       const FLOAT (*verts)[4],
	       const double *U, 
	       double *dH, 
	       double *U_vert, 
	       double t);




typedef struct FV_ELEM_ {
    int v[3];
    //double x[3], y[3];
    //double xe[3], ye[3];
    double xc, yc;
    double len[3];
    double n[3][2];		/* normal */
    double a[3][2];		/* area */
    double quad_pt[3][2][2];	/* quad point
				 * [nedge][v0|1][X|Y] */
} FV_ELEM; 


typedef struct GRID_2D_ {
    int nv;
    int ne;
    FV_ELEM *tri;		/* [ne] */
    double *X;			/* [nv] */
    double *Y;
    double *node_area;		/* [nv] */
    int *bdry_mark;		/* [nv] */
    DOF_USER_FUNC_T func_f;	/* source function */
} GRID_2D; 


//static GRID_2D g2D;

static int edge2vert[3][2] = {
    {0, 1}, {1, 2}, {2, 0}
};

/* static double *H; */
/* static double *dH; */
/* static double *H_old; */


/* boundary condition */
static double
func_h(double x, double y, double t)
{
    return 0.;
}

static double
get_tri_area(double x0, double x1, double x2, 
	     double y0, double y1, double y2, 
	     double *xc) 
{
    double a0 = x1 - x0;
    double a1 = x2 - x1;
    double a2 = x0 - x2;

    double b0 = y1 - y0;
    double b1 = y2 - y1;
    double b2 = y0 - y2;

    double l0 = sqrt(a0*a0 + b0*b0);
    double l1 = sqrt(a1*a1 + b1*b1);
    double l2 = sqrt(a2*a2 + b2*b2);

    double s = (l0+l1+l2)/2; 
    double area = sqrt(fabs(s*(s-l0)*(s-l1)*(s-l2)));

    if (xc != NULL) {
	xc[0]= (x0 + x1 + x2) / 3.;
	xc[1]= (y0 + y1 + y2) / 3.;
    }

    return area;
}


static void
set_solution_on_bdry(GRID_2D *g2D, double *solu, double time, int set_bdry_only)
/*
 *  Set solution on bdry,
 *     set_bdry_only: false, all
 *                    true , only bdry
 *
 *
 * */
{
    FV_ELEM *fv_tri = g2D->tri;
    int nv = g2D->nv;
    int ne = g2D->ne;
    double *X = g2D->X, *Y = g2D->Y;
    int i, j, k; 


#if 0
#   warning h bdry proj 
    if (set_bdry_only) {
	for (i = 0; i < nv; i++)
	    if (g2D->bdry_mark[i] == -1)
		solu[i] = 0;	/* clear bdry */
    } else {
	bzero(solu, nv * sizeof(double));
    }


    for (i = 0; i < ne; i++) {
	FV_ELEM *tri = &fv_tri[i];

	int V[3] = {tri.v[0], tri.v[1], tri.v[2]};
	double x[3] = {X[V[0]], X[V[1]], X[V[2]]};
	double y[3] = {Y[V[0]], Y[V[1]], Y[V[2]]};

	if (set_bdry_only && 
	    !(g2D->bdry_mark[V[0]] == -1 ||
	      g2D->bdry_mark[V[1]] == -1 ||
	      g2D->bdry_mark[V[2]] == -1))
	    continue;

	for (k = 0; k < 3; k++) {
	    int v0 = edge2vert[k][0];
	    int v1 = edge2vert[k][1];
	    int V0 = tri.v[v0];
	    int V1 = tri.v[v1];

	    double xe = (x[v0] + x[v1]) / 2.;
	    double ye = (y[v0] + y[v1]) / 2.;

	    double xc = tri.xc;
	    double yc = tri.yc;

	    double a0 = tri.a[k][0];
	    double a1 = tri.a[k][1];

	    // x[v0], y[v0];
	    // x[v1], y[v1];

	    double xx, yy;
	    if (!set_bdry_only || g2D->bdry_mark[V0] == -1) {
		xx = (xc + xe + x[v0]) / 3;
		yy = (yc + ye + y[v0]) / 3;
		solu[V0] += a0 * func_h(xx, yy, time);
	    }

	    if (!set_bdry_only || g2D->bdry_mark[V1] == -1) {
		xx = (xc + xe + x[v1]) / 3;
		yy = (yc + ye + y[v1]) / 3;
		solu[V1] += a1 * func_h(xx, yy, time);
	    }
	}
    }

    for (i = 0; i < nv; i++) {
	if (!set_bdry_only || g2D->bdry_mark[i] == -1)
	    solu[i] /= g2D->node_area[i];
    }
    
#else
#   warning h bdry interp pointwise
    for (i = 0; i < nv; i++) {
	if (!set_bdry_only || g2D->bdry_mark[i] == -1)
	    solu[i] = func_h(X[i], Y[i], time);
    }
#endif
}



static void
set_periodic_bdry(GRID_2D *g2D, double *dH)
{
    int i, id;
    int nv = g2D->nv;
    
    /* First, collect to master */
    for (i = 0; i < nv; i++) {
	if (g2D->bdry_mark[i] == -1
	    || g2D->bdry_mark[i] == i)
	    continue;

	id = g2D->bdry_mark[i];
	assert(g2D->bdry_mark[id] == id);
	dH[id] += dH[i];
    }

    /* Second, distribue to slaves */
    for (i = 0; i < nv; i++) {
	if (g2D->bdry_mark[i] == -1
	    || g2D->bdry_mark[i] == i)
	    continue;

	id = g2D->bdry_mark[i];
	dH[i] = dH[id];
    }
}








//#define SQUARE(x) (x)*(x)

static double
det2x2_by_formula(const double a00,  const double a01,
          const double a10,  const double a11)
{
    // First compute the det2x2
    const double m01 = a00*a11 - a10*a01;
    return m01;
}


static double
det3x3_by_formula(const double a00,  const double a01,  const double a02,
          const double a10,  const double a11,  const double a12,
          const double a20,  const double a21,  const double a22)
{
    // First compute the det2x2
    const double m01 = a00*a11 - a10*a01;
    const double m02 = a00*a21 - a20*a01;
    const double m12 = a10*a21 - a20*a11;

    // Now compute the minors of rank 3
    const double m012 = m01*a22 - m02*a12 + m12*a02;
    return m012;
}

static double
get_dir_area(const double *p, const double *q, const double *s)
{
    double psx = p[0]-s[0];
    double psy = p[1]-s[1];
    double qsx = q[0]-s[0];
    double qsy = q[1]-s[1];

    return psx * qsy - psy * qsx;
}


static void 
get_circumcenter(double *x, double *y, double *pxc, double *pyc) 
{
    double p[3] = {x[0], y[0], 0.};
    double q[3] = {x[1], y[1], 0.};
    double s[3] = {x[2], y[2], 0.};

    double psx = p[0]-s[0];
    double psy = p[1]-s[1];
    double psz = p[2]-s[2];
    double ps2 = SQUARE(psx) + SQUARE(psy) + SQUARE(psz);
    double qsx = q[0]-s[0];
    double qsy = q[1]-s[1];
    double qsz = q[2]-s[2];
    double qs2 = SQUARE(qsx) + SQUARE(qsy) + SQUARE(qsz);
    double rsx = psy*qsz-psz*qsy;
    double rsy = psz*qsx-psx*qsz;
    double rsz = psx*qsy-psy*qsx;

    double num_x = ps2 * det2x2_by_formula(qsy,qsz,rsy,rsz)
        - qs2 * det2x2_by_formula(psy,psz,rsy,rsz);
    double num_y = ps2 * det2x2_by_formula(qsx,qsz,rsx,rsz)
        - qs2 * det2x2_by_formula(psx,psz,rsx,rsz);
    double num_z = ps2 * det2x2_by_formula(qsx,qsy,rsx,rsy)
        - qs2 * det2x2_by_formula(psx,psy,rsx,rsy);

    double den   = det3x3_by_formula(psx,psy,psz,
                    qsx,qsy,qsz,
                    rsx,rsy,rsz);

    assert( fabs(den) > 1e-10 );
    double inv = 1 / (2 * den);

    double xc = s[0] + num_x*inv;
    double yc = s[1] - num_y*inv;
    double zc = s[2] + num_z*inv;

    /* 2D check */
    assert(fabs(zc) < 1e-15);
    double a[3], a0, X[2] = {xc, yc};
    a0 = get_dir_area(p, q, s);
    a[0] = get_dir_area(p, q, X);
    a[1] = get_dir_area(q, s, X);
    a[2] = get_dir_area(s, p, X);

    if (a0 * a[0] >= 0
	&& a0 * a[1] >= 0
	&& a0 * a[2] >= 0) {			
	/* inside */
	*pxc = xc;
	*pyc = yc;
    } else {
	double len_edge[3], max_len = -1;
	double xe[3],ye[3];
	int k, ilong = -1;

	for (k = 0; k < 3; k++) {
	    int v0 = edge2vert[k][0];
	    int v1 = edge2vert[k][1];

	    double vec[2] = {x[v0]-x[v1], 
			    y[v0]-y[v1]};
	    len_edge[k] = sqrt(vec[0]*vec[0] + vec[1]*vec[1]);
	    xe[k] = (x[v0] + x[v1]) / 2.;
	    ye[k] = (y[v0] + y[v1]) / 2.;

	    if (len_edge[k] > max_len) {
		max_len = len_edge[k];
		ilong = k;
	    }
	}

	*pxc = xe[ilong];
	*pyc = ye[ilong];
    }

    return;
}



void
fv_solver_init(void **fv_data_ptr,
	       const char *mesh_file,
	       double (*verts)[4],
	       DOF_USER_FUNC_T func_f)
{
    FILE *fp, *fp1;
    int i, j, k, id;
    size_t nv, ne;
    double *X, *Y;
    int nbdry = 0;
    int rank;
    GRID_2D *g2D;    
    int ncid, dimid, varid, dim2;
    int retval;

    g2D = calloc(1, sizeof(*g2D));
    *fv_data_ptr = (void *) g2D;
    g2D->func_f = func_f;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    /* Open the file. */
    if ((retval = nc_open(mesh_file, NC_NOWRITE, &ncid)))
	NC_ERR(retval);
    

    /* nodes info */
    NC_READ_DIM("nvert2d", nv);
    NC_CHECK_DIM_1D("verts2d_info", nv);
    int *vinfo = calloc(nv, sizeof(float));
    NC_READ_INT("verts2d_info", vinfo);

    X = calloc(nv, sizeof(*X));
    Y = calloc(nv, sizeof(*Y));
    g2D->nv = nv;
    g2D->X = X;
    g2D->Y = Y;
    g2D->bdry_mark = calloc(nv, sizeof(*g2D->bdry_mark));
    for (i = 0; i < nv; i++) {
#  warning Fixme: vert 2d mark 
	int mark = 0;		/* vinfo[i] */

	/* Use mapped coordinates */
	X[i] = verts[i][0];
	Y[i] = verts[i][1];
	
	/* mark: 
	 *     negtive:  bdry
	 *     positive: periodic id
	 * */
	if (mark > 0) 
	    g2D->bdry_mark[i] = mark - 1; /* periodic id */
	else if (mark < 0)
	    g2D->bdry_mark[i] = -1;       /* bdry */
	else 
	    g2D->bdry_mark[i] = i;        /* self */

	phgInfo(3, "g2D bdry mark: %5d -> %5d\n",
		mark, g2D->bdry_mark[i]);
		
	if (g2D->bdry_mark[i] == -1)
	    nbdry ++;
    }
    free(vinfo);


    /* triangles */
    NC_READ_DIM("nelem2d", ne);
    dim2 = 3;			/* 3 verts of one triangle */
    int *elem_data = calloc(ne * dim2, sizeof(int));
    NC_READ_INT("elems2d", elem_data);
    NC_CHECK_DIM_2D("elems2d", ne, dim2);

    FV_ELEM *fv_tri = (FV_ELEM *) calloc(ne, sizeof(*fv_tri));
    g2D->ne = ne;
    g2D->tri = fv_tri;
    for (i = 0; i < ne; i++) {
	int v0, v1, v2;

	v0 = elem_data[i * dim2 + 0];
	v1 = elem_data[i * dim2 + 1];
	v2 = elem_data[i * dim2 + 2];
	/* fscanf(fp, "%d %d %d %d",  */
	/*        &v0, &v1, &v2, &id); */

	fv_tri[i].v[0] = v0;
	fv_tri[i].v[1] = v1;
	fv_tri[i].v[2] = v2;
    }
    free(elem_data);

    if (rank == 0)
	printf("   Triangle mesh nv:%d ne:%d, nb:%d\n", 
	       (int) nv, (int) ne, nbdry);



    /*
     *
     * Bulid FV solver
     *
     * */
    double *node_area = calloc(nv, sizeof(double));
    g2D->node_area = node_area;
    double min_h = 1e20, max_h = -1e20;

	if (FALSE
		&& rank == 0) 
		fp = fopen("fv.dat", "w");
    else 
	fp = NULL;

    for (i = 0; i < ne; i++) {
	FV_ELEM *tri = &fv_tri[i];
	
	int V[3] = {tri->v[0], tri->v[1], tri->v[2]};
	assert(V[0] < nv);
	assert(V[1] < nv);
	assert(V[2] < nv);

	/* printf("elem %5d %5d %5d %5d\n", i, */
	/*        V[0], V[1], V[2]); */

	/* center */
	double x[3] = {X[V[0]], X[V[1]], X[V[2]]};
	double y[3] = {Y[V[0]], Y[V[1]], Y[V[2]]};

#if DUAL_GRID == 0
	double xc = (x[0] + x[1] + x[2]) / 3.;
	double yc = (y[0] + y[1] + y[2]) / 3.;
#else
	double xc, yc;
	get_circumcenter(x, y, &xc, &yc);
#endif
	if (fp) {
#define DRAW_LINE(n0, n1)				\
	    fprintf(fp, "%e %e\n", x[n0], y[n0]);	\
	    fprintf(fp, "%e %e\n\n\n", x[n1], y[n1]);
	    
	    DRAW_LINE(0, 1);
	    DRAW_LINE(1, 2);
	    DRAW_LINE(2, 0);
	}
	

	tri->xc = xc;
	tri->yc = yc;
	
	{
	    double a0 = x[1] - x[0];
	    double a1 = x[2] - x[1];
	    double a2 = x[0] - x[2];

	    double b0 = y[1] - y[0];
	    double b1 = y[2] - y[1];
	    double b2 = y[0] - y[2];

	    double l0 = sqrt(a0*a0 + b0*b0);
	    double l1 = sqrt(a1*a1 + b1*b1);
	    double l2 = sqrt(a2*a2 + b2*b2);
	    double l[3] = {l0, l1, l2};

	    for (j = 0; j < 3; j++) {
		if (max_h < l[j])
		    max_h = l[j];
		if (min_h > l[j])
		    min_h = l[j];
	    }
	}

	/* for edges */
	for (k = 0; k < 3; k++) {
	    int v0 = edge2vert[k][0];
	    int v1 = edge2vert[k][1];
	    int V0 = tri->v[v0];
	    int V1 = tri->v[v1];

	    double xe = (x[v0] + x[v1]) / 2.;
	    double ye = (y[v0] + y[v1]) / 2.;

	    double xx = xe - xc;
	    double yy = ye - yc;

	    double len = sqrt(xx*xx + yy*yy);

#if DUAL_GRID == 0
	    double n[2] = {yy, -xx};
	    n[0] /= len; n[1] /= len;
#else
	    double n[2] = {x[v0] - x[v1], 
			   y[v0] - y[v1]};
	    double nn = sqrt(n[0]*n[0] + n[1]*n[1]);
	    n[0] /= nn; n[1] /= nn;
#endif

	    double t[2] = {x[v1] - x[v0], y[v1] - y[v0]};
	    if (n[0]*t[0] + n[1]*t[1] < 0) {
		n[0] *= -1; n[1] *= -1;
	    }
	    
	    double a0 = get_tri_area(xc, x[v0], xe,
				     yc, y[v0], ye,
				     tri->quad_pt[k][0]);
	    double a1 = get_tri_area(xc, x[v1], xe,
				     yc, y[v1], ye, 
				     tri->quad_pt[k][1]);

	    tri->len[k] = len;
	    tri->n[k][0] = n[0];
	    tri->n[k][1] = n[1];
	    tri->a[k][0] = a0;
	    tri->a[k][1] = a1;

	    node_area[V0] += a0;
	    node_area[V1] += a1;
	}
    }
    set_periodic_bdry(g2D, node_area);
    
    if (fp)
	fclose(fp);

    if (0 && rank == 0) {
	fp = fopen("fv1.dat", "w");
	for (i = 0; i < ne; i++) {
	    FV_ELEM *tri = &fv_tri[i];
	    fprintf(fp, "%e %e\n", tri->xc, tri->yc);
	}
	fclose(fp);
    }

    if (rank == 0)
	printf("   Triangle mesh h: [%e %e]\n", min_h, max_h);

    if ((retval = nc_close(ncid)))
	NC_ERR(retval);

    return;
}



/*
 * FV solver update
 * Input: U [ne][3][2], 
 *        H [nv]
 * Output: dH [nv]
 * */
void
fv_update(void *g2D_,
	  const double *H, 	/* km */
	  const FLOAT (*verts)[4],  /* km */
	  const double *U, 	/* km/a */
	  double *dH, 		/* km */
	  double *U_vert,	/* km/a */
	  double t)		/* a */
{
    GRID_2D *g2D = (GRID_2D *) g2D_; /* cast */
    FV_ELEM *fv_tri = g2D->tri;
    int nv = g2D->nv;
    int ne = g2D->ne;
    double *X = g2D->X, *Y = g2D->Y;
    double *node_area = g2D->node_area;
    int i, j, k; 
    DOF_USER_FUNC_T func_f = g2D->func_f;
    FILE *fp;
	BOOLEAN output = ns_params->output_fv_vert;

    double *H_old = calloc(nv, sizeof(double));
    memcpy(H_old, H, nv * sizeof(double));
    bzero(dH, nv * sizeof(double));

    if (U_vert != NULL) {
	bzero(U_vert, nv * sizeof(double));
    }

	double *dH0_vert = NULL, *a_vert = NULL;

	if (output) {
		dH0_vert = calloc(2 * nv, sizeof(double)); /* dH caused by velocity only */
		a_vert = dH0_vert + nv;						/* accu on vert */
	}

    for (i = 0; i < ne; i++) {
	FV_ELEM *tri = &fv_tri[i];
	
	if (dbg) printf("\n\nelem: %d\n", i);
	int V[3] = {tri->v[0], tri->v[1], tri->v[2]};

	/* Compute flux */
	for (k = 0; k < 3; k++) {
	    if (dbg) printf("edge: %d\n", k);

	    int v0 = edge2vert[k][0];
	    int v1 = edge2vert[k][1];
	    int V0 = tri->v[v0];
	    int V1 = tri->v[v1];

	    double *a = tri->a[k];
	    double *n = tri->n[k];
	    double len = tri->len[k];
	    double h = 0;
	    double *qpt = &tri->quad_pt[k][0][0];

	    double vel_u = U[i*6 + 2*k    ];
	    double vel_v = U[i*6 + 2*k + 1];

#if 0
#   warning Flux central
	    /* Central */
	    h = .5 * (H_old[V0] + H_old[V1]);
#else
#   warning Flux upwind
	    /*
	     * 
	     * Hpwind: v0 -> v1
	     * 
	     * */
	    if (n[0] * vel_u + n[1] * vel_v > 0) {
		h = H_old[V0];
	    } else {
		h = H_old[V1];
	    }
#endif

	    double flux = len * (n[0] * vel_u + n[1] * vel_v) * h;

	    dH[V0] -= flux / node_area[V0]; /* unit: km */
	    dH[V1] += flux / node_area[V1];

	    if (output) {
		dH0_vert[V0] -= flux / node_area[V0]; /* unit: km */
		dH0_vert[V1] += flux / node_area[V1];
	    }

	    /* source term */
	    if (func_f != NULL) {
		double f0, f1;
		func_f(qpt[0], qpt[1], verts[V0][3], t, &f0); /* unit: m */
		func_f(qpt[2], qpt[3], verts[V1][3], t, &f1);

		f0 /= 1e3; /* unit: km */
		f1 /= 1e3;

		dH[V0] += f0 * a[0] / node_area[V0];
		dH[V1] += f1 * a[1] / node_area[V1];

		if (output) {
		    a_vert[V0] += f0 * a[0] / node_area[V0];
		    a_vert[V1] += f1 * a[1] / node_area[V1];
		}
	    }
	    

	    /* velocity at vert */
	    if (U_vert != NULL) {
		U_vert[2*V0   ] +=  vel_u *  a[0] / node_area[V0];
		U_vert[2*V0 +1] +=  vel_v *  a[0] / node_area[V0];
		U_vert[2*V1   ] +=  vel_u *  a[1] / node_area[V1];
		U_vert[2*V1 +1] +=  vel_v *  a[1] / node_area[V1];
		/* printf("vert %4d vel %15.7e %15.7e, area: %15.7e %15.7e\n",  */
		/*        V0, vel_u, vel_v, a[0] / node_area[V0], U_vert[2*V0]); */
		/* printf("vert %4d vel %15.7e %15.7e, area: %15.7e %15.7e\n",  */
		/*        V1, vel_u, vel_v, a[1] / node_area[V1], U_vert[2*V1]); */
	    }

	    if (dbg) printf("%20.10lf %20.10lf, [%20.10lf %20.10lf] [%20.10lf %20.10lf]\n",
			    len, h, n[0], n[1], dH[V0], dH[V1]);
	} /* end flux */
    }     /* end elem */

    /* Set boundary */
    set_solution_on_bdry(g2D, dH, 0., 1);

    
    /* Periodic correct */
    set_periodic_bdry(g2D, dH);


    if (output) {
	fp = fopen("output/fv-vert.vtk", "w");
	fprintf(fp, "# vtk DataFile Version 2.0\n"
		"Triangle grid created by WaveHdiv.\n"
		//"BINARY\nDATASET UNSTRUCTURED_GRID\n"
		"ASCII\nDATASET UNSTRUCTURED_GRID\n");

	fprintf(fp, "\nPOINTS %d double\n", g2D->nv);
	for (i = 0; i < g2D->nv; i++) {
	    fprintf(fp, "%f %f %f\n",
		    g2D->X[i], g2D->Y[i], 0.);
	}
	fprintf(fp, "\nCELLS %d %d\n", g2D->ne, g2D->ne * (3 + 1));
	for (i = 0; i < g2D->ne; i++) {
	    fprintf(fp, "3 %d %d %d\n",
		    g2D->tri[i].v[0],
		    g2D->tri[i].v[1],
		    g2D->tri[i].v[2]);
	}
	fprintf(fp, "\nCELL_TYPES %d\n", g2D->ne);
	for (i = 0; i < g2D->ne; i++) {
	    fprintf(fp, "%d\n", 5);	 // VTK_TRIANGLE
	}
	fprintf(fp, "\nPOINT_DATA %d\n", g2D->nv);

	fprintf(fp, "\nVECTORS %s double\n", "Vel_UV");
	for (i = 0; i < g2D->nv; i++) {
	    fprintf(fp, "%f %f %f\n",
		    U_vert[2*i + 0], 
		    U_vert[2*i + 1], 
		    0.);
	}

	fprintf(fp, "\nSCALARS %s double %d\n", "dH_w", 1);
	fprintf(fp, "LOOKUP_TABLE default\n");
	for (i = 0; i < g2D->nv; i++) {
	    fprintf(fp, "%f\n", dH0_vert[i]);
	}

	fprintf(fp, "\nSCALARS %s double %d\n", "a", 1);
	fprintf(fp, "LOOKUP_TABLE default\n");
	for (i = 0; i < g2D->nv; i++) {
	    fprintf(fp, "%f\n", a_vert[i]);
	}
		
	fprintf(fp, "\nSCALARS %s double %d\n", "dH_final", 1);
	fprintf(fp, "LOOKUP_TABLE default\n");
	for (i = 0; i < g2D->nv; i++) {
	    fprintf(fp, "%f\n", dH[i]);
	}

	fclose(fp);

	free(dH0_vert);
    }
    
    
#if 0
    fp = fopen("H_.m", "w");
    fprintf(fp, "H=[\n");
    for (i = 0; i < nv; i++)
	fprintf(fp, "%lf\n", H[i]);
    fprintf(fp, "];\n");
    fclose(fp);
#endif

    free(H_old);
    return;
}
