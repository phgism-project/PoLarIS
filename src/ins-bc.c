/*
 *  Boundary map & boundary condition.
 *
 *
 *
 *
 *  */
#include "ins.h"
#include "parameters.h"
#include "mat_op3.h"

/*********************/
/* static parameters */
/*********************/
static FLOAT Time, Re, nu;
static FLOAT Re = 1.0;
static FLOAT nu = 1.0;
static FLOAT Time0;
FLOAT nu_max;
FLOAT nu_min;



/* 
 * Wrapper for func w.r.t time.
 * */
#if USE_MOC
# define FUNC_T_WRAPPER(func_xyz)					\
    void								\
    func_xyz##_t(FLOAT x, FLOAT y, FLOAT z, FLOAT t, FLOAT *values)	\
    {									\
	FLOAT t_save = Time;						\
	setFuncTime(t);							\
	func_xyz(x, y, z, values);					\
	setFuncTime(t_save);						\
    }

FUNC_T_WRAPPER(func_u)
FUNC_T_WRAPPER(func_gradu)
FUNC_T_WRAPPER(func_p)
FUNC_T_WRAPPER(func_f)
//FUNC_T_WRAPPER(func_wind)
#undef FUNC_T_WRAPPER
#endif	/* USE_MOC */


void setFlowParameter(FLOAT Re_in, FLOAT nu_in, FLOAT time_in)
{
    Re = Re_in;
    nu = nu_in;
    Time = time_in;
}

void setFuncTime(FLOAT time_in)
{
    Time = time_in;
}

void adjust_time(FLOAT delta_time)
{
    Time0 = Time; 
    Time += delta_time;
    phgInfo(2, "* Adjust static Time: %E\n", Time);
}

void restore_time()
{
    Time = Time0;
    phgInfo(2, "* Restore static Time: %E\n", Time);
}

void									
func_xyz_(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)	   
{
    *(values++) = x;
    *(values++) = y;
    *(values++) = z; 
}






/* ------------------------------------------------------------
 *    
 *    Effective viscosity
 *    
 * ------------------------------------------------------------ */

FLOAT 
get_effective_viscosity(const FLOAT *gu, FLOAT T, FLOAT p,
			int viscosity_type)
{
    FLOAT A = 1e-16;
    const FLOAT n = POWER_N;
    const FLOAT a = SEC_PER_YEAR;
    const FLOAT L = 100;	/* 100km */

    FLOAT yeta, eps = 0.;
    FLOAT eu[Dim*Dim];
    int k;

    Unused(a);

#  if 0
#  warning const vis -------------------------------------------------
    nu_max = nu_min = 1e7;
    return 1e7;
#  endif

    const FLOAT T0 = ARRHENIUS_T;
    const FLOAT Q0 = ARRHENIUS_Q0;
    const FLOAT a0 = ARRHENIUS_A0;
    const FLOAT Q1 = ARRHENIUS_Q1;
    const FLOAT a1 = ARRHENIUS_A1;
    const FLOAT R  = GAS_CONSTANT;

    if (ns_params->temp_viscosity == 2) {
	/* Use temperature in vis */
	if (T < T0)
	    A = a0 * exp(-Q0 / R / T);
	else
	    A = a1 * exp(-Q1 / R / T);
	A *= SEC_PER_YEAR;
    }
    else if (ns_params->temp_viscosity == 1) {
	T = TEMP_WATER - 10.;
	A = a1 * exp(-Q1 / R / T);
	A *= SEC_PER_YEAR;
    }
    else {
	A = ns_params->constant_A;
	//A = 3.1536e-18; 
    }
    

    if (viscosity_type == VIS_CONST) {
	eps = (100.) / (L * 1000.);	/* initial guess */
    } else if (viscosity_type == VIS_STRAIN) {
	if (ns_params->core_type == SIA) {
	    abort();
	} 
	else if (ns_params->core_type == FIRST_ORDER) {
	    FLOAT ux, uy, uz, vx, vy, vz;
	    
	    ux = gu[0] / LEN_SCALING; 
	    uy = gu[1] / LEN_SCALING; 
	    uz = gu[2] / LEN_SCALING;
	    vx = gu[3] / LEN_SCALING; 
	    vy = gu[4] / LEN_SCALING; 
	    vz = gu[5] / LEN_SCALING;

	    eps = sqrt(ux*ux + vy*vy + ux*vy
		       + 0.25 * pow(uy + vx, 2)
		       + 0.25 * (uz*uz + vz*vz));
	} 
	else if (ns_params->core_type == STOKES) {
	    /* Stokes */
	    MAT3_SYM(gu, eu);
	    for (k = 0; k < DDim; k++)
		eu[k] /= LEN_SCALING;
	    eps = sqrt(.5) * MAT3_NORM2(eu);
	} 
	else if (ns_params->core_type == DEBUG_CORE1) {
	    eps = (100.) / (L * 1000.);	/* const */
	}
	else {
	    phgError(1, "Unknown core type\n!!!");
	}
    } else {
	phgError(1, "Unknown viscosity type\n!!!");
    }

    if (eps < MIN_EFFECTIVE_STRAIN) 
	eps = MIN_EFFECTIVE_STRAIN;

    yeta = pow(A, -1./n) * pow(eps, (1.-n)/n);

    nu_max = MAX(nu_max, yeta);
    nu_min = MIN(nu_min, yeta);

    return yeta;
}



/* ------------------------------------------------------------
 *    
 *    Fraction coef
 *    
 * ------------------------------------------------------------ */


void 
func_beta(FLOAT x, FLOAT y, FLOAT z, FLOAT *beta) 
/* Beta square */
/* Input: real coord (km).
 * Output: beta2
 * */
{
    /*
     * Suppose the coordinate is already shifted when the topo is inited.
     *
     * */
    *beta = interp_topo_data('B', x, y);

    return;
}




/* ------------------------------------------------------------
 *    
 *    B.C. funcs
 *    
 * ------------------------------------------------------------ */




#if SIMPLE_TEST
/* ------------------------------------------------------------
 *
 * Simple test for Stokes eqn convergence test
 *
 * ------------------------------------------------------------ */

/* slop */
void func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *u)
{
    SIMPLE_INIT_SLOP;
    FLOAT z0 = z + x * tan_alpha;

#if 0    
    u[0] = 0;
    u[1] = z0;
    u[2] = 0;
#else
    u[0] = z0;
    u[1] = 0;
    u[2] = -z0 * tan_alpha;
#endif
    
}

void func_p(FLOAT x, FLOAT y, FLOAT z, FLOAT *p)
{
    SIMPLE_INIT_SLOP;
    FLOAT z0 = z + x * tan_alpha;

    //p[0] = 10;
    //p[0] = y;
    p[0] = 0.01 * z0;
}

void func_gradu(FLOAT x, FLOAT y, FLOAT z, FLOAT *gradu) {
    gradu[0] = 0;
    gradu[1] = 0;
    gradu[2] = 0;
    gradu[3] = 0;
    gradu[4] = 0;
    gradu[5] = 0;
    gradu[6] = 0;
    gradu[7] = 0;
    gradu[8] = 0;
}

void func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *f)
{
    SIMPLE_INIT_SLOP;
    FLOAT z0 = z + x * tan_alpha;
    
    /* f[0] = 0; */
    /* f[1] = 0; */
    /* f[2] = 0; */

    f[0] = 0.01 * tan_alpha;
    f[1] = 0.01 * 0;
    f[2] = 0.01 * 1;
}

void func_g1(FLOAT x, FLOAT y, FLOAT z, FLOAT *g1) {
    g1[0] = 0;
    g1[1] = 0;
    g1[2] = 0;
}

void func_g2(FLOAT x, FLOAT y, FLOAT z, FLOAT *g2) {
    g2[0] = 0;
    g2[1] = 0;
    g2[2] = 0;
}

void func_g3(FLOAT x, FLOAT y, FLOAT z, FLOAT *g3) {
    g3[0] = 0;
    g3[1] = 0;
    g3[2] = 0;
}

void func_T(FLOAT x, FLOAT y, FLOAT z, FLOAT t, FLOAT *T) 
{
    *T = interp_topo_dataT('T', x, y, t);
    return;
}

void func_fT(FLOAT x, FLOAT y, FLOAT z, FLOAT *fT) {
    fT[0] = 0;
}

void func_q(FLOAT x, FLOAT y, FLOAT z, FLOAT t, FLOAT *q) 
{
    *q = interp_topo_dataT('q', x, y, t);
    return;
}
#else
/* ------------------------------------------------------------
 *
 * Ice sheet dynamics
 *
 * ------------------------------------------------------------ */
void func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *u) {
    u[0] = 0;
    u[1] = 0;
    u[2] = 0;
}

void func_p(FLOAT x, FLOAT y, FLOAT z, FLOAT *p) {
    p[0] = 0;
}

void func_gradu(FLOAT x, FLOAT y, FLOAT z, FLOAT *gradu) {
    gradu[0] = 0;
    gradu[1] = 0;
    gradu[2] = 0;
    gradu[3] = 0;
    gradu[4] = 0;
    gradu[5] = 0;
    gradu[6] = 0;
    gradu[7] = 0;
    gradu[8] = 0;
}

void func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *f) {
    f[0] = 0;
    f[1] = 0;
    f[2] = 0;
}

void func_g1(FLOAT x, FLOAT y, FLOAT z, FLOAT *g1) {
    g1[0] = 0;
    g1[1] = 0;
    g1[2] = 0;
}

void func_g2(FLOAT x, FLOAT y, FLOAT z, FLOAT *g2) {
    g2[0] = 0;
    g2[1] = 0;
    g2[2] = 0;
}

void func_g3(FLOAT x, FLOAT y, FLOAT z, FLOAT *g3) {
    g3[0] = 0;
    g3[1] = 0;
    g3[2] = 0;
}

void func_T(FLOAT x, FLOAT y, FLOAT z, FLOAT *T) 
{
    *T = interp_topo_data('T', x, y);
    return;
}

void func_T_t(FLOAT x, FLOAT y, FLOAT z, FLOAT t, FLOAT *T) 
{
    *T = interp_topo_dataT('T', x, y, t);
    return;
}

void func_fT(FLOAT x, FLOAT y, FLOAT z, FLOAT *fT) {
    fT[0] = 0;
}

void func_q(FLOAT x, FLOAT y, FLOAT z, FLOAT *q) 
{
    *q = interp_topo_data('q', x, y);
    return;
}

void func_q_t(FLOAT x, FLOAT y, FLOAT z, FLOAT t, FLOAT *q) 
{
    *q = interp_topo_dataT('q', x, y, t);
    return;
}
#endif


/* ------------------------------------------------------------
 *    
 *    b.C. map
 *    
 * ------------------------------------------------------------ */

int
my_bc_map(int bctype)
{
    switch (bctype) {
    case 3:
	return BC_LATERAL;
    case 4:
	return BC_DIVIDE;
    case 5:
	return BC_FRONT;
    case 13:
	return BC_TOP;
    case 17:
	return SLIP_BDRY;
    default:
	return DIRICHLET;
    }
}



void 
set_boundary_mask(NSSolver *ns)    
{
    DOF **u = ns->u;
    DOF **p = ns->p;
    DOF **T = ns->T;

    BTYPE DB_masks[3] = {0, 0, 0};


    if (ns_params->use_slide) {
	/*
	 * first direction: bottom out normal
	 * second direcion: lateral out normal
	 * */
	DB_masks[UN_DIR] = SLIP_BDRY;
	DB_masks[LN_DIR] = BC_LATERAL;

	/* divide is fixed */
	DB_masks[X_DIR] |= BC_DIVIDE;
	DB_masks[Y_DIR] |= BC_DIVIDE;
	DB_masks[Z_DIR] |= BC_DIVIDE;
    }
    else {
	DB_masks[0] = SLIP_BDRY | BC_DIVIDE; 
	DB_masks[1] = SLIP_BDRY | BC_DIVIDE;
	DB_masks[2] = SLIP_BDRY | BC_DIVIDE;

	DB_masks[LN_DIR] |= BC_LATERAL;
    }

#if SIMPLE_TEST
    DB_masks[0] = SLIP_BDRY | BC_TOP;
    DB_masks[1] = SLIP_BDRY | BC_TOP;
    DB_masks[2] = SLIP_BDRY | BC_TOP;
#endif
    

    phgDofSetDirichletBoundaryMasks(u[1], DB_masks);
    phgDofSetDirichletBoundaryMask(p[1], BC_PIN_NODE);
    phgDofSetDirichletBoundaryMask(T[1], BC_TOP);

    return;
}


/* Unused */
void
load_case_options()
{

}

void iceSetUserFunc(NSSolver *ns)
{
    /* Do nothing */
}

void func_a(FLOAT x, FLOAT y, FLOAT *q)
{
}

void func_s(FLOAT x, FLOAT y, FLOAT z, FLOAT *q)
{
}

void func_b(FLOAT x, FLOAT y, FLOAT z, FLOAT *q)
{
}

