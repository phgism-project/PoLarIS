#ifndef INS_H
#include "phg.h"
#include <stdlib.h>
#include <string.h>
#include <strings.h>	/* bzero() */
#include <math.h>
#include <stdarg.h>

//#include "oct-search.h"

#if USE_MG 
#  include "multi-grid.h"
#endif /* USE_MG */
#include "layers.h"

/*
 * ================================================================================
 * 
 *      Control macro
 * 
 * ================================================================================
 * */

/*
 * Note:
 * 1. face quad for curved element is NOT implemented.
 *
 *  */

/* Control Macros */
#  define STEADY_STATE 0	  /* Steaty or time-dependent */
#  define TIME_DEP_NON 1          /* Time-dependent case, using non-linear solver  */


/* Grid */
#define USE_ISOP 0		/* Iso-parametric elements */

/* Stokes pin node */
#define PIN_AT_ROOT 0			/* Pin node options */

/* Linear solver */
#define USE_QP_ONLY 0			/* Use Qp in PCD precondtioner */
#define REUSE_MAT 0                     /* Reuse matrix, if matrix is fixed for each time step */
#define MAT_HANDLE_BDRY_EQNS FALSE
#ifndef DUMP_MAT_VEC			/* Dump mat vec for debug */
#   define DUMP_MAT_VEC 0
#endif


/* Sliding B.C. */
#define USE_SLIDING_BC  1	        /* Use sliding boundary condition */
#define ZERO_TRACT_ZONE 0
#define BC_BOTTOM (SLIP_BDRY | BC_FLOAT)


/* Temp solver */
#define USE_TEMP_SDPG 1
#define USE_TEMP_TIME 1
#define USE_TEMP_CONV 1
#define USE_TEMP_DIFF 1
#define USE_TEMP_HEAT 1
#define USE_TEMP_GEOFLUX 1


/* Mulitgrid solver */
#define USE_MG 0
#define MG_SOLVER_TYPE_U  0
#define MG_SOLVER_TYPE_P  1

#define MG_SOLVER_F  0
#define MG_SOLVER_Ap 1
#define MG_SOLVER_Qp 2

#define PC_MAX_LEVEL 10

#define NS_PROBLEM "ice"


/* Scaling */
#  define EQU_SCALING 1e-8
#  define LEN_SCALING 1e3
#  define PRES_SCALING 1e5
#  define EQU_T_SCALING 1e0
#  define LEN_SCALING2 (LEN_SCALING * LEN_SCALING)
#  define LEN_SCALING3 (LEN_SCALING * LEN_SCALING * LEN_SCALING)

/* Const parameters */
#include "parameters.h"
#include "phg/netcdf-utils.h"


/*
 * ================================================================================
 * 
 *          discretization options
 * 
 * ================================================================================
 * */

//#warning ----- FIXME on DOF_G1 -------
#define DOF_G1 DOF_P8

/* Bdry type */
#define BC_TOP       (BDRY_USER1)
#define SLIP_BDRY    (BDRY_USER2)
#define BC_FLOAT     (BDRY_USER3)
#define BC_LATERAL   (BDRY_USER4)
#define BC_DIVIDE    (BDRY_USER5)
#define BC_FRONT     (BDRY_USER6)
#define BC_PIN_NODE  (BDRY_USER7)
#define BC_INACTIVE  (BDRY_USER8)

#define SETFLOW      (BC_BOTTOM | BC_LATERAL)

/* #define INFLOW  BC_BOTTOM */
/* #define OUTFLOW BC_TOP */
//#define PINNEDNODE   (BDRY_USER6)

typedef enum { PICARD, NEWTON } LTYPE;	/* linearization type */
typedef enum { STOKES, SIA, SSA, FIRST_ORDER, DEBUG_CORE1, DEBUG_CORE2} CORE_TYPE;	/* linearization type */

/*
 * ================================================================================
 * 
 *          Sturcture
 * 
 * ================================================================================
 * */

typedef struct NSSolver_ NSSolver;


/* Parameters */
typedef struct NSParams_ {
    int core_type;	 /* Ice sheet core type */

    DOF_TYPE *utype;             /* DOF type for velocity */
    DOF_TYPE *ptype;		  /* DOF type for pressure */
    DOF_TYPE *T_type;           /* DOF type for Temperature */
    char *utype_name;		  /* DOF type name for velocity */
    char *ptype_name;		  /* DOF type name for pressure */
    char *T_type_name;	          /* DOF type name for coordinates */

    DOF_TYPE *mapping_type;           /* DOF type for Isop mapping */
    char *isop_type_name;
    //char *test_name;

    DOF_TYPE *utype_prv;             /* DOF type for velocity */
    DOF_TYPE *ptype_prv;		  /* DOF type for pressure */
    DOF_TYPE *T_type_prv;           /* DOF type for Temperature */
    char *utype_prv_name;		  /* DOF type name for velocity */
    char *ptype_prv_name;		  /* DOF type name for pressure */
    char *T_type_prv_name;	          /* DOF type name for coordinates */
    
    BOOLEAN reduce_mesh;	  /* Reduce mesh where h = 0 */
    BOOLEAN layered_mesh;	  /* Use layered mesh */
    BOOLEAN struct_mesh;	  /* Use struct mesh */
    BOOLEAN solve_temp;	          /* Solve temperature */
    BOOLEAN solve_height;	  /* Solve height */
    BOOLEAN update_bdry_type;	  /* Update boundary btype by 2D mesh */
    BOOLEAN resume;	          /* Resume step from previous step */
    BOOLEAN record;	          /* Save step data to resume */
    BOOLEAN enclosed_flow;        /* Enclosed flow */
    BOOLEAN pin_node;		  /* Pin a node */
    BOOLEAN curved_bdry;	  /* Use curved boundary */
    BOOLEAN curved_refine;	  /* Use curved refinement */
    BOOLEAN start_const_vis;	  /* Viscosity start at constant */
    BOOLEAN compensate_equ;	  /* Compensate euqations */
    
    //int Nx, Ny, Nz;		  /* Structured mesh size */

    INT periodicity;		  /* periodicity */
    INT period_repeat;		  /* periodical repeat */
    FLOAT Re;			  /* Reynolds num */
    FLOAT nu;			  /* Viscosity */
    FLOAT time_start;		  /* Time start */
    FLOAT time_end;		  /* Time end */
    FLOAT dt0;			  /* Inintial dt */
    FLOAT eps_diagP;		  /* Diagnal entries for pressure matrix */
    FLOAT fp_scale;		  /* Dirichlet scale for Fp */

    /* Discrete scheme */
    INT height_scheme;		  /* Height solver scheme */
    FLOAT Theta;		  /* 0.5: Crank-Nicolson
				   * 1. : Backward Euler */
    BOOLEAN use_PCD;		  /* Use Pressure Convetion-Diffusion preconditioner
				   *   default: True */
    int init_temp_type;		  /* Init temperature field
				   * 0: diff, 1: interp, 2: read data */
    BOOLEAN use_prism_elem;	  /* Use prism elems */

    int stab_scheme;		  /* Pressure stab scheme
				   * -1: none,
				   * 0: Proj,
				   * 1: Grad, P=\int x_i x_j
				   * 2: Grad, P=(dhx^2, dhz^2)
				   * */
    BOOLEAN use_slide;		
    int sliding_bdry_scheme; 	  /* 0: Dirich
				   * 1: Lagrangian Multiplier
				   * 2: Penalty
				   *  */
    FLOAT sliding_penalty;

    FLOAT stab_nu;		  /* Pressure stab, negtive means nonlinear */
    FLOAT stab_alpha;		  /* Pressure stab alpha */
    int stab_remove_static;	  /* Remove Pressure stab hydro static pressure
				   *   0: not removed, defult 
				   *   1: solve p, with ((p - ps) - \Pi (p - p))
				   *   2: solve p' = p - p_static (P1)
				   *   3: solve p' = p - p_static (P2)
				   *
				   * */
    FLOAT   dhx;		  /* dh x */
    FLOAT   dhz;		  /* dh z */
    int     temp_viscosity;       /* viscosity relates temperate type:
				   *   0: const A
				   *   1: T = -10
				   *   2: use T
				   * */
    FLOAT constant_A;

    
#if USE_MG 
    BOOLEAN use_mg_F;		  /* Use Multigrid solver */
    BOOLEAN use_mg_Ap;		  /* Use Multigrid solver */
    BOOLEAN use_mg_Qp;		  /* Use Multigrid solver */
#endif /* USE_MG */
    BOOLEAN use_moc;		  /* Use method of characterics */
    BOOLEAN use_Fu;		  /* Use Fu & Qu in the preconditioner, default: True */
    BOOLEAN use_Zp;		  /* Use Zp in the preconditioner, default: True */
    BOOLEAN implicit_convect;     /* Implicit scheme for convetion term,
				   *   default True */
    //BOOLEAN use_symetric;	  /* Mat and PC is symetric, for symetric check */

    /* Faster code */
    BOOLEAN extern_force;	  /* Extern force, default ture! */

    /* Solver options */
    BOOLEAN non_linear;		  /* Nonlinear iteration */
    BOOLEAN noniter_temp;	  /* nonlinear iterative of velocity and temperature */
    FLOAT non_tol;		  /* Nonlinear iteration tolerance */
    FLOAT non_sub_tol;		  /* Nonlinear iter: sub linear problem tolerance */

    INT pre_refines;		  /* Grid pre-refine level */
    INT max_tstep;		  /* Max time step */
    INT max_nonstep;		  /* Max nonlinear interation step,
				   * -1 means P2P1 is skiped */
    INT min_nonstep;		  /* Min nonlinear interation step */
    INT max_nonstep0;		  /* Max nonlinear interation step for first step,
				   * negtive means using max_nonstep instead. */
    INT newton_start;		  /* Newton start step */
    INT newton_start0;		  /* Newton start step for first step,
				   * negtive means using max_nonstep instead. */
    INT step_span;		  /* Step span to output geo file */
    INT mem_max;		  /* Max memory per process */
    INT n_bdry_layer;		  /* # of boundary layers */
    /* INT moc_quad_order;		  /\* MOC quad order *\/ */
    /* INT moc_quad_nelem;		  /\* MOC quad nelem *\/ */
    BOOLEAN compute_error;	  /* Compute error */

    FLOAT grid_coord_unit;	/* Grid coord unit [m] */

    char *fn;			  /* Mesh file */
    char *resume_mesh;		  /* Resume mesh file */
    char *resume_data;		  /* Resume data file */
    char *Stokes_opts;		  /* Solver Stokes options */
    char *fo_opts;		  /* Solver Stokes options */
    char *proj_opts;		  /* Solver Stokes options */
    char *F_opts;		  /* Solver F options*/
    char *Fu_opts;		  /* Solver Fu options*/
    char *Ap_opts;		  /* Solver Ap options*/
    char *Qp_opts;		  /* Solver Qp options*/
    char *Fp_opts;		  /* Solver Fp options*/
    char *T_opts;		  /* Solver temperature opts */

    char *topo_file;		/* Topo and beta2 file  */
    double topo_shift_x;	/* Topo shift x direction, unit: km */
    double topo_shift_y;

    /* Ref solution file */
    BOOLEAN compute_error_ref;	  /* Compute error to ref */
    char *ref_sol_file;
    
    /* PC levels */
    int n_pc_level;
    char *fn_L[PC_MAX_LEVEL];			  /* Mesh file */
    char *tria_file_L[PC_MAX_LEVEL];
    char *vert_file_L[PC_MAX_LEVEL];
    char *layer_file_L[PC_MAX_LEVEL];
    char *nodeZ_file_L[PC_MAX_LEVEL];
    char *dual_file_L[PC_MAX_LEVEL];
    char *nc_file_L[PC_MAX_LEVEL];

    BOOLEAN pc_reuse_L[PC_MAX_LEVEL];
    int pc_type_L[PC_MAX_LEVEL];

    char *pc_smoother_Lp;
    char *pc_smoother_L[PC_MAX_LEVEL];


	/* debugging controls */
	BOOLEAN output_non_iter;
	BOOLEAN output_temp_init;
	BOOLEAN output_parted;
	BOOLEAN output_beta;
	BOOLEAN output_melt;
	BOOLEAN output_fv_vert;
	
} NSParams;





/* Surface bases */
typedef struct SURF_BAS_ {
    DOF_TYPE *type;
    DOF *dof;
    BOOLEAN *rotated;
    char dim;
} SURF_BAS;



typedef void (USER_FUNC_ICE)(NSSolver *ns);

/* Main solver */
struct NSSolver_ {
    GRID *g;			  /* Grid */
    DOF **u;			  /* velocity ptr */
    DOF *u_edge;		  /* velocity on edge, discontinous */
    DOF **p;			  /* pressure ptr */
    DOF **T;			  /* conformation tensor */
#if STEADY_STATE || TIME_DEP_NON 
    DOF *du;			  /* delta u in non-linear iteration */
    DOF *dp;			  /* delta p in non-linear iteration */
    DOF *dT;			  /* delta C in non-linear iteration */
#endif /* STEADY_STATE || TIME_DEP_NON */
    DOF **gradu;		  /* gradient velocity ptr */
    DOF **lapu;			  /* laplace velocity ptr */
    DOF **gradp;		  /* gradient pressure ptr */
    DOF *Gradu;		  /* gradient temperature ptr */
    DOF *f;			  /* source term momentum  */
    DOF *u_queue[3];		  /* velocity at different time */
    DOF *p_queue[3];		  /* pressure at different time */
    DOF *T_queue[3];		  
    DOF *gradu_queue[3];	  /* gradient velocity at different time */
    DOF *div0;			  /* velocity divergence on element */
    DOF *u_shape;		  /* velocity shape DOF */
    DOF *p_shape;		  /* pressure shape DOF */
    DOF *un_shape;		  /* velocity shape DOF */
    DOF *bub_shape;		  /* velocity shape DOF */
    DOF *T_shape;		  /* proj Gradient vel shape DOF */
    DOF *gn[3];			  /* Outflow bdry condition for velocity & pressure*/
    DOF *wind;			  /* predicted velocity */
    DOF *dH;			  /* coord change */
    DOF *dHt;			  /* coord change rate */
    DOF *beta;			  /* slip coef */
    DOF *bottom_normal;		  /* bedrock normal: outwards */
    DOF *top_normal;		  /* top surf normal */
    DOF *lateral_normal;	  /* lateral normal */
    DOF *coord;			  /* coord(P1) */
    DOF *coord_sigma;		  /* coord(P1), sigma */
    DOF *p_static;		  /* hydro static pressure */
    DOF *nu;			  /* Viscosity for interp */

    INT pinned_node_id;	          /* Vert index of pinned node at rank 0
				   * -1 for no pinned. */
    MAP *Vmap;			  /* velocity dof map */
    MAP *Pmap;			  /* pressure dof map */
    MAP *T_map;

    /*          | F  Bt|
     *  matNS = |      |
     *          | B  C |
     *  */
    MAT *matNS;			  /* matrix of the coupled problem  */
    MAT *matF;			  /* matNS[0][0] */
    MAT *matBt;			  /* matNS[0][1] */
    MAT *matB;			  /* matNS[1][0] */
    MAT *matC;			  /* matNS[1][1] */

#if USE_SLIDING_BC
    DOF *Un;			  /* Langrangian Mutiplier of u dot n */
    DOF *dUn;			  /* d Un */
    
    MAP *UNmap;			  /* pressure dof map */

    MAT *matUn;
    MAT *matUnT;
    MAT *matUnD;
#endif


    MAT *matT;			  /* Mat T, bottom free */
    VEC *rhsT;
    
    SOLVER *solver_u;		  /* solver of the coupled problem */
    SOLVER *solver_T;		  /* solver of conformation tensor */
    SOLVER *pc;			  /* preconditioner */

    //OCT_TREE *og;		/* Oct tree for search */

#if USE_MG 
    MULTI_GRID *mg;		  /* Multi grid solver */
#endif /* USE_MG */

    /* Temp solver, Constrains */
    BOOLEAN *T_mask;		/* Temperature constrains mask */
    BOOLEAN *T_actc;		/* Temperature active constrain */
    FLOAT *T_cntr;		/* Temperature constrains value */

    /* GradS */
    DOF *dof_gradS;		/* grad S */

    /* Depth */
    DOF *depth_P1;		/* Depth P1 */
    DOF *depth_T; 		/* Depth of T fe type */
	DOF *height;		/* Height P1 */
	DOF *sigma_z;

    /* Variables */
    FLOAT non_res;		  /* nonlinear residual */
    FLOAT *time;		  /* time */
    FLOAT time_queue[3];	  /* time queue */
    FLOAT *dt;			  /* time step ptr
				   * dt_{n} = t_{n+1} - t_{n} */
    FLOAT dt_queue[2];		  /* time step at different time */
    INT tstep;			  /* current time step */
    int viscosity_type;	          /* Viscosity types: const, T-indep, T-dep, ... */
    LTYPE ltype;	          /* Picard or Newton */

    SURF_BAS *surf_bas;
    SURF_BAS *surf_bas_2d;

    /* Layred mesh */
    LAYERED_MESH *gL;		  /* Layered mesh */
    /* MG_BLOCK_DOFS *bk;		  /\* Block dofs for line smoothing *\/ */
    /* DOF *coord; */

    NSParams *ns_params;	  /* Parameter list */
    USER_FUNC_ICE *init_func;	  /* init function */
    USER_FUNC_ICE *mark_bdry_temp; /* Mark bdry by tempertature */
    USER_FUNC_ICE *check_div;      /* Check divergence */
} ;


typedef struct NSPCSolver_ {
    
} NSPCSolver;


/*
 * ================================================================================
 * 
 *                 Global parameters
 * 
 * ================================================================================
 * */
extern NSParams *ns_params;

enum { VIS_CONST  = 0,
       VIS_STRAIN = 1,
       VIS_TEMP	  = 2};

extern FLOAT nu_max;
extern FLOAT nu_min;

extern FLOAT eps_height;

typedef struct TIME_LOCK_ {
    BOOLEAN locked;
    FLOAT time;
} TIME_LOCK;

extern TIME_LOCK *time_lock; 


#define SIMPLE_TEST 0  /* simple test for convergence */
#if SIMPLE_TEST
#  define SIMPLE_INIT_SLOP				\
    FLOAT _alpha_ = 0.;					\
    FLOAT _length_ = 5.;				\
    FLOAT tan_alpha = tan(M_PI * _alpha_ / 180);
#endif


/*
 * ================================================================================
 * 
 *                 subroutines
 * 
 * ================================================================================
 * */
# define FUNC_T_DECLARE(func_xyz)					\
    void func_xyz##_t(FLOAT x, FLOAT y, FLOAT z, FLOAT t, FLOAT *values);    

/* ins-solver.c */
NSParams *phgParametersCreate();     
NSSolver *phgNSCreate(GRID *g, LAYERED_MESH *gL, NSParams *ns_params);
void phgNSFinalize(NSSolver **ns);
void phgNSTimeAdvance(NSSolver *ns, FLOAT time, int tstep);
INT phgNSPinNode(NSSolver *ns);

void phgNSInitSolverU(NSSolver *ns);
void phgNSReInitSolverU(NSSolver *ns);
void phgNSBuildSolverUMat(NSSolver *ns);
void phgNSBuildSolverURHS(NSSolver *ns);
void phgNSSolverUAssemble(NSSolver *ns);
void phgNSDestroySolverU(NSSolver *ns);

FLOAT *get_gbas_product_stokes(const FLOAT *gi, const FLOAT *gj,
			       const FLOAT *gu, LTYPE ltype); 
FLOAT *get_gbas_product_fo(const FLOAT *gi, const FLOAT *gj,
			       const FLOAT *gu, LTYPE ltype); 


void build_prism_elems(GRID *g, LAYERED_MESH *gL);
void phgNSBuildSolverUMatPrism(NSSolver *ns);
void phgNSBuildSolverURHSPrism(NSSolver *ns);



void buildFOMatPrism(NSSolver *ns, DOF *dof_uv, DOF *dof_du, SOLVER *solver_uv);
void buildFORHSPrism(NSSolver *ns, DOF *dof_uv, DOF *dof_du, SOLVER *solver_uv);
void buildProjGuPrism(NSSolver *ns, SOLVER *solver, VEC **vec, const BYTE *components);

void getPecletNum(GRID *g, DOF *u, FLOAT nu, int order);
void phgDofMaxMin(DOF *u, FLOAT *umax, FLOAT *umin);
void estimate_error(NSSolver *ns, DOF *error);
void phgResumeLogUpdate(GRID *g, FLOAT *time, int *tstep, char *mesh_file, char *data_file);
void phgResumeStage(GRID *g, FLOAT *time, int *tstep, char *mesh_file, char *data_file);

/* PC */
NSPCSolver *nsPCCreate(NSSolver *ns, NSParams *ns_params);
void nsPCBuildMat(NSPCSolver *pcSolver, int non_step);

/* save load */
void save_dof_data(GRID *g, DOF *dof, const char *file);
void load_dof_data(GRID *g, DOF *dof, const char *data_file, const char *mesh_file);
void load_dof_data2(GRID *g, DOF *dof, const char *data_file, const char *mesh_file);
void save_dof_data3(GRID *g, DOF *dof, const char *file);
void load_dof_data3(GRID *g, DOF *dof, const char *data_file, const char *mesh_file);
void ns_dof_copy(NSSolver *ns, DOF *u, DOF *p);

void save_state(int tstep, double time, int ndof, DOF **dofs, char **dof_names);
void load_state(int tstep, double *time, int ndof, DOF **dofs, char **dof_names);



/* ins-bc.c */
int my_bc_map(int bctype);
void setFlowParameter(FLOAT Re_in, FLOAT nu_in, FLOAT time_in);
void setFuncTime(FLOAT time_in);
void adjust_time(FLOAT delta_time);
void restore_time(void);
void set_boundary_mask(NSSolver *ns);
void func_init_params(double Lx0, double alpha0);

void func_init(FLOAT x, FLOAT y, FLOAT z, FLOAT *u);
void func_u(FLOAT x, FLOAT y, FLOAT z, FLOAT *u);
void func_p(FLOAT x, FLOAT y, FLOAT z, FLOAT *p);
void func_gradp(FLOAT x, FLOAT y, FLOAT z, FLOAT *gradp);
void func_gradu(FLOAT x, FLOAT y, FLOAT z, FLOAT *gradu);
void func_f(FLOAT x, FLOAT y, FLOAT z, FLOAT *f);
void func_g1(FLOAT x, FLOAT y, FLOAT z, FLOAT *g);
void func_g2(FLOAT x, FLOAT y, FLOAT z, FLOAT *g);
void func_g3(FLOAT x, FLOAT y, FLOAT z, FLOAT *g);
void func_T(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);
void func_fT(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);
void func_beta(FLOAT x, FLOAT y, FLOAT z, FLOAT *beta);
void func_q(FLOAT x, FLOAT y, FLOAT z, FLOAT *value);
void func_xyz_(FLOAT x, FLOAT y, FLOAT z, FLOAT *coord);
void func_s(FLOAT x, FLOAT y, FLOAT z, FLOAT *s);
void func_b(FLOAT x, FLOAT y, FLOAT z, FLOAT *b);
void func_a(FLOAT x, FLOAT y, FLOAT *a);
void func_normal(FLOAT x, FLOAT y, FLOAT z, FLOAT *a);

FUNC_T_DECLARE(func_u);
FUNC_T_DECLARE(func_p);
FUNC_T_DECLARE(func_gradu);
FUNC_T_DECLARE(func_gradp);
FUNC_T_DECLARE(func_f);
FUNC_T_DECLARE(func_T);
FUNC_T_DECLARE(func_q);


/* temp-solver.c */
void phgNSInitSolverT(NSSolver *ns);
void phgNSBuildSolverT(NSSolver *ns);
void phgNSDestroySolverT(NSSolver *ns);
void phgNSBuildSolverTMat(NSSolver *ns, BOOLEAN init_T);
void phgNSBuildSolverTRHS(NSSolver *ns, BOOLEAN init_T);
void phgNSSolverTBuildConstrain(NSSolver *ns);
void find_melt_region(NSSolver *ns);
void phgNSSolverTSolve(NSSolver *ns, BOOLEAN init_T);
void phgNSTempInit(NSSolver *ns);
void proj_gradu(NSSolver *ns, DOF *gradu, const BYTE *components, int type);


/* ins-pcd.c */
void phgNSInitPc(NSSolver *ns);
void phgNSBuildPc(NSSolver *ns);
void pc_proc(void *ctx, VEC *b0, VEC **x0);
void phgNSDestroyPc(NSSolver *ns);
void estimate_error(NSSolver *ns, DOF *error);

void phgDofSetName(DOF *dof, const char *name);
double elapsed_time(GRID *g, BOOLEAN flag, double mflops);
void phgDofRemoveConstMod(GRID *g, DOF *u);

/* ins-utils.c */
void NsSolver_Options();
FLOAT dofNormL2(DOF *dof);
void dof_norm_L2(DOF *dof);
void checkBdry(GRID *g);
int my_bc_map(int bctype);
void dof_range(DOF *u);
void output_bottom_dofs(NSSolver *ns, int tstep);
DOF *proj_G1(NSSolver *ns, DOF *p0, DOF **p1_ptr);
void check_div(DOF *gradu, DOF **divu);
void plot_residual(NSSolver *ns, VEC *residual, int non_step);
void plot_surf_force(NSSolver *ns);
void mark_inactive(NSSolver *ns);
void save_solution_dat(NSSolver *ns);
void compute_error_ref(NSSolver *ns);
void init_by_ref_solution(NSSolver *ns);
void install_error_log();
void get_p_static(NSSolver *ns);
void get_p_static_prism(NSSolver *ns);
int phgCompEdge(const void *p0, const void *p1);
int phgCompEdgeMark(const void *p0, const void *p1);
void get_avg_n(GRID *g, DOF *sur_normal);


void load_ref_solution();
const FLOAT *interp_ref_solution(FLOAT *X, FLOAT sigma);


/* ins-mg.c */
void build_mg_levels(NSSolver *ns, int i_slv);

/* upwind.c */
void cr1_upwind(SIMPLEX *e, DOF *wind, FLOAT nu, 
		FLOAT Samarskij_alpha, int order, FLOAT *values);
void fv_upwind(SIMPLEX *e, DOF *wind, FLOAT nu, 
		FLOAT Samarskij_alpha, int order, FLOAT *values);

/* ice-grid.c */
void ice_grid(GRID *g);
void ice_monitor(NSSolver *ns, int nonstep);
BOOLEAN iceParter(GRID *g, int nprocs, DOF *weights, FLOAT power);

/* slip-bdry.c */
#define get_surface_bases_2d(ns, u_type) get_surface_bases(ns, u_type, 2)
#define get_surface_bases_3d(ns, u_type) get_surface_bases(ns, u_type, 3)
SURF_BAS *get_surface_bases(NSSolver *ns, DOF_TYPE *u_type, int rdim);

#define trans_left_3d(A, ncol, lda, T) trans_left(A, ncol, lda, 3, T) 
#define trans_leftT_3d(A, ncol, lda, T) trans_leftT(A, ncol, lda, 3, T) 
#define trans_rightT_3d(A, ncol, lda, T) trans_rightT(A, ncol, lda, 3, T) 

#define trans_left_2d(A, ncol, lda, T) trans_left(A, ncol, lda, 2, T) 
#define trans_leftT_2d(A, ncol, lda, T) trans_leftT(A, ncol, lda, 2, T) 
#define trans_rightT_2d(A, ncol, lda, T) trans_rightT(A, ncol, lda, 2, T) 
    
void trans_left(FLOAT *A, int ncol, int lda, int dim, const FLOAT *Trans); 
void trans_leftT(FLOAT *A, int ncol, int lda, int dim, const FLOAT *Trans); 
void trans_rightT(FLOAT *A, int ncol, int lda, int dim, const FLOAT *Trans); 

void rotate_dof_bases(DOF *u, SURF_BAS *surf_bas, BOOLEAN forward);
void dof_set_normal_data(DOF *u_h, SURF_BAS *surf_bas);
extern FLOAT trans_eye[Dim*Dim];
DOF *get_bottom_normal(NSSolver *ns);
DOF *get_lateral_normal(NSSolver *ns);
void get_orth_bases(FLOAT *normal, FLOAT *T);
void update_floating(NSSolver *ns);


/* moving-mesh.c */
void get_surf_dH(NSSolver *ns);
void get_moved_coord(NSSolver *ns, int tstep);
void move_dof(GRID *g, DOF *dz, DOF *u);
void move_mesh(NSSolver *ns);
void get_layer_height(FLOAT *H, int nv, const FLOAT *ratio, FLOAT h0, FLOAT h1);
void get_height_depth(NSSolver *ns);

/* fv-solver.c */
typedef void (*DOF_USER_FUNC_T)(double x, double y, double z, 
				double t, double *values);

void fv_solver_init(void **fv_data,
		    const char *mesh_file,
		    FLOAT (*verts)[4],
		    const DOF_USER_FUNC_T func_f);
void fv_update(void *fv_data,
	       const double *H, 
		   FLOAT (*verts)[4],
	       const double *U, 
	       double *dH, 
	       double *U_vert, 
	       double t);

/* update_surf.c */
void struct_mesh_init(GRID *g);
void struct_mesh_reinit(GRID *g);
void struct_mesh_update(NSSolver *ns, int tstep, double t);

/* core-shallow.c */
void core_SIA(NSSolver *ns);
void get_moved_coord_SIA(NSSolver *ns, int tstep);

/* core-FO.c */
void core_FO(NSSolver *ns, int tstep);
void reconstruct_velocity_w(NSSolver *ns);

/* Netcdf tools */
void load_topo_data(const char *filename);
#define interp_topo_data(v, x, y) \
    interp_topo_dataT(v, x, y, 0);
double interp_topo_dataT(char var_type, double x, double y, double t);


/* Problem dependent */
void iceInit(GRID *g, LAYERED_MESH **gL);
void iceInitPC(GRID *g, int iLevel, LAYERED_MESH **gL);
void iceSetUserFunc(NSSolver *ns);
void func_ice_slab(FLOAT x, FLOAT y, FLOAT z, FLOAT *coord);
FLOAT func_ice_topg(FLOAT x, FLOAT y);
void map_coord_inv(FLOAT x, FLOAT y, FLOAT *x_, FLOAT *y_);
FLOAT get_effective_viscosity(const FLOAT *gu, FLOAT T, FLOAT p,
			      BOOLEAN initialiszed);
void load_case_options();



/* Missing */
#define phgExportTecplot(...)
#define phgPartUserSetFunc(...)

/* redundent sizeof */
#define PHG_ALLOC(p, n) p = phgAlloc(n * sizeof(*p));
#define PHG_CALLOC(p, n) p = phgCalloc(n, sizeof(*p));
#define PHG_REALLOC(p, n1, n0) p = phgRealloc_(p, n1 * sizeof(*p), \
                                               n0 * sizeof(*p));


/*
 * ================================================================================
 * 
 *                 Utils
 * 
 * ================================================================================
 * */

#define OUTPUT_DIR "output/"
#define X_DIR 0
#define Y_DIR 1
#define Z_DIR 2
#define UN_DIR X_DIR
#define LN_DIR Y_DIR
#define DDim (Dim*Dim)

typedef INT EdgeVerts[2];
typedef INT EdgeMarks[3];


/* Quad macros */
#define NbasFace(u) (3 * (u->type->np_vert + u->type->np_edge)	\
		     + u->type->np_face)
#define SQUARE(x) ((x)*(x))
#define INNER_PRODUCT2(p, q)			\
    (*(p    ) * *(q    ) +			\
     *(p + 1) * *(q + 1))
#define INNER_PRODUCT(p, q)			\
    (*(p    ) * *(q    ) +			\
     *(p + 1) * *(q + 1) +			\
     *(p + 2) * *(q + 2))
#define CROSS_PRODUCT(a, b, n) {                      \
        n[0] =  (a[1] * b[2] - b[1] * a[2]);        \
        n[1] = -(a[0] * b[2] - b[0] * a[2]);        \
        n[2] =  (a[0] * b[1] - b[0] * a[1]);        \
    }


#define BasisOrder(u, e, i) (!DofIsHP(u) ? (u)->type->order :		\
			      (u)->hp->info->types[(u)->hp->elem_order[e->index]]->order)
#define Bzero(v) bzero(v, sizeof(v));
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

#define DATA_FILE_SURFIX     {			\
  sprintf(data_u, "%s.u", data_file);		\
  sprintf(data_p, "%s.p", data_file);		\
  sprintf(data_T, "%s.T", data_file);		\
  }

#define GET_DOF_TYPE(dof_type, dof_name) {	\
    char s[128];				\
    phgOptionsPush();				\
    sprintf(s, "-dof_type %s", dof_name);	\
    phgOptionsSetOptions(s);			\
    dof_type = DOF_DEFAULT;			\
    phgOptionsPop();				\
  }

#define phgDofBY(b, y) {			\
    	INT i_, n_ = DofGetDataCount(*y);	\
	FLOAT *p_ = DofData(*y);		\
	for (i_ = 0; i_ < n_; i_++)		\
	    *(p_++) *= b;			\
    }



/* Periodic sync */
#define PERIODIC_SYNC(dof)					\
    if (g->period != NULL) {					\
	MAP *map_ = phgMapCreate(dof, NULL);			\
	VEC *vec_ = phgMapCreateVec(map_, 1);			\
	phgMapDofToLocalData(map_, 1, &dof, vec_->data);	\
	phgMapLocalDataToDof(map_, 1, &dof, vec_->data);	\
	phgVecDestroy(&vec_);					\
	phgMapDestroy(&map_);					\
    }


/* TODO: merge isop */
#define phgGeomGetCurvedJacobianAtLambda(g, e, p, det)	\
    *det = phgGeomGetVolume(g, e) * 6.;			\

#define phgQuadGetBasisCurvedGradient(e, u, i, quad, q) \
    phgQuadGetBasisGradient(e, u, i, quad) + q*Dim


#define SORT_TWO_VERT(v0, v1) {			\
	int vv;					\
	if (v0 > v1) {				\
	    vv = v0; v0 = v1; v1 = vv;		\
	}					\
    }



/*
 * ================================================================================
 * 
 *                 Check macros
 * 
 * ================================================================================
 * */

#if 1
#define SHOW_M(matN, mat_m, mat_n) {					\
	phgInfo(2, "\n### rank: %d\n", g->rank);			\
	int i, j;							\
	phgInfo(2, " --- "#matN":(%3d * %3d)\n", mat_m, mat_n);		\
	for (i = 0; i < mat_m; i++) {					\
	    for (j = 0; j < mat_n; j++){				\
		phgInfo(2, "%30.15e, ", *(matN + i * (mat_n) + j));	\
	    }								\
	    phgInfo(2, "\n");						\
	}								\
    }

#define SHOW_M3(matN, mat_k, mat_m, mat_n) {				\
	phgInfo(2, "\n### rank: %d\n", g->rank);			\
	int i, j, k;							\
	phgInfo(2, "--- --- %15s :(%3d * %3d)\n", #matN, mat_m, mat_n);	\
	for (k = 0; k < mat_k; k++) {					\
	    phgInfo(2, "  comp: %d\n", mat_k);				\
	    for (i = 0; i < mat_m; i++) {				\
		phgInfo(2, "    ");					\
		for (j = 0; j < mat_n; j++){				\
		    phgInfo(2, "%10f, ", *(matN + k * mat_m * mat_n +	\
					   i * mat_n + j));		\
		    if (mat_n > 10 && j%5 == 4)				\
			phgInfo(2, "\n    ");				\
		}							\
		phgInfo(2, "\n");					\
	    }								\
	}								\
    }

#define SHOW_V(vec, vec_n) { int _i;		\
	phgInfo(2, " --- "#vec":(%3d)\n", vec_n);	\
	for (_i = 0; _i < vec_n; _i++) {	\
	    phgInfo(2, "%10f, ", *(vec + _i));	\
	}					\
	phgInfo(2, "\n");				\
    }

#define SHOW_iV_(verb_, vec, vec_n) { int _i;		\
	phgInfo(verb_, " --- "#vec":(%3d)\n", vec_n);	\
	for (_i = 0; _i < vec_n; _i++) {	\
	    phgInfo(verb_, "%4d, ", *(vec + _i));	\
	}					\
	phgInfo(verb_, "\n");				\
    }			

#define SHOW_iM(matN, mat_m, mat_n) {				\
	int i, j;						\
	phgInfo(2, " --- "#matN":(%3d * %3d)\n", mat_m, mat_n);	\
	for (i = 0; i < mat_m; i++) {				\
	    for (j = 0; j < mat_n; j++){			\
		phgInfo(2, "%5d, ", *(matN + i * (mat_n) + j));	\
	    }							\
	    phgInfo(2, "\n");					\
	}							\
    }


#endif

#define PRINT_ELAPSED_TIME(str, ...)	\
    phgPrintf(str);				\
    elapsed_time(__VA_ARGS__);
    

#if 0
#define DOF_SCALE(u, desp) {}
#elif 0
# define DOF_SCALE(u, desp)					\
    dof_range(u);
#elif 1
# define DOF_SCALE(u, desp)					\
    phgPrintf("    %s: [%16.8e, %16.8e]\n",			\
	      (u)->name,					\
	      phgDofMinValVec(u),				\
	      phgDofMaxValVec(u)				\
	      );						
#elif 0
# define DOF_SCALE(u, description) {				\
	char trimed[100];					\
	strncpy(trimed, __FUNCTION__, 8);			\
	trimed[8]='\0';						\
	phgPrintf("   ------------------------------\n"		\
		  "   %-10s  /* %s */\n"			\
		  "     func: %-10s, line: %03d\n"		\
		  ,#u":", description,				\
		  trimed, __LINE__);				\
	phgPrintf("     [%16.8e, %16.8e] (max,min); \n"		\
		  "     [%16.8e, %16.8e] ( L1, L2); \n",	\
		  phgDofMaxValVec(u), phgDofMinValVec(u),	\
		  phgDofNormL1(u), phgDofNormL2(u)		\
		  );						\
    }
#else
# define DOF_SCALE(u, description) {				\
	char trimed[100];					\
	DOF *_tmp_grad = phgDofGradient(u, NULL, NULL, NULL);	\
	strncpy(trimed, __FUNCTION__, 8);			\
	trimed[8]='\0';						\
	phgPrintf("   ------------------------------\n"		\
		  "   %-10s  /* %s */\n"			\
		  "     func: %-10s, line: %03d\n"		\
		  ,#u":", description,				\
		  trimed, __LINE__);				\
	phgPrintf("     [%16.8e, %16.8e] (max,min); \n"		\
		  "     [%16.8e, %16.8e] ( L1, L2); \n",	\
		  phgDofMaxValVec(u), phgDofMinValVec(u),	\
		  phgDofNormL1(u), phgDofNormL2(u)		\
		  );						\
	phgPrintf("   %-10s \n"					\
		  ,"grad("#u"):");				\
	phgPrintf("     [%16.8e, %16.8e] (max,min);\n"		\
		  "     [%16.8e, %16.8e] ( L1, L2);\n\n",	\
		  phgDofMaxValVec(_tmp_grad),			\
		  phgDofMinValVec(_tmp_grad),			\
		  phgDofNormL1(_tmp_grad),			\
		  phgDofNormL2(_tmp_grad)			\
		  );						\
	phgDofFree(&_tmp_grad);					\
    }
#endif

#if 1
#define sayHello(desp)							\
    {									\
	char hostname_[256];						\
	gethostname(hostname_, sizeof(hostname_));			\
	phgInfo(1, "Host %-15s, PID %10d, "				\
		"get line:%5d, mem: %0.4lfMB, /* %-20s */\n",		\
		hostname_, getpid(), __LINE__,				\
		phgMemoryUsage(g, NULL) / (1024.0 * 1024.0), desp);	\
	fflush(stdout);							\
	MPI_Barrier(g->comm);						\
    }
#else
#define sayHello(desp)
#endif


#define ERROR_MSG(desp) \
    phgError(1, "Error: file: %s, func: %s, line: %s\n"	\
	     "   %s not avaliable", __FILE__,		\
	     __FUNCTION__, __LINE__, desp);


#define INS_H
#endif
