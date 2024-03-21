/*
 *************************************************
 *                                               *
 * Incompressible Naiver-Stokes euqation solver. *
 *                                               *
 *************************************************
 *
 *
 * Usage & Options:
 *
 * 1. General Utils
 * 1.1 Save & resume: save & load mesh and DOF of type Pn
 * 1.2 Mesh file
 * 1.3 Max mem usage limit
 * 1.4 step span: step interval to save solution
 *
 *
 * 2. Physical flow configuration
 * 2.1 Reynolds num:
 *     When variables are dimensionless, use Reynolds num instead of
 *     viscosity num.
 * 
 * 2.2 Enclosed flow & pin node:
 *     For enclosed flow, all boundary is to set to be Dirichlet type, thus 
 *     generate a singular linear system with null space {1}. To avoid this 
 *     singularity, a vert is pinned as Dirichlet type.
 *
 * 3. Discretization
 * 3.1 Basic
 *     * time end:
 *     * dt: time step
 *     * max time step:
 *     * theta: fractional time step, 1. for Euler backward, .5 for CN
 *     * [up]type: dof type of Velocity and pressure
 *     
 * 3.2 convection treatment:
 *     (i)   explicit 
 *     (ii)  implicit
 *     (iii) MOC 
 * 
 * 3.3 Misc
 *     * symetric: use symetric preconditioner and MINRES iter
 *     * no external force: no need to compute quadrature of external force
 *     * eps on diag of pressure mat: small pertubation on diag of C to statisfy
 *       the requirement that diag should not be zero in some solver, like PETSC
 *       
 *     
 * 4. Solver option
 *
 * 4.1 non-linear solver
 *     * non linear tolerance
 *     * sub linear solver tolerance
 *     * tolerance of picard to newton: the non-linear iter type change from
 *       picard to newton if non residual is smaller than given tol.
 * 4.2 PCD preconditioner
 *     * use Pressure convection-diffusion preconditioner
 *     * Use Fu to solve the sub mat
 * 4.3 sub problems solver:
 *     * solver F : vector convection-diffusion solver
 *     * solver Fu: scalar convection-diffusion solver
 *     * solver Fp: Dirich scale 
 *     * solver Ap: pressure poisson solver
 *     * solver Qp: pressure mass solver
 *
 * Notes:
 * 1. boundary condition for PCD PC.
 *    Different choice of bdry condition could be made for pressure convc-diffu
 *    problem, for enclose flow, the best may be all Neumann with a pinned node,
 *    for open flow, the best may be (i) Robin at inflow, (ii) scaled Dirich
 *    at outflow, (iii) Neumann at wall.
 *    The difference of pressure flow boundary condition lies in three place:
 *    (1) initializing mat: choice different boundary to be Dirichlet type.
 *    (2) building matrix:
 *        Dirichlet condition needs taken care of when adding entries to mat; 
 *        Robin condition needs an extra quad on face.
 *    (3) preconditioning process: scaling for Dirichlet condition.
 *    
 *
 * $Id: */


#include "ins.h"
/* #include "io.h" */
/* #include "periodic.h" */
#if USE_PETSC
#  include <petscsys.h>
#  include <petscviewer.h>
#endif


#define _p new_params
#define _nsp (ns->ns_params)
#define _pcd (ns->pcd)

#define SUB_SOLVER_VERB -1

/* non linear solver type: picard or newton */
static const char *noniter_name[] = {
  "picard", "newton", NULL};
static const char *core_names[] = {
    "stokes", "sia", "ssa", "fo", "debug1", "debug2", NULL
};
static const char *pc_type_names[] = {
    "stokes", "schur", "block11", "fo", "debug1", NULL
};



/*******************************************************/
/* Get default parameters and user defined parameters. */
/*******************************************************/
NSParams *
phgParametersCreate()
{
    NSParams *new_params = phgCalloc(1, sizeof(*new_params));
    int iLevel;

    /* default settings */
    _p->fn = "../test/cube.dat";

    /* Get user defined parameter from file*/
    /* geom & mesh files */
    phgOptionsRegisterFilename("mesh_file", "Mesh file", (char **)&_p->fn);
    phgOptionsRegisterFilename("ref_sol_file", "Ref solution file", (char **)&_p->ref_sol_file);

    phgOptionsRegisterNoArg("reduce_mesh", "Reduce mesh", &_p->reduce_mesh); 
    phgOptionsRegisterNoArg("layered_mesh", "Layered mesh", &_p->layered_mesh);
    phgOptionsRegisterFloat("grid_coord_unit", "Grid coord unit", &_p->grid_coord_unit);
    phgOptionsRegisterInt("periodicity", "Set periodicity", &_p->periodicity);
    phgOptionsRegisterInt("period_repeat", "Set periodicity repeat", &_p->period_repeat);
    phgOptionsRegisterInt("pre_refines", "Pre refines", &_p->pre_refines);
    phgOptionsRegisterInt("n_bdry_layer", "# of boundary layers", &_p->n_bdry_layer);
    phgOptionsRegisterNoArg("curved_bdry", "Use curved boundary", &_p->curved_bdry);
    phgOptionsRegisterNoArg("curved_refine", "curved refinement", &_p->curved_refine);
    phgOptionsRegisterNoArg("update_bdry_type", "update bdry type", &_p->update_bdry_type);


    /* Topo file */
    phgOptionsRegisterFilename("topo_file", "Topo file", (char **)&_p->topo_file);
    phgOptionsRegisterFloat("topo_shift_x", "Topo shift x", &_p->topo_shift_x);
    phgOptionsRegisterFloat("topo_shift_y", "Topo shift x", &_p->topo_shift_y);


    /*
     * Structure mesh, disabled
     * */
    /* phgOptionsRegisterNoArg("struct_mesh", "Structure mesh", &_p->struct_mesh); */
    /* phgOptionsRegisterInt("Nx", "Nx", &_p->Nx); */
    /* phgOptionsRegisterInt("Ny", "Ny", &_p->Ny); */
    /* phgOptionsRegisterInt("Nz", "Nz", &_p->Nz); */


    /*
     * PC options, disabled
     * */
    /* phgOptionsRegisterInt("n_pc_level", "# of PC levels", &_p->n_pc_level); */
    /* for (iLevel = 0; iLevel < PC_MAX_LEVEL; iLevel++) { */
    /* 	char str[1000]; */
    /* 	char desp[1000]; */

    /* 	sprintf(desp, "Mesh file Level %d", iLevel); */
	
    /* 	sprintf(str, "mesh_file_L%d", iLevel);  */
    /* 	phgOptionsRegisterFilename(str, desp, (char **)&_p->fn_L[iLevel]); */
    /* 	sprintf(str, "tria_file_L%d", iLevel);  */
    /* 	phgOptionsRegisterFilename(str, desp, (char **)&_p->tria_file_L[iLevel]); */
    /* 	sprintf(str, "vert_file_L%d", iLevel);  */
    /* 	phgOptionsRegisterFilename(str, desp, (char **)&_p->vert_file_L[iLevel]); */
    /* 	sprintf(str, "layer_file_L%d", iLevel);  */
    /* 	phgOptionsRegisterFilename(str, desp, (char **)&_p->layer_file_L[iLevel]); */
    /* 	sprintf(str, "dual_file_L%d", iLevel);  */
    /* 	phgOptionsRegisterFilename(str, desp, (char **)&_p->dual_file_L[iLevel]); */
    /* 	sprintf(str, "nodeZ_file_L%d", iLevel);  */
    /* 	phgOptionsRegisterFilename(str, desp, (char **)&_p->nodeZ_file_L[iLevel]); */
    /* 	sprintf(str, "netcdf_file_L%d", iLevel);  */
    /* 	phgOptionsRegisterFilename(str, desp, (char **)&_p->nc_file_L[iLevel]); */

    /* 	sprintf(str, "pc_reuse_L%d", iLevel);  */
    /* 	sprintf(desp, "Reuse grid in PC Level %d", iLevel); */
    /* 	phgOptionsRegisterNoArg(str, desp, &_p->pc_reuse_L[iLevel]);  */

    /* 	sprintf(str, "pc_type_L%d", iLevel);  */
    /* 	sprintf(desp, "PC type Level %d", iLevel); */
    /* 	phgOptionsRegisterKeyword(str, desp, pc_type_names, &_p->pc_type_L[iLevel]); */

    /* 	sprintf(str, "pc_smoother_L%d", iLevel);  */
    /* 	sprintf(desp, "PC smoother Level %d", iLevel); */
    /* 	phgOptionsRegisterString(str, desp, (char **)&_p->pc_smoother_L[iLevel]); */
    /* } */
    /* phgOptionsRegisterString("pc_smoother_Lp", "PC smoother Level -1",  */
    /* 			     (char **)&_p->pc_smoother_Lp); */


    /* FEM discretization */
    phgOptionsRegisterKeyword("-core_type", "Ice sheet core type",
				core_names, &_p->core_type); /* default Stokes */
    phgOptionsRegisterString("utype", "DOF type for velocity", &_p->utype_name);
    phgOptionsRegisterString("ptype", "DOF type for pressure", &_p->ptype_name);
    phgOptionsRegisterString("T_type", "DOF type for Temp", &_p->T_type_name);
    phgOptionsRegisterString("isop_type", "DOF type for Isop mapping", &_p->isop_type_name);
    phgOptionsRegisterString("utype_prv", "previous DOF type for velocity", &_p->utype_prv_name);
    phgOptionsRegisterString("ptype_prv", "previous DOF type for pressure", &_p->ptype_prv_name);
    phgOptionsRegisterString("T_type_prv", "previous DOF type for Temp", &_p->T_type_prv_name);
    phgOptionsRegisterFloat("nu", "viscosity number", &_p->nu);
    phgOptionsRegisterNoArg("enclosed_flow", "Enclosed flow", &_p->enclosed_flow);
    phgOptionsRegisterNoArg("pin_node", "Pin a node for enclosed flow", &_p->pin_node);
    phgOptionsRegisterNoArg("use_moc", "Use Method of Characterics", &_p->use_moc);
    phgOptionsRegisterNoArg("use_prism_elem", "Use Prism elements", &_p->use_prism_elem);
    phgOptionsRegisterNoArg("implicit_convect", "Implicit scheme for convetion term", 
			    &_p->implicit_convect); 
    phgOptionsRegisterInt("stab_scheme", "Pressure stab scheme", &_p->stab_scheme);
    phgOptionsRegisterInt("stab_remove_static", "Pressure stab remvoe p static", &_p->stab_remove_static);
    phgOptionsRegisterFloat("stab_nu", "Stab visocity", &_p->stab_nu);
    phgOptionsRegisterFloat("stab_alpha", "Stab alpha", &_p->stab_alpha);

    phgOptionsRegisterNoArg("use_slide", "Use sliding brdy", &_p->use_slide);
    phgOptionsRegisterInt("sliding_bdry_scheme", "Sliding bdry scheme", &_p->sliding_bdry_scheme);
    phgOptionsRegisterFloat("sliding_penalty", "Sliding penalty", &_p->sliding_penalty);

    phgOptionsRegisterInt("height_scheme", "Height solver scheme", &_p->height_scheme);
    phgOptionsRegisterInt("init_temp_type", "Init temperature type", &_p->init_temp_type);
    phgOptionsRegisterInt("temp_viscosity", "Temp viscosity type", &_p->temp_viscosity);
    phgOptionsRegisterFloat("constant_A", "constant_A", &_p->constant_A);
    phgOptionsRegisterNoArg("solve_temp", "Solve temperature", &_p->solve_temp);
    phgOptionsRegisterNoArg("solve_height", "Solve height", &_p->solve_height);
    phgOptionsRegisterNoArg("start_const_vis", "Start viscosity as constant", &_p->start_const_vis);
    phgOptionsRegisterNoArg("compensate_equ", "Compensate equations", &_p->compensate_equ);
    phgOptionsRegisterNoArg("compute_error", "Compute error", &_p->compute_error);
    phgOptionsRegisterNoArg("compute_error_ref", "Compute error to ref", &_p->compute_error_ref);
    phgOptionsRegisterFloat("dhx", "dh x", &_p->dhx);
    phgOptionsRegisterFloat("dhz", "dh z", &_p->dhz);


    /* Time discretization */
    phgOptionsRegisterFloat("dt", "Time step", &_p->dt0);
    phgOptionsRegisterFloat("time_start", "Time start", &_p->time_start);
    phgOptionsRegisterFloat("time_end", "Time end", &_p->time_end);
    phgOptionsRegisterFloat("theta", "Time fraction coef", &_p->Theta);
    phgOptionsRegisterInt("step_span", "Step span to output geo file", &_p->step_span);
    phgOptionsRegisterInt("max_time_step", "Max time step", &_p->max_tstep);

    /* Nonlinear system */
    phgOptionsRegisterNoArg("non_linear", "non linear problem", &_p->non_linear);
    phgOptionsRegisterFloat("non_tol", "Nonlinear iteration tolerance", &_p->non_tol);
    phgOptionsRegisterFloat("non_sub_tol", 
			    "Nonlinear iteration sub linear solver tolerance", &_p->non_sub_tol);
    phgOptionsRegisterInt("max_non_step", "Max nonlinear iteration step", &_p->max_nonstep);
    phgOptionsRegisterInt("min_non_step", "Min nonlinear iteration step", &_p->min_nonstep);
    phgOptionsRegisterInt("max_non_step0", "Max nonlinear iteration step for 1st step", &_p->max_nonstep0);
    phgOptionsRegisterInt("newton_start", "newton start", &_p->newton_start);
    phgOptionsRegisterInt("newton_start0", "newton start for 1st step", &_p->newton_start0);
    phgOptionsRegisterNoArg("noniter_temp", "non linear Temperature", &_p->noniter_temp);


    /* Linear system */
    phgOptionsRegisterFloat("eps_diagP", "Diagnal entries for pressure mat", &_p->eps_diagP);
    phgOptionsRegisterFloat("fp_scale", "Dirich scale for Fp", &_p->fp_scale);
    phgOptionsRegisterNoArg("use_PCD", "Use PCD preconditioner", &_p->use_PCD);
    phgOptionsRegisterNoArg("use_Fu", "Use Fu in the preconditioner", &_p->use_Fu);
#if USE_MG 
    phgOptionsRegisterNoArg("use_mg_F", "Use MG preconditioner", &_p->use_mg_F);
    phgOptionsRegisterNoArg("use_mg_Ap", "Use MG preconditioner", &_p->use_mg_Ap);
    phgOptionsRegisterNoArg("use_mg_Qp", "Use MG preconditioner", &_p->use_mg_Qp);
#endif /* USE_MG */

    /* Utils */
    phgOptionsRegisterInt("mem_max", "Max memory per process(MB)", &_p->mem_max);
    phgOptionsRegisterNoArg("resume", "resume from previous step", &_p->resume);
    phgOptionsRegisterNoArg("record", "save data to resume", &_p->record);

    phgOptionsRegisterNoArg("output_non_iter", "Output non[0123...].vtk", &_p->output_non_iter); 
    phgOptionsRegisterNoArg("output_temp_init", "Output  temp_init.vtk", &_p->output_temp_init); 
    phgOptionsRegisterNoArg("output_beta", "Output  beta.vtk", &_p->output_beta); 
    phgOptionsRegisterNoArg("output_melt", "Output  melt.vtk", &_p->output_melt); 
    phgOptionsRegisterNoArg("output_parted", "Output  parted.vtk", &_p->output_parted); 
    phgOptionsRegisterNoArg("output_fv_vert", "Output fv_vert.vtk", &_p->output_fv_vert); 


	//phgOptionsRegisterNoArg("use_symetric", "Use symetric Mat and PC", &_p->use_symetric);
    phgOptionsRegisterNoArg("extern_force", "Extern force", &_p->extern_force); 

    /* Solver options */
    phgOptionsRegisterString("Stokes_opts", "Solver Stokes options", &_p->Stokes_opts);
    phgOptionsRegisterString("fo_opts", "First order options", &_p->fo_opts);
    phgOptionsRegisterString("F_opts", "Solver F options", &_p->F_opts);
    phgOptionsRegisterString("Fu_opts", "Solver Fu options", &_p->Fu_opts);
    phgOptionsRegisterString("Ap_opts", "Solver Ap options", &_p->Ap_opts);
    phgOptionsRegisterString("Qp_opts", "Solver Qp options", &_p->Qp_opts);
    phgOptionsRegisterString("T_opts", "temperature solver options", &_p->T_opts);
    phgOptionsRegisterString("proj_opts", "projection options", &_p->proj_opts);


    phgOptionsPreset("-mesh_file \"../test/cube.dat\" "
		     "-grid_coord_unit 1000 "
		     "+reduce_mesh "
		     "-utype P2 -ptype P1 " 
		     "-T_type P2 "
		     "-isop_type P1 "
		     "-nu 1.0 "
		     "-dt 1e-5 "
		     "-time_end 1.0 "
		     "-max_time_step 10 "
		     "-eps_diagP 0. "
		     "-theta 0.5 "
		     "-mem_max 2000 "
		     "-step_span 1 "
		     "-pre_refines 0 "
		     "+enclosed_flow "
		     "-non_linear "
		     "-non_sub_tol 1e-3 "
		     "-periodicity 0 "
		     "+resume "
		     "+record "
		     "+pin_node "
		     "+curved_bdry "
		     "-n_bdry_layer 2 "
		     "-use_PCD "
		     "-fp_scale 1 "
		     "+use_moc "
		     "+use_Fu "
		     "-implicit_convect "
		     "+extern_force "
		     "-max_non_step 10 " 
		     "-newton_start 100 " 
		     "-max_non_step0 -1 " 
		     "-newton_start0 -1 "
		     "-stab_scheme -1 "
		     "-stab_remove_static -1 "
#if USE_SLIDING_BC    
		     "-sliding_bdry_scheme 0 "
		     "-sliding_penalty 1 "
#endif
		     "-constant_A 1e-16"		     
		     );		     


    load_case_options();
    
    return new_params;
}

void phgDofSetName(DOF *dof, const char *name)
{
    if (dof->name != NULL) {
	phgFree(dof->name);
	dof->name = NULL;
    }

    if (name != NULL) {
	dof->name = phgAlloc(strlen(name) + 1);
	strcpy(dof->name, name);
    }

    return;
}

void
func_u_xyz(FLOAT x, FLOAT y, FLOAT z, FLOAT *v)
{
    *(v++) = x;
    *(v++) = y;
    *(v++) = z;
    *(v++) = 0;

    *(v++) = x;
    *(v++) = y;
    *(v++) = z;
    *(v++) = 1;

    *(v++) = x;
    *(v++) = y;
    *(v++) = z;
    *(v++) = 2;
}

void
func_p_xyz(FLOAT x, FLOAT y, FLOAT z, FLOAT *v)
{
    *(v++) = x;
    *(v++) = y;
    *(v++) = z;
    *(v++) = 3;
}


/***********************************/
/* Navier-stokes equations:	   */
/* Set initial values at time = 0. */
/***********************************/
NSSolver *
phgNSCreate(GRID *g, LAYERED_MESH *gL, NSParams *ns_params0)
{
    DOF **u , **p, **T, **gradu;
    DOF_TYPE *utype, *ptype;

    NSSolver *ns = (NSSolver *) phgCalloc(1, sizeof(*ns));    
    ns->g = g;
    ns->gL = gL;
    ns->ns_params = ns_params0;

    iceSetUserFunc(ns);		/* Setup user funcs */

    /* set utype and ptype */
    {
	char s[128];
	phgOptionsPush();

	sprintf(s, "-dof_type %s", _nsp->utype_name);
	phgOptionsSetOptions(s);
	_nsp->utype = DOF_DEFAULT;

	sprintf(s, "-dof_type %s", _nsp->ptype_name);
	phgOptionsSetOptions(s);
	_nsp->ptype = DOF_DEFAULT;

	sprintf(s, "-dof_type %s", _nsp->T_type_name);
	phgOptionsSetOptions(s);
	_nsp->T_type = DOF_DEFAULT;

	phgOptionsPop();

	assert(_nsp->utype->fe_space == FE_H1);
	       // && _nsp->ptype->fe_space == FE_H1 
	       // && _nsp->utype->order >= _nsp->ptype->order);
    }

   
    /* set previous utype and ptype */
    if (_nsp->utype_prv_name != NULL) {
	char s[128];
	phgOptionsPush();

	sprintf(s, "-dof_type %s", _nsp->utype_prv_name);
	phgOptionsSetOptions(s);
	_nsp->utype_prv = DOF_DEFAULT;

	sprintf(s, "-dof_type %s", _nsp->ptype_prv_name);
	phgOptionsSetOptions(s);
	_nsp->ptype_prv = DOF_DEFAULT;

	sprintf(s, "-dof_type %s", _nsp->T_type_prv_name);
	phgOptionsSetOptions(s);
	_nsp->T_type_prv = DOF_DEFAULT;

	phgOptionsPop();

	assert(_nsp->utype_prv->fe_space == FE_H1);
	       // && _nsp->ptype->fe_space == FE_H1 
	       // && _nsp->utype->order >= _nsp->ptype->order);
    } else {
	_nsp->utype_prv = _nsp->utype;
	_nsp->ptype_prv = _nsp->ptype;
	_nsp->T_type_prv = _nsp->T_type;
    }


    utype = _nsp->utype;
    ptype = _nsp->ptype;

    u = ns->u = ns->u_queue + 1;
    p = ns->p = ns->p_queue + 1;
    T = ns->T = ns->T_queue + 1;
    gradu = ns->gradu = ns->gradu_queue + 1;
    
    ns->f = phgDofNew(g, DOF_ANALYTIC, 3, "f_u", func_f);
    ns->gn[0] = phgDofNew(g, DOF_ANALYTIC, 3, "gxbc", func_g1);
    ns->gn[1] = phgDofNew(g, DOF_ANALYTIC, 3, "gybc", func_g2);
    ns->gn[2] = phgDofNew(g, DOF_ANALYTIC, 3, "gzbc", func_g3);

    /* Set initial value */
    u[1] = phgDofNew(g, utype, Dim, "u_{n+1}", DofInterpolation);
    p[1] = phgDofNew(g, ptype, 1, "p_{n+1}", DofInterpolation);
    T[1] = phgDofNew(g, _nsp->T_type, 1, "T_{n+1}", DofInterpolation);
    set_boundary_mask(ns);	/* Problem depenedent */

#if USE_SLIDING_BC    
    /* Init surface bases */
    ns->surf_bas = get_surface_bases_3d(ns, utype);
    ns->surf_bas_2d = get_surface_bases_2d(ns, utype);
#else
    ns->surf_bas = NULL;
    ns->surf_bas_2d = NULL:
#endif    
    
    if (_nsp->pin_node) {
	phgNSPinNode(ns);  /* pin node before set data by func */
    }

#if STEADY_STATE
    /* Set initial values only at boundary.
     * For enclosed flow, pinned node has not been decided,
     *   so it is useless to set pinned node value here. */
    phgDofSetBdryDataByFunction(u[1], func_u, SETFLOW);
#else
    /* Set initial values in whole domain */
    phgDofSetDataByFunction(u[1], func_u);
    phgDofSetDataByFunction(p[1], func_p);
    phgDofSetDataByFunction(T[1], func_T);
    PERIODIC_SYNC(u[1]);
    PERIODIC_SYNC(p[1]);
    PERIODIC_SYNC(T[1]);
#endif /* STEADY_STATE */
    gradu[1] = phgDofGradient(u[1], NULL, NULL, "gradu_{n+1}");
    phgDofSetDataByFunction(T[1], func_T);
    PERIODIC_SYNC(T[1]);


    u[0] = NULL;
    p[0] = NULL;
    T[0] = NULL;
    gradu[0] = NULL;

    u[-1] = NULL;
    p[-1] = NULL;
    T[-1] = NULL;
    gradu[-1] = NULL;

    /* Init height change varible. */
    if (ns->dH == NULL) {
	DOF *mapping = GetGeomMapping(g);
	ns->dH = phgDofNew(g, mapping->type, 1, "DH", DofNoAction);
	ns->dHt = phgDofNew(g, mapping->type, 1, "DHt", DofNoAction);
	phgDofSetDirichletBoundaryMask(ns->dH, 0);
	phgDofSetDirichletBoundaryMask(ns->dHt, 0);
    }

    /* Time and dt */
    ns->time = ns->time_queue + 1;
    ns->dt = ns->dt_queue + 1;
    Bzero(ns->time_queue);
    Bzero(ns->dt_queue);

    /* Init nonlinear iteration DOF */
#if STEADY_STATE || TIME_DEP_NON
    ns->du = phgDofCopy(u[1], NULL, NULL, "du");
    ns->dp = phgDofCopy(p[1], NULL, NULL, "dp");
    ns->dT = phgDofCopy(T[1], NULL, NULL, "dT");
    phgDofSetDataByValue(ns->du, 0.);
    phgDofSetDataByValue(ns->dp, 0.);
    phgDofSetDataByValue(ns->dT, 0.);

#  if USE_SLIDING_BC    
    if (ns_params->sliding_bdry_scheme == 1) {/* surf L2 constraint */
	ns->Un = phgDofNew(g, ns_params->utype, 1, "Un", DofNoAction);
	ns->dUn = phgDofNew(g, ns_params->utype, 1, "dUn", DofNoAction);
    }
#  endif
    
#else
    TIME_DEP_LINEAR_ENTRY; /* Unavailable */
#endif /* STEADY_STATE || TIME_DEP_NON*/
    ns->non_res = 1e10;
 
    /* Coordinates */
    ns->coord = phgDofNew(g, DOF_P1, Dim, "Coord P1", func_xyz_);
    ns->coord_sigma
	= phgDofNew(g, DOF_P1, 1, "Coord sigma P1", DofNoAction);

    /* Beta */
    ns->beta = phgDofNew(g, ns_params->utype, 1, "beta", func_beta);
	if (ns_params->output_beta)
		phgExportVTK(ns->g, "beta.vtk", ns->beta, NULL);
 
    /* Div0 */
    ns->div0 = phgDofNew(g, DOF_P0, 1, "div0", DofNoAction);

    /* Viscosity */
    ns->nu = phgDofNew(g, ns_params->utype, 1, "nu", DofNoAction);


#if STUDY_ILU    
    /* Dof coord */
    {
	FILE *fp;
	DOF *coord_u, *coord_p, *dofs[2];
	MAP *map_;
	VEC *vec_;
	INT i, j, k;

	coord_u = phgDofNew(g, ns->du->type, Dim * (Dim+1),
			    "Coord u", func_u_xyz);
	coord_p = phgDofNew(g, ns->dp->type, Dim+1,
			    "Coord p", func_p_xyz);

	dofs[0]= coord_u;
	dofs[1]= coord_p;
	map_ = phgMapCreate(coord_u, coord_p, NULL);
	vec_ = phgMapCreateVec(map_, 1);

	
	
#if ICE_BENCH_TEST && TEST_CASE != ICE_BENCH_E
	/* Handle periodic:
	 *
	 * phgMapDofToLocalData:
	 *    multi dof data ==--> vec data,
	 *    when dof data diffs, (e.g. coord),
	 *    vec data is not uniq.
	 * */
	FLOAT *vec_dat = vec_->data;
	SIMPLEX *e;
	ForAllElements(g, e) {
	    const FLOAT *vcrd[NVert];
	    FLOAT *dof_dat;
	    INT dof_ofs, dim;

#define SET_VERT_DATA(DOF, NO)				\
	    dof_dat = DofVertexData(DOF, e->verts[i]);	\
	    dof_ofs = dof_dat - DOF->data;		\
	    dim = DOF->count_vert;			\
	    for (k = 0; k < dim; k++) {			\
		j = phgMapD2L(map_, NO, dof_ofs + k);	\
		vec_dat[j] = dof_dat[k];		\
	    }

#define SET_EDGE_DATA(DOF, NO)				\
	    dof_dat = DofEdgeData(DOF, e->edges[i]);	\
	    dof_ofs = dof_dat - DOF->data;		\
	    dim = DOF->count_edge;			\
	    for (k = 0; k < dim; k++) {			\
		j = phgMapD2L(map_, NO, dof_ofs + k);	\
		vec_dat[j] = dof_dat[k];		\
	    }
	    
	    /* Vert */
	    for (i = 0; i < NVert; i++) {
		vcrd[i] = g->verts[e->verts[i]];

		if (vcrd[i][0] > _Length_ - 1e-6
		    || vcrd[i][1] > _Length_ - 1e-6)
		    continue;

		SET_VERT_DATA(coord_u, 0);
		SET_VERT_DATA(coord_p, 1);
	    }

	    for (i = 0; i < NEdge; i++) {
		int v0, v1;
		FLOAT ecrd[Dim];
		GetEdgeVertices(e, i, v0, v1);

		for (k = 0; k < Dim; k++)
		    ecrd[k] = .5 * (vcrd[v0][k] + vcrd[v1][k]);

		if (ecrd[0] > _Length_ - 1e-6
		    || ecrd[1] > _Length_ - 1e-6)
		    continue;
		
		SET_EDGE_DATA(coord_u, 0);
		SET_EDGE_DATA(coord_p, 1);
	    }
	}
#else
	phgMapDofToLocalData(map_, 2, dofs, vec_->data);
#endif

	
	/* fp = fopen("coord_u.dat", "w"); */
	/* fclose(fp); */
	phgVecDumpMATLAB(vec_, "xyz", "dof_xyz_.m");
	
	phgVecDestroy(&vec_);
	phgMapDestroy(&map_);
	phgDofFree(&coord_u);
	phgDofFree(&coord_p);
    }
#endif	/* STUDY_ILU */

    /*
     *
     * Init coord sigma
     *
     * */
    {
	const FLOAT *ratio = gL->layer_ratio;
	int ii, i, j;
	for (ii = 0; ii < gL->nvert_bot; ii++) {
	    i = gL->vert_bot_Lidx[ii];
	    assert(gL->vert_local_lists[i] != NULL);

	    int nv = gL->vert_local_lists[i][0];
	    int *iL = &gL->vert_local_lists[i][1];
	    int iv = gL->vert_L2S[i];
	    assert(nv > 0);

	    assert(gL->max_nlayer + 1 == nv);
	    for (j = 0; j < nv; j++) {
		*DofVertexData(ns->coord_sigma, iL[j]) = ratio[j];
	    }
	}
    }

    
    /* Init Other things */
    if (ns->init_func)
	ns->init_func(ns);

    //phgFinalize();
    //exit(0);
    
    return ns;
}

void phgNSFinalize(NSSolver **ns_ptr)
{
    NSSolver *ns = *ns_ptr;
    DOF **u = ns->u, **p = ns->p, **T = ns->T,
	**gradu = ns->gradu;

#if STEADY_STATE || TIME_DEP_NON 
    phgDofFree(&ns->du);
    phgDofFree(&ns->dp);
    phgDofFree(&ns->dT);

#  if USE_SLIDING_BC    
    if (ns_params->sliding_bdry_scheme == 1) {/* surf L2 constraint */
	phgDofFree(&ns->Un);
	phgDofFree(&ns->dUn);
    }
#  endif    
#else
    TIME_DEP_LINEAR_ENTRY; /* Unavailable */
#endif /* STEADY_STATE || TIME_DEP_NON */
    
    if (u[-1] != NULL) {
	phgDofFree(&u[-1]);
	phgDofFree(&p[-1]);
	phgDofFree(&T[-1]);
	phgDofFree(&gradu[-1]);
    }

    if (u[0] != NULL) {
	phgDofFree(&u[0]);
	phgDofFree(&p[0]);
	phgDofFree(&T[0]); 
	phgDofFree(&gradu[0]);
    }

    phgDofFree(&u[1]);
    phgDofFree(&p[1]);
    phgDofFree(&T[1]);
    phgDofFree(&gradu[1]);

    phgDofFree(&ns->f);
    phgDofFree(&ns->gn[0]);
    phgDofFree(&ns->gn[1]);
    phgDofFree(&ns->gn[2]);
    phgDofFree(&ns->dH);
    phgDofFree(&ns->dHt);
    phgDofFree(&ns->coord);
    phgDofFree(&ns->depth_P1);
    phgDofFree(&ns->depth_T);

    if (ns->Gradu) phgDofFree(&ns->Gradu);
    if (ns->div0) phgDofFree(&ns->div0);
    if (ns->beta) phgDofFree(&ns->beta);
    if (ns->p_static) phgDofFree(&ns->p_static);
    if (ns->surf_bas) {
	phgDofFree(&ns->surf_bas->dof);
	phgFree(ns->surf_bas);
    }

    
    phgFree(ns);
    *ns_ptr = NULL;

    return;
}



/*******************************************************************************/
/* Navier-stokes equations:						       */
/* Time advance: shift solutions at different time.			       */
/* Note: 								       */
/* 1. this routine MUST be called AFTER time update,			       */
/*     since DOF f is evauled using t_{n+1}.				       */
/* 2. this routine MUST be called AFTER grid change (refine/coercen, balance), */
/*     since gradu is not interpolated during grid change, instead,	       */
/*     they are recaculated using the interpolated u.			       */
/*******************************************************************************/
void
phgNSTimeAdvance(NSSolver *ns, FLOAT time, int tstep)
{
    GRID *g = ns->g;
    DOF **u = ns->u, **p = ns->p, **T = ns->T,
	**gradu = ns->gradu;

    elapsed_time(g, FALSE, 0.);	/* reset timer */

    ns->tstep = tstep;
    ns->time[0] = ns->time[1];
    ns->time[1] = time;
    /* phgPrintf("Time advance: "); */

    /* --------------
     * 1. Last level:
     *    u_{n-1} -> discard,
     *    u_{n}   -> u{n-1}.
     *  */
    phgDofFree(&u[-1]);
    phgDofFree(&p[-1]);
    phgDofFree(&T[-1]);
    if (u[0] != NULL) {
	/* phgPrintf("Free u_{n-1}; "); */

	assert (!strcmp(u[0]->name, "u_{n}"));
	assert (!strcmp(p[0]->name, "p_{n}"));
	assert (!strcmp(T[0]->name, "T_{n}"));
	assert (!strcmp(gradu[0]->name, "gradu_{n}"));

	u[-1] = u[0];
	p[-1] = p[0];
	T[-1] = T[0];
	/* grad u_{n-1} unused */
	phgDofFree(&gradu[0]);

	phgDofSetName(u[-1], "u_{n-1}");
	phgDofSetName(p[-1], "p_{n-1}");
	phgDofSetName(T[-1], "T_{n-1}");

    } 
    
    /* --------------
     * 2. mid level:
     *    copy its upper level's value sequentially from 2nd last level to 2nd top level.
     * */

    assert (!strcmp(u[1]->name, "u_{n+1}"));
    assert (!strcmp(p[1]->name, "p_{n+1}"));
    assert (!strcmp(T[1]->name, "T_{n+1}"));
    u[0] = u[1];
    p[0] = p[1];
    T[0] = T[1];
    gradu[0] = gradu[1];
    phgDofSetName(u[0], "u_{n}");
    phgDofSetName(p[0], "p_{n}");
    phgDofSetName(T[0], "T_{n}");
    phgDofSetName(gradu[0], "gradu_{n}");


    /* --------------
     * 3. top level:
     *    created by copying 2nd top level.
     * */
    phgPrintf("      new u");
    u[1] = phgDofCopy(u[0], NULL, NULL, "u_{n+1}");
    p[1] = phgDofCopy(p[0], NULL, NULL, "p_{n+1}");
    T[1] = phgDofCopy(T[0], NULL, NULL, "T_{n+1}");

    /* Set Dirich B.C */
    phgPrintf("   Set init flow at bdry, ");
    phgDofSetBdryDataByFunction(u[1], func_u, SETFLOW); /* No inflow */


    /* Pin_node value */
    if (_nsp->pin_node) {
	phgNSPinNode(ns);  /* pin node before set data by func */
	
	adjust_time(- (1. - _nsp->Theta) * ns->dt[0]);
	phgDofSetBdryDataByFunction(p[1], func_p, BC_PIN_NODE);
	restore_time();
    }
    PERIODIC_SYNC(u[1]);
    PERIODIC_SYNC(p[1]);
    PERIODIC_SYNC(T[1]);
    elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));


    /* Recompute grad u_{n+1},
     * grad T_{n+1} not include. */
    phgPrintf("      gradu");
#if STEADY_STATE
    gradu[1] = phgDofCopy(gradu[0], NULL, NULL, "gradu_{n+1}");
#elif TIME_DEP_NON
    gradu[1] = phgDofGradient(u[1], NULL, NULL, "gradu_{n+1}");
#else
    TIME_DEP_LINEAR_ENTRY; /* Unavailable */
#endif
    elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));
    ns->non_res = 1e10;

    /* dt change not here*/

    //phgExportVTK(g, "ns_init.vtk", u[1], p[1], NULL); /* Output to check grid and DOFs. */


    /* Compute hydro static pressure */
    get_p_static(ns);

    /* Update floating mask */
    update_floating(ns);
	
    return;
}


/* Pin a node to remove the singularity of the system.
 *
 * Note:
 *   1. g->types_vert[pinned_node] will be changed when refine, coarsen, balance grid,
 *      so after calling this routine, grid should be fixed until linear solver finishes.
 *
 *   2. the pinned vertex need to be local to avoid global change,
 *      TODO: fixed this for the case when all verts are remote.
 *      
 * */
INT 
phgNSPinNode(NSSolver *ns)
{
    GRID *g = ns->g;
    SIMPLEX *e;
    int i, j, id, pinned_node_id = -1;

    if (!_nsp->enclosed_flow)
    	phgError(1, "Pining node is only needed for enclosed flow.\n");

    phgUpdateBoundaryTypes(g); 	/* clean pinned node */
    MPI_Barrier(g->comm);
    pinned_node_id = -1;

#if PIN_AT_ROOT
    if (g->rank == 0) {
	ForAllElements(g, e) {
	    for (j = 0; j < NVert; j++) {
		if (pinned_node_id == -1) {
		    FLOAT vp;
		    i = e->verts[j];
		    if((g->types_vert[i] & REMOTE) || !(g->types_vert[i] & BDRY_MASK))
		    	continue;

		    assert(g->verts != NULL && g->verts[j] != NULL);
		    g->types_vert[i] |= PINNEDNODE;
		    pinned_node_id = i;
		    printf("   Pin a node at[%d] in e[%d]: (%f, %f, %f), at rank: %d\n",
			   i, e->index,
			   g->verts[i][0], g->verts[i][1], g->verts[i][2], 
			   g->rank);
		    func_p(g->verts[i][0], g->verts[i][1], g->verts[i][2], &vp);
		    printf("      value: %16.8E\n", vp);
		}
	    }
	}
	if (pinned_node_id == -1)
	    phgError(-1, "All verts in rank 0 is either remote or interior, "
		     "leaves no one to pin!\n");
    }
#else
    /* pick one node in each proc as candidates to be pinned */
    ForAllElements(g, e) {
	for (j = 0; j < NVert; j++) {
	    if (pinned_node_id == -1) {
		//FLOAT vp;
		i = e->verts[j];
		if(!(g->types_vert[i] & BDRY_MASK))
		    continue;
		pinned_node_id = GlobalVertex(g, i);

		assert(g->verts != NULL && g->verts[i] != NULL);
#if 0
		/* check */
		printf("   [%2d] possible pin node at vert local[%d] global[%d]\n",
		       g->rank, i, pinned_node_id);
		printf("       in e[%d]: (%f, %f, %f)\n",
			e->index, g->verts[i][0], g->verts[i][1], g->verts[i][2]);
		func_p(g->verts[i][0], g->verts[i][1], g->verts[i][2], &vp);
		printf("       value: %16.8E\n", vp);
#endif
	    }
	}
    }

    /* Pin the node with MAX global index */
    id = pinned_node_id;
    MPI_Allreduce(&id, &pinned_node_id, 1, MPI_INT, MPI_MAX, g->comm);
    phgPrintf("   Pin node at vert global[%d]\n", pinned_node_id);
    if (pinned_node_id == -1)
	phgError(-1, "Pin verts err!\n");

    /* Set pinned node type */
    for (i = 0; i < g->nvert; i++) {
	if (!(g->types_vert[i] & BDRY_MASK))
	    continue;
	if (GlobalVertex(g, i) == pinned_node_id) {
	    g->types_vert[i] |= BC_PIN_NODE;
	    printf("   [%2d] pin node at vert[%d]\n",
		    g->rank, i);
	    break;
	}
    }
    MPI_Barrier(g->comm);
#endif	/* PIN_AT_ROOT */

    /* pinned_node_id = -1 on non-root rank */
    ns->pinned_node_id = pinned_node_id; 


    /*
     * Set dof value after pin node.
     * */
    /* if (_nsp->pin_node) { */
    /* 	adjust_time(- (1. - _nsp->Theta) * ns->dt[0]); */
    /* 	phgDofSetBdryDataByFunction(ns->p[1], func_p, PINNEDNODE); */
    /* 	restore_time(); */
    /* } */

    MPI_Barrier(g->comm);
    return ns->pinned_node_id;
}


/******************/
/* Init NS solver */
/******************/
void phgNSInitSolverU(NSSolver *ns)
{
    GRID *g = ns->g;
    MAP *Vmap, *Pmap;
    //MAT *pmat[4];
    MAT *pmat[16];		/* Un penalty */
    DOF **u = ns->u; 
    FLOAT *dt = ns->dt;
    INT verb;
#if USE_SLIDING_BC
    MAP *UNmap;
    BOOLEAN use_un_constraint = (ns_params->sliding_bdry_scheme == 1) ? TRUE : FALSE;
#endif
    MAP *bubmap;
    
    
    Unused(dt);  
    Unused(verb);  
    /* dof copy */
    ns->u_shape = phgDofCopy(ns->du, NULL, NULL, "u shape");
    ns->p_shape = phgDofCopy(ns->dp, NULL, NULL, "p shape");
  
    /* dof map */
    ns->Vmap = phgMapCreate(ns->u_shape, NULL);
    ns->Pmap = phgMapCreate(ns->p_shape, NULL);

    phgPrintf("   solver_u size [%d, %d]\n", ns->Vmap->nglobal, ns->Vmap->bdry_nglobal);
    phgPrintf("   solver_p size [%d, %d]\n", ns->Pmap->nglobal, ns->Pmap->bdry_nglobal);

    Vmap = ns->Vmap;
    Pmap = ns->Pmap;

    /* matrices */
    ns->matF = phgMapCreateMat(Vmap, Vmap);
    ns->matBt = phgMapCreateMat(Vmap, Pmap);
    ns->matB = phgMapCreateMat(Pmap, Vmap);
    ns->matC = phgMapCreateMat(Pmap, Pmap);


    ns->matF->handle_bdry_eqns = MAT_HANDLE_BDRY_EQNS;
    if (ns->matC != NULL)
	ns->matC->handle_bdry_eqns = MAT_HANDLE_BDRY_EQNS;

    ns->matB->handle_bdry_eqns = FALSE;
    ns->matBt->handle_bdry_eqns = FALSE;


    
#if USE_SLIDING_BC    
    if (use_un_constraint) {
	ns->un_shape = phgDofCopy(ns->dUn, NULL, NULL, "un shape");

	ns->UNmap = phgMapCreate(ns->un_shape, NULL);
	UNmap = ns->UNmap;
    
	ns->matUn = phgMapCreateMat(UNmap, Vmap);
	ns->matUnT = phgMapCreateMat(Vmap, UNmap);
	ns->matUnD = phgMapCreateMat(UNmap, UNmap);

	ns->matUn->handle_bdry_eqns = FALSE;
	ns->matUnT->handle_bdry_eqns = FALSE;
	ns->matUnD->handle_bdry_eqns = FALSE;
    }
#endif

    
#if USE_SLIDING_BC    
    if (use_un_constraint) {
	phgPrintf("--- Use 3x3 mat ---\n");
	    
	/*          | F   Bt  UnT|
	 *  matNS = | B   C      |
	 *          | Un         |
	 *  */
	pmat[0] = ns->matF;
	pmat[1] = ns->matBt;
	pmat[2] = ns->matUnT;

	pmat[3] = ns->matB;
	pmat[4] = ns->matC;
	pmat[5] = NULL;

	pmat[6] = ns->matUn;
	pmat[7] = NULL;
	pmat[8] = ns->matUnD;

	ns->matNS = phgMatCreateBlockMatrix(g->comm, 3, 3, pmat, NULL, NULL);
    }
    else {
	phgPrintf("--- Use 2x2 mat ---\n");
	/*          | F  Bt|
	 *  matNS = |      |
	 *          | B  C |
	 *  */
	pmat[0] = ns->matF;
	pmat[1] = ns->matBt;
	pmat[2] = ns->matB;
	pmat[3] = ns->matC;

	ns->matNS = phgMatCreateBlockMatrix(g->comm, 2, 2, pmat, NULL, NULL);
    }
#else    
    /*          | F  Bt|
     *  matNS = |      |
     *          | B  C |
     *  */
    pmat[0] = ns->matF;
    pmat[1] = ns->matBt;
    pmat[2] = ns->matB;
    pmat[3] = ns->matC;

    ns->matNS = phgMatCreateBlockMatrix(g->comm, 2, 2, pmat, NULL, NULL);
#endif    


    ns->matNS->mv_data = phgAlloc(sizeof(*ns->matNS->mv_data));
    ns->matNS->mv_data[0] = (void *) ns;
    ns->matNS->handle_bdry_eqns = MAT_HANDLE_BDRY_EQNS;

    /* solver_u */
    /* Note: can't use phgMat2Solver here because build_rhs()
     * requires solver->rhs->map to map U into vector indices */
    phgOptionsPush();
    phgOptionsSetOptions(_nsp->Stokes_opts);

#if USE_SLIDING_BC    
    if (use_un_constraint) 
	ns->solver_u = phgSolverCreate(SOLVER_DEFAULT, ns->u_shape, ns->p_shape, ns->un_shape,
				       NULL);
    else
	ns->solver_u = phgSolverCreate(SOLVER_DEFAULT, ns->u_shape, ns->p_shape,
				       NULL);
#else    
	ns->solver_u = phgSolverCreate(SOLVER_DEFAULT, ns->u_shape, ns->p_shape,
				       NULL);
#endif    

    phgOptionsPop();

    {
	phgInfo(0, "nlocal1: %d\n", ns->solver_u->rhs->map->nlocal);
	phgInfo(0, "nlocals: %d\n", ns->solver_u->rhs->map->localsize);
    }

    phgMatDestroy(&ns->solver_u->mat);
    ns->solver_u->mat = ns->matNS;
    ns->solver_u->rhs->mat = ns->solver_u->mat;
    if (_nsp->non_linear)
	ns->solver_u->rtol = _nsp->non_sub_tol;

#if STEADY_STATE || TIME_DEP_NON
    phgDofSetDataByValue(ns->du, 0.);
    phgDofSetDataByValue(ns->dp, 0.);
    
#   if USE_SLIDING_BC    
    if (use_un_constraint) 
	phgDofSetDataByValue(ns->dUn, 0.);
#   endif    
    
    phgDofFree(&ns->wind);
    ns->wind = phgDofCopy(u[1], NULL, NULL, "wind");
#else
    TIME_DEP_LINEAR; /* Unavailable */
#endif /* STEADY_STATE || TIME_DEP_NON */

#if GIVEN_WIND 
    phgDofSetDataByFunction(ns->wind, func_wind);
#endif /* GIVEN_WIND */

    /* Set boundary types of ice-sheet according to temprature. */
    if (ns->mark_bdry_temp)
	ns->mark_bdry_temp(ns);

    return;
}

void phgNSReInitSolverU(NSSolver *ns)
{
    DOF **u = ns->u; 

    if (_nsp->non_linear)
	ns->solver_u->rtol = _nsp->non_sub_tol;

#if STEADY_STATE || TIME_DEP_NON
    phgDofSetDataByValue(ns->du, 0.);
    phgDofSetDataByValue(ns->dp, 0.);
    phgDofFree(&ns->wind);
    ns->wind = phgDofCopy(u[1], NULL, NULL, "wind");
#else
    TIME_DEP_LINEAR; /* Unavailable */
#endif /* STEADY_STATE || TIME_DEP_NON */

#if GIVEN_WIND 
    phgDofSetDataByFunction(ns->wind, func_wind);
#endif /* GIVEN_WIND */

    return;
}


/**********************/
/* Destroy NS solver  */
/**********************/
void phgNSDestroySolverU(NSSolver *ns)
{
    phgInfo(2, "   Destroy solver U\n");
    phgDofFree(&ns->wind);
    ns->wind = NULL;

    phgMatDestroy(&ns->matF);
    phgMatDestroy(&ns->matB);
    phgMatDestroy(&ns->matBt);
    phgMatDestroy(&ns->matC);


#if USE_SLIDING_BC    
    if (ns_params->sliding_bdry_scheme == 1) {
	phgMatDestroy(&ns->matUn);
	phgMatDestroy(&ns->matUnT);
	phgMatDestroy(&ns->matUnD);
    }
#endif
    
    phgSolverDestroy(&ns->solver_u);
    ns->solver_u = NULL;
    phgMapDestroy(&ns->Vmap);
    phgMapDestroy(&ns->Pmap);
    phgDofFree(&ns->u_shape);
    phgDofFree(&ns->p_shape);

#if USE_SLIDING_BC    
    if (ns_params->sliding_bdry_scheme == 1) {
	phgMapDestroy(&ns->UNmap);
	phgDofFree(&ns->un_shape);
    }
#endif
}










