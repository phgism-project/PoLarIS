
#include "ins.h"
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#if USE_PETSC
#  include <petscksp.h>
#endif 




/***************/
/* GLOBAL vars */
/***************/
NSParams *ns_params = NULL;
FLOAT eps_height = 1e-3;	/* 1m */

char vtk_file[1000];

/* Functions to check DOF */
/* Velocity */
static
void func_u0(FLOAT x, FLOAT y, FLOAT z, FLOAT *u) {
    /* const FLOAT x0 = 750e3 / LEN_SCALING; */
    /* const FLOAT y0 = 750e3 / LEN_SCALING; */

    FLOAT d;
    //u[0] = .5 * z * z ;
    u[0] = 0.;
    u[1] = 0.;
    u[2] = 0.;
    
    /* x -= x0; */
    /* y -= y0; */
    /* z /= 3.2; */

    /* u[0] = x / 5  * (1-pow(1-z, 4)); */
    /* u[1] = y / 10 * (1-pow(1-z, 4)); */
    /* u[2] = - z; */

    /* u[0] *= 1+sin(_t_*0.1); */
    /* u[1] *= 1+sin(_t_*0.3); */
    /* u[2] *= 1+sin(_t_*0.2); */

    /* u[0] = 100; */
    /* u[1] = 0; */
    /* u[2] = -1; */

    return;
}

/* Pressure */
static
void func_p0(FLOAT x, FLOAT y, FLOAT z, FLOAT *p) {
    /* *p =  z * (1+sin(_t_*0.1)); */
    *p = (x > 40) ? -z : 0;
    return;
}


/* Accumulation
 *  unit: km/a
 * */
static
void func_Q(FLOAT x, FLOAT y, FLOAT z, FLOAT *q) {
    func_q(x, y, z, q);
    *q /= LEN_SCALING;
    return;
}


/****************/
/* Main program */
/****************/
int
main(int argc, char *argv[])
{
    GRID *g;
    SIMPLEX *e;
    DOF **u, **p, **T, **gradu, 
	*u_exact, *p_exact, *gradu_exact, *T_exact,
	*dp1 = NULL, *p1 = NULL,
	*eu, *ep, *egradu, *ediv, *eT,
	*dH = NULL;
    FLOAT Time, *dt, res, non_du, non_dp, non_dT;
    INT tstep = 0, nelem;
    char mesh_file[100], hostname[256],
	data_file[100], data_u[100], data_p[100], data_T[100], data_Crd[100];
    size_t mem, mem_peak;
    int verb;
    double tt[3], tt1[3];
    /* ---------- NS ---------- */
    NSSolver *ns = NULL;
    SURF_BAS *surf_bas = NULL;
    LAYERED_MESH *gL = NULL;
    /* MG_BLOCK_DOFS *bk = NULL; */
    BOOLEAN debug = FALSE;
    NSPCSolver *pc = NULL;

    /* ================================================================================
     *
     *         Initialize Grid & parameters
     *
     * ================================================================================
     */
    phgOptionsRegisterNoArg("debug", "mpi debug", &debug);

    /* Global (static) options */
    Unused(verb);  
    ns_params = phgParametersCreate();     
    phgInit(&argc, &argv);
    phgOptionsShowUsed();

    if (debug) {
        int _i_ = 0;
        unsigned int t = 5;
        int pid = getpid();

        gethostname(hostname, sizeof(hostname));
        printf("#### Lengweee debug "
               "PID %d on %s ready for attach\n", getpid(), hostname);
        fflush(stdout);
        while (0 == _i_) {
            MPI_Barrier(MPI_COMM_WORLD);
            printf("%d after bar\n", pid);
            fflush(stdout);
            for (t=0; t<50; t++)
                sleep(1);
            printf("%d after sleep\n", pid);
            fflush(stdout);
        }
        printf("### PID %d, ready!\n", getpid());
    }


    install_error_log();


    g = phgNewGrid(-1);
    phgSetPeriodicity(g, ns_params->periodicity);

    phgImportSetBdryMapFunc(my_bc_map);
    if (ns_params->resume) {
	phgResumeStage(g, &Time, &tstep, mesh_file, data_file);
	phgPrintf("================================\n\n");
	phgPrintf("* RESUME from time:%E, tstep:%d\n", Time, tstep);
	phgPrintf("*             mesh:%s\n", mesh_file);
	phgPrintf("*             data:%s\n", data_file);
	phgPrintf("================================\n");

	if (!phgImport(g, mesh_file, FALSE))
	    phgError(1, "can't read file \"%s\".\n", ns_params->fn);
    } else {
	phgPrintf("Using mesh: %s\n", ns_params->fn);
	if (!phgImport(g, ns_params->fn, FALSE))
	    phgError(1, "can't read file \"%s\".\n", ns_params->fn);
    }

    phgPrintf("Check boundary before mapping to height.\n");
    checkBdry(g);
    elapsed_time(g, FALSE, 0.);	/* reset timer */
    gethostname(hostname, sizeof(hostname));
    printf("#%5d# runing PID %5d on %s \n", phgRank, getpid(), hostname);


    NsSolver_Options();
    phgPrintf("  Pre-refine & repartition ");
    phgRefineAllElements(g, ns_params->pre_refines);

    /* Set Reynolds number */
    Time = ns_params->time_start;	/* default: 0 */
    setFuncTime(Time);
    setFlowParameter(ns_params->Re, ns_params->nu, Time);
    


    /* ================================================================================
     *
     *         build ice grid
     *
     * ================================================================================
     */
    iceInit(g, &gL);


    phgPrintf("Check boundary after mapping to height.\n");
    checkBdry(g);
    //phgExportVTK(g, "ice_domain.vtk", NULL, NULL);


    
    /* ================================================================================
     *
     *         Create INS solver
     *
     * ================================================================================
     */

    /* Note: pointers u, p, gradu, dt
     *       DIRECTLY access private member of INS solver.
     * */
    phgPrintf("  Create INS solver");
    tstep = 1;			/* time step start at 1 */
    setFuncTime(Time);          /* in file ins-bc.c: static */
    ns = phgNSCreate(g, gL, ns_params);
    ns->time[1] = Time;
    u = ns->u; 
    p = ns->p;
    T = ns->T;
    gradu = ns->gradu;
    dH = ns->dH;
    dt = ns->dt;		/* direct accses ns  */
    dt[0] = ns_params->dt0;
    //ns->bk = bk = NULL; //init_line_block(T[1], gL);   /* Use line block */
    elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));
    phgExportTecplot(g, "init_beta.plt", ns->beta, NULL);

    /* Removed */
    //build_prism_elems(g, gL);    

    /* Init height & depth */
    if (gL != NULL) {
	get_height_depth(ns);
	//phgExportVTK(g, OUTPUT_DIR "depth.vtk", ns->depth_P1, NULL);
    }

    /* Export exact solution to check flow pattern */
    if (0) {
	phgDofSetDataByFunction(u[1], func_u);
	phgDofSetDataByFunction(p[1], func_p);
	if (ns_params->ptype == DOF_G1)
	    memset(p[1]->data_elem, 0, g->nelem * sizeof(FLOAT));
	phgExportVTK(g, OUTPUT_DIR "/ins_" NS_PROBLEM "_init.vtk", u[1], p[1], T[1], NULL);
	DOF_SCALE(u[1], "init u");
	DOF_SCALE(p[1], "init p");
	phgDofDump(u[1]);
	phgDofDump(p[1]);
	MPI_Barrier(MPI_COMM_WORLD);
	//phgDofDump(u[1]);
    }

    /* Exact solutions */
    u_exact = phgDofNew(g, DOF_ANALYTIC, 3, "exact u", func_u);
    p_exact = phgDofNew(g, DOF_ANALYTIC, 1, "exact p", func_p);
    T_exact = phgDofNew(g, DOF_ANALYTIC, 1, "exact T", func_T);
    gradu_exact = phgDofNew(g, DOF_ANALYTIC, 9, "exact grad u", func_gradu);

    if (ns_params->ptype == DOF_G1) {
	p1 = phgDofNew(g, DOF_P1, 1, "pp1", DofNoAction);
	dp1 = phgDofNew(g, DOF_P1, 1, "dpp1", DofNoAction);
    }

    /* surf bases */
    surf_bas = ns->surf_bas;

#if 0
    DOF *beta = phgDofNew(g, DOF_P1, 1, "beta", func_beta);
    phgExportEnsight(g, "check", beta, NULL);
    phgFinalize();
    exit(1);
#endif



    /* ------------------------------------------------------------
     * 
     * Resume dof data.
     * 
     * ------------------------------------------------------------ */
    if (ns_params->resume) {
	FILE *fp = NULL;
	char fname[1000];
	static DOF *u_prv, *p_prv, *T_prv;
	phgResumeStage(g, &Time, &tstep, mesh_file, data_file);

	u_prv = phgDofNew(g, ns_params->utype_prv, Dim, "u_prv", DofNoAction);
	p_prv = phgDofNew(g, ns_params->ptype_prv, 1, "p_prv", DofNoAction);
	T_prv = phgDofNew(g, ns_params->T_type_prv, 1, "T_prv", DofNoAction);

	/* resume coord */
	{
	    const FLOAT *v = DofData(ns->coord);
	    int i, k;
	    sprintf(data_Crd, OUTPUT_DIR "/ins_" NS_PROBLEM "_%05d.dat.Crd", tstep - 1);
	    assert(ns->coord->type == DOF_P1);
	    load_dof_data3(g, ns->coord, data_Crd, mesh_file);
	    for (i = 0; i < g->nvert; i++) 
		for (k = 0; k < Dim; k++)
		    g->verts[i][k] = *(v++);
	    phgGeomInit_(g, TRUE);
	}

	/* resmue u_{n-1} */
	sprintf(data_file, OUTPUT_DIR "/ins_" NS_PROBLEM "_%05d.dat", tstep - 2);
	DATA_FILE_SURFIX;
	sprintf(fname, "%s.p%03d", data_u, g->rank);
	if ((fp = fopen(fname, "r")) == NULL) {
	    phgPrintf("*  u_{%d} unavailable.\n", tstep - 2);
	} 
	else {
	    fclose(fp);

	    load_dof_data3(g, u_prv, data_u, mesh_file);
	    load_dof_data3(g, p_prv, data_p, mesh_file);
	    load_dof_data3(g, T_prv, data_T, mesh_file);
 
	    phgDofCopy(u_prv, &u[1], NULL, "u_{n+1}");
	    phgDofCopy(p_prv, &p[1], NULL, "p_{n+1}");
	    phgDofCopy(T_prv, &T[1], NULL, "T_{n+1}");
	    set_boundary_mask(ns); /* set DB for n+1 */
	    phgDofCopy(u[1], &u[0], NULL, "u_{n}");
	    phgDofCopy(u[1], &p[0], NULL, "p_{n}");
	    phgDofCopy(u[1], &T[0], NULL, "T_{n}");

	    phgPrintf("   Resume u_ {%5d}[%8d]:%24.12E p_ {%5d}[%8d]:%24.12E\n", 
		      tstep - 2, DofGetDataCountGlobal(u[0]), phgDofNormL2(u[0]), 
		      tstep - 2, DofGetDataCountGlobal(p[0]), phgDofNormL2(p[0]));
	    phgPrintf("   Resume T_{%5d}[%8d]:%24.12E\n", 
		      tstep - 2, DofGetDataCountGlobal(T[0]), phgDofNormL2(T[0]));
      
	    phgDofGradient(u[0], &gradu[0], NULL, "gradu_{n}");
	    phgDofSetFunction(u[0], DofInterpolation);
	    phgDofSetFunction(p[0], DofInterpolation);
	    //phgDofSetBdryDataByFunction(u[0], func_u, SETFLOW);
	    DOF_SCALE(u[0], "resume");
	    DOF_SCALE(p[0], "resume");
	    DOF_SCALE(T[0], "resume");
	    DOF_SCALE(gradu[0], "resume");

	    elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));
	}

	/* resmue u_{n} */
	sprintf(data_file, OUTPUT_DIR "/ins_" NS_PROBLEM "_%05d.dat", tstep - 1);
	DATA_FILE_SURFIX;
	sprintf(fname, "%s.p%03d", data_u, g->rank);
	if ((fp = fopen(fname, "r")) == NULL) {
	    phgError(1, "read Dof data %s failed!\n", data_file);
	} else {
	    fclose(fp);

	    load_dof_data3(g, u_prv, data_u, mesh_file);
	    load_dof_data3(g, p_prv, data_p, mesh_file);
	    load_dof_data3(g, T_prv, data_T, mesh_file);
 
	    phgDofCopy(u_prv, &u[1], NULL, "u_{n+1}");
	    phgDofCopy(p_prv, &p[1], NULL, "p_{n+1}");
	    phgDofCopy(T_prv, &T[1], NULL, "T_{n+1}");
	    set_boundary_mask(ns);

	    phgPrintf("   Resume u_ {%5d}[%8d]:%24.12E p_ {%5d}[%8d]:%24.12E\n", 
		      tstep - 1, DofGetDataCountGlobal(u[1]), phgDofNormL2(u[1]), 
		      tstep - 1, DofGetDataCountGlobal(p[1]), phgDofNormL2(p[1]));
	    phgPrintf("   Resume T_{%5d}[%8d]:%24.12E\n", 
		      tstep - 1, DofGetDataCountGlobal(T[1]), phgDofNormL2(T[1])); 
      
	    phgDofGradient(u[1], &gradu[1], NULL, "gradu_{n+1}");
	    phgDofSetFunction(u[1], DofInterpolation);
	    phgDofSetFunction(p[1], DofInterpolation);
	    //phgDofSetBdryDataByFunction(u[1], func_u, SETFLOW);
	    DOF_SCALE(u[1], "resume");
	    DOF_SCALE(p[1], "resume");
	    DOF_SCALE(T[1], "resume");
	    DOF_SCALE(gradu[1], "resume");

	    elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));
	}

	/* Re init height & depth */
	if (gL != NULL) {
	    get_height_depth(ns);
	    build_layered_mesh_height(g, gL);
	    check_height(g, gL);
	}

	/* reconstruct last time step */
	Time -= dt[0];
	ns->time[1] = Time;
	ns->time[0] = Time;
	setFuncTime(Time);
#if 0
	phgExportEnsightT(g, OUTPUT_DIR "/ins_" NS_PROBLEM , tstep, tstep, u[1], p[1], NULL);
	phgFinalize();
	return 0;
#endif	/* debug exit */

    }	/* end of resume */


    /* Init temp field */
    if (ns_params->init_temp_type >= 0
	&& tstep == 1) {
	phgNSTempInit(ns);
	/* DOF *temp_diff = phgDofCopy(ns->T[1], NULL, NULL, "Td"); */
	/* { */
	/* 	FLOAT *vt = temp_diff->data; */
	/* 	const FLOAT *vh = ns->depth_T->data; */

	/* 	INT i, n = DofGetDataCount(temp_diff); */
	/* 	for (i = 0; i < n; i++, vh++, vt++) */
	/* 	    *vt = TEMP_WATER - BETA_MELT * (*vh) * LEN_SCALING - (*vt); */
	/* } */
	/* phgExportVTK(g, OUTPUT_DIR "temp_init.vtk", ns->T[1], temp_diff, NULL); */
	/* phgDofFree(&temp_diff); */
    }


    //pc = nsPCCreate(ns, ns_params);

    if (ns_params->compute_error_ref) /* Use ref solution as init values */
	init_by_ref_solution(ns);
	
    /* ================================================================================
     *
     * 
     *    Main loop:
     *       1. Steady state:   adaptive refinement.
     *       2. Time dependent: time advance.
     *
     * ================================================================================
     * */
    while (TRUE) {
	static BOOLEAN initialized = FALSE;
	FLOAT time_end = ns_params->time_end;

	elapsed_time(g, FALSE, 0.);	/* reset timer */
	phgGetTime(tt);

	if (Fabs(time_end - Time) < 1e-12) {
	    phgPrintf("\n=======\nTime reach end: %lf, exit.\n", Time);
	    break;
	}

	if (tstep > ns_params->max_tstep) {
	    phgPrintf("\n=======\nTime step reach end: %d, exit.\n", 
		      tstep);
	    break;
	}

#if 0
	/* use time t^{n+1} */
	dt[-1] = dt[0];
	if (Time + dt[0] > time_end)
	    dt[0] = time_end - Time;

	Time += dt[0];
	setFuncTime(Time);
#endif

	phgPrintf("\n==========\ntime: %lf, step:%d\n", (double)Time, tstep);
	phgPrintf("    %d DOF (u:%d, p:%d), %d elements, %d submesh%s, load imbalance: %lg\n",
		  DofGetDataCountGlobal(u[1]) + DofGetDataCountGlobal(p[1]), 
		  DofGetDataCountGlobal(u[1]), DofGetDataCountGlobal(p[1]), 
		  g->nleaf_global, g->nprocs,
		  g->nprocs > 1 ? "es" : "", (double)g->lif);

	/* save mesh */
	if (ns_params->record
	    && tstep % ns_params->step_span == 0) {			
	    phgResumeLogUpdate(g, &Time, &tstep, ns_params->fn, NULL);
	}

	if (!initialized) {
	    /* reset mem_peak */
	    phgMemoryUsageReset();
	    initialized = TRUE;
	}




	/* ------------------------------------------------------------
	 * 
	 *   Time marching 
	 * 
	 * ------------------------------------------------------------ */
	/* update variale,  
	 *   call this routine after time update and/or grid change.*/
	phgNSTimeAdvance(ns, Time, tstep);
	phgPrintf("    update solution");
	elapsed_time(g, TRUE, 0.);
	//load_ref_solution();	/* load ref for interp nu*/

	
	/* --------------------------------------------------------------------------------
	 * 
	 *  Step 3.
	 *
	 *  Momentum Equations.
	 *
	 * -------------------------------------------------------------------------------- */
	elapsed_time(g, FALSE, 0.);	/* reset timer */


	if (ns_params->core_type == STOKES) {
	    phgPrintf("   * Using Full Stokes core\n");

	    /* -------------------------------------------------
	     *
	     * 
	     * Full Stokes core
	     *     non-linear iteration.
	     *
	     * -------------------------------------------------*/

	    int max_nonstep = 0, newton_start = 0;
	    //assert(ns_params->utype == DOF_P2);

	    /* For nonlinear iter */
	    int nonstep = 0; 
	    non_du = non_dp = non_dT = 1e+10;
	    DOF *u_last = phgDofCopy(u[1], NULL, NULL, "u_last");
	    DOF *p_last = phgDofCopy(p[1], NULL, NULL, "p_last");
	    FLOAT non_res_last = 1.;
	    LTYPE ltype_last = PICARD;

	    /* First step, change max non step.  */
	    if (tstep == 1) {	
		if (ns_params->max_nonstep0 > 0)
		    max_nonstep = ns_params->max_nonstep0;
		else
		    max_nonstep = ns_params->max_nonstep;
		if (ns_params->newton_start0 > 0)
		    newton_start = ns_params->newton_start0;
		else
		    newton_start = ns_params->newton_start;
		phgPrintf("   * Set max nonstep to %d for first step.\n", max_nonstep);
		phgPrintf("   * Set Newton start to %d for first step.\n", newton_start);
	    } else {
		max_nonstep = ns_params->max_nonstep;
		newton_start = ns_params->newton_start;
	    }

	    static DOF *gradu_nonstep[100];
	    Bzero(gradu_nonstep);
	    
	    while (TRUE) {
		//break;		/* debug other solver */
		phgPrintf("\n   ==================\n");
		phgPrintf("   Non-linear interation step: %3d (%3d)\n", nonstep, max_nonstep);
	    
		/* Init const viscosity */
		if (ns_params->start_const_vis &&
		    tstep == 1 && nonstep == 0) {
		    phgPrintf("* vis: const\n");
		    ns->viscosity_type = VIS_CONST;
		} else {		
		    phgPrintf("* vis: strain\n");
		    ns->viscosity_type = VIS_STRAIN;		
		}

		if (ns_params->reduce_mesh && gL != NULL)
		    mark_inactive(ns);

		sayHello("Non linear solve begin");
		phgNSInitSolverU(ns);

		if (nonstep < newton_start)
		    ns->ltype = PICARD;
		else 
		    ns->ltype = NEWTON;

		phgPrintf("   Build RHS: ");
		phgUseIsop = USE_ISOP;		/* Note: Isop !!! */
		phgNSBuildSolverURHS(ns);
		phgUseIsop = FALSE;
		elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));

		ns->non_res = 
		    res = phgVecNorm2(ns->solver_u->rhs, 0, NULL);
		if (0)
		    plot_residual(ns, ns->solver_u->rhs, nonstep);
		phgPrintf("   nonlinear residual: %24.12E\n", res);

		/* Restore Picard if no improvement */
		if (ltype_last == NEWTON &&
		    res > non_res_last * .75) {
		    phgPrintf("   !!! Newton step failed, use Picard to run again\n");
		    ns->ltype = PICARD;
		    max_nonstep += 5; /* Add more Picard steps  */

		    /* resotre dofs:
		     * Fix me: temprature */
		    phgDofCopy(u_last, &u[1], NULL, "u_{n+1}");
		    phgDofCopy(p_last, &p[1], NULL, "p_{n+1}");
		    phgDofGradient(u[1], &gradu[1], NULL, "gradu_{n+1}");
		    phgUseIsop = USE_ISOP;		/* Note: Isop !!! */
		    phgNSBuildSolverURHS(ns);
		    phgUseIsop = FALSE;
		    ns->non_res = 
			res = phgVecNorm2(ns->solver_u->rhs, 0, NULL);
		    phgPrintf("   nonlinear residual: %24.12E\n", res);
		} 

		/* save non res */
		non_res_last = res;
		ltype_last = ns->ltype;
	    
		/* build matrices */
		/* if (ns_params->use_PCD) */
		/*     phgNSInitPc(ns); */

		phgPrintf("   Build matrices:\n");
		phgUseIsop = USE_ISOP; /* Note: Isop !!! */
		phgNSBuildSolverUMat(ns);
		phgUseIsop = FALSE;
		phgPrintf("      done ");
		elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));

		/* /\* PC Mat *\/ */
		/* nsPCBuildMat(pc, nonstep); */
		
		/* if (ns_params->use_PCD) { */
		/*     phgPrintf("   Build Pc: \n"); */
		/*     phgNSBuildPc(ns); */
		/*     phgPrintf("      done "); */
		/*     elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL)); */
		/* } */


		/*
		 * solve equation and update (u, p)
		 * */
		phgPrintf("solver tol: %E\n", ns->solver_u->rtol);
		phgSolverSolve(ns->solver_u, FALSE,
			       ns->du, ns->dp, 
#if USE_SLIDING_BC
			       ns->dUn, /* is NULL if sliding_bdry_scheme > 0 */
#endif
			       NULL);
		


#if USE_SLIDING_BC
		if (ns_params->sliding_bdry_scheme == 0) /* Direct */
		    rotate_dof_bases(ns->du, surf_bas, FALSE);
#endif
		phgPrintf("      solver_u: nits = %d, resid = %0.4lg ",
			  ns->solver_u->nits, ns->solver_u->residual);
		elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));

	
		/* save dofs */
		phgDofCopy(u[1], &u_last, NULL, "u_last");
		phgDofCopy(p[1], &p_last, NULL, "p_last");

		/* nonlinear correction */
		phgDofAXPY(1.0, ns->du, &u[1]);
		phgDofAXPY(1.0, ns->dp, &p[1]);
		
		assert(u[1]->type == ns_params->utype);
		assert(p[1]->type == ns_params->ptype);
#if USE_SLIDING_BC
		if (ns_params->sliding_bdry_scheme == 0) /* Direct */
		    dof_set_normal_data(u[1], surf_bas);
		else if (ns_params->sliding_bdry_scheme == 1) /* Surf bdry Proj */
		    phgDofAXPY(1.0, ns->dUn, &ns->Un);
#endif
		PERIODIC_SYNC(u[1]);
		PERIODIC_SYNC(p[1]);

		/* non_du = phgDofNormL2(ns->du); */
		/* non_dp = phgDofNormL2(ns->dp); */
		non_du = phgDofNormInftyVec(ns->du);
		if (ns_params->ptype == DOF_G1) {
		    proj_G1(ns, ns->dp, &dp1);
		    non_dp = phgDofNormInftyVec(dp1);
		} else {
		    non_dp = phgDofNormInftyVec(ns->dp);
		}
#if USE_SLIDING_BC
		if (ns_params->sliding_bdry_scheme == 1) {
		    FLOAT non_dUn = phgDofNormInftyVec(ns->dUn);
		    phgPrintf("   du: %24.12E dp: %24.12E dUn: %24.12E\n",
			      non_du, non_dp, non_dUn);
		} else {
		    phgPrintf("   du: %24.12E dp: %24.12E\n", non_du, non_dp);
		}
#else		
		phgPrintf("   du: %24.12E dp: %24.12E\n", non_du, non_dp);
#endif		
		phgPrintf("   u: [%24.12E, %24.12E]\n", 
			  phgDofMinValVec(u[1]), 
			  phgDofMaxValVec(u[1]));
		phgPrintf("   p: [%24.12E, %24.12E]\n", 
			  phgDofMinValVec(p[1]), 
			  phgDofMaxValVec(p[1]));
		phgDofGradient(u[1], &gradu[1], NULL, "gradu_{n+1}");

		dof_range(ns->u[1]);
		dof_range(ns->p[1]);
		
		/* if (ns_params->use_PCD) */
		/*     phgNSDestroyPc(ns); */
		phgNSDestroySolverU(ns);

		/* evolution of u */
		//DOF_SCALE(u[1], "after solve");
		//DOF_SCALE(p[1], "after solve");
#if 0
#  warning check solution U,p
		sprintf(vtk_file, OUTPUT_DIR "non_%02d_u.vtk", nonstep);
		phgExportVTK(g, vtk_file, u[1], p[1], NULL);
		phgExportEnsightT(g, OUTPUT_DIR "ins_" NS_PROBLEM ,
				  nonstep, nonstep, u[1], p[1], T[1],
				  ns->du, ns->dp, ns->dT, NULL);
#endif



		if (ns_params->output_non_iter
		    && nonstep % ns_params->step_span == 0) {
		    phgPrintf("   Output solution to ensight ");
		    /* phgExportEnsightT(g, OUTPUT_DIR "/ins_" NS_PROBLEM , nonstep, nonstep, */
		    /* 		      u[1], p[1], T[1], NULL);  /\* ensight *\/ */
		    sprintf(vtk_file, OUTPUT_DIR "non_%02d_T.vtk", nonstep);
		    phgExportVTK(g, vtk_file , 
				 u[1], p[1], T[1], ns->du, NULL);
		    elapsed_time(g, TRUE, 0.);
		    //ice_monitor(ns, nonstep);
		}


		/* Linearized */
		if (!ns_params->non_linear
		    && nonstep >= 0) {
		    phgPrintf("   Linearized iteration converges.\n");
		    break;
		}

		phgGetTime(tt1);
		phgPrintf("    time usage of current non step: %lfs\n",
			  (double)(tt1[2] - tt[2]));

		nonstep++;

		/*
		 * Nonliner iteration break, 
		 *   converge for characteristic value.
		 * Velocity: 100 m/a
		 * Pressure: 1e8 Pa
		 *
		 *  */
		const FLOAT U0 = 100;
		const FLOAT P0 = 1e8;

		if (//(res < ns_params->non_tol) 
		    nonstep >= ns_params->min_nonstep 
		    && ((ns->viscosity_type != VIS_CONST && 
			 non_du < ns_params->non_tol * U0
			 && non_dp * PRES_SCALING < ns_params->non_tol * P0
			 )
			|| nonstep > max_nonstep)
		    ) {

		    if (nonstep > max_nonstep) 
			phgWarning("   Non-linear iteration reach max step,"
				   " results may be inaccrate!\n");
		    else
			phgPrintf("   Non-linear iteration converges.\n");
		    break;
		}
	    } /* solve */


#if 0
#  warning  -------- Full Stokesreconstruct velocity w !!! -----------
	    if (ns_params->utype == DOF_P1) {
		//BOOLEAN use_prism_save = ns_params->use_prism_elem;
		//ns_params->use_prism_elem = FALSE;

		/* reconstruct using linear isop */
		reconstruct_velocity_w(ns);
		phgDofGradient(u[1], &gradu[1], NULL, "gradu_{n+1}"); /* update gradu */

		// ns_params->use_prism_elem = use_prism_save;
	    }
#endif
		
	    
	    phgDofFree(&u_last);
	    phgDofFree(&p_last);
	    if (ns_params->stab_scheme >= 1) {
		int istep;
		for (istep = 0; istep < nonstep; istep++)
		    phgDofFree(&gradu_nonstep[istep]);
	    }

	    phgPrintf("Save Dofs\n");
	    save_dof_data3(g, u[1], OUTPUT_DIR"u.dat");
	    save_dof_data3(g, p[1], OUTPUT_DIR"p.dat");
	    DOF_SCALE(u[1], "save");
	    DOF_SCALE(p[1], "save");
	}
	else if (ns_params->core_type == SIA) {

	    /* -------------------------------------------------
	     * 
	     * SIA core.
	     * 
	     * -------------------------------------------------*/
	    phgPrintf("   * Using SIA core\n");

	    ns->viscosity_type = VIS_STRAIN;		
	    core_SIA(ns);
	    phgDofGradient(u[1], &gradu[1], NULL, "gradu_{n+1}");
	}
	else if (ns_params->core_type == FIRST_ORDER) {

	    /* -------------------------------------------------
	     * 
	     * FO core.
	     * 
	     * -------------------------------------------------*/
	    phgPrintf("   * Using FO core\n");

	    core_FO(ns, tstep);
	    phgDofGradient(u[1], &gradu[1], NULL, "gradu_{n+1}");
	}
	else if (ns_params->core_type == DEBUG_CORE1) {

	    /* -------------------------------------------------
	     * 
	     * Debug core. set velocity by given function
	     * 
	     * -------------------------------------------------*/

	    ns->viscosity_type = VIS_STRAIN;
	    phgDofSetDataByFunction(u[1], func_u0);
	    phgDofSetDataByFunction(p[1], func_p0);
	    phgDofGradient(u[1], &gradu[1], NULL, "gradu_{n+1}");
	}
	else if (ns_params->core_type == DEBUG_CORE2) {

	    /* -------------------------------------------------
	     * 
	     * Debug core. Load previous velocity
	     * 
	     * -------------------------------------------------*/

	    phgPrintf("Load Dofs\n");

	    load_dof_data3(g, u[1], OUTPUT_DIR"u.dat", NULL);
	    load_dof_data3(g, p[1], OUTPUT_DIR"p.dat", NULL);
	    DOF_SCALE(u[1], "load");
	    DOF_SCALE(p[1], "load");
	    
	    phgDofGradient(u[1], &gradu[1], NULL, "gradu_{n+1}");
	}


	dof_range(ns->u[1]);
	dof_range(ns->p[1]);

	
	/* Project to continuous gradu */
	if (ns_params->solve_temp) {
	    phgPrintf("\n* Project gradu.\n");
	    proj_gradu(ns, ns->gradu[1], NULL, 0);
	    DOF_SCALE(gradu[1], "grad u");
	    DOF_SCALE(ns->Gradu, "Grad u");
	}

	/* Export gradu */
	if (0) {
	    int i, k;
	    DOF *Gu[DDim], *gu[DDim], *guDG0;
	    guDG0 = phgDofCopy(ns->gradu[1], NULL, DOF_P0, NULL);

	    for (k = 0; k < DDim; k++) {
		FLOAT *vGu;
		INT n;
		char name[1000];

		sprintf(name, "Gu%d", k);
		Gu[k] = phgDofNew(g, DOF_P1, 1, name, DofNoAction);
		vGu = ns->Gradu->data; /* DOF_P1 */
		n = DofGetDataCount(Gu[k]);
		for (i = 0; i < n; i++)
		    Gu[k]->data[i] = vGu[i * DDim + k];

		
		sprintf(name, "gu%d", k);
		gu[k] = phgDofNew(g, DOF_P0, 1, name, DofNoAction);
		vGu = guDG0->data;   /* DOF_P0 */
		n = DofGetDataCount(gu[k]);
		for (i = 0; i < n; i++)
		    gu[k]->data[i] = vGu[i * DDim + k];
	    }

	    
	    phgExportVTK(g, "Gu.vtk",
			 Gu[0], Gu[1], Gu[2],
			 Gu[3], Gu[4], Gu[5],
			 Gu[6], Gu[7], Gu[8],
			 gu[0], gu[1], gu[2],
			 gu[3], gu[4], gu[5],
			 gu[6], gu[7], gu[8],
			 NULL);
	    phgFinalize();
	}





	/* --------------------------------------------------------------------------------
	 *
	 *  Step 4.
	 *
	 *   Solve temperature.
	 *
	 * -------------------------------------------------------------------------------- */


	if (ns_params->solve_temp) {
	    phgPrintf("\n   ==================\n");
	    phgPrintf("   Temperature solve \n");
	    phgPrintf("   ==================\n\n");
	    phgPrintf("   T type: %s\n", T[1]->type->name);

	    elapsed_time(g, FALSE, 0.);	/* reset timer */

	    phgNSInitSolverT(ns);

	    phgPrintf("   Build Mat: ");
	    phgNSBuildSolverTMat(ns, FALSE);
	    elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));
	    phgPrintf("   Build RHS: ");
	    phgNSBuildSolverTRHS(ns, FALSE);
	    elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));

	    phgNSSolverTBuildConstrain(ns);
	    phgDofCopy(ns->T[1], &ns->dT, NULL, "dT");

	    phgNSSolverTSolve(ns, FALSE);
	    elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));

	    phgNSDestroySolverT(ns);

	    find_melt_region(ns);
	    DOF_SCALE(ns->T[1], "after solve");
	    phgDofAXPY(-1.0, ns->T[1], &ns->dT);
	    non_dT = phgDofNormInftyVec(ns->dT);
	    phgPrintf("   dT: %24.12E\n", non_dT);


	    DOF *temp_diff = phgDofCopy(ns->T[1], NULL, NULL, "Td");
	    {
		FLOAT *vt = temp_diff->data;
		const FLOAT *vh = ns->depth_T->data;

		INT i, n = DofGetDataCount(temp_diff);
		for (i = 0; i < n; i++, vh++, vt++)
		    *vt = TEMP_WATER - BETA_MELT * (*vh) *LEN_SCALING  - (*vt);
	    }
	    phgDofFree(&temp_diff);
	} else {
	    phgPrintf("Temp not updated.\n");
	    non_dT = 0.;
	}












	


	/* ------------------------------------------------------------
	 * 
	 *   Error check
	 * 
	 * ------------------------------------------------------------ */
	if (!ns_params->compute_error) {
	    eu = ep = egradu = eT = NULL;
	    ediv = phgDofDivergence(u[1], NULL, NULL, "err div u");
	    phgPrintf(            "            normL2(u, p) = (%20.12E, %20.12E)\n"
				  "            normH1(u)    = (%20.12E)\n"
				  "            normDiv(u)   = (%20.12E)\n",
				  dofNormL2(u[1]), dofNormL2(p[1]),
				  dofNormL2(gradu[1]), dofNormL2(ediv));
	    elapsed_time(g, TRUE, 0.);
	} else {
	    /* ----Error check------- */
	    phgPrintf("    Errors: \n");
	    eu = phgDofCopy(u_exact, NULL, ns_params->utype, "erru");
	    ep = phgDofCopy(p_exact, NULL, DOF_P1, "errp");
	    eT = phgDofCopy(T_exact, NULL, ns_params->T_type, "errT");
	    egradu = phgDofCopy(gradu_exact, NULL, gradu[1]->type, "err grad u");
	    ediv = phgDofDivergence(u[1], NULL, NULL, "err div u");

	    phgExportVTK(g, "diff.vtk",
			 u[1], p[1], T[1], 
			 eu, ep, eT,
			 NULL);
	    
	    phgDofAXPY(-1.0, u[1], &eu);
	    phgDofAXPY(-1.0, gradu[1], &egradu);
	    phgDofAXPY(-1.0, T[1], &eT);
	    if (ns_params->ptype == DOF_G1) {
		proj_G1(ns, p[1], &p1);
		phgDofAXPY(-1.0, p1, &ep);
	    } else {
		phgDofAXPY(-1.0, p[1], &ep);
	    }

	    phgPrintf("            errL2(u, p) = (%20.12E, %20.12E)\n"
		      "            errH1(u)    = (%20.12E)\n"
		      "            errDiv(u)   = (%20.12E)\n",
		      dofNormL2(eu), dofNormL2(ep),
		      //dofDiffNormL2(p[1], p_exact),
		      dofNormL2(egradu), dofNormL2(ediv));
	    elapsed_time(g, TRUE, 0.);

	    phgExportVTK(g, "error.vtk",
			 u[1], p[1], T[1], 
			 eu, ep, eT,
			 NULL);
	    
	    phgDofFree(&egradu);
	}


	//DOF_SCALE(ediv, "check div");
	//check_div(gradu[1], &ns->div0);
	

	getPecletNum(g, u[1], ns_params->nu, 6);	
	mem = phgMemoryUsage(g, &mem_peak);
	phgPrintf("    Memory usage: current %0.4lgMB, peak %0.4lgMB\n",
		  (double)mem / (1024.0 * 1024.0),
		  (double)mem_peak / (1024.0 * 1024.0));



	/* ------------------------------------------------------------
	 * 
	 *   Move mesh
	 * 
	 * ------------------------------------------------------------ */
	if (ns_params->solve_height) {
	    /* Unstructed layered mesh */

	    assert (gL != NULL) ;
	    phgPrintf("Move mesh on UNSTRUCT mesh.\n");
	    get_moved_coord(ns, tstep);
	    move_mesh(ns);
	    check_height(g, gL);

	    phgDofGradient(u[1], &gradu[1], NULL, "gradu_{n+1}");
	    phgDofGradient(u[0], &gradu[0], NULL, "gradu_{n}");
	} else {
	    phgPrintf("Mesh not moved.\n");
	}


	
	if (0)
	    plot_surf_force(ns);



	/* ------------------------------------------------------------
	 * 
	 *   Output solution
	 * 
	 * ------------------------------------------------------------ */
	if (tstep % ns_params->step_span == 0) { 
	    //ice_monitor(ns, tstep);

#if 1
	    /* DOF *p1 = NULL, *p0 = NULL; */
	    if (ns_params->ptype == DOF_G1) 
		proj_G1(ns, p[1], &p1);

	    if (ns_params->compute_error) {
		phgPrintf("    Output solution and error to ensight ");
		if (ns_params->ptype == DOF_G1) 
		    proj_G1(ns, ep, &dp1);
		/* phgExportEnsightT(g, OUTPUT_DIR "/ins_" NS_PROBLEM , Time, tstep, */
		/* 		  u[1], p[1], T[1], eu, ep, eT, ns->div0, p1, dp1, NULL); */
		sprintf(vtk_file, OUTPUT_DIR "ice_%05d.plt", tstep);
		phgExportTecplot(g, vtk_file ,
				 u[1], p[1], T[1], 				 
				 ns->dHt, 
				 eu, ep, eT, 
				 /* ns->div0,  */
				 NULL);
	    } else {
		/* div0 proj */
		/* DOF *div1 = phgDofNew(g, DOF_P1, 1, "div1", DofNoAction); */
		/* proj_G1(ns, ns->div0, &div1); */

		phgPrintf("    Output solution to ensight ");
		/* phgExportEnsightT(g, OUTPUT_DIR "/ins_" NS_PROBLEM , Time, tstep, */
		/* 		  u[1], p[1], T[1], ns->div0, /\* div1， *\/ ns->beta, p1, NULL); */
		/* sprintf(vtk_file, OUTPUT_DIR "ice_%05d.plt", tstep); */
		/* phgExportTecplot(g, vtk_file , */
		/* 		 u[1], p[1], T[1], ns->dHt, ns->div0, ns->beta, p1, NULL); */
		sprintf(vtk_file, OUTPUT_DIR "ice_%05d.vtk", tstep);
		phgExportVTK(g, vtk_file ,
			     u[1], p[1], T[1], 
			     ns->dHt, ns->beta, 
			     ns->sigma_z,
			     //ns->depth_P1,
			     ns->height,
			     /* ns->div0,  */
			     NULL);
	    }
	    
	    elapsed_time(g, TRUE, 0.);
#endif

	    /* Save coord data */
	    {
		sprintf(data_Crd, OUTPUT_DIR "/ins_" NS_PROBLEM "_%05d.dat.Crd", tstep);
		assert(ns->coord->type == DOF_P1);
		save_dof_data3(g, ns->coord, data_Crd);
	    }

	    if (ns_params->record) {			
		/* save dof data for time step {n}, {n+1} */
		sprintf(data_file, OUTPUT_DIR "/ins_" NS_PROBLEM "_%05d.dat", tstep - 1);
		DATA_FILE_SURFIX;
		save_dof_data3(g, u[0], data_u);
		save_dof_data3(g, p[0], data_p);
		save_dof_data3(g, T[0], data_T);
		phgPrintf("   Save u_ {%5d}[%8d]:%24.12E p_ {%5d}[%8d]:%24.12E\n", 
			  tstep - 1, DofGetDataCountGlobal(u[0]), phgDofNormL2(u[0]), 
			  tstep - 1, DofGetDataCountGlobal(p[0]), phgDofNormL2(p[0]));
		phgPrintf("   Save T_{%5d}[%8d]:%24.12E\n", 
			  tstep - 1, DofGetDataCountGlobal(T[0]), phgDofNormL2(T[0])); 
		DOF_SCALE(u[0], "save");
		DOF_SCALE(p[0], "save");
		DOF_SCALE(T[0], "save");
		DOF_SCALE(gradu[0], "save");
		

		sprintf(data_file, OUTPUT_DIR "/ins_" NS_PROBLEM "_%05d.dat", tstep);
		DATA_FILE_SURFIX;
		save_dof_data3(g, u[1], data_u);
		save_dof_data3(g, p[1], data_p);
		save_dof_data3(g, T[1], data_T);
		phgPrintf("   Save u_ {%5d}[%8d]:%24.12E p_ {%5d}[%8d]:%24.12E\n", 
			  tstep, DofGetDataCountGlobal(u[1]), phgDofNormL2(u[1]), 
			  tstep, DofGetDataCountGlobal(p[1]), phgDofNormL2(p[1]));
		phgPrintf("   Save T_{%5d}[%8d]:%24.12E\n", 
			  tstep, DofGetDataCountGlobal(T[1]), phgDofNormL2(T[1]));
		DOF_SCALE(u[1], "save");
		DOF_SCALE(p[1], "save");
		DOF_SCALE(T[1], "save");
		DOF_SCALE(gradu[1], "save");
		phgResumeLogUpdate(g, NULL, NULL, NULL, data_file);
	    } /* end of record */
	    if (gL != NULL) 
		check_height(g, gL);
	    sayHello("After record final solution data");

	}
	if (0)
	    plot_surf_force(ns);


	/* Test on new save load */
	{
	    DOF *dofs[] = {u[1], p[1], T[1], NULL};
	    char *names[] = {"velocity", "pressure", "temperature", NULL};
	    save_state(tstep, Time, 3, dofs, names);

	    phgDofSetDataByValue(u[1], 0.);
	    phgDofSetDataByValue(p[1], 0.);
	    phgDofSetDataByValue(T[1], 0.);

	    load_state(tstep, &Time, 3, dofs, names);
	    phgFinalize();
	    exit(0);
	}

	

	if (0) //(ns_params->save_solution_dat > 0)
	    save_solution_dat(ns);
	if (ns_params->compute_error_ref)
	    compute_error_ref(ns);


	/* clean up */
	phgDofFree(&ediv);
	phgDofFree(&eu);
	phgDofFree(&ep);
	phgDofFree(&eT);



	/* ----------------------------------------------------------------------
	 * 
	 * Compute drag force FD
	 *
	 * ---------------------------------------------------------------------- 
	 * */






        phgGetTime(tt1);
        phgPrintf("    total time usage of current time step: %lfs\n",
		  (double)(tt1[2] - tt[2]));

	if (mem_peak > 1024 * (size_t)ns_params->mem_max * 1024) {
	    phgPrintf("\n=======\nMem usage reach max, exit.\n");
	    break;
	}
	tstep++;


#if 1
	/* use time t^{n} */
	dt[-1] = dt[0];
	if (Time + dt[0] > time_end)
	    dt[0] = time_end - Time;

	Time += dt[0];
	setFuncTime(Time);
#endif
    }				/* end of time advaning */



    /* destroy line block */
    //destroy_line_block(&ns->bk);


    if (ns_params->ptype == DOF_G1) {
	phgDofFree(&p1);
	phgDofFree(&dp1);
    }

    /* destroy reused solver */
    if (ns->solver_u != NULL) {
	/* if (ns_params->use_PCD) */
	/*     phgNSDestroyPc(ns); */
	phgNSDestroySolverU(ns);
    }


    phgNSFinalize(&ns);
    phgDofFree(&u_exact);
    phgDofFree(&p_exact);
    phgDofFree(&gradu_exact);
    phgDofFree(&T_exact);
	    
    phgFreeGrid(&g);
    phgFinalize();
    phgFree(ns_params);

    return 0;
}
