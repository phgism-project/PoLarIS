/* ================================================================================
 *
 *  Nonlinear solver for Temp.
 *  Equations:
 *
 *  rho c dT/dt + rho c u \nabla T - K \Delta T = f    ... (1)
 *    T = T_0        when z = s                        ... (2)
 *    K dT/dn = G    when z = b                        ... (3)
 *    T <= T'                                          ... (*)
 *  
 *  We use a nonlinear iterative solver to solve this problem
 *  For time step k,
 *    1) Build system with out constrain (*),
 *       A * x = b
 *    2) Get marker
 *    
 *    3) Change system using marker,
 *       A ==> A'
 *       b ==> b'
 *    4) solve A' x' = b',
 *    5) if x' does not changes much, end;
 *       else
 *    6) compute residual r = b - A*x'
 *    7) Update marker:
 *         if T_i > T'_i, mark_i = 1 (constrained)
 *         if r_i < 0,    mark_i = 0 (free)
 *    8) Goto 3
 * 
 *  In our implementation, both constrain and marker is recorded in
 *  vector, since it will be more convenient to change matrix and vector.
 *
 *  Routines:
 *  
 *  1. InitSolver
 *  2. Build constrain
 *  3. Build mat
 *  4. Build rhs
 *  5. non solve
 *  6. Destroy solver
 *
 *  
 *
 *
 * ================================================================================
 */

#include "ins.h"
#include "mat_op3.h"
#if USE_PETSC
#  include <petscsys.h>
#  include <petscviewer.h>
#endif
#define _nsp (ns->ns_params)
#define _pcd (ns->pcd)
#define DDim (Dim*Dim)

enum { TEMP_FREE = 0,
       CONSTRAINED = 1};


void phgNSInitSolverT(NSSolver *ns)
{
    MAP *T_map;

    /* dof copy */
    ns->T_shape = phgDofCopy(ns->T[1], NULL, NULL, "T shape");

    /* dof map */
    T_map = ns->T_map = phgMapCreate(ns->T_shape, NULL);

    /* matrices */
    ns->matT = phgMapCreateMat(T_map, T_map);
    ns->matT->handle_bdry_eqns = FALSE; //MAT_HANDLE_BDRY_EQNS;

    /* rhs */
    ns->rhsT = phgMapCreateVec(T_map, 1);

    /* constrain */
    INT nlocal = T_map->nlocal;
    PHG_CALLOC(ns->T_mask, nlocal); /* all free */
    PHG_CALLOC(ns->T_cntr, nlocal);

    int i;
    for (i = 0; i < nlocal; i++)
	ns->T_mask[i] = TEMP_FREE; /* mark free when init */

    return;
}





void phgNSDestroySolverT(NSSolver *ns)
{
    phgInfo(2, "   Destroy solver T\n");

    /* Constrain */
    phgFree(ns->T_mask);
    phgFree(ns->T_cntr);

    /* Mat Rhs */
    phgMatDestroy(&ns->matT);
    phgVecDestroy(&ns->rhsT);


    /* Map, DOF */
    phgMapDestroy(&ns->T_map);
    phgDofFree(&ns->T_shape);
}



void 
phgNSBuildSolverTMat(NSSolver *ns, BOOLEAN init_T)
{
    GRID *g = ns->g;
    SIMPLEX *e;
    DOF **T = ns->T;
    int i, j, q;
    FLOAT Theta = _nsp->Theta, nu = _nsp->nu, Thet1;
    FLOAT *dt = ns->dt;
    int viscosity_type = ns->viscosity_type;

    
    Unused(dt);
    ForAllElements(g, e) {
	int M = T[1]->type->nbas;	/* num of bases of Velocity */
	int order = DofTypeOrder(T[1], e) * 2;
	FLOAT A[M][M], bufT[M], rhsT[M];
	INT I_[M];
	QUAD *quad;
	FLOAT vol, det, diam;
	const FLOAT *w, *p, *gu, *vu, *vT, *vT0, *vTe;

	Bzero(A); Bzero(bufT); Bzero(rhsT); 

	quad = phgQuadGetQuad3D(order);
	if (init_T) {
	    /* Initial temperature:
	     *    Use steady diffusion model.
	     * */
	    diam = phgGeomGetDiameter(g, e);
	    vT = phgQuadGetDofValues(e, ns->T[1], quad); 

	    const FLOAT rho = RHO_ICE;
	    const FLOAT c = HEAT_CAPACITY;
	    const FLOAT K = THEM_CONDUCT;
	    const FLOAT a = SEC_PER_YEAR;

	    p = quad->points;
	    w = quad->weights;
	    for (q = 0; q < quad->npoints; q++) {
		phgGeomGetCurvedJacobianAtLambda(g, e, p, &det);
		vol = fabs(det / 6.);

		for (i = 0; i < M; i++) {
		    const FLOAT *gi = phgQuadGetBasisValues(e, T[1], i, quad) + q;    
		    const FLOAT *ggi = phgQuadGetBasisCurvedGradient(e, T[1], i, quad, q); 
		    for (j = 0; j < M; j++) {
			const FLOAT *gj = phgQuadGetBasisValues(e, T[1], j, quad) + q;       
			const FLOAT *ggj = phgQuadGetBasisCurvedGradient(e, T[1], j, quad, q); 

			FLOAT diff_T = INNER_PRODUCT(ggi, ggj) / LEN_SCALING2;
			A[i][j] += vol*(*w) * K * diff_T * a * EQU_T_SCALING;
		    }
		}

		w++; p += Dim + 1;
	    }
	} else {
	    /* 
	     * Convection diffusion model.
	     * */
	    diam = phgGeomGetDiameter(g, e);
	    vT = phgQuadGetDofValues(e, ns->T[1], quad); 
	    vu = phgQuadGetDofValues(e, ns->u[1], quad); 
	    gu = phgQuadGetDofValues(e, ns->Gradu, quad); 
	    vT0 = phgQuadGetDofValues(e, ns->T[0], quad); 

	    const FLOAT rho = RHO_ICE;
	    const FLOAT c = HEAT_CAPACITY;
	    const FLOAT K = THEM_CONDUCT;
	    const FLOAT a = SEC_PER_YEAR;


	    p = quad->points;
	    w = quad->weights;
	    for (q = 0; q < quad->npoints; q++) {
		phgGeomGetCurvedJacobianAtLambda(g, e, p, &det);
		vol = fabs(det / 6.);

		nu = get_effective_viscosity(gu, *vT0, 0, viscosity_type);

#if USE_TEMP_SDPG
		/* Stream line diffusion PG */
		FLOAT U = sqrt(INNER_PRODUCT(vu, vu));
		FLOAT h = diam;
		FLOAT P_h = U * h / (K / LEN_SCALING) * (rho * c / a);
		FLOAT SD_delta = 0.;
		if (P_h > 1)
		    SD_delta = h / (2*U) * (1. - 1./P_h);
#endif
	
		for (i = 0; i < M; i++) {
		    const FLOAT *gi = phgQuadGetBasisValues(e, T[1], i, quad) + q;    
		    const FLOAT *ggi = phgQuadGetBasisCurvedGradient(e, T[1], i, quad, q); 
		    for (j = 0; j < M; j++) {
			const FLOAT *gj = phgQuadGetBasisValues(e, T[1], j, quad) + q;       
			const FLOAT *ggj = phgQuadGetBasisCurvedGradient(e, T[1], j, quad, q); 

			FLOAT qmass = (*gj) * (*gi);
			//FLOAT conv = INNER_PRODUCT(vu, ggj) * (*gi) / LEN_SCALING;
			FLOAT diff_T = INNER_PRODUCT(ggi, ggj) / LEN_SCALING2;
			FLOAT test_phi_i = *gi + SD_delta * (INNER_PRODUCT(vu, ggi));

			A[i][j] += vol*(*w) * (
#if USE_TEMP_TIME
					       rho * c * (*gj) * (*gi) / dt[0]    
#endif
#if USE_TEMP_CONV
					       + rho * c * (INNER_PRODUCT(vu, ggj) 
							    * test_phi_i ) / LEN_SCALING
#endif
#if USE_TEMP_DIFF
					       + K * diff_T * a        
					       //+ K * ggj[Z_DIR] * ggi[Z_DIR] / LEN_SCALING2 * a
#endif
					       ) * EQU_T_SCALING;
		    }
		}

		vu += Dim;
		gu += DDim;
		vT++;
		vT0++;
		w++; p += Dim + 1;
	    }

	}

	/* Map: Element -> system */
	for (i = 0; i < M; i++)
	    I_[i] = phgMapE2L(ns->matT->cmap, 0, e, i);

	/* Global res */
	for (i = 0; i < M; i++) {
	    if (phgDofDirichletBC_(ns->T[1], e, i, NULL, bufT, &rhsT[i],
				   DOF_PROJ_NONE)) {
		phgMatAddEntries(ns->matT, 1, I_ + i, M, I_, bufT);
	    } else {
		phgMatAddEntries(ns->matT, 1, I_ + i, M, I_, &(A[i][0])); 
	    }
	}
    }				/* end element */
 
    if (DUMP_MAT_VEC) {
	static int cntr = 0;
	static char fname[1000];
	phgPrintf("Dumping MatT\n");
	sprintf(fname, "MT%d_.m", cntr++);
	phgMatDumpMATLAB(ns->matT, "T", fname);
    }
 
    return;
}


static double _t1 = 0;
void 
func_T_wrapper(FLOAT x, FLOAT y, FLOAT z, FLOAT *value)
{
    func_T_t(x, y, z, _t1, value);
}


void 
phgNSBuildSolverTRHS(NSSolver *ns, BOOLEAN init_T)
{
    GRID *g = ns->g;
    SIMPLEX *e;
    DOF **T = ns->T;
    int i, j, q;
    FLOAT *dt = ns->dt;
    FLOAT Theta = _nsp->Theta, nu = _nsp->nu, Thet1;
    MAP *map = ns->T_map;
    VEC *rhsT = ns->rhsT;
    int viscosity_type = ns->viscosity_type;

    _t1 = ns->time[1];

    Unused(dt);
    Unused(j);
    phgVecDisassemble(rhsT);
    ForAllElements(g, e) {
	int s, M = T[1]->type->nbas;	/* num of bases of Velocity */
	int order = DofTypeOrder(T[1], e) * 2;
	FLOAT rhs[M], buf[M], T_value[M];
	INT I_[M];
	QUAD *quad;
	FLOAT vol, det, area, diam;
	const FLOAT *w, *p, *vu, *gu, *vT, *vT0, *normal;

	Bzero(rhs);

	quad = phgQuadGetQuad3D(order);
	phgDofGetElementData(ns->T[1], e, T_value);

	if (init_T) {
	} else {
	    diam = phgGeomGetDiameter(g, e);
	    vT = phgQuadGetDofValues(e, ns->T[1], quad); 
	    vu = phgQuadGetDofValues(e, ns->u[1], quad); 
	    gu = phgQuadGetDofValues(e, ns->Gradu, quad); 
	    vT0 = phgQuadGetDofValues(e, ns->T[0], quad); 

	    const FLOAT rho = RHO_ICE;
	    const FLOAT c = HEAT_CAPACITY;
	    const FLOAT K = THEM_CONDUCT;
	    const FLOAT a = SEC_PER_YEAR;

	    p = quad->points;
	    w = quad->weights;
	    for (q = 0; q < quad->npoints; q++) {
		phgGeomGetCurvedJacobianAtLambda(g, e, p, &det);
		vol = fabs(det / 6.);

		FLOAT gT[Dim] = {0, 0, 0};
		FLOAT eu[Dim*Dim];
		MAT3_SYM(gu, eu);

		for (i = 0; i < M; i++) {
		    const FLOAT *ggi = phgQuadGetBasisGradient(e, T[1], i, quad) + q*Dim;
		    int k;
		    for (k = 0; k < Dim; k++)
			gT[k] += T_value[i] * ggi[k];
		}

		nu = get_effective_viscosity(gu, *vT0, 0, viscosity_type);
#if 1
		/* Stream line PG */
		FLOAT U = sqrt(INNER_PRODUCT(vu, vu));
		FLOAT h = diam;
		FLOAT P_h = U * h / (K / LEN_SCALING) * (rho * c / a);
		FLOAT SD_delta = 0.;
		if (P_h > 1)
		    SD_delta = h / (2*U) * (1. - 1./P_h);
#endif

		for (i = 0; i < M; i++) {
		    const FLOAT *gi = phgQuadGetBasisValues(e, T[1], i, quad) + q;    
		    const FLOAT *ggi = phgQuadGetBasisCurvedGradient(e, T[1], i, quad, q); 
		    FLOAT test_phi_i = *gi + SD_delta * (INNER_PRODUCT(vu, ggi));

		    FLOAT conv = INNER_PRODUCT(vu, gT) * (*gi);
		    FLOAT diff_T = INNER_PRODUCT(gT, ggi);
		    FLOAT diff_u = MAT3_NORM2_2(eu) * (*gi);
		    
		    const FLOAT rho = RHO_ICE;
		    const FLOAT c = HEAT_CAPACITY;
		    const FLOAT K = THEM_CONDUCT;
		    const FLOAT a = SEC_PER_YEAR;
		    //FLOAT q = HEAT_SOURCE;
		    FLOAT strain_heat = 0;

		    if (0) {
			FLOAT ux, uy, uz, vx, vy, vz;
	    
			ux = gu[0] / LEN_SCALING; 
			uy = gu[1] / LEN_SCALING; 
			uz = gu[2] / LEN_SCALING;
			vx = gu[3] / LEN_SCALING; 
			vy = gu[4] / LEN_SCALING; 
			vz = gu[5] / LEN_SCALING;

			phgInfo(0, "   strain: %e %e %e\n", 
				(.5*gu[2]*gu[2] + .5*gu[5]*gu[5]) * (*gi)  / LEN_SCALING2,
				2* (ux*ux + vy*vy + ux*vy
				    + 0.25 * pow(uy + vx, 2)
				    + 0.25 * (uz*uz + vz*vz)) * (*gi),
				MAT3_NORM2_2(eu) * (*gi) / LEN_SCALING2
				);
		    }

		    /* Strain heat */
		    if (ns_params->core_type == SIA) {

#if 0
#   warning SIA Strain heat: P1(Gu)
			/* Strain heat, from P1(Gu)  */
			strain_heat = nu * (.5*gu[2]*gu[2] + .5*gu[5]*gu[5]) * (*gi)  / LEN_SCALING2;
#else
#   warning SIA Strain heat: nabla u		
			/* Strain heat, from nabla u
			 * 
			 * PHI = 2 * A(T) * (rho g (s-z))^{n+1} |\nabla s|^{n+1}
			 *
			 * Rutt et al.: GLIMMER Ice sheet model, P3-4
			 * */

			const FLOAT n = POWER_N;
			const FLOAT a = SEC_PER_YEAR;
			const FLOAT T0 = ARRHENIUS_T;
			const FLOAT Q0 = ARRHENIUS_Q0;
			const FLOAT a0 = ARRHENIUS_A0;
			const FLOAT Q1 = ARRHENIUS_Q1;
			const FLOAT a1 = ARRHENIUS_A1;
			const FLOAT R  = GAS_CONSTANT;

			const FLOAT *lambda = p;
			FLOAT vT, vgS[2], vdepth, norm_gS, A;

			phgDofEval(ns->T[1], e, lambda, &vT);
			phgDofEval(ns->depth_P1, e, lambda, &vdepth);
			phgDofEval(ns->dof_gradS, e, lambda, vgS); 
			norm_gS = sqrt(vgS[0]*vgS[0] + vgS[1]*vgS[1]);

			if (vT < T0)
			    A = a0 * exp(-Q0 / R / vT);
			else
			    A = a1 * exp(-Q1 / R / vT);
			A *= SEC_PER_YEAR;
#       if CONST_A
			A = 1e-16;
#       endif
			strain_heat = 2*A
			    * pow(RHO_ICE * GRAVITY * (vdepth) * LEN_SCALING, n+1)
			    * pow(norm_gS, n+1)
			    * (*gi) / LEN_SCALING2;
#endif

		    } 
		    else if (ns_params->core_type == FIRST_ORDER) {
			FLOAT ux, uy, uz, vx, vy, vz;
	    
			ux = gu[0] / LEN_SCALING; 
			uy = gu[1] / LEN_SCALING; 
			uz = gu[2] / LEN_SCALING;
			vx = gu[3] / LEN_SCALING; 
			vy = gu[4] / LEN_SCALING; 
			vz = gu[5] / LEN_SCALING;

			strain_heat = nu * 2* (ux*ux + vy*vy + ux*vy
					       + 0.25 * pow(uy + vx, 2)
					       + 0.25 * (uz*uz + vz*vz)) * (*gi);
		    } 
		    else if (ns_params->core_type == STOKES) {
			strain_heat = nu * MAT3_NORM2_2(eu) * (*gi) / LEN_SCALING2;
		    } 
		    else if (ns_params->core_type == DEBUG_CORE1) {
			strain_heat = 0.;
		    }
		    else {
			phgError(0, "Unknown core type!!!\n");
		    }

		    rhs[i] += vol*(*w) * (
#if USE_TEMP_TIME
					  rho * c * (*vT0)*(*gi) / dt[0]
#endif
#if USE_TEMP_HEAT
					  //+ nu * MAT3_NORM2_2(gu) * (test_phi_i) / LEN_SCALING2
					  //+ q * a
					  + strain_heat
#endif
					  +0
					  ) * EQU_T_SCALING;
		}

		vu += Dim;
		gu += DDim;
		vT++;
		vT0++;
		w++; p += Dim + 1;
	    }
	}

	/* Geo thermal flux */
	for (s = 0; s < NFace; s++) {
	    if (e->bound_type[s] & BC_BOTTOM) {
		int v0, v1, v2;
		int nbas_face = NbasFace(ns->T[1]);
		SHORT bases[nbas_face];
		FLOAT lambda[Dim + 1];
		order = DofTypeOrder(ns->T[1], e) * 3 - 1;

		phgDofGetBasesOnFace(ns->T[1], e, s, bases);
		v0 = GetFaceVertex(s, 0);
		v1 = GetFaceVertex(s, 1);
		v2 = GetFaceVertex(s, 2);
		lambda[s] = 0.;
		    
		area = phgGeomGetFaceArea(g, e, s);
		normal = phgGeomGetFaceOutNormal(g, e, s);
		quad = phgQuadGetQuad2D(order);
		
		p = quad->points;
		w = quad->weights;
		for (q = 0; q < quad->npoints; q++) {
		    FLOAT vu[Dim];
		    lambda[v0] = *(p++);
		    lambda[v1] = *(p++);
		    lambda[v2] = *(p++);

		    const FLOAT *gi_T = 
			ns->T[1]->type->BasFuncs(ns->T[1], e, 0, -1, lambda);

		    for (i = 0; i < nbas_face; i++) {
			int ii = bases[i];
			const FLOAT a = SEC_PER_YEAR;
#if USE_TEMP_GEOFLUX
			rhs[ii] += area*(*w) * (GEOTHE_FLUX) * (gi_T[ii]) * a
			    / LEN_SCALING * EQU_T_SCALING ;
#endif
		    } /* end of bas_i */

		    w++;
		} /* end of quad point */
	    }	  /* end of face outflow */
	}	  /* end of faces */

	for (i = 0; i < M; i++) 
	    if (phgDofDirichletBC_(ns->T[1], e, i, func_T_wrapper, buf, rhs + i,
				   DOF_PROJ_NONE)) {
		;
	    }



	/* Map: Element -> system */
	for (i = 0; i < M; i++)
	    I_[i] = phgMapE2L(ns->matT->cmap, 0, e, i);

	phgVecAddEntries(rhsT, 0, M, I_, &rhs[0]);
    }				/* end element */
    phgVecAssemble(rhsT);

    if (DUMP_MAT_VEC) {
	static int cntr = 0;
	static char fname[1000];
	phgPrintf("Dumping rhsT\n");
	sprintf(fname, "bT%d_.m", cntr++);
	phgVecDumpMATLAB(ns->solver_T->rhs, "bT", fname);
    }

    return;
}


/*
 * See EISMINT Benchmark.
 * At the ice-bedrock interface, mixed boundary conditions are emplyerd:
 *   if T < T', dT/dz = -G/k;
 *   else,      T = T',
 *   where T' = T_0 - beta H.
 * Here we build this constrain for vector, rather than DOF,
 *   since we may use this in the solvering interation.
 *  
 * */
void 
phgNSSolverTBuildConstrain(NSSolver *ns)
{
    GRID *g = ns->g;
    LAYERED_MESH *gL = ns->gL;
    SIMPLEX *e;
    DOF *T = ns->T[1];
    MAP *T_map = ns->T_map;
    BOOLEAN *T_mask;
    INT nlocal = T_map->nlocal;
    int ii, i, j;
    DOF *dof_depth = NULL;


    /* depth and T are all P2 */
    if (ns->T[1]->type == DOF_P1)
	dof_depth = ns->depth_P1;
    else 
	dof_depth = ns->depth_T;
    INT dof_size = DofGetDataCount(dof_depth);

    FLOAT *vH = dof_depth->data;
    FLOAT *T_cntr = ns->T_cntr;
    ForAllElements(g, e) {
	int N = T->type->nbas;
	INT iD, iL, iV;

	for (i = 0; i < N; i++) {
	    /* Dof idx */
	    iD = phgDofMapE2D(dof_depth, e, i);
	    assert(iD >= 0 && iD < dof_size);
	    /* local idx */
	    iL = phgMapE2L(T_map, 0, e, i); /* index in Dof data */
	    assert(iL >= 0);
	    /* vector idx */
	    iV = phgMapL2V(T_map, iL);      /* index in Vec data */
	    assert(iV < T_map->localsize && iV >= 0);
	    
	    if (iV < nlocal) {
		T_cntr[iV] = TEMP_WATER - BETA_MELT * vH[iD] * LEN_SCALING;
		/* Error: Scaling not right !!! */
	    }
	}
    }

    /* Output check */
    if (0) {
	char name[1000];
	sprintf(name, "cntrT.plt");
	DOF *dof_cntr = phgDofCopy(ns->T[1], NULL, NULL, "T_cntr");
	phgMapLocalDataToDof(T_map, 1, &dof_cntr, T_cntr);
	phgExportTecplot(g, name, dof_cntr, NULL);
	phgDofFree(&dof_cntr);
    }

    return;
}

#define EPS_T 1e-10

void 
phgNSSolverTSolve(NSSolver *ns, BOOLEAN init_T)
{
    GRID *g = ns->g;
    BOOLEAN *T_mask = ns->T_mask;
    FLOAT *T_cntr = ns->T_cntr;
    MAP *T_map = ns->T_map;
    MAT *matT = ns->matT;
    VEC *rhsT = ns->rhsT;
    VEC *vecT= phgMapCreateVec(T_map, 1);
    VEC *Res = NULL;
    INT nlocal = T_map->nlocal;
    FLOAT res0, res;
    BOOLEAN remote_data = FALSE;
    INT niter = 0, ncntr;
    int i, j, k;


    phgPrintf("Solver T: constrained solve\n");

    /* ------------------------------------------------------------
     * Step 1. Set Masker,
     *   since the height is updated, the masker need to be updated too.
     * ------------------------------------------------------------ */
    if (init_T) 
	phgMapDofToLocalData(T_map, 1, &ns->T[1], vecT->data);
    else
	phgMapDofToLocalData(T_map, 1, &ns->T[0], vecT->data);
    for (i = 0; i < nlocal; i++) {
	if (vecT->data[i] >= T_cntr[i] - EPS_T)
	    T_mask[i] = CONSTRAINED;
    }
    if (T_map->nprocs > 1)
	remote_data = TRUE;


    niter = 0;
    while (TRUE) {

	phgPrintf("\n\n\n   --------------------------\n");
	phgPrintf("   * Temp non iter step: %d\n", niter);
	ncntr = 0;
	for (i = 0; i < nlocal; i++)
	    if (T_mask[i] == CONSTRAINED)
		ncntr++;
	INT n = ncntr;
	MPI_Reduce(&n, &ncntr, 1, PHG_MPI_INT, MPI_SUM, 0, g->comm);
	phgPrintf("   # constrained %d\n", ncntr);

	/* Output check */
	if (0) {
	    char name[1000];
	    sprintf(name, "Tmask_%03d.plt", niter);
	    DOF *dof_mask = phgDofCopy(ns->T[1], NULL, NULL, "T_mask");
	    FLOAT *vec_mask = phgCalloc(nlocal, sizeof(FLOAT));
	    for (i = 0; i < nlocal; i++)
		vec_mask[i] = T_mask[i];
	    phgMapLocalDataToDof(T_map, 1, &dof_mask, vec_mask);
	    phgExportTecplot(g, name, dof_mask, NULL);
	    phgDofFree(&dof_mask);
	}



	/* ------------------------------------------------------------
	 * Step 2. Change matrix
	 * ------------------------------------------------------------ */
	MAT *A_fix = NULL;
	VEC *b_fix = NULL;
	VEC *oldT = NULL;
	phgMatCopy(matT, &A_fix);
	phgVecCopy(rhsT, &b_fix);
	phgVecCopy(vecT, &oldT);
	phgMatPack(A_fix);

	{
	    int n, jcol, *pc, *pc0, *pc_offp, nlocal;
	    size_t *ps = A_fix->packed_ind;
	    FLOAT *pd, *pd0, *pd_offp;

	    nlocal = A_fix->rmap->nlocal;
	    pc0 = A_fix->packed_cols;
	    pd0 = A_fix->packed_data;

	    for (i = 0; i < nlocal; i++)
		if (T_mask[i] == CONSTRAINED) { 
		    pc = pc0 + ps[i];
		    pd = pd0 + ps[i];
		    n = ps[i+1] - ps[i];
			
		    /* mat row: local columes */
		    for (j = 0; j < n; j++) {
			jcol = *(pc++);
			if (jcol != i) {
			    *(pd++) = 0.;
			} else {
			    *(pd++) = 1.;
			}
		    }
	    
		    /* remote mat: remote columes  */
		    if (remote_data) {
			pc_offp = pc0 + ps[nlocal + i];
			pd_offp = pd0 + ps[nlocal + i];
			n = ps[nlocal + i + 1] - ps[nlocal + i];

			for (j = 0; j < n; j++) {
			    jcol = *(pc_offp++);
			    *(pd_offp++) = 0.;
			}
		    } /* end remote */

		    /* set vector */
		    b_fix->data[i] = T_cntr[i];
		}	  /* end constrian */
	}
	A_fix->handle_bdry_eqns = FALSE;


	/* ------------------------------------------------------------
	 * Step 3. Linear Solver
	 * ------------------------------------------------------------ */
	/* solver_T */
	elapsed_time(g, FALSE, 0.);	/* reset timer */
	phgPrintf("   Linear Solve: \n");
	phgOptionsPush();
	phgOptionsSetOptions("-solver gmres "
			     "-solver_maxit 1000 "
			     "-solver_rtol 1e-10");
	phgOptionsSetOptions(_nsp->T_opts);
	ns->solver_T = phgMat2Solver(SOLVER_DEFAULT, A_fix);
	phgSolverAssemble(ns->solver_T);
	phgOptionsPop();
	ns->solver_T->rhs_updated = TRUE;

	/* linear solve */
	phgVecCopy(b_fix, &ns->solver_T->rhs);
	phgSolverVecSolve(ns->solver_T, TRUE, vecT);
	phgPrintf("      solver_T: nits = %d, resid = %0.4lg ",
		  ns->solver_T->nits, ns->solver_T->residual);
	phgMapLocalDataToDof(T_map, 1, &ns->T[1], vecT->data);
	elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));

	phgSolverDestroy(&ns->solver_T);
	phgMatDestroy(&A_fix);	
	phgVecDestroy(&b_fix);	

	/* Output check */
	if (0) {
	    char name[1000];
	    sprintf(name, "nonT_%03d.plt", niter);
	    phgMapLocalDataToDof(T_map, 1, &ns->T[1], vecT->data);
	    phgExportTecplot(g, name, ns->T[1], NULL);
	}


	FLOAT dT = 0;
	phgVecAXPBY(-1., vecT, 1., &oldT);
	phgVecNormInfty(oldT, 0, &dT);
	phgVecDestroy(&oldT);
	phgPrintf("   dT %e\n", dT);
	if (dT < 1e-6)
	    break;
	
	
	/* ------------------------------------------------------------
	 * Step 4. Linear Solver
	 * ------------------------------------------------------------ */
	/* Residual */
	phgVecCopy(rhsT, &Res);
	phgMatVec(MAT_OP_N, -1., matT, vecT, 1., &Res);
	
	for (i = 0; i < nlocal; i++) {

	    /*
	     * Note:
	     * 
	     * The linear solver may solve the constrained value with
	     * some small error, then constrained Temp vecT[i] maybe not
	     * equal to T_cntr[i].
	     *
	     * To avoid set marker back and forth,
	     * set a tolerance for reaching constrain.
	     *
	     * */
	    if (T_mask[i] == TEMP_FREE) {
		if (vecT->data[i] >= T_cntr[i] - EPS_T)
		    T_mask[i] = CONSTRAINED;
	    } else {
		//T_mask[i] == CONSTRAINED 
		if (Res->data[i] < 0)
		    T_mask[i] = TEMP_FREE;
	    }
	}
	phgVecDestroy(&Res);



	/* Output T_mask to check */
	if (1) {
	    char name[1000];
	    sprintf(name, "Tmask_%03d.plt", niter);
	    DOF *dof_mask = phgDofCopy(ns->T[1], NULL, NULL, "T_mask");
	    FLOAT *vec_mask = phgCalloc(nlocal, sizeof(FLOAT));
	    for (i = 0; i < nlocal; i++)
		vec_mask[i] = T_mask[i];
	    phgMapLocalDataToDof(T_map, 1, &dof_mask, vec_mask);
	    phgExportTecplot(g, name, dof_mask, NULL);
	    phgDofFree(&dof_mask);
	}

	niter++;
    }


    phgVecDestroy(&vecT);
    return;
}




/*
 *
 * Init temperature by diffusion equation.
 *
 * */
static void
init_temp_diff(NSSolver *ns) 
{
    GRID *g = ns->g;
    DOF **T = ns->T;

    elapsed_time(g, FALSE, 0.);	/* reset timer */
    phgNSInitSolverT(ns);

    phgPrintf("   Build Mat: ");
    phgNSBuildSolverTMat(ns, TRUE);
    elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));
    phgPrintf("   Build RHS: ");
    phgNSBuildSolverTRHS(ns, TRUE);
    elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));

    /* Constrain */
    phgNSSolverTBuildConstrain(ns);

    /* Solve */
    phgNSSolverTSolve(ns, TRUE);

    elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));
    phgNSDestroySolverT(ns);

    //find_melt_region(ns);
    //phgDofCopy(T[1], &T[0], NULL, "T_{n}");
    
    if (0)
		phgExportVTK(g, OUTPUT_DIR "Temp_init.vtk", T[1], NULL);
}


/*
 *
 * Init temperature by interpolation.
 *
 * */
static void
init_temp_interp(NSSolver *ns) 
{
    GRID *g = ns->g;
    DOF **T = ns->T, *T_P1;
    LAYERED_MESH *gL = ns->gL;
    int i, j, ii, iG, nly = gL->max_nlayer;
    BTYPE btype = T[1]->DB_mask;
    FILE *fp;

    phgPrintf("\n   ==================\n");
    phgPrintf("   Temperature init interp, from top\n");
    phgPrintf("   ==================\n\n");
    phgPrintf("   T type: %s\n", T[1]->type->name);

    T_P1 = phgDofNew(g, DOF_P1, 1, "temp P1", func_T);
    
    for (ii = 0; ii < gL->nvert_bot; ii++) {
	i = gL->vert_bot_Lidx[ii];
	assert(gL->vert_local_lists[i] != NULL);

	iG = gL->vert_bot_Gidx[ii];
	assert(iG < gL->nvert);

	int nv = gL->vert_local_lists[i][0];
	int *iL = &gL->vert_local_lists[i][1];
	assert(nv > 0);

	/* interp from top to bottom,
	 *   top: surface temp observed
	 *   bot: melting point
	 * */
	FLOAT Tn = *DofVertexData(T_P1, iL[nv-1]);
	FLOAT T0 = TEMP_WATER;
	FLOAT dT = (Tn - T0) / (nv - 1);

	for (j = 0; j < nv; j++) {
	    *DofVertexData(T_P1, iL[j])
		= T0 + j * dT;
	}
    }

    /* interp */
    phgDofCopy(T_P1, &T[1], NULL, "T_{n+1}");
    T[1]->DB_mask = btype;	/* restore */
    //phgDofGradient(T[1], &ns->gradT[1], NULL, "gradT_{n+1}");
    phgExportTecplot(g, "init_T_interp.plt", T[1], NULL);
    
    phgDofFree(&T_P1);
    DOF_SCALE(ns->T[1], "interp");
}



/*
 *
 * Init temperature by read data.
 *   temp data: temperature.txt
 *   beta data: beta.txt
 *   
 * */
static void
init_temp_data(NSSolver *ns) 
{
    GRID *g = ns->g;
    DOF **T = ns->T, *T_P1;
    LAYERED_MESH *gL = ns->gL;
    int i, j, ii, iG, nly = gL->max_nlayer;
    FILE *fp;

    /* Temp */
    if (0) {
	phgPrintf("\n   ==================\n");
	phgPrintf("   Temperature init interp, from 2D data\n");
	phgPrintf("   ==================\n\n");
	phgPrintf("   T type: %s\n", T[1]->type->name);

	FLOAT *data;
	data = phgCalloc(gL->nvert * (nly+1), sizeof(data));

	T_P1 = phgDofNew(g, DOF_P1, 1, "temp P1", DofNoAction);
    
	fp = fopen("./temperature.txt", "r");

	for (j = 0; j < nly + 1; j++)
	    for (i = 0; i < gL->nvert; i++)
		fscanf(fp, "%lf", data + j*gL->nvert + i);

	for (ii = 0; ii < gL->nvert_bot; ii++) {
	    i = gL->vert_bot_Lidx[ii];
	    assert(gL->vert_local_lists[i] != NULL);

	    iG = gL->vert_bot_Gidx[ii];
	    assert(iG < gL->nvert);

	    int nv = gL->vert_local_lists[i][0];
	    int *iL = &gL->vert_local_lists[i][1];
	    assert(nv > 0);

#if 0
	    /* interp from top to bottom */
	    FLOAT Tn = data[iG];
	    FLOAT T0 = TEMP_WATER;
	    FLOAT dT = (Tn - T0) / nv;

	    for (j = 0; j < nv; j++) {
		*DofVertexData(T_P1, iL[j])
		    = T0 + j * dT;
	    }
#else
	    /* use read-in data */
	    for (j = 0; j < nv; j++) {
		*DofVertexData(T_P1, iL[j])
		    = data[j * gL->nvert + iG];
	    }	    
#endif
	}

	/* interp */
	phgDofCopy(T_P1, &T[1], NULL, "T_{n+1}");
	//phgDofGradient(T[1], &ns->gradT[1], NULL, "gradT_{n+1}");
	phgExportTecplot(g, "init_T_interp.plt", T[1], NULL);
    
	phgDofFree(&T_P1);
	DOF_SCALE(ns->T[1], "interp");
	fclose(fp);
    }




    return;
}



void 
phgNSTempInit(NSSolver *ns)
{
    phgPrintf("\n   ==================\n");
    phgPrintf("   Temperature init \n");
    phgPrintf("   ==================\n\n");

    if (_nsp->init_temp_type == 0) {
	init_temp_diff(ns);
    }
    else if (_nsp->init_temp_type == 1) {
	init_temp_interp(ns);
    }
    else if (_nsp->init_temp_type == 2) {
	init_temp_data(ns);
    } 
    else {
	phgError(0, "init_temp_type %d unknown!\n", _nsp->init_temp_type);
    }

}



/*
 * Find melting region.
 *
 *  */
void
find_melt_region(NSSolver *ns)
{
    GRID *g = ns->g;
    SIMPLEX *e;
    DOF *T = ns->T[1];
    DOF *melt = NULL;
    int i, s, n;
    FLOAT Area = 0.;
    FLOAT *vH = ns->depth_T->data, *vM;

    //ERROR_MSG("Disable melt region.");
    assert(ns->depth_T->type == ns->T[1]->type); 
    melt = phgDofNew(g, ns->T[1]->type, 1, "melt mark", DofNoAction);
    n = DofGetDataCount(melt);
    vM = melt->data;
    for (i = 0; i < n; i++, vH++, vM++) {
		*vM = TEMP_WATER - BETA_MELT * (*vH) * LEN_SCALING;
    }
    //phgDofDump(ns->depth_T); 
    //phgDofDump(melt); 
    //phgDofDump(T); 
    phgDofAXPBY(-1., T, 1, &melt);
    if (ns_params->output_melt)
		phgExportVTK(g, "output/melt.vtk", melt, NULL);

    ForAllElements(g, e) {
		int v0, v1, v2;
		FLOAT area, vmark, lambda[Dim+1];

		e->mark = 0;
		for (s = 0; s < NFace; s++) {
			if (e->bound_type[s] & BC_BOTTOM) {
				v0 = GetFaceVertex(s, 0);
				v1 = GetFaceVertex(s, 1);
				v2 = GetFaceVertex(s, 2);
				lambda[v0] = 1./3.;
				lambda[v1] = 1./3.;
				lambda[v2] = 1./3.;
				lambda[s] = 0.;
		    
				phgDofEval(melt, e, lambda, &vmark);
				phgInfo(0, "vmark %e\n", vmark);
				if (vmark < 1e-3) { /* fixed */
					area = phgGeomGetFaceArea(g, e, s);
					Area += area;
					e->mark = 1;
				}
			}
		}
    }
    phgInfo(0, "   Local melting regain: %e\n", Area);

    FLOAT tmp = Area;
    MPI_Allreduce(&tmp, &Area, 1, MPI_DOUBLE, MPI_SUM, g->comm);
    phgPrintf("   Melting regain: %e\n", Area);

    phgDofFree(&melt);
    return;
}


















/* -------------------------------------------------------------------- 
 *
 *
 *  L2 projection of gradient 
 * 
 *                                                                      
 * -------------------------------------------------------------------- */

    

void
proj_gradu(NSSolver *ns, DOF *gradu,
	   const BYTE *components, int type)
{
    GRID *g = ns->g;
    SIMPLEX *e;
    int i, j, k, q;
    DOF *Gradu, *Gradu0;
    SOLVER *solver_Gu;
    const BYTE all_components[DDim] = {1, 1, 1,
				       1, 1, 1,
				       1, 1, 1};
    if (components == NULL)
	components = all_components;


    if (type == 0) {
	/* ------------------------------------------------------------
	 *
	 *
	 * Proj gradu by L2
	 *
	 *
	 * ------------------------------------------------------------ */

	/* L2 projection of gradu */
	Gradu = phgDofNew(g, DOF_P1, DDim, "Gradu", DofNoAction);
	Gradu0 = phgDofNew(g, DOF_P1, 1, "Gradu", DofNoAction);

	/* solver_Gu */
	phgOptionsPush();
	phgOptionsSetOptions(ns_params->proj_opts);
	solver_Gu = phgSolverCreate(SOLVER_DEFAULT, Gradu0, NULL);
	solver_Gu->verb = 0;
	phgOptionsPop();

	VEC *vec[DDim];
	for (k = 0; k < DDim; k++) {
	    if (components[k]) {
		vec[k] = phgMapCreateVec(solver_Gu->rhs->map, 1) ;
		phgVecDisassemble(vec[k]);
	    }
	}

	/* Build linear system */
	if (0) {
	    /* Removed */
	    /* ns_params->utype == DOF_P1 */
	    /* && ns_params->use_prism_elem) { */
	    /* buildProjGuPrism(ns, solver_Gu, &vec[0], components); */
	}
	else {
	    ForAllElements(g, e) {
		int M = Gradu0->type->nbas;	/* num of bases of Velocity */
		int order = DofTypeOrder(Gradu0, e) * 2;
		FLOAT A[M][M], rhs[DDim][M];
		INT I_[M];
		QUAD *quad;
		FLOAT vol, det;
		const FLOAT *w, *p, *gu;

		Bzero(A); Bzero(rhs);
		quad = phgQuadGetQuad3D(order);
		gu = phgQuadGetDofValues(e, gradu, quad); 

		p = quad->points;
		w = quad->weights;
		for (q = 0; q < quad->npoints; q++) {
		    phgGeomGetCurvedJacobianAtLambda(g, e, p, &det);
		    vol = fabs(det / 6.);

		    for (i = 0; i < M; i++) {
			const FLOAT *gi = phgQuadGetBasisValues(e, Gradu0, i, quad) + q;    
			for (j = 0; j < M; j++) {
			    const FLOAT *gj = phgQuadGetBasisValues(e, Gradu0, j, quad) + q;       
			    FLOAT qmass = vol*(*w) * (*gj) * (*gi);
			    A[i][j] += qmass;
			}
		    
			for (k = 0; k < DDim; k++) {
			    if (components[k]) {
				rhs[k][i] += vol*(*w) * gu[k] * (*gi); 
			    }
			}
		    }
		    gu += DDim;
		    w++; p += Dim + 1;
		}

		/* Map: Element -> system */
		for (i = 0; i < M; i++)
		    I_[i] = phgMapE2L(solver_Gu->mat->cmap, 0, e, i);

		/* Global res */
		for (i = 0; i < M; i++)
		    phgMatAddEntries(solver_Gu->mat, 1, I_ + i, M, I_,
				     &(A[i][0])); 

		for (k = 0; k < DDim; k++)
		    if (components[k]) 
			phgVecAddEntries(vec[k], 0, M, I_, &rhs[k][0]);
	    }/* end element */
	}
	
    

	for (k = 0; k < DDim; k++)
	    if (components[k]) 
		phgVecAssemble(vec[k]);
	solver_Gu->rhs_updated = TRUE;

	INT n = DofGetDataCount(Gradu0);
	for (k = 0; k < DDim; k++) {
	    if (components[k]) {
		phgVecCopy(vec[k], &solver_Gu->rhs);
		phgDofSetDataByValue(Gradu0, 0.);
		phgSolverSolve(solver_Gu, FALSE, Gradu0, NULL);
		phgPrintf("      solver_Gu: nits = %d, resid = %0.4lg ",
			  solver_Gu->nits, solver_Gu->residual);
		elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));

		FLOAT *vg = Gradu->data, *vg0 = Gradu0->data;
		for (i = 0; i < n; i++, vg0++, vg += DDim)
		    vg[k] = *vg0;
	    }
	}

	if (DUMP_MAT_VEC) {
	    phgPrintf("Dumping MatGu, rhsGu\n");
	    phgMatDumpMATLAB(solver_Gu->mat, "A_gu", "mat_gu_.m");
	    phgVecDumpMATLAB(solver_Gu->rhs, "b_gu", "rhs_gu_.m");
	}

	phgInfo(2, "   Destroy solver Gu\n");
	for (k = 0; k < DDim; k++)
	    if (components[k]) 
		phgVecDestroy(&vec[k]);
	phgSolverDestroy(&solver_Gu);
	phgDofFree(&Gradu0);

	phgExportTecplot(g, "gradu0.plt", Gradu, NULL);
    } 
    else {
	/* ------------------------------------------------------------
	 *
	 *
	 * Proj gradu by average
	 *
	 *
	 * ------------------------------------------------------------ */
	Gradu = phgDofNew(g, DOF_P1, DDim, "Gradu", DofNoAction);
	DOF *scale = phgDofNew(g, DOF_P1, 1, "scale", DofNoAction);
	phgDofSetDataByValue(scale, 0.);

	MAP *map0 = phgMapCreate(scale, NULL);
	VEC *vec0 = phgMapCreateVec(map0, 1);

	VEC *vec[DDim];
	phgVecDisassemble(vec0);
	for (k = 0; k < DDim; k++) {
	    if (components[k]) {
		vec[k] = phgMapCreateVec(map0, 1) ;
		phgVecDisassemble(vec[k]);
	    }
	}
    
	ForAllElements(g, e) {
	    int M = scale->type->nbas;	/* num of bases of Velocity */	
	    INT I_[M];
	    FLOAT a[M], rhs[DDim][M], vol, gu[DDim];
	    const FLOAT lambda0[4] = {.25, .25, .25, .25};

#if 1	
	    vol = phgGeomGetVolume(g, e);
#else
	    vol = 1.; 
#endif

	    Bzero(a); Bzero(rhs);
	    for (i = 0; i < M; i++) {
		a[i] += vol;		/* scale */

		phgDofEval(gradu, e, lambda0, gu);
		for (k = 0; k < DDim; k++) {
		    if (components[k]) {
			rhs[k][i] += vol * gu[k];
		    }
		}
	    }

	    for (i = 0; i < M; i++) 
		I_[i] = phgMapE2L(map0, 0, e, i);

	    phgVecAddEntries(vec0, 0, M, I_, &a[0]);
	    for (k = 0; k < DDim; k++)
		if (components[k]) 
		    phgVecAddEntries(vec[k], 0, M, I_, &rhs[k][0]);
	}

	/* assemble */
	phgVecAssemble(vec0);
	for (k = 0; k < DDim; k++)
	    if (components[k]) 
		phgVecAssemble(vec[k]);

	/* average */
	for (k = 0; k < DDim; k++)
	    if (components[k]) {
		const FLOAT *s = vec0->data;
		FLOAT *v = vec[k]->data;

		for (i = 0; i < map0->nlocal; i++) {
		    v[i] /= s[i];
		}
	    }

	/* assign dofs */
	for (k = 0; k < DDim; k++)
	    if (components[k]) {
		phgMapLocalDataToDof(map0, 1, &scale, vec[k]->data);
		FLOAT *v0 = scale->data;
		FLOAT *v = Gradu->data;

		for (i = 0; i < g->nvert; i++) {
		    v[k] = *v0;

		    v0++;
		    v += DDim;
		}
	    }    
    

	phgVecDestroy(&vec0);
	for (k = 0; k < DDim; k++)
	    if (components[k]) 
		phgVecDestroy(&vec[k]);
	phgMapDestroy(&map0);
	phgDofFree(&scale);

	phgExportTecplot(g, "gradu1.plt", Gradu, NULL);
    }



    if (ns->Gradu == NULL)
	ns->Gradu = Gradu;
    else {
	phgDofCopy(Gradu, &ns->Gradu, NULL, "Gradu");
	phgDofFree(&Gradu);
    }
}




/*
 * Output bottom temperature
 *
 *  */
void 
output_bottom_dofs(NSSolver *ns, int tstep)
{
    GRID *g = ns->g;
    SIMPLEX *e;

    phgPrintf("Output bottom temp.\n");

    DOF *temp = phgDofCopy(ns->T[1], NULL, DOF_P3,"temp");
    DOF *coord = phgDofNew(g, temp->type, Dim, "coord", func_xyz_);
    int ndof = DofGetDataCount(temp);
    int ndof_global = DofGetDataCountGlobal(temp);
    MAP *map = phgMapCreate(temp, NULL);

    /* send buf:
     * 
     * 1. Gidx    1
     * 2. Coord   3
     * 3. temp    1
     *
     * */
    FLOAT *sbuf, *p;		
    INT *mark;
    sbuf = phgCalloc(5 * ndof,  sizeof(*sbuf));
    mark = phgCalloc(ndof, sizeof(*mark));

    INT ndat = 0;
    p = sbuf;

    ForAllElements(g, e) {
	int N = temp->type->nbas;
	INT i, s, iD, iG;
	int nbas_face = NbasFace(temp);
	SHORT bases[nbas_face];

	for (s = 0; s < NFace; s++) {
	    if (e->bound_type[s] & BC_BOTTOM) {
		phgDofGetBasesOnFace(temp, e, s, bases);

		for (i = 0; i < nbas_face; i++) {
		    int ii = bases[i];

		    /* Dof idx */
		    iD = phgDofMapE2D(temp, e, ii);
		    iG = phgMapE2G(map, 0, e, ii);
		    
		    if (mark[iD] == 0) {

			p[0] = iG;
			p[1] = coord->data[iD*3 + 0];
			p[2] = coord->data[iD*3 + 1];
			p[3] = coord->data[iD*3 + 2];
			p[4] = temp->data[iD];
			p += 5;
			ndat += 5;

			mark[iD] = 1;
		    }
		    
		}	
	    } /* bottom face */
	} /* end of face */
    } /* end of elems */
    
    /* Gather to rank0 */
    INT *counts, *displs;
    counts = phgAlloc(2 * g->nprocs * sizeof(*counts));
    displs = counts + g->nprocs;


    MPI_Allgather(&ndat, 1, MPI_INT, counts, 1, MPI_INT, g->comm);

    INT i, m = 0, ntotal;
    for (i = 0; i < g->nprocs; i++) {
        displs[i] = m;
        m += counts[i];
    }
    ntotal = m;

    FLOAT *rbuf = NULL;
    if (g->rank == 0)
	rbuf = phgAlloc(ntotal * sizeof(*rbuf));

    MPI_Gatherv(sbuf, ndat, PHG_MPI_FLOAT,
		rbuf, counts, displs, PHG_MPI_FLOAT, 0, g->comm);
    
    if (g->rank == 0) {
	char fname[1000];
	sprintf(fname, OUTPUT_DIR "/ins_" NS_PROBLEM "_bot_%05d.dat", tstep - 1);

	FILE *fp;
	if ((fp = fopen(fname, "w")) == NULL) {
	    phgError(0, "%s can not be open!!!\n", fname);
	} 

	p = rbuf;
	for (i = 0; i < ntotal / 5; i++, p+=5) {
	    fprintf(fp, "%5d %20.10e %20.10e %20.10e %20.10e\n", 
		   (int) p[0], p[1], p[2], p[3], p[4]);
	}

	fclose(fp);
	phgFree(rbuf);
    }

    phgFree(counts);
    phgFree(sbuf);
    phgFree(mark);

    phgMapDestroy(&map);
    phgDofFree(&temp);
    phgDofFree(&coord);
}












