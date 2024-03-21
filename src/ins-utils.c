/*
 *  Usefull subroutines.
 *
 *
 *
 *  */
#include "ins.h"
#include "mat_op3.h"
#if USE_PETSC
#  include <petscsys.h>
#  include <petscviewer.h>
#endif


/* Segv */
#include <signal.h>
#include <execinfo.h> // backtrace(), etc
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdio.h>


/* Print Macro defination */
void 
NsSolver_Options() 
{
    /* Information of solver */
    phgPrintf("\n===============================================\n");
    phgPrintf(" Ice sheet solver options:  \n\n");
    phgPrintf("===============================================\n\n");

    phgPrintf("Scaling:\n");
    phgPrintf("   equ : %e\n", EQU_SCALING);
    phgPrintf("   len : %e\n", LEN_SCALING);
    phgPrintf("   pres: %e\n", PRES_SCALING);

    if (ns_params->use_slide) {
	phgPrintf("Sliding BC.\n");
    }
    else {
	phgPrintf("No sliding BC.\n");
    }

    return;
}



static
void handler(int sig) {
    void *array[100];
    size_t size;

    // get void*'s for all entries on the stack
    size = backtrace(array, 100);

    // print out all the frames to stderr
    fprintf(stderr, "Error: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, STDERR_FILENO);

    MPI_Finalize();             /* Cause other procs waiting,
                                 * then manually interupte mpirun.
                                 *  */

    fprintf(stderr, "MPI finalized.\n");
    exit(1);
}

void install_error_log()
{
    /* Redirect stderr */
    {
	char name[1000];
	FILE *fp;

	sprintf(name, "err/%04d.log", phgRank);
	fp = fopen(name, "w");
	if (fp == NULL) {
	    phgError(0, "Error creating error file: %s, check that the directory ``err'' has been created.\n", name);
	}
	assert(fp != NULL);
	dup2(fileno(fp), fileno(stderr));
	fclose(fp);
	fprintf(stderr, "redirect stderr.\n");
    }


    // install segv handler
    signal(SIGSEGV, handler);
    signal(SIGABRT, handler);
    signal(SIGINT, handler);
    //signal(SIGFPE, handler);  // 1/0 is not trapped
    signal(SIGTERM, handler);
}



double
elapsed_time(GRID *g, BOOLEAN flag, double mflops)
/* returns and prints elapsed time since last call to this function */
{
    static double t0 = 0.0;
    double et, tt[3];
    size_t mem;

    phgGetTime(tt);
    et = tt[2] - t0;
    t0 = tt[2];

    mem = phgMemoryUsage(g, NULL);

    if (flag) {
	if (mflops <= 0)
	    phgPrintf("[%0.4lgMB %0.4lfs]\n", mem / (1024.0 * 1024.0), et);
	else
	    phgPrintf("[%0.4lgMB %0.4lfs %0.4lgGF]\n",
		      mem / (1024.0 * 1024.0), et, mflops * 1e-3);
    }

    return et;
}




/* --------------------------------------------------------------------------------
 *
 *  Solution utils 
 *
 *
 * -------------------------------------------------------------------------------- */

/* Dof norm for curved elements */
FLOAT
dofNormL2(DOF *dof)
{
    GRID *g = dof->g;
    SIMPLEX *e;
    FLOAT norm = 0;
    FLOAT a = 0, b = 0;
    int i, dim = DofDim(dof);
    
    ForAllElements(g, e) {
	QUAD *quad = phgQuadGetQuad3D(2*DofTypeOrder(dof, e));
	FLOAT vol, det, v[dim];
	const FLOAT *p, *w;
	int q;
	
	p = quad->points;
 	w = quad->weights;
	for (q = 0; q < quad->npoints; q++) {
	    phgGeomGetCurvedJacobianAtLambda(g, e, p, &det);
	    vol = fabs(det / 6.);
	    phgDofEval(dof, e, p, v);
	    for (i = 0; i < dim; i++)
		norm += vol*(*w) * (v[i] * v[i]);
	    w++; p += Dim+1;
	}
    }

    a = norm;
    MPI_Allreduce(&a, &b, 1, MPI_DOUBLE, MPI_SUM, g->comm);
    norm = sqrt(b);
    
    return norm;
}

/* Dof norm for curved elements */
FLOAT
dofDiffNormL2(DOF *dof1, DOF *dof2)
{
    GRID *g = dof1->g;
    SIMPLEX *e;
    FLOAT norm = 0;
    FLOAT a = 0, b = 0;
    int i, dim = DofDim(dof1);
    
    ForAllElements(g, e) {
	QUAD *quad = phgQuadGetQuad3D(2*DofTypeOrder(dof1, e));
	FLOAT vol, det, v1[dim], v2[dim];
	const FLOAT *p, *w;
	int q;
	
	p = quad->points;
 	w = quad->weights;
	for (q = 0; q < quad->npoints; q++) {
	    phgGeomGetCurvedJacobianAtLambda(g, e, p, &det);
	    vol = fabs(det / 6.);
	    phgDofEval(dof1, e, p, v1);
	    phgDofEval(dof2, e, p, v2);
	    for (i = 0; i < dim; i++)
		norm += vol*(*w) * (v1[i] - v2[i]) * (v1[i] - v2[i]);
	    w++; p += Dim+1;
	}
    }

    a = norm;
    MPI_Allreduce(&a, &b, 1, MPI_DOUBLE, MPI_SUM, g->comm);
    norm = sqrt(b);
    
    return norm;
}


#define MAX_QUAD_ORDER 12
void
dof_norm_L2(DOF *dof)
{
    GRID *g = dof->g;
    SIMPLEX *e;
    int i, dim = DofDim(dof);
    FLOAT a[dim], b[dim], norm[dim];

    Bzero(norm);
    ForAllElements(g, e) {
	QUAD *quad;
	FLOAT vol, det, v[dim];
	const FLOAT *p, *w;
	int q, order = DofTypeOrder(dof, e) * 2;
	
	if (order > MAX_QUAD_ORDER)
	    order = MAX_QUAD_ORDER;
	quad = phgQuadGetQuad3D(order);
	p = quad->points;
 	w = quad->weights;
	for (q = 0; q < quad->npoints; q++) {
	    phgGeomGetCurvedJacobianAtLambda(g, e, p, &det);
	    vol = fabs(det / 6.);
	    phgDofEval(dof, e, p, v);
	    for (i = 0; i < dim; i++)
		norm[i] += vol*(*w) * (v[i] * v[i]);
	    w++; p += Dim+1;
	}
    }

    memcpy(a, norm, dim * sizeof(FLOAT));
    MPI_Allreduce(&a, &b, dim, MPI_DOUBLE, MPI_SUM, g->comm);

    phgPrintf("   Dof %s normL2\n", dof->name);
    for (i = 0; i < dim; i++) {
	norm[i] = sqrt(b[i]);
	phgPrintf("      %12.4E\n", norm[i]);
    }
    
    return;
}

void dof_range(DOF *u)
{
    GRID *g = u->g;
    int dim = u->dim;
    INT i, k, N;
    FLOAT max[dim*2], min[dim], tmp[2*dim];
    const FLOAT *v;

    for (k = 0; k < dim; k++) {
	max[k] = -1e32;
	min[k] = 1e32;
    }
    
    v = DofData(u);
    N = DofGetDataCount(u);
    for (i = 0; i < N / dim; i++) {
	for (k = 0; k < dim; k++) {
	    max[k] = MAX(max[k], *v);
	    min[k] = MIN(min[k], *v);
	    v++;
	}
    }

    for (k = 0; k < dim; k++) 
	max[k + dim] = - min[k];
    memcpy(tmp, max, 2 * dim * sizeof(FLOAT));
    MPI_Allreduce(&tmp, &max, 2 * dim, MPI_DOUBLE, MPI_MAX, g->comm);
    for (k = 0; k < dim; k++) 
	min[k]  = - max[k + dim];
    

    if (dim == 1) {
 	phgPrintf("   %-5s: [%24.12e, %24.12e]\n", 
		  u->name, min[0], max[0]); 
    } else {
	for (k = 0; k < dim; k++) {
	    if (k == 0)
		phgPrintf("   %-5s[0]: [%24.12e, %24.12e]\n", 
			  u->name, min[k], max[k]); 
	    else
		phgPrintf("          [%d]: [%24.12e, %24.12e]\n", 
			  k, min[k], max[k]); 
	}
    }
}



/* Symetric check:
 * Use random vector x, y check whether Mat and PC is symetric,
 *   if they are, then x'*A*y == y'*A*x.
 *   This is NOT a good way, since rounding error may be not small.
 * */
void sym_check(NSSolver *ns)
{
#if 0
    if (0 && ns_params->use_symetric) {
	VEC *xup = phgMapCreateVec(ns->solver_u->rhs->map, 1);
	VEC *yup, *zup;
	FLOAT product[2], p0, test;
		
	yup = phgVecCopy(xup, NULL);
	zup = phgVecCopy(xup, NULL);
	phgVecRandomize(xup, 0);
	phgVecRandomize(yup, 0);
		
	phgMatVec(MAT_OP_N, 1.0, ns->solver_u->mat, xup, 0., &zup);
	phgVecDot(yup, 0, zup, 0, product);
	phgMatVec(MAT_OP_N, 1.0, ns->solver_u->mat, yup, 0., &zup);
	phgVecDot(xup, 0, zup, 0, product+1);
		
	if (fabs(product[0]) > 1e-10) 
	    p0 = product[0];
	else if (fabs(product[1]) > 1e-10) 
	    p0 = product[1];
	else
	    p0 = 1.;
	test = fabs((product[0] - product[1]) / p0);
	if (test > 1e-10) {
	    phgError(-1, "Solver matrix may be unsymetric!!!\n random vec test:(%E, %E), diff:%E\n", 
		     product[0], product[1], test);
	}
		
	if (ns_params->use_PCD) {
	    pc_proc(ns->solver_u, xup, &zup);
	    phgVecDot(yup, 0, zup, 0, product);
	    pc_proc(ns->solver_u, yup, &zup);
	    phgVecDot(xup, 0, zup, 0, product+1);
		
	    if (fabs(product[0]) > 1e-10) 
		p0 = product[0];
	    else if (fabs(product[1]) > 1e-10) 
		p0 = product[1];
	    else
		p0 = 1.;
	    test = fabs((product[0] - product[1]) / p0);
	    if (test > 1e-10) {
		phgPrintf("Dumping sym mat vec\n");
		printf("p:%E, %E\n", product[0], product[1]);
		phgVecDumpMATLAB(xup, "vx", "vx_.m");
		phgVecDumpMATLAB(yup, "vy", "vy_.m");
		phgVecDumpMATLAB(ns->pcd->rhsScale, "vs", "vs_.m");
		phgMatDumpMATLAB(ns->matF, "F", "F_.m");
		phgMatDumpMATLAB(ns->pcd->matQp, "Qp", "Qp_.m");
		phgMatDumpMATLAB(ns->pcd->matQp, "Qp", "Qp_.m");
		phgMatDumpMATLAB(ns->pcd->matAp, "Ap", "Ap_.m");
		phgMatDumpMATLAB(ns->pcd->matFp, "Fp", "Fp_.m");
		if (ns_params->use_symetric)
		    phgMatDumpMATLAB(ns->pcd->matFu, "Fu", "Fu_.m");
		phgError(-1, "PC matrix may be unsymetric!!!\n random vec test:(%E, %E), diff:%E\n", 
			 product[0], product[1], test);
	    }
	}
    }
#endif
}


/***************************/
/* Check all kinds of bdry */
/***************************/

void
checkBdry(GRID *g)
{
    SIMPLEX *e;
    int s;
    double a[8];

    Bzero(a);
    ForAllElements(g, e) {
        a[7] += phgGeomGetVolume(g, e);
        for (s = 0; s < NFace; s++) {
            FLOAT area = phgGeomGetFaceArea(g, e, s);
            if (e->bound_type[s] & BDRY_MASK) {
                if (e->bound_type[s] & BC_TOP)
                    a[0] += area;
                else if (e->bound_type[s] & BC_BOTTOM)
                    a[1] += area;
                else if (e->bound_type[s] & BC_LATERAL)
                    a[2] += area;
                else if (e->bound_type[s] & BC_DIVIDE)
                    a[3] += area;
                else if (e->bound_type[s] & BC_FRONT)
                    a[4] += area;
		else
		    a[6] += area;
            }
            else {
                a[7] += area;
            }
        }

	if (isnan(a[7])) {
	    phgPrintf("voL: %d %e\n", e->index, phgGeomGetVolume(g, e));
	    phgDumpElement(g, e);
	    phgError(1, "Grid problem.\n");
	}
    }

    phgInfo(0, "\tTop     : %20.10e\n"
	    "\tBottom  : %20.10e\n"
	    "\tLateral : %20.10e\n"
	    "\tDivide  : %20.10e\n"
	    "\tFront   : %20.10e\n"
	    "\tOther   : %20.10e\n",
	    a[0], a[1], a[2], a[3], a[4], a[6]);


#if USE_MPI
    {
        double b[8];
        MPI_Reduce(a, b, 8, MPI_DOUBLE, MPI_SUM, 0, g->comm);
        memcpy(a, b, sizeof(b));
    }
#endif

    phgPrintf("Boundary types check:\n");
    phgPrintf( "\tTop     : %20.10e\n"
	       "\tBottom  : %20.10e\n"
	       "\tLateral : %20.10e\n"
	       "\tDivide  : %20.10e\n"
	       "\tFront   : %20.10e\n"
	       "\tOther   : %20.10e\n",
	       a[0], a[1], a[2], a[3], a[4], a[6]);


    return;
}


void 
getPecletNum(GRID *g, DOF *u, FLOAT nu, int order)	
/* 
 * Mesh Peclet number:
 * Pe := U * h / nu,
 * where, U = ||u||, h = diam of circumsphere of tetra as h.
 * This suggest that the wind direction coincides with the direction
 * of max length of the element.
 * 
 * */
{
    SIMPLEX *e;
    int q;
    FLOAT pe, pe0, diam, maxu, tmp;

    assert(u->dim == Dim);
    pe0 = -1;
    ForAllElements(g, e) {
	QUAD *quad;
	const FLOAT *vu;
	
	diam = phgGeomGetDiameter(g, e);
	quad = phgQuadGetQuad3D(order);
	vu = phgQuadGetDofValues(e, u, quad);
	maxu = -1;
	for (q = 0; q < quad->npoints; q++) {
	    tmp = INNER_PRODUCT(vu, vu);
	    maxu = MAX(maxu, tmp);
	    vu += Dim;
	}
	tmp = sqrt(maxu) * diam;
	pe0 = MAX(pe0, tmp);
    }

    MPI_Allreduce(&pe0, &pe, 1, MPI_DOUBLE, MPI_MAX, g->comm);
    pe /= nu;
    phgPrintf("* Mesh Peclet num: %16.8E\n", pe);
    if (pe > 2.)
	phgPrintf("* Pe num %16.8E > 2, FEM may get unstable!\n", pe);

    return;
}






void 
get_p_static(NSSolver *ns)
{
    GRID *g = ns->g;
    LAYERED_MESH *gL = ns->gL;
    int i, j, k;
    DOF *dofp_static;

    /* Removed */
    /* if (ns_params->use_prism_elem) { */
    /* 	get_p_static_prism(ns); */
    /* 	return; */
    /* } */

    if (ns_params->stab_remove_static == 3) {
	dofp_static = phgDofNew(g, DOF_P2, 1, "dofp_static", DofNoAction);
	phgError(0, "Unimplemented!!!\n");
    }
    else {
	dofp_static = phgDofNew(g, DOF_P1, 1, "dofp_static", DofNoAction);
    }
	
    for (i = 0; i < g->nvert; i++) {
	if (g->types_vert[i] & BC_BOTTOM) {
	    assert(gL->vert_local_lists[i] != NULL);

	    FLOAT ztop;
	    int nv = gL->vert_local_lists[i][0];
	    int *iL = &gL->vert_local_lists[i][1];

	    /* Pressure SIA:
	     *  p =  rho g ( s - z )
	     * */
	    ztop = g->verts[iL[nv-1]][Z_DIR];
	    for (j = 0; j < nv; j++) {
		const FLOAT rho = RHO_ICE;
		const FLOAT grav = GRAVITY;
#if 1
		*DofVertexData(dofp_static, iL[j]) =
		    (rho * grav * (ztop - g->verts[iL[j]][Z_DIR]) * LEN_SCALING) / PRES_SCALING;
#else
#  warning ----- GLS only pressure -----		
		/* GLS: only (0, 0, -rho g z) */
		*DofVertexData(dofp_static, iL[j]) =
		    (rho * grav * (0 - g->verts[iL[j]][Z_DIR]) * LEN_SCALING) / PRES_SCALING;
#endif

		phgInfo(3, "Vert: %10.5e %10.5e %10.5e: %10.5e %10.5e\n",
			g->verts[iL[j]][0],
			g->verts[iL[j]][1],
			g->verts[iL[j]][2],
			ztop, *DofVertexData(dofp_static, iL[j]));
	    }
	}
    }

    
    if (ns->p_static == NULL) 
	ns->p_static = dofp_static;
    else 
	phgDofCopy(dofp_static, &ns->p_static, NULL, "p_static");
	
    return;
}



/* projection of p0 to p1 */
DOF *
proj_G1(NSSolver *ns, DOF *p0, DOF **p1_ptr)
{
    GRID *g = ns->g;
    SIMPLEX *e;
    int i, j, q;
    DOF *p1 = *p1_ptr;
    VEC *vec_rhs;
    static SOLVER *solver_p1 = NULL;
    static BOOLEAN init = TRUE;
    
    /* L2 projection of pressure */
    //p1 = phgDofNew(g, DOF_P1, 1, "pp1", DofNoAction);

    /* solver_p */
    if (init) {
	phgOptionsPush();
	phgOptionsSetOptions("-solver gmres "
			     "-solver_maxit 100 "
			     "-solver_rtol 1e-10");
	//phgOptionsSetOptions(Gu_opts);
	solver_p1 = phgSolverCreate(SOLVER_DEFAULT, p1, NULL);
	solver_p1->verb = 2;
	phgOptionsPop();
    }

    vec_rhs = phgMapCreateVec(solver_p1->rhs->map, 1);
    phgVecDisassemble(vec_rhs);

    /* Build linear system */
    ForAllElements(g, e) {
	int M = p1->type->nbas;	/* num of bases of Velocity */
	int order = DofTypeOrder(p1, e) * 2;
	FLOAT A[M][M], rhs[M];
	INT I_[M];
	QUAD *quad;
	FLOAT vol, det;
	const FLOAT *w, *p, *vp;

	Bzero(A); Bzero(rhs);
	quad = phgQuadGetQuad3D(order);
	vp = phgQuadGetDofValues(e, p0, quad); 

	p = quad->points;
	w = quad->weights;
	for (q = 0; q < quad->npoints; q++) {
	    phgGeomGetCurvedJacobianAtLambda(g, e, p, &det);
	    vol = fabs(det / 6.);

	    for (i = 0; i < M; i++) {
		const FLOAT *gi = phgQuadGetBasisValues(e, p1, i, quad) + q;    

		if (init) {
		    for (j = 0; j < M; j++) {
			const FLOAT *gj = phgQuadGetBasisValues(e, p1, j, quad) + q;       
			FLOAT qmass = vol*(*w) * (*gj) * (*gi);
			A[i][j] += qmass;
		    }
		}
		    
		rhs[i] += vol*(*w) * (*vp) * (*gi); 
	    }
	    vp++;
	    w++; p += Dim + 1;
	}

	/* Map: Element -> system */
	for (i = 0; i < M; i++)
	    I_[i] = phgMapE2L(solver_p1->mat->cmap, 0, e, i);

	/* Global res */
	if (init)
	    for (i = 0; i < M; i++)
		phgMatAddEntries(solver_p1->mat, 1, I_ + i, M, I_,
				 &(A[i][0])); 

	phgVecAddEntries(vec_rhs, 0, M, I_, &rhs[0]);
    }				/* end element */

    phgVecAssemble(vec_rhs);
    phgVecAssemble(solver_p1->rhs);
    phgVecAXPBY(1., vec_rhs, 0, &solver_p1->rhs);
    phgVecDestroy(&vec_rhs);
    solver_p1->rhs_updated = FALSE;

    phgDofSetDataByValue(p1, 0.);
    phgSolverSolve(solver_p1, FALSE, p1, NULL);
    phgPrintf("      solver_p1: nits = %d, resid = %0.4lg\n",
	      solver_p1->nits, solver_p1->residual);

    init = FALSE;


    return NULL;
}



void
check_div(DOF *gradu, DOF **divu_ptr)
{
    GRID *g = gradu->g;
    SIMPLEX *e;
    int i;

    /*
     * div0 is a P0 DOF,
     * in each element div0(e) = \int_e div(u) / \int_e.
     *
     * */
    //DOF *div0 = phgDofNew(g, DOF_P0, 1, "div0", DofNoAction);
    DOF *div0 = *divu_ptr;

    ForAllElements(g, e) {
	QUAD *quad;
	FLOAT vol, det, dive;
	const FLOAT *p, *w, *gu;
	int q, order = DofTypeOrder(gradu, e) * 2;
	
	dive = 0;
	quad = phgQuadGetQuad3D(order);
	gu = phgQuadGetDofValues(e, gradu, quad);  /* grad u */

	p = quad->points;
 	w = quad->weights;
	for (q = 0; q < quad->npoints; q++) {
	    phgGeomGetCurvedJacobianAtLambda(g, e, p, &det);
	    vol = fabs(det / 6.);

	    dive += (*w) * (gu[0] + gu[4] + gu[8]);
	    w++; p += Dim+1;
	    gu += DDim;
	}
	*DofElementData(div0, e->index) = dive;
    }

    DOF_SCALE(div0, "div0");
    //phgDofFree(&div0);

    return;
}




void
plot_residual(NSSolver *ns, VEC *residual, int non_step)
{
    GRID *g = ns->g;
    SIMPLEX *e;
    int i, j, k, q;
    static DOF *resu, *resp;
    static SOLVER *solver_u = NULL, *solver_p = NULL;

    if (solver_u == NULL) {
	resu = phgDofCopy(ns->du, NULL, NULL, "resu");
	resp = phgDofCopy(ns->dp, NULL, NULL, "resp");
	resu->DB_masks[0] = 0;
	resu->DB_masks[1] = 0;
	resu->DB_masks[2] = 0;
	resp->DB_mask = 0;

	solver_u = phgSolverCreate(SOLVER_DEFAULT, resu, NULL);
	solver_p = phgSolverCreate(SOLVER_DEFAULT, resp, NULL);
	solver_u->verb = 0;
	solver_p->verb = 0;
    
	ForAllElements(g, e) {
	    int M = resu->type->nbas;	/* num of bases of Velocity */
	    int N = resp->type->nbas;	/* num of bases of Pressure */
	    int order = DofTypeOrder(resu, e) * 2;
	    FLOAT Mu[M][Dim][M][Dim], Mp[N][N];
	    QUAD *quad;
	    FLOAT vol, det;
	    INT Iu[M][Dim], Ju[Dim][M], Ip[N];
	    const FLOAT *w, *p;

	    for (i = 0; i < M; i++)
		for (k = 0; k < Dim; k++)
		    Ju[k][i] = Iu[i][k] = phgMapE2L(ns->matF->cmap, 0, e, i * Dim + k);
	    for (i = 0; i < N; i++)
		Ip[i] = phgMapE2L(ns->matC->cmap, 0, e, i);

	    quad = phgQuadGetQuad3D(order);
	
	    Bzero(Mu); Bzero(Mp); 
	    p = quad->points;
	    w = quad->weights;
	    for (q = 0; q < quad->npoints; q++) {
		phgGeomGetCurvedJacobianAtLambda(g, e, p, &det);
		vol = fabs(det / 6.);

		for (i = 0; i < M; i++) {
		    const FLOAT *gi_u = phgQuadGetBasisValues(e, resu, i, quad) + q;       /* phi_i */
		    for (j = 0; j < M; j++) {
			const FLOAT *gj_u = phgQuadGetBasisValues(e, resu, j, quad) + q;       /* phi_j */

			for (k = 0; k < Dim; k++) 
			    Mu[i][k][j][k] += vol*(*w) * (*gi_u)*(*gj_u);
		    }
		}

		for (i = 0; i < N; i++) {
		    const FLOAT *gi_p = phgQuadGetBasisValues(e, resp, i, quad) + q;       /* phi_i */
		    for (j = 0; j < N; j++) {
			const FLOAT *gj_p = phgQuadGetBasisValues(e, resp, j, quad) + q;       /* phi_j */

			Mp[i][j] += vol*(*w) * (*gi_p)*(*gj_p);
		    }
		}
		w++; p += Dim+1;
	    } /* end quad pts */
	

	    for (i = 0; i < M; i++) 
		for (k = 0; k < Dim; k++) 
		    phgMatAddEntries(solver_u->mat, 1, Iu[i] + k, M*Dim, Iu[0],
				     &(Mu[i][k][0][0]));
	    for (i = 0; i < N; i++) 
		phgMatAddEntries(solver_p->mat, 1, Ip + i, N, Ip,
				 &Mp[i][0]);
	
	} /* end elements */

	phgVecAssemble(solver_u->rhs);
	phgVecAssemble(solver_p->rhs);
	solver_u->rhs_updated = FALSE;
	solver_p->rhs_updated = FALSE;
    }

    INT nu = solver_u->rhs->map->nlocal,
	np = solver_p->rhs->map->nlocal;

    memcpy(solver_u->rhs->data, residual->data, nu * sizeof(FLOAT));
    memcpy(solver_p->rhs->data, residual->data + nu, np * sizeof(FLOAT));

    phgSolverSolve(solver_u, FALSE, resu, NULL);
    phgSolverSolve(solver_p, FALSE, resp, NULL);

    char name[100];
    /* sprintf(name, "res_%02d.plt", non_step); */
    /* phgExportTecplot(g, name, ns->u[1], ns->p[1], resu, resp, NULL); */
    sprintf(name, "res_%02d.vtk", non_step);
    phgExportVTK(g, name, ns->u[1], ns->p[1], resu, resp, NULL);

    /* phgDofFree(&resu); */
    /* phgDofFree(&resp); */

    phgFinalize();
    exit(1);
    return;
}


void
plot_surf_force(NSSolver *ns)
/* Plot surface force
 *   f = p n - n \dot \tau 
 * */
{
    GRID *g = ns->g;
    INT i, j;
    SURF_BAS *surf_bas = ns->surf_bas;
    DOF *surf_dof = surf_bas->dof;
    DOF *sforce_dof = phgDofNew(g, DOF_P1, 1, "sforce", DofNoAction);

    phgDofSetDataByValue(sforce_dof, 0.);
    
    proj_gradu(ns, ns->gradu[1], NULL, 0);


    FLOAT *dat_Gu = ns->Gradu->data;
    FLOAT *dat_p  = ns->p[1]->data;
    FLOAT *dat_f  = sforce_dof->data;
    for (i = 0; i < g->nvert; i++) {
	const FLOAT *normal, *gu;
	FLOAT nu, p, en[Dim], f;

	if (g->types_vert[i] & (BC_BOTTOM | BC_LATERAL)) { 

	    if (g->types_vert[i] & BC_LATERAL)
		normal = surf_dof->data + i * (Dim*Dim) + LN_DIR * Dim;

	    if (g->types_vert[i] & BC_BOTTOM)
		normal = surf_dof->data + i * (Dim*Dim) + UN_DIR * Dim;

	    /* The case of bottom and latereal uses bottom normal */

	    gu = dat_Gu + i * (Dim*Dim);
	    nu = get_effective_viscosity(gu, 0, 0, ns->viscosity_type);
	    p = dat_p[i];

	    FLOAT eu[DDim];
	    /* eu = .5 * gu + .5 * gu^T */
	    MAT3_SYM(gu, eu);

	    /* eu = eu / 1000 */
	    MAT3_AXPBY(0., gu, 1./LEN_SCALING, eu);

	    /* eu * n */
	    MAT3_MV(eu, normal, en);
	
	    dat_f[i] = p * PRES_SCALING - nu * INNER_PRODUCT(en, normal);
	}
	else if (g->types_vert[i] & BC_FRONT) {
	    /*
	     * Note: ONLY for MISMIP front !!!
	     * */
	    const FLOAT x_dir[Dim] = {1, 0, 0};
	    normal = x_dir;

	    gu = dat_Gu + i * (Dim*Dim);
	    nu = get_effective_viscosity(gu, 0, 0, ns->viscosity_type);
	    p = dat_p[i];

	    FLOAT eu[DDim];
	    /* eu = .5 * gu + .5 * gu^T */
	    MAT3_SYM(gu, eu);

	    /* eu = eu / 1000 */
	    MAT3_AXPBY(0., gu, 1./LEN_SCALING, eu);

	    /* eu * n */
	    MAT3_MV(eu, normal, en);
	
	    dat_f[i] = p * PRES_SCALING - nu * INNER_PRODUCT(en, normal);
	}
	    
    }

    phgExportVTK(g, "surf_force.vtk", sforce_dof, ns->p[1], ns->u[1], NULL);

    /* phgDofFree(&resu); */
    /* phgDofFree(&resp); */

    phgFinalize();
    exit(1);

    return;
}



void 
mark_inactive(NSSolver *ns)
{
    GRID *g = ns->g;
    LAYERED_MESH *gL = ns->gL;
    SIMPLEX **elist, *e;
    int i, ii, k, n_inactive = 0; 
    


    /* first clear mark */
    ForAllElements(g, e) {
	for (k = 0; k < NFace; k++)
	    e->bound_type[k] &= ~BC_INACTIVE;
    }


    BOTTOM_FACE *fb = gL->face_bot;
    for (ii = 0; ii < gL->nface_bot; ii++, fb++) {
	TRIA *t = gL->trias + fb->tria;
	INT vert_tri[3] = {t->verts[0],
			   t->verts[1],
			   t->verts[2]};
	FLOAT *x_tri[3] = {gL->verts[vert_tri[0]], 
			   gL->verts[vert_tri[1]], 
			   gL->verts[vert_tri[2]]};

	/* any of height over height_eps  */
	for (k = 0; k < 3; k++)
	    if (x_tri[k][3] > HEIGHT_EPS * 1.01)
		break;
	if (k < 3)
	    continue;

	n_inactive++;
	for (i = 0; i < fb->ne; i++) { /* tets */
	    SIMPLEX *e = fb->elems[i];
	    for (k = 0; k < NFace; k++)
		e->bound_type[k] |= BC_INACTIVE;
	}	
    }

    phgUpdateBoundaryTypes(g);


    phgPrintf("   Mark inactive Elements %d.\n", n_inactive);

    return;
}








/* --------------------------------------------------------------------------------
 *
 *  Solution errors 
 *
 *
 * -------------------------------------------------------------------------------- */


static int
comp_coord(const void *p0, const void *p1)
/* compares two 'INT's (used in qsort and bsearch) */
{
    int i;
    FLOAT *Xa = (FLOAT *) p0;
    FLOAT *Xb = (FLOAT *) p1;
    const FLOAT eps = 1e-12;

    /* X */
    if (Xa[0] > Xb[0] + eps) 
	return 1;
    else if (Xa[0] < Xb[0] - eps)
	return -1;

    /* Y */
    if (Xa[1] > Xb[1] + eps) 
	return 1;
    else if (Xa[1] < Xb[1] - eps)
	return -1;

    /* /\* Z *\/ */
    /* if (Xa[2] > Xb[2] + eps)  */
    /* 	return 1; */
    /* else if (Xa[2] < Xb[2] - eps) */
    /* 	return -1; */

    /* Sigma */
    if (Xa[3] > Xb[3] + eps)
    	return 1;
    else if (Xa[3] < Xb[3] - eps)
    	return -1;
    

    return 0;			/* same */
}



void
save_solution_dat(NSSolver *ns)
/*
 * Save solution on vertex,
 * Reduce # of procs to ONE.
 *
 * */
{

    assert(!ns_params->solve_temp);

#if 0
    /* ------------------------------
     *     Use parallel interp
     * ------------------------------ */
    
    /* GIRD *g2; */
    /* g2 = phgNewGrid(-1); */
    /* phgSetPeriodicity(g2, ns_params->periodicity); */
    /* phgImportSetBdryMapFunc(my_bc_map); */
    /* if (!phgImport(g2, ns_params->fn, FALSE)) */
    /* 	phgError(1, "can't read file \"%s\".\n", ns_params->fn); */
    /* OCT_TREE *og = phgOctTreeBuild(ns->g); */
    /* phgFreeGrid(&g2); */
#endif



#if 0    
    /* ------------------------------
     *     Save all data
     * ------------------------------ */
    
    /* GRID *g = ns->g; */
    /* char data_file[1000], fname[1000]; */
    /* DOF *dof_coord, *dof_Gu, *dof_nu; */
    
    /* sprintf(data_file, OUTPUT_DIR "/result.dat"); */


    /* sprintf(fname, "%s.Crd", data_file);		 */
    /* save_dof_data3(g, dof_coord, fname); */

    /* sprintf(fname, "%s.u", data_file);		 */
    /* save_dof_data3(g, ns->u[1], fname); */

    /* sprintf(fname, "%s.p", data_file);		 */
    /* save_dof_data3(g, ns->p[1], fname); */

    /* sprintf(fname, "%s.Gu", data_file);		 */
    /* save_dof_data3(g, dof_Gu, fname); */

    /* sprintf(fname, "%s.nu", data_file);		 */
    /* save_dof_data3(g, dof_nu, fname); */
#endif    




    /* ------------------------------
     *     Save only vertex data
     * ------------------------------ */
#if 0    
    FILE *fp = fopen("result.dat", "w");
    int rank, i, k;

    for (rank = 0; rank < g->nprocs; rank++) {
	if (g->rank != rank)
	    continue;

	/* for (i = 0; i < g->nvert; i++) { */
	/*     if (!(g->types_vert & OWNER)) */
	/* 	continue; */
	/* } */

	/* 3+3+1+9+1 = 17 */

	/* Coord */
	fwrite(&g->verts[0][0], g->nvert * Dim, sizeof(FLOAT), fp);

	/* u */
	fwrite(ns->u[1]->data, g->nvert * Dim, sizeof(FLOAT), fp);

	/* p */
	fwrite(ns->p[1]->data, g->nvert, sizeof(FLOAT), fp);

	/* Gradu */
	fwrite(dof_Gu, g->nvert * Dim*Dim, sizeof(FLOAT), fp);

	/* nu */
	fwrite(dof_nv, g->nvert, sizeof(FLOAT), fp);
	
    }
    fclose(fp);
#endif


#define NDATA 19  // 3 + 1 + 3 + 1 + 9 + 1 + 1
    
    /* Gather to root */
    GRID *g = ns->g;
    INT ntotal, nlocal, *counts, *displs, ntotal0;
    FLOAT *send_buf, *recv_buf = NULL, *send_ptr;
    DOF *dof_nu = phgDofNew(g, DOF_P1, 1, "nu", DofNoAction);
    int i, rank;
    FILE *fp;

    phgPrintf("   *Save solution data.\n");
    counts = phgAlloc(2 * g->nprocs * sizeof(*counts));
    displs = counts + g->nprocs;

    nlocal = g->nvert * NDATA;
    MPI_Allgather(&nlocal, 1, MPI_INT,
		  counts, 1, MPI_INT, g->comm);

    ntotal = 0;
    for (rank = 0; rank < g->nprocs; rank++) {
	displs[rank] = ntotal;
	ntotal += counts[rank];
    }

    send_buf = phgCalloc(nlocal, NDATA * sizeof(FLOAT));
    if (phgRank == 0)
	recv_buf = phgCalloc(ntotal, NDATA * sizeof(FLOAT));

    phgInfo(0, "nlocal: %d, ntotal: %d\n", nlocal / NDATA, ntotal / NDATA);
    
    proj_gradu(ns, ns->gradu[1], NULL, 0);
    
    FLOAT *dat_Crd = &g->verts[0][0];
    FLOAT *dat_sigma = ns->coord_sigma->data;
    FLOAT *dat_u   = ns->u[1]->data;
    FLOAT *dat_p   = ns->p[1]->data;;
    FLOAT *dat_p_static  = ns->p_static->data;;
    FLOAT *dat_Gu  = ns->Gradu->data;
    FLOAT *dat_nu  = dof_nu->data;
    FLOAT *dat_dht = ns->dHt->data;

    //phgDofDump(ns->dHt);

    send_ptr = send_buf;
    for (i = 0; i < g->nvert; i++) {
	FLOAT nu;
	
	*(send_ptr++) = *(dat_Crd++);
	*(send_ptr++) = *(dat_Crd++);
	*(send_ptr++) = *(dat_Crd++);

	*(send_ptr++) = *(dat_sigma++);
	
	*(send_ptr++) = *(dat_u++);
	*(send_ptr++) = *(dat_u++);
	*(send_ptr++) = *(dat_u++);

    	*(send_ptr++) = *(dat_p++) - *(dat_p_static++);

	nu = get_effective_viscosity(dat_Gu, 0, 0, ns->viscosity_type);
	*(dat_nu++) = nu;
	
	*(send_ptr++) = *(dat_Gu++);
	*(send_ptr++) = *(dat_Gu++);
	*(send_ptr++) = *(dat_Gu++);
	*(send_ptr++) = *(dat_Gu++);
	*(send_ptr++) = *(dat_Gu++);
	*(send_ptr++) = *(dat_Gu++);
	*(send_ptr++) = *(dat_Gu++);
	*(send_ptr++) = *(dat_Gu++);
	*(send_ptr++) = *(dat_Gu++);

	*(send_ptr++) = nu;

	*(send_ptr++) = *(dat_dht++);
    }

    /* phgExportVTK(g, "output/all.vtk", */
    /* 		 ns->u[1], ns->p[1], ns->Gradu, dof_nu, ns->dHt, NULL); */
    
    
    MPI_Gatherv(send_buf, nlocal, PHG_MPI_FLOAT,
		recv_buf, counts, displs, 
		PHG_MPI_FLOAT, 0, g->comm);


    phgInfo(0, "   Gather soltuion done.\n");

    
    if (phgRank == 0) {
	/* for (i = 0; i < ntotal / NDATA; i++) { */
	/*     const FLOAT *v = recv_buf + i*NDATA; */
	/*     phgInfo(0, "Unsorted %20.10e %20.10e %20.10e    %20.10e %20.10e %20.10e %20.10e\n", */
	/* 	    v[0], v[1], v[2], v[3], v[4], v[5], v[6]); */
	/* } */

	
	/* Sort & Uniq */
        qsort(recv_buf, ntotal / NDATA, NDATA * sizeof(FLOAT), comp_coord);
	
	{
            INT i0 = 0;
            for (i = i0+1; i < ntotal / NDATA; i++) {
                int cmp = comp_coord(recv_buf + i0 * NDATA,
				     recv_buf + i  * NDATA);
                if (cmp < 0) {
                    i0++;
                    memmove(recv_buf + i0 * NDATA,
			    recv_buf + i  * NDATA,
			    NDATA * sizeof(FLOAT));
                }
            }
	    ntotal0 = i0 + 1;
        }
	phgInfo(0, "   ntotal: %d\n", ntotal0);

	/* Check */
	for (i = 0; i < ntotal0; i++) {
	    const FLOAT *v = recv_buf + i*NDATA;
	    phgInfo(3, "sorted %20.10e %20.10e %20.10e %20.10e     %20.10e %20.10e %20.10e ... %20.10e \n",
		    v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[NDATA-1]);
	}
	

	/* Save to file */
	fp = fopen("result.dat", "w");
	fwrite(&ntotal0, 1, sizeof(INT), fp);
	fwrite(recv_buf, ntotal0, NDATA * sizeof(FLOAT), fp);
	fclose(fp);
	phgFree(recv_buf);
    }
    
    phgFree(send_buf);

    phgDofFree(&dof_nu);
    return;
}

static BOOLEAN ref_loaded = FALSE;
static int n_ref_dat = 0;
static FLOAT *ref_dat = NULL;

void
load_ref_solution()
{
    FILE *fp;


    if (ref_loaded == FALSE) {
	fp = fopen(ns_params->ref_sol_file, "r");
	assert(fp != NULL);
	assert(n_ref_dat == 0);
	fread(&n_ref_dat, 1, sizeof(INT), fp);
	ref_dat = phgCalloc(n_ref_dat, NDATA * sizeof(FLOAT));
	fread(ref_dat, n_ref_dat, NDATA * sizeof(FLOAT), fp);
	fclose(fp);
	ref_loaded = TRUE;
    }

}


const FLOAT *
interp_ref_solution(FLOAT *X, FLOAT sigma)
{
    FLOAT Xext[NDATA];
    FLOAT *found;
    
    Xext[0] = X[0];
    Xext[1] = X[1];
    Xext[2] = X[2];           /* z, not used  */
    Xext[3] = sigma;	      /* Use sigma */

    found = bsearch(Xext, ref_dat,
		    n_ref_dat, NDATA * sizeof(FLOAT), comp_coord);
    assert(found);

    return found + (Dim+1);
}



void
compute_error_ref(NSSolver *ns)
/* Compute error w.r.t ref  */
{
    int ndof = 15;		/* 3+1+9+1+1 */
    GRID *g = ns->g;
    DOF *dofs[ndof];
    DOF *errs[ndof];
    FLOAT *ddats[ndof];
    FLOAT *edats[ndof];
    int i, k;
    
    load_ref_solution();

    for (k = 0; k < ndof; k++) {
	dofs[k] = phgDofNew(g, DOF_P1, 1, NULL, DofNoAction);
	errs[k] = phgDofNew(g, DOF_P1, 1, NULL, DofNoAction);
	ddats[k] = dofs[k]->data;
	edats[k] = errs[k]->data;
    }

    proj_gradu(ns, ns->gradu[1], NULL, 0);

    DOF *uP1 = phgDofNew(g, DOF_P1, Dim, NULL, DofNoAction);
    DOF *pP1 = phgDofNew(g, DOF_P1, 1, NULL, DofNoAction);
    phgDofCopy(ns->u[1], &uP1, NULL, "uP1");
    phgDofCopy(ns->p[1], &pP1, NULL, "pP1");

    FLOAT *dat_Crd = &g->verts[0][0];
    FLOAT *dat_u   = uP1->data;
    FLOAT *dat_p   = pP1->data;;
    FLOAT *dat_p_static   = ns->p_static->data;;
    FLOAT *dat_Gu  = ns->Gradu->data;
    FLOAT *dat_dht = ns->dHt->data;
    FLOAT *coord_sigma = ns->coord_sigma->data;

    for (i = 0; i < g->nvert; i++) {
	const FLOAT *vref;
	FLOAT nu;

	vref = interp_ref_solution(g->verts[i], coord_sigma[i]);
	phgInfo(3, "found %20.10e %20.10e %20.10e %20.10e    %20.10e   %20.10e %20.10e %20.10e\n",
		vref[0], vref[1], vref[2], vref[3], vref[4], vref[5], vref[6], vref[7]);

	k = 0;

#define DIFF_DAT(solu_dat) {				\
	    *(edats[k] + i) = *(solu_dat++) - vref[k];	\
	    *(ddats[k] + i) = vref[k];			\
	    k++;					\
	}

	DIFF_DAT(dat_u);
	DIFF_DAT(dat_u);
	DIFF_DAT(dat_u);

	//DIFF_DAT(dat_p);
	{
	    *(edats[k] + i) = *(dat_p++) - *(dat_p_static++) - vref[k];	
	    *(ddats[k] + i) = vref[k];			
	    k++;					
	}

	nu = get_effective_viscosity(dat_Gu, 0, 0, ns->viscosity_type);
	
	DIFF_DAT(dat_Gu);
	DIFF_DAT(dat_Gu);
	DIFF_DAT(dat_Gu);
	DIFF_DAT(dat_Gu);
	DIFF_DAT(dat_Gu);
	DIFF_DAT(dat_Gu);
	DIFF_DAT(dat_Gu);
	DIFF_DAT(dat_Gu);
	DIFF_DAT(dat_Gu);

	{
	    *(edats[k] + i) = nu - vref[k];	
	    *(ddats[k] + i) = vref[k];			
	    k++;
	}

	DIFF_DAT(dat_dht);
	    
	assert(k == ndof);
    }
	

#define PRINT_ERR(name) {						\
	FLOAT uL2 = phgDofNormL2(dofs[k]);				\
	FLOAT eL2 = phgDofNormL2(errs[k]);				\
	FLOAT eMax = phgDofMaxValVec(errs[k]);				\
	FLOAT eMin = phgDofMinValVec(errs[k]);				\
	FLOAT uMax = phgDofMaxValVec(dofs[k]);				\
	FLOAT uMin = phgDofMinValVec(dofs[k]);				\
	phgPrintf("   Err " name					\
		  ": %20.10e %20.10e   %20.10e %20.10e   %20.10e %20.10e \n", \
		  eL2, eL2 / uL2, eMax, eMin, uMax, uMin		\
		  );							\
	k++;								\
    }
    
    k = 0;
    PRINT_ERR("u_");
    PRINT_ERR("v_");
    PRINT_ERR("w_");

    PRINT_ERR("p_");

    PRINT_ERR("ux");
    PRINT_ERR("uy");
    PRINT_ERR("uz");
    PRINT_ERR("vx");
    PRINT_ERR("vy");
    PRINT_ERR("vz");
    PRINT_ERR("wx");
    PRINT_ERR("wy");
    PRINT_ERR("wz");

    PRINT_ERR("nu");

    PRINT_ERR("dht");
    
    assert(k == ndof);
    

    for (k = 0; k < ndof; k++) {
	phgDofFree(&dofs[k]);
	phgDofFree(&errs[k]);
    }

    phgDofFree(&uP1);
    phgDofFree(&pP1);

    //dump_dof_ref(ns);
    return;
}


void
init_by_ref_solution(NSSolver *ns)
{
    GRID *g = ns->g;
    int i, k;
    
    phgPrintf("   Init by ref...\n");
    
    load_ref_solution();

    DOF *uP1 = phgDofNew(g, DOF_P1, Dim, NULL, DofNoAction);
    DOF *pP1 = phgDofNew(g, DOF_P1, 1, NULL, DofNoAction);

    FLOAT *dat_Crd = &g->verts[0][0];
    FLOAT *dat_u   = uP1->data;
    FLOAT *dat_p   = pP1->data;;
    FLOAT *coord_sigma = ns->coord_sigma->data;
    FLOAT *dat_nu   = ns->nu->data;;

    phgPrintf("Init by interp ref soltuion.\n");
    
    for (i = 0; i < g->nvert; i++) {
	const FLOAT *vref;
	FLOAT nu;

	vref = interp_ref_solution(g->verts[i], coord_sigma[i]);
	phgInfo(3, "found %20.10e %20.10e %20.10e %20.10e    %20.10e   %20.10e %20.10e %20.10e\n",
		vref[0], vref[1], vref[2], vref[3], vref[4], vref[5], vref[6], vref[7]);

	dat_u[i * Dim + 0] = vref[0];
	dat_u[i * Dim + 1] = vref[1];
	dat_u[i * Dim + 2] = vref[2];
	dat_p[i] = vref[3];
	dat_nu[i] = vref[13];
    }

    get_p_static(ns);
    phgDofAXPBY(1., ns->p_static, 1., &ns->p[1]); /* restore */

    /* #warning --- Debug set pressure static ---     */
    /*     phgDofAXPBY(1., ns->p_static, 0., &ns->p[1]); /\* restore *\/ */
    //phgExportVTK(g, "dof_nu.vtk", ns->u[1], ns->p[1], ns->nu, NULL);


    phgDofCopy(uP1, &ns->u[1], NULL, "u_{n+1}");
    phgDofCopy(pP1, &ns->p[1], NULL, "p_{n+1}");

    phgDofFree(&uP1);
    phgDofFree(&pP1);

    phgPrintf("   Init by ref Done.\n");

    return;
}










/* void */
/* dump_dof_ref(NSSolver *ns) */
/* /\* Compute error w.r.t ref  *\/ */
/* { */
/*     int ndof = 15;		/\* 3+1+9+1 *\/ */
/*     GRID *g = ns->g; */
/*     DOF *dofs[ndof]; */
/*     DOF *errs[ndof]; */
/*     FLOAT *ddats[ndof]; */
/*     FLOAT *edats[ndof]; */
/*     int i, k; */
    
/*     for (k = 0; k < ndof; k++) { */
/* 	dofs[k] = phgDofNew(g, DOF_P1, 1, NULL, DofNoAction); */
/* 	errs[k] = phgDofNew(g, DOF_P1, 1, NULL, DofNoAction); */
/* 	ddats[k] = dofs[k]->data; */
/* 	edats[k] = errs[k]->data; */
/*     } */

/*     FLOAT *dat_Crd = &g->verts[0][0]; */
/*     FLOAT *dat_u   = ns->u[1]->data; */
/*     FLOAT *dat_p   = ns->p[1]->data;; */
/*     FLOAT *dat_Gu  = ns->Gradu->data; */
/*     FLOAT *dat_dht = ns->dHt->data; */

/*     for (i = 0; i < g->nvert; i++) { */
/* 	const FLOAT *vref; */
/* 	FLOAT nu; */

/* 	vref = interp_ref_solution(g->verts[i]); */


/* 	k = 0; */

/* #define DIFF_DAT(solu_dat) {				\ */
/* 	    *(edats[k] + i) = *(solu_dat++);		\ */
/* 	    *(ddats[k] + i) = vref[k];			\ */
/* 	    k++;					\ */
/* 	} */

/* 	DIFF_DAT(dat_u); */
/* 	DIFF_DAT(dat_u); */
/* 	DIFF_DAT(dat_u); */

/* 	DIFF_DAT(dat_p); */

/* 	nu = get_effective_viscosity(dat_Gu, 0, 0, ns->viscosity_type); */
	
/* 	DIFF_DAT(dat_Gu); */
/* 	DIFF_DAT(dat_Gu); */
/* 	DIFF_DAT(dat_Gu); */
/* 	DIFF_DAT(dat_Gu); */
/* 	DIFF_DAT(dat_Gu); */
/* 	DIFF_DAT(dat_Gu); */
/* 	DIFF_DAT(dat_Gu); */
/* 	DIFF_DAT(dat_Gu); */
/* 	DIFF_DAT(dat_Gu); */

/* 	{ */
/* 	    *(edats[k] + i) = nu; */
/* 	    *(ddats[k] + i) = vref[k];			 */
/* 	    k++; */
/* 	} */
	
/* 	assert(k == ndof); */
/*     } */
	

/* #define PRINT_VEC(name) {					\ */
/* 	MAP *map = phgMapCreate(dofs[k], NULL);			\ */
/* 	VEC *vec = phgMapCreateVec(map, 1);			\ */
/* 								\ */
/* 	phgMapDofToLocalData(map, 1, &dofs[k], vec->data);	\ */
/* 	phgVecDumpMATLAB(vec, name, name"r_.m");		\ */
/* 								\ */
/* 	phgMapDofToLocalData(map, 1, &errs[k], vec->data);	\ */
/* 	phgVecDumpMATLAB(vec, name, name"_.m");			\ */
/* 								\ */
/* 	phgMapDestroy(&map);					\ */
/* 	k++;							\ */
/*     } */

    
/*     k = 0; */
/*     PRINT_VEC("u"); */
/*     PRINT_VEC("v"); */
/*     PRINT_VEC("w"); */

/*     PRINT_VEC("p"); */


/*     for (k = 0; k < ndof; k++) { */
/* 	phgDofFree(&dofs[k]); */
/* 	phgDofFree(&errs[k]); */
/*     } */

/*     return; */
/* } */


/* --------------------------------------------------------------------------------
 *
 *  Misc utils
 *
 *
 * -------------------------------------------------------------------------------- */

const char *
bdrymask2bit(short k)
{
    static char string[1024];

#if 0
    for (int c = 16; c >= 0; c--) {
	k = n >> c;

	if (k & 1)
	    printf("1");
	else
	    printf("0");
    }
#else

    int i = 0;

    if (k == 0) {
	string[i++] = 'U';
	string[i++] = 'r';
	string[i] = '\0';
	return string;
    }

    if (k & OWNER)
	string[i++] = 'O';

    if (k & INTERIOR)
	string[i++] = 'I';

    if (k & REMOTE)
	string[i++] = 'R';

    if (k & DIRICHLET)
	string[i++] = 'D';

    if (k & NEUMANN)
	string[i++] = 'N';

    if (k & BDRY_USER0) {
	string[i++] = 'u';
	string[i++] = '0';
    }
    if (k & BDRY_USER1) {
	string[i++] = 'u';
	string[i++] = '1';
    }
    if (k & BDRY_USER2) {
	string[i++] = 'u';
	string[i++] = '2';
    }
    if (k & BDRY_USER3) {
	string[i++] = 'u';
	string[i++] = '3';
    }
    if (k & BDRY_USER4) {
	string[i++] = 'u';
	string[i++] = '4';
    }
    if (k & BDRY_USER5) {
	string[i++] = 'u';
	string[i++] = '5';
    }
    if (k & BDRY_USER6) {
	string[i++] = 'u';
	string[i++] = '6';
    }
    if (k & BDRY_USER7) {
	string[i++] = 'u';
	string[i++] = '7';
    }
    if (k & BDRY_USER8) {
	string[i++] = 'u';
	string[i++] = '8';
    }
    if (k & BDRY_USER9) {
	string[i++] = 'u';
	string[i++] = '9';
    }

    if (k & UNDEFINED) {
	string[i++] = 'U';
	string[i++] = 'd';
    }

    string[i] = '\0';
    return string;

#endif
    
}




int
phgCompEdge(const void *p0, const void *p1)
/* compares two 'INT's (used in qsort and bsearch) */
{
    int i;
    EdgeVerts *vts0 = (EdgeVerts *) p0;
    EdgeVerts *vts1 = (EdgeVerts *) p1;

    i = (*vts0)[0] - (*vts1)[0];

    if (i < 0) {
	return -1;
    } else if (i > 0) {
	return 1;
    } else {
	i = (*vts0)[1] - (*vts1)[1];
	return (i < 0) ? -1 : ((i == 0) ? 0 : 1);
    }

}


int
phgCompEdgeMark(const void *p0, const void *p1)
/* compares edges with mark (used in qsort and bsearch) */
{
    int i;
    EdgeMarks *vts0 = (EdgeMarks *) p0;
    EdgeMarks *vts1 = (EdgeMarks *) p1;

    i = (*vts0)[0] - (*vts1)[0];

    if (i < 0) {
	return -1;
    } else if (i > 0) {
	return 1;
    } else {
	i = (*vts0)[1] - (*vts1)[1];
	return (i < 0) ? -1 : ((i == 0) ? 0 : 1);
    }

}


void get_avg_n(GRID *g, DOF *sur_normal)
{
    SIMPLEX *e;
    int i, j, k, l, s, q;
    int i_face, j_face;

    MAP *map_n = phgMapCreate(sur_normal, NULL);
    VEC *vec_n = phgMapCreateVec(map_n, Dim);
    phgVecDisassemble(vec_n);

    ForAllElements(g, e)
    {
        int M = sur_normal->type->nbas;	
        FLOAT rhs[M][Dim];
        BOOLEAN btype[M][Dim];
        INT local_map_idx[M][Dim];

        memset(rhs, 0, sizeof rhs); 
        memset(btype, 0, sizeof btype); 

        for (s = 0; s < NFace; s++)
        {
            int M_face = 3*(sur_normal->type->np_vert + sur_normal->type->np_edge);
            SHORT bases[M_face];
            const FLOAT *normal;

            if (e->bound_type[s] & BC_TOP || e->bound_type[s] & BC_BOTTOM)
            {
                phgDofGetBasesOnFace(sur_normal, e, s, bases);

                normal = phgGeomGetFaceOutNormal(g, e, s);

                for (i = 0; i < M_face; i++)
                {
                    i_face = bases[i];

                    for (k = 0; k < Dim; k++)
                    {
                        rhs[i_face][k] += normal[k]; 
                    }

                }

            }
        }

        for (i = 0; i < M; i++) 
        {
            for (k = 0; k < Dim; k++)
            {
            btype[i][k] = phgDofGetElementBoundaryType(sur_normal, e, i*Dim+k);
            if (!(btype[i][k] & BC_TOP) && !(btype[i][k] & BC_BOTTOM)) 
            {
                rhs[i][k] = 0.;
            } 
            }
        }


        for (i = 0; i < M; i++)
        for (k = 0; k < Dim; k++)
        {
            local_map_idx[i][k] = phgMapE2L(map_n, 0, e, i*Dim+k);
        }
        phgVecAddEntries(vec_n, 0, M * Dim, local_map_idx[0], &rhs[0][0]);
    }

    phgVecAssemble(vec_n);
    phgMapLocalDataToDof(map_n, 1, &sur_normal, vec_n->data);

    
    FLOAT *data = DofData(sur_normal);
    INT len = DofGetDataCount(sur_normal)/3.;
    FLOAT nx, ny, nz, val;

    for (i = 0; i < len; i++)
    {

        nx = data[Dim*i];
        ny = data[Dim*i+1];
        nz = data[Dim*i+2];

        val = sqrt(nx*nx+ny*ny+nz*nz);

        if (val != 0)
        {
        nx = nx/val;ny = ny/val; nz = nz/val;
        }
        else
        {
            nx=0;ny=0;nz=0;
        }

        data[Dim*i+0] = nx;
        data[Dim*i+1] = ny;
        data[Dim*i+2] = nz;
    }
   

    phgMapDestroy(&map_n);
    phgVecDestroy(&vec_n);


    //phgExportVTK(g, "avg_n.vtk", sur_normal, NULL);
}

