/*
 *  Subroutines of Building Mat & RHS of solver
 *
 *
 *  */
#include "ins.h"
#include "mat_op3.h"
#define _nsp (ns->ns_params)


#define USE_FS 1

#define SIGN_FRICTION -1
#define SIGN_STAB -1.

#define QUAD_EXTRA_ORDER 0





static void
proj_pres_coef(GRID *g, SIMPLEX *e, QUAD *quad, FLOAT *A)
/*
 * Find the coef of (p - \pi p, q - \pi q) expressed as (A \nabla p, \nabla q)
 * In the case to tetrahedron, this is possible because \nabla p is constant in element,
 *   p - \pi p = p_x (x - x_0) + p_y (y - y_0) + p_z (z - z_0),
 *   thus A_ij = \int (x_i - x_i_0) * (x_j - x_j_0) / \int
 *
 * 
 * */
{
    const FLOAT *w, *p;
    FLOAT X[Dim], X0[Dim];
    FLOAT *VX[NVert];
    FLOAT lambda0[] = {.25, .25, .25, .25};
    int q, i, j, k;

#define EVAL_COORD(X, lambda) {			\
	Bzero(X);				\
	for (i = 0; i < NVert; i++)		\
	    for (k = 0; k < Dim; k++)		\
		X[k] += VX[i][k] * lambda[i];	\
    }
    

    if (ns_params->stab_scheme == 0) {
	/* (p - \Pi p, q - \Pi q) */
    }
    else if (ns_params->stab_scheme == 1) {
	/*
	 * P_ij = \int x_i x_j
	 *
	 *  */
	for (i = 0; i < NVert; i++) {
	    VX[i] = g->verts[e->verts[i]];
	    phgInfo(3, "Coord %30.15e %30.15e %30.15e \n",
		    VX[i][0], VX[i][1], VX[i][2]);
	}

	EVAL_COORD(X0, lambda0);

	bzero(A, Dim*Dim * sizeof(*A));
    
	p = quad->points;
	w = quad->weights;
	for (q = 0; q < quad->npoints; q++) {
	    EVAL_COORD(X, p);

	    for (i = 0; i < Dim; i++)
		for (j = 0; j < Dim; j++)
		    A[i*Dim + j] +=
			(*w) * (X[i] - X0[i]) * (X[j]- X0[j]);
	
	    w++; p += Dim+1;
	}
    }
    else if (ns_params->stab_scheme == 2) {
	bzero(A, Dim*Dim * sizeof(*A));

	FLOAT dx = ns_params->dhx;
	FLOAT dz = ns_params->dhz;
	FLOAT dx2 = dx*dx, dz2 = dz*dz;

	A[0] = dx2;    A[1] = 0;         A[2] = 0;
	A[3] = 0;      A[4] = dx2;       A[5] = 0;
	A[6] = 0;      A[7] = 0;         A[8] = dz2;
    }

    return;
}



#define MAT3_SYM_FO(a, b) {			\
	FLOAT ux = a[0];			\
	FLOAT uy = a[1];			\
	FLOAT uz = a[2];			\
	FLOAT vx = a[3];			\
	FLOAT vy = a[4];			\
	FLOAT vz = a[5];			\
	FLOAT wx = a[6];			\
	FLOAT wy = a[7];			\
	FLOAT wz = a[8];			\
						\
	b[0] = a[0];				\
	b[4] = a[4];				\
	b[8] = a[8];				\
						\
	b[1] = b[3] = .5 * (a[1] + a[3]);	\
	b[2] = .5 * (a[2]);			\
	b[5] = .5 * (a[5]);			\
	b[6] = b[7] = 0.;			\
    }



/* static void */
/* dof_eval(DOF *dof, const FLOAT *elem_dat, , FLOAT *dof_value) */
/* { */
/*     int i, j; */
/*     int dim = dof->dim; */
/*     FLOAT *p; */

/*     for (j = 0; j < dim; j++) */
/* 	dof_value[j] = 0.; */
    
/*     for (i = 0; i < NBAS_PRISM; i++) { */
/* 	for (j = 0; j < dim; j++) */
/* 	    dof_value[j] += bas_value[i] * elem_dat[i * dim + j];  */
/*     } */
/* } */


static void
grad_dof_eval(DOF *dof, const FLOAT *elem_dat, SIMPLEX *e, QUAD *quad, int q, FLOAT *grad_value)
{
    int i, j, k;
    int N = dof->type->nbas;
    int dim = dof->dim;
    FLOAT *p;

    for (j = 0; j < dim; j++)
	for (k = 0; k < Dim; k++)
	    grad_value[j*Dim + k] = 0.;
    
    for (i = 0; i < N; i++) {
	const FLOAT *ggi =
	    phgQuadGetBasisGradientIsop(e, dof, i, quad, q);   
	
	for (j = 0; j < dim; j++)
	    for (k = 0; k < Dim; k++)
		grad_value[j * Dim + k]
		    += ggi[k] * elem_dat[i * dim + j]; 
    }
}



static FLOAT *
get_bas_dot_normal(const FLOAT *bas, const FLOAT *normal, INT i_e, INT j_e)
{
	static FLOAT basdotn[Dim][Dim];

	basdotn[0][0] = normal[0]*normal[0]*bas[i_e]*bas[j_e];
	basdotn[0][1] = normal[1]*normal[0]*bas[i_e]*bas[j_e];
	basdotn[0][2] = normal[2]*normal[0]*bas[i_e]*bas[j_e];
	basdotn[1][0] = normal[0]*normal[1]*bas[i_e]*bas[j_e];
	basdotn[1][1] = normal[1]*normal[1]*bas[i_e]*bas[j_e];
	basdotn[1][2] = normal[2]*normal[1]*bas[i_e]*bas[j_e];
	basdotn[2][0] = normal[0]*normal[2]*bas[i_e]*bas[j_e];
	basdotn[2][1] = normal[1]*normal[2]*bas[i_e]*bas[j_e];
	basdotn[2][2] = normal[2]*normal[2]*bas[i_e]*bas[j_e];

	return basdotn[0];
}


/****************************************************************
 * Build RHS which is the residual of the nonlinear system.
 ***************************************************************/
void 
phgNSBuildSolverURHS(NSSolver *ns)
{
    GRID *g = ns->g;
    SIMPLEX *e;
    SOLVER *solver_u = ns->solver_u;
    int i, k, l, q, s;
    FLOAT *dt = ns->dt;
    BOOLEAN tstep_minus = (ns->u[-1] != NULL);
    VEC *vec_rhs = phgMapCreateVec(solver_u->rhs->map, 1);
    FLOAT Theta = _nsp->Theta, nu = _nsp->nu, Thet1;
    int viscosity_type = ns->viscosity_type;
    FLOAT dhx = ns_params->dhx;
    FLOAT dhz = ns_params->dhz;
    FLOAT Vol = 0., Vol0;
    DOF *gradp, *gradp_static;


    SURF_BAS *surf_bas = ns->surf_bas;
    DOF *surf_dof = surf_bas->dof;
    BOOLEAN *rotated = surf_bas->rotated;
    const FLOAT *Trans = DofData(surf_dof);
    assert(surf_dof->type == ns->du->type);

 
#if STEADY_STATE
    assert(fabs(Theta - 1) < 1e-12);
    Thet1 = 0; Unused(Thet1);
    Unused(dt);
#else
    Thet1 = 1 - Theta;
    Unused(dt);
#endif /* STEADY_STATE */

    if (phgUseIsop)
	phgPrintf("  With Isop  ");
    else
	phgPrintf("  WithOut Isop  ");

    /* Removed */
    /* if (_nsp->use_prism_elem) { */
    /* 	phgNSBuildSolverURHSPrism(ns); */
    /* 	return; */
    /* } */

    if (ns_params->stab_scheme >= 1) 
	gradp = phgDofGradient(ns->p[1], NULL, NULL, NULL);


    if (ns_params->stab_scheme >= 0) {
	assert(ns_params->utype == DOF_P1);
    }

    if (ns_params->stab_remove_static > 1) {
	/* Note:
	 *    solve p' = p - p_static,
	 *    and p_static might be DOF_P2,
	 *    thus the follow is only on vert.
	 * */    
	phgDofAXPBY(-1, ns->p_static, 1., &ns->p[1]);
    }
    if (ns_params->stab_remove_static) {
	gradp_static = phgDofGradient(ns->p_static, NULL, NULL, NULL);
    }

    
    phgPrintf("   DB_mask: [");
    for (k = 0; k < Dim; k++)
	phgPrintf("%d ", ns->du->DB_masks[k]);
    phgPrintf("]   ");
    phgPrintf("  p mask: %d\n", ns->dp->DB_mask);

    nu_max = -1e100;
    nu_min = +1e100;

	//DOF *avg_n = phgDofNew(g, DOF_P2, 3, "avg n", DofNoAction);
    //get_avg_n(g, avg_n);

    phgVecDisassemble(vec_rhs);
    ForAllElements(g, e) {
	int M = ns->du->type->nbas;	/* num of bases of Velocity */
	int N = ns->dp->type->nbas;	/* num of bases of Pressure */
	int order = 2 * DofTypeOrder(ns->du, e) + QUAD_EXTRA_ORDER; /* Note:
								     *   quad order is really high here,
								     *   highest order term (u \nabla u, phi)  */
	FLOAT bufu[M], bufp[N], rhsu[M][Dim], rhsp[N], resUn[M];
	INT Iu[M][Dim], Ip[N], IUn[M];
	QUAD *quad;
	FLOAT vol, area, det;
	FLOAT dof_u_n1_values[M][Dim];
	FLOAT dof_u_n0_values[M][Dim];
	FLOAT dof_p_n1_values[N];
	FLOAT dof_p_st_values[N];
	FLOAT *normal;
	const FLOAT *w, *p, 
		**vu, *vu_queue[3],
		*vf[2], *gu[2], *vp[2], *gp, *vp_static, *gp_static,
		*vw, *vTe;
	FLOAT *vf_cache[2];
	FLOAT vp_mean = 0., vp_static_mean, P[Dim][Dim];

	vu = vu_queue + 1;

	quad = phgQuadGetQuad3D(order);
	vu[0] = phgQuadGetDofValues(e, ns->u[0], quad);	                /* u^{n} */
	vp[0] = phgQuadGetDofValues(e, ns->p[0], quad);	                /* p^{n} */
	gu[0] = phgQuadGetDofValues(e, ns->gradu[0], quad);             /* grad u^{n} */
	phgDofGetElementDatas(ns->u[0], e, dof_u_n0_values[0]);

	if (ns_params->stab_scheme >= 1) {
	    gp = phgQuadGetDofValues(e, gradp, quad);                   /* grad p */
	    phgDofGetElementDatas(ns->p[1], e, dof_p_n1_values);
	}
	if (ns_params->stab_remove_static) {
	    vp_static = phgQuadGetDofValues(e, ns->p_static, quad);     /* p_staic */
	    gp_static = phgQuadGetDofValues(e, gradp_static, quad);     /* grad p_static */
	    phgDofGetElementDatas(ns->p_static, e, dof_p_st_values);
	}

	
	//vw = phgQuadGetDofValues(e, ns->wind, quad);                  /* wind */
	if (ns_params->noniter_temp)			                /* nonlinear temp iter */
	    vTe = phgQuadGetDofValues(e, ns->T[1], quad);	        /* T^{n+1} */
	else
	    vTe = phgQuadGetDofValues(e, ns->T[0], quad);	        /* T^{n} */

	if (tstep_minus) { 
	    vu[-1] = phgQuadGetDofValues(e, ns->u[-1], quad);            /* u^{n-1} */
	} else {
	    vu[-1] = vu[0];
	}

#if STEADY_STATE || TIME_DEP_NON
	vu[1] = phgQuadGetDofValues(e, ns->u[1], quad);          /* u^{n+1} */
	gu[1] = phgQuadGetDofValues(e, ns->gradu[1], quad);      /* grad u^{n} */
	vp[1] = phgQuadGetDofValues(e, ns->p[1], quad);          /* p^{n+1} */
	phgDofGetElementDatas(ns->u[1], e, dof_u_n1_values[0]);
#else
	TIME_DEP_LINEAR_ENTRY; /* Unavailable */
#endif /* STEADY_STATE || TIME_DEP_NON */
	
	if (ns_params->stab_scheme >= 1) {
	    proj_pres_coef(g, e, quad, P[0]);
	}

	Unused(l);
	Bzero(vf); Bzero(vf_cache); 

	if (_nsp->extern_force) {
	    /* cache f values,
	     * only time_{n} */
	    for (l = 1; l < 2; l++) {
		const FLOAT *cache;
		size_t cache_size;
		setFuncTime(ns->time[l]); /* set static time in ins-test.c */

		/* cache f */
		cache_size = Dim * quad->npoints * sizeof(FLOAT);
		cache = phgQuadGetFuncValues(g, e, Dim, func_f, quad);
		vf[l] = vf_cache[l] = phgAlloc(cache_size);
		memcpy(vf_cache[l], cache, cache_size);

		phgQuadGetFuncValues(NULL, NULL, 0, NULL, NULL); /* clear cache */
	    }
	}

	/* Global Matrix */
	Bzero(rhsu); Bzero(rhsp);
	Bzero(resUn);
    
	/* average pressure */
	if (ns_params->stab_scheme == 0) {
	    vp_mean = 0.;
	    vp_static_mean = 0.;
	    w = quad->weights;
	    for (q = 0; q < quad->npoints; q++, w++) {
		vp_mean += (*w) * vp[1][q];
		if (ns_params->stab_remove_static == 1) 
		    vp_static_mean += (*w) * vp_static[q];
	    }
	}

	p = quad->points;
	w = quad->weights;
	for (q = 0; q < quad->npoints; q++) {
	    const FLOAT *J;
	    if (phgUseIsop) {
		J = phgGeomGetJacobian_(g, e, p, &det);
		vol = fabs(det) / 6;
	    }
	    else {
		vol = phgGeomGetVolume(g, e);
	    }
	    Vol += vol*(*w);

	    if (phgUseIsop) {
		static FLOAT gradu_n0_values[Dim*Dim];
		static FLOAT gradu_n1_values[Dim*Dim];

		grad_dof_eval(ns->du, dof_u_n0_values[0], e, quad, q, gradu_n0_values);
		grad_dof_eval(ns->du, dof_u_n1_values[0], e, quad, q, gradu_n1_values);
		
		gu[0] = gradu_n0_values;
		gu[1] = gradu_n1_values;
	    }

	    nu = get_effective_viscosity(gu[1], *vTe, 0, viscosity_type);


#if SIMPLE_TEST
	    /* ------------------------------
	     *
	     * Simple test for Stokes eqn convergence test
	     *
	     * ------------------------------ */
	    FLOAT hmax = phgGeomGetDiameter(g, e);
	    
	    /* rhs u */
	    for (i = 0; i < M; i++) {
		/* interior node or Neumann */
		const FLOAT *gi_u = phgQuadGetBasisValues(e, ns->du, i, quad) + q;       /* phi_i */
		const FLOAT *ggi_u = phgQuadGetBasisGradientIsop(e, ns->du, i, quad, q);    /* grad phi_i */

		for (k = 0; k < Dim; k++) {
		    
		    rhsu[i][k] += vol*(*w) * (- INNER_PRODUCT(gu[1] + k*Dim, ggi_u) 
					      + (*vp[1]) * *(ggi_u+k) 
					      );

		    rhsu[i][k] += vol*(*w) * ( vf[1][k] * (*gi_u) );
		}
	    }

	    /* rhs p */
	    for (i = 0; i < N; i++) {
		const FLOAT *gi_p = phgQuadGetBasisValues(e, ns->dp, i, quad) + q;       /* psi_i */
		const FLOAT *ggi_p = phgQuadGetBasisGradientIsop(e, ns->dp, i, quad, q);    /* grad phi_i */
		FLOAT divu1 = gu[1][0] + gu[1][4] + gu[1][8];
		//FLOAT divu0 = gu[0][0] + gu[0][4] + gu[0][8];
		rhsp[i] += vol*(*w) * (divu1 * (*gi_p)
				       );

		if (ns_params->stab_scheme >= 1) {
		    rhsp[i] -= SIGN_STAB * vol*(*w) * hmax * hmax * INNER_PRODUCT(gp, ggi_p);
		}
	    }	  /* end res p */
#else
	    /* ------------------------------
	     *
	     * Ice sheet dynamics
	     *
	     * ------------------------------ */
	    
	    /* rhs u */
	    for (i = 0; i < M; i++) {
		/* interior node or Neumann */
		const FLOAT *gi_u = phgQuadGetBasisValues(e, ns->du, i, quad) + q;       /* phi_i */
		const FLOAT *ggi_u = phgQuadGetBasisGradientIsop(e, ns->du, i, quad, q);    /* grad phi_i */

		for (k = 0; k < Dim; k++) {
		    FLOAT eu[DDim];

#   if USE_FS
		    MAT3_SYM(gu[1], eu);
#   else
#   warning ---------- Yet another FO -----------	
		    MAT3_SYM_FO(gu[1], eu);
#   endif

		    rhsu[i][k] += vol*(*w) * (- nu * INNER_PRODUCT(eu+k*Dim, ggi_u) 
					      + (*vp[1]) * *(ggi_u+k) * LEN_SCALING * PRES_SCALING
					      ) * EQU_SCALING;     /* left */

		    if (ns_params->stab_remove_static > 1) {
			rhsu[i][k] += vol*(*w) * (
			    + (*vp_static) * *(ggi_u+k) * LEN_SCALING * PRES_SCALING
						  ) * EQU_SCALING;     /* left */
		    }

		    if (k == Z_DIR) { 
			const FLOAT rho = RHO_ICE;
			const FLOAT grav = GRAVITY;
			const FLOAT a = SEC_PER_YEAR;
			const FLOAT f = rho*grav * LEN_SCALING2 * EQU_SCALING; 

			Unused(a);
			rhsu[i][k] += vol*(*w) * (-f * (*gi_u) 
						  ); /* right */
		    }

		    if (_nsp->compensate_equ)
			rhsu[i][k] += vol*(*w) * EQU_SCALING * (- vf[1][k] * (*gi_u) * LEN_SCALING2);
		}
	    }
	    /* 
	     * Notes on scaling:
	     *   The momentum equation is scaled with EQU_SCALING * LEN_SCALING * LEN_SCALING,
	     *       1. (nu eu, eu) term has two gradients thus has scaling EQU_SCALING
	     *       2. grad p term has one gradients thus has scaling EQU_SCALING * LEN_SCALING
	     *       3. f term has none gradients thus has scaling EQU_SCALING * LEN_SCALING * LEN_SCALING
	     *
	     * */
		    

	    

	    /* rhs p */
	    for (i = 0; i < N; i++) {
		const FLOAT *gi_p = phgQuadGetBasisValues(e, ns->dp, i, quad) + q;       /* psi_i */
		const FLOAT *ggi_p = phgQuadGetBasisGradientIsop(e, ns->dp, i, quad, q);    /* grad phi_i */
		FLOAT divu1 = gu[1][0] + gu[1][4] + gu[1][8];
		//FLOAT divu0 = gu[0][0] + gu[0][4] + gu[0][8];
		rhsp[i] += vol*(*w) * (divu1 * (*gi_p)
				       );

		const FLOAT alpha_stab = ns_params->stab_alpha;
		FLOAT nu_stab = 1.;
		if (ns_params->stab_nu < 0) {
		    nu_stab = nu;
		}
		else {
		    nu_stab = ns_params->stab_nu;
		}

	    
		if (ns_params->stab_scheme == 0) {
		    /* Proj */
		    if (ns_params->stab_remove_static == 1)
			rhsp[i] -= SIGN_STAB * vol*(*w) * 1./(EQU_SCALING * nu_stab) * alpha_stab
				* ( ((*vp[1] - *vp_static) - (vp_mean - vp_static_mean) ) * ((*gi_p) -  1./4.) ) ;
		    else
			rhsp[i] -= SIGN_STAB * vol*(*w) * 1./(EQU_SCALING * nu_stab) * alpha_stab
				* ( ((*vp[1]) - vp_mean) * ((*gi_p) -  1./4.) ) ;
		}
		else if (ns_params->stab_scheme >= 1) {

		    FLOAT Jgi[3], gradp[3];
		    //FLOAT f3[3] = {0, 0, -f};
		    const FLOAT *iJ_ = P[0];
		
		    if (ns_params->stab_remove_static == 1) {
			for (k = 0; k < Dim; k++) 
			    gradp[k] = gp[k] - gp_static[k];
		    }
		    else {
			for (k = 0; k < Dim; k++) 
			    gradp[k] = gp[k];
		    }
		    MAT3_MV(iJ_, ggi_p, Jgi);
		
		    rhsp[i] -= SIGN_STAB * vol*(*w) * 1./(EQU_SCALING * nu_stab) * alpha_stab 
			    * (INNER_PRODUCT(gradp, Jgi)
			       );
		} /* end stab scheme */
		
	    }	  /* end res p */
	    
#endif	/* SIMPLE_TEST */

   

	    
	    if (tstep_minus) 
		vu[-1] += Dim;

#if STEADY_STATE || TIME_DEP_NON
	    vu[1] += Dim;
	    gu[1] += Dim*Dim;
	    vp[1]++;
#else
	    TIME_DEP_LINEAR; /* Unavailable */
#endif /* STEADY_STATE || TIME_DEP_NON */
	    vu[0] += Dim;
	    gu[0] += Dim * Dim;
	    vp[0]++; 
	    //vw += Dim;
	    if (_nsp->extern_force) {
		vf[0] += Dim; vf[1] += Dim;
	    }
	    if (ns_params->stab_scheme >= 1) 
		gp += Dim;
	    if (ns_params->stab_remove_static) {
		vp_static++;
		gp_static += Dim;
	    }
	    vTe++;
	    w++; p += Dim + 1;
	}

	if (_nsp->extern_force) {
	    phgFree(vf_cache[0]);
	    phgFree(vf_cache[1]);
	}

	normal = NULL; Unused(normal);
	area = 0; Unused(area);


	/* slip and lateral  boundary */
	for (s = 0; s < NFace; s++) {
	    if (e->bound_type[s] & SLIP_BDRY) {
		/* ------------------------------
		 *
		 * 
		 * slip boundary
		 *
		 *
		 * ------------------------------ */
		int v0, v1, v2;
		int nbas_face = NbasFace(ns->du);
		SHORT bases[nbas_face];
		FLOAT lambda[Dim + 1], x,y,z, beta;
		order = 2 * DofTypeOrder(ns->du, e) + QUAD_EXTRA_ORDER;

		phgDofGetBasesOnFace(ns->du, e, s, bases);
		v0 = GetFaceVertex(s, 0);
		v1 = GetFaceVertex(s, 1);
		v2 = GetFaceVertex(s, 2);
		lambda[s] = 0.;
		    
		area = phgGeomGetFaceArea(g, e, s);
		normal = phgGeomGetFaceOutNormal(g, e, s);
		quad = phgQuadGetQuad2D(order);

		/* phgInfo(0, "\n\nBot face   : %20.10e [%20.10e %20.10e %20.10e]\n", */
		/* 	area, normal[0], normal[1], normal[2]); */
		
		p = quad->points;
		w = quad->weights;
		for (q = 0; q < quad->npoints; q++) {
		    FLOAT vu[Dim], un_multpiler;
		    lambda[v0] = *(p++);
		    lambda[v1] = *(p++);
		    lambda[v2] = *(p++);

		    if (phgUseIsop) {
			area = phgGeomGetFaceArea_(g, e, s, lambda, normal);
			/* phgInfo(0, "\nbot face  %d: %20.10e [%20.10e %20.10e %20.10e]\n", q, */
			/* 	area, normal[0], normal[1], normal[2]); */

			/* FLOAT X[Dim], normal_analy[Dim]; */
			/* phgGeomLambda2XYZ(g, e, lambda, X, X+1, X+2); */
			/* func_normal(X[0], X[1], X[2], normal_analy); */
			/* phgInfo(0, "coord     %d: %20.10e [%20.10e %20.10e %20.10e]\n", q, */
			/* 	area, X[0], X[1], X[2]); */
			/* phgInfo(0, "bot face  %d: %20.10e [%20.10e %20.10e %20.10e]\n", q, */
			/* 	area, normal_analy[0], normal_analy[1], normal_analy[2]); */
		    }
		    
		    phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);
		    phgDofEval(ns->u[1], e, lambda, vu);

		    if (ns_params->sliding_bdry_scheme == 1) /* Surf L2 proj */
			phgDofEval(ns->Un, e, lambda, &un_multpiler);


#if 0
		    func_beta(x, y, z, &beta);
#else
		    phgDofEval(ns->beta, e, lambda, &beta);
#endif

		    const FLOAT *gi_u = 
			    ns->du->type->BasFuncs(ns->du, e, 0, -1, lambda);

		    for (i = 0; i < nbas_face; i++) {
			int ii = bases[i];
			for (k = 0; k < Dim; k++) {
			    rhsu[ii][k] += SIGN_FRICTION * area*(*w) * beta * vu[k] * (gi_u[ii])
				    * EQU_SCALING * LEN_SCALING;
			    /* This should only be in tangential direction */
			}
			/* 
			 * Note on scaling: The area integral term is for area, thus compared to volume
			 *   intergral terms, it has scaling EQU_SCALING * LEN_SCALING,
			 *
			 * */

			if (ns_params->sliding_bdry_scheme == 1) {/* Surf L2 proj */
			    for (k = 0; k < Dim; k++) 
				resUn[ii] -= area*(*w) * normal[k] * vu[k] * (gi_u[ii]);

			    /* rhsu */
			    for (k = 0; k < Dim; k++) 
				rhsu[ii][k] -= area*(*w) * normal[k] * (un_multpiler) * (gi_u[ii]);
			}

		    }     /* end of bas_i */
		    
		    w++;
		} /* end of quad point */
	    }	  /* end of slip boundary (bottom)  */
	    else if (e->bound_type[s] & BC_FRONT ) {
		/* ------------------------------
		 *
		 *
		 * Calving front boundary
		 * 
		 *
		 * ------------------------------ */
		int v0, v1, v2;
		int nbas_face = NbasFace(ns->du);
		SHORT bases[nbas_face];
		FLOAT lambda[Dim + 1], x,y,z, beta;
		order = 2 * DofTypeOrder(ns->du, e) + QUAD_EXTRA_ORDER;

		phgDofGetBasesOnFace(ns->du, e, s, bases);
		v0 = GetFaceVertex(s, 0);
		v1 = GetFaceVertex(s, 1);
		v2 = GetFaceVertex(s, 2);
		lambda[s] = 0.;
		    
		area = phgGeomGetFaceArea(g, e, s);
		normal = phgGeomGetFaceOutNormal(g, e, s);
		quad = phgQuadGetQuad2D(order);

		/* phgInfo(0, "\n\nLateral face   : %20.10e [%20.10e %20.10e %20.10e]\n", */
		/* 	area, normal[0], normal[1], normal[2]); */
		
		p = quad->points;
		w = quad->weights;
		for (q = 0; q < quad->npoints; q++) {
		    FLOAT vp_water;
		    lambda[v0] = *(p++);
		    lambda[v1] = *(p++);
		    lambda[v2] = *(p++);

		    if (phgUseIsop) {
			area = phgGeomGetFaceArea_(g, e, s, lambda, normal);
			/* phgInfo(0, "\nbot face  %d: %20.10e [%20.10e %20.10e %20.10e]\n", q, */
			/* 	area, normal[0], normal[1], normal[2]); */

			/* FLOAT X[Dim], normal_analy[Dim]; */
			/* phgGeomLambda2XYZ(g, e, lambda, X, X+1, X+2); */
			/* func_normal(X[0], X[1], X[2], normal_analy); */
			/* phgInfo(0, "coord     %d: %20.10e [%20.10e %20.10e %20.10e]\n", q, */
			/* 	area, X[0], X[1], X[2]); */
			/* phgInfo(0, "bot face  %d: %20.10e [%20.10e %20.10e %20.10e]\n", q, */
			/* 	area, normal_analy[0], normal_analy[1], normal_analy[2]); */
		    }
		    
		    phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);
		    //phgDofEval(ns->p_static, e, lambda, &vp_static);

		    const FLOAT rho = RHO_WATER;
		    const FLOAT grav = GRAVITY;
		
		    if (z > 0)
			vp_water = 0.;
		    else
			vp_water = rho * grav * (-z) * LEN_SCALING;

		    const FLOAT *gi_u = 
			    ns->du->type->BasFuncs(ns->du, e, 0, -1, lambda);

		    for (i = 0; i < nbas_face; i++) {
			int ii = bases[i];
			for (k = 0; k < Dim; k++) {
			    rhsu[ii][k] += area*(*w) * (-1) * vp_water * normal[k] * (gi_u[ii])
				    * EQU_SCALING * LEN_SCALING;
			    /* 
			     * Note on scaling: The area integral term is for area, thus compared to volume
			     *   intergral terms, it has scaling EQU_SCALING * LEN_SCALING. Also note that the
			     *   coord z is in unit km, thus has scaling LEN_SCALING.
			     *
			     * */
			}

		    }     /* end of bas_i */
		    
		    w++;
		}		/* end of quad point */
	    }			/* end calving front */
            else if (e->bound_type[s] & BC_FLOAT)
            {
		/* ------------------------------
		 *
		 *
		 * Ice shelf boundary, by Tong
		 *
		 *
		 * ------------------------------ */
		
                FLOAT lambda[Dim+1], lambda_x, lambda_y, lambda_z;
                INT quad_order = 5;
                FLOAT Ns, vu[Dim];
                QUAD *quad = phgQuadGetQuad2D(quad_order);
                FLOAT *w = quad->weights;
                FLOAT *p = quad->points;

                FLOAT avg_n_v[Dim];
                FLOAT nx, ny, nz, x, y, z;

                FLOAT area = phgGeomGetFaceArea(g, e, s);
                FLOAT *normal = phgGeomGetFaceOutNormal(g, e, s);
                //FLOAT normal[Dim];

                INT v0 = GetFaceVertex(s, 0);
                INT v1 = GetFaceVertex(s, 1);
                INT v2 = GetFaceVertex(s, 2);
                lambda[s] = 0;

                FLOAT dt1 = ns->dt[0];
                INT M_face = NbasFace(ns->du);
                SHORT bas_idx_e[M_face];

                phgDofGetBasesOnFace(ns->du, e, s, bas_idx_e);

                for (q = 0; q < quad->npoints; q++) {
                    lambda[v0] = *(p++);
                    lambda[v1] = *(p++);
                    lambda[v2] = *(p++);
                    
		    phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);

		    phgDofEval(ns->u[1], e, lambda, vu);  /* OLD u value */

		    /* Use continuous normal */
		    assert(UN_DIR == X_DIR);
		    assert(surf_dof->type->continuity >= 0);
                    //phgDofEval(avg_n, e, lambda, avg_n_v);
		    //normal[0] = avg_n_v[0];
             //       normal[1] = avg_n_v[1];
             //       normal[2] = avg_n_v[2];
		    
		    if (fabs(normal[2]) < 1.0e-8) {
			Ns = 1.0e50;
			phgWarning("Base normal direction nearly horizontal!!!\n");
		    }
		    else {
			Ns = sqrt(1 + SQUARE(normal[0]/normal[2])+SQUARE(normal[1]/normal[2]));
		    }
                    
                    const FLOAT *bas = ns->du->type->BasFuncs(ns->du, e, 0, -1, lambda);
                    
		    if (z < 0)
			for (i = 0; i < M_face; i++)  {
			    //if (phgDofGetElementBoundaryType(ns->du, e, i) & BC_BOTTOM_GRD)
			    //    continue;

			    INT i_e = bas_idx_e[i];
			    const FLOAT *bas_i = phgQuadGetBasisValues(e, ns->du, i_e, quad);

			    for (k = 0; k < Dim; k++) {
				rhsu[i_e][k] += area * w[q] * RHO_WATER * GRAVITY * 
						( z * LEN_SCALING - INNER_PRODUCT(vu, normal) * Ns * dt1) * bas[i_e]
						* normal[k] * EQU_SCALING * LEN_SCALING;
			    /* 
			     * Note on scaling: see calving front
			     * */
			    }
			}  /* end bas face i */
                }	   /* end quad point */
            }		   /* end ice shelf face */


	    
	} /* end of all outflow face in element */


	    
	if (_nsp->compensate_equ) {
	    /* Compensate surface boundary condition */
	    for (s = 0; s < NFace; s++) {
		if (e->bound_type[s] & BC_TOP) { /* Note: only top surf is Neumann */
		    int v0, v1, v2;
		    int nbas_face = NbasFace(ns->u[1]);
		    SHORT bases[nbas_face];
		    FLOAT lambda[Dim + 1], x,y,z;
		    order = 2 * DofTypeOrder(ns->u[1], e) + QUAD_EXTRA_ORDER;

		    phgDofGetBasesOnFace(ns->u[1], e, s, bases);
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
			FLOAT gn[Dim];
			lambda[v0] = *(p++);
			lambda[v1] = *(p++);
			lambda[v2] = *(p++);
			
			phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);
			//if (e->bound_type[s] & BC_TOP)
			func_g1(x, y, z, gn); /* Surf */

			for (i = 0; i < nbas_face; i++) {
			    int ii = bases[i];
			    FLOAT gi_u = 
				    *ns->u[1]->type->BasFuncs(ns->u[1], e, ii, ii + 1, lambda);

			    for (k = 0; k < Dim; k++) {
				rhsu[ii][k] += area*(*w) * gn[k] * (gi_u) * LEN_SCALING
					* EQU_SCALING;
			    }
			}     /* end of bas_i */
			w++;
		    } /* end of quad point */
		}	  /* end of face outflow */
	    }	  /* end of all outflow face in element */
	}		  /* end of compensate equations */



	/* Rotate bases */
	if (ns_params->sliding_bdry_scheme == 0) /* Dirich */
	    for (i = 0; i < M; i++) {
		INT id = phgDofMapE2D(surf_dof, e, i * (Dim*Dim)) / (Dim*Dim);
		if (!rotated[id])
		    continue;	
		const FLOAT *trans = Trans + id*(Dim*Dim);

		trans_left_3d(&rhsu[i][0], 1, 1, trans);
	    }



	/* Map: Element -> system */
	for (i = 0; i < M; i++)
	    for (k = 0; k < Dim; k++)
		Iu[i][k] = phgMapE2L(solver_u->rhs->map, 0, e, i * Dim + k);
	for (i = 0; i < N; i++)
	    Ip[i] = phgMapE2L(solver_u->rhs->map, 1, e, i);

	if (ns_params->sliding_bdry_scheme == 1) /* surf L2 constraint */
	    for (i = 0; i < M; i++)
		IUn[i] = phgMapE2L(solver_u->rhs->map, 2, e, i);

	/* set velocity dirichlet bdry */
	FLOAT tmp[Dim];
	for (i = 0; i < M; i++)
	    for (k = 0; k < Dim; k++) {
		if (phgDofDirichletBC_(ns->du, e, i*Dim+k, NULL, bufu, tmp,
				       DOF_PROJ_NONE)) {
		    rhsu[i][k] = 0.;
		}
	    }


#if STEADY_STATE || TIME_DEP_NON
	/* set pressure dirichlet bdry for pinned point */
	for (i = 0; i < N; i++)
	    if (phgDofDirichletBC(ns->dp, e, i, NULL, bufp, &rhsp[i],
				  DOF_PROJ_NONE)) {
		if (!_nsp->enclosed_flow)
		    phgError(1, "No dirichlet bc for Unenclosed flow!\n");
		if (_nsp->pin_node) {
# if PIN_AT_ROOT 
		    if (g->rank != 0)
		    	phgError(1, "Pinned node only on rank 0!\n");
		    if (g, e->verts[i] != ns->pinned_node_id)
			phgError(1, "Build rhs: pinned node e:%d, bas:%d, [%d] and [%d] "
				 "doesn't coincide when build RHS!\n", 
				 e->index, i, e->verts[i], ns->pinned_node_id);
# else
		    if (GlobalVertex(g, e->verts[i]) != ns->pinned_node_id)
			phgError(1, "Build rhs: pinned node e:%d, bas:%d, [%d] and [%d] "
				 "doesn't coincide when build RHS!\n", 
				 e->index, i, e->verts[i], ns->pinned_node_id);
# endif /* PIN_AT_ROOT */
		}
	    } 
#else
	TIME_DEP_LINEAR; /* Unavailable */
#endif /* STEADY_STATE || TIME_DEP_NON */


	
	phgVecAddEntries(vec_rhs, 0, M * Dim, Iu[0], &rhsu[0][0]);
	phgVecAddEntries(vec_rhs, 0, N, Ip, rhsp);

	if (ns_params->sliding_bdry_scheme == 1) /* surf L2 constraint */
	    phgVecAddEntries(vec_rhs, 0, M, IUn, resUn);
    }				/* end element */
    
	//phgDofFree(&avg_n);

    phgVecAssemble(vec_rhs);
    phgVecAssemble(solver_u->rhs);
    phgVecAXPBY(1., vec_rhs, 0, &solver_u->rhs);
    solver_u->rhs_updated = FALSE;

    if (DUMP_MAT_VEC) {
	static int count = 0;
	char name[1000];
	phgPrintf("Dumping rhs\n");
	sprintf(name, "rhs%d_.m", count++);
	phgVecDumpMATLAB(ns->solver_u->rhs, "rhs", name);
    }

    phgVecDestroy(&vec_rhs);

    {
	FLOAT a[2] = {nu_max, -nu_min}, b[2]; 
	MPI_Allreduce(&a, &b, 2, PHG_MPI_FLOAT, MPI_MAX, g->comm);
	nu_max = b[0];
	nu_min = -b[1];
	phgPrintf("  vis: [%8.4e, %8.4e] ", nu_min, nu_max);
    }

    Vol0 = Vol;
    MPI_Allreduce(&Vol0, &Vol, 1, PHG_MPI_FLOAT, MPI_SUM, g->comm);
    //    phgPrintf("\nVol: %14.8e, ref: %14.8e\n", Vol * LEN_SCALING3, func_vol());
    phgPrintf("\nVol: %14.8e\n", Vol * LEN_SCALING3);


    if (ns_params->stab_scheme >= 1) 
	phgDofFree(&gradp);


    if (ns_params->stab_remove_static > 1) {
	/* recover p */
	phgDofAXPBY(1., ns->p_static, 1., &ns->p[1]); 
    } 
    if (ns_params->stab_remove_static) {
	phgDofFree(&gradp_static);
    }

    return;
}




#define eu_xx eu[0]
#define eu_xy eu[1]
#define eu_xz eu[2]
#define eu_yx eu[3]
#define eu_yy eu[4]
#define eu_yz eu[5]
#define eu_zx eu[6]
#define eu_zy eu[7]
#define eu_zz eu[8]

FLOAT * 
get_gbas_product_stokes(const FLOAT *gi, const FLOAT *gj,
		 const FLOAT *gu, LTYPE ltype) 
{
    static FLOAT prod[Dim][Dim];
    FLOAT Gi[Dim], Gj[Dim];
    FLOAT eu[DDim];
    FLOAT eps, eps2;
    int k;

    /* Picard term */
#if USE_FS    
    prod[0][0] = gi[0] * gj[0] + .5 * (gi[1] * gj[1] + gi[2] * gj[2]);
    prod[0][1] = .5 * gi[0] * gj[1];
    prod[0][2] = .5 * gi[0] * gj[2];

    prod[1][0] = .5 * gi[1] * gj[0]; 
    prod[1][1] = gi[1] * gj[1] + .5 * (gi[0] * gj[0] + gi[2] * gj[2]);
    prod[1][2] = .5 * gi[1] * gj[2];

    prod[2][0] = .5 * gi[2] * gj[0]; 
    prod[2][1] = .5 * gi[2] * gj[1]; 
    prod[2][2] = gi[2] * gj[2] + .5 * (gi[0] * gj[0] + gi[1] * gj[1]);
#else
#   warning ---------- Yet another FO -----------	
    prod[0][0] = gi[0] * gj[0] + .5 * (gi[1] * gj[1] + gi[2] * gj[2]);
    prod[0][1] = .5 * gi[0] * gj[1];
    prod[0][2] = 0; //.5 * gi[0] * gj[2];

    prod[1][0] = .5 * gi[1] * gj[0]; 
    prod[1][1] = gi[1] * gj[1] + .5 * (gi[0] * gj[0] + gi[2] * gj[2]);
    prod[1][2] = 0; //.5 * gi[1] * gj[2];

    prod[2][0] = 0; //.5 * gi[2] * gj[0]; 
    prod[2][1] = 0; //.5 * gi[2] * gj[1]; 
    prod[2][2] = gi[2] * gj[2]; //+ .5 * (gi[0] * gj[0] + gi[1] * gj[1]);
#endif

    if (ltype == PICARD) {
	return prod[0];
    }

    /* Newton term */
    MAT3_SYM(gu, eu);
    for (k = 0; k < DDim; k++)
	eu[k] /= LEN_SCALING;
    eps = sqrt(.5) * MAT3_NORM2(eu);
    
    if (eps < MIN_EFFECTIVE_STRAIN) 
	eps = MIN_EFFECTIVE_STRAIN;

    eps2 = - (1./3.) / (eps*eps);

    Gi[0] = eu_xx * gi[0] + eu_xy * gi[1] + eu_xz * gi[2];
    Gi[1] = eu_yx * gi[0] + eu_yy * gi[1] + eu_yz * gi[2];
    Gi[2] = eu_zx * gi[0] + eu_zy * gi[1] + eu_zz * gi[2];

    Gj[0] = eu_xx * gj[0] + eu_xy * gj[1] + eu_xz * gj[2];
    Gj[1] = eu_yx * gj[0] + eu_yy * gj[1] + eu_yz * gj[2];
    Gj[2] = eu_zx * gj[0] + eu_zy * gj[1] + eu_zz * gj[2];
    
    prod[0][0] += Gi[0] * Gj[0] * eps2;
    prod[0][1] += Gi[1] * Gj[0] * eps2;
    prod[0][2] += Gi[2] * Gj[0] * eps2;

    prod[1][0] += Gi[0] * Gj[1] * eps2;
    prod[1][1] += Gi[1] * Gj[1] * eps2;
    prod[1][2] += Gi[2] * Gj[1] * eps2;

    prod[2][0] += Gi[0] * Gj[2] * eps2;
    prod[2][1] += Gi[1] * Gj[2] * eps2;
    prod[2][2] += Gi[2] * Gj[2] * eps2;

    return prod[0];
}

#undef eu_xx
#undef eu_xy
#undef eu_xz
#undef eu_yx
#undef eu_yy
#undef eu_yz
#undef eu_zx
#undef eu_zy
#undef eu_zz




/***************************************************************************
 * Build matrices *F, *Fu, *B, *Bt, *Fp, *Ap, and *Qp used * by the .
 **************************************************************************/
void
phgNSBuildSolverUMat(NSSolver *ns)
{
    GRID *g = ns->g;
    SIMPLEX *e;
    int i, j, k, l, q, s;
    FLOAT *dt = ns->dt;
    FLOAT Theta = _nsp->Theta, nu = _nsp->nu, Thet1;
    BOOLEAN use_Fu = _nsp->use_Fu;
    int viscosity_type = ns->viscosity_type;
    FLOAT dhx = ns_params->dhx;
    FLOAT dhz = ns_params->dhz;
    LTYPE ltype = ns->ltype;

    SURF_BAS *surf_bas = ns->surf_bas;
    DOF *surf_dof = surf_bas->dof;
    BOOLEAN *rotated = surf_bas->rotated;
    const FLOAT *Trans = DofData(surf_dof);
    assert(surf_dof->type == ns->du->type);

    l = 0; Unused(l); 
    Unused(nu);
#if STEADY_STATE
    assert(fabs(Theta - 1) < 1e-12);
    Thet1 = 0; Unused(Thet1);
    Unused(dt);
#else
    Thet1 = 1 - Theta;
    Unused(dt);
#endif /* STEADY_STATE */

    if (phgUseIsop)
      phgPrintf("  With Isop  ");
    else
      phgPrintf("  WithOut Isop  ");

    /* Removed */
    /* if (_nsp->use_prism_elem) { */
    /* 	phgNSBuildSolverUMatPrism(ns); */
    /* 	return; */
    /* } */

    if (ltype == PICARD)
	phgPrintf("   LinearType: Picard");
    else 
	phgPrintf("   LinearType: Newton");

	//DOF *avg_n = phgDofNew(ns->g, DOF_P2, 3, "avg_n", DofNoAction);
    //get_avg_n(ns->g, avg_n);
    // get averaged normals for ice shelf bottom

    ForAllElements(g, e) {
	int M = ns->du->type->nbas;	/* num of bases of Velocity */
	int N = ns->dp->type->nbas;	/* num of bases of Pressure */
	int order = 2 * DofTypeOrder(ns->du, e) + QUAD_EXTRA_ORDER;
	FLOAT F[M][Dim][M][Dim], Fu[M][M],
	    B[N][M][Dim], Bt[M][Dim][N], C[N][N],
	    bufu[M], bufp[N],
	    rhsu[M][Dim], CM[M][M][Dim], CMt[M][Dim][M], CMD[M][M];
	INT Iu[M][Dim], Iu1[M], Ju[Dim][M], Ip[N], IUn[M];
	QUAD *quad;
	FLOAT dof_u_n1_values[M][Dim];
	FLOAT dof_u_n0_values[M][Dim];
	FLOAT dof_p_n1_values[N];
	FLOAT dof_p_st_values[N];
	FLOAT vol, det, P[Dim][Dim];
	FLOAT *normal;
	const FLOAT *w, *p, *vw, *gu, *vTe; 

	Unused(Iu1);
	for (i = 0; i < M; i++)
	    for (k = 0; k < Dim; k++)
		Ju[k][i] = Iu[i][k] = phgMapE2L(ns->matF->cmap, 0, e, i * Dim + k);
	for (i = 0; i < N; i++)
	    Ip[i] = phgMapE2L(ns->matC->cmap, 0, e, i);

	if (ns_params->sliding_bdry_scheme == 1) /* surf L2 constraint */
	    for (i = 0; i < M; i++)
		IUn[i] = phgMapE2L(ns->matUn->rmap, 0, e, i);


	quad = phgQuadGetQuad3D(order);
	//vw = phgQuadGetDofValues(e, ns->wind, quad);      /* value wind */
	gu = phgQuadGetDofValues(e, ns->gradu[1], quad);  /* grad u */
	phgDofGetElementDatas(ns->u[1], e, dof_u_n1_values[0]);
	if (ns_params->noniter_temp)
	    vTe = phgQuadGetDofValues(e, ns->T[1], quad);	  /* T^{n+1} */
	else
	    vTe = phgQuadGetDofValues(e, ns->T[0], quad);	  /* T^{n} */

	if (ns_params->stab_scheme >= 1) {
	    proj_pres_coef(g, e, quad, P[0]);
	    //SHOW_M(P[0], Dim, Dim);
	}

	
	
	Bzero(F); Bzero(Fu); 
	Bzero(B); Bzero(Bt); Bzero(C);
	Bzero(CM); Bzero(CMt); Bzero(CMD);


	p = quad->points;
	w = quad->weights;
	for (q = 0; q < quad->npoints; q++) {
	    const FLOAT *J;
	    if (phgUseIsop) {
		J = phgGeomGetJacobian_(g, e, p, &det);
		vol = fabs(det) / 6;
	    }
	    else {
		vol = phgGeomGetVolume(g, e);
	    }

	    if (vol == 0.) {
		phgDumpElement(g, e);
		phgError(1, "Element volume is zero!\n");
	    }
	    
	    if (phgUseIsop) {
		static FLOAT gradu_n1_values[Dim*Dim];

		grad_dof_eval(ns->du, dof_u_n1_values[0], e, quad, q, gradu_n1_values);
		
		gu = gradu_n1_values;
	    }

	    nu = get_effective_viscosity(gu, *vTe, 0, viscosity_type);



#if SIMPLE_TEST
	    /* ------------------------------
	     *
	     * Simple test for Stokes eqn convergence test
	     *
	     * ------------------------------ */

	    
	    /* Test func vel type */
	    for (i = 0; i < M; i++) {
		const FLOAT *gi_u = phgQuadGetBasisValues(e, ns->du, i, quad) + q;       /* phi_i */
		const FLOAT *ggi_u = phgQuadGetBasisGradientIsop(e, ns->du, i, quad, q);    /* grad phi_i */

		/* Mat F */
		for (j = 0; j < M; j++) {
		    const FLOAT *gj_u = phgQuadGetBasisValues(e, ns->du, j, quad) + q;       /* phi_j */
		    const FLOAT *ggj_u = phgQuadGetBasisGradientIsop(e, ns->du, j, quad, q);    /* grad phi_i */
		    FLOAT mass = (*gj_u) * (*gi_u);
		    FLOAT diffu = INNER_PRODUCT(ggj_u, ggi_u);

		    Unused(mass);

		    for (k = 0; k < Dim; k++) 
			F[j][k][i][k] += vol*(*w) * diffu;
		}

		/* Mat B & Bt */
		for (j = 0; j < N; j++) {
		    const FLOAT *gj_p = phgQuadGetBasisValues(e, ns->dp, j, quad) + q;       /* psi_j */
		    for (k = 0; k < Dim; k++) {
			FLOAT b = vol*(*w) * (*gj_p) * ggi_u[k];
			B[j][i][k]  -= b; /* divergence */
			Bt[i][k][j] -= b; /* momentum */
		    }
		}
	    } /* end phi_i */


	    /*
	     * P1/P1 pressure local projection stabilization.
	     * See Elman P243.
	     *
	     *  */
	    FLOAT hmax = phgGeomGetDiameter(g, e);
	    
	    for (i = 0; i < N; i++) {
		const FLOAT *gi_p = phgQuadGetBasisValues(e, ns->dp, i, quad) + q;       /* psi_i */
		const FLOAT *ggi_p = phgQuadGetBasisGradientIsop(e, ns->dp, i, quad, q);    /* grad phi_i */
		for (j = 0; j < N; j++) {
		    const FLOAT *gj_p = phgQuadGetBasisValues(e, ns->dp, j, quad) + q;       /* psi_j */
		    const FLOAT *ggj_p = phgQuadGetBasisGradientIsop(e, ns->dp, j, quad, q);    /* grad phi_i */

		    if (ns_params->stab_scheme >= 1) {
			C[i][j] += SIGN_STAB * vol*(*w) * hmax * hmax * INNER_PRODUCT(ggi_p, ggj_p);
		    }
		} /* end bas j */
	    }     /* end bas i */
	    
#else
	    /* ------------------------------
	     *
	     * Ice sheet dynamics
	     *
	     * ------------------------------ */

	    
	    /* Test func vel type */
	    for (i = 0; i < M; i++) {
		const FLOAT *gi_u = phgQuadGetBasisValues(e, ns->du, i, quad) + q;       /* phi_i */
		const FLOAT *ggi_u = phgQuadGetBasisGradientIsop(e, ns->du, i, quad, q);    /* grad phi_i */


		/* Mat F */
		for (j = 0; j < M; j++) {
		    const FLOAT *gj_u = phgQuadGetBasisValues(e, ns->du, j, quad) + q;       /* phi_j */
		    const FLOAT *ggj_u = phgQuadGetBasisGradientIsop(e, ns->du, j, quad, q);    /* grad phi_i */
		    FLOAT mass = (*gj_u) * (*gi_u);
		    FLOAT diffu = INNER_PRODUCT(ggj_u, ggi_u);

		    Unused(mass);


		    const FLOAT *tp = get_gbas_product_stokes(ggi_u, ggj_u, gu, ltype);

		    for (k = 0; k < Dim; k++) 
			for (l = 0; l < Dim; l++) 
			    F[j][l][i][k] += vol*(*w) * EQU_SCALING * nu * tp[k+l*Dim];
		}

		/* Mat B & Bt */
		for (j = 0; j < N; j++) {
		    const FLOAT *gj_p = phgQuadGetBasisValues(e, ns->dp, j, quad) + q;       /* psi_j */
		    for (k = 0; k < Dim; k++) {
			FLOAT b = vol*(*w) * (*gj_p) * ggi_u[k];
			B[j][i][k]  -= b; /* divergence */
			Bt[i][k][j] -= EQU_SCALING * LEN_SCALING * PRES_SCALING * b; /* momentum */
		    }
		}
	    } /* end phi_i */

	    /* Test func pressure type
	     * Some solver doesn't allow zero diagnal entries in Matrix, in that case,
	     *   we perturb diagnal by a small number.
	     * */
	    for (i = 0; i < N; i++) {
		const FLOAT *gi_p = phgQuadGetBasisValues(e, ns->dp, i, quad) + q;       /* psi_i */
		C[i][i] += vol*(*w) * (_nsp->eps_diagP * (*gi_p) * (*gi_p)
				       );
	    }

	    if (ns_params->sliding_bdry_scheme == 1) /* surf L2 constraint */
		for (i = 0; i < M; i++) {
		    CMD[i][i] = 1e-15;
		}

	    /*
	     * P1/P1 pressure local projection stabilization.
	     * See Elman P243.
	     *
	     *  */
	    for (i = 0; i < N; i++) {
		const FLOAT *gi_p = phgQuadGetBasisValues(e, ns->dp, i, quad) + q;       /* psi_i */
		const FLOAT *ggi_p = phgQuadGetBasisGradientIsop(e, ns->dp, i, quad, q);    /* grad phi_i */
		for (j = 0; j < N; j++) {
		    const FLOAT *gj_p = phgQuadGetBasisValues(e, ns->dp, j, quad) + q;       /* psi_j */
		    const FLOAT *ggj_p = phgQuadGetBasisGradientIsop(e, ns->dp, j, quad, q);    /* grad phi_i */

		    const FLOAT alpha_stab = ns_params->stab_alpha;
		    FLOAT nu_stab = 1.;
		    if (ns_params->stab_nu < 0) {
			nu_stab = nu;
		    }
		    else {
			nu_stab = ns_params->stab_nu;
		    }
		    
		    if (ns_params->stab_scheme == 0) {
			C[i][j] += SIGN_STAB * vol*(*w) * 1./(EQU_SCALING * nu_stab) * alpha_stab
			    * ((*gj_p - 1./4.) * (*gi_p - 1./4.));
			
		    }
		    else if (ns_params->stab_scheme >= 1) {
			FLOAT Jgi[3];

			const FLOAT *iJ_ = P[0];
			MAT3_MV(iJ_, ggi_p, Jgi);
			    
			C[i][j] += SIGN_STAB * vol*(*w) * 1./(EQU_SCALING * nu_stab) * alpha_stab
			    * INNER_PRODUCT(Jgi, ggj_p);
		    }

		} /* end bas j */
	    }     /* end bas i */
	    
#endif	/* SIMPLE_TEST */
		

	    /* Next quad point */
	    //vw += Dim;
	    gu += Dim*Dim;
	    vTe++;
	    w++; p += Dim+1;
	}


	/* slip boundary */
	for (s = 0; s < NFace; s++) {
	    if (e->bound_type[s] & SLIP_BDRY) {
		/* ------------------------------
		 *
		 *
		 * Slip boundary
		 *
		 *
		 * ------------------------------ */
		int v0, v1, v2;
		int nbas_face = NbasFace(ns->du);
		SHORT bases[nbas_face];
		FLOAT lambda[Dim + 1], area, x, y, z, beta;
		order = 2 * DofTypeOrder(ns->du, e) + QUAD_EXTRA_ORDER;
		    
		phgDofGetBasesOnFace(ns->du, e, s, bases);
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
		    lambda[v0] = *(p++);
		    lambda[v1] = *(p++);
		    lambda[v2] = *(p++);

		    if (phgUseIsop) {
			area = phgGeomGetFaceArea_(g, e, s, lambda, normal);
		    }
		    
		    phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);
#if 0
		    func_beta(x, y, z, &beta);
#else
		    phgDofEval(ns->beta, e, lambda, &beta);
#endif

		    const FLOAT *gi_u = 
			ns->du->type->BasFuncs(ns->du, e, 0, -1, lambda);
		    for (i = 0; i < nbas_face; i++) {
			int ii = bases[i];
			for (j = 0; j < nbas_face; j++) { 
			    int jj = bases[j];
			    FLOAT mass_face = area*(*w) * beta * (gi_u[jj])*(gi_u[ii])
				* EQU_SCALING * LEN_SCALING;
			    /* See scaling note for slip bdry in build RHS. */

			    for (k = 0; k < Dim; k++) {
				F[ii][k][jj][k] -= SIGN_FRICTION * mass_face;
			    }

			    /* Penalty */
			    if (ns_params->sliding_bdry_scheme == 1) /* surf L2 constraint */
				for (k = 0; k < Dim; k++) {
				    CM[ii][jj][k]  += area*(*w) * normal[k] * (gi_u[jj])*(gi_u[ii]);
				    CMt[jj][k][ii] += area*(*w) * normal[k] * (gi_u[jj])*(gi_u[ii]);
			    }

				
			} /* end of bas_j */
		    }     /* end of bas_i */
		    w++;
		}         /* end of quad point */
	    }		/* end of face outflow */
            else if (e->bound_type[s] & BC_FLOAT)
            {
		/* ------------------------------
		 *
		 *
		 * Ice shelf boundary, by Tong
		 *
		 *
		 * ------------------------------ */
		
                FLOAT lambda[Dim+1], lambda_x, lambda_y, lambda_z;
                INT quad_order = 5;
                FLOAT Ns;
                QUAD *quad = phgQuadGetQuad2D(quad_order);
                FLOAT *w = quad->weights;
                FLOAT *p = quad->points;
		FLOAT *normal = phgGeomGetFaceOutNormal(g, e, s);

                FLOAT avg_n_v[Dim];

                FLOAT nx, ny, nz, x, y, z;

                FLOAT area = phgGeomGetFaceArea(g, e, s);
                //FLOAT normal[Dim];

                INT v0 = GetFaceVertex(s, 0);
                INT v1 = GetFaceVertex(s, 1);
                INT v2 = GetFaceVertex(s, 2);
                lambda[s] = 0;

                FLOAT dt1 = ns->dt[0];
                INT M_face = NbasFace(ns->du);
                SHORT bas_idx_e[M_face];

                phgDofGetBasesOnFace(ns->du, e, s, bas_idx_e);

                for (q = 0; q < quad->npoints; q++) {
                    lambda[v0] = *(p++);
                    lambda[v1] = *(p++);
                    lambda[v2] = *(p++);
                    
		    phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);

                    //phgDofEval(avg_n, e, lambda, avg_n_v);
                    //normal[0] = avg_n_v[0];
                    //normal[1] = avg_n_v[1];
                    //normal[2] = avg_n_v[2];
		    if (fabs(normal[2]) < 1.0e-8) {
			Ns = 1.0e50;
			phgWarning("Base normal direction nearly horizontal!!!\n");
		    }
		    else {
			Ns = sqrt(1 + SQUARE(normal[0]/normal[2])+SQUARE(normal[1]/normal[2]));
		    }
                    
                    const FLOAT *bas = ns->du->type->BasFuncs(ns->du, e, 0, -1, lambda);

		    if (z < 0)
			for (i = 0; i < M_face; i++)  {
			    //if (phgDofGetElementBoundaryType(ns->du, e, i) & BC_BOTTOM_GRD)
			    //    continue;

			    INT i_e = bas_idx_e[i];
			    const FLOAT *bas_i = phgQuadGetBasisValues(e, ns->du, i_e, quad);
			    for (j = 0; j < M_face; j++) {
				//if (phgDofGetElementBoundaryType(ns->du, e, j) & BC_BOTTOM_GRD)
				//    continue;

				INT j_e = bas_idx_e[j];
				const FLOAT *bas_j = phgQuadGetBasisValues(e, ns->du, j_e, quad);
                            
				const FLOAT *bas_dot_n = get_bas_dot_normal(bas, normal, i_e, j_e);
				for (k = 0; k < Dim; k++) {
				    //F[j_e][k][i_e][k] += area*w[q]*RHO_WAT*GRAVITY*bas[i_e]*bas[j_e]*normal[k]*normal[k]*Ns*dt1*EQU_SCALING;
#  if 1
				    for (l = 0; l < Dim; l++)  {
					F[j_e][l][i_e][k] += area  * w[q] * RHO_WATER * GRAVITY *
						bas_dot_n[k+l*Dim] * Ns * dt1 * EQU_SCALING * LEN_SCALING;
					/* See scaling note for calving front in build RHS. */
				    }
#  endif

				}  /* end dim */
			    }      /* end bas face j */
			}	   /* end bas face i */
                }	       /* end quad point */
            }		       /* end ice shelf face */

	    
	}		/* end of all outflow face in element */


	/* Rotate bases */
	if (ns_params->sliding_bdry_scheme == 0) /* Dirich */
	    for (i = 0; i < M; i++) {
		INT id = phgDofMapE2D(surf_dof, e, i * (Dim*Dim)) / (Dim*Dim);
		if (!rotated[id])
		    continue;	
		const FLOAT *trans = Trans + id*(Dim*Dim);


		trans_left_3d(&F[i][0][0][0], Dim*M, Dim*M, trans);
		trans_rightT_3d(&F[0][0][i][0], Dim*M, Dim*M, trans);

		trans_left_3d(&Bt[i][0][0], N, N, trans);
		trans_rightT_3d(&B[0][i][0], N, Dim*M, trans);
	    }


	/* Global Matrix */
	/* Mat u-p Block (1, *) */
	for (i = 0; i < M; i++) {
	    /* du = 0 at Dirichlet boundary */
	    for (k = 0; k < Dim; k++) {
		if (phgDofDirichletBC_(ns->du, e, i*Dim+k, NULL, bufu, &rhsu[i][0],
				       DOF_PROJ_NONE)) {

		    phgMatAddEntries(ns->matF, 1, Iu[i] + k, M, Ju[k], bufu);

		}
		else {
		    phgMatAddEntries(ns->matF, 1, Iu[i] + k, M*Dim, Iu[0],
				     &(F[i][k][0][0]));
		    phgMatAddEntries(ns->matBt, 1, &Iu[i][k], N, Ip,
				     &Bt[i][k][0]);

		    if (ns_params->sliding_bdry_scheme == 1) /* surf L2 constraint */
			phgMatAddEntries(ns->matUnT, 1, &Iu[i][k], M, IUn,
					 &CMt[i][k][0]);
		}
	    }
	}

	/* Mat u-p (1, *) */
	for (i = 0; i < N; i++) {
	    if (phgDofDirichletBC(ns->dp, e, i, NULL, bufp, NULL, DOF_PROJ_NONE)) {
		phgMatAddEntries(ns->matC, 1, Ip + i, N, Ip, bufp);
	    } else {
		phgMatAddEntries(ns->matB, 1, Ip + i, M * Dim, Iu[0],
				      &B[i][0][0]);
		phgMatAddEntries(ns->matC, 1, Ip + i, N, Ip,
				      &C[i][0]);
	    }
	}

	if (ns_params->sliding_bdry_scheme == 1) /* surf L2 constraint */
	    for (i = 0; i < M; i++) {
		phgMatAddEntries(ns->matUn, 1, IUn + i, M * Dim, Iu[0],
				 &CM[i][0][0]);
		phgMatAddEntries(ns->matUnD, 1, IUn + i, M, IUn,
				 &CMD[i][0]);
	    }
    }				/* end element */

	//phgDofFree(&avg_n);

    /* mat check */
#define MAT_CHECK_DUP(mat)    {					\
	MAT_ROW *row = mat->rows;				\
	for (i = 0; i < mat->rmap->nlocal; i++, row++) {	\
	    int k_, ncol = row->ncols;				\
	    INT cols[ncol];					\
	    for (j = 0; j < ncol; j++) {			\
		cols[j] = row->cols[j];				\
		for (k_ = 0; k_ < j; k_++)			\
		    assert(cols[k_] != cols[j]);		\
	    }							\
	}							\
    }
    /* MAT_CHECK_DUP(ns->matF); */
    /* MAT_CHECK_DUP(ns->matB); */
    /* MAT_CHECK_DUP(ns->matBt); */
    /* MAT_CHECK_DUP(ns->matC); */

#define MAT_CHECK_ZERO(mat)    {					\
	MAT_ROW *row = mat->rows;					\
	for (i = 0; i < mat->rmap->nlocal; i++, row++) {		\
	    int k_, ncol = row->ncols;					\
	    if (ncol <= 0) {						\
		phgInfo(0, "mat %s row %d is zero.\n", #mat, i);	\
	    }								\
	}								\
    }
    MAT_CHECK_ZERO(ns->matF);
    /* MAT_CHECK_ZERO(ns->matB); */
    /* MAT_CHECK_ZERO(ns->matBt); */
    /* MAT_CHECK_ZERO(ns->matC); */


    if ( DUMP_MAT_VEC) {
	phgPrintf("dumping F,B,Bt,C\n");
	phgMatDumpMATLAB(ns->matF, "F", "F_.m");
	phgMatDumpMATLAB(ns->matB, "B", "B_.m");
	phgMatDumpMATLAB(ns->matBt, "Bt", "Bt_.m");
	phgMatDumpMATLAB(ns->matC, "C", "C_.m");

	if (ns_params->sliding_bdry_scheme == 1) {
	    phgMatDumpMATLAB(ns->matUn, "Un", "Un_.m");
	    phgMatDumpMATLAB(ns->matUnT, "UnT", "UnT_.m");
	}
    }


    phgMatAssemble(ns->matNS);
    {
	MAT *mat = ns->matNS;
	phgInfo(0, "local rows: %d\n", mat->rmap->nlocal);
	for (i = 0; i < mat->rmap->nlocal; i++, j++) {
	    const MAT_ROW *row;
	    row = phgMatGetRow(mat, i);
	    if (row->ncols <= 0) {
		phgInfo(0, "zero rows %d\n", i);
	    }
	}
    }

    /* Exit on checking matrix. */
    if (0) { //&& viscosity_type) {
    	phgFinalize();
    	exit(1);
    }
    return;
}



