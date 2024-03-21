/*
 *
 * First-Order 
 *
 * 
 *  */
#include "ins.h"
#include <unistd.h>

#define SIGN_FRICTION -1
#define SMALL_BETA  1e-8

/*
 *
 *  Note: SIA change to be done for z_b != 0, 
 *    1. gL->verts[3] = zbot and gL->verts[4] = ztop;
 *    2. Where height is updated, the 2D verts coord need to be updated
 *    3. grad S ==> grad height ?
 *
 *  */


static int verb_core = 0;
#define CHECK_INTU 0


/*	 Paterson & budd (1982) */
#  define	 ARRHENIUS_T	263.15
#  define	 ARRHENIUS_A0   3.61E-13     
#  define	 ARRHENIUS_A1	1.73E3
#  define	 ARRHENIUS_Q0	6.0E4
#  define	 ARRHENIUS_Q1	13.9E4


#define GRAD_SURF_SCHEME 0
#if GRAD_SURF_SCHEME == 0
/* gradS, triangle wise */
static FLOAT *gradS_tri = NULL, *gradS_tri0 = NULL;
static FLOAT *gradS_vert = NULL;
#elif GRAD_SURF_SCHEME == 1
/* use (S, grad phi)  */
static DOF *dofS = NULL;
#elif GRAD_SURF_SCHEME == 2
/* gradS, point wise */
static FLOAT *gradS_tri = NULL, *gradS_tri0 = NULL;
static FLOAT *gradS_vert = NULL;
static DOF *dof_gradS = NULL;
#endif
static DOF *dof_uv = NULL, *dof_du = NULL, *dof_w = NULL, *grad_uv = NULL;
static SOLVER *solver_uv = NULL;



FLOAT * 
get_gbas_product_fo(const FLOAT *gi, const FLOAT *gj,
		 const FLOAT *gu, LTYPE ltype) 
{
    /*
     *            phi_j, Dofs
     *  phi_i [                  ]
     *  test  [                  ]
     *
     *  */
    static FLOAT prod[2][2], jac[2][2];
    FLOAT sxx, sxy, sxz, syx, syy, syz, szx, szy, szz;
    FLOAT eps, eps2, a0, a1, b, c0, c1;
    FLOAT ux, uy, uz, vx, vy, vz;
    int k;

    sxx = gi[0] * gj[0];
    syy = gi[1] * gj[1];
    szz = gi[2] * gj[2];
    sxy = gi[0] * gj[1];
    syx = gi[1] * gj[0];
    sxz = gi[0] * gj[2];
    syz = gi[1] * gj[2];
    szx = gi[2] * gj[0];
    szy = gi[2] * gj[1];

    /* Picard term */
    prod[0][0] = 2 * sxx + 0.5 * (syy + szz); 
    prod[1][0] = sxy + 0.5 * syx; 
    prod[0][1] = syx + 0.5 * sxy; 
    prod[1][1] = 2 * syy + 0.5 * (sxx + szz); 

    if (ltype == PICARD) {
	return prod[0];
    }

    /* Newton term */
    ux = gu[0] / LEN_SCALING; 
    uy = gu[1] / LEN_SCALING; 
    uz = gu[2] / LEN_SCALING;
    vx = gu[3] / LEN_SCALING; 
    vy = gu[4] / LEN_SCALING; 
    vz = gu[5] / LEN_SCALING;

    eps = sqrt(ux*ux + vy*vy + ux*vy
	       + 0.25 * pow(uy + vx, 2)
	       + 0.25 * (uz*uz + vz*vz));
    
    if (eps < MIN_EFFECTIVE_STRAIN) 
	eps = MIN_EFFECTIVE_STRAIN;

    eps2 = - (1./3.) / (eps*eps);


    a0 = 2*ux + vy;
    a1 = 2*vy + ux;
    b  = .5 * (uy + vx);
    c0 = .5 * uz;	
    c1 = .5 * vz;		 

    prod[0][0] += (  a0 * (a0 * sxx + b * sxy + c0 * sxz)
    		    + b  * (a0 * syx + b * syy + c0 * syz)
    		    + c0 * (a0 * szx + b * szy + c0 * szz)) * eps2;
    prod[1][0] += (  a0 * (a1 * sxy + b * sxx + c1 * sxz)
    		    + b  * (a1 * syy + b * syx + c1 * syz)
    		    + c0 * (a1 * szy + b * szx + c1 * szz)) * eps2;
    prod[0][1] += (  a1 * (a0 * syx + b * syy + c0 * syz)
    		    + b  * (a0 * sxx + b * sxy + c0 * sxz)
    		    + c1 * (a0 * szx + b * szy + c0 * szz)) * eps2;
    prod[1][1] += (  a1 * (a1 * syy + b * syx + c1 * syz)
    		    + b  * (a1 * sxy + b * sxx + c1 * sxz)
    		    + c1 * (a1 * szy + b * szx + c1 * szz)) * eps2;

    return prod[0];
}


static void 
build_matrix(NSSolver *ns)
{
    GRID *g = ns->g;
    SIMPLEX *e;
    int i, j, k, l, q, s;
    FLOAT *dt = ns->dt;
    FLOAT Theta = ns_params->Theta, nu = ns_params->nu, Thet1;
    int viscosity_type = ns->viscosity_type;
    LTYPE ltype = ns->ltype;

    /* Lateral BC */
    SURF_BAS *surf_bas = ns->surf_bas_2d;
    DOF *surf_dof = surf_bas->dof;
    BOOLEAN *rotated = surf_bas->rotated;
    const FLOAT *Trans = DofData(surf_dof);
    assert(surf_dof->type == ns->du->type);


    /* Removed */
    /* if (ns_params->use_prism_elem) { */
    /* 	buildFOMatPrism(ns, dof_uv, dof_du, solver_uv); */
    /* 	return; */
    /* } */
    
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

    phgPrintf("   DB_mask: %d [", __LINE__);
    for (k = 0; k < 2; k++)
	phgPrintf("%d(%d) ", dof_du->DB_masks[k], dof_uv->DB_masks[k]);
    phgPrintf("]   \n");

    if (ltype == PICARD)
	phgPrintf("   LinearType: Picard");
    else 
	phgPrintf("   LinearType: Newton");

    ForAllElements(g, e) {
	int M = dof_du->type->nbas;	/* num of bases of Velocity */
	int order = 2 * DofTypeOrder(dof_du, e);
	FLOAT F[M][2][M][2], 
	    bufu[M], rhsu[M][2];
	INT Iu[M][2], Ju[2][M];
	QUAD *quad;
	FLOAT vol, det;
	const FLOAT *w, *p, *gu, *vTe;

	for (i = 0; i < M; i++)
	    for (k = 0; k < 2; k++)
		Ju[k][i] = Iu[i][k] = phgMapE2L(solver_uv->rhs->map, 0, e, i * 2 + k);

	quad = phgQuadGetQuad3D(order);
	gu = phgQuadGetDofValues(e, grad_uv, quad);  /* grad u */
	if (ns_params->noniter_temp)
	    vTe = phgQuadGetDofValues(e, ns->T[1], quad);	  /* T^{n+1} */
	else
	    vTe = phgQuadGetDofValues(e, ns->T[0], quad);	  /* T^{n} */

	Bzero(F); 

	p = quad->points;
	w = quad->weights;
	for (q = 0; q < quad->npoints; q++) {
	    phgGeomGetCurvedJacobianAtLambda(g, e, p, &det);
	    vol = fabs(det / 6.);


	    /* Test func vel type */
	    for (i = 0; i < M; i++) {
		const FLOAT *gi_u = phgQuadGetBasisValues(e, ns->du, i, quad) + q;       /* phi_i */
		const FLOAT *ggi_u = phgQuadGetBasisCurvedGradient(e, ns->du, i, quad, q);    /* grad phi_i */

		/* Mat F */
		for (j = 0; j < M; j++) {
		    const FLOAT *gj_u = phgQuadGetBasisValues(e, ns->du, j, quad) + q;       /* phi_j */
		    const FLOAT *ggj_u = phgQuadGetBasisCurvedGradient(e, ns->du, j, quad, q);    /* grad phi_i */
		    FLOAT mass = (*gj_u) * (*gi_u);
		    FLOAT diffu = INNER_PRODUCT(ggj_u, ggi_u);


		    Unused(mass);
		    nu = get_effective_viscosity(gu, *vTe, 0, viscosity_type);

		    const FLOAT *tp = get_gbas_product_fo(ggi_u, ggj_u, gu, ltype);

		    for (k = 0; k < 2; k++) 
			for (l = 0; l < 2; l++) 
			    F[j][l][i][k] += vol*(*w) * EQU_SCALING * nu * tp[k + l*2];

		}
	    } /* end phi_i */

	    gu += Dim*2;
	    vTe++;
	    w++; p += Dim+1;
	}


	/* Sliding boundary */
	for (s = 0; s < NFace; s++) {
	    if (e->bound_type[s] & (SLIP_BDRY | BC_FLOAT) ) {
		int v0, v1, v2;
		int nbas_face = NbasFace(ns->du);
		SHORT bases[nbas_face];
		FLOAT lambda[Dim + 1], area, x, y, z, beta;
		order = 2 * DofTypeOrder(ns->du, e);
		    
		phgDofGetBasesOnFace(ns->du, e, s, bases);
		v0 = GetFaceVertex(s, 0);
		v1 = GetFaceVertex(s, 1);
		v2 = GetFaceVertex(s, 2);
		lambda[s] = 0.;

		area = phgGeomGetFaceArea(g, e, s);
		quad = phgQuadGetQuad2D(order);

		p = quad->points;
		w = quad->weights;
		for (q = 0; q < quad->npoints; q++) {
		    lambda[v0] = *(p++);
		    lambda[v1] = *(p++);
		    lambda[v2] = *(p++);
			
		    phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);

		    if (e->bound_type[s] & SLIP_BDRY) {
#if 0
			func_beta(x, y, z, &beta);
#else
			phgDofEval(ns->beta, e, lambda, &beta);
#endif
		    }
		    else {
			beta = SMALL_BETA;
		    }

		    const FLOAT *gi_u = 
			    ns->du->type->BasFuncs(ns->du, e, 0, -1, lambda);
		    for (i = 0; i < nbas_face; i++) {
			int ii = bases[i];
			for (j = 0; j < nbas_face; j++) { 
			    int jj = bases[j];
			    FLOAT mass_face = area*(*w) * beta * (gi_u[jj])*(gi_u[ii])
					      * EQU_SCALING * LEN_SCALING;

			    for (k = 0; k < 2; k++) {
				F[ii][k][jj][k] -= SIGN_FRICTION * mass_face;
			    }
			} /* end of bas_j */
		    }     /* end of bas_i */
		    w++;
		}		/* end of quad point */
	    }		/* end of face outflow */
	}


	/*
	 * Lateral rotation
	 *  */
	for (i = 0; i < M; i++) {
	    const int rdim = 2;
	    INT id = phgDofMapE2D(surf_dof, e, i * (rdim*rdim)) / (rdim*rdim);
	    if (!rotated[id])
		continue;	
	    const FLOAT *trans = Trans + id*(rdim*rdim);

	    trans_left_2d(&F[i][0][0][0], rdim*M, rdim*M, trans);
	    trans_rightT_2d(&F[0][0][i][0], rdim*M, rdim*M, trans);
	}


	
	/* Global Matrix */
	/* Mat u-p Block (1, *) */
	for (i = 0; i < M; i++) {
	    /* du = 0 at Dirichlet boundary */
	    for (k = 0; k < 2; k++) {
		if (phgDofDirichletBC_(dof_du, e, i*2+k, NULL, bufu, &rhsu[i][0],
				       DOF_PROJ_NONE)) {
		    phgMatAddEntries(solver_uv->mat, 1, Iu[i] + k, M, Ju[k], bufu);
		}
		else {
		    phgMatAddEntries(solver_uv->mat, 1, Iu[i] + k, M*2, Iu[0],
				     &(F[i][k][0][0]));
		}
	    }
	}

    }				/* end element */



    if (DUMP_MAT_VEC) {
	phgPrintf("dumping FO mat\n");
	phgMatDumpMATLAB(solver_uv->mat, "A", "A_.m");
    }

    /* Exit on checking matrix. */
    if (0) { //&& viscosity_type) {
    	phgFinalize();
    	exit(1);
    }
    return;
    
}



static void
build_rhs(NSSolver *ns)
{
    GRID *g = ns->g;
    LAYERED_MESH *gL = ns->gL;
    int i, ii, jj, k, l, q, s;
    FLOAT *dt = ns->dt;
    BOOLEAN tstep_minus = (ns->u[-1] != NULL);
    VEC *vec_rhs = phgMapCreateVec(solver_uv->rhs->map, 1);
    FLOAT Theta = ns_params->Theta, nu = ns_params->nu, Thet1;
    int viscosity_type = ns->viscosity_type;
    FLOAT Vol = 0., Vol0;

    // Lateral BC
    SURF_BAS *surf_bas = ns->surf_bas_2d;
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

    assert(dof_du->type == ns->u[1]->type);
    nu_max = -1e100;
    nu_min = +1e100;

    /* Removed */
    /* if (ns_params->use_prism_elem) { */
    /* 	buildFORHSPrism(ns, dof_uv, dof_du, solver_uv); */
    /* 	return; */
    /* } */
    
    phgVecDisassemble(vec_rhs);
    BOTTOM_FACE *fb = gL->face_bot;
    for (ii = 0; ii < gL->nface_bot; ii++, fb++) {
	TRIA *t = gL->trias + fb->tria;

	for (jj = 0; jj < fb->ne; jj++) { /* tets */
	    SIMPLEX *e = fb->elems[jj];
	    int M = dof_du->type->nbas;	/* num of bases of Velocity */
	    int order = 2 * DofTypeOrder(dof_du, e);
	    FLOAT bufu[M], rhsu[M][2];
	    INT Iu[M][2];
	    QUAD *quad;
	    FLOAT vol, area, det, *vf; 
	    const FLOAT *w, *p, *normal, *vS, *gS,
		*gu, *vTe;


	    quad = phgQuadGetQuad3D(order);
	    if (ns_params->noniter_temp)			                /* nonlinear temp iter */
		vTe = phgQuadGetDofValues(e, ns->T[1], quad);	        /* T^{n+1} */
	    else
		vTe = phgQuadGetDofValues(e, ns->T[0], quad);	        /* T^{n} */

#if STEADY_STATE || TIME_DEP_NON
	    gu = phgQuadGetDofValues(e, grad_uv, quad);      /* grad u^{n} */
#else
	    TIME_DEP_LINEAR_ENTRY; /* Unavailable */
#endif /* STEADY_STATE || TIME_DEP_NON */

#if GRAD_SURF_SCHEME == 1
	    vS = phgQuadGetDofValues(e, dofS, quad);
#elif GRAD_SURF_SCHEME == 2
	    gS = phgQuadGetDofValues(e, dof_gradS, quad);
#endif
	    

	    Unused(l);

	    if (ns_params->extern_force) {
		/* cache f values,
		 * only time_{n} */
		for (l = 1; l < 2; l++) {
		    const FLOAT *cache;
		    size_t cache_size;
		    setFuncTime(ns->time[l]); /* set static time in ins-test.c */

		    /* cache f */
		    cache_size = Dim * quad->npoints * sizeof(FLOAT);
		    cache = phgQuadGetFuncValues(g, e, Dim, func_f, quad);
		    vf = phgAlloc(cache_size);
		    memcpy(vf, cache, cache_size);

		    phgQuadGetFuncValues(NULL, NULL, 0, NULL, NULL); /* clear cache */
		}
	    }

	    /* Global Matrix */
	    Bzero(rhsu); 
    
	    p = quad->points;
	    w = quad->weights;
	    for (q = 0; q < quad->npoints; q++) {
		phgGeomGetCurvedJacobianAtLambda(g, e, p, &det);
		vol = fabs(det / 6.);
		Vol += vol*(*w);

		/* rhs u */
		for (i = 0; i < M; i++) {
		    /* interior node or Neumann */
		    const FLOAT *gi = phgQuadGetBasisValues(e, ns->du, i, quad) + q;       /* phi_i */
		    const FLOAT *ggi = phgQuadGetBasisCurvedGradient(e, ns->du, i, quad, q);    /* grad phi_i */

		    for (k = 0; k < 2; k++) {

			nu = get_effective_viscosity(gu, *vTe, 0, viscosity_type);

			FLOAT ux, uy, uz, vx, vy, vz;
			ux = gu[0]; // LEN_SCALING; 
			uy = gu[1]; // LEN_SCALING; 
			uz = gu[2]; // LEN_SCALING;
			vx = gu[3]; // LEN_SCALING; 
			vy = gu[4]; // LEN_SCALING; 
			vz = gu[5]; // LEN_SCALING;

			if (k == 0)
			    rhsu[i][k] += -vol*(*w) * nu * (ggi[0]*(2*ux+vy) + 0.5*ggi[1]*(uy+vx) + 0.5*ggi[2]*uz
							    ) * EQU_SCALING;
			else
			    rhsu[i][k] += -vol*(*w) * nu * (0.5*ggi[0]*(uy+vx) + ggi[1]*(ux+2*vy) + 0.5*ggi[2]*vz
							    ) * EQU_SCALING;

			{ 
			    const FLOAT rho = RHO_ICE;
			    const FLOAT grav = GRAVITY;
			    const FLOAT a = SEC_PER_YEAR;
#if GRAD_SURF_SCHEME == 0
			    const FLOAT f = rho*grav * gradS_tri[2*fb->tria + k]
				* LEN_SCALING2 * EQU_SCALING; 
			    Unused(a);
			    rhsu[i][k] += vol*(*w) * (- f * (*gi) 
						      ); /* right */
#elif GRAD_SURF_SCHEME == 1
			    rhsu[i][k] += vol*(*w) 
				* rho*grav * (*vS) * ggi[k]
				* LEN_SCALING2 * EQU_SCALING;
#elif GRAD_SURF_SCHEME == 2
			    rhsu[i][k] += vol*(*w) 
				* rho*grav * (- gS[k]) * (*gi)
				* LEN_SCALING2 * EQU_SCALING;
#endif

			}

			if (ns_params->extern_force)
			    rhsu[i][k] += vol*(*w) * (*(vf+k) * (*gi)
						      ); /* right */

		    }
		}

		gu += Dim*2;
		//vw += Dim;
		if (ns_params->extern_force) {
		    vf += Dim;
		}
		vTe++;
#if GRAD_SURF_SCHEME == 1
		vS++;
#elif GRAD_SURF_SCHEME == 2
		gS += 2;
#endif
		w++; p += Dim + 1;
	    }

	    if (ns_params->extern_force) {
		phgFree(vf);
	    }


	    /* Sliding boundary */
	    for (s = 0; s < NFace; s++) {
		if (e->bound_type[s] & (SLIP_BDRY | BC_FLOAT) ) {
		    int v0, v1, v2;
		    int nbas_face = NbasFace(ns->du);
		    SHORT bases[nbas_face];
		    FLOAT lambda[Dim + 1], x,y,z, beta;
		    order = 2 * DofTypeOrder(ns->du, e);

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
			FLOAT vu[2];
			lambda[v0] = *(p++);
			lambda[v1] = *(p++);
			lambda[v2] = *(p++);
			
			phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);
			phgDofEval(dof_uv, e, lambda, vu);

			if (e->bound_type[s] & SLIP_BDRY) {
#if 0
			    func_beta(x, y, z, &beta);
#else
			    phgDofEval(ns->beta, e, lambda, &beta);
#endif
			}
			else {
			    beta = SMALL_BETA;
			}

			const FLOAT *gi_u = 
				ns->du->type->BasFuncs(ns->du, e, 0, -1, lambda);

			for (i = 0; i < nbas_face; i++) {
			    int ii = bases[i];
			    for (k = 0; k < 2; k++) {
				rhsu[ii][k] += SIGN_FRICTION * area*(*w) * beta * vu[k] * (gi_u[ii])
					       * EQU_SCALING * LEN_SCALING;
			    }
			}     /* end of bas_i */
			w++;
		    }   /* end of quad point */
		}	/* end face type */
		else if (e->bound_type[s] & BC_FRONT ) {
		    int v0, v1, v2;
		    int nbas_face = NbasFace(ns->du);
		    SHORT bases[nbas_face];
		    FLOAT lambda[Dim + 1], x,y,z;
		    order = 2 * DofTypeOrder(ns->du, e);

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
			FLOAT vp_water, vp_ice, vdepth;
			lambda[v0] = *(p++);
			lambda[v1] = *(p++);
			lambda[v2] = *(p++);

			phgDofEval(ns->depth_P1, e, lambda, &vdepth);
			phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);

			const FLOAT rho_w = RHO_WATER;
			const FLOAT rho_i = RHO_ICE;
			const FLOAT grav = GRAVITY;

			if (z > 0)
			    vp_water = 0.;
			else
			    vp_water = rho_w * grav * (-z) * LEN_SCALING;
			vp_ice = rho_i * grav * (vdepth) * LEN_SCALING;

			const FLOAT *gi_u =
				ns->du->type->BasFuncs(ns->du, e, 0, -1, lambda);

			for (i = 0; i < nbas_face; i++) {
			    int ii = bases[i];
			    for (k = 0; k < 2; k++) {
				rhsu[ii][k] += area*(*w) * (vp_ice - vp_water) * normal[k] * (gi_u[ii])
					       * EQU_SCALING * LEN_SCALING;
			    }
			}     /* end of bas_i */
			w++;
		    }   /* end of quad point */
		}	/* end face type */
	    }		/* end face */

	    /* Lateral rotation */
	    for (i = 0; i < M; i++) {
		const int rdim = 2 ;
		INT id = phgDofMapE2D(surf_dof, e, i * (rdim*rdim)) / (rdim*rdim);
		if (!rotated[id])
		    continue;	
		const FLOAT *trans = Trans + id*(rdim*rdim);

		trans_left_2d(&rhsu[i][0], 1, 1, trans);
	    }

	    
	    /* Map: Element -> system */
	    for (i = 0; i < M; i++)
		for (k = 0; k < 2; k++)
		    Iu[i][k] = phgMapE2L(solver_uv->rhs->map, 0, e, i * 2 + k);

	    /* set velocity dirichlet bdry */
	    FLOAT tmp[2];
	    for (i = 0; i < M; i++)
		for (k = 0; k < 2; k++) {
		    if (phgDofDirichletBC_(dof_du, e, i*2+k, NULL, bufu, tmp,
					   DOF_PROJ_NONE)) {
			rhsu[i][k] = 0.;
		    }
		}

	    /* Global res */
	    phgVecAddEntries(vec_rhs, 0, M * 2, Iu[0], &rhsu[0][0]);
	}				/* end element */
    }
    
    phgVecAssemble(vec_rhs);
    phgVecAssemble(solver_uv->rhs);
    phgVecAXPBY(1., vec_rhs, 0, &solver_uv->rhs);
    solver_uv->rhs_updated = FALSE;

    if (DUMP_MAT_VEC) {
	phgPrintf("Dumping rhs\n");
	phgVecDumpMATLAB(solver_uv->rhs, "b", "b_.m");
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

    return;    
}



static void 
copy_solution(DOF *u3, DOF *uv, /* DOF *w, */ 
	      BOOLEAN forward)
{
    INT i, ndof = DofGetDataCount(u3) / Dim;
    FLOAT *vU, *vu; //*vw;

    vU = u3->data;
    vu = uv->data;
    //vw = w->data;
    if (forward)
	for (i = 0; i < ndof; i++) {
	    *(vu++) = *(vU++);
	    *(vu++) = *(vU++);
	    vU++;
	    //*(vw++) = *(vU++);
	}
    else 
	for (i = 0; i < ndof; i++) {
	    *(vU++) = *(vu++) ;
	    *(vU++) = *(vu++) ;
	    vU++;
	    //*(vU++) = *(vw++) ;
	}

    return;
}

#if 0
#  define CHECK_DBMASK							\
    {									\
	phgPrintf("   DB_mask: %d [", __LINE__);			\
	for (k = 0; k < 2; k++)						\
	    phgPrintf("%d(%d) ", dof_du->DB_masks[k], dof_uv->DB_masks[k]); \
	phgPrintf("]   \n");						\
    }
#else
# define CHECK_DBMASK
#endif


void 
core_FO(NSSolver *ns, int tstep)
{
    GRID *g = ns->g;
    LAYERED_MESH *gL = ns->gL;
    int i, ii, j, k;
    FILE *fp;
    double tt[3], tt1[3];
    char vtk_file[1000];

#if USE_SLIDING_BC    
    SURF_BAS *surf_bas = ns->surf_bas_2d;
    DOF *surf_dof = surf_bas->dof;
    const FLOAT *Trans = DofData(surf_dof);
    FLOAT *vU = DofData(ns->u[1]);
#endif    



#if GRAD_SURF_SCHEME == 0
    /* ------------------------------------------------------------
     * 
     *
     * Use grad S, triangle wize
     *
     *
     * ------------------------------------------------------------ */
    phgPrintf("Use grad S(tri).\n");
    if (verb_core) {
	char name[100];
	sprintf(name, "grads_%02d.dat", phgRank);
	fp = fopen(name, "w");
    }

    if (gradS_tri == NULL) {
	PHG_CALLOC(gradS_tri, 2*gL->ntria);
	PHG_CALLOC(gradS_tri0, 2*gL->ntria);
	PHG_CALLOC(gradS_vert, 2*gL->nvert);
    }
    /* clear */
    for (i = 0; i < 2*gL->ntria; i++)
	gradS_tri[i] = -1e30;

    BOTTOM_FACE *fb = gL->face_bot;
    for (ii = 0; ii < gL->nface_bot; ii++, fb++) {
	TRIA *t = gL->trias + fb->tria;
	FLOAT X[3][3];
	int nv, *iL;
	for (k = 0; k < 3; k++) {
	    i = t->verts[k];
	    X[k][0] = gL->verts[i][0];
	    X[k][1] = gL->verts[i][1];

	    j = gL->vert_S2L[i];
	    assert(gL->vert_local_lists[j] != NULL);
	    nv = gL->vert_local_lists[j][0];
	    iL = &gL->vert_local_lists[j][1];
	    X[k][2] = g->verts[iL[nv-1]][Z_DIR];
	}
	/* printf(" %4d %4d %4d\n", */
	/*        t->verts[0], t->verts[1], t->verts[2]); */
	/* printf("%12.5e %12.5e %12.5e\n%12.5e %12.5e %12.5e\n%12.5e %12.5e %12.5e\n", */
	/*        X[0][0],  X[1][0],  X[2][0], */
	/*        X[0][1],  X[1][1], X[2][1], */
	/*        X[0][2],  X[1][2], X[2][2]); */

	FLOAT A[2][2] = {
	    {X[0][0] - X[2][0], X[0][1] - X[2][1]},
	    {X[1][0] - X[2][0], X[1][1] - X[2][1]}
	};
	FLOAT b[2] = {X[0][2] - X[2][2],
		      X[1][2] - X[2][2]};
	FLOAT det = A[0][0] * A[1][1] - A[0][1] * A[1][0]; 
	FLOAT iA[2][2] = {
	    {A[1][1], -A[0][1]},
	    {-A[1][0], A[0][0]},
	};
	FLOAT gx = (iA[0][0] * b[0] + iA[0][1] * b[1]) / det;
	FLOAT gy = (iA[1][0] * b[0] + iA[1][1] * b[1]) / det;


	gradS_tri[2*fb->tria  ] = gx;
	gradS_tri[2*fb->tria+1] = gy;
	if (verb_core)
	    fprintf(fp, " %e %e %e %e %e \n",
		    1./3. * (X[0][0] + X[1][0] + X[2][0]),
		    1./3. * (X[0][1] + X[1][1] + X[2][1]),
		    1./3. * (X[0][2] + X[1][2] + X[2][2]),
		    gx, gy
		    );
    }
    if (verb_core)
	fclose(fp);


    /* Gather to root */
    MPI_Allreduce(gradS_tri, gradS_tri0, 2 * gL->ntria,
    		  PHG_MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    if (verb_core && phgRank == 0)
	for (ii = 0; ii < gL->ntria; ii++) {
	    printf("grads: %5d %e %e\n", ii, gradS_tri0[2*ii], gradS_tri0[2*ii+1]);
	}

#elif GRAD_SURF_SCHEME == 1
    /* ------------------------------------------------------------
     * 
     *
     * Use (S, grad phi)
     *
     *
     * ------------------------------------------------------------ */
    phgPrintf("Use (S, grad phi).\n");


    if (dofS != NULL)
	phgDofFree(&dofS);
    dofS = phgDofNew(g, DOF_P1, 1, "Surf", DofNoAction);

    BOTTOM_FACE *fb = gL->face_bot;
    FLOAT *vS = dofS->data;
    for (i = 0; i < g->nvert; i++) {
	if (g->types_vert[i] & BC_BOTTOM) {
	    assert(gL->vert_local_lists[i] != NULL);

	    FLOAT z;
	    int nv = gL->vert_local_lists[i][0];
	    int *iL = &gL->vert_local_lists[i][1];

	    z = g->verts[iL[nv-1]][Z_DIR];
	    for (j = 0; j < nv; j++)
		vS[iL[j]] = z;
	}
    }
    phgExportVTK(g, "surf.vtk", dofS, NULL);

#elif GRAD_SURF_SCHEME == 2
    /* ------------------------------------------------------------
     * 
     *
     * Use grad S, point wize
     *
     *
     * ------------------------------------------------------------ */
    phgPrintf("Use grad S(vert).\n");

    if (gradS_tri == NULL) {
	PHG_CALLOC(gradS_tri, 2*gL->ntria);
	PHG_CALLOC(gradS_tri0, 2*gL->ntria);
	PHG_CALLOC(gradS_vert, 2*gL->nvert);
    }
    /* clear */
    for (i = 0; i < 2*gL->ntria; i++)
	gradS_tri[i] = -1e30;

    BOTTOM_FACE *fb = gL->face_bot;
    for (ii = 0; ii < gL->nface_bot; ii++, fb++) {
	TRIA *t = gL->trias + fb->tria;
	FLOAT X[3][3];
	int nv, *iL;
	for (k = 0; k < 3; k++) {
	    i = t->verts[k];
	    X[k][0] = gL->verts[i][0];
	    X[k][1] = gL->verts[i][1];

	    j = gL->vert_S2L[i];
	    assert(gL->vert_local_lists[j] != NULL);
	    nv = gL->vert_local_lists[j][0];
	    iL = &gL->vert_local_lists[j][1];
	    X[k][2] = g->verts[iL[nv-1]][Z_DIR];
	}
	/* printf(" %4d %4d %4d\n", */
	/*        t->verts[0], t->verts[1], t->verts[2]); */
	/* printf("%12.5e %12.5e %12.5e\n%12.5e %12.5e %12.5e\n%12.5e %12.5e %12.5e\n", */
	/*        X[0][0],  X[1][0],  X[2][0], */
	/*        X[0][1],  X[1][1], X[2][1], */
	/*        X[0][2],  X[1][2], X[2][2]); */

	FLOAT A[2][2] = {
	    {X[0][0] - X[2][0], X[0][1] - X[2][1]},
	    {X[1][0] - X[2][0], X[1][1] - X[2][1]}
	};
	FLOAT b[2] = {X[0][2] - X[2][2],
		      X[1][2] - X[2][2]};
	FLOAT det = A[0][0] * A[1][1] - A[0][1] * A[1][0]; 
	FLOAT iA[2][2] = {
	    {A[1][1], -A[0][1]},
	    {-A[1][0], A[0][0]},
	};
	FLOAT gx = (iA[0][0] * b[0] + iA[0][1] * b[1]) / det;
	FLOAT gy = (iA[1][0] * b[0] + iA[1][1] * b[1]) / det;


	gradS_tri[2*fb->tria  ] = gx;
	gradS_tri[2*fb->tria+1] = gy;
	if (verb_core)
	    fprintf(fp, " %e %e %e %e %e \n",
		    1./3. * (X[0][0] + X[1][0] + X[2][0]),
		    1./3. * (X[0][1] + X[1][1] + X[2][1]),
		    1./3. * (X[0][2] + X[1][2] + X[2][2]),
		    gx, gy
		    );
    }
    if (verb_core)
	fclose(fp);


    /* Gather to root */
    MPI_Allreduce(gradS_tri, gradS_tri0, 2 * gL->ntria,
    		  PHG_MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    if (verb_core && phgRank == 0)
	for (ii = 0; ii < gL->ntria; ii++) {
	    printf("grads: %5d %e %e\n", ii, gradS_tri0[2*ii], gradS_tri0[2*ii+1]);
	}

    /* Projection */
    if (phgRank == 0) {
	static MAP *map = NULL;
	static MAT *mat = NULL;
	static SOLVER *solver = NULL;

	if (solver == NULL) {
	    phgPrintf("Create grad S proj solver\n");

	    /* Create solver */
	    map = phgMapCreateSimpleMap(MPI_COMM_SELF, gL->nvert, gL->nvert);
	    mat = phgMapCreateMat(map, map);
	
	    TRIA *t = gL->trias;
	    for (ii = 0; ii < gL->ntria; ii++, t++) {
		int I_[3] = {t->verts[0],
			    t->verts[1],
 			    t->verts[2]};
		FLOAT area = t->area;
		FLOAT A[3][3];
		for (i = 0; i < 3; i++)
		    for (j = 0; j < 3; j++)
			A[i][j] = (i == j) ? (area * 1./6.) : (area * 1./12.);

		for (i = 0; i < 3; i++)
		    phgMatAddEntries(mat, 1, I_+i, 3, I_, A[i]);
	    }
	    phgMatAssemble(mat);


	    phgOptionsPush();
	    phgOptionsSetOptions("-solver pcg -pcg_pc_type jacobi");
	    solver = phgMat2Solver(SOLVER_DEFAULT, mat);
	    phgOptionsPop();
	    phgMatDestroy(&mat);
	    phgSolverAssemble(solver);

	    if (0)
		phgMatDumpMATLAB(solver->mat, "M", "M_.m");
	}

	
	/* assemble rhs */
	VEC *rhsx = phgMapCreateVec(map, 1);
	VEC *rhsy = phgMapCreateVec(map, 1);
	VEC *vec_Sx = phgMapCreateVec(map, 1);
	VEC *vec_Sy = phgMapCreateVec(map, 1);
	FLOAT *rhs_dat = solver->rhs->data;
	phgVecDisassemble(rhsx);
	phgVecDisassemble(rhsy);

	TRIA *t = gL->trias;
	for (ii = 0; ii < gL->ntria; ii++, t++) {
	    int I_[3] = {t->verts[0],
			t->verts[1],
			t->verts[2]};
	    FLOAT area = t->area;
	    FLOAT bx[3], by[3];

	    for (i = 0; i < 3; i++) {
		bx[i] = area * 1./3. * gradS_tri0[2*ii];
		by[i] = area * 1./3. * gradS_tri0[2*ii+1];
	    }

	    phgVecAddEntries(rhsx, 0, 3, I_, bx);
	    phgVecAddEntries(rhsy, 0, 3, I_, by);
	}
	phgVecAssemble(rhsx);
	phgVecAssemble(rhsy);


	/* solve */
	solver->rhs->assembled = TRUE;
	solver->rhs->data = rhsx->data;
	phgPrintf("* Solver gradS[x].\n");
	phgSolverVecSolve(solver, FALSE, vec_Sx);
	memcpy(gradS_vert, vec_Sx->data, gL->nvert * sizeof(*gradS_vert));

	solver->rhs->assembled = TRUE;
	solver->rhs->data = rhsy->data;
	phgPrintf("* Solver gradS[y].\n");
	phgSolverVecSolve(solver, FALSE, vec_Sy);
	memcpy(gradS_vert + gL->nvert,
	       vec_Sy->data, gL->nvert * sizeof(*gradS_vert));

	/* Dump to check */
	if (verb_core) {
	    fp = fopen("grads1.dat", "w");
	    for (i = 0; i < gL->nvert; i++) 
		fprintf(fp, "%e %e %e %e %e\n", 
			gL->verts[i][0], 
			gL->verts[i][1], 
			0.,
			vec_Sx->data[i],
			vec_Sy->data[i]
			);
	    fclose(fp);
	}

	/* restorte */
	solver->rhs->data = rhs_dat;
	phgVecDestroy(&vec_Sx);
	phgVecDestroy(&vec_Sy);
    }
    
    /* Broadcast to other */
    MPI_Bcast(gradS_vert, 2 * gL->nvert,
	      PHG_MPI_FLOAT, 0, MPI_COMM_WORLD);
    
    /* interp to 3D */
    if (dof_gradS == NULL) {
	phgDofFree(&dof_gradS);
	dof_gradS = phgDofNew(g, DOF_P1, 2, "gradS", DofNoAction);
    }

    for (ii = 0; ii < gL->nvert_bot; ii++) {
	i = gL->vert_bot_Lidx[ii];
	assert(gL->vert_local_lists[i] != NULL);

	INT iG = gL->vert_bot_Gidx[ii];
	assert(iG < gL->nvert);

	int nv = gL->vert_local_lists[i][0];
	int *iL = &gL->vert_local_lists[i][1];
	    
	assert(nv > 0);
	for (j = 0; j < nv; j++) {
	    FLOAT *vg = DofVertexData(dof_gradS, iL[j]);
	    vg[0] = gradS_vert[iG];
	    vg[1] = gradS_vert[iG + gL->nvert];
	}
    }
    //phgExportVTK(g, "grad_surf.vtk", dof_gradS, NULL);

#endif



    /* -------------------------------------------------
     *
     * 
     * First order core
     *     non-linear iteration.
     *
     * -------------------------------------------------*/
    BTYPE DB_masks[3];
    dof_uv = phgDofNew(g, ns_params->utype, 2, "dof_uv", DofNoAction);
    phgPrintf("   DB_mask u: [");
    for (k = 0; k < 3; k++)
    	phgPrintf("%d ", ns->u[1]->DB_masks[k]);
    phgPrintf("]   ");

#if 0    
    memcpy(DB_masks, ns->u[1]->DB_masks, sizeof(DB_masks));
#else
    /* Shift DB masks */
    assert(UN_DIR == X_DIR);
    DB_masks[0] = ns->u[1]->DB_masks[1];
    DB_masks[1] = ns->u[1]->DB_masks[2];
    DB_masks[2] = ns->u[1]->DB_masks[0];
#endif    

/* #if USE_SLIDING_BC */
/*     Bzero(DB_masks);		/\* Sliding: all robin *\/ */
/* #endif     */

    phgDofSetDirichletBoundaryMasks(dof_uv, DB_masks);
    dof_du = phgDofCopy(dof_uv, NULL, NULL, "dof_du2");

    int max_nonstep = 0, newton_start = 0;
    //assert(ns_params->utype == DOF_P2);

    /* Copy solutions */
    copy_solution(ns->u[1], dof_uv, TRUE);
    grad_uv = phgDofGradient(dof_uv, NULL, NULL, "grad_uv");


    /* For nonlinear iter */
    int nonstep = 0; 
    FLOAT non_du = 1e+10, res;
    DOF *u_last = phgDofCopy(dof_uv, NULL, NULL, "u_last");
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
    CHECK_DBMASK;
    

    while (TRUE) {
	phgPrintf("\n   ==================\n");
	phgPrintf("   Non-linear interation step: %d\n", nonstep);
	    
	elapsed_time(g, FALSE, 0.);	/* reset timer */
	phgGetTime(tt);

	/* Init const viscosity */
	if (ns_params->start_const_vis &&
	    tstep == 1 && nonstep == 0) {
	    phgPrintf("* vis: const\n");
	    ns->viscosity_type = VIS_CONST;
	} else {		
	    phgPrintf("* vis: strain\n");
	    ns->viscosity_type = VIS_STRAIN;		
	}
    CHECK_DBMASK;
	    
	if (ns_params->reduce_mesh && gL != NULL)
	    mark_inactive(ns);

	sayHello("Non linear solve begin");
	{
	    phgOptionsPush();
	    phgOptionsSetOptions(ns_params->fo_opts);
	    
	    solver_uv = phgSolverCreate(SOLVER_DEFAULT, dof_uv, NULL);
	    
	    phgOptionsPop();
	}
    CHECK_DBMASK;

	if (nonstep < newton_start)
	    ns->ltype = PICARD;
	else 
	    ns->ltype = NEWTON;

	phgPrintf("   Build RHS: ");
	build_rhs(ns);
	elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));

	ns->non_res = 
	    res = phgVecNorm2(solver_uv->rhs, 0, NULL);
	if (0)
	    plot_residual(ns, solver_uv->rhs, nonstep);
	phgPrintf("   nonlinear residual: %24.12E\n", res);
    CHECK_DBMASK;

	/* Restore Picard if no improvement */
	if (ltype_last == NEWTON &&
	    res > non_res_last * .75) {
	    phgPrintf("   !!! Newton step failed, use Picard to run again\n");
	    ns->ltype = PICARD;
	    max_nonstep += 5; /* Add more Picard steps  */

	    /* resotre dofs:
	     * Fix me: temprature */
	    phgDofCopy(u_last, &dof_uv, NULL, "dof_uv");
	    phgDofGradient(dof_uv, &grad_uv, NULL, "grad_uv");
	    build_rhs(ns);
	    ns->non_res = 
		res = phgVecNorm2(solver_uv->rhs, 0, NULL);
	    phgPrintf("   nonlinear residual: %24.12E\n", res);
	} 
    CHECK_DBMASK;

	/* save non res */
	non_res_last = res;
	ltype_last = ns->ltype;
	    
	/* build matrices */
	/* if (ns_params->use_PCD) */
	/*     phgNSInitPc(ns); */
	phgPrintf("   Build matrices:\n");
	build_matrix(ns);
	phgPrintf("      done ");
	elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));
    CHECK_DBMASK;


	/*
	 * solve equation and update (u, p)
	 * */
	phgPrintf("solver tol: %E\n", solver_uv->rtol);
	phgSolverSolve(solver_uv, FALSE, dof_du, NULL);
#if USE_SLIDING_BC
	rotate_dof_bases(dof_du, surf_bas, FALSE);
#endif
	CHECK_DBMASK;
	phgPrintf("      solver_uv: nits = %d, resid = %0.4lg ",
		  solver_uv->nits, solver_uv->residual);
	elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));

	/* save dofs */
	phgDofCopy(dof_uv, &u_last, NULL, "u_last");
	CHECK_DBMASK;

	/* nonlinear correction */
	phgDofAXPY(1.0, dof_du, &dof_uv);
	CHECK_DBMASK;
	//phgDofDump(dof_du);

	assert(dof_uv->type == ns_params->utype);
	//assert(p[1]->type == ns_params->ptype);
#if USE_SLIDING_BC
	dof_set_normal_data(dof_uv, surf_bas);
#endif
	PERIODIC_SYNC(dof_uv);
	CHECK_DBMASK;


	/* non_du = phgDofNormL2(dof_uv); */
	non_du = phgDofNormInftyVec(dof_du);

	phgPrintf("   du: %24.12E \n", non_du);
	phgPrintf("   u: [%24.12E, %24.12E]\n", 
		  phgDofMinValVec(dof_uv), 
		  phgDofMaxValVec(dof_uv));
	phgDofGradient(dof_uv, &grad_uv, NULL, "gradu_{n+1}");

	phgSolverDestroy(&solver_uv);
	CHECK_DBMASK;
	

	/* evolution of u */
	//DOF_SCALE(u[1], "after solve");
	//DOF_SCALE(p[1], "after solve");


    
	if (ns_params->output_non_iter
	    && nonstep % ns_params->step_span == 0) {
	    phgPrintf("   Output solution to ensight ");
	    sprintf(vtk_file, OUTPUT_DIR "non_%02d_u.vtk", nonstep);
	    phgExportVTK(g, vtk_file , 
			 dof_uv, dof_du, NULL);
	    elapsed_time(g, TRUE, 0.);
	    //ice_monitor(ns, nonstep);
	}
	CHECK_DBMASK;


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
	if (//(res < ns_params->non_tol) 
	    nonstep >= ns_params->min_nonstep 
	    && ((ns->viscosity_type != VIS_CONST && 
		 non_du < ns_params->non_tol * U0
		 //&& non_dp * PRES_SCALING < ns_params->non_tol * P0
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

    phgDofFree(&u_last);



    /* Copy solutions */
    copy_solution(ns->u[1], dof_uv, FALSE);

    /* Reconstruct w and pressure */
    reconstruct_velocity_w(ns);

    dof_range(ns->u[1]);
    dof_range(ns->p[1]);


    /* temp save */
    phgPrintf("Save Dofs\n");
    save_dof_data3(g, ns->u[1], OUTPUT_DIR"u.dat");
    DOF_SCALE(ns->u[1], "save");

    
    /* free Dofs */
    phgDofFree(&dof_uv);
    phgDofFree(&dof_du);
    phgDofFree(&dof_w);
    phgDofFree(&grad_uv);

    return;
}
















void
reconstruct_velocity_w(NSSolver *ns)
/*
 * Reconstruct velocity w, by grad u, v
 * Only P1
 *
 * */
{
    GRID *g = ns->g;
    LAYERED_MESH *gL = ns->gL;
    BYTE components[DDim] = {1, 1, 1,
			     1, 1, 1,
			     0, 0, 0};
    int i, j, k; 
    assert(ns->u[1]->type == DOF_P1); /* ONLY P1 */

    
    phgPrintf("Reconstruct velocity w by incompressibility.\n");
    

    phgDofGradient(ns->u[1], &ns->gradu[1], NULL, "gradu_{n+1}");
    proj_gradu(ns, ns->gradu[1], components, 0);
    QUAD *quad = phgQuadGetQuad1D(4);

    

#if USE_SLIDING_BC    
    SURF_BAS *surf_bas = ns->surf_bas; /* 3D rotation */
    DOF *surf_dof = surf_bas->dof;
    BOOLEAN *rotated = surf_bas->rotated;
    const FLOAT *Trans = DofData(surf_dof);
    FLOAT *vU = DofData(ns->u[1]);
#endif    


    
    for (i = 0; i < g->nvert; i++) {
	if (g->types_vert[i] & BC_BOTTOM) {
	    assert(gL->vert_local_lists[i] != NULL);

	    FLOAT z0, z1, za, zb, ztop;
	    FLOAT gu[DDim], x, y, z, *vu;
	    int nv = gL->vert_local_lists[i][0];
	    int *iL = &gL->vert_local_lists[i][1];
	    FLOAT intz, len;
	    SIMPLEX *e, **elist;
	    int q, l;
	    const FLOAT *w, *p, *lambda;

	    elist = gL->vert_elements[i];
	    x = g->verts[iL[0]][X_DIR];
	    y = g->verts[iL[0]][Y_DIR];
	    intz = 0;

	    /*
	     * Non penetrate on bottom:
	     * 
	     * u * nx + v * ny + w * nz = 0
	     *
	     * */
	    if (ns_params->use_slide) {
		const FLOAT *trans = Trans + i*(Dim*Dim) + Dim*UN_DIR;
		intz = - (trans[X_DIR] * vU[i*Dim] + trans[Y_DIR] * vU[i*Dim+1]) / trans[Z_DIR];
		vU[i*Dim + Z_DIR] = intz;
	    }
	    
	    /* Pressure:
	     *  p = - DS_xx - DS_yy + rho g ( s - z )
	     * */
	    ztop = g->verts[iL[nv-1]][Z_DIR];
	    for (j = 0; j < nv; j++) {
		const FLOAT rho = RHO_ICE;
		const FLOAT grav = GRAVITY;
		const FLOAT *Gu = DofVertexData(ns->Gradu, iL[j]);
		const FLOAT *vT = DofVertexData(ns->T[1], iL[j]);

		FLOAT nu = get_effective_viscosity(Gu, *vT, 0, ns->viscosity_type);
		
		*DofVertexData(ns->p[1], iL[j]) =
		    (- nu * (Gu[0] + Gu[4]) / LEN_SCALING 
		       + rho * grav * (ztop - g->verts[iL[j]][Z_DIR]) * LEN_SCALING) / PRES_SCALING;
		/* *DofVertexData(dofp_1, iL[j]) = */
		/*     (+ nu * (Gu[0] + Gu[4]) / LEN_SCALING  */
		/*        + rho * grav * (ztop - g->verts[iL[j]][Z_DIR]) * LEN_SCALING) / PRES_SCALING; */
		/* *DofVertexData(dofp_2, iL[j]) = */
		/*     (+ 0 */
		/*        + rho * grav * (ztop - g->verts[iL[j]][Z_DIR]) * LEN_SCALING) / PRES_SCALING; */
		/* *DofVertexData(dofp_3, iL[j]) = */
		/*     (- nu * (Gu[0] + Gu[4]) / LEN_SCALING  */
		/*        + rho * grav * (ztop - g->verts[iL[j]][Z_DIR]) * LEN_SCALING) / PRES_SCALING; */
	    }
	    
	    for (j = 0; j < nv-1; j++) {
		z0 = g->verts[iL[j]  ][Z_DIR];
		z1 = g->verts[iL[j+1]][Z_DIR];
		e = elist[j];

		int edge;
		for (edge = 0; edge < NEdge; edge++) {
		    INT v0 = e->verts[GetEdgeVertex(edge, 0)];
		    INT v1 = e->verts[GetEdgeVertex(edge, 1)];
		    if ((v0 == iL[j] && v1 == iL[j+1])
			|| (v0 == iL[j+1] && v1 == iL[j])
			)
			break;
		}
		assert(edge < NEdge);


		/* segment */
		for (k = 0; k < 2; k++) {

		    if (k == 0) {
			za = z0;
			zb = .5 * (z0 + z1);
		    } else {
			za = .5 * (z0 + z1);
			zb = z1;
		    }
		    len = zb - za;

		    p = quad->points;
		    w = quad->weights;
		    for (q = 0; q < quad->npoints; q++) {
			z = za * p[0] + zb * p[1];
			lambda = phgGeomXYZ2Lambda(g, e, x, y, z);
			
			for (l = 0; l < 4; l++)
			    assert(-1e-10 < lambda[l] && lambda[l] < 1 + 1e-10);

			phgDofEval(ns->Gradu, e, lambda, gu);

			intz += (*w) * (len) * (-gu[0] - gu[4]);
			p += 2;
			w++;
		    } /* end quad pts */

		    /* Evaluation of horizontal velocity, part I */
		    if (k == 0) {
			if (ns->u[1]->type == DOF_P2) {
			    vu = DofEdgeData(ns->u[1], e->edges[edge]);
			    vu[Z_DIR] = intz;
			}
		    } else {
			vu = DofVertexData(ns->u[1], iL[j+1]);
			vu[Z_DIR] = intz;
		    } /* end int */

		}     /* end sub segment */
	    }	      /* end segment */
	}
    }

    //phgDofDump(ns->u[1]); 
    /* phgDofDump(ns->p[1]); */
    /* phgDofDump(dofp_1); */
    /* phgDofDump(dofp_2); */
    /* phgDofDump(dofp_3); */


    BOTTOM_FACE *fb = gL->face_bot;
    //e = NULL;
    for (i = 0; i < gL->nface_bot; i++, fb++) {
	TRIA *t = gL->trias + fb->tria;

	FLOAT intz[3] = {0, 0, 0}, len;
	int e2v_tri[3][2] = {
	    {0, 1}, {1, 2}, {2, 0}
	};
	INT vert_tri[3] = {t->verts[0],
			   t->verts[1],
			   t->verts[2]};
	FLOAT *x_tri[3] = {gL->verts[vert_tri[0]], 
			   gL->verts[vert_tri[1]], 
			   gL->verts[vert_tri[2]]};
	FLOAT mid_tet[NEdge][Dim];
	FLOAT mid_x[3], mid_y[3], z_top[3],
	    x, y, gu[DDim];	/* triangle edge coord */
	const FLOAT *w, *p, *lambda;
	int ii, l, q;
	FLOAT *vu;

	/* edge mid-pt */
	for (k = 0; k < 3; k++) {
	    int v0 = e2v_tri[k][0];
	    int v1 = e2v_tri[k][1];
	    mid_x[k] = .5 * (x_tri[v0][0] + x_tri[v1][0]); /* x */
	    mid_y[k] = .5 * (x_tri[v0][1] + x_tri[v1][1]); /* y */
	    z_top[k] = .5 * (x_tri[v0][3] + x_tri[v1][3]); /* z top */
	}
	
	for (ii = 0; ii < fb->ne; ii++) { /* tets */
	    SIMPLEX *e = fb->elems[ii];
	    int *edge2to3 = fb->edge2to3 + ii*6; /* [3][2] */
	    FLOAT za, zb, zz, z;
	    
	    /* tetra edge */
	    for (j = 0; j < NEdge; j++) {
		int v0 = GetEdgeVertex(j, 0); /* 3D */
		int v1 = GetEdgeVertex(j, 1);
		const FLOAT *x0 = g->verts[e->verts[v0]];
		const FLOAT *x1 = g->verts[e->verts[v1]];
		
		mid_tet[j][0] = .5 * (x0[0] + x1[0]);
		mid_tet[j][1] = .5 * (x0[1] + x1[1]);
		mid_tet[j][2] = .5 * (x0[2] + x1[2]);
	    }
	    
	    for (k = 0; k < 3; k++) { /* tri edge */
		int e0 = edge2to3[k*2];
		int e1 = edge2to3[k*2+1];
		int etmp;
		if (e1 == -1)
		    continue;
		
		x = mid_x[k];
		y = mid_y[k];
		
		za = mid_tet[e0][2];
		zb = mid_tet[e1][2];
		if (za > zb) {	/* za < zb, and edge e0 below edge e1 */
		    zz = za; za = zb; zb = zz;
		    etmp = e0; e0 = e1; e1 = etmp;
		}
		len = zb - za;
		if (fabs(zb - za) < 1e-8) {
		    printf("elem: %d\n", e->index);
		    printf("edge: %d %d\n", e0, e1);
		    printf("vert: %d %d %d %d\n", e->verts[0], e->verts[1], e->verts[2], e->verts[3]);
		    printf("vert: %15.7e %15.7e %15.7e\n", g->verts[e->verts[0]][0], g->verts[e->verts[0]][1], g->verts[e->verts[0]][2]);
		    printf("vert: %15.7e %15.7e %15.7e\n", g->verts[e->verts[1]][0], g->verts[e->verts[1]][1], g->verts[e->verts[1]][2]);
		    printf("vert: %15.7e %15.7e %15.7e\n", g->verts[e->verts[2]][0], g->verts[e->verts[2]][1], g->verts[e->verts[2]][2]);
		    printf("vert: %15.7e %15.7e %15.7e\n", g->verts[e->verts[3]][0], g->verts[e->verts[3]][1], g->verts[e->verts[3]][2]);
		    abort();
		}

		p = quad->points;
		w = quad->weights; 
		for (q = 0; q < quad->npoints; q++) {
		    z = za * p[0] + zb * p[1];
		    lambda = phgGeomXYZ2Lambda(g, e, x, y, z);
		    
		    for (l = 0; l < 4; l++)
			assert(-1e-10 < lambda[l] && lambda[l] < 1 + 1e-10);

		    phgDofEval(ns->Gradu, e, lambda, gu);

		    intz[k] += (*w) * (len) * (-gu[0] - gu[4]);

		    p += 2;
		    w++;
		} /* end quad pts */

		/* Evaluation of horizontal velocity, part II */
		if (ns->u[1]->type == DOF_P2) {
		    vu = DofEdgeData(ns->u[1], e->edges[e1]);
		    vu[Z_DIR] = intz[k];
		}


	    }	  /* end tri edge */
	}	  /* end tet */
    }		  /* end bot face */

}



