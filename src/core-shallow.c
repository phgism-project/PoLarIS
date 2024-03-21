/*
 *
 * Shallow ice approximation.
 *
 * * velocity formulation
 *    u = - 2 (rho g)^n |grad(s)|^(n-1) grad(s) \int_{b}^{z} A(T) (s-z)^n dz
 * * implementation
 *    //1. compute u element wise, including grad(s) and int from bottom
 *    //2. L2 projection to continuous u, e.g. P2 finite element 
 *   1. compute grads L2 projection
 *   2. compute int_{b}^{z} A(T) (s-z)^2 dz
 * 
 *  */
#include "ins.h"

#define CONST_A 0

/*
 *
 *  Note: SIA change to be done for z_b != 0, 
 *    1. gL->verts[3] = zbot and gL->verts[4] = ztop;
 *    2. Where height is updated, the 2D verts coord need to be updated
 *    3. grad S ==> grad height ?
 *
 *  */


int verb_core = 0;
#define CHECK_INTU 0


/*	 Paterson & budd (1982) */
#  define	 ARRHENIUS_T	263.15
#  define	 ARRHENIUS_A0   3.61E-13     
#  define	 ARRHENIUS_A1	1.73E3
#  define	 ARRHENIUS_Q0	6.0E4
#  define	 ARRHENIUS_Q1	13.9E4


static FLOAT *gradS_tri = NULL, *gradS_tri0 = NULL;
static FLOAT *gradS_vert = NULL;



void 
core_SIA(NSSolver *ns)
{
    GRID *g = ns->g;
    LAYERED_MESH *gL = ns->gL;
    SIMPLEX **elist, *e;
    int i, ii, j, k, l, q;
    FLOAT x, y, z, A, vT, vgS[2], norm_gS, *vu;
    const FLOAT *w, *p, *lambda;
    TRIA *t;
    FILE *fp;


    if (ns->u_edge == NULL) 
	ns->u_edge = phgDofNew(g, DOF_DG2, Dim, "u_edge", DofNoAction);


    /* ------------------------------------------------------------
     *
     * 1. Compute grads L2 projection
     *
     * ------------------------------------------------------------ */
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


    /*
     *
     *
     * Projection gradS to vert
     *
     *
     * */
    if (phgRank == 0) {
#if 0
#  warning proj gradS by L2
	/*
	 *
	 * Projection gradS to vert by global L2 
	 *
	 *  */
	static MAP *map = NULL;
	static MAT *mat = NULL;
	static SOLVER *solver = NULL;

	if (solver == NULL) {
	    phgPrintf("Create grad S proj solver\n");

	    /* Create solver */
	    map = phgMapCreateSimpleMap(MPI_COMM_SELF, gL->nvert, gL->nvert);
	    mat = phgMapCreateMat(map, map);
	
	    t = gL->trias;
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
	    phgOptionsSetOptions("-solver pcg -pcg_pc_type jacobi -solver_rtol 1e-10");
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

	t = gL->trias;
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
#else
#  warning proj gradS by averaging
	/*
	 *
	 * Projection gradS to vert by averaging
	 *
	 *  */
	static FLOAT *scale = NULL;
	scale = phgCalloc(gL->nvert, sizeof(*scale));

	for (i = 0; i < gL->nvert; i++) {
	    scale[i] = 0.;
	    gradS_vert[2*i  ] = 0.;
	    gradS_vert[2*i+1] = 0.;
	}

	t = gL->trias;
	for (ii = 0; ii < gL->ntria; ii++, t++) {
	    int I_[3] = {t->verts[0],
			t->verts[1],
			t->verts[2]};
	    FLOAT area = t->area;
	    
	    for (i = 0; i < 3; i++) {
		phgInfo(3, "vert node %3d %d %4d\n", ii, i, I_[i]);
		scale[I_[i]] += 1.;
		gradS_vert[I_[i]   ] += gradS_tri0[2*ii];
		gradS_vert[I_[i] + gL->nvert] 
		    += gradS_tri0[2*ii+1];
	    }
	}	

	for (i = 0; i < gL->nvert; i++) {
	    gradS_vert[i] /= scale[i];
	    gradS_vert[i + gL->nvert] /= scale[i];
	}

	phgFree(scale);
#endif
    }
    
    /* Broadcast to other */
    MPI_Bcast(gradS_vert, 2 * gL->nvert,
	      PHG_MPI_FLOAT, 0, MPI_COMM_WORLD);
    
    /* interp to 3D */
    if (ns->dof_gradS == NULL) 
	ns->dof_gradS = phgDofNew(g, DOF_P1, 2, "gradS", DofNoAction);
    DOF *dof_gradS = ns->dof_gradS;

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
 
    if (1)
	phgExportTecplot(g, "gradS.plt", dof_gradS, NULL);




    /*
     *
     *
     * Use gradS on edge, which is the average of neighs.
     * int the case of structure mesh, average of four points.
     *
     *
     * */
    static int nedge = 0;
    static int *elem_edges = NULL; /* nelem * 3 */
    static int *edge_elems = NULL; /* nedge * 2 */
    static int *elem_neigh = NULL; /* nelem * 3 */

    if (elem_edges == NULL) {
	EdgeVerts *edge2vert;
	INT count, size;

	count = 0;
	edge2vert = phgCalloc((size = 100), sizeof(EdgeVerts));

	t = gL->trias;
	for (ii = 0; ii < gL->ntria; ii++, t++) {
	    int ev[3][2] = {{0, 1}, {1, 2}, {2, 0}};

	    for (k = 0; k < 3; k++) {
		INT v0 = t->verts[ev[k][0]];
		INT v1 = t->verts[ev[k][1]];
		INT vv;

		SORT_TWO_VERT(v0, v1);

		if (count >= size) {
		    edge2vert = phgRealloc_(edge2vert,
					    2*size * sizeof(EdgeVerts),
					    size * sizeof(EdgeVerts));
		    size *= 2;
		}
		edge2vert[count][0] = v0;
		edge2vert[count][1] = v1;
		count++;
	    }
	}

	/* remove duplicate */
	qsort(edge2vert, count, sizeof(EdgeVerts), phgCompEdge);

	{
	    INT i0 = 0;
	    for (i = i0+1; i < count; i++) {
		int cmp = phgCompEdge(edge2vert + i0,
				      edge2vert + i);
		if (cmp < 0) {
		    i0++;
		    memcpy(edge2vert + i0,
			   edge2vert + i,
			   sizeof(EdgeVerts));
		}
	    }
	    count = i0 + 1;
	}
	nedge = count;


	elem_edges = phgCalloc(gL->ntria * 3, sizeof(*elem_edges));
	elem_neigh = phgCalloc(gL->ntria * 3, sizeof(*elem_edges));
	edge_elems = phgCalloc(nedge * 2, sizeof(*elem_edges));
	for (i = 0; i < nedge*2; i++) 
	    edge_elems[i] = -1;
	for (i = 0; i < gL->ntria*3; i++) 
	    elem_neigh[i] = -1;


	t = gL->trias;
	for (ii = 0; ii < gL->ntria; ii++, t++) {
	    int ev[3][2] = {{0, 1}, {1, 2}, {2, 0}};

	    for (k = 0; k < 3; k++) {
		INT v0 = t->verts[ev[k][0]];
		INT v1 = t->verts[ev[k][1]];
		INT vv;

		SORT_TWO_VERT(v0, v1);

		EdgeVerts vts, *p;
		vts[0] = v0;
		vts[1] = v1;
		p = bsearch(&vts, edge2vert,
			    count, sizeof(EdgeVerts),
			    phgCompEdge);
		assert(p != NULL);

		int edge_idx = p - edge2vert;
		elem_edges[3 * ii + k] = edge_idx;
		if (edge_elems[2 * edge_idx] == -1) {
		    edge_elems[2 * edge_idx] = ii;
		} else {
		    assert(edge_elems[2 * edge_idx + 1] == -1);
		    edge_elems[2 * edge_idx + 1] = ii;
		} 
	    }
	}

	for (i = 0; i < gL->ntria; i++) {
	    for (k = 0; k < 3; k++) {
		int edge_no = elem_edges[3*i + k];
		int e0 = edge_elems[2*edge_no + 0];
		int e1 = edge_elems[2*edge_no + 1];
		if (i == e0) {
		    elem_neigh[3*i + k] = e1; /* possible -1 */
		} else if (i == e1) {
		    elem_neigh[3*i + k] = e0;
		    assert(e0 != -1);
		} else {
		    abort();
		}
	    }

	    phgInfo(3, "Elem edges %d: %d %d %d\n", 
		    i, 
		    elem_edges[3 * i], 
		    elem_edges[3 * i+1], 
		    elem_edges[3 * i+2] 
		    );
	    phgInfo(3, "Elem neigh %d: %d %d %d\n", 
		    i, 
		    elem_neigh[3 * i], 
		    elem_neigh[3 * i+1], 
		    elem_neigh[3 * i+2] 
		    );
	}

	for (i = 0; i < nedge; i++) {
	    phgInfo(3, "Edge elems: %d %d\n", 
		    edge_elems[2 * i], 
		    edge_elems[2 * i+1] 
		    );
	}
    }







    /* ------------------------------------------------------------
     *
     * 2. compute int_{b}^{z} A(T) (s-z)^2 dz
     *
     * ------------------------------------------------------------ */
    const FLOAT n = POWER_N;
    const FLOAT a = SEC_PER_YEAR;
    const FLOAT T0 = ARRHENIUS_T;
    const FLOAT Q0 = ARRHENIUS_Q0;
    const FLOAT a0 = ARRHENIUS_A0;
    const FLOAT Q1 = ARRHENIUS_Q1;
    const FLOAT a1 = ARRHENIUS_A1;
    const FLOAT R  = GAS_CONSTANT;

    /*
     *
     * U = int_{b}^{z} A(T) (s-z)^2 dz
     *
     * Value on the vertex and veritcal edge mid-pt. 
     * Quad from triange vertices along z direction.
     * 
     * */
    QUAD *quad = phgQuadGetQuad1D(4);
    for (i = 0; i < g->nvert; i++) {
	if (g->types_vert[i] & BC_BOTTOM) {
	    assert(gL->vert_local_lists[i] != NULL);

	    FLOAT z0, z1, za, zb, z_top;
	    int nv = gL->vert_local_lists[i][0];
	    int *iL = &gL->vert_local_lists[i][1];
	    FLOAT intz, len;

	    elist = gL->vert_elements[i];
	    x = g->verts[iL[0]][X_DIR];
	    y = g->verts[iL[0]][Y_DIR];
	    z_top = g->verts[iL[nv-1]][Z_DIR];
	    intz = 0;

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

			phgDofEval(ns->T[1], e, lambda, &vT);
			if (q == 0)
			    phgDofEval(dof_gradS, e, lambda, vgS); /* same for all qpoints */
			//printf("q: %2d, vgS: %e %e %e %e\n", q, x, y, vgS[0], vgS[1]);

			if (vT < T0)
			    A = a0 * exp(-Q0 / R / vT);
			else
			    A = a1 * exp(-Q1 / R / vT);
			A *= SEC_PER_YEAR;
#if CONST_A
			A = 1e-16;
#endif

			intz += (*w) * (len * LEN_SCALING)
			    * A * pow((z_top - z) * LEN_SCALING, n);
			p += 2;
			w++;
		    } /* end quad pts */
		    norm_gS = sqrt(vgS[0]*vgS[0] + vgS[1]*vgS[1]);

		    /* Evaluation of horizontal velocity, part I */
		    if (k == 0) {
			if (ns->u[1]->type == DOF_P2) {
			    vu = DofEdgeData(ns->u[1], e->edges[edge]);

			    FLOAT tmp = (-2.) * pow(RHO_ICE * GRAVITY, n)
				* pow(norm_gS, n-1.) * intz;
			    vu[X_DIR] = vgS[X_DIR] * tmp;
			    vu[Y_DIR] = vgS[Y_DIR] * tmp;
			}
		    } else {
			vu = DofVertexData(ns->u[1], iL[j+1]);

			FLOAT tmp = (-2.) * pow(RHO_ICE * GRAVITY, n)
			    * pow(norm_gS, n-1.) * intz;
			vu[X_DIR] = tmp * vgS[X_DIR];
			vu[Y_DIR] = tmp * vgS[Y_DIR];
		    } /* end int */

		}     /* end sub segment */
	    }	      /* end segment */
	}
    }

    /*
     *
     * Value on the horizontal edge mid-pt.
     * Quad from triange edge mid-pt along z direction.
     * 
     * */
    fb = gL->face_bot;
    //e = NULL;
    for (i = 0; i < gL->nface_bot; i++, fb++) {
	TRIA *t = gL->trias + fb->tria;

	FLOAT intz[3] = {0, 0, 0}, len;
	FLOAT intz0[3] = {0, 0, 0};
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
	    x, y, tmp;	/* triangle edge coord */

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
	    FLOAT za, zb, zz;
	    
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


		/*
		 *
		 * gradS on edge: average
		 *
		 * */
		int t_neigh;
		if ((t_neigh = elem_neigh[3*fb->tria + k]) >= 0) {
		    vgS[X_DIR] = .5 * (gradS_tri0[2*t_neigh + X_DIR]
				       + gradS_tri0[2*fb->tria + X_DIR]
				       );
		    vgS[Y_DIR] = .5 * (gradS_tri0[2*t_neigh + Y_DIR]
				       + gradS_tri0[2*fb->tria + Y_DIR]
				       );
		} else {
		    vgS[X_DIR] = gradS_tri0[2*fb->tria + X_DIR];
		    vgS[Y_DIR] = gradS_tri0[2*fb->tria + Y_DIR];
		}

		if (ii == 0) {
		    phgInfo(3, "_gradS_ %e %e %e %e\n", 
			    mid_x[k], mid_y[k], vgS[0], vgS[1]);
		}

		
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

		    phgDofEval(ns->T[1], e, lambda, &vT);

#if 0
		    /* Use point wise grad S */
		    if (q == 0)
			phgDofEval(dof_gradS, e, lambda, vgS); /* same for all qpoints */
#endif

		    if (vT < T0)
			A = a0 * exp(-Q0 / R / vT);
		    else
			A = a1 * exp(-Q1 / R / vT);
		    A *= SEC_PER_YEAR;
#if CONST_A
		    A = 1e-16;
#endif

		    intz[k] += (*w) * (len * LEN_SCALING) * A
			* pow((z_top[k] - z) * LEN_SCALING, n);

		    p += 2;
		    w++;
		} /* end quad pts */
		norm_gS = sqrt(vgS[0]*vgS[0] + vgS[1]*vgS[1]);

		/* Evaluation of horizontal velocity, part II */
		if (ns->u[1]->type == DOF_P2) {
		    vu = DofEdgeData(ns->u[1], e->edges[e1]);

		    tmp = (-2.) * pow(RHO_ICE * GRAVITY, n)
			* pow(norm_gS, n-1.) * intz[k];
		    vu[X_DIR] = vgS[X_DIR] * tmp;
		    vu[Y_DIR] = vgS[Y_DIR] * tmp;
		}



#if 0
		/*
		 * Discontinouse edge velocity
		 * ??? why ???
		 *
		 * */
		vgS[0] = gradS_tri0[2*fb->tria + X_DIR];
		vgS[1] = gradS_tri0[2*fb->tria + Y_DIR];
		norm_gS = sqrt(vgS[0]*vgS[0] + vgS[1]*vgS[1]);

		/* upper edge */
		vu = DofElementData(ns->u_edge, e->index)
		    + (NVert + e1) * Dim;
		vu[X_DIR] = gradS_tri0[2*fb->tria + X_DIR] * tmp;
		vu[Y_DIR] = gradS_tri0[2*fb->tria + Y_DIR] * tmp;
		/* lower edge */
		tmp = (-2.) * pow(RHO_ICE * GRAVITY, n)
		    * pow(norm_gS, n-1.) * intz0[k];
		vu = DofElementData(ns->u_edge, e->index)
		    + (NVert + e0) * Dim;
		vu[X_DIR] = gradS_tri0[2*fb->tria + X_DIR] * tmp;
		vu[Y_DIR] = gradS_tri0[2*fb->tria + Y_DIR] * tmp;

		intz0[k] = intz[k];
#endif
	    }	  /* end tri edge */
	}	  /* end tet */

    }	      /* end bot face */




#if 0
    /* ------------------------------------------------------------
     *
     * 3. Compute w, by dw/dz = - du/dx - dv/dy
     *    solve (dw/dz, \phi) = (-du/dx - dv/dy, \phi)
     *  
     * ------------------------------------------------------------ */
    DOF *dof_w;
    SOLVER *solver_w;

    dof_w = phgDofNew(g, ns_params->utype, 1, "vel w", DofNoAction);
    phgDofSetDirichletBoundaryMask(dof_w, BC_BOTTOM);
    phgOptionsPush();
    phgOptionsSetOptions("-solver gmres "
			 "-solver_maxit 100 "
			 "-solver_rtol 1e-10");
    solver_w = phgSolverCreate(SOLVER_DEFAULT, dof_w,
			       NULL);
    solver_w->verb = 1;
    phgOptionsPop();
    
    phgDofGradient(ns->u[1], &ns->gradu[1], NULL, "gradu_{n+1}");

    ForAllElements(g, e) {
	int M = dof_w->type->nbas;	/* num of bases of Velocity */
	int order = DofTypeOrder(dof_w, e) * 2;
	FLOAT A[M][M], buf[M], rhs[M], values_u[M][Dim], ux, vy;
	INT I_[M];
	QUAD *quad;
	FLOAT vol, det;
	const FLOAT *w, *p, *gu;
	
	quad = phgQuadGetQuad3D(order);
	phgDofGetElementData(ns->u[1], e, values_u[0]);
	gu = phgQuadGetDofValues(e, ns->gradu[1], quad);  /* grad u */

	Bzero(A); Bzero(rhs); 

	p = quad->points;
	w = quad->weights;
	for (q = 0; q < quad->npoints; q++) {
	    phgGeomGetCurvedJacobianAtLambda(g, e, p, &det);
	    vol = fabs(det / 6.);

	    /* ux = vy = 0.; */
	    /* for (i = 0; i < M; i++) { */
	    /* 	const FLOAT *ggi = phgQuadGetBasisCurvedGradient(e, dof_w, i, quad, q);  */
	    /* 	ux += ggi[0] * values_u[i][0];  */
	    /* 	vy += ggi[1] * values_u[i][1];  */
	    /* } */
	
	    for (i = 0; i < M; i++) {
		const FLOAT *gi_u = phgQuadGetBasisValues(e, dof_w, i, quad) + q;    
		for (j = 0; j < M; j++) {
		    const FLOAT *ggj_u = phgQuadGetBasisCurvedGradient(e, dof_w, j, quad, q); 
		    
		    A[j][i] += vol*(*w) * (ggj_u[Z_DIR]) * (*gi_u);
		}

		//rhs[i] += vol*(*w) * (-ux - vy)  * (*gi_u);
		rhs[i] += vol*(*w) * (-gu[0] - gu[4])  * (*gi_u);
	    }
   
	    w++; p += Dim+1;
	    gu += DDim;
	} /* end quad pts */


	/* Map: Element -> system */
	for (i = 0; i < M; i++)
	    I_[i] = phgMapE2L(solver_w->mat->cmap, 0, e, i);

	for (i = 0; i < M; i++) {
	    if (phgDofDirichletBC_(dof_w, e, i, NULL, buf, &rhs[i],
				   DOF_PROJ_NONE)) {
		phgMatAddEntries(solver_w->mat, 1, I_ + i, M, I_, buf);
	    } else {
		phgMatAddEntries(solver_w->mat, 1, I_ + i, M, I_, &(A[i][0])); 
	    }
	}

	phgVecAddEntries(solver_w->rhs, 0, M, I_, &rhs[0]);
    }

    phgPrintf("* Solver vertical velocity.\n");
    phgSolverSolve(solver_w, TRUE, dof_w, NULL);

    phgExportVTK(g, "dofw.vtk", dof_w, NULL);

    /* copy to u^{n+1} */
    INT N = DofGetDataCount(dof_w);
    vu = ns->u[1]->data + Z_DIR;
    vw = dof_w->data;
    for (i = 0; i < N; i++,
	     vu += Dim, vw++) 
	*vu = *vw;
    phgDofFree(&dof_w);
#else
    /* ------------------------------------------------------------
     *
     * 3. Compute w, by w = \int_z (- du/dx - dv/dy)
     *    
     *  
     * ------------------------------------------------------------ */
    BYTE components[DDim] = {1, 0, 0,
			     0, 1, 0,
			     0, 0, 0};
    phgDofGradient(ns->u[1], &ns->gradu[1], NULL, "gradu_{n+1}");
    proj_gradu(ns, ns->gradu[1], components, 1); /* average */

    /* int ns->Gradu */
    for (i = 0; i < g->nvert; i++) {
	if (g->types_vert[i] & BC_BOTTOM) {
	    assert(gL->vert_local_lists[i] != NULL);

	    FLOAT z0, z1, za, zb, gu[DDim];
	    int nv = gL->vert_local_lists[i][0];
	    int *iL = &gL->vert_local_lists[i][1];
	    FLOAT intz, len;

	    elist = gL->vert_elements[i];
	    x = g->verts[iL[0]][X_DIR];
	    y = g->verts[iL[0]][Y_DIR];
	    intz = 0;

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

    fb = gL->face_bot;
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
	    FLOAT za, zb, zz;
	    
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



#endif	/* solve w */














    /*
     *
     * Check
     *
     *  */


#if CHECK_INTU 
    phgDofDump(ns->u[1]);
#endif /* CHECK_INTU */


    phgPrintf("   u: [%24.12E, %24.12E]\n", 
	      phgDofMinValVec(ns->u[1]), 
	      phgDofMaxValVec(ns->u[1]));
    phgPrintf("   p: [%24.12E, %24.12E]\n", 
	      phgDofMinValVec(ns->p[1]), 
	      phgDofMaxValVec(ns->p[1]));

    if (0)
	phgExportVTK(g, "test-SIA.vtk", ns->u[1], NULL);
    //phgDofFree(&dof_gradS);

    return;
}
























