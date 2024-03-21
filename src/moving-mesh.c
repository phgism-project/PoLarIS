#include "ins.h"
#if USE_PETSC
#  include <petscksp.h>
#endif


#define _nsp (ns->ns_params)


static void get_moved_coord_FEM(NSSolver *ns, int nstep);
static void get_moved_coord_FV_cell(NSSolver *ns, int nstep);
static void get_moved_coord_FV_point(NSSolver *ns, int nstep);
static void get_moved_coord_FEM_P2(NSSolver *ns, int tstep);


/*
 * Scheme:
 * 1. FEM implicit
 * 2. FEM explicit
 * 3. FV cell based
 * 4. FV point based
 *
 *  */
void 
get_moved_coord(NSSolver *ns, int tstep)
{
    /* if (_nsp->core_type == SIA)  */
    /* 	get_moved_coord_SIA(ns, tstep); */

    //assert(_nsp->height_scheme == 3);

    if (_nsp->height_scheme == 0) {
	phgPrintf("Height solver scheme: FEM implicit.\n");
	get_moved_coord_FEM(ns, tstep);
    } else if (_nsp->height_scheme == 1) {
	phgPrintf("Height solver scheme: FEM explicit.\n");
	get_moved_coord_FEM(ns, tstep);
    } else if (_nsp->height_scheme == 2) {
	phgPrintf("Height solver scheme: FV cell.\n");
	get_moved_coord_FV_cell(ns, tstep);
    } else if (_nsp->height_scheme == 3) {
	phgPrintf("Height solver scheme: FV point.\n");
	get_moved_coord_FV_point(ns, tstep);
    } else if (_nsp->height_scheme == 4) {
	phgPrintf("Height solver scheme: FEM P2.\n");
	get_moved_coord_FEM_P2(ns, tstep);
    } else {
	phgError(0, "Wrong scheme!\n");
    }
}



static void 
get_grad_bas(const FLOAT *J, const FLOAT *G,
		  FLOAT *p)
{
    *(p++) = G[0] * J[0 * (Dim + 1) + 0] +
	G[1] * J[1 * (Dim + 1) + 0] +
	G[2] * J[2 * (Dim + 1) + 0] +
	G[3] * J[3 * (Dim + 1) + 0];
    *(p++) = G[0] * J[0 * (Dim + 1) + 1] +
	G[1] * J[1 * (Dim + 1) + 1] +
	G[2] * J[2 * (Dim + 1) + 1] +
	G[3] * J[3 * (Dim + 1) + 1];
    *(p++) = G[0] * J[0 * (Dim + 1) + 2] +
	G[1] * J[1 * (Dim + 1) + 2] +
	G[2] * J[2 * (Dim + 1) + 2] +
	G[3] * J[3 * (Dim + 1) + 2];
}

static void
func_z(FLOAT x, FLOAT y, FLOAT z, FLOAT *v) 
{
    *v = z;
}



/* ================================================================================
 *
 *   update height using FEM 
 *   
 * ================================================================================
 */

void inv3(FLOAT a[3][3], FLOAT ia[3][3])
{
    int i;
    FLOAT det = a[0][0]*a[1][1]*a[2][2] - a[0][0]*a[1][2]*a[2][1] - a[0][1]*a[1][0]*a[2][2] 
	+ a[0][1]*a[1][2]*a[2][0] + a[0][2]*a[1][0]*a[2][1] - a[0][2]*a[1][1]*a[2][0];
    FLOAT b[3][3] = {
	{ a[1][1]*a[2][2] - a[1][2]*a[2][1], a[0][2]*a[2][1] - a[0][1]*a[2][2], a[0][1]*a[1][2] - a[0][2]*a[1][1]},
	{ a[1][2]*a[2][0] - a[1][0]*a[2][2], a[0][0]*a[2][2] - a[0][2]*a[2][0], a[0][2]*a[1][0] - a[0][0]*a[1][2]},
	{ a[1][0]*a[2][1] - a[1][1]*a[2][0], a[0][1]*a[2][0] - a[0][0]*a[2][1], a[0][0]*a[1][1] - a[0][1]*a[1][0]}
    };
    for (i = 0; i < 9; i++)
	ia[0][i] = b[0][i] / det;
}

/*
 *  
 *  Get Surface height change
 *  Output dH |_{bdry surf}
 *
 *  */
void
get_surf_dH(NSSolver *ns)
{
    GRID *g = ns->g;
    LAYERED_MESH *gL = ns->gL;
    DOF *u = ns->u[1];
    DOF *dH = ns->dH;
    FLOAT *dt = ns->dt;
    int height_scheme = _nsp->height_scheme;
    SOLVER *solver;
    SIMPLEX *e;
    INT i, j, ii, s;
    int verb;
    

    static FLOAT *surf_velo, *surf_velo0, *surf_elev, *surf_elev0;
    if (surf_velo == NULL) {
	PHG_CALLOC(surf_velo, 3*gL->nvert);
	PHG_CALLOC(surf_velo0, 3*gL->nvert);
	PHG_CALLOC(surf_elev, gL->nvert);
	PHG_CALLOC(surf_elev0, gL->nvert);
    }
    /* clear */
    for (i = 0; i < gL->nvert; i++) {
	surf_velo[3*i + 0] = -1e30;
	surf_velo[3*i + 1] = -1e30;
	surf_velo[3*i + 2] = -1e30;
	surf_elev[i] = -1e30;
    }

    /* Collect vertex velocity, elevation */
    for (ii = 0; ii < gL->nvert_bot; ii++) {
        i = gL->vert_bot_Lidx[ii];
        assert(gL->vert_local_lists[i] != NULL);

        INT iG = gL->vert_bot_Gidx[ii];
        assert(iG < gL->nvert);

        int nv = gL->vert_local_lists[i][0];
        int *iL = &gL->vert_local_lists[i][1];

        assert(nv > 0);
        //for (j = 0; j < nv; j++) {
	j = nv-1;
	{
            FLOAT *vu = DofVertexData(ns->u[1], iL[j]);
	    surf_velo[3*ii  ] = vu[0];
	    surf_velo[3*ii+1] = vu[1];
	    surf_velo[3*ii+2] = vu[2];
	    surf_elev[ii] = g->verts[iL[j]][Z_DIR];
        }
    }

    /* Gather to root */
    MPI_Allreduce(surf_velo, surf_velo0, 3 * gL->nvert,
    		  PHG_MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(surf_elev, surf_elev0, gL->nvert,
    		  PHG_MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    if (1 || phgRank == 0)
	for (ii = 0; ii < gL->nvert; ii++) {
	    printf("velocity: %5d %e %e %e %e\n", ii, 
		   surf_velo0[3*ii], 
		   surf_velo0[3*ii+1],
		   surf_velo0[3*ii+2], 
		   surf_elev0[ii]);
	}


    /* FEM 2D Solve */
    if (phgRank == 0) {
	MAP *map = phgMapCreateSimpleMap(MPI_COMM_SELF, gL->nvert, gL->nvert);
	MAT *mat = phgMapCreateMat(map, map);
	SOLVER *solver = NULL;
	VEC *rhs = phgMapCreateVec(map, 1);
	VEC *vec = phgMapCreateVec(map, 1);
	
	TRIA *t = gL->trias;
	for (ii = 0; ii < gL->ntria; ii++, t++) {
	    QUAD *quad;
	    int q;
	    const FLOAT *w, *p;
	    int I_[3] = {t->verts[0],
			t->verts[1],
			t->verts[2]};
	    INT vert_tri[3] = {t->verts[0],
			       t->verts[1],
			       t->verts[2]};
	    FLOAT *x[3] = {gL->verts[vert_tri[0]], 
			   gL->verts[vert_tri[1]], 
			   gL->verts[vert_tri[2]]};
	    FLOAT iJ[3][3] = {
		{x[0][0], x[1][0], x[2][0]},
		{x[0][1], x[1][1], x[2][1]},
		{1, 1, 1},
	    }, J[3][3];
	    FLOAT area = t->area;
	    FLOAT A[3][3], b[3];

	    bzero(A, sizeof(A));
	    bzero(b, sizeof(b));
	    inv3(iJ, J);

	    quad = phgQuadGetQuad2D(2);
	    p = quad->points;
	    w = quad->weights;
	    for (q = 0; q < quad->npoints; q++) {

		if (height_scheme == 1) {
		    /*
		     *
		     * Implicit
		     *
		     * */
		    FLOAT vu[3] = {0, 0, 0}, z = 0;
		    for (i = 0; i < 3; i++) {
			int k;
			FLOAT gi = p[i];
			FLOAT *ggi = J[i];
			for (k = 0; k < Dim; k++)
			    vu[k] += surf_velo0[3*I_[i] + k] * gi;
			z += surf_elev0[I_[i]] * gi;
		    }

		    for (i = 0; i < 3; i++) {
			FLOAT gi = p[i];
			FLOAT *ggi = J[i];
			for (j = 0; j < 3; j++) {
			    FLOAT gj = p[i];
			    FLOAT *ggj = J[j];

			    A[i][j] += area*(*w) * (gi*gj*LEN_SCALING / dt[0]
						    + (vu[0] * ggj[0] + vu[1] * ggj[1]) * gi
						    );
			    
			}

			b[i] += area*(*w) * (z * LEN_SCALING / dt[0]
					     + vu[Z_DIR]) * gi;
		    }
		} 
		else {
		    /*
		     *
		     * Explicit
		     *
		     * */
		    FLOAT vu[3] = {0, 0, 0}, dSx = 0, dSy = 0;
		    for (i = 0; i < 3; i++) {
			int k;
			FLOAT gi = p[i];
			FLOAT *ggi = J[i];
			for (k = 0; k < Dim; k++)
			    vu[k] += surf_velo0[3*I_[i] + k] * gi;
			dSx += surf_elev0[I_[i]] * ggi[0];
			dSy += surf_elev0[I_[i]] * ggi[1];
		    }

		    for (i = 0; i < 3; i++) {
			FLOAT gi = p[i];
			FLOAT *ggi = J[i];
			for (j = 0; j < 3; j++) {
			    FLOAT gj = p[i];
			    FLOAT *ggj = J[j];

			    A[i][j] += area*(*w) * (gi*gj);
			}

			FLOAT dh = vu[Z_DIR] 
			    - vu[Y_DIR] * (-dSx)
			    - vu[Z_DIR] * (-dSy);
			dh *= dt[0];
			b[i] += area*(*w) * (dh) * gi / LEN_SCALING;
		    }
		}
		
		p += 3;
		w++;
	    } /* end quad */

	    for (i = 0; i < 3; i++)
		phgMatAddEntries(mat, 1, I_+i, 3, I_, A[i]);
	    phgVecAddEntries(rhs, 0, 3, I_, b);
	}
	phgMatAssemble(mat);
	phgVecAssemble(rhs);


	phgOptionsPush();
	phgOptionsSetOptions("-solver pcg -pcg_pc_type jacobi");
	solver = phgMat2Solver(SOLVER_DEFAULT, mat);
	phgOptionsPop();
	phgMatDestroy(&mat);
	phgSolverAssemble(solver);

	if (0)
	    phgMatDumpMATLAB(solver->mat, "M", "M_.m");


	/* solve */
	solver->rhs->assembled = TRUE;
	phgVecAXPBY(1., rhs, 0, &solver->rhs);
	phgPrintf("* Solver 2D.\n");
	phgSolverVecSolve(solver, TRUE, vec);
	memcpy(surf_elev,
	       vec->data, gL->nvert * sizeof(*surf_elev));

	/* Dump to check */
	if (1) {
	    FILE *fp;
	    fp = fopen("ds.dat", "w");
	    for (i = 0; i < gL->nvert; i++) 
		fprintf(fp, "%e %e %e %e %e\n", 
			gL->verts[i][0], 
			gL->verts[i][1], 
			0.,
			vec->data[i],
			0.
			);
	    fclose(fp);
	}

	phgVecDestroy(&vec);
	phgMapDestroy(&map);
	exit(0);
    }
    
    /* Broadcast to other */
    MPI_Bcast(surf_elev, gL->nvert,
	      PHG_MPI_FLOAT, 0, MPI_COMM_WORLD);
    
    /* for (ii = 0; ii < gL->nvert_bot; ii++) { */
    /* 	i = gL->vert_bot_Lidx[ii]; */
    /* 	assert(gL->vert_local_lists[i] != NULL); */

    /* 	INT iG = gL->vert_bot_Gidx[ii]; */
    /* 	assert(iG < gL->nvert); */

    /* 	int nv = gL->vert_local_lists[i][0]; */
    /* 	int *iL = &gL->vert_local_lists[i][1]; */
	    
    /* 	assert(nv > 0); */
    /* 	for (j = 0; j < nv; j++) { */
    /* 	    FLOAT *vg = DofVertexData(dof_gradS, iL[j]); */
    /* 	    vg[0] = surf_velo[iG]; */
    /* 	    vg[1] = surf_velo[iG + gL->nvert]; */
    /* 	} */
    /* } */
    /* //phgExportVTK(g, "grad_surf.vtk", dof_gradS, NULL); */

    return;
}





/*
 * Get grid point z coord change
 * Output dH.
 * */
static void 
get_moved_coord_FEM(NSSolver *ns, int nstep)
{
    GRID *g = ns->g;
    DOF *u = ns->u[0];
    DOF *dH = ns->dH;
    LAYERED_MESH *gL = ns->gL;
    INT i, j;
    int ii, rank;
    FLOAT *sbuf, *rbuf, *dHS, *vH;
    FLOAT dH_max = -1e10, dH_min = 1e10;
    int verb;

    /* Use layer info */
    for (ii = 0; ii < gL->nvert_bot; ii++) {
	i = gL->vert_bot_Lidx[ii];
	assert(gL->vert_local_lists[i] != NULL);

	int nv = gL->vert_local_lists[i][0];
	int *iL = &gL->vert_local_lists[i][1];
	assert(nv > 0);
	
	FLOAT dh, dh0, dh1;
	dh0 = 0;
	dh1 = *DofVertexData(dH, iL[nv-1]);

	/* Set zero height to eps,
	 * dH is not used in this case. */
	if (TRUE) {
	    FLOAT *z = &g->verts[iL[nv-1]][Z_DIR];
	    if (*z + dh1  < HEIGHT_EPS) 
		dh1 = HEIGHT_EPS - *z; /* use H_eps */
	}

	dh = (dh1 - dh0) / nv;

	dH_max = MAX(dH_max, dh1);
	dH_min = MIN(dH_min, dh1);

	for (j = 0; j < nv; j++) {
	    *DofVertexData(dH, iL[j]) = j * dh;
	}
    }
    PERIODIC_SYNC(dH);	 
    //phgExportVTK(g, "output/dH.vtk", dH, NULL);

    {
	FLOAT a[2] = {dH_max, -dH_min}, b[2]; 
	MPI_Allreduce(&a, &b, 2, PHG_MPI_FLOAT, MPI_MAX, g->comm);
	dH_max = b[0];
	dH_min = -b[1];
	phgPrintf("   dH: [%8.4e, %8.4e]\n", dH_min, dH_max);
    }

    return;
}






/* ================================================================================
 *
 *   update height using FV (cell based)
 *      ONLY support NO FLOATING case !!!
 *   
 * ================================================================================
 */
static void 
get_moved_coord_FV_cell(NSSolver *ns, int nstep)
{
    GRID *g = ns->g;
    DOF *u = ns->u[1];
    LAYERED_MESH *gL = ns->gL;
    int i, j, k;
    DOF *dH = ns->dH;
    DOF *dHt = ns->dHt;
    BOTTOM_FACE *fb;
    FLOAT dt = ns_params->dt0;

    /* ------------------------------------------------------------
     * 
     * Step 1. compute flux
     *
     * ------------------------------------------------------------ */
    //FLOAT dh_vert[gL->nvert];
    FLOAT dh_gather[gL->ntria], dh_global[gL->ntria],
	dh_local[gL->nface_bot], max_dht = -1e10;
    DOF *dof_dh = phgDofNew(ns->g, DOF_P1, 1, "dof_dh", DofNoAction);

    if (0) {
	/* Removed */
	/* if (ns_params->use_prism_elem) { */
	/* 	get_cell_dh_prism(ns, dh_local); */
	/* } */
    }
    else {
	/*
	 *
	 * Use Tet  
	 *
	 *  */
	fb = gL->face_bot;
	//e = NULL;
	for (i = 0; i < gL->nface_bot; i++, fb++) {
	    FLOAT flux = 0.;
	    TRIA *t = gL->trias + fb->tria;
	    FLOAT area_tria = t->area;
	    int ii;
	
	    for (ii = 0; ii < fb->ne; ii++) {
		SIMPLEX *e = fb->elems[ii];
		INT *faces = fb->faces + ii*3;

		for (j = 0; j < 3; j++) { /* tria edges */
		    FLOAT *n = t->normals + j*2;
		    FLOAT nt[3] = {n[0], n[1], 0.};
		    FLOAT lambda[Dim+1];
		    int q, v0, v1, v2;
		    int s = faces[j];

		    if (s >= 0) {
			assert(s < NFace);
		    } else {
			continue;
		    }
		
		    QUAD *quad = phgQuadGetQuad2D(2);
		    FLOAT area = phgGeomGetFaceArea(g, e, s);
		    lambda[s] = 0;
		    v0 = GetFaceVertex(s, 0);
		    v1 = GetFaceVertex(s, 1);
		    v2 = GetFaceVertex(s, 2);
		
		    const FLOAT *p, *w;
		    p = quad->points;
		    w = quad->weights;
		    for (q = 0; q < quad->npoints; q++) {
			FLOAT vu[Dim];
			FLOAT x, y, z;
			lambda[v0] = *(p++);
			lambda[v1] = *(p++);
			lambda[v2] = *(p++);

			phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z);
			phgDofEval(u, e, lambda, vu);
			flux -= area*(*w) * INNER_PRODUCT(nt, vu);
			w++;
		    } /* end quad */
		}	  /* end tri edges */
	    }	  /* end elements */


	    dh_local[i] = flux / LEN_SCALING / area_tria; /* dh per year */
	    max_dht = MAX(fabs(dh_local[i]), max_dht);
	} /* end bot face */
    }




    /* ------------------------------------------------------------
     * 
     * Step 2.  gather h_bar
     *
     * ------------------------------------------------------------ */
    INT id_tria_global[gL->ntria],
	id_tria_local[gL->nface_bot];
    int rank, *tria_cnts, *tria_dsps, *tria_idxs, ntotal;

    tria_cnts = phgCalloc(2 * g->nprocs, sizeof(*tria_cnts));
    tria_dsps = tria_cnts + g->nprocs;

    MPI_Allgather(&gL->nface_bot, 1, MPI_INT, 
		  tria_cnts, 1, MPI_INT, g->comm);
    ntotal = 0;
    for (rank = 0; rank < g->nprocs; rank++) {
	tria_dsps[rank] = ntotal;
	ntotal += tria_cnts[rank];
    }

    fb = gL->face_bot;
    //e = NULL;
    for (i = 0; i < gL->nface_bot; i++, fb++) {
	id_tria_local[i] = fb->tria;
    }

    MPI_Allgatherv(id_tria_local,  gL->nface_bot, MPI_INT,
		   id_tria_global, tria_cnts, tria_dsps, 
		   MPI_INT, g->comm);
    MPI_Allgatherv(dh_local,    gL->nface_bot, PHG_MPI_FLOAT,
		   dh_gather,   tria_cnts, tria_dsps, 
		   PHG_MPI_FLOAT, g->comm);

    for (i = 0; i < gL->ntria; i++) {
	INT id_tria = id_tria_global[i];
	dh_global[id_tria] = dh_gather[i];
    }

    /* ------------------------------------------------------------
     * 
     * Step 3.  project from trigs to verts
     *
     * ------------------------------------------------------------ */
    /* \sum_T (dh_vert, \phi)_T = \sum_T (dh_tri, \phi)_T  */
    {
	SOLVER *solver;
	int ii, s, q, ib;
	
	phgOptionsPush();
	phgOptionsSetOptions(ns_params->proj_opts);
	solver = phgSolverCreate(SOLVER_DEFAULT, dof_dh, NULL);
	solver->verb = 2;
	phgOptionsPop();

	fb = gL->face_bot;
	for (ib = 0; ib < gL->nface_bot; ib++, fb++) {
	    TRIA *t = gL->trias + fb->tria;
	    FLOAT area_tri = t->area;
	
	    for (ii = 0; ii < fb->ne; ii++) { /* tets */
		SIMPLEX *e = fb->elems[ii];

		int M = dof_dh->type->nbas;	/* num of bases of Velocity */
		int order = DofTypeOrder(dof_dh, e) * 2;
		FLOAT A[M][M], rhs[M];
		INT I_[M];
		QUAD *quad;
		FLOAT vol, det, area;
		const FLOAT *w, *p, *gS;

	    
		Bzero(A); Bzero(rhs);
		/* ------------------------------
		 *
		 *   Face proj
		 * 
		 * ------------------------------ */
		for (s = 0; s < NFace; s++) {
		    if (e->bound_type[s] & BC_TOP) {
			int v0, v1, v2;
			int nbas_face = NbasFace(dof_dh);
			SHORT bases[nbas_face];
			FLOAT lambda[Dim + 1], x, y, z, gS[Dim];
			order = 2 * DofTypeOrder(dof_dh, e);

			phgDofGetBasesOnFace(dof_dh, e, s, bases);
			v0 = GetFaceVertex(s, 0);
			v1 = GetFaceVertex(s, 1);
			v2 = GetFaceVertex(s, 2);
			lambda[s] = 0.;


			/* Proj on manifold */
			// area = phgGeomGetFaceArea(g, e, s);
			area = area_tri;
			quad = phgQuadGetQuad2D(order);
		    
			p = quad->points;
			w = quad->weights;
			for (q = 0; q < quad->npoints; q++) {
			    lambda[v0] = *(p++);
			    lambda[v1] = *(p++);
			    lambda[v2] = *(p++);

			    const FLOAT *gi = 
				dof_dh->type->BasFuncs(dof_dh, e, 0, -1, lambda);
			
			    for (i = 0; i < nbas_face; i++) {
				int ii = bases[i];
				for (j = 0; j < nbas_face; j++) { 
				    int jj = bases[j];
				
				    FLOAT mass_face = area*(*w) * (gi[jj])*(gi[ii]);
				    A[ii][jj] += mass_face;
				}
		    
				rhs[ii] += area*(*w) * dh_global[fb->tria] * (gi[ii]); 
			    }

			    w++;
			} /* end qp */

			/* 		    SHOW_M(&A[0][0], M, M); */
			/* 		    SHOW_M(&rhs[0][0], 2, M); */
		    }     /* end face type*/
		}	  /* end face */

		for (i = 0; i < M; i++) {
		    BTYPE btype =
			phgDofGetElementBasisInfo(dof_dh, e, i,
						  NULL, NULL, NULL);
		    if (btype & BC_TOP)
			continue;

		    A[i][i] += 1;
		}

		/* Map: Element -> system */
		for (i = 0; i < M; i++)
		    I_[i] = phgMapE2L(solver->mat->cmap, 0, e, i);

		/* Global res */
		for (i = 0; i < M; i++)
		    phgMatAddEntries(solver->mat, 1, I_ + i, M, I_,
				     &(A[i][0])); 

		phgVecAddEntries(solver->rhs, 0, M, I_, rhs);
		
	    } /* end elems */
	}     /* end trias */


	phgSolverSolve(solver, FALSE, dof_dh, NULL);
	phgPrintf("      solver_dh: nits = %d, resid = %0.4lg ",
		  solver->nits, solver->residual);

	phgSolverDestroy(&solver);
    }

    /* ------------------------------------------------------------
     * 
     * Step 4.  update height
     *
     * ------------------------------------------------------------ */
    int ii;
    const FLOAT *ratio = gL->layer_ratio;
    FLOAT dH_max = -1e20, dH_min = 1e20;
    FLOAT *dh_vert = dof_dh->data;
    for (ii = 0; ii < gL->nvert_bot; ii++) {
	i = gL->vert_bot_Lidx[ii];
	assert(gL->vert_local_lists[i] != NULL);

	int nv = gL->vert_local_lists[i][0];
	int *iL = &gL->vert_local_lists[i][1];
	int iv = gL->vert_L2S[i];
	assert(nv > 0);
	
	FLOAT h0, h1;
	gL->height[iv] += dh_vert[iL[nv-1]] * dt;

# warning set zero height to eps
	/* Set zero height to eps,
	 * dH is not used in this case. */
	if (TRUE) {
	    if (gL->height[iv] < HEIGHT_EPS) 
		gL->height[iv] = HEIGHT_EPS;
	}

	h0 = gL->topg[iv];
	h1 = gL->height[iv] + h0;

	FLOAT H[nv];
	get_layer_height(H, nv, ratio, h0, h1);

	assert(gL->max_nlayer + 1 == nv);
	for (j = 0; j < nv; j++) {
	    g->verts[iL[j]][Z_DIR] = H[j];
	    *DofVertexData(ns->dHt, iL[j]) = dh_vert[iL[nv-1]];
	}
    }
    
    /* { */
    /*     FLOAT a[2] = {dH_max, -dH_min}, b[2]; */
    /*     MPI_Allreduce(&a, &b, 2, PHG_MPI_FLOAT, MPI_MAX, g->comm); */
    /*     dH_max = b[0]; */
    /*     dH_min = -b[1]; */
    /*     phgPrintf("   dH: [%8.4e, %8.4e]\n", dH_min, dH_max); */
    /* } */

    /* phgDofCopy(dof_dh, &ns->dHt, NULL, NULL); */
    
    phgDofFree(&dof_dh);
    phgFree(tria_cnts);

#warning FIXME: gL->height ???    
}

/* static void  */
/* get_moved_coord_FV_cell(NSSolver *ns, int nstep) */
/* { */
/*     GRID *g = ns->g; */
/*     DOF *u = ns->u[1]; */
/*     LAYERED_MESH *gL = ns->gL; */
/*     int i, j, k; */
/*     DOF *dH = ns->dH; */
/*     BOTTOM_FACE *fb; */
/*     FLOAT dt = ns_params->dt0; */

/*     /\* ------------------------------------------------------------ */
/*      *  */
/*      * Step 1. compute flux */
/*      * */
/*      * ------------------------------------------------------------ *\/ */
/*     FLOAT h_bar[gL->ntria],  */
/* 	h_bar_local[gL->nface_bot], */
/* 	h_bar_global[gL->ntria], */
/* 	height[gL->nvert]; */
/*     FLOAT dh_local[gL->nface_bot], max_dht = -1e10; */
    
/*     fb = gL->face_bot; */
/*     //e = NULL; */
/*     for (i = 0; i < gL->nface_bot; i++, fb++) { */
/* 	FLOAT flux = 0.; */
/* 	TRIA *t = gL->trias + fb->tria; */
/* 	FLOAT area_tria = t->area; */
/* 	int ii; */
	
/* 	for (ii = 0; ii < fb->ne; ii++) { */
/* 	    SIMPLEX *e = fb->elems[ii]; */
/* 	    INT *faces = fb->faces + ii*3; */

/* 	    for (j = 0; j < 3; j++) { /\* tria edges *\/ */
/* 		FLOAT *n = t->normals + j*2; */
/* 		FLOAT nt[3] = {n[0], n[1], 0.}; */
/* 		FLOAT lambda[Dim+1]; */
/* 		int q, v0, v1, v2; */
/* 		int s = faces[j]; */

/* 		if (s >= 0) { */
/* 		    assert(s < NFace); */
/* 		} else { */
/* 		    continue; */
/* 		} */
		
/* 		QUAD *quad = phgQuadGetQuad2D(1); */
/* 		FLOAT area = phgGeomGetFaceArea(g, e, s); */
/* 		lambda[s] = 0; */
/* 		v0 = GetFaceVertex(s, 0); */
/* 		v1 = GetFaceVertex(s, 1); */
/* 		v2 = GetFaceVertex(s, 2); */
		
/* 		const FLOAT *p, *w; */
/* 		p = quad->points; */
/* 		w = quad->weights; */
/* 		for (q = 0; q < quad->npoints; q++) { */
/* 		    FLOAT vu[Dim]; */
/* 		    FLOAT x, y, z; */
/* 		    lambda[v0] = *(p++); */
/* 		    lambda[v1] = *(p++); */
/* 		    lambda[v2] = *(p++); */

/* 		    phgGeomLambda2XYZ(g, e, lambda, &x, &y, &z); */
/* 		    phgDofEval(u, e, lambda, vu); */
/* 		    flux += area*(*w) * INNER_PRODUCT(nt, vu); */
/* 		    w++; */
/* 		} /\* end quad *\/ */
/* 	    }	  /\* end tri edges *\/ */
/* 	}	  /\* end elements *\/ */


/* 	dh_local[i] = flux / LEN_SCALING / area_tria; /\* dh per year *\/ */
/* 	max_dht = MAX(fabs(dh_local[i]), max_dht); */
/*     }		  /\* end bot face *\/ */

/* #   warning CFL limit */
/*     /\* */
/*      *  */
/*      * CFL limit: dh < 5e-3 */
/*      * */
/*      * *\/ */
/*     FLOAT dht = max_dht; */
/*     MPI_Allreduce(&dht, &max_dht, 1, PHG_MPI_FLOAT, MPI_MAX, g->comm); */
/*     while (dht * dt > 5e-3) */
/* 	dt  /= 2.; */
/*     phgPrintf("  time step change to %f\n", dt); */
/*     ns->dt[0] = dt; */

/*     fb = gL->face_bot; */
/*     for (i = 0; i < gL->nface_bot; i++, fb++) { */
/* 	FLOAT flux = 0.; */
/* 	TRIA *t = gL->trias + fb->tria; */
/* 	FLOAT area_tria = t->area; */
/* 	int ii; */
	
/* 	FLOAT vol = gL->volumes[fb->tria] - dh_local[i] * dt * area_tria; */
/* 	phgInfo(0, " volume[%3d]: %e, flux: %e\n", */
/* 		fb->tria, */
/* 		gL->volumes[fb->tria], */
/* 		flux * dt / LEN_SCALING); */

/* 	if (vol < HEIGHT_EPS * area_tria) { */
/* 	    vol = HEIGHT_EPS * area_tria; /\* Fix height *\/ */
/* 	    phgInfo(0, "*  fix height[%3d] e: %e\n", fb->tria, vol / area_tria); */
/* 	} */
/* 	gL->volumes[fb->tria] = vol; /\* only local tria */
/* 				      * volume *\/ */
/* 	h_bar_local[i] = vol / area_tria;  */
/* 	phgInfo(0, "  h_bar_local[%3d]: %e\n", i, h_bar_local[i]); */
/*     } */
/*     //printf("%3d %e\n", i, flux); */





/*     /\* ------------------------------------------------------------ */
/*      *  */
/*      * Step 2.  gather h_bar */
/*      * */
/*      * ------------------------------------------------------------ *\/ */
/*     INT   id_tria_global[gL->ntria], */
/* 	id_tria_local[gL->nface_bot]; */
/*     int rank, *tria_cnts, *tria_dsps, *tria_idxs, ntotal; */

/*     tria_cnts = phgCalloc(2 * g->nprocs, sizeof(*tria_cnts)); */
/*     tria_dsps = tria_cnts + g->nprocs; */

/*     MPI_Allgather(&gL->nface_bot, 1, MPI_INT,  */
/* 		  tria_cnts, 1, MPI_INT, g->comm); */
/*     ntotal = 0; */
/*     for (rank = 0; rank < g->nprocs; rank++) { */
/* 	tria_dsps[rank] = ntotal; */
/* 	ntotal += tria_cnts[rank]; */
/*     } */

/*     fb = gL->face_bot; */
/*     //e = NULL; */
/*     for (i = 0; i < gL->nface_bot; i++, fb++) { */
/* 	id_tria_local[i] = fb->tria; */
/*     } */

/*     MPI_Allgatherv(id_tria_local,  gL->nface_bot, MPI_INT, */
/* 		   id_tria_global, tria_cnts, tria_dsps,  */
/* 		   MPI_INT, g->comm); */
/*     MPI_Allgatherv(h_bar_local,    gL->nface_bot, PHG_MPI_FLOAT, */
/* 		   h_bar_global,   tria_cnts, tria_dsps,  */
/* 		   PHG_MPI_FLOAT, g->comm); */

/*     for (i = 0; i < gL->ntria; i++) { */
/* 	INT id_tria = id_tria_global[i]; */
/* 	h_bar[id_tria] = h_bar_global[i]; */
/*     } */

/*     /\* ------------------------------------------------------------ */
/*      *  */
/*      * Step 3.  project from trigs to verts */
/*      * */
/*      * ------------------------------------------------------------ *\/ */
    
/* #if USE_PETSC */
/*     /\* \sum_T (h, \phi)_T = \sum_T (h^bar, \phi)_T  *\/ */
/*     { */
/* 	Vec x, b, u; */
/* 	Mat A; */
/* 	KSP ksp; */
/* 	PC pc; */
/* 	PetscReal norm; */
/* 	PetscErrorCode ierr; */
/* 	PetscInt i, j, k, n = gL->nvert, col[3], its; */
/* 	PetscMPIInt size; */
/* 	PetscScalar one = 1.0, value[3][3], r[3]; */
/* 	int rank = 0; */

/* 	//PetscInitialize(&argc, &args, (char *) 0, help); */
/* 	MPI_Group orig_group, new_group;  */
/* 	MPI_Comm new_comm;  */
/* 	//MPI_Comm_rank(MPI_COMM_WORLD, &rank);  */
/* 	MPI_Comm_group(MPI_COMM_WORLD, &orig_group);  */
/* 	MPI_Group_incl(orig_group, 1, &rank, &new_group); */
/* 	MPI_Comm_create(MPI_COMM_WORLD, new_group, &new_comm);  */

/* 	if (g->rank == 0) { */

/* 	    MPI_Comm_size(new_comm, &size); */
/* 	    if (size != 1) */
/* #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 1 */
/* 		SETERRQ(1, "This is a uniprocessor example only!"); */
/* #else */
/* 		SETERRQ(MPI_COMM_WORLD, 1, "This is a uniprocessor example only!"); */
/* #endif */

/* 	    /\* Create vectors. *\/ */
/* 	    VecCreate(new_comm, &x);	 */
/* 	    PetscObjectSetName((PetscObject) x, "Solution");	 */
/* 	    VecSetSizes(x, PETSC_DECIDE, n);	 */
/* 	    VecSetFromOptions(x);	 */
/* 	    VecDuplicate(x, &b);	 */
/* 	    VecDuplicate(x, &u);	 */

/* 	    /\* Create matrix. *\/ */
/* 	    MatCreate(new_comm, &A);	 */
/* 	    MatSetType(A,MATSEQAIJ);      */
/* 	    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n);	 */
/* 	    //MatSetFromOptions(A);	 */

/* 	    /\* Assemble matrix *\/ */
/* 	    TRIA *t = gL->trias; */
/* 	    for (k = 0; k < gL->ntria; k++, t++) { */
/* 	    	/\* nodes of tria *\/ */
/* 	    	PetscInt I_[3] = {t->verts[0], */
/* 	    			 t->verts[1], */
/* 	    			 t->verts[2]}; */
/* 		FLOAT *coord[3]; */
/* 	    	PetscScalar area = t->area; */

/* 		for (i = 0; i < 3; i++) { */
/* 		    coord[i] = gL->verts[I_[i]]; */
/* 		} */
/* 		FLOAT x, y; */
/* 		x = (coord[0][0]+ coord[1][0] + coord[2][0]) / 3.; */
/* 		y = (coord[0][1]+ coord[1][1] + coord[2][1]) / 3.; */

/* 		/\* New height */
/* 		 * Vialov */
/* 		 * *\/ */
/* 		if (0) { */
/* 		    x -= 750.; */
/* 		    y -= 750.; */
/* 		    FLOAT d = sqrt(x*x + y*y); */
/* 		    FLOAT h0 = 300.; */
/* 		    FLOAT L = 500.; */
/* 		    FLOAT n = POWER_N; */

/* 		    if (d < L) */
/* 			h_bar[k] = h0 * pow(fabs(1 - pow(d/L, (n+1)/n)),  */
/* 					    n/(2*n+2)); */
/* 		    else */
/* 			h_bar[k] = HEIGHT_EPS; */
/* 		} */

/* 	    	for (i = 0; i < 3; i++) { */
/* 	    	    for (j = 0; j < 3; j++) */
/* 	    		value[i][j] = area * ((i == j) ? (1./6.) : (1./12.)); */
/* 	    	    r[i] = area * h_bar[k] * (1./3.); */
/* 	    	} */

/* 	    	for (i = 0; i < 3; i++) { */
/* 	    	    MatSetValues(A, 1, I_+i, 3, I_, */
/* 	    				value[i], ADD_VALUES);   */
/* 	    	} */

/* 	    	VecSetValues(b, 3, I_, r, ADD_VALUES);   */
/* 	    } */

/* 	    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);	 */
/* 	    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);	 */
/* 	    VecAssemblyBegin(b); */
/* 	    VecAssemblyEnd(b); */

/* 	    /\* PetscViewer viewer; *\/ */
/* 	    /\* PetscViewerASCIIOpen(PETSC_COMM_WORLD,"A_.m",&viewer); *\/ */
/* 	    /\* PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); *\/ */
/* 	    /\* MatView(A, viewer);  *\/ */
/* 	    /\* PetscViewerASCIIOpen(PETSC_COMM_WORLD,"b_.m",&viewer); *\/ */
/* 	    /\* PetscViewerSetFormat(viewer, PETSC_VIEWER_ASCII_MATLAB); *\/ */
/* 	    /\* VecView(b, viewer); *\/ */

/* 	    VecSet(u, 0);	  */

/* 	    /\* KSP *\/ */
/* 	    KSPCreate(new_comm, &ksp);	 */
/* #if PETSC_VERSION_MAJOR<3 || (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<5) */
/* 	    KSPSetOperators(ksp, A, A, DIFFERENT_NONZERO_PATTERN);	 */
/* #else	/\* PETSc < 3.5 *\/ */
/* 	    KSPSetOperators(ksp, A, A); */
/* #endif	/\* PETSc < 3.5 *\/ */
/* 	    KSPGetPC(ksp, &pc);	 */
/* 	    PCSetType(pc, PCJACOBI);	 */
/* 	    ierr = */
/* 		KSPSetTolerances(ksp, 1.e-12, PETSC_DEFAULT, PETSC_DEFAULT, */
/* 				 PETSC_DEFAULT);	 */
/* 	    KSPSetFromOptions(ksp);	 */
/* 	    KSPSolve(ksp, b, x);	 */

/* 	    PetscScalar *xx; */
/* 	    phgInfo(0, "New height\n"); */
/* 	    VecGetArray(x, &xx);  */
/* 	    for (i = 0; i < gL->nvert; i++) { */
/* 		height[i] = xx[i]; */
/* 		if (height[i] < HEIGHT_EPS) {/\* Fix height *\/ */
/* 		    phgInfo(0, "*   fix height[%3d] v %e\n", i, height[i]); */
/* 		    height[i] = HEIGHT_EPS; */
/* 		} */
/* 		phgInfo(0, "  %3d %e\n", i, height[i]); */
/* 	    } */

/* 	    /\* Free *\/ */
/* #if PETSC_VERSION_MAJOR == 3 && PETSC_VERSION_MINOR == 1 */
/* 	    VecDestroy(x);	 */
/* 	    VecDestroy(u);	 */
/* 	    VecDestroy(b);	 */
/* 	    MatDestroy(A);	 */
/* 	    KSPDestroy(ksp);	 */
/* #else */
/* 	    VecDestroy(&x);	 */
/* 	    VecDestroy(&u);	 */
/* 	    VecDestroy(&b);	 */
/* 	    MatDestroy(&A);	 */
/* 	    KSPDestroy(&ksp);	 */
/* #endif */
/* 	} */
/*     } */
/* #else */
/*     phgError(0, "Petsc solver disabled!!!\n"); */
/* #endif */

/*     /\* ------------------------------------------------------------ */
/*      *  */
/*      * Step 4.  update height */
/*      * */
/*      * ------------------------------------------------------------ *\/ */
/*     MPI_Bcast(height, gL->nvert, PHG_MPI_FLOAT, 0, MPI_COMM_WORLD);     */

/*     int ii; */
/*     FLOAT dH_max = -1e20, dH_min = 1e20; */
/*     for (ii = 0; ii < gL->nvert_bot; ii++) { */
/* 	i = gL->vert_bot_Lidx[ii]; */
/* 	assert(gL->vert_local_lists[i] != NULL); */

/* 	int nv = gL->vert_local_lists[i][0]; */
/* 	int *iL = &gL->vert_local_lists[i][1]; */
/* 	int iv = gL->vert_L2S[i]; */
/* 	assert(nv > 0); */
	
/* 	FLOAT dh, dh0, dh1; */
/* 	dh0 = 0; */
/* 	//dh1 = *DofVertexData(dH, iL[nv-1]); */
/* 	dh1 = height[iv] - g->verts[iL[nv-1]][Z_DIR]; */

/* # warning set zero height to eps */
/* 	/\* Set zero height to eps, */
/* 	 * dH is not used in this case. *\/ */
/* 	if (TRUE) { */
/* 	    FLOAT *z = &g->verts[iL[nv-1]][Z_DIR]; */
/* 	    if (*z + dh1  < HEIGHT_EPS)  */
/* 		dh1 = HEIGHT_EPS - *z; /\* use H_eps *\/ */
/* 	} */

/* 	dh = (dh1 - dh0) / nv; */

/* 	dH_max = MAX(dH_max, dh1); */
/* 	dH_min = MIN(dH_min, dh1); */

/* 	for (j = 0; j < nv; j++) { */
/* 	    *DofVertexData(dH, iL[j]) = j * dh; */
/* 	} */
/*     } */
    
/*     { */
/*         FLOAT a[2] = {dH_max, -dH_min}, b[2]; */
/*         MPI_Allreduce(&a, &b, 2, PHG_MPI_FLOAT, MPI_MAX, g->comm); */
/*         dH_max = b[0]; */
/*         dH_min = -b[1]; */
/*         phgPrintf("   dH: [%8.4e, %8.4e]\n", dH_min, dH_max); */
/*     } */
/*     phgFree(tria_cnts); */

/* } */








/* ================================================================================
 *
 *   update height using FV  (point based),
 *      ONLY support NO FLOATING case !!!
 *   
 * ================================================================================
 */
static void
get_moved_coord_FV_point(NSSolver *ns, int tstep)
{
    GRID *g = ns->g;
    DOF *u = ns->u[1];
    LAYERED_MESH *gL = ns->gL;
    int i, j, k;
    DOF *dH = ns->dH;
    DOF *dHt = ns->dHt;
    BOTTOM_FACE *fb;
    FLOAT dt = ns_params->dt0;
    INT nL = gL->max_nlayer;
    FLOAT *dht;

    dht = phgCalloc(gL->nvert, sizeof(*dht)); 
    

    /* ------------------------------------------------------------
     * 
     * Step 1. compute flux
     *
     * ------------------------------------------------------------ */
    FLOAT max_dht = -1e10;
    FLOAT *U_local, *U_vert, *U; 
    U       = phgCalloc(gL->ntria * 6, sizeof(*U)); /* [ne][3][2] */
    U_local = phgCalloc(gL->ntria * 6, sizeof(*U)); /* [ne][3][2] */
    U_vert  = phgCalloc(gL->nvert * 2, sizeof(*U)); /* [nv][2] */

    fb = gL->face_bot;
    //e = NULL;
    for (i = 0; i < gL->nface_bot; i++, fb++) {
	TRIA *t = gL->trias + fb->tria;
	FLOAT area_tri = t->area;
	int ii;

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
	FLOAT lens[3] = {0, 0, 0};
	FLOAT *u_local = U_local + fb->tria * 6; /* [3][2] */

	for (ii = 0; ii < fb->ne; ii++) { /* tets */
	    SIMPLEX *e = fb->elems[ii];
	    static FLOAT vu_edge[NEdge][Dim];
	    const FLOAT *vu[NEdge];
	    int *edge2to3 = fb->edge2to3 + ii*6; /* [3][2] */
	    
	    /* tetra edge */
	    for (j = 0; j < NEdge; j++) {
		if (ns_params->core_type == STOKES) {
		    if (ns->u[1]->type == DOF_P2) {
			vu[j] = DofEdgeData(ns->u[1], e->edges[j]);
		    } else if (ns->u[1]->type == DOF_P1){
			/* P1, edge middle */
			vu[j] = vu_edge[j];
			int v0 = e->verts[GetEdgeVertex(j, 0)]; 
			int v1 = e->verts[GetEdgeVertex(j, 1)]; 
			
			const FLOAT *u0 = &ns->u[1]->data[v0*Dim];
			const FLOAT *u1 = &ns->u[1]->data[v1*Dim];

			for (k = 0; k < Dim; k++)
			    vu_edge[j][k] = .5 * (u0[k] + u1[k]);
		    } else {
			abort();
		    }
		} else if (ns_params->core_type == FIRST_ORDER) {
		    if (ns->u[1]->type == DOF_P2) {
			vu[j] = DofEdgeData(ns->u[1], e->edges[j]);
		    } else if (ns->u[1]->type == DOF_P1){
			/* P1, edge middle */
			vu[j] = vu_edge[j];
			int v0 = e->verts[GetEdgeVertex(j, 0)]; 
			int v1 = e->verts[GetEdgeVertex(j, 1)]; 
			
			const FLOAT *u0 = &ns->u[1]->data[v0*Dim];
			const FLOAT *u1 = &ns->u[1]->data[v1*Dim];

			for (k = 0; k < Dim; k++)
			    vu_edge[j][k] = .5 * (u0[k] + u1[k]);

		    } else {
			abort();
		    }
		} 
		else if (ns_params->core_type == SIA || 
			 ns_params->core_type == DEBUG_CORE1) {
		    if (ns->u[1]->type == DOF_P2) {
			vu[j] = DofEdgeData(ns->u[1], e->edges[j]);
		    } else {
			/* P1, edge middle */
			vu[j] = vu_edge[j];
			int v0 = e->verts[GetEdgeVertex(j, 0)]; 
			int v1 = e->verts[GetEdgeVertex(j, 1)]; 
			
			const FLOAT *u0 = &ns->u[1]->data[v0*Dim];
			const FLOAT *u1 = &ns->u[1]->data[v1*Dim];

			for (k = 0; k < Dim; k++)
			    vu_edge[j][k] = .5 * (u0[k] + u1[k]);
		    }
		}
		else {
		    phgError(0, "Unknown core!!!\n");
		}
		
		int v0 = GetEdgeVertex(j, 0); /* 3D */
		int v1 = GetEdgeVertex(j, 1);
		FLOAT *x0 = g->verts[e->verts[v0]];
		FLOAT *x1 = g->verts[e->verts[v1]];
		
		mid_tet[j][0] = .5 * (x0[0] + x1[0]);
		mid_tet[j][1] = .5 * (x0[1] + x1[1]);
		mid_tet[j][2] = .5 * (x0[2] + x1[2]);
	    }


	    const double x0 = 6.515425e+02, y0 = 3.519559e+02;

	    for (k = 0; k < 3; k++) { /* tri edge */
		int e0 = edge2to3[k*2];
		int e1 = edge2to3[k*2+1];
		if (e1 == -1)
		    continue;
		
		/* flux v0 -> v1 */
		FLOAT len = fabs(mid_tet[e0][2] - mid_tet[e1][2]); /* z dir */

		FLOAT u = .5 * (vu[e0][0] + vu[e1][0]) * len; /* \int u dz */
		FLOAT v = .5 * (vu[e0][1] + vu[e1][1]) * len;

		u_local[2*k   ] += u / LEN_SCALING;
		u_local[2*k +1] += v / LEN_SCALING;  /* unit km/a */

		/* if (fabs(x_tri[k][0] - x0) < 1 */
		/*     && fabs(x_tri[k][1] - y0) < 1) { */
		/*     printf("vel %e %e, %e %e %e, %e %e %e\n", */
		/* 	   x_tri[k][0], x_tri[k][1], */
		/* 	   mid_tet[e0][2], vu[e0][0], vu[e0][1],  */
		/* 	   mid_tet[e1][2], vu[e1][0], vu[e1][1]  */
		/* 	   );  */
		/* } */

		/* if (g->types_edge[e->edges[e0]] & BC_TOP) */
		/*     printf("vel %e %e %e %e\n", x_tri[k][0], x_tri[k][1], vu[e0][0], vu[e0][1]); */
		/* if (g->types_edge[e->edges[e1]] & BC_TOP) */
		/*     printf("vel %e %e %e %e\n", x_tri[k][0], x_tri[k][1], vu[e1][0], vu[e1][1]); */
		    

		lens[k] += len;
	    } /* end tri edge */
	}     /* end tet */

	for (k = 0; k < 3; k++) { /* velocity mean */
	    u_local[2*k   ] /= lens[k];
	    u_local[2*k +1] /= lens[k];
	    /* printf("U   %e %e %e %e\n", */
	    /* 	   x_tri[k][0], x_tri[k][1], */
	    /* 	   u_local[2*k], u_local[2*k +1]); */
	}

	phgInfo(3, "U_local[%d * 6]: %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e\n", i,
		u_local[0], u_local[1], u_local[2],
		u_local[3], u_local[4], u_local[5]);
    }	      /* end bot face */


    /* Reduce to root
     * (actually non is sumed)
     * */
    MPI_Reduce(U_local, U, gL->ntria * 6,
	       PHG_MPI_FLOAT, MPI_SUM, 0, g->comm);
    for (i = 0; i < gL->ntria; i++) {
	FLOAT *u = U + i * 6;
	phgInfo(3, "U[%d * 6]: %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e\n", i, 
		u[0], u[1], u[2], u[3], u[4], u[5]);
    }

    /*
     * FV solve
     * */
    if (g->rank == 0) {
	FLOAT max_dht = 0;
	fv_update(gL->fv_data, gL->height, gL->verts,
		  U, dht, U_vert, 
		  ns->time[1]);

	/* CFL limit */
	for (i = 0; i < gL->nvert; i++)
	    if (max_dht < fabs(dht[i]))
		max_dht = fabs(dht[i]);
	phgPrintf("   Max dht %15.7e\n", max_dht);

#warning CFL condition 1e-2
	while (max_dht * dt > 1e-1)
	    dt  /= 2.;
	phgPrintf("   Time step change to %f\n", dt);
	ns->dt[0] = dt;

	for (i = 0; i < gL->nvert; i++) {
	    FLOAT h = gL->height[i] + dht[i] * dt;

	    if (h < HEIGHT_EPS) {
		h = HEIGHT_EPS; /* Fix height */
		phgInfo(3, "*  fix height[%3d] e: %15.7e\n", i, h);
	    }

	    gL->height[i] = h;
	}
	
	/* debug */
	if (0) {
	    char name[100];
	    sprintf(name, "./output/H_%05d_.m", tstep);
	    FILE *fp = fopen(name, "w");
	    //fprintf(fp, "H = [\n");
	    fprintf(fp, "#x, y, h, dht, u, v\n");
	    for (i = 0; i < gL->nvert; i++) {
		fprintf(fp, "%15.7e %15.7e %15.7e %15.7e %15.7e %15.7e\n",
			gL->verts[i][0], gL->verts[i][1],
			gL->height[i], dht[i],
			U_vert[2*i], U_vert[2*i+1]);
	    }
	    //fprintf(fp, "];\n");
	    fclose(fp);
	}
    }
    phgFree(U_local);
    phgFree(U);
    phgFree(U_vert);


    MPI_Bcast(ns->dt, 1,
	      PHG_MPI_FLOAT, 0, MPI_COMM_WORLD);    
    MPI_Bcast(gL->height, gL->nvert,
	      PHG_MPI_FLOAT, 0, MPI_COMM_WORLD);    
    MPI_Bcast(dht, gL->nvert,
	      PHG_MPI_FLOAT, 0, MPI_COMM_WORLD);    

    /* update gL->verts */
    for (i = 0; i < gL->nvert; i++) {
	gL->verts[i][3] = gL->height[i];
    }

    /* update height */
    int ii;
    const FLOAT *ratio = gL->layer_ratio;
    for (ii = 0; ii < gL->nvert_bot; ii++) {
	i = gL->vert_bot_Lidx[ii];
	assert(gL->vert_local_lists[i] != NULL);

	int nv = gL->vert_local_lists[i][0];
	int *iL = &gL->vert_local_lists[i][1];
	int iv = gL->vert_L2S[i];
	assert(nv > 0);

	FLOAT h0, h1;
	h0 = gL->topg[iv];
	h1 = gL->height[iv] + h0;

	FLOAT H[nv];
	get_layer_height(H, nv, ratio, h0, h1);
		
	assert(gL->max_nlayer + 1 == nv);
	for (j = 0; j < nv; j++) {
	    g->verts[iL[j]][Z_DIR] = H[j];
	    *DofVertexData(dHt, iL[j]) = dht[iv];
	}
    }

    phgFree(dht);
}



void
check_height(GRID *g, LAYERED_MESH *gL)
{
    INT n = gL->nvert, i;
    FLOAT sum_h = 0., de; 
    
    phgInfo(0, "%s not implemented.\n", __FUNCTION__);
    return;

    for (i = 0; i < n; i++){
	sum_h += gL->height[i] * gL->height[i];
    }
    sum_h = sqrt(sum_h);

    {
	FLOAT a[2] = {sum_h, -sum_h}, b[2];
	MPI_Allreduce(a, b, 2, 
		      PHG_MPI_FLOAT, MPI_MAX, g->comm);
	de = (fabs(b[0] - sum_h)
	      + fabs(b[1] + sum_h) ) / 2;
    }
    phgPrintf("   *height %15.7e +- %15.7e\n", sum_h, de);
    
}




/* Note: Move mesh coord is divied to two step,
 * 1. Compute dH(dZ)
 * 2. Z += dH
 * This is because Z is not periodic, but dH is.
 *
 * */
void
move_mesh(NSSolver *ns)
{
    GRID *g = ns->g;
    SIMPLEX *e;
    FLOAT *vc, *vz;
    LAYERED_MESH *gL = ns->gL;
    INT i;

    /* g->verts changed */
    phgGeomInit_(g, TRUE);

    /* Set inactive element masks */
    ForAllElements(g, e) {
	for (i = 0; i < NVert; i++)
	    if (g->verts[i][Z_DIR] > HEIGHT_EPS)
		break;
	if (i < NVert) {	/* at least one non-zero: active */
	    e->region_mark = 0;
	} else {		/* all zero: inactive */
	    e->region_mark = 1;
	}
    }


    //phgDofDump(dz);
    phgExportMedit(g, "moved_mesh.mesh");
    phgDofSetDataByFunction(ns->coord, func_xyz_);
    get_height_depth(ns);
    

    /* Volume */
    FLOAT Vol = 0;
    ForAllElements(g, e) {
	Vol += phgGeomGetVolume(g, e);
    }
    FLOAT tmp = Vol;
    MPI_Allreduce(&tmp, &Vol, 1, MPI_DOUBLE, MPI_SUM, g->comm);
    phgPrintf("   Total volume: %24.12e\n", Vol);

    return;
}    




/* ------------------------------------------------------------
 *
 * 
 * Use FEM P2,
 *    use 3D grid to solve
 * 
 *
 * ------------------------------------------------------------ */
#define DEBUG_2D_P2 0

#if DEBUG_2D_P2
void func_vel_test(FLOAT x, FLOAT y, FLOAT z, FLOAT *values)
{
    values[0] = 1.;
    values[1] = 0;
    values[2] = 0.1;
}
#endif

static void
get_moved_coord_FEM_P2(NSSolver *ns, int tstep)
{
    GRID *g = ns->g;
    DOF *u = ns->u[1];
    LAYERED_MESH *gL = ns->gL;
    SOLVER *solver;
    DOF *mapping = GetGeomMapping(g);
    DOF *dof_dht = ns->dHt;
    BOTTOM_FACE *fb;
    int i, j, k, ib, ii, jj, s, q, itet;

    /*
     * First, 
     *   Solve 2D problem on TOP surface
     *
     * */
    assert(mapping->type == DOF_P2);


#if DEBUG_2D_P2
    phgDofSetDataByFunction(ns->u[1], func_vel_test);
#endif    

    
    //dof_dht = phgDofNew(g, mapping->type, 1, "dH", DofNoAction);
    phgDofSetDataByValue(dof_dht, 0);

    phgOptionsPush();
    phgOptionsSetOptions(ns_params->proj_opts);
    solver = phgSolverCreate(SOLVER_DEFAULT, dof_dht, NULL);
    solver->verb = 2;
    phgOptionsPop();


    /*
     *  Build system
     *  */
    fb = gL->face_bot;
    for (ib = 0; ib < gL->nface_bot; ib++, fb++) {
	TRIA *t = gL->trias + fb->tria;
	FLOAT area_tri = t->area;
	
	for (itet = 0; itet < fb->ne; itet++) { /* tets */
	    SIMPLEX *e = fb->elems[itet];

	    int M = dof_dht->type->nbas;	/* num of bases of Velocity */
	    int order = DofTypeOrder(dof_dht, e) * 2;
	    FLOAT coord_elem_values[M][Dim];
	    FLOAT A[M][M], rhs[M];
	    INT I_[M];
	    QUAD *quad;
	    FLOAT vol, det, area;
	    const FLOAT *w, *p, *gS;

	    if (g->non_period) {
		SIMPLEX *ee = g->non_period->elems[e->index];
		assert(ee->index == e->index);
		phgDofGetElementDatas(mapping, ee, coord_elem_values[0]);
	    } else {
		phgDofGetElementDatas(mapping, e, coord_elem_values[0]);
	    }
	    
	    Bzero(A); Bzero(rhs);
	    /* ------------------------------
	     *
	     *   Face proj
	     * 
	     * ------------------------------ */
	    for (s = 0; s < NFace; s++) {
		if (e->bound_type[s] & BC_TOP) {
		    int v0, v1, v2;
		    int nbas_face = NbasFace(dof_dht);
		    SHORT bases[nbas_face];
		    FLOAT lambda[Dim + 1], vu[Dim], x, y, z, gS[Dim];
		    order = 2 * DofTypeOrder(dof_dht, e);

		    phgDofGetBasesOnFace(dof_dht, e, s, bases);
		    v0 = GetFaceVertex(s, 0);
		    v1 = GetFaceVertex(s, 1);
		    v2 = GetFaceVertex(s, 2);
		    lambda[s] = 0.;

		    /* 2D jacobian */
		    FLOAT *X[3] = {g->verts[e->verts[v0]],
				   g->verts[e->verts[v1]],
				   g->verts[e->verts[v2]]};
		    FLOAT iJ[3][3] = {
			{X[0][0], X[1][0], X[2][0]},
			{X[0][1], X[1][1], X[2][1]},
			{1, 1, 1},
		    }, J[3][3];
		    FLOAT gbas[M][Dim];
		    Bzero(gbas);

		    inv3(iJ, J);
		    
		    /* Proj on manifold */
		    // area = phgGeomGetFaceArea(g, e, s);
		    area = area_tri;
		    quad = phgQuadGetQuad2D(order);
		    
		    p = quad->points;
		    w = quad->weights;
		    for (q = 0; q < quad->npoints; q++) {
			FLOAT gS[Dim] = {0, 0, 0};
			lambda[v0] = *(p++);
			lambda[v1] = *(p++);
			lambda[v2] = *(p++);

			const FLOAT *vbas = 
			    dof_dht->type->BasFuncs(dof_dht, e, 0, -1, lambda);
			phgDofEval(ns->u[1], e, lambda, vu);
			phgInfo(0, "vel: %e %e %e\n", vu[0], vu[1], vu[2]);

			/* Gbas --> gbas */
			const FLOAT *Gbas = 
			    dof_dht->type->BasGrads(dof_dht, e, 0, -1, lambda);

			for (ii = 0; ii < nbas_face; ii++) {
			    int i = bases[ii];
			    const FLOAT *G = Gbas + i * (Dim + 1);
			    
			    /* compute data[Dim] = G[Dim] * J[Dim][Dim] */
			    gbas[i][0] = G[v0] * J[0][0]
				       + G[v1] * J[1][0]
				       + G[v2] * J[2][0];
			    gbas[i][1] = G[v0] * J[0][1]
				       + G[v1] * J[1][1]
				       + G[v2] * J[2][1];
			    gbas[i][2] = G[v0] * J[0][2]
				       + G[v1] * J[1][2]
				       + G[v2] * J[2][2];

			    //phgInfo(0, "   gbas[%d] %f %f %f\n", i, gbas[i][0], gbas[i][1], gbas[i][2]);
			    
			    gS[X_DIR] += coord_elem_values[i][Z_DIR] * gbas[i][X_DIR];
			    gS[Y_DIR] += coord_elem_values[i][Z_DIR] * gbas[i][Y_DIR];
			}

			//phgInfo(0, " gS %f %f\n", gS[X_DIR], gS[Y_DIR]);
			
			
			for (ii = 0; ii < nbas_face; ii++) {
			    int i = bases[ii];
			    for (jj = 0; jj < nbas_face; jj++) { 
				int j = bases[jj];
				
				A[i][j] += area*(*w)
				    * (
				       vbas[j] * LEN_SCALING
				       ) * vbas[i];
			    }
		    
			    rhs[i] += area*(*w) *
				(
				 vu[Z_DIR]
				 - vu[X_DIR] * gS[X_DIR]
				 - vu[Y_DIR] * gS[Y_DIR]
				 ) * vbas[i];
				
			}

			w++;
		    } /* end qp */

		    /* 		    SHOW_M(&A[0][0], M, M); */
		    /* 		    SHOW_M(&rhs[0][0], 2, M); */
		}     /* end face type*/
	    }	  /* end face */

	    for (i = 0; i < M; i++) {
		BTYPE btype =
		    phgDofGetElementBasisInfo(dof_dht, e, i,
					      NULL, NULL, NULL);
		if (btype & BC_TOP)
		    continue;

		A[i][i] += 1;
	    }

	    /* Map: Element -> system */
	    for (i = 0; i < M; i++)
		I_[i] = phgMapE2L(solver->mat->cmap, 0, e, i);

	    /* SHOW_iV_(2, I_, M); */
	    /* SHOW_M(A[0], M, M); */
	    
	    /* Global res */
	    for (i = 0; i < M; i++)
		phgMatAddEntries(solver->mat, 1, I_ + i, M, I_,
				 &(A[i][0])); 

	    phgVecAddEntries(solver->rhs, 0, M, I_, rhs);
		
	} /* end elems */
    }     /* end trias */


    phgSolverSolve(solver, TRUE, dof_dht, NULL);
    phgPrintf("      solver_dh: nits = %d, resid = %0.4lg ",
	      solver->nits, solver->residual);
    phgSolverDestroy(&solver);




    /* Copy dht from top to bottom verts */
    const FLOAT *ratio = gL->layer_ratio;
    for (ii = 0; ii < gL->nvert_bot; ii++) {
	i = gL->vert_bot_Lidx[ii];
	assert(gL->vert_local_lists[i] != NULL);

	int nv = gL->vert_local_lists[i][0];
	int *iL = &gL->vert_local_lists[i][1];
	int iv = gL->vert_L2S[i];
	assert(nv > 0);

	FLOAT dht = *DofVertexData(ns->dHt, iL[nv-1]);
	
	assert(gL->max_nlayer + 1 == nv);
	for (j = 0; j < nv; j++) {
	    *DofVertexData(ns->dHt, iL[j]) = dht;
	}
    }


    
    /* FIXME  TODO:
     *
     * ??? How to update mapping ???
     *
     * */

#if DEBUG_2D_P2
    phgDofDump(dof_dh);
    phgExportVTK(g, "dh.vtk", dof_dh, NULL);
    phgFinalize();
    exit(0);
#endif
    
}




/*
 * Build Dof Height and Depth,
 *   and change NO ice BC mask
 *  
 *  */
void
get_height_depth(NSSolver *ns)
{
    GRID *g = ns->g;
    LAYERED_MESH *gL = ns->gL;
    DOF *depth_P1, *depth_T, *height_P1, *dof_mask; 
    FLOAT z0, z1;
    int i, j, ii;


    phgPrintf("   Get height depth and CHANGE NO ICE BC MASK...\n");
    depth_P1 = phgDofNew(g, DOF_P1, 1, "depth P1", DofNoAction);
    depth_T = phgDofNew(g, ns_params->T_type, 1, "depth P2", DofNoAction);
    height_P1 = phgDofNew(g, DOF_P1, 1, "height P1", DofNoAction);
    dof_mask = phgDofNew(g, DOF_P1, 1, "mask", DofNoAction);
    

    for (ii = 0; ii < gL->nvert_bot; ii++) {
	i = gL->vert_bot_Lidx[ii];
	assert(gL->vert_local_lists[i] != NULL);

	int nv = gL->vert_local_lists[i][0];
	int *iL = &gL->vert_local_lists[i][1];
	assert(nv > 0);

	z0 = g->verts[iL[nv-1]][Z_DIR];
	z1 = g->verts[iL[0]][Z_DIR];

	for (j = 0; j < nv; j++) {
	    *DofVertexData(depth_P1, iL[j])
			= z0 - g->verts[iL[j]][Z_DIR];
	    *DofVertexData(height_P1, iL[j])
			= g->verts[iL[j]][Z_DIR] - z1;
	}

	double thk = z0 - z1;
	if (thk <= HEIGHT_EPS * (1 + 1e-6) ) {
	    /* NO ice: h <= eps */
	    g->types_vert[iL[0]] |= BC_INACTIVE; /* inactive */
	    *DofVertexData(dof_mask, iL[0]) = 1.;
	}
	else {
	    /* h > eps */
	    g->types_vert[iL[0]] &= ~BC_INACTIVE; /* active */
	}
    }

    /* interp */
    phgDofCopy(depth_P1, &depth_T, NULL, "depth P2");
    //phgExportVTK(g, "depth.vtk", depth_P1, depth_T, NULL);
    if (0) {
	phgExportVTK(g, "thk_mask.vtk", depth_P1, dof_mask, NULL);
	phgFinalize();
	exit(1);
    }

    if (ns->depth_P1 != NULL)
	phgDofFree(&ns->depth_P1);
    if (ns->depth_T != NULL)
	phgDofFree(&ns->depth_T);

    ns->depth_P1 = depth_P1;
    ns->depth_T = depth_T;
    ns->height = height_P1;

	/* Init sigma coord */
	if (ns->sigma_z == NULL) {
		int ii, j;
		ns->sigma_z = phgDofNew(g, DOF_P1, 1, "sigma", DofNoAction);

		for (ii = 0; ii < gL->nvert_bot; ii++) {
			int i = gL->vert_bot_Lidx[ii];
			int nv = gL->vert_local_lists[i][0];
			int *iL = &gL->vert_local_lists[i][1];

			for (j = 0; j < nv; j++) {
				*DofVertexData(ns->sigma_z, iL[j]) =
					gL->layer_ratio[j];
			}
		}
	}
	

    phgPrintf("   Get height depth Done.\n");
    phgDofFree(&dof_mask);
    return;
}


/*
 * Update Dof value when mesh coords change.
 * NOT used in the simulation, since the ice-sheet movement is small.
 *
 * Use algorithm in Y. Di, R. Li, T. Tang, and P. Zhang. Moving mesh
 * finite element methods for the incompressible Navier-Stokes
 * equations. SIAM J. Sci. Comput., 26(3):1036–1056, 2005
 *
 *  */
void 
move_dof(GRID *g, DOF *dz, DOF *u)
{
    SIMPLEX *e;
    MAP *map;
    MAT *mat;
    SOLVER *solver;
    int dim = u->dim;
    INT i, j, k;
    DOF *gradu;
    FLOAT msl = 1.;
    BTYPE DB_mask_save[100];
    int verb;

    gradu = phgDofGradient(u, NULL, NULL, NULL);


    /* Note: save DB_mask of dofs, and currently set it to UNDEFINED
     *       to remove the effect of DIRICHLET bdry nodes. */
    if (u->DB_masks != NULL) {
	memcpy(DB_mask_save, u->DB_masks, dim * sizeof(*DB_mask_save));
	memset(u->DB_masks, dim, sizeof(*DB_mask_save));
    } else
	DB_mask_save[0] = u->DB_mask;

    /* Create dof update solver */
    map = phgMapCreate(u, NULL);
    mat = phgMapCreateMat(map, map);
    phgOptionsPush();
    phgOptionsSetOptions("-solver petsc "
			 "-solver_maxit 10000 "
			 "-solver_rtol 1e-10");
    verb = phgVerbosity;
    phgVerbosity = 1;
    solver = phgMat2Solver(SOLVER_DEFAULT, mat);
    phgVerbosity = verb;
    phgVecDisassemble(solver->rhs);
    phgOptionsPop();

    ForAllElements(g, e) {
	int M = DofGetNBas(u, e);	/* number of basises in the element */
	int order = DofTypeOrder(u, e) * 2;
	FLOAT A[M][dim][M][dim], vol, buf[M], rhs[M][dim];
	INT I_[M][dim], J[dim][M];
	QUAD *quad;
	int q;
	const FLOAT *w, *p, *gu, *vu, *vz;

	bzero(A, sizeof(A));
	bzero(rhs, sizeof(rhs));
	bzero(buf, sizeof(buf));
	vol = phgGeomGetVolume(g, e);
	quad = phgQuadGetQuad3D(order);
	gu = phgQuadGetDofValues(e, gradu, quad); 
	vu = phgQuadGetDofValues(e, u, quad);        
	vz = phgQuadGetDofValues(e, dz, quad);     
	    
	p = quad->points;
	w = quad->weights;
	for (q = 0; q < quad->npoints; q++) {
	    /* Mat */
	    for (i = 0; i < M; i++) {
		for (j = 0; j < M; j++) {
		    const FLOAT *gi = phgQuadGetBasisValues(e, u, i, quad) + q; /* phi_i */
		    const FLOAT *gj = phgQuadGetBasisValues(e, u, j, quad) + q; /* phi_j */
		    FLOAT m = vol*(*w) * (*gi) * (*gj);
		    for (k = 0; k < dim; k++)
			A[i][k][j][k] += m;
		}
	    }

	    for (i = 0; i < M; i++) {
		const FLOAT *gi = phgQuadGetBasisValues(e, u, i, quad) + q; /* phi_i */
		for (k = 0; k < dim; k++) {
		    rhs[i][k] += vol*(*w) * (*gi) * (vu[k] + msl * gu[k*Dim+Z_DIR] * *vz
						     );
		}
	    }
	    
	    vu += dim; 
	    gu += Dim*dim; 
	    vz++;
	    w++; p += Dim + 1;
	}

	for (i = 0; i < M; i++)
	    for (k = 0; k < dim; k++)
		J[k][i] = I_[i][k] =
		    phgMapE2L(solver->rhs->map, 0, e, i * dim + k);

	for (i = 0; i < M; i++) {
	    for (k = 0; k < dim; k++) {
		if (phgDofDirichletBC_(u, e, i*dim+k, NULL, buf, &rhs[i][0],
				       DOF_PROJ_NONE)) {
		    phgMatAddEntries(mat, 1, I_[i] + k, M, J[k], buf);
		} else {
		    phgMatAddEntries(mat, 1, I_[i] + k, M*dim, I_[0],
				     &(A[i][k][0][0]));
		}
	    }
	}

	phgVecAddEntries(solver->rhs, 0, M * dim, I_[0], &rhs[0][0]);
    } /* end build rhs */

    phgSolverSolve(solver, FALSE, u, NULL);
    phgPrintf("   Dof %s moving solve: nits %d, res %e\n", 
	      u->name, solver->nits, solver->residual);
    phgDofFree(&gradu);

    phgSolverDestroy(&solver);
    phgMatDestroy(&mat);
    phgMapDestroy(&map);

    /* restore */
    if (u->DB_masks != NULL)
	memcpy(u->DB_masks, DB_mask_save, dim * sizeof(*DB_mask_save));
    else
	u->DB_mask = DB_mask_save[0];
    return;
}



