#include "ins.h"

#define SURF_DOF_TYPE (ns_params->utype)


/* Get local bases for boundary vertex i,
 *  return num of constrained bases,
 *  and noramlized direction of rotated bases.
 *  */
SURF_BAS *
get_surface_bases(NSSolver *ns, DOF_TYPE *u_type, int rdim)
{
    GRID *g = ns->g;
    SIMPLEX *e;
    DOF *surf_dof = NULL;
    DOF *dof_bottom_normal, *dof_lateral_normal;
    BOOLEAN *rotated = NULL;
    INT nrot = 0;
    surf_dof = phgDofNew(g, SURF_DOF_TYPE, rdim*rdim, "Surf bases", DofNoAction);
    rotated = phgCalloc(DofGetDataCount(surf_dof) / (rdim*rdim), sizeof(*rotated));
    SURF_BAS *surf_bas;
    int i;

    surf_bas = phgCalloc(1, sizeof(*surf_bas));
    surf_bas->type = u_type;
    surf_bas->dof = surf_dof;
    surf_bas->rotated = rotated;
    surf_bas->dim = rdim;
    phgDofSetDataByValue(surf_dof, 0.);
    //phgDofSetDataByFunction(surf_dof, func_eye);



    /* ----------------------------------------
     *
     *
     * Use grad b to compute dof normal
     *
     *
     * ---------------------------------------- */

    phgPrintf("   * Proj dof normal\n");
    dof_bottom_normal = get_bottom_normal(ns);
    dof_lateral_normal = get_lateral_normal(ns);


    //phgDofDump(dof_lateral_normal);
    if (0) {
	/* phgDofDump(dof_bottom_normal); */
	/* phgDofDump(dof_lateral_normal); */
	phgExportVTK(g, "normals.vtk", dof_bottom_normal, dof_lateral_normal, NULL);
	phgFinalize();
	exit(0);
    }



    if (rdim == 3) {
	/* ------------------------------------------------------------
	 *
	 *
	 * 3D rotation
	 *
	 *
	 * ------------------------------------------------------------ */

	/*
	 *
	 * First bottom direction
	 *
	 *  */
	ForAllElements(g, e) {
	    int s, ii, i, m;
	    int N = surf_dof->type->nbas;
	    FLOAT norm;
	    FLOAT bas_value[N][DDim],  H[Dim][Dim], c[Dim],
		    bxyz[Dim][Dim] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
	    FLOAT bottom_normal_value[N][Dim]; //, lateral_normal_value[N][Dim];


	    phgDofGetElementData(surf_dof, e, &bas_value[0][0]);
	    phgDofGetElementData(dof_bottom_normal, e, &bottom_normal_value[0][0]);
	    //phgDofGetElementData(dof_lateral_normal, e, &lateral_normal_value[0][0]);

	    for (s = 0; s < NFace; s++) {
		int nbas_face = NbasFace(surf_dof);
		INT id;
		SHORT ibas_face[nbas_face];

		/* Slip and lateral both need rotate. */
		if (!(e->bound_type[s] & BC_BOTTOM ))
		    continue;

		/* SHOW_M(&bottom_normal_value[0][0], N, Dim); */
		/* SHOW_M(&lateral_normal_value[0][0], N, Dim); */

		phgDofGetBasesOnFace(surf_dof, e, s, ibas_face);
		for (ii = 0; ii < nbas_face; ii++) {
		    i = ibas_face[ii];
		    /* Use Gram Schmidt process to get orthogonal bases,
		     * one constrains */

		    id = phgDofMapE2D(surf_dof, e, i * (DDim)) / (DDim);
		    rotated[id] = TRUE;
		    nrot++;

		    /* fisrt basis:
		     *   bottom out normal if at bottom,
		     *   z dir otherwize
		     */
		    memcpy(H[0], bottom_normal_value[i], Dim * sizeof(FLOAT));

		    /* second basis
		     *   lateral out normal if on lateral face
		     *   x or y dir otherwize
		     */
		    const FLOAT *dd = NULL;
		    m = 0;
		    for (; m < Dim; m++)
			if (fabs(c[0] = INNER_PRODUCT(H[0], bxyz[m])) < 0.9)
			    break;
		    assert(m < Dim);
		    dd = bxyz[m];
		    m++;

		    H[1][0] = dd[0] - c[0] * H[0][0];
		    H[1][1] = dd[1] - c[0] * H[0][1];
		    H[1][2] = dd[2] - c[0] * H[0][2];
		    norm = sqrt(INNER_PRODUCT(H[1], H[1]));
		    assert(norm > 1e-10);
		    H[1][0] /= norm;
		    H[1][1] /= norm;
		    H[1][2] /= norm;

		    /* third basis,
		     *  cross product
		     *
		     * */
#define CROSS_PRODUCT(a, b, n) {				\
			n[0] =  (a[1] * b[2] - b[1] * a[2]);	\
			n[1] = -(a[0] * b[2] - b[0] * a[2]);	\
			n[2] =  (a[0] * b[1] - b[0] * a[1]);	\
		    }
		    CROSS_PRODUCT(H[0], H[1], H[2]);
#undef CROSS_PRODUCT
		    norm = sqrt(INNER_PRODUCT(H[2], H[2]));
		    assert(norm > 1e-10);
		    H[2][0] /= norm;
		    H[2][1] /= norm;
		    H[2][2] /= norm;

		    memcpy(bas_value[i], H[0], DDim*sizeof(FLOAT));

		} /* end bas */
	    }     /* end face */
	    phgDofSetElementData(surf_dof, e, &bas_value[0][0]);
	}	      /* end elem */


	/*
	 *
	 * Then lateral, override the second direction in the bottom case
	 *
	 *  */
	ForAllElements(g, e) {
	    int s, ii, i, m;
	    int N = surf_dof->type->nbas;
	    FLOAT norm;
	    FLOAT bas_value[N][DDim],  H[Dim][Dim], c[Dim],
		    bxyz[Dim][Dim] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
	    FLOAT bottom_normal_value[N][Dim], lateral_normal_value[N][Dim];


	    phgDofGetElementData(surf_dof, e, &bas_value[0][0]);
	    phgDofGetElementData(dof_bottom_normal, e, &bottom_normal_value[0][0]);
	    phgDofGetElementData(dof_lateral_normal, e, &lateral_normal_value[0][0]);

	    for (s = 0; s < NFace; s++) {
		int nbas_face = NbasFace(surf_dof);
		INT id;
		SHORT ibas_face[nbas_face];

		/* Slip and lateral both need rotate. */
		if (!(e->bound_type[s] & BC_LATERAL ))
		    continue;

		/* SHOW_M(&bottom_normal_value[0][0], N, Dim); */
		//SHOW_M(&lateral_normal_value[0][0], N, Dim);

		phgDofGetBasesOnFace(surf_dof, e, s, ibas_face);
		//SHOW_iV_(0, ibas_face, nbas_face);
		for (ii = 0; ii < nbas_face; ii++) {
		    i = ibas_face[ii];
		    /* Use Gram Schmidt process to get orthogonal bases,
		     * one constrains */

		    id = phgDofMapE2D(surf_dof, e, i * (DDim)) / (DDim);
		    rotated[id] = TRUE;
		    nrot++;

		    /* fisrt basis:
		     *   bottom out normal if at bottom,
		     *   z dir otherwize
		     */
		    BTYPE btype = phgDofGetElementBasisInfo(surf_dof, e,
							    i, NULL, NULL, NULL);
		    if (btype & BC_BOTTOM) {
			memcpy(H[0], bottom_normal_value[i], Dim * sizeof(FLOAT));
		    }
		    else {
			H[0][X_DIR] = 0.;
			H[0][Y_DIR] = 0.;
			H[0][Z_DIR] = -1.;
		    }
		    //SHOW_M(H[0], Dim, 1);

		    /* second basis
		     *   lateral out normal if on lateral face
		     *   x or y dir otherwize
		     */
		    const FLOAT *dd = NULL;
		    dd = lateral_normal_value[i];
		    assert(INNER_PRODUCT(dd, dd) > 0.99999);
		    c[0] = INNER_PRODUCT(H[0], dd);

		    H[1][0] = dd[0] - c[0] * H[0][0];
		    H[1][1] = dd[1] - c[0] * H[0][1];
		    H[1][2] = dd[2] - c[0] * H[0][2];
		    norm = sqrt(INNER_PRODUCT(H[1], H[1]));
		    assert(norm > 1e-10);
		    H[1][0] /= norm;
		    H[1][1] /= norm;
		    H[1][2] /= norm;

		    //SHOW_M(H[1], Dim, 1);

		    /* third basis,
		     *  cross product
		     *
		     * */
#define CROSS_PRODUCT(a, b, n) {				\
			n[0] =  (a[1] * b[2] - b[1] * a[2]);	\
			n[1] = -(a[0] * b[2] - b[0] * a[2]);	\
			n[2] =  (a[0] * b[1] - b[0] * a[1]);	\
		    }
		    CROSS_PRODUCT(H[0], H[1], H[2]);
#undef CROSS_PRODUCT
		    norm = sqrt(INNER_PRODUCT(H[2], H[2]));
		    assert(norm > 1e-10);
		    H[2][0] /= norm;
		    H[2][1] /= norm;
		    H[2][2] /= norm;

		    memcpy(bas_value[i], H[0], DDim*sizeof(FLOAT));

		} /* end bas */
	    }     /* end face */
	    phgDofSetElementData(surf_dof, e, &bas_value[0][0]);
	}
    }
    else if (rdim == 2) {
	/* ------------------------------------------------------------
	 *
	 *
	 * 2D rotation
	 *
	 *
	 * ------------------------------------------------------------ */
	/*
	 *
	 * Only lateral
	 *
	 *  */
	ForAllElements(g, e) {
	    int s, ii, i, m;
	    int N = surf_dof->type->nbas;
	    FLOAT norm;
	    FLOAT bas_value[N][rdim*rdim],  H[rdim][rdim], c[rdim],
		    bxy[2][2] = {{1, 0}, {0, 1}};
	    FLOAT lateral_normal_value[N][Dim];


	    phgDofGetElementData(surf_dof, e, &bas_value[0][0]);
	    phgDofGetElementData(dof_lateral_normal, e, &lateral_normal_value[0][0]);

	    for (s = 0; s < NFace; s++) {
		int nbas_face = NbasFace(surf_dof);
		INT id;
		SHORT ibas_face[nbas_face];

		/* Lateral need rotate. */
		if (!(e->bound_type[s] & BC_LATERAL ))
		    continue;

		/* SHOW_M(&bottom_normal_value[0][0], N, Dim); */
		//SHOW_M(&lateral_normal_value[0][0], N, Dim);

		phgDofGetBasesOnFace(surf_dof, e, s, ibas_face);
		//SHOW_iV_(0, ibas_face, nbas_face);
		for (ii = 0; ii < nbas_face; ii++) {
		    i = ibas_face[ii];
		    /* Use Gram Schmidt process to get orthogonal bases,
		     * one constrains */

		    id = phgDofMapE2D(surf_dof, e, i * (rdim*rdim)) / (rdim*rdim);
		    rotated[id] = TRUE;
		    nrot++;

		    /* fisrt basis:
		     *   bottom out normal if at bottom,
		     *   z dir otherwize
		     */
		    memcpy(H[0], lateral_normal_value[i], rdim * sizeof(FLOAT));

		    /* second basis
		     */
		    const FLOAT *dd = NULL;
		    m = 0;
		    for (; m < rdim; m++)
			if (fabs(c[0] = INNER_PRODUCT2(H[0], bxy[m])) < 0.9) /* sqrt(2)/2 ? */
			    break;
		    assert(m < rdim);
		    dd = bxy[m];
		    m++;

		    H[1][0] = dd[0] - c[0] * H[0][0];
		    H[1][1] = dd[1] - c[0] * H[0][1];
		    norm = sqrt(INNER_PRODUCT2(H[1], H[1]));
		    assert(norm > 1e-10);
		    H[1][0] /= norm;
		    H[1][1] /= norm;


		    //SHOW_M(&H[0][0], rdim, rdim);

		    memcpy(bas_value[i], H[0], rdim * rdim *sizeof(FLOAT));
		} /* end bas */
	    }     /* end face */
	    phgDofSetElementData(surf_dof, e, &bas_value[0][0]);
	}
    }
    else {
	abort();
    }


    INT nrot0 = nrot;
    MPI_Allreduce(&nrot0, &nrot, 1, PHG_MPI_INT, MPI_SUM, g->comm);
    phgPrintf("   Rotate bdry (%d) Dofs: %d, ", rdim, nrot);

    //phgDofDump(surf_dof);
    PERIODIC_SYNC(surf_dof);

    if (0 && rdim == 3) {
	/* Plot check */
	DOF *dof_un = phgDofNew(g, SURF_DOF_TYPE, Dim, "UN", DofNoAction);
	DOF *dof_ln = phgDofNew(g, SURF_DOF_TYPE, Dim, "LN", DofNoAction);
	DOF *dof_fr = phgDofNew(g, SURF_DOF_TYPE, Dim, "free", DofNoAction);

	ForAllElements(g, e) {
	    int N = surf_dof->type->nbas;
	    FLOAT bas_value[N][DDim];
	    FLOAT bas_comp_value[N][Dim];
	    int k;

	    phgDofGetElementData(surf_dof, e, &bas_value[0][0]);

	    for (i = 0; i < N; i++)
		for (k = 0; k < Dim; k++)
		    bas_comp_value[i][k] = bas_value[i][k];
	    phgDofSetElementData(dof_un, e, &bas_comp_value[0][0]);

	    for (i = 0; i < N; i++)
		for (k = 0; k < Dim; k++)
		    bas_comp_value[i][k] = bas_value[i][3 + k];
	    phgDofSetElementData(dof_ln, e, &bas_comp_value[0][0]);

	    for (i = 0; i < N; i++)
		for (k = 0; k < Dim; k++)
		    bas_comp_value[i][k] = bas_value[i][6 + k];
	    phgDofSetElementData(dof_fr, e, &bas_comp_value[0][0]);
	}


	phgExportVTK(g, "rotate_bases.vtk", dof_un, dof_ln, dof_fr, NULL);
	/* phgFinalize(); */
	/* exit(1); */
    }

    if (0 && rdim == 2) {
	/* Plot check */
	DOF *dof_ln = phgDofNew(g, SURF_DOF_TYPE, Dim, "LN", DofNoAction);
	DOF *dof_fr = phgDofNew(g, SURF_DOF_TYPE, Dim, "free", DofNoAction);

	phgDofSetDataByValue(dof_ln, 0);
	phgDofSetDataByValue(dof_fr, 0);

	ForAllElements(g, e) {
	    int N = surf_dof->type->nbas;
	    FLOAT bas_value[N][rdim*rdim];
	    FLOAT bas_comp_value[N][Dim];
	    int k;

	    phgDofGetElementData(surf_dof, e, &bas_value[0][0]);
	    Bzero(bas_comp_value);

	    for (i = 0; i < N; i++)
		for (k = 0; k < rdim; k++)
		    bas_comp_value[i][k] = bas_value[i][k];
	    phgDofSetElementData(dof_ln, e, &bas_comp_value[0][0]);

	    for (i = 0; i < N; i++)
		for (k = 0; k < rdim; k++)
		    bas_comp_value[i][k] = bas_value[i][2 + k];
	    phgDofSetElementData(dof_fr, e, &bas_comp_value[0][0]);
	}


	phgExportVTK(g, "rotate_bases_2d.vtk", dof_ln, dof_fr, NULL);
	phgDofDump(surf_dof);
	phgDofDump(dof_ln);
	phgDofDump(dof_fr);
	phgFinalize();
	exit(1);
    }

    return surf_bas;
}







/*
 * A_{3, ncol} = Trans_{3,3} * A_{3, ncol}
 *  */
void
trans_left(FLOAT *A, int ncol, int lda, int rdim, const FLOAT *Trans)
{
    int i, j, k;
    FLOAT tmp[rdim][ncol];
    bzero(tmp, sizeof(tmp));
    for (i = 0; i < rdim; i++)
	for (j = 0; j < ncol; j++)
	    for (k = 0; k < rdim; k++)
		tmp[i][j] += Trans[i*rdim + k] * *(A + lda*k + j);
    for (i = 0; i < rdim; i++)
	for (j = 0; j < ncol; j++) {
	    *(A + lda*i + j) = tmp[i][j];
	}
}


/*
 * A_{3, ncol} = Trans_{3,3}^T * A_{3, ncol}
 *  */
void
trans_leftT(FLOAT *A, int ncol, int lda, int rdim, const FLOAT *Trans)
{
    int i, j, k;
    FLOAT tmp[rdim][ncol];
    bzero(tmp, sizeof(tmp));
    for (i = 0; i < rdim; i++)
	for (j = 0; j < ncol; j++)
	    for (k = 0; k < rdim; k++)
		tmp[i][j] += Trans[i + k*rdim] * *(A + lda*k + j);
    for (i = 0; i < rdim; i++)
	for (j = 0; j < ncol; j++)
	    *(A + lda*i + j) = tmp[i][j];
}


/*
 * A_{nrow, 3} = A_{nrow, 3} * Trans_{3,3}^T
 *  */
void
trans_rightT(FLOAT *A, int nrow, int lda, int rdim, const FLOAT *Trans)
{
    int i, j, k;
    FLOAT tmp[nrow][rdim];
    bzero(tmp, sizeof(tmp));
    for (i = 0; i < nrow; i++)
	for (j = 0; j < rdim; j++)
	    for (k = 0; k < rdim; k++)
		tmp[i][j] += *(A + lda*i + k) * Trans[k + j*rdim];
    for (i = 0; i < nrow; i++)
	for (j = 0; j < rdim; j++)
	    *(A + lda*i + j) = tmp[i][j];
}









void
rotate_dof_bases(DOF *u, SURF_BAS *surf_bas, BOOLEAN forward)
{
    int rdim = surf_bas->dim;
    INT i, N = DofGetDataCount(u) / rdim;
    DOF *surf_dof = surf_bas->dof;
    BOOLEAN *rotated = surf_bas->rotated;
    FLOAT *v = DofData(u);
    const FLOAT *trans = DofData(surf_dof);

    for (i = 0; i < N; i++) {
	if (rotated[i]) {
	    if (forward)
		trans_left(v, 1, 1, rdim, trans);
	    else
		trans_leftT(v, 1, 1, rdim, trans);
	}
	trans += rdim*rdim;
	v += rdim;
    }

    return;
}



void
dof_set_normal_data(DOF *u_h, SURF_BAS *surf_bas)
{
    GRID *g = u_h->g;
    SIMPLEX *e;

    DOF *surf_dof = surf_bas->dof;
    int rdim = surf_bas->dim;
    BOOLEAN *rotated = surf_bas->rotated;
    const FLOAT *Trans = DofData(surf_dof);

    assert(u_h->dim == rdim);

    ForAllElements(g, e) {
	int i, N = u_h->type->nbas;
	FLOAT vu[N][rdim], u0[Dim];

	phgDofGetElementData(u_h, e, vu[0]);

	for (i = 0; i < N; i++) {
	    INT id = phgDofMapE2D(surf_dof, e, i * (rdim*rdim)) / (rdim*rdim);
	    if (!rotated[id])
		continue;

	    const FLOAT *trans = Trans + id*(rdim*rdim);
	    FLOAT *coord = phgDofGetElementCoordinates(u_h, e, i*rdim);
	    func_u(coord[0], coord[1], coord[2], u0); /* Note: u0 should be of [Dim] */


	    BTYPE btype = phgDofGetElementBasisInfo(u_h, e,
						    i, NULL, NULL, NULL);
	    if (rdim == 3) {
		trans_left_3d  (vu[i], 1, 1, trans);
		if (btype & SLIP_BDRY) {
		    vu[i][UN_DIR] = INNER_PRODUCT(u0, trans + UN_DIR * rdim);
		}
		if (btype & BC_LATERAL) {
		    vu[i][LN_DIR] = INNER_PRODUCT(u0, trans + LN_DIR * rdim);
		}
		trans_leftT_3d (vu[i], 1, 1, trans);
	    }
	    else if (rdim == 2) {
		trans_left_2d  (vu[i], 1, 1, trans);
		if (btype & BC_LATERAL) { /* 2d: un is ln */
		    /* Note: use inner_product2 */
		    vu[i][UN_DIR] = INNER_PRODUCT2(u0, trans + UN_DIR * rdim); 
		}
		trans_leftT_2d (vu[i], 1, 1, trans);
	    }
	    else {
		abort();
	    }
	}

	phgDofSetElementData(u_h, e, vu[0]);
    }	      /* end elem */


    return;
}



void
update_floating(NSSolver *ns)
/* Update floating mask by height and bedrock topo */
{
    GRID *g = ns->g;
    LAYERED_MESH *gL = ns->gL;
    SIMPLEX *e;
    const FLOAT EPS = 1e-6;
    BTYPE *vert_mark;
    INT ii, i, k, s;

    vert_mark = phgCalloc(gL->nvert, sizeof(*vert_mark));

    /* Fisrt get the vert marks for Dirich BC */
    for (ii = 0; ii < gL->nvert_bot; ii++) {
	i = gL->vert_bot_Lidx[ii];
	assert(gL->vert_local_lists[i] != NULL);
	int iv = gL->vert_L2S[i];

	FLOAT h_topg, h_bottom;
	h_topg = gL->topg[iv];
	h_bottom = gL->bottom[iv];


	if (h_topg < h_bottom - EPS) {
	    /* floating */
	    vert_mark[iv] = BC_FLOAT;
	    phgInfo(0, "float_vert: %5d, %f %f %e\n", iv, h_topg, h_bottom, h_bottom - h_topg);
	}
	else {
	    /* touch ground */
	    vert_mark[iv] = SLIP_BDRY;
	}
    }



    /* Next update triangle marks for Neumann BC */
    double a[3];
    Bzero(a);

    ForAllElements(g, e) {
	for (k = 0; k < NFace; k++) {
	    int v[3+1], iS[3], iL[3];
	    if (!(e->bound_type[k] & BC_BOTTOM))
		continue;

            FLOAT area = phgGeomGetFaceArea(g, e, k);

	    /* Bottom face */
	    BOOLEAN floating = TRUE;
	    GetFaceVertices(e, k, v[0], v[1], v[2], v[3]);
	    for (i = 0; i < 3; i++) {
		iL[i] = e->verts[v[i]];
		iS[i] = gL->vert_L2S[iL[i]];
		assert(iS[i] >= 0 && iS[i] < gL->nvert);

		/* All three verts floats then the face is floating */
		if (vert_mark[iS[i]] & SLIP_BDRY)
		    floating = FALSE;
	    }

	    if (floating) {
		phgInfo(0, "float_face: %5d, %d\n", e->index, k);
		e->bound_type[k] &= ~SLIP_BDRY;
		e->bound_type[k] |= BC_FLOAT;
	    }
	    else {
		e->bound_type[k] &= ~BC_FLOAT;
		e->bound_type[k] |= SLIP_BDRY;
	    }

	    if (e->bound_type[k] & SLIP_BDRY)
		a[0] += area;
	    if (e->bound_type[k] & BC_FLOAT)
		a[1] += area;
	    a[2] += area;
	}
    }
    free(vert_mark);


#if USE_MPI
    {
        double b[3];
        MPI_Reduce(a, b, 3, MPI_DOUBLE, MPI_SUM, 0, g->comm);
        memcpy(a, b, sizeof(b));
    }
#endif

    phgPrintf("\tFloating area:\n");
    phgPrintf( "\tBottom  : %20.10e (%e)\n"
	       "\tSlip    : %20.10e\n"
	       "\tFloat   : %20.10e\n",
	       a[2], a[2] - a[0] - a[1], a[0], a[1]);
    if (phgRank == 0) {
	assert( fabs(a[2] - a[0] - a[1]) < 1e-8 );
    }


    phgUpdateBoundaryTypes(g);

    return;
}
