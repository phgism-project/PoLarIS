/*
 * Ice sheet - layers module.
 * Layers info is maintained for computaion involving surface-bottom. 
 *
 *
 * Prerequsite:
 * 1. 2D vert to 3D vert maping,
 *    L2Gmap_vert, 
 *    Z direction vert connection.
 * 2. 2D triangle to vert connection.
 * 3. 3D partation is based on 2D partation.
 *
 * Layers info:
 * 1. 2D vert to vert chain in Z dir.
 * 2. 2D triangle to elem chain in Z dir.
 * 3. Bottom Dof to Dof chain in Z difr. 
 *
 *
 *  */
#include "ins.h"
#include "layers.h"
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <netcdf.h>
#include "phg/netcdf-utils.h"
#include "vtk-draw.h"
 
#ifdef VTK_DEBUG
#else
int vtk_verb = 3;
char vtk_tmp_str[1000];
#endif



#define VTK_VERB(VERB) {				\
	verb = (VERB < vtk_verb) ? vtk_verb : VERB;	\
	phgInfo(0, "layers verb: %d\n", verb);		\
    }


/******************/
/* GEO DEFINATION */
/******************/
#define phgElemInfoCount 6

typedef enum {
    QUADRILATERAL_ = 0,
    TRIANGLE_	= 1,
    TETRA_	= 2,
    PYRAMID_	= 3,
    PRISM_	= 4,
    BRICK_	= 5
} ELEM_TYPE_;

/* typedef short EDGE_INFO[2];	/\* list of vertices *\/ */
/* typedef short FACE_INFO[5];	/\* 0: # of vertices, 1..: list of vertices *\/ */
/* typedef void (*ADD_TETRA_FUNC)(int v0, int v1, int v2, int v3, int bound_type[]); */
/* typedef void (*FUNC2TET)(int verts[], ADD_TETRA_FUNC func, int bound_type[]); */

/* typedef struct { */
/*     const char  *name;		/\* name of the element type *\/ */
/*     EDGE_INFO   *edge_info;	/\* list of edges *\/ */
/*     FACE_INFO   *face_info;	/\* list of faces *\/ */
/*     short       nvert;		/\* number of vertices *\/ */
/*     short	nedge;		/\* number of edges *\/ */
/*     short       nface;		/\* number of faces *\/ */
/*     FUNC2TET    func2tet;	/\* function to covert to tet *\/ */
/* } ELEM_INFO; */


/*****************/
/* Reading Utils */
/*****************/
#define READ_NUMBER						\
    if (!get_token(fp, token)) strcpy(token, "End ALL");	\
    if (isalpha((int)(token[0]))) {				\
	phgWarning("fewer entries (%d) than expected.\n", i);	\
	break;							\
    }

#define GET_TOKEN {					\
	if (!get_token(fp, token)) {			\
	    phgPrintf("error on line: %d\n", __LINE__);	\
	    goto error;					\
	}						\
    }

#define UNUSED_LINE { 				\
	fgets(token, 1024, fp);			\
    }
//	fprintf(stderr, "Unuesd: %s", token);	
#define NEXT_LINE UNUSED_LINE

#define UNUSED_LINE_CHECK(str) {				\
	fgets(token, 1024, fp);					\
	if (!strcmp(token, str))				\
	    phgError(-1, "Unmatched phase when read mesh file,"	\
		     " line:%d!\n", __LINE__);			\
    }
//	fprintf(stderr, "Unuesd: %s", token);	


static ELEM_INFO phg2DElemInfo_[] = {
    /* name	list of edges	list of faces	nvert	nedge	nface  func2tet*/
    {"quad",	NULL,	NULL,	4,	4,	4, NULL},
    {"tria",	NULL,	NULL,	3,	3,	3, NULL},
};

static int Gambit_Elemtype[phgElemInfoCount] = {
    2,				/* Quadrilateral */
    3,				/* Triangle */
    6, 				/* Tetrahedra */
    7, 				/* Pyramid */
    5, 				/* Wedge */
    4				/* Brick */
};

static int Gambit_Vertmap[phgElemInfoCount][100] = {
    {0, 1, 2, 3},		/* Quadrilateral */
    {0, 1, 2},			/* Triangle */
    {0, 1, 2, 3},               /* Tetrahedra */
    {0, 1, 3, 2, 4}, 		/* Pyramid */
    {0, 1, 2, 3, 4, 5},		/* Wedge */
    {0, 1, 3, 2, 4, 5, 7, 6}	/* Brick */
};

static int Gambit_Facemap[phgElemInfoCount][100] = {
    {0, 1, 2, 3},		/* Quadrilateral */
    {0, 1, 2},			/* Triangle */
    {3, 2, 0, 1},               /* Tetrahedra */
    {4, 0, 1, 2, 3}, 		/* Pyramid */
    {0, 1, 2, 3, 4, 5},		/* Wedge */
    {0, 1, 2, 3, 4, 5, 6, 7}	/* Brick */
};


static BOOLEAN
get_token(FILE *fp, char *token)
{
    int c;
    char *p;

    while (TRUE) {
	if (fscanf(fp, "%s", token) != 1)
	    return FALSE;
	if (token[0] != '#')
	    break;
	/* skip to newline */
	do
	    if ((c = fgetc(fp)) == EOF)
		return FALSE;
	while (c != '\n');
    }
    if ((p = strchr(token, '#')) != NULL)
	*p = '\0';
    return TRUE;
}

FLOAT 
get_face_area(FLOAT *p0, FLOAT *p1, FLOAT *p2)
{
    FLOAT x1 = p0[0], y1 = p0[1];
    FLOAT x2 = p1[0], y2 = p1[1];
    FLOAT x3 = p2[0], y3 = p2[1];

    return .5 * fabs(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2);
}

void 
get_face_normals(FLOAT *p0, FLOAT *p1, FLOAT *p2, FLOAT *normals)
{
    FLOAT x1 = p0[0], y1 = p0[1];
    FLOAT x2 = p1[0], y2 = p1[1];
    FLOAT x3 = p2[0], y3 = p2[1];

    FLOAT det = .5 * (x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2);
    char sgn = 1.;

    if (det > 0)
	sgn = 1.;
    else
	sgn = -1;

    FLOAT *n1 = normals + 0;	/* 1,2 */
    FLOAT *n2 = normals + 2;	/* 2,3 */
    FLOAT *n3 = normals + 4;	/* 3,1 */
    
    FLOAT x, y, d;

#define SET_NORM(X, Y, N)			\
    x = (X); y = (Y);				\
    d = sqrt(x*x + y*y);			\
    x /= d;  y /= d;				\
    N[0] = sgn*y; N[1] = -sgn*x;		

    SET_NORM(x2 - x1, y2 - y1, n1);
    SET_NORM(x3 - x2, y3 - y2, n2);
    SET_NORM(x1 - x3, y1 - y3, n3);
}



/*
 * Import (global) 2D mesh, and layer, using netcdf format.
 * Global Only.
 * 
 * */
LAYERED_MESH *
import_layered_mesh(const char *filename)
{
    LAYERED_MESH *gL;
    FILE *fp;
    char token[1024]; 
    INT i, j, n, k;
    FLOAT *R;
    int retval;
    int ncid, dimid, varid;
    size_t nvert, nelem, nbdry, nlayer;
    int dim2;
    int verb;

    PHG_CALLOC(gL, 1);
 


    /* Open the file. */
    if ((retval = nc_open(filename, NC_NOWRITE, &ncid)))
	NC_ERR(retval);

    /*
     * Read in dims
     * */
    NC_READ_DIM("nvert2d", nvert);
    NC_READ_DIM("nelem2d", nelem);
    NC_READ_DIM("nbdry2d", nbdry);
    NC_READ_DIM("nlayer", nlayer);



    /*********************/
    /* Nodal Coordinates */
    /*********************/
    phgInfo(0, "number of 2d vertices: %d\n", nvert);
    float *coord = calloc(nvert * 2, sizeof(float));
    NC_READ_FLOAT("verts2d", coord);
    NC_CHECK_DIM_2D("verts2d", nvert, 2);

    n = nvert;
    phgInfo(2, "number of 2d vertices: %"dFMT"\n", n);
    gL->verts = calloc(n, sizeof(*gL->verts));
    gL->vert_type = calloc(n, sizeof(*gL->vert_type));
    for (i = 0; i < n; i++) {
	FLOAT x, y, X1[3], X0[3];
	x = coord[2*i + 0] * 0.001;
	y = coord[2*i + 1] * 0.001;

	/* map z coord [0,1] to physical [base elevation, surf elevation]*/
	func_ice_slab(x, y, 0., X0); /* base */
	func_ice_slab(x, y, 1., X1); /* surf */

	gL->verts[i][0] = X0[X_DIR];     /* mapped x, should be the same x */
	gL->verts[i][1] = X0[Y_DIR];     /* mapped y,  */
	gL->verts[i][2] = X0[Z_DIR];     /* z base */
	gL->verts[i][3] = X1[Z_DIR];     /* z surf */

	gL->vert_type[i] = 0;	/* Unused */
    }
    gL->nvert = n;
    free(coord);
    

    /*****************************/
    /* Element/Cell Connectivity */
    /*****************************/
    n = nelem;
    dim2 = 3;			/* 3 verts of one triangle */
    phgInfo(2, "number of elements: %"dFMT"\n", n);
    int *elem_data = calloc(nelem * dim2, sizeof(int));
    NC_READ_INT("elems2d", elem_data);
    NC_CHECK_DIM_2D("elems2d", nelem, dim2);

    PHG_CALLOC(gL->trias, n);
    phgInfo(0, "number of triangles: %d\n", n);
    for (i = 0; i < n; i++) {
	int *verts = gL->trias[i].verts;
	verts[0] = elem_data[i * dim2 + 0];
	verts[1] = elem_data[i * dim2 + 1];
	verts[2] = elem_data[i * dim2 + 2];
	gL->trias[i].area = 
	    get_face_area(gL->verts[verts[0]],
			  gL->verts[verts[1]],
			  gL->verts[verts[2]]);
	get_face_normals(gL->verts[verts[0]],
			 gL->verts[verts[1]],
			 gL->verts[verts[2]], 
			 gL->trias[i].normals);
    }
    gL->ntria = n;
    free(elem_data);



    /*
     * Boundaries
     * */
    n = nbdry;
    dim2 = 3;			/* 3 verts of one triangle */
    phgInfo(2, "number of bdries: %"dFMT"\n", n);
    int *bdry_data = calloc(nbdry * dim2, sizeof(int));
    NC_READ_INT("bdrys2d", bdry_data);
    NC_CHECK_DIM_2D("bdrys2d", nbdry, dim2);

#define NORM2(x) sqrt(x[0]*x[0] + x[1]*x[1])

    EdgeMarks *bdrys;
    PHG_CALLOC(bdrys, nbdry);
    for (i = 0; i < nbdry; i++) {
	bdrys[i][X_DIR] = bdry_data[3*i + X_DIR];
	bdrys[i][Y_DIR] = bdry_data[3*i + Y_DIR];
	bdrys[i][2]     = bdry_data[3*i + 2];  /* mark */
	SORT_TWO_VERT(bdrys[i][X_DIR], bdrys[i][Y_DIR]);
    }
    qsort(bdrys, nbdry, sizeof(EdgeMarks), phgCompEdgeMark);
    
    /* compute region out normal */
    gL->normals = calloc(nvert, sizeof(*gL->normals));
    for (i = 0; i < nelem; i++) {
	const int edgeVerts[3][2] = {
	    {1, 2}, {2, 0}, {0, 1},
	    /* op: 0, 1, 2 */
	};
	int *verts = gL->trias[i].verts;

	for (j = 0; j < 3; j++) {
	    int v0 = verts[edgeVerts[j][0]];
	    int v1 = verts[edgeVerts[j][1]];
	    int vv = verts[j];
	    SORT_TWO_VERT(v0, v1);

	    /* find edge. */
	    EdgeMarks vts, *p;
	    vts[0] = v0;
	    vts[1] = v1;
	    p = bsearch(&vts, bdrys,
			nbdry, sizeof(EdgeMarks),
			phgCompEdgeMark);
	    if (p == NULL) {
		continue;
	    }
	    if ((*p)[2] != 3)  /* only lateral */
		continue;

	    FLOAT normal[2];
	    normal[X_DIR] = - (gL->verts[v0][Y_DIR] - gL->verts[v1][Y_DIR]);
	    normal[Y_DIR] = + (gL->verts[v0][X_DIR] - gL->verts[v1][X_DIR]);


	    FLOAT dir[2];
	    dir[X_DIR] = gL->verts[vv][X_DIR] - .5 * (gL->verts[v0][X_DIR] + gL->verts[v1][X_DIR]) ;
	    dir[Y_DIR] = gL->verts[vv][Y_DIR] - .5 * (gL->verts[v0][Y_DIR] + gL->verts[v1][Y_DIR]) ;

	    if (normal[X_DIR] * dir[X_DIR] + normal[Y_DIR] * dir[Y_DIR] > 0) {
		normal[X_DIR] *= -1;
		normal[Y_DIR] *= -1;
	    }

	    FLOAT length = NORM2(normal);
	    normal[X_DIR] /= length;
	    normal[Y_DIR] /= length;
	    
	    if (NORM2(gL->normals[v0]) > 0) {
		/* have value, take average */
		gL->normals[v0][X_DIR] = .5 * (gL->normals[v0][X_DIR] + normal[X_DIR]) ;
		gL->normals[v0][Y_DIR] = .5 * (gL->normals[v0][Y_DIR] + normal[Y_DIR]) ;
	    } else {
		gL->normals[v0][X_DIR] = normal[X_DIR];
		gL->normals[v0][Y_DIR] = normal[Y_DIR];
	    }

	    if (NORM2(gL->normals[v1]) > 0) {
		/* have value, take average */
		gL->normals[v1][X_DIR] = .5 * (gL->normals[v1][X_DIR] + normal[X_DIR]) ;
		gL->normals[v1][Y_DIR] = .5 * (gL->normals[v1][Y_DIR] + normal[Y_DIR]) ;
	    } else {
		gL->normals[v1][X_DIR] = normal[X_DIR];
		gL->normals[v1][Y_DIR] = normal[Y_DIR];
	    }
	    
	}
    }
    free(bdry_data);

    /* Re normalize */
    for (i = 0; i < nvert; i++) {
	FLOAT nn = NORM2(gL->normals[i]);
	if (NORM2(gL->normals[i]) > 0) {
	    gL->normals[i][X_DIR] /= nn;
	    gL->normals[i][Y_DIR] /= nn;
	}
    }

    
    /* plot check */
    if (0) {
	for (i = 0; i < nvert; i++) {
	    if (NORM2(gL->normals[i]) > 0) {
		printf("normal: %f %f\n%f %f\n\n", 
		       gL->verts[i][X_DIR], gL->verts[i][Y_DIR],
		       gL->verts[i][X_DIR] + 2 * gL->normals[i][X_DIR],
		       gL->verts[i][Y_DIR] + 2 * gL->normals[i][Y_DIR]);
	    }
	}
	for (i = 0; i < nbdry; i++) {
	    printf("bdrys: %f %f\n%f %f\n\n",
		   gL->verts[bdrys[i][0]][X_DIR],
		   gL->verts[bdrys[i][0]][Y_DIR],
		   gL->verts[bdrys[i][1]][X_DIR],
		   gL->verts[bdrys[i][1]][Y_DIR]
		   );
	}
    }

#undef NORM2

    

    /**************/
    /* Layer Info */
    /**************/
    phgInfo(2, "number of vert layers: %"dFMT"\n", nlayer);
    float *layers = calloc(nlayer, sizeof(float));
    NC_READ_FLOAT("layers", layers);

    gL->max_nlayer = nlayer - 1;	/* vert layers -> elem layers */
    gL->layer_ratio = phgCalloc(nlayer,
				sizeof(*gL->layer_ratio));
    for (i = 0; i < nlayer; i++) {
	gL->layer_ratio[i] = layers[i];
    }
    free(layers);


    /******************/
    /* Vert List Info */
    /******************/
    PHG_CALLOC(gL->vert_global_lists, nvert);
    for (i = 0; i < nvert; i++) {
	int nv = nlayer;	/* # vert layers */
	PHG_CALLOC(gL->vert_global_lists[i], nv + 1);
	
	/* fisrt as num of verts in list */
	gL->vert_global_lists[i][0] = nv; 

	for (j = 0; j < nv; j++) {
	    gL->vert_global_lists[i][j+1] = i * nlayer + j;
	}
    }
 
    if ((retval = nc_close(ncid)))
	NC_ERR(retval);
    
    return gL;
}



void
build_layered_mesh(GRID *g, LAYERED_MESH *gL)
/*
 *
 * Build layered mesh, globally or locally  
 *
 *  */
{
    SIMPLEX *e;
    TRIA *t;
    int i, j, k;
    int verb = 0;

    gL->g = g;

    /* ------------------------------------------------------------
     *
     *  Step 1:
     *    Build vertex relation: 
     *
     * ------------------------------------------------------------ */

    /* Step 1.1:
     *   3D mesh vertex connection */
    VEFMAP_3D *vef_3d;
    PHG_CALLOC(vef_3d, g->nvert); 
    gL->vef_3d = vef_3d;
    for (i = 0; i < g->nvert; i++) {
	vef_3d[i].size = 20;
	PHG_CALLOC(vef_3d[i].v, vef_3d[i].size);
    }

    ForAllElements(g, e) {
	VEFMAP_3D *vf;
	int v[NVert];
	for (i = 0; i < NVert; i++)
	    v[i] = e->verts[i];
	for (i = 0; i < NVert; i++) {
	    vf = vef_3d + v[i];
	    for (j = 0; j < NVert; j++) {
		for (k = 0; k < vf->nv; k++) 
		    if (v[j] == vf->v[k])
			break;
		if (k == vf->nv) {	/* new point connection
					 * v[i] ~ v[j]
					 * */
		    if (vf->nv >= vf->size) {
			int old_size = vf->size;
			vf->size += 20;
			PHG_REALLOC(vf->v,
				    vf->size, old_size);
		    }
		    vf->v[vf->nv++] = v[j];
		}
	    }
	}
    }

    /* /\* debug *\/ */
    /* for (i = 0; i < g->nvert; i++) { */
    /* 	phgInfo(0, "vef3d [%d]\n", i); */
    /* 	SHOW_iV_(0, vef_3d[i].v, vef_3d[i].nv); */
    /* } */


    /*
     * Step 1.2:
     *   Build Vert L2S
     * */
    int *vert_list0, *found;
    /* map bottom verts (2D) to global (3D) */
    PHG_CALLOC(vert_list0, gL->nvert); /* only first vert,
					* the one on bottom */
    for (i = 0; i < gL->nvert; i++) {
	vert_list0[i] = gL->vert_global_lists[i][1];
	if (i > 0)
	    assert(vert_list0[i] > vert_list0[i-1]);
    }
    //SHOW_iV_(0, vert_list0, gL->nvert);

    /*
     * Map local verts (3D) to (2D)
     *
     * gL->vert_bot_Lidx[n_local_vert_bot]: -> 3d local index
     * gL->vert_bot_Gidx[n_local_vert_bot]: -> 3d global index
     *
     * 
     * Note:
     *   assume gL->vert_global_list[.][1] is increasing !!! */
    gL->nvert_bot = 0;
    int nvb = 0;		/* # of bottom vert */
    for (i = 0; i < g->nvert; i++) 
	if (g->types_vert[i] & BC_BOTTOM) 
	    nvb++;
    gL->nvert_bot = nvb;
    PHG_CALLOC(gL->vert_bot_Lidx, gL->nvert_bot);
    PHG_CALLOC(gL->vert_bot_Gidx, gL->nvert_bot);
    nvb = 0;
    for (i = 0; i < g->nvert; i++) 
	if (g->types_vert[i] & BC_BOTTOM) 
	    gL->vert_bot_Lidx[nvb++] = i;

    /*
     *
     * gL->vert_L2S[n_local_vert_3d]:   -> 2d (global) index,  colume to bottom 1
     * gL->vert_S2L[n_lcoal_vert_bot]:  -> 3d local index,     1 to 1   
     * gL->vert_local_lists[n_lcoal_vert_bot]: ???
     *                                  -> To change
     */

    nvb = 0;
    PHG_CALLOC(gL->vert_L2S, g->nvert);
    PHG_CALLOC(gL->vert_S2L, gL->nvert);
    PHG_CALLOC(gL->vert_local_lists, g->nvert);

    for (i = 0; i < gL->nvert; i++) {
	gL->vert_S2L[i] = -1;
    }

    for (i = 0; i < g->nvert; i++) {
	if (g->types_vert[i] & BC_BOTTOM) {
	    //printf(" BC vert[%d]\n", i);
	    INT iS, iG = GlobalVertex(g, i);
	    found = bsearch(&iG, vert_list0,
			    gL->nvert, sizeof(INT), phgCompINT);
	    assert(found != NULL);
	    gL->vert_L2S[i] =
		iS = found - vert_list0;
	    gL->vert_S2L[iS] = i;

	    /* Map vG to vL */
	    INT nly, *vG, *vL;
	    gL->vert_bot_Gidx[nvb] = iS;
	    vG = gL->vert_global_lists[iS];
	    nly = vG[0];
	    //SHOW_iV_(0, vG, nly+1);
	    PHG_CALLOC(vL, nly+1);
	    gL->vert_local_lists[i] = vL;
	    vL[0] = nly;		/* first as num */
	    vL[1] = i;			/* second as self */
	    for (j = 2; j < nly+1; j++) {
		/* In the neighbor of prev vert will
		 * be the current vert. */
		VEFMAP_3D *vf = vef_3d + vL[j-1];
		for (k = 0; k < vf->nv; k++) {
		    if (vG[j] == GlobalVertex(g, vf->v[k])) 
			break;
		}
		if (k >= vf->nv) {
		    int q;
		    for (q = 0; q < vf->nv; q++) {
			phgInfo(0, "%d (%d), ", 
				vf->v[q],
				GlobalVertex(g, vf->v[q]));
		    }
		    phgInfo(0, "\nvG[%d]:%d not found to vL\n", j, vG[j]);
		}
		assert(k < vf->nv);
		vL[j] = vf->v[k];
	    }
	    //SHOW_iV_(0, vL, nly+1);
	    //phgInfo(0, "    add: %x\n", vL);
	    nvb++;
	} else {
	    //printf(" Ly vert[%d]\n", i);
	    gL->vert_L2S[i] = -1;
	    gL->vert_local_lists[i] = NULL;
	}
	/* phgInfo(0, "vert local list %d: %x\n",  */
	/* 	i, gL->vert_local_lists[i] */
	/* 	); */
    }
    

    /* 2D verts to 3D vert lists, ref element */
    {
	VEF_MAP2 *vef = NULL;
	vef = phgDofSetupVEFMap2(g, NULL, VERT_FLAG);
	PHG_CALLOC(gL->vert_elements, g->nvert);

	for (i = 0; i < g->nvert; i++) {
	    if (g->types_vert[i] & BC_BOTTOM) {
		assert(gL->vert_local_lists[i] != NULL);
		
		SIMPLEX **elist, *e;
		int nv = gL->vert_local_lists[i][0];
		int *iL = &gL->vert_local_lists[i][1];
		PHG_CALLOC(elist, nv);
		gL->vert_elements[i] = elist;
		
		for (j = 0; j < nv-1; j++) {
		    INT v0 = iL[j];
		    INT v1 = iL[j+1];
		    int j0, j1;
		    
		    e = NULL;
		    for (j0 = 0; j0 < vef->Vsize[v0]; j0++) {
			for (j1 = 0; j1 < vef->Vsize[v1]; j1++) {
			    if (vef->Vmap[v0][j0] == vef->Vmap[v1][j1]) 
				e = vef->Vmap[v0][j0];
			}
		    }
		    assert(e != NULL);

		    elist[j] = e;
		} /* segments */

	    }
	}

	phgDofFreeVEFMap2(&vef);
    }

    /*
     * Create gather communicator
     *
     *  */
    if (1) {
	static int pass = 0;
	int rank, *vert_bot_cnts, *vert_bot_dsps, *vert_bot_idxs, ntotal;

	phgInfo(0, "[%d] Build vert 2D gather comm: pass %d\n", 
		phgRank, pass++);
        PHG_CALLOC(vert_bot_cnts, g->nprocs);
        PHG_CALLOC(vert_bot_dsps, g->nprocs);
	MPI_Allgather(&gL->nvert_bot, 1, MPI_INT, 
		      vert_bot_cnts, 1, MPI_INT, g->comm);
	ntotal = 0;
	for (rank = 0; rank < g->nprocs; rank++) {
	    vert_bot_dsps[rank] = ntotal;
	    ntotal += vert_bot_cnts[rank];
	}
	PHG_CALLOC(vert_bot_idxs, ntotal);
	MPI_Allgatherv(gL->vert_bot_Gidx, gL->nvert_bot, PHG_MPI_INT,
		       vert_bot_idxs, vert_bot_cnts, vert_bot_dsps, 
		       MPI_INT, g->comm);

	
	gL->vert_bot_total = ntotal;
	gL->vert_bot_cnts = vert_bot_cnts;
	gL->vert_bot_dsps = vert_bot_dsps;
	gL->vert_bot_idxs = vert_bot_idxs;

	/* SHOW_iV_(0, gL->vert_bot_cnts, g->nprocs); */
	/* SHOW_iV_(0, gL->vert_bot_dsps, g->nprocs); */
	/* SHOW_iV_(0, gL->vert_bot_idxs, */
	/* 	 gL->vert_bot_total); */
    }


    /* ------------------------------------------------------------
     *
     *  Step 1:
     *    Build triangle relation: VEFmap_2d, mapE3to2
     *
     * ------------------------------------------------------------ */


    /* Build VEFmap_2d */
    VEFMAP_2D *vef_2d;
    PHG_CALLOC(vef_2d, gL->nvert);
    gL->vef_2d = vef_2d;
    t = gL->trias;
    for (j = 0; j < gL->ntria; j++, t++) {
	for (k = 0; k < 3; k++) {
	    i = t->verts[k];
	    assert(i < gL->nvert);
	    int ne = vef_2d[i].ne;
	    assert(ne < 12);
	    vef_2d[i].idx[ne] = j;
	    vef_2d[i].ne++;
	    //SHOW_iV(vef_2d[i].idx, vef_2d[i].ne);
	}
    }

    /* Build mapE3to2 */
    int nface_bot = 0;
    BOTTOM_FACE *face_bot, *fb;

    ForAllElements(g, e) 
	for (k = 0; k < NFace; k++) 
	    if (e->bound_type[k] & BC_BOTTOM)
		nface_bot++;
    PHG_CALLOC(face_bot, nface_bot);
    phgInfo(0, "Bottom face: %d\n", nface_bot);
    gL->nface_bot = nface_bot;
    gL->face_bot = face_bot;

    VTK_VERB(3);
    fb = face_bot;
    ForAllElements(g, e) {
	for (k = 0; k < NFace; k++) {
	    int v[3+1], iS[3], iL[3], ii;
	    if (!(e->bound_type[k] & BC_BOTTOM))
		continue;

	    /* Bottom face */
	    GetFaceVertices(e, k, v[0], v[1], v[2], v[3]);
	    for (i = 0; i < 3; i++) {
		iL[i] = e->verts[v[i]];
		iS[i] = gL->vert_L2S[iL[i]];
		assert(iS[i] >= 0 && iS[i] < gL->nvert);
	    }

	    /* find tria contains v[012] */
	    int n0, n1, n2, *t0, *t1, *t2,
		i0, i1, i2;

	    //printf("v: [%d %d %d]\n", iS[0], iS[1], iS[2]);
	    n0 = vef_2d[iS[0]].ne;
	    n1 = vef_2d[iS[1]].ne;
	    n2 = vef_2d[iS[2]].ne;
	    t0 = vef_2d[iS[0]].idx;
	    t1 = vef_2d[iS[1]].idx;
	    t2 = vef_2d[iS[2]].idx;
	    //SHOW_iV(t0, n0);
	    //SHOW_iV(t1, n1);
	    //SHOW_iV(t2, n2);
#if 1
	    /* faster */
	    t0 = vef_2d[iS[0]].idx;
	    for (i0 = 0; i0 < n0; i0++, t0++) {
		t1 = vef_2d[iS[1]].idx;
		for (i1 = 0; i1 < n1; i1++, t1++) {
		    if (*t0 == *t1) {
			t2 = vef_2d[iS[2]].idx;
			for (i2 = 0; i2 < n2; i2++, t2++) {
			    if (*t0 == *t2) {
				//printf("   found: %d\n", *t0);
				break;
			    }
			}
			if (i2 < n2)	/* found */
			    break;	
		    }
		}		
		if (i1 < n1)		/* found */
		    break;
	    }
#else
	    t0 = vef_2d[iS[0]].idx;
	    for (i0 = 0; i0 < n0; i0++, t0++) {
		t1 = vef_2d[iS[1]].idx;
		for (i1 = 0; i1 < n1; i1++, t1++) {
		    t2 = vef_2d[iS[2]].idx;
		    for (i2 = 0; i2 < n2; i2++, t2++) {
			if (*t0 == *t2 && *t0 == *t1) {
			    printf("   found: %d\n", *t0);
			    break;
			}
		    }
		    if (i2 < n2)	/* found */
			break;
		}		
		if (i1 < n1)		/* found */
		    break;
	    }
#endif

	    if (i0 >= n0) {
		printf("can't find triangle in 2D mesh shared by verts: [%d %d %d]\n",
		       iS[0], iS[1], iS[2]);
	    }
	    assert (i0 < n0);		/* found */
	    fb->e = e;
	    fb->face = k;
	    for (j = 0; j < 3; j++)
		fb->vert[j] = iL[j];
	    fb->tria = *t0;

	    /* debug */
	    {
		FLOAT c0[3], c1[3];
		Bzero(c0);
		Bzero(c1);
		for (i = 0; i < 3; i++)
		    for (j = 0; j < 3; j++) {
			c0[j] += g->verts[iL[i]][j] / 3.;
			c1[j] += gL->verts[iS[i]][j] / 3.; 
		    }
		vtkSetColor(verb, "yellow");
		vtkDrawLine(verb, c0, c1);
	    }

	    fb++;
	}
    }


    /* ------------------------------------------------------------
     *
     *  Step 2:
     *   Build tetra relation:    VEFmap_3d, mapE2to3
     *
     * ------------------------------------------------------------ */
    phgInfo(0, "Tria to tetra relation\n");
    VTK_VERB(3);
    //assert(nface_bot == gL->ntria); /* Fixme */
    fb = face_bot;
    e = NULL;
    for (i = 0; i < nface_bot; i++, fb++) {
	SIMPLEX *e1 = fb->e, *e0 = NULL;
	int f1 = fb->face, f0 = -1;
	int *vL[3], n[3], iv[3], ne = 0;
	
	//phgInfo(0, "face bot %d\n", i);

	for (k = 0; k < 3; k++) {
	    assert(g->types_vert[fb->vert[k]] & BC_BOTTOM);
	    //printf(" BC vert[%d]\n", fb->vert[k]);
	    vL[k] = gL->vert_local_lists[fb->vert[k]];
	    //printf("    add: %x\n", vL[k]);
	    assert(vL[k] != NULL);
	    n[k] = vL[k][0];
	    iv[k] = 1;		/* start */
	    //SHOW_iV(vL[k], n[k]+1);
	    if (n[k] > ne)
		ne = n[k];
	}
	
	/* tets in layers */
	ne *= 3;
	PHG_CALLOC(fb->elems, ne);
	PHG_CALLOC(fb->faces, 3*ne);
	fb->ne = 0;

	/* Algorithm desp:
	 * Given tet(e0) and face(f1), find the next face(f1)
	 *   and next tet(e1) in layer.
	 * 1. find the vert not on the face, it must belongs
	 *    to one of three vertical line which contains,
	 *    say, v(j).
	 * 2. find the face opsite to v(j), the next tet
	 *    is the neigh of current tet on face v(j)
	 *
	 * */

	while (TRUE) {
	    int v[4], ff;

	    e0 = e1;
	    f0 = f1;
	    fb->elems[fb->ne++] = e0;

	    vtkSetColor(verb, "yellow");
	    vtkDrawElement(verb, e0);

	    for (k = 0; k < 3; k++) 
		v[k] = GetFaceVertex(f0, k);
	    v[3] = f0;
	    for (k = 0; k < 3; k++)
		if (iv[k]+1 < n[k]+1
		    && vL[k][iv[k]+1] == e0->verts[f0])
		    break;

	    assert(k < 3);	/* advance on line[k] */
	    for (j = 0; j < 3; j++)
		if (e0->verts[v[j]] == vL[k][iv[k]])
		    break;
	    assert(j < 3);	/* elem v[j] on line vL[k],
				 * and the next face */
	    iv[k]++;
	    ff = v[j];	

	    vtkSetColor(verb, "red");
	    vtkDrawFace(verb, e0, ff);
	    //vtkPause(0);

	    /* reach top */
	    if (e0->neighbours[ff] == NULL) {
		break;
	    } else {
		assert(e0->bound_type[ff] & INTERIOR);
		e1 = e0->neighbours[ff];
		for (k = 0; k < NFace; k++)
		    if (e1->faces[k] == e0->faces[ff])
			break;
		assert(k < NFace);
		f1 = k;
	    }
	} /* end layers */
	vtkPause(verb);
	assert(fb->ne == gL->max_nlayer * 3);


	/* { */
	/*     int i_; */
	/*     phgInfo(0, "Debug face bot [%d]\n", fb->ne); */
	/*     for (i_ = 0; i_ < fb->ne / 3; i_++) { */
	/* 	phgInfo(0, "%d_%d, ", fb->elems[i_*3 + 2]->index, fb->elems[i_*3 + 2]->index - (fb->tria * gL->max_nlayer * 3 + i_*3 + 0)); */
	/* 	phgInfo(0, "%d_%d, ", fb->elems[i_*3 + 1]->index, fb->elems[i_*3 + 1]->index - (fb->tria * gL->max_nlayer * 3 + i_*3 + 1)); */
	/* 	phgInfo(0, "%d_%d, ", fb->elems[i_*3 + 0]->index, fb->elems[i_*3 + 0]->index - (fb->tria * gL->max_nlayer * 3 + i_*3 + 2)); */
	/* 	/\* assert(fb->elems[i_*3 + 2]->index == fb->tria * gL->max_nlayer * 3 + i_*3 + 0); *\/ */
	/* 	/\* assert(fb->elems[i_*3 + 1]->index == fb->tria * gL->max_nlayer * 3 + i_*3 + 1); *\/ */
	/* 	/\* assert(fb->elems[i_*3 + 0]->index == fb->tria * gL->max_nlayer * 3 + i_*3 + 2); *\/ */
	/*     } */
	/*     phgInfo(0, "\n"); */
	/* } */

	/* Flux face */
	if (1) {
	    /* Normals */
	    TRIA *t = gL->trias + fb->tria;
	    FLOAT *normals = t->normals;
	    int ii;

	    for (ii = 0; ii < fb->ne; ii++) { /* elems */
		SIMPLEX *e = fb->elems[ii];
		
		for (j = 0; j < 3; j++) { /* tri faces */
		    FLOAT *n = normals + j * 2;
		    for (k = 0; k < NFace; k++) {
			const FLOAT *nn
			    = phgGeomGetFaceOutNormal(g, e, k);
			FLOAT nt[3] = {n[0], n[1], 0.};
			FLOAT r = INNER_PRODUCT(nn, nt);
			if (1. - r < 1e-7)
			    break;
		    }
		    if (k < NFace) {
			fb->faces[ii*3 + j] = k;

			/* check */
			int kk = 0;
			for (kk = 0; kk < j-1; kk++)
			    assert(k != fb->faces[ii*3 + kk]);
			    
		    } else {
			fb->faces[ii*3 + j] = -99;
		    }
		} /* end tria face */
	    }	  /* end elements */
	}	  /* end flux */
    }

    phgFree(vert_list0);



#if 1
    /* ------------------------------------------------------------
     * 
     *    Dual mesh
     *
     *  Build fb->edge2to3[ne][3][2]
     *    triagle edges [3] to tetra edges [ne][2]
     *    if tri_edge -> only one tet edge
     *    then edge2to3[ii][k][1] = -1 
     *
     * ------------------------------------------------------------
     * */
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
	FLOAT mid_tri[3][2];
	FLOAT mid_tet[NEdge][3];
	int verb_edge = 3;

	PHG_CALLOC(fb->edge2to3, 6 * fb->ne);

	for (k = 0; k < 3; k++) { /* edge */
	    mid_tri[k][0] = .5*(  x_tri[ e2v_tri[k][0] ][0]
				+ x_tri[ e2v_tri[k][1] ][0]);
	    mid_tri[k][1] = .5*(  x_tri[ e2v_tri[k][0] ][1]
				+ x_tri[ e2v_tri[k][1] ][1]);
	    phgInfo(verb_edge, "mid_tri: %e %e\n",
		    mid_tri[k][0], mid_tri[k][1]);
	}


	
	for (ii = 0; ii < fb->ne; ii++) { /* tets */
	    SIMPLEX *e = fb->elems[ii];
	    int edge3to2[NEdge];
	    int edge2to3[3][2];
	    
	    /* tetra edge */
	    for (j = 0; j < NEdge; j++) {
		int v0 = GetEdgeVertex(j, 0);
		int v1 = GetEdgeVertex(j, 1);
		FLOAT *x0 = g->verts[e->verts[v0]];
		FLOAT *x1 = g->verts[e->verts[v1]];
		
		mid_tet[j][0] = .5 * (x0[0] + x1[0]);
		mid_tet[j][1] = .5 * (x0[1] + x1[1]);
		mid_tet[j][2] = .5 * (x0[2] + x1[2]);

		phgInfo(verb_edge, "mid_tet: %e %e\n",
			mid_tet[j][0], mid_tet[j][1]);
		
		for (k = 0; k < 3; k++) {
		    FLOAT d[2] = {mid_tet[j][0] - mid_tri[k][0], 
				  mid_tet[j][1] - mid_tri[k][1]};
		    FLOAT dist = sqrt(d[0]*d[0] + d[1]*d[1]);
		    phgInfo(verb_edge, "dist: %e \n", dist);
		    if (dist < 1e-7) {
			break;
		    }
		}
		if (k < 3) {
		    edge3to2[j] = k; /* found */
		} else {
		    edge3to2[j] = -1;
		}
	    }

	    phgInfo(verb_edge, "edge3to2: %d %d %d %d %d %d\n",
		    edge3to2[0], edge3to2[1], edge3to2[2],
		    edge3to2[3], edge3to2[4], edge3to2[5]
		    );
	    
	    {
	    	/*
	    	 * check: edge3to2
	    	 *     tet edge map to tri edge
	    	 *     e.g. 2 k1, 2 k2, only 1 k3
	    	 * */
	    	int nk[3] = {0, 0, 0}, k0;
	    	for (j = 0; j < NEdge; j++)
		    if (edge3to2[j] >= 0)
			nk[edge3to2[j]]++;
		phgInfo(verb_edge, "nk %d %d %d\n", nk[0], nk[1], nk[2]);
	    	for (k = 0; k < 3; k++)
	    	    if (nk[k] == 1)
	    		break;
	    	assert(k < 3);
	    	k0 = k;
	    	for (k = 0; k < 3; k++)
	    	    if (k != k0)
	    		assert(nk[k] == 2);
	    }
	 
	    for (k = 0; k < 3; k++) {
		int v0 = -1, v1 = -1;
		for (j = 0; j < NEdge; j++)
		    if (edge3to2[j] == k) {
			if (v0 == -1)
			    v0 = j;
			else 
			    v1 = j;
		    }
		assert(v0 != -1);
		if (v1 == -1) {
		    edge2to3[k][0] = v0;
		    edge2to3[k][1] = -1; /* note 2nd one */
		} else {
		    edge2to3[k][0] = v0;
		    edge2to3[k][1] = v1;
		}
	    }

	    for (k = 0; k < 3; k++) {
		fb->edge2to3[ii * 6 + k*2  ] = edge2to3[k][0];
		fb->edge2to3[ii * 6 + k*2+1] = edge2to3[k][1];
	    }
	}
    }		  /* end bot face */


    /* Volume of prism colume */
    phgInfo(0, "* Initial voulme\n");
    FLOAT grid_coord_unit = ns_params->grid_coord_unit;
    gL->volumes = phgCalloc(gL->ntria, sizeof(FLOAT));
    t = gL->trias;
    for (j = 0; j < gL->ntria; j++, t++) {
	int I_[3] = {t->verts[0],
		    t->verts[1],
		    t->verts[2]};
	FLOAT *coord[3];

	for (i = 0; i < 3; i++) {
	    coord[i] = gL->verts[I_[i]];
	}
	FLOAT x, y, z, X[3];
	x = (coord[0][0]+ coord[1][0] + coord[2][0]) / 3.; /* [km] */
	y = (coord[0][1]+ coord[1][1] + coord[2][1]) / 3.; /* [km] */
	z = 1.;			/* Note: height assume to be 1 */

	/* Note:
	 *  this map workes only with [grid_coord_unit],
	 *  here x, y is already in [grid_coord_unit] (km).
	 *  */
	if (0) {
#warning -- FIXME !!! ----	    
	    //x /= grid_coord_unit;
	    //y /= grid_coord_unit;
	    func_ice_slab(x, y, z, X);
	    // input:  km 
	    // output: km
	    //x *= grid_coord_unit;
	    //y *= grid_coord_unit;
	}
	else {
	    X[2] = 1.;
	}
	gL->volumes[j] = t->area * X[2];
	phgInfo(3, " X[%e %e %e], height %e\n", x, y, z, X[2]);
	phgInfo(3, "vol[%3d] %e\n", j, gL->volumes[j]);
    }    
#else
#  warning Daul mesh disabled!
#endif

#if 1
    /* Height of points */
    build_layered_mesh_height(g, gL);
    check_height(g, gL);
#else
    /* Use ice-grid. */
#endif

    return;
}


void 
build_layered_mesh_height(GRID *g, LAYERED_MESH *gL)
{
    SIMPLEX *e;
    TRIA *t;
    int i, j, k;
    int verb = 0;

    gL->height =
	phgCalloc(gL->nvert, sizeof(*gL->height));
    gL->height_local =
	phgCalloc(gL->nvert_bot, sizeof(*gL->height_local));
    gL->height_gather =
	phgCalloc(gL->vert_bot_total, sizeof(*gL->height_gather));

    gL->topg =
	phgCalloc(gL->nvert, sizeof(*gL->topg));
    gL->topg_local =
	phgCalloc(gL->nvert_bot, sizeof(*gL->topg_local));
    gL->topg_gather =
	phgCalloc(gL->vert_bot_total, sizeof(*gL->topg_gather));

    gL->bottom =
	phgCalloc(gL->nvert, sizeof(*gL->bottom));
    gL->bottom_local =
	phgCalloc(gL->nvert_bot, sizeof(*gL->bottom_local));
    gL->bottom_gather =
	phgCalloc(gL->vert_bot_total, sizeof(*gL->bottom_gather));
    
    
    INT nvb = 0;
    for (i = 0; i < g->nvert; i++) {
	if (g->types_vert[i] & BC_BOTTOM) {
	    assert(gL->vert_local_lists[i] != NULL);

	    int nv = gL->vert_local_lists[i][0];
	    int *iL = &gL->vert_local_lists[i][1];
	    gL->height_local[nvb] = 
		g->verts[iL[nv-1]][Z_DIR];  /* surf */
	    gL->bottom_local[nvb] = 
		g->verts[iL[0]][Z_DIR];	 /* base */

	    gL->topg_local[nvb] = func_ice_topg(g->verts[iL[0]][X_DIR],
						g->verts[iL[0]][Y_DIR]);
	    
	    phgInfo(3, "  Asgn Bidx: %5d, Gidx: %5d height: %e\n",
		    gL->vert_bot_Gidx[nvb],
		    GlobalVertex(g, iL[nv-1]),
		    g->verts[iL[nv-1]][Z_DIR]);
	    
	    nvb++;
	}
    }

    MPI_Allgatherv(gL->height_local, gL->nvert_bot, PHG_MPI_FLOAT,
		   gL->height_gather, gL->vert_bot_cnts, gL->vert_bot_dsps, 
		   PHG_MPI_FLOAT, g->comm);
    MPI_Allgatherv(gL->bottom_local, gL->nvert_bot, PHG_MPI_FLOAT,
		   gL->bottom_gather, gL->vert_bot_cnts, gL->vert_bot_dsps, 
		   PHG_MPI_FLOAT, g->comm);
    MPI_Allgatherv(gL->topg_local, gL->nvert_bot, PHG_MPI_FLOAT,
		   gL->topg_gather, gL->vert_bot_cnts, gL->vert_bot_dsps, 
		   PHG_MPI_FLOAT, g->comm);


    for (i = 0; i < gL->vert_bot_total; i++) {
	INT idx = gL->vert_bot_idxs[i];
	gL->height[idx] = gL->height_gather[i];
	gL->verts[idx][3] = gL->height[idx];

    	gL->bottom[idx] = gL->bottom_gather[i];
	gL->verts[idx][2] = gL->bottom[idx];

	gL->height[idx] -= gL->bottom[idx]; /* height = surf - base */

    	gL->topg[idx] = gL->topg_gather[i];
    }    

    /* debug */
    if (0 && g->rank == 0) {
	FILE *fp = fopen("./output/H0_.m", "w");
	fprintf(fp, "H0 = [\n");
	for (i = 0; i < gL->nvert; i++) {
	    fprintf(fp, "%e\n", gL->height[i]);
	}
	fprintf(fp, "];\n");
	fclose(fp);
    }
    if (0 && g->rank == 0) {
	FILE *fp = fopen("./output/B0_.m", "w");
	fprintf(fp, "H0 = [\n");
	for (i = 0; i < gL->nvert; i++) {
	    fprintf(fp, "%e\n", gL->bottom[i]);
	}
	fprintf(fp, "];\n");
	fclose(fp);
    }


#if 0    
    /* Mesh remap has been done in ice_grid(). */
    phgPrintf("   Re Build layers using height ratio ?\n");
    const FLOAT *ratio = gL->layer_ratio;
    int ii;
    for (ii = 0; ii < gL->nvert_bot; ii++) {
	i = gL->vert_bot_Lidx[ii];
	assert(gL->vert_local_lists[i] != NULL);

	int nv = gL->vert_local_lists[i][0];
	int *iL = &gL->vert_local_lists[i][1];
	int iv = gL->vert_L2S[i];
	assert(nv > 0);

	FLOAT h0, h1;
	h0 = gL->bottom[iv];
	h1 = gL->height[iv] + h0;

	FLOAT H[nv];
	get_layer_height(H, nv, ratio, h0, h1);
		
	assert(gL->max_nlayer + 1 == nv);
	for (j = 0; j < nv; j++) {
	    g->verts[iL[j]][Z_DIR] = H[j];
	}
    }

    phgGeomInit_(g, TRUE);
#endif
    
    //phgExportVTK(g, OUTPUT_DIR "maped.vtk", NULL);
}

#include <metis.h>
#include <parmetis.h>

#if PARMETIS_MAJOR_VERSION == 4
# define ID_TYPE idx_t
#else
# define ID_TYPE idxtype
#endif

static ID_TYPE *idx_array;

void
part_layered_mesh(GRID *g, LAYERED_MESH *gL)
{
    /* Fix me: part on the fly */
    assert(g->nleaf == g->nleaf_global);
    if (phgRank > 0)
	return;

    int ne = gL->ntria, nn = gL->nvert, etype = 1; /* triangle */
    int numflag = 0, edgecut = 0;		   /* c-style */
    ID_TYPE *epart = calloc(ne, sizeof(ID_TYPE));
    ID_TYPE *npart = calloc(nn, sizeof(ID_TYPE));
    int i, j;
    static int *part = NULL; 

#if PARMETIS_MAJOR_VERSION == 4
    ID_TYPE nprocs = phgNProcs;
    ID_TYPE *eptr = calloc(ne+1, sizeof(ID_TYPE));
    for (i = 0; i < ne+1; i++) {
	eptr[i] = 3 * i;
    }
#else
    int nprocs = phgNProcs;
#endif

    if (nprocs <= 1) {
        part = calloc(ne, sizeof(*part));
	return;
    }

    ID_TYPE *elmnts = calloc(3*ne, sizeof(ID_TYPE));
    idx_array = elmnts;
    /* bottom faces */
    for (i = 0; i < ne; i++) 
	for (j = 0; j < 3; j++) 
	    *(idx_array++) = gL->trias[i].verts[j];


#if 0
    /* Metis adaptive */
    phgPrintf("Use metis partitaion.\n");
    ParMETIS_V3_AdaptiveRepart(vtxdist, xadj, adjncy, vwgt, vsize, adjwgt,
			       &wgtflag, &numflag, &ncon, &nparts, tpwgts, ubvec, &ipc2redist, 
			       options, &edgecut, part, &comm);
#elif 1
    phgPrintf("Use metis daul partitaion.\n");
    /* Metis dual */
#  if PARMETIS_MAJOR_VERSION == 3
    METIS_PartMeshDual(&ne, &nn, elmnts, &etype, &numflag, &nprocs,
    		       &edgecut, epart, npart);
#  elif PARMETIS_MAJOR_VERSION == 4
    idx_t ncommon = 2, objval;
    METIS_PartMeshDual(&ne, &nn, eptr, elmnts,
		       NULL, NULL, &ncommon, &nprocs,
		       NULL, NULL, &objval, epart, npart);
#  endif
#else
    /* Read in part info */
    phgPrintf("Use user partitaion.\n");
    FILE *fp = fopen("tri-part.dat", "r");
    int ne0;
    fread(&ne0, sizeof(int), 1, fp);
    assert(ne0 == ne);
    fread(epart, sizeof(*epart), ne, fp);
    fclose(fp);
#endif

    /* statistics */
    int *nep, max_ne, min_ne;
    part = calloc(ne, sizeof(*part));
    nep = calloc(nprocs, sizeof(*nep));
    //printf("\n --part: \n");
    for (i = 0; i < ne; i++) {
	part[i] = epart[i];
	//printf("   elem[%5d]: %2d\n", i, part[i]);
	nep[part[i]]++;
    }
    
    max_ne = -1e6;
    min_ne = 1e6;
    for (i = 0; i < nprocs; i++) {
	if (max_ne < nep[i])
	    max_ne = nep[i];
	if (min_ne > nep[i])
	    min_ne = nep[i];
    }
    phgPrintf("Ne per rank: [%d %d]\n", min_ne, max_ne);

#if 0
    /* output */
    FILE *fp = fopen("tri-part.dat", "w");
    fwrite(&ne, sizeof(ne), 1, fp);
    fwrite(part, sizeof(*part), ne, fp);
    fclose(fp);
#endif


    /* Bring to 3D */
    assert(gL->nface_bot == gL->ntria);
    BOTTOM_FACE *fb;
    fb = gL->face_bot;
    for (i = 0; i < gL->nface_bot; i++, fb++) {
	for (j = 0; j < fb->ne; j++) {
	    SIMPLEX *e = fb->elems[j];
	    e->mark = part[fb->tria];
	}
    }

    
    free(part);
    free(npart);
    free(epart);
    free(elmnts);
    return;
}


void 
destory_layerd_mesh(LAYERED_MESH **gL_ptr)
{
    LAYERED_MESH *gL = *gL_ptr;
    GRID *g = gL->g;
    SIMPLEX ***e;
    int i, **p;

    phgInfo(0, "GRID nvert: %d\n", g->nvert);

    phgFree(gL->verts);
    phgFree(gL->trias);

    for (i = 0; i < gL->nvert; i++)
	phgFree(gL->vert_global_lists[i]);
    phgFree(gL->vert_global_lists);

    if ((p = gL->vert_local_lists) != NULL) {
	assert((e = gL->vert_elements) != NULL);
	for (i = 0; i < g->nvert; i++, p++, e++)
	    if (g->types_vert[i] & BC_BOTTOM) {
		phgFree(*p);
		phgFree(*e);
	    }
	phgFree(gL->vert_local_lists);
	phgFree(gL->vert_elements);
    }

    phgFree(gL->vert_bot_cnts);
    phgFree(gL->vert_bot_dsps);
    phgFree(gL->vert_bot_Lidx);
    phgFree(gL->vert_bot_Gidx);
    phgFree(gL->vert_bot_idxs);
    phgFree(gL->vert_L2S);

    for (i = 0; i < gL->nface_bot; i++)
	phgFree(gL->face_bot[i].elems);
    phgFree(gL->face_bot);

    for (i = 0; i < g->nvert; i++)
	phgFree(gL->vef_3d[i].v);
    phgFree(gL->vef_3d);
    phgFree(gL->vef_2d);
    phgFree(gL);

    *gL_ptr = NULL;
    return;
}







//#define SURF_DOF_TYPE DOF_P2
//#define SURF_DOF_TYPE (ns_params->mapping_type)
#define SURF_DOF_TYPE (ns_params->utype)

DOF *
get_bottom_normal(NSSolver *ns)
/*
 * Compute lower surface grad,
 * type: 0. lower
 *       1. upper
 *
 * (To be used in momentum equ)
 * To be used in slip bdry.
 *
 * */
{
    GRID *g = ns->g;
    SIMPLEX *e;
    LAYERED_MESH *gL = ns->gL;
    int i, j, k, ii, q, ib, s;
    DOF *dof_surf, *dof_gS, *dof_sx;;
    DOF *dof_normal;
    BOTTOM_FACE *fb;
    


#if 0
    /* Analytic base normal */
#  warning     ------------ Analytic base normal -------------
    phgUseIsop = TRUE;		/* Note: Isop Normal !!! */
    dof_normal = phgDofNew(g, SURF_DOF_TYPE, Dim, "dof_normal", DofNoAction);
    phgDofSetDataByFunction(dof_normal, func_normal);
    phgUseIsop = FALSE;

    if (ns->normal == NULL)
	ns->normal = dof_normal;
    else
	phgDofCopy(dof_normal, &ns->normal, NULL, NULL);

    //phgDofDump(dof_normal);
    return dof_normal;
#else    

    phgPrintf("* Compute surf grad with P1\n");

    dof_surf = phgDofNew(g, DOF_DG1, 1, "dof_surf", DofNoAction);
    
    /* for (i = 0; i < g->nvert; i++) { */
    /* 	if (g->types_vert[i] & BC_BOTTOM) { */
    /* 	    assert(gL->vert_local_lists[i] != NULL); */

    /* 	    FLOAT ztop, zbot; */
    /* 	    int nv = gL->vert_local_lists[i][0]; */
    /* 	    int *iL = &gL->vert_local_lists[i][1]; */

    /* 	    ztop = g->verts[iL[nv-1]][Z_DIR]; */
    /* 	    zbot = g->verts[iL[0]][Z_DIR]; */
    /* 	    for (j = 0; j < nv; j++) { */
    /* 		const FLOAT rho = RHO_ICE; */
    /* 		const FLOAT grav = GRAVITY; */
    /* 		FLOAT z = g->verts[iL[j]][Z_DIR]; */

    /* 		*DofVertexData(dof_surf, iL[j]) = zbot; */
    /* 	    } */
    /* 	} */
    /* } */

    phgInfo(0, "Porj grad bas\n");
    
    fb = gL->face_bot;
    //e = NULL;
    for (i = 0; i < gL->nface_bot; i++, fb++) {
	TRIA *t = gL->trias + fb->tria;
	int *vL[3], *iL[3], nL;
	FLOAT zbot[3];

	for (k = 0; k < 3; k++) {
	    assert(g->types_vert[fb->vert[k]] & BC_BOTTOM);
	    vL[k] = gL->vert_local_lists[fb->vert[k]];
	    assert(vL[k] != NULL);

	    nL = vL[k][0];
	    iL[k] = &vL[k][1];
		
	    zbot[k] = g->verts[iL[k][0]][Z_DIR];
	}
	
	for (ii = 0; ii < fb->ne; ii++) { /* tets */
	    SIMPLEX *e = fb->elems[ii];
	    FLOAT *vdat = DofElementData(dof_surf, e->index);

	    for (j = 0; j < NVert; j++) {
		for (k = 0; k < 3; k++)
		    if (e->verts[j] >= iL[k][0]
			&& e->verts[j] <= iL[k][nL-1]) {
			vdat[j] = zbot[k];
			break;
		    }
		assert(k < 3);
	    }
		
	}
    }


    
    /* 
     *
     * Proj grad surf to Dofs
     *
     * */
    SOLVER *solver;
    VEC *vec[2];
    
    dof_gS = phgDofGradient(dof_surf, NULL, NULL, "dof_gS");
    dof_sx = phgDofNew(g, SURF_DOF_TYPE, 1, "dof_sx", DofNoAction);
    dof_normal = phgDofNew(g, SURF_DOF_TYPE, Dim, "dof_bottom_normal", DofNoAction);

    phgOptionsPush();
    phgOptionsSetOptions(ns_params->proj_opts);
    solver = phgSolverCreate(SOLVER_DEFAULT, dof_sx, NULL);
    solver->verb = 0;
    phgOptionsPop();

    
    for (k = 0; k < 2; k++) {
	vec[k] = phgMapCreateVec(solver->rhs->map, 1) ;
	phgVecDisassemble(vec[k]);
    }

    fb = gL->face_bot;
    for (ib = 0; ib < gL->nface_bot; ib++, fb++) {
	TRIA *t = gL->trias + fb->tria;
	FLOAT area_tri = t->area;
	
	for (ii = 0; ii < fb->ne; ii++) { /* tets */
	    SIMPLEX *e = fb->elems[ii];

	    int M = dof_sx->type->nbas;	/* num of bases of Velocity */
	    int order = DofTypeOrder(dof_sx, e) * 2;
	    FLOAT A[M][M], rhs[2][M];
	    INT I_[M];
	    QUAD *quad;
	    FLOAT vol, det, area;
	    const FLOAT *w, *p, *gS;

	    
	    Bzero(A); Bzero(rhs);

#if 0
	    /* ------------------------------
	     *
	     *   Vol proj
	     * 
	     * ------------------------------ */
	    quad = phgQuadGetQuad3D(order);
	    gS = phgQuadGetDofValues(e, dof_gS, quad); 
	    vol = area_tri;

	    p = quad->points;
	    w = quad->weights;
	    for (q = 0; q < quad->npoints; q++) {
		phgGeomGetCurvedJacobianAtLambda(g, e, p, &det);

		for (i = 0; i < M; i++) {
		    const FLOAT *gi = phgQuadGetBasisValues(e, dof_sx, i, quad) + q;    
		    for (j = 0; j < M; j++) {
			const FLOAT *gj = phgQuadGetBasisValues(e, dof_sx, j, quad) + q;       
			FLOAT qmass = vol*(*w) * (*gj) * (*gi);
			A[i][j] += qmass;
		    }
		    
		    for (k = 0; k < 2; k++) {
			rhs[k][i] += vol*(*w) * gS[k] * (*gi); 
		    }
		}
		gS += Dim;
		w++; p += Dim + 1;
	    }
#else
	    /* ------------------------------
	     *
	     *   Face proj
	     * 
	     * ------------------------------ */
	    for (s = 0; s < NFace; s++) {
		if (e->bound_type[s] & BC_BOTTOM) {
		    int v0, v1, v2;
		    int nbas_face = NbasFace(dof_sx);
		    SHORT bases[nbas_face];
		    FLOAT lambda[Dim + 1], x, y, z, gS[Dim];
		    order = 2 * DofTypeOrder(dof_sx, e);

		    phgDofGetBasesOnFace(dof_sx, e, s, bases);
		    v0 = GetFaceVertex(s, 0);
		    v1 = GetFaceVertex(s, 1);
		    v2 = GetFaceVertex(s, 2);
		    lambda[s] = 0.;

/* 		    SHOW_iV_(0, bases, nbas_face); */

#if 0
#   warning --- Sx proj on manifold -----
		    /* Proj on manifold */
		    area = phgGeomGetFaceArea(g, e, s);
#else		    
#   warning --- Sx proj on plane -----
		    /* Proj on plane */
		    area = area_tri;
#endif		    
		    quad = phgQuadGetQuad2D(order);

		    
		    p = quad->points;
		    w = quad->weights;
		    for (q = 0; q < quad->npoints; q++) {
			lambda[v0] = *(p++);
			lambda[v1] = *(p++);
			lambda[v2] = *(p++);

			phgDofEval(dof_gS, e, lambda, gS);

			const FLOAT *gi = 
			    dof_sx->type->BasFuncs(dof_sx, e, 0, -1, lambda);
			
/* 			SHOW_V(gi, M); */
			
			for (i = 0; i < nbas_face; i++) {
			    int ii = bases[i];
			    for (j = 0; j < nbas_face; j++) { 
				int jj = bases[j];
				
				FLOAT mass_face = area*(*w) * (gi[jj])*(gi[ii]);
				A[ii][jj] += mass_face;
			    }
		    
			    for (k = 0; k < 2; k++) {
				rhs[k][ii] += area*(*w) * gS[k] * (gi[ii]); 
			    }
			}

			phgInfo(3, "w: %30.15e g: %30.15e\n", *w, gi[bases[0]]);
			
			w++;
		    } /* end qp */

/* 		    SHOW_M(&A[0][0], M, M); */
/* 		    SHOW_M(&rhs[0][0], 2, M); */
		    
		}     /* end face type*/
	    }	      /* end face */

	    
	    for (i = 0; i < M; i++) {
		BTYPE btype =
		    phgDofGetElementBasisInfo(dof_sx, e, i,
					      NULL, NULL, NULL);
		if (btype & BC_BOTTOM)
		    continue;

		A[i][i] += 1;
	    }
#endif	    

	    /* Map: Element -> system */
	    for (i = 0; i < M; i++)
		I_[i] = phgMapE2L(solver->mat->cmap, 0, e, i);

	    /* Global res */
	    for (i = 0; i < M; i++)
		phgMatAddEntries(solver->mat, 1, I_ + i, M, I_,
				 &(A[i][0])); 

	    for (k = 0; k < 2; k++)
		phgVecAddEntries(vec[k], 0, M, I_, &rhs[k][0]);
	}
    }/* end element */
    

    for (k = 0; k < 2; k++)
	phgVecAssemble(vec[k]);
    solver->rhs_updated = TRUE;
    
    INT n = DofGetDataCount(dof_sx);
    for (k = 0; k < 2; k++) {
	phgVecCopy(vec[k], &solver->rhs);
	phgDofSetDataByValue(dof_sx, 0.);
	phgSolverSolve(solver, FALSE, dof_sx, NULL);
	phgPrintf("      solver_Sx: nits = %d, resid = %0.4lg ",
		  solver->nits, solver->residual);
	elapsed_time(g, TRUE, phgPerfGetMflops(g, NULL, NULL));

	FLOAT *vg = dof_normal->data, *vg0 = dof_sx->data;
	for (i = 0; i < n; i++, vg0++, vg += Dim)
	    vg[k] = *vg0;
    }
     
    for (k = 0; k < 2; k++)
	phgVecDestroy(&vec[k]);
    phgSolverDestroy(&solver);


    FLOAT *vg = dof_normal->data;
    FLOAT norm;
    for (i = 0; i < n; i++, vg += Dim) {
	vg[Z_DIR] = -1.;
	norm = sqrt(INNER_PRODUCT(vg, vg));
	vg[0] /= norm;
	vg[1] /= norm;
	vg[2] /= norm;
    }

    //phgDofDump(dof_normal);
    //phgExportVTK(g, "gradB.vtk", dof_surf, dof_normal, NULL);
    
    phgDofFree(&dof_gS);
    phgDofFree(&dof_surf);
    phgDofFree(&dof_sx);

    if (ns->bottom_normal == NULL)
	ns->bottom_normal = dof_normal;
    else
	phgDofCopy(dof_normal, &ns->bottom_normal, NULL, NULL);
    
    return dof_normal;
#endif    
}



DOF *
get_lateral_normal(NSSolver *ns)
{
    GRID *g = ns->g;
    LAYERED_MESH *gL = ns->gL;
    INT i, j;

#if 1
    /*
     * The lateral normal is computed using 2D grid info, and then carried to 3D 
     *
     * */
    
    DOF *dof_normal_P1, *dof_normal = NULL;
    dof_normal_P1 = phgDofNew(g, DOF_P1, Dim, "dof_lateral_normal", DofNoAction);
    
    for (i = 0; i < g->nvert; i++) {
	if (g->types_vert[i] & BC_BOTTOM
	    && g->types_vert[i] & BC_LATERAL) {
	    assert(gL->vert_local_lists[i] != NULL);

	    int nv = gL->vert_local_lists[i][0];
	    int *iL = &gL->vert_local_lists[i][1];

	    const FLOAT *normal = gL->normals[gL->vert_L2S[i]];
#define NORM2(x) sqrt(x[0]*x[0] + x[1]*x[1])
	    assert(NORM2(normal) > 0);
#undef NORM2
	    
	    for (j = 0; j < nv; j++) {
		FLOAT *vd = DofVertexData(dof_normal_P1, iL[j]);
		vd[X_DIR] = normal[X_DIR];
		vd[Y_DIR] = normal[Y_DIR];
		vd[Z_DIR] = 0.;
	    }
	}
    }


    /* Copy P1 to P2(or P1) */
    dof_normal = phgDofCopy(dof_normal_P1, NULL,
			    SURF_DOF_TYPE, "dof_lateral_normal");

    FLOAT *vd =	dof_normal->data;
    for (i = 0; i < DofGetDataCount(dof_normal) / Dim; i++) {
	FLOAT nn = sqrt(INNER_PRODUCT(vd, vd));
	if (nn > 0) {
	    vd[X_DIR] /= nn;
	    vd[Y_DIR] /= nn;
	    vd[Z_DIR] /= nn;
	}
	vd += Dim;	
    }
    phgDofFree(&dof_normal_P1);
#else
    /*
     * Use 3D normal average
     *  */
    DOF *dof_normal = phgDofNew(g, SURF_DOF_TYPE, Dim, "dof_lateral_normal", DofNoAction);
    phgDofSetDataByValue(dof_normal, 0.);

    ForAllElements(g, e) {
	for (s = 0; s < NFace; s++) 
	    if (e->bound_type[s] & BC_LATERL) {
		int v[3] = {
		    GetFaceVertex(s, 0),
		    GetFaceVertex(s, 1),
		    GetFaceVertex(s, 2)
		};
		
		normal = phgGeomGetFaceOutNormal(g, e, s);

#  if 0		
		/* no weights, accumulate ? */
		for (k = 0; k < 3; k++) {
		    dof_normal->data[e->verts[v[k]] * Dim + 0] += normal[0];
		    dof_normal->data[e->verts[v[k]] * Dim + 1] += normal[1];
		    dof_normal->data[e->verts[v[k]] * Dim + 2] += normal[2];
		}
#  else
		/* weights by angle */
		
#  endif		
	    }
    }

    /* global average */
    
#endif    


    
    if (ns->lateral_normal == NULL)
	ns->lateral_normal = dof_normal;
    else
	phgDofCopy(dof_normal, &ns->lateral_normal, NULL, NULL);

   
    return dof_normal;
}
