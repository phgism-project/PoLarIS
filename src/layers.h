#ifndef LAYERED_MESH_H
#include "phg.h"


typedef struct TRIA_ {
    int verts[3];
    FLOAT area;
    FLOAT normals[2*3];
    FLOAT dual_edge[3];
    BYTE region_mark;
    BYTE bound_types[3];
} TRIA;


/* vert to tria (2D) */
typedef struct VEFMAP_2D_ {
    int ne;
    int idx[12];		/* max 12 elem */
} VEFMAP_2D; 

/* vert to vert (3D) */
typedef struct VEFMAP_3D_ {
    int nv;
    int size;
    int *v;			/* var length */
} VEFMAP_3D; 


typedef struct BOTTOM_FACE_ {
    SIMPLEX *e;			/* element  */
    int face;			/* element[face] */
    int vert[3];		/* local vert */
    int tria;			/* Global tria index */
    
    /* List of tets */
    int ne;			/* # of tets in Z dir */
    SIMPLEX **elems;		/* tets in Z dir */
    int *faces;			/* tria face[3] as a,b,c,
				 * faces[ne][3]
				 *   = {x, x, x;
				 *      x, x, x;}
				 *      a  b  c
				 * */
    int *edge2to3;		/* map 2D edge to 3D edges.
				 * [ne][3][2]
				 * */
} BOTTOM_FACE;

typedef struct LAYERED_MESH_ {
    int nvert;			/* # of verts (2D) */
    int ntria;			/* # of trias (2D) */
    int max_nlayer;		/* max # of layers */
    FLOAT *layer_ratio;	/* layer height ratio [0:...:1] */

    /* 2D */
    FLOAT (*verts)[4];		/* verts coord,
				 * (x, y, z_base, z_surf)
				 * */
    int *vert_type;		/* verts type */
    TRIA *trias;		/* tira to verts (2D) */
    int **vert_global_lists;	/* verts (2D) to vert chain (3D)
				 * [nv_2d][nlayer_i+1],
				 * Note: First num nlayer_i */
    int **vert_local_lists;	/* [nL],    [nL_2d][nlayer+1] */
    SIMPLEX ***vert_elements;	/* [nL],    [nL_2d][nlayer] */

    FLOAT *volumes;
    FLOAT *height;			/* [nv_2d] */
    FLOAT *height_local;		/* [nvert_bot] */
    FLOAT *height_gather;		/* [nvert_bot] */

    FLOAT *topg;			/* [nv_2d] */
    FLOAT *topg_local;			/* [nv_bot] */
    FLOAT *topg_gather;			/* [nv_bot] */

    FLOAT *bottom;			/* [nv_2d] */
    FLOAT *bottom_local;		/* [nv_bot] */
    FLOAT *bottom_gather;		/* [nv_bot] */

    FLOAT (*normals)[2];		/* [nv_2d], x-y plane out normal */
    
    /* 3D */
    GRID *g;
    int nvert_bot;		/* # of local verts on bottom */
    int *vert_bot_Lidx;		/* local verts on bottom */
    int *vert_bot_Gidx;		/* global verts on bottom */

    int vert_bot_total;
    int *vert_bot_cnts;
    int *vert_bot_dsps;		
    int *vert_bot_idxs;		

    int *vert_L2S;		/* map vert Local to Bottom,
				 * [g->nvert] -> [nv_2D] */
    int *vert_S2L;		/* map Bottom to vert Local,
				 * [nv_2D] -> [g->nvert] */

    int nface_bot;		/* # of local tria on bottom */
    BOTTOM_FACE *face_bot;

    /* Output */
    VEFMAP_2D *vef_2d;		/* (2D): vert to tria connection
				 * [nv_2d][...]*/
    int *mapE3to2;		/* face (3D) to tria (2D)
				 * [ne_2d] */
    VEFMAP_3D *vef_3d;		/* vert to vert (3D)
				 * [g->nvert] */

    /* Extra */
    void *fv_data;
    
    MPI_Comm comm;
} LAYERED_MESH;

typedef struct MAT_FACTOR_ {
    INT *line_poc0;
    INT *line_pocs;
    INT *line_piv0;
    INT *line_pivs;
    INT *line_mat0;
    FLOAT *line_matv;
    BOOLEAN factorized;
} MAT_FACTOR;
 
/* block dof smoother */
typedef struct MG_BLOCK_DOFS_ {
    int nds[10];
    DOF_TYPE *dof_type;

    MAP *map_cnnt;
    MAT *mat_cnnt;
    INT nline;
    INT *line_dofs;
    INT *line_dof0;
    INT *line_ndof;

    int neigh_size;
 
    MAT *mat;
    MAT *mat_fix;
    MAT_FACTOR *mat_factor;
    MAT_FACTOR *mat_fix_factor;

    BOOLEAN constrained;
} MG_BLOCK_DOFS ;


/* Main functions */
LAYERED_MESH *import_layered_mesh(const char *file_name);
void build_layered_mesh(GRID *g, LAYERED_MESH *gL);
void build_layered_mesh_height(GRID *g, LAYERED_MESH *gL);
void part_layered_mesh(GRID *g, LAYERED_MESH *gL);
void destory_layerd_mesh(LAYERED_MESH **gL_ptr);
void destroy_line_block(MG_BLOCK_DOFS **bk_ptr);
void check_height(GRID *g, LAYERED_MESH *gL);

#define check_height_(gL) {				\
	phgPrintf("   # file: %15s, line: %d\n",	\
		  __FILE__, __LINE__);			\
	check_height(g, gL);				\
    }

MG_BLOCK_DOFS *init_line_block(DOF *u0, LAYERED_MESH *gL);
void factorize_line_block(GRID *g, MG_BLOCK_DOFS *bk, MAT *A, MAT_FACTOR **mf_ptr);
void free_fact_line_block(MAT_FACTOR **mf_ptr);
void line_smooth(MG_BLOCK_DOFS *bk, MAT *A, VEC *x, VEC *b,
		 int nsmooth, void *ctx);

/* Lapack routines */
#define integer int
#define doublereal double
int dgbtrf_(integer *m, integer *n, integer *kl, integer *ku,
	    doublereal *ab, integer *ldab, integer *ipiv, integer *info);

int dgbtrs_(char *trans, integer *n, integer *kl, integer *
	    ku, integer *nrhs, doublereal *ab, integer *ldab, integer *ipiv,
	    doublereal *b, integer *ldb, integer *info);
#undef integer
#undef doublereal


#define LAYERED_MESH_H
#endif
