#ifndef PRISM_ELEM_H
#include "phg.h"





#define NVERT_PRISM 6
#define NFACE_PRISM 5
#define NDOF_P2_PRISM 18
#define NDOF_MINI_PRISM (6+3)

//#define NDOF_P0_PRISM 1
//#define NBAS_PRISM 6



typedef struct PRISM_ELEM_ {
    INT verts[NVERT_PRISM];
    //SIMPLEX *tets[3];
    INT tetIds[3];		/* mapping to tet */
    INT *dofIds[4];		/* mapping to tet, P1 P2 MINI  */
    INT *dofId2s[4];		/* mapping to tet, P1 P2 MINI, non period */
    BTYPE bound_type[NFACE_PRISM];
    COORD face_normal[NFACE_PRISM];
    FLOAT face_area[NFACE_PRISM];
    FLOAT face_diam[NFACE_PRISM];
} PRISM_ELEM;
#define PRISM_ELEM_H
#endif
