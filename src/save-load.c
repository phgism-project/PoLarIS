/*
 * Save & load DOF data, only for nodal DOF type.
 *
 * Usage:
 *
 * 1. Save: 
 * 1.1 phgExportMedit2: save mesh, vert coord,
 * 1.2 save_map       : save elem to vert edge face map
 * 1.3 save_dof_data  : save dof data
 *
 * 2. Load:
 * 2.0 Use exteranl program to generate text mesh file with binary
 *     coord file. 
 * 2.1 save_element_id: save element index, before redistribute
 * 2.2 load_dof_data2 : load dof date, after redistribute
 * 2.3 free_element_id
 *
 * Note:
 * Load is essentially a sequtial process, the time is mostly spend
 * on distribute mesh.
 * Using load_dof_data2, no dof data distribute is not needed, still,
 * further improvement could be made to read in only locate map and 
 * dof data. 
 *
 *  */

#include "phg.h"
#include "phg/io.h"
#include <string.h>
#include <math.h>

#define DofVertData(dof, i) DofVertexData(dof, i)
#define GlobalVert(g,no) GlobalVertex(g, no)
//#define GlobalElem(g,no) GlobalElement(g, no)
#define NElem 0

#define GetElementVertGMap(e, M) {					\
	int _v;								\
	GetElementVertices(e, 0, M[0], M[1], M[2], M[3], _v, _v, _v, _v); \
    }
#define GetElementFaceGMap(e, M) {					\
	int _v;								\
	GetElementVertices(e, 0, M[0], M[1], M[2], M[3], _v, _v, _v, _v); \
    }

static DOF *elem_id = NULL;

static void 
GetElementEdgeGMap(SIMPLEX *e, int *M)
{
    int i, v[2], N[NVert]; 
    
    GetElementVertGMap(e, N);
    bzero(M, NEdge * sizeof(*M));
    for (i = 0; i < NEdge; i++) {
	v[0] = GetEdgeVertex(i, 0);
	v[1] = GetEdgeVertex(i, 1);
	v[0] = N[v[0]];
	v[1] = N[v[1]];
	M[i] = GetEdgeNo(v[0], v[1]);
    }
    //SHOW_V(M, NFace);
    return;
}


void 
save_map(GRID *g, const char *file) {
    SIMPLEX *e;
    int i, j, n, nEmap = NVert + NEdge + NFace + NElem;
    FILE *fp = NULL;
    MPI_Status status;
    //FLOAT *fbuffer;
    INT *ibuffer;
# define TRUNK_SIZE 65536

    /* save E2Gmap */
    if (g->rank == 0) {
	char e2g_file[1000];
	sprintf(e2g_file, "%s.map", file);
	if ((fp = fopen(e2g_file, "w+t")) == NULL) {
	    phgError(1, "cannot open output file \"%s\".\n", e2g_file);
	    return;
	}

	fwrite(&g->nelem_global, sizeof(INT), 1, fp);
	ForAllElements(g, e) {
	    int M[6];
#define WRITE_MAP(GType, gtype)					\
	    GetElement##GType##GMap(e, M);			\
	    for (i = 0; i < N##GType; i++) {			\
		INT ind = Global##GType(g, e->gtype##s[M[i]]);	\
		fwrite(&ind, sizeof(INT), 1, fp);		\
	    }

	    WRITE_MAP(Vert, vert);
	    WRITE_MAP(Edge, edge);
	    WRITE_MAP(Face, face);
	    //WRITE_MAP(Elem, elem);
#undef WRITE_MAP
	}
	/* receive elements from other processes */
	ibuffer = phgAlloc(TRUNK_SIZE * (nEmap) * sizeof(*ibuffer));
	for (j = 1; j < g->nprocs; j++) {
	    while (TRUE) {
		MPI_Probe(j, 222, g->comm, &status);
		MPI_Get_count(&status, PHG_MPI_INT, &n);
		assert(n <= ((nEmap) * TRUNK_SIZE) && n % (nEmap) == 0);
		MPI_Recv(ibuffer, n, PHG_MPI_INT, j, 222, g->comm, &status);
		/* process received vertices */
		for (i = 0; i < n; i += nEmap) {
		    fwrite(ibuffer + i, sizeof(INT), nEmap, fp);
		}
		if (n < (nEmap) * TRUNK_SIZE)
		    break;
	    }
	}
	phgFree(ibuffer);
	fclose(fp);
    }
    else {
	/* send elements to root process */
	ibuffer = phgAlloc(TRUNK_SIZE * (nEmap) * sizeof(*ibuffer));
	n = 0;
	ForAllElements(g, e) {
	    int k = 0, M[6];
#define WRITE_MAP(GType, gtype)					\
	    GetElement##GType##GMap(e, M);			\
	    for (i = 0; i < N##GType; i++) {			\
		ibuffer[(nEmap) * n + k++] =			\
		    Global##GType(g, e->gtype##s[M[i]]);	\
	    }

	    WRITE_MAP(Vert, vert);
	    WRITE_MAP(Edge, edge);
	    WRITE_MAP(Face, face);
	    //WRITE_MAP(Elem, elem);
#undef WRITE_MAP
	    if (++n >= TRUNK_SIZE) {
		MPI_Send(ibuffer, (nEmap) * n, PHG_MPI_INT, 0, 222, g->comm);
		n = 0;
	    }
	}
	/* send the last block (may be size 0), which also marks EOD */
	MPI_Send(ibuffer, (nEmap) * n, PHG_MPI_INT, 0, 222, g->comm);
	phgFree(ibuffer);
    }

    return;
}

void 
save_dof_data(GRID *g, DOF *dof, const char *file)
{
    FLOAT *buffer = NULL;
    size_t buffer_size = 0;
    INT ndata;

    if (phgIOOpen(g, file) == FALSE)
	phgError(1, "Export Dof data file %s failed!", file);
    buffer_size += DofGetDataCount(dof);
    buffer = phgCalloc(buffer_size,  sizeof(FLOAT));

#define WRITE_DATA(dof, GType, gtype)					\
    if (dof->type->np_##gtype > 0){					\
	ndata = g->n##gtype * dof->count_##gtype;			\
	memcpy(buffer, dof->data_##gtype, ndata * sizeof(FLOAT));	\
	phgIOWrite(g, buffer, sizeof(FLOAT) * dof->count_##gtype,	\
		   g->n##gtype, g->n##gtype##_global,			\
		   g->types_##gtype, g->L2Gmap_##gtype, TRUE);		\
    }

    WRITE_DATA(dof, Vertex, vert);
    WRITE_DATA(dof, Edge, edge);
    WRITE_DATA(dof, Face, face);
    //WRITE_DATA(dof, Elem, elem);
#undef WRITE_DATA
    phgFree(buffer);
    phgIOClose();

    return;
}


void 
load_dof_data(GRID *g, DOF *dof, const char *data_file, const char *mesh_file)
/* Load dof data BEFORE mesh redistribute. The loading process is sequential.
 *
 * E2Gmap is read element wize only on root proc.
 *  */
{
    SIMPLEX *e;
    FLOAT *buffer = NULL, *buf_vert, *buf_edge,	*buf_face, *buf_elem;
    int dim = dof->dim;
    size_t buffer_size = 0;
    int i, j, k, n, count, nEmap = NVert + NEdge + NFace + NElem;
    INT ndata, nelem, *E2Gmap = NULL;
    FILE *fp;
    char e2g_file[1000];

    if (phgRank == 0) {
	if ((fp = fopen(data_file, "r")) == NULL)
	    phgError(1, "read Dof data %s failed!\n", data_file);
	buffer_size += DofGetDataCount(dof);
	buffer = phgCalloc(buffer_size, sizeof(FLOAT));
	buf_vert = buffer;
	buf_edge = buf_vert + g->nvert * dof->count_vert;
	buf_face = buf_edge + g->nedge * dof->count_edge;
	buf_elem = buf_face + g->nface * dof->count_face;

#define READ_DATA(dof, GType, gtype)				\
	if (dof->type->np_##gtype > 0){				\
	    int n##gtype = g->n##gtype;			\
	    ndata = n##gtype * dof->count_##gtype;		\
	    n = fread(buf_##gtype, ndata, sizeof(FLOAT), fp);	\
	}

	READ_DATA(dof, Vertex, vert);
	READ_DATA(dof, Edge, edge);
	READ_DATA(dof, Face, face);
	//READ_DATA(dof, Elem, elem);
#undef READ_DATA
	fclose(fp);

	/* read E2Gmap */
	sprintf(e2g_file, "%s.map", mesh_file);
	if ((fp = fopen(e2g_file, "r+t")) == NULL) {
	    phgError(1, "cannot open output file \"%s\".\n", e2g_file);
	    return;
	}

	E2Gmap = phgCalloc(nEmap, sizeof(*E2Gmap));
	n = fread(&nelem, sizeof(INT), 1, fp);
	assert(g->nelem == nelem);

	/* Note:
	 * element order of E2Gmap is the same as vector g->roots, but NOT
	 * the same as g->elements, as well as ForAllElements.
	 */
	e = g->roots;
	for (j = 0; j < g->nelem; j++, e++) {
	    FLOAT *val;
	    int M[6], shift = 0;
	    
	    n = fread(E2Gmap, sizeof(INT), nEmap, fp);

	    /* Vertex data */
	    GetElementVertGMap(e, M);
	    for (i = 0; i < NVert; i++) {
		val = DofVertexData(dof, GlobalVertex(g, e->verts[M[i]]));
		for (k = 0; k < dim; k++)
		    val[k] = buf_vert[E2Gmap[i]*dim + k];
		//printf("E:%d, vert[%d]:%d\n", e->index, i, e->verts[i]);
	    }

	    /* Edge data */
	    shift += NVert;
	    count = dof->count_edge;
	    GetElementEdgeGMap(e, M);
	    for (i = 0; i < NEdge; i++) {
		val = DofEdgeData(dof, GlobalEdge(g, e->edges[M[i]]));
		for (k = 0; k < count; k++)
		    val[k] = buf_edge[E2Gmap[i+shift]*count + k];
	    }
		
	    /* Face data */
	    shift += NEdge;
	    count = dof->count_face;
	    GetElementFaceGMap(e, M);
	    for (i = 0; i < NFace; i++) {
		val = DofFaceData(dof, GlobalFace(g, e->faces[M[i]]));
		for (k = 0; k < count; k++) 
		    val[k] = buf_face[E2Gmap[i+shift]*count + k];
	    }
	}

	fclose(fp);
	phgFree(buffer);
    }
    return;
}

void 
load_dof_data2(GRID *g, DOF *dof, const char *data_file, const char *mesh_file)
/* Load dof data AFTER mesh redistribute. The loading process is parallel.
 *
 * E2Gmap is read as a whole on each proc.
 * TODO: read only the local part of E2Gmap and local part of dof data.
 *  */
{
    SIMPLEX *e;
    FLOAT *buffer = NULL, *buf_vert, *buf_edge,	*buf_face, *buf_elem;
    int dim = dof->dim;
    size_t buffer_size = 0;
    int i, k, n, count, nEmap = NVert + NEdge + NFace + NElem;
    INT ndata, nelem_global, *E2Gmap = NULL;
    FILE *fp;
    char e2g_file[1000];

    if ((fp = fopen(data_file, "r")) == NULL)
	phgError(1, "read Dof data %s failed!\n", data_file);
    buffer_size += DofGetDataCountGlobal(dof);
    buffer = phgCalloc(buffer_size, sizeof(FLOAT));
    buf_vert = buffer;
    buf_edge = buf_vert + g->nvert_global * dof->count_vert;
    buf_face = buf_edge + g->nedge_global * dof->count_edge;
    buf_elem = buf_face + g->nface_global * dof->count_face;
    phgPrintf("* Load dof: load datafile: %0.4lfMB\n", 
	      phgMemoryUsage(g, NULL) / (1024.0 * 1024.0));

#define READ_DATA(dof, GType, gtype)				\
    if (dof->type->np_##gtype > 0){				\
	int n##gtype_global = g->n##gtype##_global;		\
	ndata = n##gtype_global * dof->count_##gtype;		\
	n = fread(buf_##gtype, ndata, sizeof(FLOAT), fp);	\
    }

    READ_DATA(dof, Vertex, vert);
    READ_DATA(dof, Edge, edge);
    READ_DATA(dof, Face, face);
    //READ_DATA(dof, Elem, elem);
#undef READ_DATA
    fclose(fp);

    /* read E2Gmap */
    sprintf(e2g_file, "%s.map", mesh_file);
    if ((fp = fopen(e2g_file, "r+t")) == NULL) {
	phgError(1, "cannot open output file \"%s\".\n", e2g_file);
	return;
    }
    E2Gmap = phgCalloc(nEmap * g->nelem_global, sizeof(*E2Gmap));
    phgPrintf("* Load dof: load E2Gmap: %0.4lfMB\n", 
	      phgMemoryUsage(g, NULL) / (1024.0 * 1024.0));
    n = fread(&nelem_global, sizeof(INT), 1, fp);
    n = fread(E2Gmap, sizeof(INT), nEmap * g->nelem_global, fp);
    assert(g->nelem_global == nelem_global);
    ForAllElements(g, e) {
	FLOAT *val;
	INT *e2g_map, M[6], shift = 0, e_id;

	e_id = (INT)(*DofElementData(elem_id, e->index));
	e2g_map = E2Gmap + e_id * nEmap;

	/* Vertex data */
	GetElementVertGMap(e, M);
	count = dof->count_vert;
	for (i = 0; i < NVert; i++) {
	    val = DofVertexData(dof, e->verts[M[i]]);
	    for (k = 0; k < dim; k++)
		val[k] = buf_vert[e2g_map[i]*count + k];
	    //printf("E:%d, vert[%d]:%d\n", e->index, i, e->verts[i]);
	}

	/* Edge data */
	shift += NVert;
	count = dof->count_edge;
	GetElementEdgeGMap(e, M);
	for (i = 0; i < NEdge; i++) {
	    val = DofEdgeData(dof, e->edges[M[i]]);
	    for (k = 0; k < count; k++)
		val[k] = buf_edge[e2g_map[i+shift]*count + k];
	}
		
	/* Face data */
	shift += NEdge;
	count = dof->count_face;
	GetElementFaceGMap(e, M);
	for (i = 0; i < NFace; i++) {
	    val = DofFaceData(dof, e->faces[M[i]]);
	    for (k = 0; k < count; k++) 
		val[k] = buf_face[e2g_map[i+shift]*count + k];
	}
    }

    fclose(fp);
    phgFree(buffer);
    phgFree(E2Gmap);
    return;
}

void 
save_dof_data3(GRID *g, DOF *dof, const char *file)
/*
 * Each proc save it's own Dof data
 *  */
{
    FILE *fp = NULL;
    char fname[100];
    INT n = DofGetDataCount(dof);

    sprintf(fname, "%s.p%03d", file, g->rank);
    phgInfo(1, "* save data to %s\n", fname);
    if ((fp = fopen(fname, "w")) == NULL) {
	phgError(1, "read Dof data %s failed!\n", fname);
    } else {
	fwrite(&n, sizeof(INT), 1, fp);
	fwrite(dof->data, sizeof(FLOAT), n, fp);
	fclose(fp);
    }

    return;
}

void 
load_dof_data3(GRID *g, DOF *dof, const char *file, const char *mesh_file)
/*
 * Each proc read it's own Dof data
 *  */
{
    FILE *fp = NULL;
    char fname[100];
    INT n0 = DofGetDataCount(dof), n;
    int nr;

    Unused(mesh_file);
    sprintf(fname, "%s.p%03d", file, g->rank);
    phgInfo(1, "* read data from %s\n", fname);
    if ((fp = fopen(fname, "r")) == NULL) {
	phgError(1, "read Dof data %s failed!\n", fname);
    } else {
	nr = fread(&n, sizeof(INT), 1, fp);
	if (n0 != n) 
	    phgError(1, "read [%d] but saved [%d] !\n", n0, n);
	nr = fread(dof->data, sizeof(FLOAT), n, fp);
	fclose(fp);
    }

    return;
}


void
save_element_id(GRID *g)
{
    SIMPLEX *e;
    INT i;
    FLOAT *v;

    elem_id = phgDofNew(g, DOF_P0, 1, "element_id", DofInterpolation);
    if (phgRank == 0) {
	e = g->roots;
	for (i = 0; i < g->nelem; i++, e++) {
	    v = DofElementData(elem_id, e->index);
	    *v = i;
	}
    }
    return;
}

void 
free_element_id(GRID *g)
{
    phgDofFree(&elem_id);
    elem_id = NULL;
}


void
phgResumeStage(GRID *g, FLOAT *time, int *tstep, char *mesh_file, char *data_file)
/* Load resume info (time, tstep, mesh file name, data file name)
 * from file "resume.log".
 *
 * Note: resume should be execute by ALL ranks.
 * */
{

    FILE *fp;
    FLOAT time0;
    int tstep0, n;
    char mesh0[100];
    char data0[100];

    if ((fp = fopen("resume.log", "r")) == NULL)
	phgError(1, "read resume.log failed!\n");
    n = fscanf(fp, "%lf", &time0);
    n = fscanf(fp, "%d", &tstep0);
    n = fscanf(fp, "%s", mesh0);
    n = fscanf(fp, "%s", data0);
    n = fclose(fp);
    if (time != NULL)
	*time = time0;
    if (tstep != NULL)
	*tstep = tstep0;
    if (mesh_file != NULL)
	strcpy(mesh_file, mesh0);
    if (data_file != NULL)
	strcpy(data_file, data0);

    return;
}


/*******************/
/* Resume routines */
/*******************/
void
phgResumeLogUpdate(GRID *g, FLOAT *time, int *tstep, char *mesh_file, char *data_file)
/* Save resume info (time, tstep, mesh file name, data file name)
 * to file "resume.log".
 *
 * Resume info is stored as static variable, component is not updated
 * if input points is NULL.
 * */
{
    static FLOAT time0;
    static int tstep0;
    static char mesh0[100];
    static char data0[100];
    FILE *fp;

    if (g->rank != 0)
	return;

    if (time != NULL)
	time0 = *time;
    if (tstep != NULL)
	tstep0 = *tstep;
    if (mesh_file != NULL)
	strcpy(mesh0, mesh_file);
    if (data_file != NULL)
	strcpy(data0, data_file);

    if ((fp = fopen("resume.log", "w")) == NULL)
	phgError(1, "open resume log file failed!\n");
    fprintf(fp, "%24.12E\n%d\n%s\n%s\n", time0, tstep0, mesh0, data0);
    fclose(fp);

    return;
}




/* --------------------------------------------------------------------------------
 *
 *
 * Save load using pnetcdf
 *
 * 
 * -------------------------------------------------------------------------------- */
#include <pnetcdf.h>
#define NAME_MAX_LENGTH 256


static void
handle_error(int status, int lineno)
{
    fprintf(stderr, "Error at line %d: %s\n", lineno, ncmpi_strerror(status));
    MPI_Abort(MPI_COMM_WORLD, 1);
}

static char *save_prefix = "output/step";
static char *load_prefix = "output/step";


void
save_state(int tstep, double time, int ndof, DOF **dofs, char **dof_names)
//
// Create netcdf files to store the current state:
//   step#tstep.nc:
//   containig: time, dofs
//   ...
//
{
    int nprocs = phgNProcs, rank = phgRank;
    int ret, ncfile, dimid, ndims=1;
    MPI_Offset start, count=1;
    char file_name[NAME_MAX_LENGTH], tmp_str[NAME_MAX_LENGTH];

    MPI_Offset ndat[ndof];
    MPI_Offset cnts[nprocs][ndof];
    MPI_Offset offs[nprocs+1][ndof];
    int cnts_int[ndof][nprocs];	 // copy to int type, and transpose
    int var_id[nprocs], var_part_id[nprocs];

    phgInfo(0, "\n* Save state\n");
    for (int idof = 0; idof < ndof; idof++) {
	ndat[idof] = DofGetDataCount(dofs[idof]);
	phgInfo(0, "dof[%d]: %d\n", idof, ndat[idof]);
    }
    MPI_Allgather(&ndat[0], ndof, MPI_OFFSET,
		  &cnts[0][0], ndof, MPI_OFFSET, MPI_COMM_WORLD);


    for (int idof = 0; idof < ndof; idof++) {
	phgInfo(0, "all dof[%d]: ", idof);
	offs[0][idof] = 0;
	for (int p = 0; p < nprocs; p++) {
	    offs[p+1][idof] = offs[p][idof] + cnts[p][idof];
	    cnts_int[idof][p] = cnts[p][idof];

	    phgInfo(0, "%d, ", cnts[p][idof]);
	}
	phgInfo(0, "\n");
    }


    //
    //
    //
    // Header
    //
    //
    //
    snprintf(file_name, NAME_MAX_LENGTH, "%s_%05d.nc", save_prefix, tstep);

    ret = ncmpi_create(MPI_COMM_WORLD, file_name,
		       NC_CLOBBER|NC_64BIT_OFFSET, MPI_INFO_NULL, &ncfile);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    // Dims
    int dimid_np, dimid_ndof, dimid_dof[ndof];
    ret = ncmpi_def_dim(ncfile, "nprocs", nprocs, &dimid_np);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    ret = ncmpi_def_dim(ncfile, "ndof", ndof, &dimid_ndof);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);


    // Time
    ret = ncmpi_put_att_double(ncfile, NC_GLOBAL, "Time",
			       NC_DOUBLE, 1, &time);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    
    // Dofs dims
    for (int idof = 0; idof < ndof; idof++) {

	// dof dim
	sprintf(tmp_str, "length_%s", dof_names[idof]);
	phgInfo(0, "len[%d] %d\n", idof, offs[nprocs][idof]);
	ret = ncmpi_def_dim(ncfile, tmp_str, offs[nprocs][idof], &dimid_dof[idof]);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);

	ret = ncmpi_def_var(ncfile, dof_names[idof], NC_DOUBLE, 1, &dimid_dof[idof], &var_id[idof]);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);

	ret = ncmpi_put_att_text(ncfile, var_id[idof], "long_name",
				 strlen(dof_names[idof]), dof_names[idof]);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);

    	ret = ncmpi_put_att_text(ncfile, var_id[idof], "FE_space",
				 strlen(dofs[idof]->type->name),
				 dofs[idof]->type->name);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);

	const int dim = dofs[idof]->dim;
    	ret = ncmpi_put_att_int(ncfile, var_id[idof], "dim",
				NC_INT, 1, &dim);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);

	FLOAT norm_L2 = phgDofNormL2(dofs[idof]);
    	ret = ncmpi_put_att_double(ncfile, var_id[idof], "L2_norm",
				NC_DOUBLE, 1, &norm_L2);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);


    	// ret = ncmpi_put_att_int(ncfile, var_id[idof], "localsize",
	// 			NC_INT, nprocs, cnts_int[idof]);
	// if (ret != NC_NOERR) handle_error(ret, __LINE__);

	//
	// localsize as var
	//
	sprintf(tmp_str, "%s%s", dof_names[idof], "_part");
	ret = ncmpi_def_var(ncfile, tmp_str,
			    NC_INT, 1,
			    &dimid_np, &var_part_id[idof]);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);

	sprintf(tmp_str, "localsize of dof %s", dof_names[idof]);
	ret = ncmpi_put_att_text(ncfile, var_part_id[idof], "long_name",
				 strlen(tmp_str),
				 tmp_str);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);
    }

    ret = ncmpi_enddef(ncfile);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);


    //
    //
    // Data
    //
    //
    for (int idof = 0; idof < ndof; idof++) {

	//
	// Part
	// 
	/* entering independent data mode */
	ret = ncmpi_begin_indep_data(ncfile);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);

	/* independently write values into netCDF variable */
	if (rank == 0) {
	    const MPI_Offset start = 0;
	    const MPI_Offset count = nprocs;
	    ret = ncmpi_put_vara_int(ncfile, var_part_id[idof],
				     &start, &count, cnts_int[idof]);
	    if (ret != NC_NOERR) handle_error(ret, __LINE__);
	}

	/* exiting independent data mode */
	ret = ncmpi_end_indep_data(ncfile);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);


	//
	// Data
	// 
	ret = ncmpi_put_vara_double_all(ncfile, var_id[idof],
					&offs[rank][idof], &cnts[rank][idof], dofs[idof]->data);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);
    }


    ret = ncmpi_close(ncfile);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    return;
}






void
load_state(int tstep, double *time, int ndof, DOF **dofs, char **dof_names)
//
//
// Note: The dofs are store and loaded with order, and the poisition is the dof id.
//       On the other hand, the name could be changes.
//
//
{
    int nprocs = phgNProcs, rank = phgRank;
    int ret, ncfile, dimid, ndims=1;
    MPI_Offset start, count=1;
    char file_name[NAME_MAX_LENGTH], tmp_str[NAME_MAX_LENGTH];

    MPI_Offset ndat[ndof];
    MPI_Offset cnts[nprocs][ndof];
    MPI_Offset offs[nprocs+1][ndof];
    int cnts_int[ndof][nprocs];	 // copy to int type, and transpose

    phgInfo(0, "\n* Load state\n");
    for (int idof = 0; idof < ndof; idof++) {
	ndat[idof] = DofGetDataCount(dofs[idof]);
	phgInfo(0, "dof[%d]: %d\n", idof, ndat[idof]);
    }
    MPI_Allgather(&ndat[0], ndof, MPI_OFFSET,
		  &cnts[0][0], ndof, MPI_OFFSET, MPI_COMM_WORLD);


    for (int idof = 0; idof < ndof; idof++) {
	phgInfo(0, "all dof[%d]: ", idof);
	offs[0][idof] = 0;
	for (int p = 0; p < nprocs; p++) {
	    offs[p+1][idof] = offs[p][idof] + cnts[p][idof];
	    cnts_int[idof][p] = cnts[p][idof];

	    phgInfo(0, "%d, ", cnts[p][idof]);
	}
	phgInfo(0, "\n");
    }


    //
    //
    //
    // Header
    //
    //
    //
    snprintf(file_name, NAME_MAX_LENGTH, "%s_%05d.nc", save_prefix, tstep);

    ret = ncmpi_open(MPI_COMM_WORLD, file_name,
		     NC_NOWRITE, MPI_INFO_NULL, &ncfile);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    ret = ncmpi_get_att_double(ncfile, NC_GLOBAL, "Time", time);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);
    phgInfo(0, "time: %f\n", *time);
    
    // Dims
#define GET_SIZE(name, dim)  {				\
	int dimid;					\
	ret = ncmpi_inq_dimid(ncfile, name, &dimid);	\
	if (ret != NC_NOERR)				\
	    handle_error(ret, __LINE__);		\
	ret = ncmpi_inq_dimlen(ncfile, dimid, &size);	\
	if (ret != NC_NOERR)				\
	    handle_error(ret, __LINE__);		\
    }

    MPI_Offset size;
    GET_SIZE("nprocs", size);
    if (size != nprocs) {
	phgError(1, "Saved with %d procs but read with %d procs !!!\n", size, nprocs);
    }

    GET_SIZE("ndof", size);
    if (size != ndof) {
	phgError(1, "Saved %d dofs but read %d dofs !!!\n", size, ndof);
    }

    // Dofs dims
    for (int idof = 0; idof < ndof; idof++) {
	int var_id, var_part_id;
    	ret = ncmpi_inq_varid(ncfile, dof_names[idof], &var_id);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);

	sprintf(tmp_str, "%s%s", dof_names[idof], "_part");
    	ret = ncmpi_inq_varid(ncfile, tmp_str,
			      &var_part_id);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);
	
	// dof dim
	sprintf(tmp_str, "length_%s", dof_names[idof]);
	GET_SIZE(tmp_str, size);
	if (size != offs[nprocs][idof]) {
	    phgError(1, "Saved dof %s size %d but read size %d !!!\n",
		     dof_names[idof], size, offs[nprocs][idof]);
	}


	// check fe space
	char fe_space_name[NAME_MAX_LENGTH];
	MPI_Offset str_len;

	ret = ncmpi_inq_attlen(ncfile, var_id, "FE_space", &str_len);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);
    	ret = ncmpi_get_att_text(ncfile, var_id, "FE_space", fe_space_name);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);
	fe_space_name[str_len] = '\0';	// Note nc get text has no ending !!!

	if (strcmp(fe_space_name, dofs[idof]->type->name)) {
	    phgPrintf("Saved dof %s space =%s= but read space =%s=  %d!!!\n",
		     dof_names[idof], dofs[idof]->type->name,
		      fe_space_name, strcmp(fe_space_name, dofs[idof]->type->name));
	    phgError(1, "Saved dof %s space %s but read space %s !!!\n",
		     dof_names[idof], dofs[idof]->type->name,
		     fe_space_name);
	}


	// check dof dim
	int dof_dim;
    	ret = ncmpi_get_att_int(ncfile, var_id, "dim", &dof_dim);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);
	if (dof_dim != dofs[idof]->dim) {
	    phgError(1, "Saved dof %s dim %d but read dim %d !!!\n",
		     dof_names[idof], dofs[idof]->dim, dof_dim);
	}


	// check local size
	int cnts_read[nprocs];
    	// ret = ncmpi_get_att_int(ncfile, var_id, "localsize", cnts_read);
	// if (ret != NC_NOERR) handle_error(ret, __LINE__);
	{
	    const MPI_Offset start = 0;
	    const MPI_Offset count = nprocs;
	    ret = ncmpi_get_vara_int_all(ncfile, var_part_id,
					 &start, &count, cnts_read); // read by root then bcast ???
	}
	for (int p = 0; p < nprocs; p++) {
	    if (cnts_int[idof][p] != cnts_read[p])
		phgError(1, "Saved dof %s cnts on proc %d is %d but read %d !!!\n",
			 dof_names[idof], p, cnts_int[idof][p], cnts_read[p]);
	}


	//
	// Load data
	//
    	ret = ncmpi_get_vara_double_all(ncfile, var_id,
    					&offs[rank][idof], &cnts[rank][idof], dofs[idof]->data);
    	if (ret != NC_NOERR) handle_error(ret, __LINE__);


	// check dof norm
	FLOAT norm_L2 = phgDofNormL2(dofs[idof]), norm_read;
    	ret = ncmpi_get_att_double(ncfile, var_id, "L2_norm", &norm_read);
	if (ret != NC_NOERR) handle_error(ret, __LINE__);
	if (norm_L2 > 0 && fabs(norm_L2 - norm_read) > 1e-8 * norm_L2) {
	    phgError(1, "Saved dof %s norm %30.15f but read %30.15f !!!\n",
		     dof_names[idof], norm_L2, norm_read);
	}
	phgInfo(0, "Load %s, norm %30.15f %30.15f, err %.10e\n", dof_names[idof],
		norm_L2, norm_read, fabs(norm_L2 - norm_read)
		);
    }




    ret = ncmpi_close(ncfile);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);

    MPI_Barrier(MPI_COMM_WORLD);
    phgPrintf("Pnetcdf read done.\n");

    return;
}
