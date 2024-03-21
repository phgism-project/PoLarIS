/*
 *
 * Read Netcdf data of topo info and beta square,
 *   single test file: interp-data.c .
 *
 *  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <netcdf.h>
#include "ins.h"
#include "phg/netcdf-utils.h" 




typedef struct GEO_DATA_ {
    char *filename;
    int nx, ny, nt;
    float *x, *y;  /* unit: meter */
    float dx, dy;
    float x0, y0;
    char flip_x, flip_y;
    /* [ny, nx */
    float *bed, *thickness, *surface;  /* unit: meter */
    float *beta2;		       /* beta square */

    /* [ny, nx] or [nt, ny, nx] */
    float *temp;		       
    float *accu;

    BOOLEAN temp_time_dependent;
    BOOLEAN accu_time_dependent;
} GEO_DATA;


static GEO_DATA *geo_data = NULL;
static int ncid = -1;



/* ---------------------------------------------------------------------
 *
 *   Load 
 *
 * --------------------------------------------------------------------- */


void 
load_topo_data(const char *filename)
/* Load netcdf file */
{
    int retval;
    int dimid, varid;
    char name[1000];

    phgPrintf("Loading topo file: %s\n", filename);

    assert(geo_data == NULL);

    if ((retval = nc_open(filename, NC_NOWRITE, &ncid)))
	NC_ERR(retval);

    
    geo_data = calloc(1, sizeof(GEO_DATA));

    
    /*
     * Read in dims
     * */
    size_t nx, ny;
    size_t nn;
    size_t nt = 1; 
    
    NC_READ_DIM("x", nx);
    NC_READ_DIM("y", ny);
    nn = nx * ny;

    {
	char name[1024];					
	if ((retval = nc_inq_dimid(ncid, "time", &dimid))) {
	    phgInfo(0, "nc doesn't has dim time\n");
	    //nt = 1;
	}
	else {
	    if ((retval = nc_inq_dim(ncid, dimid, name,		
				     &nt))) {
		NC_ERR(retval);
	    } 
	    phgInfo(0, "nc has dim time, %d\n", nt);
	}
    }

    
    /* Read in coords */
    float *coordx = calloc(nx + ny, sizeof(float));
    float *coordy = coordx + nx;


    NC_CHECK_DIM_1D("x", nx);    
    NC_CHECK_DIM_1D("y", ny);    

    nc_type type_coord;
    if ((retval = nc_inq_varid(ncid, "x", &varid)))
	NC_ERR(retval);
    if ((retval = nc_inq_vartype(ncid, varid,
				 &type_coord)))
	NC_ERR(retval);

    if (type_coord == NC_INT) {
	/* int to float */
	int i;
	int *coordx_int = calloc(nx + ny, sizeof(int));
	int *coordy_int = coordx_int + nx;

	NC_READ_INT("x", coordx_int);
	NC_READ_INT("y", coordy_int);

	for (i = 0; i < nx; i++)
	    coordx[i] = coordx_int[i];

	for (i = 0; i < ny; i++)
	    coordy[i] = coordy_int[i];

	free(coordx_int);
    }
    else if (type_coord == NC_FLOAT) {
	NC_READ_FLOAT("x", coordx);
	NC_READ_FLOAT("y", coordy);
    }
    else {
	phgError(-1, "Unknown coord type %d!!!\n", type_coord);
    }

        
    float *bed = calloc(nn, sizeof(float));
    NC_READ_FLOAT("bed", bed);
    NC_CHECK_DIM_2D("bed", ny, nx);    


    float *surface = calloc(nn, sizeof(float));
    NC_READ_FLOAT("surface", surface);
    NC_CHECK_DIM_2D("surface", ny, nx);    


    float *thickness = calloc(nn, sizeof(float));
    NC_READ_FLOAT("thickness", thickness);
    NC_CHECK_DIM_2D("thickness", ny, nx);    

    float *beta2 = calloc(nn, sizeof(float));
    retval = nc_inq_varid(ncid, "beta2", &varid);
    if (retval == NC_ENOTVAR) {
	int i;
	phgPrintf("Beta2 not found in %s, use beta2 = 1000000\n",
		   filename);
	for (i = 0; i < nn; i++)
	    beta2[i] = 1000000;
    }
    else {
	phgPrintf("Reading beta\n");
	NC_READ_FLOAT("beta2", beta2);
	NC_CHECK_DIM_2D("beta2", ny, nx);    
    }


    float *temp = calloc(nn, sizeof(float));
    retval = nc_inq_varid(ncid, "surface_temperature", &varid);
    if (retval == NC_ENOTVAR) {
	int i;
	phgPrintf("Surface temperature not found in %s, use T_surf = TEMP_WATER - 10\n",
		   filename);
	for (i = 0; i < nn; i++)
	    temp[i] = TEMP_WATER - 10;
    }
    else {
	phgPrintf("Reading surface temperature\n");

	int dim;
	NC_GET_NDIM("surface_temperature", dim);

	if (dim == 2) {
	    NC_CHECK_DIM_2D("surface_temperature", ny, nx);
	}
	else {
	    assert(nt > 0);
	    free(temp);
	    temp = calloc(nt * nn, sizeof(float));
	    NC_CHECK_DIM_3D("surface_temperature", nt, ny, nx);
	    if (nt > 1) {
		geo_data->temp_time_dependent = TRUE;
		phgInfo(0, "time dependent surface temp\n");
	    }
	}

	NC_READ_FLOAT("surface_temperature", temp);
    }


    float *accu = calloc(nn, sizeof(float));
    retval = nc_inq_varid(ncid, "accumulation_rate", &varid);
    if (retval == NC_ENOTVAR) {
	phgPrintf("Accumulation rate not found in %s, use accu = 0\n",
		   filename);
    }
    else {
	phgPrintf("Reading accumulation rate\n");

	int dim;
	NC_GET_NDIM("accumulation_rate", dim);

	if (dim == 2) {
	    NC_CHECK_DIM_2D("accumulation_rate", ny, nx);
	}
	else {
	    assert(nt > 0);
	    free(accu);
	    accu = calloc(nt * nn, sizeof(float));
	    NC_CHECK_DIM_3D("accumulation_rate", nt, ny, nx);
	    if (nt > 1) {
		geo_data->accu_time_dependent = TRUE;
		phgInfo(0, "time dependent accumulation\n");
	    }
	}

	NC_READ_FLOAT("accumulation_rate", accu);
    }




    float dx = coordx[1] - coordx[0];
    float dy = coordy[1] - coordy[0];
    float x0 = coordx[0];
    float y0 = coordy[0];
    char flip_x = 0;
    char flip_y = 0;
    if (dx < 0) {
	dx *= -1;
	flip_x = 1;
	x0 = coordx[nx-1];
	phgInfo(3, "nc topo flip x\n");
    }
    if (dy < 0) {
	dy *= -1;
	flip_y = 1;
	y0 = coordy[ny-1];
	phgInfo(3, "nc topo flip y\n");
    }

    geo_data->nx = nx;
    geo_data->ny = ny;
    geo_data->nt = nt;
    geo_data->x = coordx;
    geo_data->y = coordy;
    
    geo_data->dx = dx;
    geo_data->dy = dy;
    geo_data->flip_x = flip_x;
    geo_data->flip_y = flip_y;
    geo_data->x0 = x0;
    geo_data->y0 = y0;
    
    geo_data->bed = bed;
    geo_data->surface = surface;
    geo_data->thickness = thickness;
    geo_data->beta2 = beta2;

    geo_data->temp = temp;
    geo_data->accu = accu;

    if ((retval = nc_close(ncid)))
	NC_ERR(retval);

    return;
}








/* ---------------------------------------------------------------------
 *
 *   Interpolate
 *
 * --------------------------------------------------------------------- */

double
interp_topo_dataT(char var_type,
		 double x, double y, double t)
/*
 * Interp data from a[ny][nx]
 * Input: x, y (unit km)
 *        t (year)
 *
 *
 *  */
{
    int i, j;
    double a00, a01, a10, a11;
    double wx, wy, wt;

    assert(geo_data != NULL);

    int nx = geo_data->nx;
    int ny = geo_data->ny;
    int nt = geo_data->nt;
    float dx = geo_data->dx;
    float dy = geo_data->dy;
    float *X = geo_data->x;
    float *Y = geo_data->y;
    float X0 = geo_data->x0;
    float Y0 = geo_data->y0;
    int nn = nx * ny;
    float scale = 1.;

    
    /* relative position */
    i = (x * 1000. - X0) / dx;
    j = (y * 1000. - Y0) / dy;


    /* The outside is interp as border */
    if (i < 0) i = 0;
    if (j < 0) j = 0;
    if (i >= nx-1) i = nx-2;
    if (j >= ny-1) j = ny-2;  /* range [0, n-2] */
    /* if (ii < 0) ii = 0; */
    /* if (jj < 0) jj = 0; */
    /* if (ii > nx-1) ii = nx-1; */
    /* if (jj > ny-1) jj = ny-1; */
    
    int i0 = (!geo_data->flip_x) ? i : nx-1 - i;
    int i1 = (!geo_data->flip_x) ? i + 1: nx-2 - i;
    /* [i0, i1] range:
     *  [0, 1] -> [nx-2, nx-1]
     *  [nx-2, nx-1] -> [1, 0] */

    int j0 = (!geo_data->flip_y) ? j : ny-1 - j;
    int j1 = (!geo_data->flip_y) ? j + 1: ny-2 - j;
   

#if NC_DEBUG
    printf("\n");
    printf("interp x %d %d\n", x * 1000 - X0, y * 1000 - Y0);
    printf("interp I [%d %d] [%d %d]\n", i0, i1, j0, j1);
    printf("interp B [%d %d] [%d %d]\n",
	   X[i0] - X0,
	   X[i1] - X0,
	   Y[j0] - Y0,
	   Y[j1] - Y0);
#endif	/* NC_DEBUG */


#if NC_DEBUG
    printf("interp a [%f %f, %f %f]\n", a00, a01, a10, a11);
#endif	/* NC_DEBUG */

    wx = (x * 1000 - X[i0]) / ((double) dx);
    wy = (y * 1000 - Y[j0]) / ((double) dy);


#if NC_DEBUG
    printf("interp w [%f %f]\n", wx, wy);
#endif	/* NC_DEBUG */

    if (wx < 0.) wx = 0.;
    if (wy < 0.) wy = 0.;
    if (wx > 1.) wx = 1.;
    if (wy > 1.) wy = 1.;

    int it = floor(t);
    if (it >= nt-1) it = nt-2;
    if (it < 0) it = 0;
	    
    wt = (t - it);
    if (wt < 0.) wt = 0.;
    if (wt > 1.) wt = 1.;

    
#if NC_DEBUG
    printf("interp w [%f %f]\n", wx, wy);
#endif	/* NC_DEBUG */
    
    float *var_data = NULL;
    float *var_data2 = NULL;

    if (var_type == 't') {
	/* top */
	var_data = geo_data->surface;
	scale = 0.001;		/* m to km */
    } 
    else if (var_type == 'b') {
	var_data = geo_data->bed;
	scale = 0.001;		/* m to km */
    }
    else if (var_type == 'h') {
	var_data = geo_data->thickness;
	scale = 0.001;		/* m to km */
    }
    else if (var_type == 'B') {
	/* beta square */
	var_data = geo_data->beta2;
    }
    else if (var_type == 'T') {
	/* surface temperature */
	if (!geo_data->temp_time_dependent) {
	    var_data = geo_data->temp ;
	}
	else {
	    var_data = geo_data->temp + nn * it;
	    var_data2 = geo_data->temp + nn * (it + 1);
	}
    }
    else if (var_type == 'q') {
	/* accumulation rate */
	if (!geo_data->accu_time_dependent) {
	    var_data = geo_data->accu;
	}
	else {
	    var_data = geo_data->accu + nn * it;
	    var_data2 = geo_data->accu + nn * (it + 1);
	}
    }
    else {
	phgError(-1, "Topo data interp unknown var type %c!!!\n", 
		 var_type);
    }

    /* Bilinear interpolation */
    a00 = var_data[i0 + j0*nx];
    a01 = var_data[i0 + j1*nx];
    a10 = var_data[i1 + j0*nx];
    a11 = var_data[i1 + j1*nx];


    double a = a00 * (1-wx) * (1-wy)
	+ a01 * (1-wx) * wy
	+ a10 * wx * (1-wy)
	+ a11 * wx * wy;

    if (var_data2 != NULL) {
	a00 = var_data2[i0 + j0*nx];
	a01 = var_data2[i0 + j1*nx];
	a10 = var_data2[i1 + j0*nx];
	a11 = var_data2[i1 + j1*nx];
	
	double a2 = a00 * (1-wx) * (1-wy)
		+ a01 * (1-wx) * wy
		+ a10 * wx * (1-wy)
		+ a11 * wx * wy;

	a = a * (1-wt) + a2 * wt;

	/* if (var_type == 'q') */
	/* 	abort(); */
    }

    /* if (var_type == 'T') */
    /* 	abort(); */
    
    a *= scale;


#if NC_DEBUG
    printf("interp %f\n", a);
#endif	/* NC_DEBUG */

    return a;
}







