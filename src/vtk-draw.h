#ifndef VTK_DRAW
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
//#include <bitset>
#include <stdarg.h>
#include "phg.h"

/*
 * Vtk debug control.
 * */
//#define VTK_DEBUG 

#define POINT_SIZE 0.005
#define TEXT_SIZE 0.005

typedef struct VTK_HANDLER_ {
    void (*hi)();
    void (*Init)(GRID *g);
    void (*ReInit)();
    void (*Pause)();
    void (*Finalize)();

    void (*DrawPoint)(double *p, char *words);
    void (*DrawLine)(double *x1, double *x2);
    void (*DrawEdge)(SIMPLEX *e, int edge);
    void (*DrawTria)(double *v0, double *v1, double *v2);
    void (*DrawFace)(SIMPLEX *e, int face);
    void (*DrawElement)(SIMPLEX *e); 
    void (*DrawRemoteFace)();
    void (*DrawBbox)();
    void (*DrawMesh)();

    void (*SetViewstop)(int stop);
    void (*SetColor)(char *colorname);
    void (*SetTransparent)(float a);

    void (*TmpActorsBegin)();
    void (*TmpActorsClear)();
} VTK_HANDLER;

extern VTK_HANDLER vtk;
extern char vtk_tmp_str[1000];
extern int vtk_verb;


#ifdef VTK_DEBUG 
# define vtkDrawEdge(verb, ...)        {if ((verb) <= phgVerbosity) { vtk.DrawEdge(__VA_ARGS__);}}
# define vtkDrawElement(verb, ...)     {if ((verb) <= phgVerbosity) { vtk.DrawElement(__VA_ARGS__);}}
# define vtkDrawTria(verb, ...)        {if ((verb) <= phgVerbosity) { vtk.DrawTria(__VA_ARGS__);}}
# define vtkDrawFace(verb, ...)        {if ((verb) <= phgVerbosity) { vtk.DrawFace(__VA_ARGS__);}}
# define vtkDrawLine(verb, ...)        {if ((verb) <= phgVerbosity) { vtk.DrawLine(__VA_ARGS__);}}
# define vtkDrawRemoteFace(verb, ...)  {if ((verb) <= phgVerbosity) { vtk.DrawRemoteFace(__VA_ARGS__);}}
# define vtkDrawBbox(verb, ...)        {if ((verb) <= phgVerbosity) { vtk.DrawBbox(__VA_ARGS__);}}
# define vtkDrawMesh(verb, ...)        {if ((verb) <= phgVerbosity) { vtk.DrawMesh(__VA_ARGS__);}}
# define vtkDrawPoint(verb, ...)       {if ((verb) <= phgVerbosity) { vtk.DrawPoint(__VA_ARGS__);}}
# define vtkhi(verb, ...)              {if ((verb) <= phgVerbosity) { vtk.hi(__VA_ARGS__);}}
# define vtkInit(verb, ...)            {if ((verb) <= phgVerbosity) { vtk.Init(__VA_ARGS__);}}
# define vtkPause(verb, ...)           {				\
	if ((verb) <= phgVerbosity) {					\
	    fprintf(stderr, "*%2d* vtk pause \tfile:%-25s line:%5d\tfunc:%-25s\n", \
		    phgRank, __FILE__, __LINE__, __FUNCTION__);		\
	    vtk.Pause(__VA_ARGS__);					\
	}								\
    }
# define vtkStop(verb, ...)            {if ((verb) <= phgVerbosity) {while(1) {vtkPause(verb)}}};
# define vtkSetColor(verb, ...)        {if ((verb) <= phgVerbosity) { vtk.SetColor(__VA_ARGS__);}}
# define vtkSetTransparent(verb, ...)  {if ((verb) <= phgVerbosity) { vtk.SetTransparent( __VA_ARGS__ );}}
# define vtkTmpActorsBegin(verb, ...)  {if ((verb) <= phgVerbosity) { vtk.TmpActorsBegin( __VA_ARGS__ );}}
# define vtkTmpActorsClear(verb, ...)  {if ((verb) <= phgVerbosity) { vtk.TmpActorsClear( __VA_ARGS__ );}}
# define vtkSetString(verb, ...)       {if ((verb) <= phgVerbosity) { sprintf(vtk_tmp_str, __VA_ARGS__);}}
#else
# define vtkDrawEdge(verb, ...)        
# define vtkDrawElement(verb, ...)     
# define vtkDrawTria(verb, ...)        
# define vtkDrawFace(verb, ...)        
# define vtkDrawLine(verb, ...)        
# define vtkDrawRemoteFace(verb, ...)  
# define vtkDrawBbox(verb, ...)  
# define vtkDrawMesh(verb, ...)        
# define vtkDrawPoint(verb, ...)       
# define vtkhi(verb, ...)              
# define vtkInit(verb, ...) 
# define vtkPause(verb, ...)           
# define vtkStop(verb, ...)
# define vtkSetColor(verb, ...)        
# define vtkSetTransparent(verb, ...)  
# define vtkTmpActorsBegin(verb, ...)  
# define vtkTmpActorsClear(verb, ...)  
# define vtkSetString(verb, ...)
#endif	/* VTK_DEBUG */

/* Utils */
#define SET_MULTI_COLOR(i)			\
    switch(i % 3) {				\
    case 0: vtkSetColor(0, "red");    break;	\
    case 1: vtkSetColor(0, "yellow"); break;	\
    case 2: vtkSetColor(0, "green");  break;	\
    }


#define VTK_DRAW
#endif
