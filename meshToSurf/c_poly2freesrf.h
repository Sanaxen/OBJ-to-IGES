#ifndef _C_POLY2FREESRF_H
#define _C_POLY2FREESRF_H

void ReadMeshUV( char* filename, int* vtxnum, int* trinum, double** vertex, double** uv, int** tri );
void ReadMeshOBJUV(char* filename, int* vtxnum, int* trinum, double** vertex, double** uv, int** tri);

int Polygon2SbezFitC(int* chkgridA, double *ctrlpnt, int upatch, int vpatch, double rectA[][2], 
					 int vtxnum, int trinum, double* vertex, double* uv, int* tri, double* normal, 
					double* start, double* step, double* (*fit)(int, double*) );

void SrBez(double cp[4][4][3], double* uv, double p[3]);

int	SBezFitCP0(double pnt[4][4][3], double patch[4][4][3]);
int	Conv16_3PtrTo4_4_3Ary(double *ptr, double arry[4][4][3]);
int	Conv4_4_3AryTo16_3Ptr(double arry[4][4][3], double *ptr);

#endif
