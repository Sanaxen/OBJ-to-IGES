#ifndef _MAKEMSHPNT_H
#define _MAKEMSHPNT_H

#ifdef __cplusplus
extern "C" {
#endif

int MakeMeshFromPoints(int* numpnt, double* points, double* normal, 
					   int gridnum[3], int sampling, int grid, 
					   int *vtxnumP, int *facenumP, double** vtxP, int** faceP,
					    void* (*alc)(int), void* (*realc)(void*, int),
						char* rootname, int debug_out_mesh);

int DeleteMakeMeshFromPoints( double* vtxP, int* faceP);
#ifdef __cplusplus
};
#endif

#endif