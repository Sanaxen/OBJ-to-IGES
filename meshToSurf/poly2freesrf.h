#ifndef _POLY2FREESRF_H
#define _POLY2FREESRF_H

#include <vector>
#include "mesh_points.h"



class PolygonSurface {
	double* ctrlpnt_;
public:
	int	vtxnum_;
	int trinum_;
	double* vertex_;
	double* uv_;
	int* tri_;
	double* normal_;
	int chkgrid_[3];
	double start_[3];
	double step_[3];
	int		alcflg_;

	PolygonSurface(){
		ctrlpnt_ = NULL;
		alcflg_ = 0;
	};
	virtual ~PolygonSurface(){
	};

	void SetCPArea(double *cp ){
		ctrlpnt_ = cp;
	};
	double* GetCPArea(){
		return ctrlpnt_;
	};


	meshGrdid mgrid_;
	
	void Create(int chkgrid[3], int upatch, int vpatch, double rect[2][2],
			int vtxnum, int trinum, double* vertex, double* uv, int* tri, double* normal, 
			double* start, double* step);

	int Evalue(double* uvparam, double point[3] );
	void DumpUVLines(FILE* fp, int un, int vn);


};


#endif

	
