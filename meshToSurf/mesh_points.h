#ifndef _MESH_POINTS_H
#define _MESH_POINTS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <vector>


class meshGrdidCell {

public:
	float minA[3];
	float maxA[3];
	std::vector<int> facet;

	meshGrdidCell(){
		minA[0] = (float)99999999999999.0;
		minA[1] = (float)99999999999999.0;
		minA[2] = (float)99999999999999.0;

		maxA[0] = -(float)99999999999999.0;
		maxA[1] = -(float)99999999999999.0;
		maxA[2] = -(float)99999999999999.0;
		facet.clear();
	};

	virtual ~meshGrdidCell(){
		facet.clear();
	};

	void Add(float vmin[3], float vmax[3], int fid){
		if ( minA[0] > vmin[0] ) minA[0] = vmin[0];
		if ( minA[1] > vmin[1] ) minA[1] = vmin[1];
		if ( minA[2] > vmin[2] ) minA[2] = vmin[2];
		if ( maxA[0] < vmax[0] ) maxA[0] = vmax[0];
		if ( maxA[1] < vmax[1] ) maxA[1] = vmax[1];
		if ( maxA[2] < vmax[2] ) maxA[2] = vmax[2];
		facet.push_back( fid );
	};

	void DumpTri(FILE* fp, double* vtx, int* triid );

	int OnTri( double* point, double* vtx, int* triid, double& dist, double* on_point  );
};

class meshGrdid {
	float minA_[3];
	float maxA_[3];
	float len_[3];
	int gridsize_[3];
	std::vector< std::vector< std::vector<meshGrdidCell> > > cell_;
	float org_[3];
	int vtxnum_;
	int trinum_;
	double* vtx_;
	int* triid_;

	void Add(float vp[3][3], int fid);
public:
	double overlap_coef_;

	meshGrdid(){;};

	meshGrdid(int vtxnum, int trinum, double* vtx, int* triid){
		Init(vtxnum, trinum, vtx, triid);
	};

	void Init(int vtxnum, int trinum, double* vtx, int* triid){
		vtxnum_ = vtxnum;
		trinum_ = trinum;
		vtx_ = vtx;
		triid_ = triid;
		overlap_coef_ = 0.05;
	};

	void Create(int gridsize[3], float minA[3], float maxA[3] ){
		gridsize_[0] = gridsize[0];
		gridsize_[1] = gridsize[1];
		gridsize_[2] = gridsize[2];

		int i, j;
		cell_.resize(gridsize_[0]);
		for ( i = 0; i < gridsize_[0]; i++ ){
			cell_[i].resize(gridsize_[1] );

			for ( j = 0; j < gridsize_[1]; j++ ){
				cell_[i][j].resize(gridsize_[2] );
			}
		}
		len_[0] = (maxA[0] - minA[0])/gridsize_[0];
		len_[1] = (maxA[1] - minA[1])/gridsize_[1];
		len_[2] = (maxA[2] - minA[2])/gridsize_[2];
		minA_[0] = (float)minA[0]-len_[0]*0.01;
		minA_[1] = (float)minA[1]-len_[1]*0.01;
		minA_[2] = (float)minA[2]-len_[2]*0.01;
		maxA_[0] = (float)maxA[0]+len_[0]*0.01;
		maxA_[1] = (float)maxA[1]+len_[1]*0.01;
		maxA_[2] = (float)maxA[2]+len_[2]*0.01;
		
		len_[0] = (maxA_[0] - minA_[0])/gridsize_[0];
		len_[1] = (maxA_[1] - minA_[1])/gridsize_[1];
		len_[2] = (maxA_[2] - minA_[2])/gridsize_[2];



		org_[0] = minA_[0];
		org_[1] = minA_[1];
		org_[2] = minA_[2];
	};

	virtual ~meshGrdid(){
	};


	void GetTri(int id, double vtx[3][3] ){
		int vid[3];

		vid[0] = triid_[3*id  ];
		vid[1] = triid_[3*id+1];
		vid[2] = triid_[3*id+2];

		for (int i = 0; i < 3; i++ ){
			vtx[i][0] = vtx_[3*vid[i]];
			vtx[i][1] = vtx_[3*vid[i]+1];
			vtx[i][2] = vtx_[3*vid[i]+2];
		}
	};


	void mesh_to_grid(int gridsize[3]);
	void Dump( char* filename);

	int OnTri( double* point, double& dist, double* on_point  );
	void DumpOnTri(char* filename, double* point);

};


class PointToNormal{
	meshGrdid mgrid_;

public:
	PointToNormal(int vtxnum, int trinum, double* vtx, int* triid){
		mgrid_.Init(vtxnum, trinum, vtx, triid);
	};

	PointToNormal(){;};
	virtual ~PointToNormal(){;};


	void DumpGrig( char* filename ){
		if ( filename != NULL ){
			mgrid_.Dump( filename );
		}else{
			mgrid_.Dump( "PointToNormal_DumpGrig.txt" );
		}
	};
	void MakeNormal(int numpoint, double* points, double* normal, int grid[3], int* errnum, int** errid );
};


double* PointToNormal_PointsFromFile( char* filename, int* numpoint );
double* PointToNormal_PointsFromFile2( char* filename, int* numpoint );

#endif

