#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <vector>

#include "mesh_points.h"
#include "makemshpnts.h"
#include "TriToPoint.hpp"

double* PointToNormal_PointsFromFile( char* filename, int* numpoint )
{
   FILE* fp = fopen( filename, "r");

	char buf[256];
    *numpoint = 0;
    while( fgets(buf, 256, fp) != NULL ){
            if ( buf[0] != 'p' ) continue;
            (*numpoint)++;
    }
    fclose( fp );
	if ( *numpoint == 0 ) return NULL;

    double x, y, z;
    double* pnts = NULL;

    pnts = (double*)malloc( 3*(*numpoint)*sizeof(double));
    if ( pnts == NULL ){
            return NULL;
    }
    fp = fopen( filename, "r");


    int i = 0;
    while( fgets(buf, 256, fp) != NULL ){
            if ( buf[0] != 'p' ) continue;
            sscanf(buf, "p %lf %lf %lf", &x, &y, &z);
            pnts[3*i+0] = x;
            pnts[3*i+1] = y;
            pnts[3*i+2] = z;
            i++;
    }
    fclose( fp );

    return pnts;
}

double* PointToNormal_PointsFromFile2( char* filename, int* numpoint )
{
    FILE* fp = fopen( filename, "r");
    double* pnts = NULL;
	int i;

	char buf[256];
    *numpoint = 0;
	if ( fgets(buf, 256, fp) != NULL ){
		sscanf( buf, "%d", numpoint );
	}else{
		goto err;
	}

	double x, y, z;
    pnts = (double*)malloc( 3*(*numpoint)*sizeof(double));
    if ( pnts == NULL ){
            return NULL;
    }
	for (i = 0; i < *numpoint; i++ ){
        if( fgets(buf, 256, fp) != NULL ){
			sscanf( buf, "%lf %lf %lf", &x, &y, &z );
                pnts[3*i+0] = x;
                pnts[3*i+1] = y;
                pnts[3*i+2] = z;
		}else{
			goto err;
		}
	}
    fclose( fp );
	fp = NULL;
	printf( "%d points\n",*numpoint);

	goto rtn;
err:	;
	if ( fp != NULL ) fclose( fp );
	if ( pnts != NULL ) free( pnts );

rtn:	;
        return pnts;
}

void meshGrdidCell::DumpTri(FILE* fp, double* vtx, int* triid )
{
	std::vector<int>::iterator it = facet.begin();

	int fid;
	int	vid[3];
	double p[3][3];
	for ( ; it != facet.end(); it++ ){
		fid = *it;
		vid[0] = triid[3*fid  ];
		vid[1] = triid[3*fid+1];
		vid[2] = triid[3*fid+2];

		for ( int i = 0; i < 3; i++ ){
			p[i][0] = vtx[3*vid[i]  ];
			p[i][1] = vtx[3*vid[i]+1];
			p[i][2] = vtx[3*vid[i]+2];
		}
		fprintf(fp, "Plane%d=Triangle([%f,%f,%f],[%f,%f,%f],[%f,%f,%f]);\n", fid,
			p[0][0], p[0][1], p[0][2],
			p[1][0], p[1][1], p[1][2],
			p[2][0], p[2][1], p[2][2]);
	}
}
int meshGrdidCell::OnTri( double* point, double* vtx, int* triid, double& dist, double* on_point  )
{
	std::vector<int>::iterator it = facet.begin();

	int fid;
	int	vid[3];
	double p[3][3];

	double dist_min = dist;
	double wknpnt[3];
	int id = -1;
	
	for ( ; it != facet.end(); it++ ){
		fid = *it;
		vid[0] = triid[3*fid  ];
		vid[1] = triid[3*fid+1];
		vid[2] = triid[3*fid+2];

		for ( int i = 0; i < 3; i++ ){
			p[i][0] = vtx[3*vid[i]  ];
			p[i][1] = vtx[3*vid[i]+1];
			p[i][2] = vtx[3*vid[i]+2];
		}

		//int stat = TriToPoint(3, point, p[0], p[1], p[2], 0, wknpnt, &dist); 
		int stat = TriToPoint(point, p[0], p[1], p[2], wknpnt, dist);
		if ( stat != 0 ){
			continue;
		}
		if( dist < dist_min ){
			dist_min = dist;
			on_point[0] = wknpnt[0];
			on_point[1] = wknpnt[1];
			on_point[2] = wknpnt[2];
			id  = fid;
		}
	}
	return id;
}

void meshGrdid::Add(float vp[3][3], int fid)
{
	float vmin[3], vmax[3];
	vmin[0] = (float)99999999999999.0;
	vmin[1] = (float)99999999999999.0;
	vmin[2] = (float)99999999999999.0;

	vmax[0] = -(float)99999999999999.0;
	vmax[1] = -(float)99999999999999.0;
	vmax[2] = -(float)99999999999999.0;

	int i, j, k;
	for ( i = 0; i < 3; i++ ){
		if (vmin[0] > vp[i][0] ) vmin[0] = vp[i][0];		
		if (vmin[1] > vp[i][1] ) vmin[1] = vp[i][1];		
		if (vmin[2] > vp[i][2] ) vmin[2] = vp[i][2];

		if (vmax[0] < vp[i][0] ) vmax[0] = vp[i][0];		
		if (vmax[1] < vp[i][1] ) vmax[1] = vp[i][1];		
		if (vmax[2] < vp[i][2] ) vmax[2] = vp[i][2];
	}

	double x[2], y[2], z[2];
	for ( i = 0; i < gridsize_[0]; i++ ){
		x[0] = org_[0] + len_[0]*((float)i    );
		x[1] = org_[0] + len_[0]*((float)(i+1));

		if ( x[0] > vmax[0] )continue; 
		if ( x[1] < vmin[0] )continue; 
		for ( j = 0; j < gridsize_[1]; j++ ){
				y[0] = org_[1] + len_[1]*((float)j    );
				y[1] = org_[1] + len_[1]*((float)(j+1));

				if ( y[0] > vmax[1] )continue; 
				if ( y[1] < vmin[1] )continue;
				
				for ( k = 0; k < gridsize_[2]; k++ ){
					z[0] = org_[2] + len_[2]*((float)k    );
					z[1] = org_[2] + len_[2]*((float)(k+1));

					if ( z[0] > vmax[2] )continue; 
					if ( z[1] < vmin[2] )continue;
					
					cell_[i][j][k].Add( vmin, vmax, fid );
				}
		}
	}
}

void meshGrdid::mesh_to_grid(int gridsize[3])
{
	int i;
	int vid[3];
	float	vp[3][3];
	float	minA[3], maxA[3];

	minA[0] = (float)99999999999999.0;
	minA[1] = (float)99999999999999.0;
	minA[2] = (float)99999999999999.0;

	maxA[0] = -(float)99999999999999.0;
	maxA[1] = -(float)99999999999999.0;
	maxA[2] = -(float)99999999999999.0;

	for ( i = 0; i < vtxnum_; i++ ){
		if ( minA[0] > vtx_[3*i  ] ) minA[0] = vtx_[3*i  ];
		if ( minA[1] > vtx_[3*i+1] ) minA[1] = vtx_[3*i+1];
		if ( minA[2] > vtx_[3*i+2] ) minA[2] = vtx_[3*i+2];
		if ( maxA[0] < vtx_[3*i  ] ) maxA[0] = vtx_[3*i  ];
		if ( maxA[1] < vtx_[3*i+1] ) maxA[1] = vtx_[3*i+1];
		if ( maxA[2] < vtx_[3*i+2] ) maxA[2] = vtx_[3*i+2];
	}

	Create(gridsize, minA, maxA);


	for ( i = 0; i < trinum_; i++ ){
		vid[0] = triid_[3*i  ];
		vid[1] = triid_[3*i+1];
		vid[2] = triid_[3*i+2];

		vp[0][0] = vtx_[3*vid[0]  ];
		vp[0][1] = vtx_[3*vid[0]+1];
		vp[0][2] = vtx_[3*vid[0]+2];

		vp[1][0] = vtx_[3*vid[1]  ];
		vp[1][1] = vtx_[3*vid[1]+1];
		vp[1][2] = vtx_[3*vid[1]+2];

		vp[2][0] = vtx_[3*vid[2]  ];
		vp[2][1] = vtx_[3*vid[2]+1];
		vp[2][2] = vtx_[3*vid[2]+2];

		Add(vp, i);
	}
}

void meshGrdid::Dump(char* filename)
{
	FILE* fp = fopen(filename, "w");
	if ( fp == NULL ){
		return;
	}

	double x[2], y[2], z[2];
	int	color[3];
	int i, j, k;

	for ( i = 0; i < gridsize_[0]; i++ ){
		x[0] = org_[0] + len_[0]*((float)i    );
		x[1] = org_[0] + len_[0]*((float)(i+1));

			for ( j = 0; j < gridsize_[1]; j++ ){
				y[0] = org_[1] + len_[1]*((float)j    );
				y[1] = org_[1] + len_[1]*((float)(j+1));

					for ( k = 0; k < gridsize_[2]; k++ ){
						z[0] = org_[2] + len_[2]*((float)k    );
						z[1] = org_[2] + len_[2]*((float)(k+1));

						color[0] = rand()%255;
						color[1] = rand()%255;
						color[2] = rand()%255;

							fprintf(fp, "modalcolor(%d,%d,%d);\n", color[0], color[1], color[2] );

							fprintf(fp, "Line=[%f,%f,%f],[%f,%f,%f];\n",
										x[0],y[0],z[0],x[1],y[0],z[0]);
							fprintf(fp, "Line=[%f,%f,%f],[%f,%f,%f];\n",
										x[0],y[0],z[0],x[0],y[1],z[0]);
							fprintf(fp, "Line=[%f,%f,%f],[%f,%f,%f];\n",
										x[0],y[0],z[0],x[0],y[0],z[1]);

							fprintf(fp, "Line=[%f,%f,%f],[%f,%f,%f];\n",
										x[0],y[0],z[1],x[1],y[0],z[1]);
							fprintf(fp, "Line=[%f,%f,%f],[%f,%f,%f];\n",
										x[0],y[0],z[1],x[0],y[1],z[1]);

							fprintf(fp, "Line=[%f,%f,%f],[%f,%f,%f];\n",
										x[1],y[1],z[1],x[0],y[1],z[1]);
							fprintf(fp, "Line=[%f,%f,%f],[%f,%f,%f];\n",
										x[1],y[1],z[1],x[1],y[0],z[1]);

							fprintf(fp, "Line=[%f,%f,%f],[%f,%f,%f];\n",
										x[1],y[1],z[1],x[1],y[1],z[0]);

							fprintf(fp, "Line=[%f,%f,%f],[%f,%f,%f];\n",
										x[1],y[0],z[1],x[1],y[0],z[0]);

							fprintf(fp, "Line=[%f,%f,%f],[%f,%f,%f];\n",
										x[0],y[1],z[1],x[0],y[1],z[0]);

							fprintf(fp, "Line=[%f,%f,%f],[%f,%f,%f];\n",
										x[1],y[0],z[0],x[1],y[1],z[0]);

							fprintf(fp, "Line=[%f,%f,%f],[%f,%f,%f];\n",
										x[1],y[1],z[0],x[0],y[1],z[0]);
							cell_[i][j][k].DumpTri(fp,  vtx_, triid_);
					}
			}
	}
	fprintf(fp, "end;\n");
	fclose(fp);
}

int meshGrdid::OnTri( double* point, double& dist, double* on_point  )
{
	double x[2], y[2], z[2];
	int i, j, k;

	double dist_min = 9999999999999999.0;
	int		id = -1;
	int		wkid;
	double wkdist;
	double wkpnt[3];

	wkdist = dist_min;
	for ( i = 0; i < gridsize_[0]; i++ ){
		x[0] = org_[0] + len_[0]*((float)i    )-overlap_coef_*len_[0];
		x[1] = org_[0] + len_[0]*((float)(i+1))+overlap_coef_*len_[0];

		if ( x[0] > point[0] )continue; 
		if ( x[1] < point[0] )continue; 
	for ( j = 0; j < gridsize_[1]; j++ ){
			y[0] = org_[1] + len_[1]*((float)j    )-overlap_coef_*len_[1];
			y[1] = org_[1] + len_[1]*((float)(j+1))+overlap_coef_*len_[1];

			if ( y[0] > point[1] )continue; 
			if ( y[1] < point[1] )continue;
				
			for ( k = 0; k < gridsize_[2]; k++ )
			{
				z[0] = org_[2] + len_[2]*((float)k    )-overlap_coef_*len_[2];
				z[1] = org_[2] + len_[2]*((float)(k+1))+overlap_coef_*len_[2];

				if ( z[0] > point[2] )continue; 
				if ( z[1] < point[2] )continue;
					
				wkdist = dist_min;
				wkid = cell_[i][j][k].OnTri( point, vtx_, triid_, wkdist, wkpnt);
				if ( wkdist < dist_min ){
					dist_min = wkdist;
					on_point[0] = wkpnt[0];
					on_point[1] = wkpnt[1];
					on_point[2] = wkpnt[2];
					id = wkid;
				}
			}
		}
	}

	dist = dist_min;
	return id;
}

void meshGrdid::DumpOnTri(char* filename, double* point)
{
	FILE* fp = fopen(filename, "w");
	if ( fp == NULL ){
		return;
	}

	double x[2], y[2], z[2];
	int	color[3];
	int i, j, k;

	for ( i = 0; i < gridsize_[0]; i++ ){
		x[0] = org_[0] + len_[0]*((float)i    );
		x[1] = org_[0] + len_[0]*((float)(i+1));

			for ( j = 0; j < gridsize_[1]; j++ ){
				y[0] = org_[1] + len_[1]*((float)j    );
				y[1] = org_[1] + len_[1]*((float)(j+1));

					for ( k = 0; k < gridsize_[2]; k++ ){
						z[0] = org_[2] + len_[2]*((float)k    );
						z[1] = org_[2] + len_[2]*((float)(k+1));

						color[0] = rand()%255;
						color[1] = rand()%255;
						color[2] = rand()%255;

							fprintf(fp, "modalcolor(%d,%d,%d);\n", color[0], color[1], color[2] );

							fprintf(fp, "Line=[%f,%f,%f],[%f,%f,%f];\n",
										x[0],y[0],z[0],x[1],y[0],z[0]);
							fprintf(fp, "Line=[%f,%f,%f],[%f,%f,%f];\n",
										x[0],y[0],z[0],x[0],y[1],z[0]);
							fprintf(fp, "Line=[%f,%f,%f],[%f,%f,%f];\n",
										x[0],y[0],z[0],x[0],y[0],z[1]);

							fprintf(fp, "Line=[%f,%f,%f],[%f,%f,%f];\n",
										x[0],y[0],z[1],x[1],y[0],z[1]);
							fprintf(fp, "Line=[%f,%f,%f],[%f,%f,%f];\n",
										x[0],y[0],z[1],x[0],y[1],z[1]);

							fprintf(fp, "Line=[%f,%f,%f],[%f,%f,%f];\n",
										x[1],y[1],z[1],x[0],y[1],z[1]);
							fprintf(fp, "Line=[%f,%f,%f],[%f,%f,%f];\n",
										x[1],y[1],z[1],x[1],y[0],z[1]);

							fprintf(fp, "Line=[%f,%f,%f],[%f,%f,%f];\n",
										x[1],y[1],z[1],x[1],y[1],z[0]);

							fprintf(fp, "Line=[%f,%f,%f],[%f,%f,%f];\n",
										x[1],y[0],z[1],x[1],y[0],z[0]);

							fprintf(fp, "Line=[%f,%f,%f],[%f,%f,%f];\n",
										x[0],y[1],z[1],x[0],y[1],z[0]);

							fprintf(fp, "Line=[%f,%f,%f],[%f,%f,%f];\n",
										x[1],y[0],z[0],x[1],y[1],z[0]);

							fprintf(fp, "Line=[%f,%f,%f],[%f,%f,%f];\n",
										x[1],y[1],z[0],x[0],y[1],z[0]);
							cell_[i][j][k].DumpTri(fp,  vtx_, triid_);
					}
			}
	}
	double dist;
	double on_point[3];

	OnTri(point, dist, on_point );

	fprintf(fp, "modalcolor(10, 200, 255);modalsize(5);\n");
	fprintf(fp, "P1=[%f,%f,%f];\n", point[0], point[1], point[2] );
	fprintf(fp, "P2=[%f,%f,%f];\n", on_point[0], on_point[1], on_point[2] );
	fprintf(fp, "LineP1P2=P1,P2;\n");

	fprintf(fp, "end;\n");
	fclose(fp);
}


void PointToNormal::MakeNormal(int numpoint, double* points, double* normal, int grid[3], int* errnum, int** errid )
{
	int id;
	double vtx[3][3];
	double	dist;
	double on_point[3];

	*errnum = 0;
	*errid = NULL;
	mgrid_.mesh_to_grid(grid);
	for ( int i = 0; i < numpoint; i++ ){
		this->mgrid_.overlap_coef_ = 0.001;
		do{
			id = this->mgrid_.OnTri(&(points[3*i]), dist, on_point);
			if ( id < 0 ){
				normal[3*i  ] = 0.0;
				normal[3*i+1] = 0.0;
				normal[3*i+2] = 0.0;
				this->mgrid_.overlap_coef_ *= 1.005;
				continue;
			}
			break;
		}while(this->mgrid_.overlap_coef_ < 2.0);

		if (id < 0 ){
			printf("normal vector error.:pointID:%d\n", i);
			if ( *errid == NULL ){
				*errid = (int*)malloc(numpoint*sizeof(int));
				if ( *errid == NULL ) throw("no enough memoery.");
				memset(*errid, '\0', numpoint*sizeof(int));
			}
			(*errid)[i] = 1;
			(*errnum)++;
			continue;
		}

#if 10
		this->mgrid_.GetTri(id, vtx);

		double	drc1[3];
		double	drc2[3];
		double	drc3[3];
		double	norm[3];
		double	lnglng;

		drc1[0] = vtx[1][0] - vtx[0][0];
		drc1[1] = vtx[1][1] - vtx[0][1];
		drc1[2] = vtx[1][2] - vtx[0][2];
		drc2[0] = vtx[2][0] - vtx[1][0];
		drc2[1] = vtx[2][1] - vtx[1][1];
		drc2[2] = vtx[2][2] - vtx[1][2];
		drc3[0] = vtx[0][0] - vtx[2][0];
		drc3[1] = vtx[0][1] - vtx[2][1];
		drc3[2] = vtx[0][2] - vtx[2][2];
		norm[0] = drc1[2] * drc3[1] - drc1[1] * drc3[2];
		norm[1] = drc1[0] * drc3[2] - drc1[2] * drc3[0];
		norm[2] = drc1[1] * drc3[0] - drc1[0] * drc3[1];
		lnglng = norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2];
		double lng = sqrt( lnglng );
		//lng = 1.0;

		normal[3*i  ] = norm[0]/lng;
		normal[3*i+1] = norm[1]/lng;
		normal[3*i+2] = norm[2]/lng;
#else
		double	norm[3];
		double	lnglng;

		norm[0] = on_point[0] - points[3*i  ];
		norm[1] = on_point[1] - points[3*i+1];
		norm[2] = on_point[2] - points[3*i+2];

		lnglng = norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2];
		double lng = sqrt( lnglng );

		normal[3*i  ] = norm[0]/lng;
		normal[3*i+1] = norm[1]/lng;
		normal[3*i+2] = norm[2]/lng;
#endif
	}
}