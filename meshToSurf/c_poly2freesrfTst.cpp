#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>

#include <random>
#include <omp.h>
#include <concurrent_vector.h>

#include "c_poly2freesrf.h"
#include "TriToPoint.hpp"

using namespace concurrency;

class CtrLpnt
{
public:
	int error;
	double ctrlpnt[48];
	CtrLpnt()
	{
		error = 0;
	}
	CtrLpnt(double pnt[48])
	{
		error = 0;
#pragma omp parallel for
		for (int i = 0; i < 48; i++) ctrlpnt[i] = pnt[i];
	}
};
class MinMax3d
{
public:
	double min[3];
	double max[3];

	MinMax3d(double min_[3], double max_[3])
	{
		min[0] = min_[0];
		min[1] = min_[1];
		min[2] = min_[2];
		max[0] = max_[0];
		max[1] = max_[1];
		max[2] = max_[2];
	}
};

class SamplePoint
{
public:
	double error_dist;
	int error;
	double p[3];
	SamplePoint() {
		error = 0; error_dist = -99999.0;
	}
	SamplePoint(double q[3])
	{
		p[0] = q[0];
		p[1] = q[1];
		p[2] = q[2];
		error = 0;
		error_dist = -99999.0;
	}
};

static FILE* DumpSrfStart(char* filename)
{
	FILE* fp = fopen(filename, "w");
	if ( fp == NULL ){
		return NULL;
	}
	fprintf(fp, "pconst_num(10,10);\n");

	return fp;
}
static void DumpSrfEnd(FILE* fp)
{
	if ( fp == NULL ){
		return;
	}
	fprintf(fp, "end;\n");
	fclose(fp);
}

static void DumpSrf(FILE*fp, int id, int un, int vn, double* cp)
{	


	int i, j, k;
/*
	k = 0;
	for ( i = 0; i < 3*vn+1; i++ ){
		for ( j = 0; j < 3*un+1; j++ ){
			fprintf(fp, "Pcp=[%f,%f,%f];\n", cp[3*k], cp[3*k+1], cp[3*k+2]);
			k++;
		}
	}
*/
	k = 0;
	fprintf(fp, "wavfsave(\"mesh_%d.obj\");\n", id);
	fprintf(fp, "Sbez.%d = %d %d(\n", id, un, vn);
	for ( i = 0; i < 3*vn+1; i++ ){
		for ( j = 0; j < 3*un+1; j++ ){
			fprintf(fp, "[%f,%f,%f]", cp[3*k], cp[3*k+1], cp[3*k+2]);
			if ( i < 3*vn || j < 3*un){
				fprintf(fp, ",\n");
			}else{
				fprintf(fp, "\n");
			}
			k++;
		}
	}

	fprintf(fp, ");\n");
	fprintf(fp, "uvmesh(20,20);\n");
	fprintf(fp, "etrgmesh(Sbez.%d);\n", id);
	fprintf(fp, "wavfsave();\n");
}

void ExportBezier0(int error, FILE* fp, double cp[4][4][3], int divnum, int* vcount)
{
	double r = (double)rand() / (double)RAND_MAX;
	double g = (double)rand() / (double)RAND_MAX;
	double b = (double)rand() / (double)RAND_MAX;

	r = (r > 0.8) ? 0.8 : r;
	g = (g > 0.8) ? 0.8 : g;
	b = (b > 0.8) ? 0.8 : b;
	if (error != 0)
	{
		r = 1.0;
		g = 1.0;
		b = 1.0;
	}
	double u, v;
	double n = divnum;
	for (u = 0.0; u <= 1.0; u += 1.0 / n)
	{
		for (v = 0.0; v <= 1.0; v += 1.0 / n)
		{
			double uv1[2] = { u, v };
			double uv2[2] = { u + 1.0 / n, v };
			double uv3[2] = { u + 1.0 / n, v + 1.0 / n };
			double uv4[2] = { u, v + 1.0 / n };
			double p[4][3];
			SrBez(cp, uv1, p[0]);
			SrBez(cp, uv2, p[1]);
			SrBez(cp, uv3, p[2]);
			SrBez(cp, uv4, p[3]);
			for (int k = 0; k < 4; k++)
			{
				fprintf(fp, "v %.3f %.3f %.3f %.2f %.2f %.2f\n", p[k][0], p[k][1], p[k][2], r, g, b);
			}
			fprintf(fp, "f %d %d %d\n", *vcount, *vcount + 2, *vcount + 1);
			fprintf(fp, "f %d %d %d\n", *vcount + 2, *vcount, *vcount + 3);
			*vcount += 4;
		}
	}
}

void ExportBezier(char* fname, double cp[4][4][3], int divnum)
{
	int vcount = 1;
	FILE* fp = fopen(fname, "w");
	ExportBezier0(0, fp, cp, divnum, &vcount);
	fclose(fp);
}

int vtx_count = 1;
FILE* fp_patch;
double Tol = -1.0;

int PatchFit(double rectA[2][2], int chkgrid[3], double* ctrlpnt, int upatch, int vpatch, int vtxnum, int trinum, double* vertex, double* uv, int* tri, double* normal, double* start, double* step, int depth);
int PatchFit0(double rectA[2][2], int chkgrid[3], double ctrlpnt[4][48], int upatch, int vpatch, int vtxnum, int trinum, double* vertex, double* uv, int* tri, double* normal, double* start, double* step, int depth);

//std::vector<CtrLpnt> PatchList;
concurrent_vector<CtrLpnt> PatchList;

int main_bezierPatch(char* file, char* out, int un, int vn);
int main_bezierPatch_BSpline(char* file, char* out, int un, int vn, int type);

int main(int argc, char** argv)
{
	if (argc < 2)
	{
		printf("%s meshFile -o CPFile [-t (0..2) -p N M]\n");
		printf("type 0 BezierPatch initial patch=MxN\n");
		printf("type 1 BezierPatch MxN\n");
		printf("type 2 One B-SplineSurface MxN points\n");
	}

	char* infile = NULL;
	char* outfile = NULL;
	int type = 2;
	int upatch = 0;
	int vpatch = 0;

	infile = argv[1];
	for (int i = 2; i < argc; i++)
	{
		if (strcmp(argv[i], "-o") == 0)
		{
			outfile = argv[i + 1];
			i++;
			continue;
		}
		if (strcmp(argv[i], "-tol") == 0)
		{
			Tol = atof(argv[i + 1]);
			i++;
			continue;
		}
		if (strcmp(argv[i], "-t") == 0)
		{
			type = atoi(argv[i + 1]);
			i++;
			continue;
		}
		if (strcmp(argv[i], "-p") == 0)
		{
			upatch = atoi(argv[i + 1]);
			vpatch = atoi(argv[i + 2]);
			i++;
			i++;
			continue;
		}
	}
	if (outfile == NULL)
	{
		printf("error output file ?\n");
		return -1;
	}
	if (type < 0 || type >= 5)
	{
		printf("error type = 0,1,2..4!!\n");
		return -2;
	}

	if (type == 0 && (upatch == 0 || vpatch == 0))
	{
		upatch = 5;
		vpatch = 5;
	}
	if (type == 1 && (upatch == 0 || vpatch == 0))
	{
		upatch = 55;
		vpatch = 55;
	}
	if (type >= 2 && (upatch == 0 || vpatch == 0))
	{
		upatch = 85;
		vpatch = 75;
	}

	printf("in :%s\n", infile);
	printf("out:%s\n", outfile);
	printf("type:%d\n", type);
	printf("upatch:%d\n", upatch);
	printf("vpatch:%d\n", vpatch);
	if (type == 0)
	{
		main_bezierPatch(infile, outfile, upatch, vpatch);
	}
	if (type == 1 || type >= 2)
	{
		main_bezierPatch_BSpline(infile, outfile, upatch, vpatch, type);
	}

	return 0;
}

double MinMax[2][3];
double Tlength;
void getMinMax(int vtxnum, int trinum, double* vertex)
{
	MinMax[0][0] = 99999.0;
	MinMax[1][0] = 99999.0;
	MinMax[2][0] = 99999.0;
	MinMax[0][1] = -99999.0;
	MinMax[1][1] = -99999.0;
	MinMax[2][1] = -99999.0;

	for (int i = 0; i < vtxnum; i++)
	{
		if (vertex[3 * i + 0] < MinMax[0][0]) MinMax[0][0] = vertex[3 * i];
		if (vertex[3 * i + 1] < MinMax[1][0]) MinMax[1][0] = vertex[3 * i + 1];
		if (vertex[3 * i + 2] < MinMax[2][0]) MinMax[2][0] = vertex[3 * i + 2];

		if (vertex[3 * i + 0] > MinMax[0][1]) MinMax[0][1] = vertex[3 * i];
		if (vertex[3 * i + 1] > MinMax[1][1]) MinMax[1][1] = vertex[3 * i + 1];
		if (vertex[3 * i + 2] > MinMax[2][1]) MinMax[2][1] = vertex[3 * i + 2];
	}

	Tlength =
		(MinMax[0][1] - MinMax[0][0])*(MinMax[0][1] - MinMax[0][0]) +
		(MinMax[1][1] - MinMax[1][0])*(MinMax[1][1] - MinMax[1][0]) +
		(MinMax[2][1] - MinMax[2][0])*(MinMax[2][1] - MinMax[2][0]);
	Tlength = sqrt(Tlength);

	printf("MinMax (%f %f %f) - (%f %f %f)\n",
		MinMax[0][0], MinMax[0][1], MinMax[0][2],
		MinMax[1][0], MinMax[1][1], MinMax[1][2]);
	printf("T=%f\n", Tlength);
}

int main_bezierPatch_BSpline(char* file, char* out, int un, int vn, int type)
{
	int vtxnum;
	int trinum;
	double* vertex;
	double* uv;
	double* normal;
	int* tri;
	double* step;
	double* start;
	int upatch;
	int vpatch;
	int chkgrid[3];
	int grid[3];
	double* ctrlpnt;
	int stat;
	FILE* fp;

	ReadMeshOBJUV(file, &vtxnum, &trinum, &vertex, &uv, &tri);
	getMinMax(vtxnum, trinum, vertex);

	normal = NULL;
	start = NULL;
	step = NULL;

	chkgrid[0] = 20;
	chkgrid[1] = 20;
	chkgrid[2] = 2;

	upatch = un;
	vpatch = vn;

	grid[0] = 3 * upatch + 1;
	grid[1] = 3 * vpatch + 1;
	grid[2] = 1;

	ctrlpnt = new double[grid[0] * grid[1] * 3];

	fp = DumpSrfStart("sbezer.txt");

	double rectA[2][2];

	rectA[0][0] = 0.0;
	rectA[0][1] = 1.0;
	rectA[1][0] = 0.0;
	rectA[1][1] = 1.0;

	stat = Polygon2SbezFitC(chkgrid, ctrlpnt, upatch, vpatch, rectA,
		vtxnum, trinum, vertex, uv, tri, normal, start, step, NULL);
	DumpSrf(fp, 0, upatch, vpatch, ctrlpnt);
	DumpSrfEnd(fp);

	char fname[256];
	sprintf(fname, out);
	fp_patch = fopen(fname, "w");

	if (fp_patch == NULL)
	{
		printf("output[%s] open error.\n");
		return -1;
	}

	fprintf(fp_patch, "%d\n", type);
	fprintf(fp_patch, "%d %d\n", upatch, vpatch);

	if (type >= 2)
	{
		int k = 0;
		for (int j = 0; j < grid[1]; j++)
		{
			for (int i = 0; i < grid[0]; i++)
			{
				fprintf(fp_patch, "%.3f %.3f %.3f\n", ctrlpnt[3 * k], ctrlpnt[3 * k + 1], ctrlpnt[3 * k + 2]);
				k++;
			}
		}
		return 0;
	}

	//Bezier Patch
	for (int iv = 1; iv <= vpatch; iv++)
	{
		for (int iu = 1; iu <= upatch; iu++)
		{
			double pcp[16 * 3];
			/*   Get Patch's Control   */
			int ucp = 9 * upatch + 3;
			int icp = 3 * ucp * (iv - 1) + 9 * iu - 9;
			for (int i = 0; i<4; i++) 
			{
				for (int j = 0; j<4; j++) 
				{
					*(pcp + 12 * i + 3 * j) = *(ctrlpnt + icp + 3 * j);
					*(pcp + 12 * i + 3 * j + 1) = *(ctrlpnt + icp + 3 * j + 1);
					*(pcp + 12 * i + 3 * j + 2) = *(ctrlpnt + icp + 3 * j + 2);
				}
				icp += ucp;
			}

			double cp[4][4][3];
			Conv16_3PtrTo4_4_3Ary(pcp, cp);
#if 10
			double cp_tmp[4][4][3];
			//パッチ内のフィッティング
			stat = SBezFitCP0(cp, cp_tmp);

			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					fprintf(fp_patch, "%.3f %.3f %.3f\n", cp_tmp[i][j][0], cp_tmp[i][j][1], cp_tmp[i][j][2]);
				}
			}
#else
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					fprintf(fp_patch, "%.3f %.3f %.3f\n", cp[i][j][0], cp[i][j][1], cp[i][j][2]);
				}
			}
#endif
		}
	}
	return 0;
}

//std::vector<MinMax3d> MinMax3dList;
concurrent_vector<MinMax3d> MinMax3dList;

int main_bezierPatch(char* file, char* out, int un, int vn)
{
	int vtxnum;
	int trinum;
	double* vertex;
	double* uv;
	double* normal;
	int* tri;
	double* step;
	double* start;
	int upatch;
	int vpatch;
	int chkgrid[3];
	int grid[3];
	double ctrlpnt[4][48];
    int stat;
	FILE* fp;

 	ReadMeshOBJUV(file, &vtxnum, &trinum, &vertex, &uv, &tri );
	getMinMax(vtxnum, trinum, vertex);


	normal = NULL;
	start = NULL;
	step = NULL;

	chkgrid[0] = 20;
	chkgrid[1] = 20;
	chkgrid[2] = 2;

	upatch = 1;
	vpatch = 1;

	grid[0] = 3*upatch+1;
	grid[1] = 3*vpatch+1;
	grid[2] = 1;

	fp = DumpSrfStart("sbezer.txt");

	double rectA[2][2];

	rectA[0][0] = 0.0;
	rectA[0][1] = 1.0;
	rectA[1][0] = 0.0;
	rectA[1][1] = 1.0;

	int depth = 0;

	if (MinMax3dList.size() == 0)
	{
#pragma omp for
		for (int k = 0; k < trinum; k++)
		{
			int vid[3];
			vid[0] = tri[3 * k];
			vid[1] = tri[3 * k + 1];
			vid[2] = tri[3 * k + 2];

			double minvtx[3] = { 99.0, 99.0, 99.0 };
			double maxvtx[3] = { -99.0, -99.0, -99.0 };
			double minu = 99.0, minv = 99.0;
			double maxu = -99.0, maxv = -99.0;
			double vtx[3][3], vtx_uv[3][2];
			for (int ki = 0; ki < 3; ki++)
			{
				vtx[ki][0] = vertex[3 * vid[ki]];
				vtx[ki][1] = vertex[3 * vid[ki] + 1];
				vtx[ki][2] = vertex[3 * vid[ki] + 2];

				vtx_uv[ki][0] = uv[3 * vid[ki]];
				vtx_uv[ki][1] = uv[3 * vid[ki] + 1];

				if (vtx[ki][0] < minvtx[0]) minvtx[0] = vtx[ki][0];
				if (vtx[ki][1] < minvtx[1]) minvtx[1] = vtx[ki][1];
				if (vtx[ki][2] < minvtx[2]) minvtx[2] = vtx[ki][2];

				if (vtx[ki][0] > maxvtx[0]) maxvtx[0] = vtx[ki][0];
				if (vtx[ki][1] > maxvtx[1]) maxvtx[1] = vtx[ki][1];
				if (vtx[ki][2] > maxvtx[2]) maxvtx[2] = vtx[ki][2];

				if (vtx_uv[ki][0] < minu) minu = vtx_uv[ki][0];
				if (vtx_uv[ki][1] < minv) minv = vtx_uv[ki][1];
				if (vtx_uv[ki][0] > maxu) maxu = vtx_uv[ki][0];
				if (vtx_uv[ki][1] > maxv) maxv = vtx_uv[ki][1];
			}
			;
			MinMax3dList.push_back(MinMax3d(minvtx, maxvtx));
		}
	}

#if 0
	upatch = 1;
	vpatch = 1;
	PatchFit0(rectA, chkgrid, ctrlpnt, upatch, vpatch, vtxnum, trinum, vertex, uv, tri, normal, start, step, depth);
#else	
	int n = std::max(un,vn);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			rectA[0][0] = (double)i / (double)n;
			rectA[0][1] = (double)(i + 1) / (double)n;
			rectA[1][0] = (double)j / (double)n;
			rectA[1][1] = (double)(j + 1) / (double)n;
			{
				double ctrlpnt_wrk[4][48];
				stat = PatchFit0(rectA, chkgrid, ctrlpnt_wrk, upatch, vpatch, vtxnum, trinum, vertex, uv, tri, normal, start, step, depth);
			}
		}
	}
#endif

	int np = 0;
	char fname[256];
#if 0
	sprintf(fname, "aaa_%04d.obj", np);
	fp_patch = fopen(fname, "w");
	for (int i = 0; i < PatchList.size(); i++)
	{
		if (i && i % 200 == 0)
		{
			np++;
			fclose(fp_patch);
			sprintf(fname, "aaa_%04d.obj", np);
			fp_patch = fopen(fname, "w");
			vtx_count = 1;
		}
		double cp[4][4][3];
		Conv16_3PtrTo4_4_3Ary(PatchList[i].ctrlpnt, cp);

		ExportBezier0(PatchList[i].error, fp_patch, cp, 80, &vtx_count);
	}
	fclose(fp_patch);
#endif

	sprintf(fname, out);
	fp_patch = fopen(fname, "w");
	if (fp_patch == NULL)
	{
		printf("output[%s] open error.\n");
		return -1;
	}

	fprintf(fp_patch, "%d\n", 0);
	fprintf(fp_patch, "%d %d\n", 1, PatchList.size());
	for (int i = 0; i < PatchList.size(); i++)
	{
		double cp[4][4][3];
		Conv16_3PtrTo4_4_3Ary(PatchList[i].ctrlpnt, cp);
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				fprintf(fp_patch, "%.3f %.3f %.3f\n", cp[i][j][0], cp[i][j][1], cp[i][j][2]);
			}
		}
	}
	fclose(fp_patch);

	printf("Patch=%d\n", PatchList.size());
	return 0;
}


int NumPatch = 0;
int PatchFit00(int i, int j, double rectA[2][2], int chkgrid[3], double ctrlpnt[4][48], int upatch, int vpatch, int vtxnum, int trinum, double* vertex, double* uv, int* tri, double* normal, double* start, double* step, int depth)
{
	const int n = 2;

	const double deltau = (rectA[0][1] - rectA[0][0]) / (double)n;
	const double deltav = (rectA[1][1] - rectA[1][0]) / (double)n;
	
	double rect2A[2][2];
	
	rect2A[0][0] = rectA[0][0] + (double)i*deltau;
	rect2A[0][1] = rectA[0][0] + (double)(i + 1)*deltau;

	rect2A[1][0] = rectA[1][0] + (double)j*deltav;
	rect2A[1][1] = rectA[1][0] + (double)(j + 1)*deltav;

	int stat = PatchFit(rect2A, chkgrid, ctrlpnt[2 * i + j], upatch, vpatch, vtxnum, trinum, vertex, uv, tri, normal, start, step, depth);
	int error = stat;
	if (depth > 10)
	{
		printf("@@@@@delta %f %f\n", (rect2A[0][1] - rect2A[0][0]), (rect2A[1][1] - rect2A[1][0]));
		stat = 0;
	}
	if (stat != 0)
	{
		double ctrlpnt_wrk[4][48];
		stat = PatchFit0(rect2A, chkgrid, ctrlpnt_wrk, upatch, vpatch, vtxnum, trinum, vertex, uv, tri, normal, start, step, depth);
	}
	else
	{
		CtrLpnt& p = CtrLpnt(ctrlpnt[2 * i + j]);
		p.error = error;
		PatchList.push_back(p);
		//double cp[4][4][3];
		//Conv16_3PtrTo4_4_3Ary(ctrlpnt[2 * i + j], cp);

		//ExportBezier0(fp_patch, cp, 30, &vtx_count);
		//NumPatch++;
		//printf("%d NumPatch=%d\n", depth, NumPatch);
	}
	return stat;
}

int PatchFit0(double rectA[2][2], int chkgrid[3], double ctrlpnt[4][48], int upatch, int vpatch, int vtxnum, int trinum, double* vertex, double* uv, int* tri, double* normal, double* start, double* step, int depth)
{
	depth += 1;

	const int n = 2;

	const double deltau = (rectA[0][1] - rectA[0][0]) / (double)n;
	const double deltav = (rectA[1][1] - rectA[1][0]) / (double)n;

#if 10
	double ctrlpnt2[4][4][48];
#pragma omp parallel sections
	{
#pragma omp section
		{
			PatchFit00(0, 0, rectA, chkgrid, ctrlpnt2[0], upatch, vpatch, vtxnum, trinum, vertex, uv, tri, normal, start, step, depth);
		}
#pragma omp section
		{
			PatchFit00(0, 1, rectA, chkgrid, ctrlpnt2[1], upatch, vpatch, vtxnum, trinum, vertex, uv, tri, normal, start, step, depth);
		}
#pragma omp section
		{
			PatchFit00(1, 0, rectA, chkgrid, ctrlpnt2[2], upatch, vpatch, vtxnum, trinum, vertex, uv, tri, normal, start, step, depth);
		}
#pragma omp section
		{
			PatchFit00(1, 1, rectA, chkgrid, ctrlpnt2[3], upatch, vpatch, vtxnum, trinum, vertex, uv, tri, normal, start, step, depth);
		}
	}
#else
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			double rect2A[2][2];
			rect2A[0][0] = rectA[0][0] + (double)i*deltau;
			rect2A[0][1] = rectA[0][0] + (double)(i + 1)*deltau;

			rect2A[1][0] = rectA[1][0] + (double)j*deltav;
			rect2A[1][1] = rectA[1][0] + (double)(j + 1)*deltav;

			int stat = PatchFit(rect2A, chkgrid, ctrlpnt[2 * i + j], upatch, vpatch, vtxnum, trinum, vertex, uv, tri, normal, start, step, depth);
			if (depth > 8)
			{
				printf("@@@@@delta %f %f\n", (rect2A[0][1] - rect2A[0][0]), (rect2A[1][1] - rect2A[1][0]));
				stat = 0;
			}
			if (stat != 0)
			{
				double ctrlpnt_wrk[4][48];
				stat = PatchFit0(rect2A, chkgrid, ctrlpnt_wrk, upatch, vpatch, vtxnum, trinum, vertex, uv, tri, normal, start, step, depth);
			}
			else
			{
				PatchList.push_back(CtrLpnt(ctrlpnt[2 * i + j]));
				//double cp[4][4][3];
				//Conv16_3PtrTo4_4_3Ary(ctrlpnt[2 * i + j], cp);

				//ExportBezier0(fp_patch, cp, 30, &vtx_count);
				NumPatch++;
				printf("%d NumPatch=%d\n", depth, NumPatch);
			}
		}
	}
#endif

	return 0;
}



int PatchFit(double rectA[2][2], int chkgrid[3], double* ctrlpnt, int upatch, int vpatch, int vtxnum, int trinum, double* vertex, double* uv, int* tri, double* normal, double* start, double* step, int depth)
{
	double defalut_tol = 0.0015;

	if (Tol > 0.0) defalut_tol = Tol;

	double tol = defalut_tol * pow(0.9, (double)(depth + 1));
	//tol = 0.01;
	tol = tol*Tlength/300.0;
	if (Tol < 0.0)
	{
		if (tol > 0.1) tol = 0.1;
	}
	if (tol < 1.0e-6) tol = 1.0e-6;

	int stat = 0;
	stat = Polygon2SbezFitC(chkgrid, ctrlpnt, upatch, vpatch, rectA,
		vtxnum, trinum, vertex, uv, tri, normal, start, step, NULL);

	double deltau = (rectA[0][1] - rectA[0][0]);
	double deltav = (rectA[1][1] - rectA[1][0]);


	double dist_max = -1.0;
	double m = 9;
	double n = 8;
	double deltau_s = (rectA[0][1] - rectA[0][0]) / (n - 1.0);
	double deltav_s = (rectA[1][1] - rectA[1][0]) / (n - 1.0);


	std::mt19937 mt;
	mt = std::mt19937(1);
	std::uniform_real_distribution<double> u_rand(0.0, 1.0);
	std::uniform_real_distribution<double> v_rand(0.0, 1.0);

	std::vector<SamplePoint> SamplePointList;
	SamplePointList.resize(2*n+2*n+m*m);

#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < (int)n; i++)
		{
			//外周 V一定線
			for (int j = 0; j < 2; j++)
			{
				double cp[4][4][3];
				Conv16_3PtrTo4_4_3Ary(ctrlpnt, cp);

				double uv1[2];
				uv1[0] = ((double)i)*deltau_s;
				uv1[1] = (j == 0) ? -0.0000001 : 1.0000001;
				SrBez(cp, uv1, SamplePointList[2 * i + j].p);
			}
			//外周 U一定線
			for (int j = 0; j < 2; j++)
			{
				double cp[4][4][3];
				Conv16_3PtrTo4_4_3Ary(ctrlpnt, cp);

				double uv1[2];
				uv1[0] = (j == 0) ? -0.0000001 : 1.0000001;	//<==  j で正しい
				uv1[1] = ((double)i)*deltav_s;	//<==	i で正しい
				SrBez(cp, uv1, SamplePointList[2 * n + 2 * i + j].p);
			}
		}

		deltau_s = (rectA[0][1] - rectA[0][0]) / (m - 1.0);
		deltav_s = (rectA[1][1] - rectA[1][0]) / (m - 1.0);

#pragma omp for
		for (int i = 0; i < (int)m; i++)
		{
			for (int j = 0; j < (int)m; j++)
			{
				double cp[4][4][3];
				Conv16_3PtrTo4_4_3Ary(ctrlpnt, cp);

				double uv1[2];
				//uv1[0] = u_rand(mt);
				//uv1[1] = v_rand(mt);
				uv1[0] = (double)i*deltau_s;
				uv1[1] = (double)j*deltav_s;
				SrBez(cp, uv1, SamplePointList[2 * n + 2 * n + m*i + j].p);
			}
		}
		//printf("%d == %d\n", (int)(2 * n + 2 * n + m), SamplePointList.size());
	}

	const int sz = SamplePointList.size();
#pragma omp parallel for
	for (int i = 0; i < sz; i++)
	{
		double dist_min = 9999.0;
		{
			double* p = SamplePointList[i].p;
			const int thn = omp_get_max_threads();
			std::vector<double> dist_min_th(thn + 1, 9999.0);
			for (int k = 0; k < trinum; k++)
			{
				const int th_id = omp_get_thread_num();
				const MinMax3d& minmax = MinMax3dList[k];
				if (p[0] < minmax.min[0] - tol) continue;
				if (p[1] < minmax.min[1] - tol) continue;
				if (p[2] < minmax.min[2] - tol) continue;
				if (p[0] > minmax.max[0] + tol) continue;
				if (p[1] > minmax.max[1] + tol) continue;
				if (p[2] > minmax.max[2] + tol) continue;

				double vtx[3][3];
				int vid[3];
				vid[0] = tri[3 * k];
				vid[1] = tri[3 * k + 1];
				vid[2] = tri[3 * k + 2];
				for (int ki = 0; ki < 3; ki++)
				{
					vtx[ki][0] = vertex[3 * vid[ki]];
					vtx[ki][1] = vertex[3 * vid[ki] + 1];
					vtx[ki][2] = vertex[3 * vid[ki] + 2];
				}


				double wknpnt[3];
				double dist;
				//int stat = TriToPoint(3, p, vtx[0], vtx[1], vtx[2], 0, wknpnt, &dist);
				int stat = TriToPoint( p, vtx[0], vtx[1], vtx[2], wknpnt, dist);
				if (stat != 0) {
					continue;
				}
				if (dist < dist_min_th[th_id]) {
					dist_min_th[th_id] = dist;
				}
			}

			for (int ki = 0; ki < thn; ki++)
			{
				if (dist_min_th[ki] < dist_min) dist_min = dist_min_th[ki];
			}
		}
		SamplePointList[i].error_dist = dist_min;
	}

	int on_edge_error = 0;
	double dist_s = 0.0;
	for (int i = 0; i < SamplePointList.size(); i++)
	{
		if (i < 4 * n)
		{
			if (SamplePointList[i].error_dist > tol)
			{
				on_edge_error++;
			}
		}
		dist_s += SamplePointList[i].error_dist;
		if (dist_max < SamplePointList[i].error_dist)
		{
			dist_max = SamplePointList[i].error_dist;
		}
	}
	dist_s /= (double)SamplePointList.size();

	if (dist_max < tol)
	{
		return 0;
	}

	printf("%d dist max=%.3f(%.3f) (%f,%f)-(%f,%f)\n", depth, dist_max, tol, rectA[0][0], rectA[1][0], rectA[0][1], rectA[1][1]);
	printf("dist s=%.3f delta %f %f\n", dist_s, deltau, deltav);

	if (deltau < 1.0e-6 || deltav < 1.0e-6)
	{
		return 0;
	}
	if (on_edge_error) return -1;

	if (dist_s > tol)
	{
		return -1;
	}
	if (dist_s > tol*1.1)
	{
		return -1;
	}
	//if (dist_max*dist_max > tol*tol)
	//{
	//	return -1;
	//}

	//double cp[4][4][3];
	//Conv16_3PtrTo4_4_3Ary(ctrlpnt, cp);

	//ExportBezier0(fp_patch, cp, 10, &vtx_count);

	return 0;
}