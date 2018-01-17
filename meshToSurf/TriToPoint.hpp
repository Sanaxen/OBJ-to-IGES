#include <math.h>

inline void sub_vector(const double a[3], const double b[3], double c[3])
{
	c[0] = a[0] - b[0];
	c[1] = a[1] - b[1];
	c[2] = a[2] - b[2];
}

inline void cross_product(const double vl[3], const double vr[3], double c[3])
{
	c[0] = vl[1] * vr[2] - vl[2] * vr[1];
	c[1] = vl[2] * vr[0] - vl[0] * vr[2];
	c[2] = vl[0] * vr[1] - vl[1] * vr[0];
}

inline double dot_product(const double vl[3], const double  vr[3])
{
	return vl[0] * vr[0] + vl[1] * vr[1] + vl[2] * vr[2];
}

// Triangle and point hit determination
inline int hittest_point_polygon_3d(double A[3], double B[3], double C[3], double P[3])
{

	double AB[3], BP[3], BC[3], CP[3], CA[3], AP[3];
	sub_vector(B, A, AB);
	sub_vector(P, B, BP);

	sub_vector(C, B, BC);
	sub_vector(P, C, CP);

	sub_vector(A, C, CA);
	sub_vector(P, A, AP);

	double c1[3], c2[3], c3[3];
	cross_product(AB, BP, c1);
	cross_product(BC, CP, c2);
	cross_product(CA, AP, c3);

	double dot_12 = dot_product(c1, c2);
	double dot_13 = dot_product(c1, c3);

	if (dot_12 > 0.0 && dot_13 > 0.0) 
	{
		return 0;
	}

	return 1;
}

inline double DistPointToPoint(double p1[3], double p2[3])
{
	double v[3];
	sub_vector(p1, p2, v);

	return sqrt(dot_product(v, v));
}


//Find the closest point on the line from point P
inline int NearPosOnLine(double P[3], double A[3], double B[3], double NearPoint[3])
{
	double AB[3], AP[3];

	AB[0] = B[0] - A[0];
	AB[1] = B[1] - A[1];
	AB[2] = B[2] - A[2];
	AP[0] = P[0] - A[0];
	AP[1] = P[1] - A[1];
	AP[2] = P[2] - A[2];

	double len = sqrt(dot_product(AB, AB));
	if (fabs(len) < 1.0E-16)
	{
		return -1;
	}

	double nAB[3];
	nAB[0] = AB[0] / len;
	nAB[1] = AB[1] / len;
	nAB[2] = AB[2] / len;

	double dist_AX = dot_product(nAB, AP);

	NearPoint[0] = A[0] + (nAB[0] * dist_AX);
	NearPoint[1] = A[1] + (nAB[1] * dist_AX);
	NearPoint[2] = A[2] + (nAB[2] * dist_AX);

	if (dist_AX >= 0.0 && dist_AX <= 1.0)
	{
		return 0;
	}
	return 1;
}

inline double clamp(double x, double minv, double maxv)
{
	if (x < minv) return minv;
	if (x > maxv) return maxv;
	return x;
}

inline int closesPointOnTriangle(const double triangle[3][3], const double sourcePosition[3], double trianglePosition[3])
{
	double edge0[3];
	sub_vector(triangle[1], triangle[0], edge0);

	double edge1[3];
	sub_vector(triangle[2], triangle[0], edge1);

	double v0[3];
	sub_vector(triangle[0], sourcePosition, v0);

	const double a = dot_product(edge0, edge0);
	const double b = dot_product(edge0, edge1);
	const double c = dot_product(edge1, edge1);
	const double d = dot_product(edge0, v0);
	const double e = dot_product(edge1, v0);

	const double det = a*c - b*b;
	double s = b*e - c*d;
	double t = b*d - a*e;

	int stat = 1;
	if (s + t <= 1.0 && s >= 0.0 && t >= 0.0)
	{
		stat = 0;
	}

	if (s + t < det)
	{
		if (s < 0.0)
		{
			if (t < 0.0)
			{
				if (d < 0.0)
				{
					s = clamp(-d / a, 0.0, 1.0);
					t = 0.0;
				}
				else
				{
					s = 0.0;
					t = clamp(-e / c, 0.0, 1.0);
				}
			}
			else
			{
				s = 0.0;
				t = clamp(-e / c, 0.0, 1.0);
			}
		}
		else if (t < 0.0)
		{
			s = clamp(-d / a, 0.0, 1.0);
			t = 0.0;
		}
		else
		{
			float invDet = 1.0 / det;
			s *= invDet;
			t *= invDet;
			if (s + t <= 1.0 && s >= 0.0 && t >= 0.0)
			{
				stat = 0;
			}
		}
	}
	else
	{
		if (s < 0.0)
		{
			const double tmp0 = b + d;
			const double tmp1 = c + e;
			if (tmp1 > tmp0)
			{
				const double numer = tmp1 - tmp0;
				const double denom = a - 2 * b + c;
				s = clamp(numer / denom, 0.0, 1.0);
				t = 1 - s;
			}
			else
			{
				t = clamp(-e / c, 0.0, 1.0);
				s = 0.0;
			}
		}
		else if (t < 0.0)
		{
			if (a + d > b + e)
			{
				const double numer = c + e - b - d;
				const double denom = a - 2 * b + c;
				s = clamp(numer / denom, 0.0, 1.0);
				t = 1 - s;
			}
			else
			{
				s = clamp(-e / c, 0.0, 1.0);
				t = 0.0;
			}
		}
		else
		{
			const double numer = c + e - b - d;
			const double denom = a - 2 * b + c;
			s = clamp(numer / denom, 0.0, 1.0);
			t = 1.0 - s;
			if (s + t <= 1.0 && s >= 0.0 && t >= 0.0)
			{
				stat = 0;
			}
		}
	}
	trianglePosition[0] = triangle[0][0] + s * edge0[0] + t * edge1[0];
	trianglePosition[1] = triangle[0][1] + s * edge0[1] + t * edge1[1];
	trianglePosition[2] = triangle[0][2] + s * edge0[2] + t * edge1[2];

	return stat;
}


inline int	TriToPoint(double pnt[3], double vtx1[3], double vtx2[3], double vtx3[3], double NearPos[3], double& dist)
{
	double triangle[3][3];
	triangle[0][0] = vtx1[0];
	triangle[0][1] = vtx1[1];
	triangle[0][2] = vtx1[2];
	triangle[1][0] = vtx2[0];
	triangle[1][1] = vtx2[1];
	triangle[1][2] = vtx2[2];
	triangle[2][0] = vtx3[0];
	triangle[2][1] = vtx3[1];
	triangle[2][2] = vtx3[2];

	int stat = closesPointOnTriangle(triangle, pnt, NearPos);
	dist = DistPointToPoint(NearPos, pnt);

	//stat = hittest_point_polygon_3d(vtx1, vtx2, vtx3, NearPos);

#if 0
	if (stat != 0)
	{
		double d;
		double nearP[3];
		if (NearPosOnLine(pnt, vtx1, vtx2, nearP) == 0)
		{
			d = DistPointToPoint(nearP, pnt);
			if (d < dist)
			{
				dist = d;
				stat = 0;
			}
		}
		else
			if (NearPosOnLine(pnt, vtx2, vtx3, nearP) == 0)
			{
				d = DistPointToPoint(nearP, pnt);
				if (d < dist)
				{
					dist = d;
					stat = 0;
				}
			}
			else
				if (NearPosOnLine(pnt, vtx3, vtx1, nearP) == 0)
				{
					d = DistPointToPoint(nearP, pnt);
					if (d < dist)
					{
						dist = d;
						stat = 0;
					}
				}
	}

	if (stat != 0)
	{
		double d = DistPointToPoint(vtx1, pnt);
		if (d < dist)
		{
			dist = d;
			stat = 0;
		}
		else
		{
			d = DistPointToPoint(vtx2, pnt);
			if (d < dist)
			{
				dist = d;
				stat = 0;
			}
			else
			{
				d = DistPointToPoint(vtx3, pnt);
				if (d < dist)
				{
					dist = d;
					stat = 0;
				}
			}
		}
	}
#endif
	return 0;
}
