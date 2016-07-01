/* Molecular Surface Package */
/* General Math Routines */

/* Written by Michael L. Connolly */
/* January 7, 2002 */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

/* low level mathematical functions */


/* double-precision routines */


/* change vector length to 1 (i.e. normalize) */
int normalize (double vect[3])
{
	int k;
	double n;

	n = norm (vect);
	if (n <= 0.0) {
		return (0);
	}

	for (k = 0; k < 3; k++)
		vect[k] /= n;
	return (1);
}

/* return norm of vector */
double norm(double f[3])
{
	double n;

	n = norm_squared (f);
	if (n <= 0.0)
		return (0.0);
	n = sqrt (n);
	return (n);
}

/* return square of norm of vector */
double norm_squared (double f[3])
{
	int k;
	double n;

	n = 0.0;
	for (k = 0; k < 3; k++)
		n += f[k] * f[k];
	return (n);
}

/* return dot product of two vectors */
double dot_product (double f[3], double g[3])
{
	int k;
	double n;

	n = 0.0;
	for (k = 0; k < 3; k++)
		n += f[k] * g[k];
	return (n);
}

/* return squared distance between points */
double distance_squared (double f[3], double g[3])
{
	int k;
	double fg[3];

	for (k = 0; k < 3; k++)
		fg[k] = f[k] - g[k];
	return (norm_squared (fg));
}

int negdot (double vect1[3], double vect2[3])
{
	int isneg;
	double dt;

	dt = dot_product (vect1, vect2);
	isneg = ((int) (dt < 0.0));
	return (isneg);
}

/* reverse direction of vector */
void reverse (double vect[3])
{
	int k;

	for (k = 0; k < 3; k++) vect[k] = (-vect[k]);
}

/* return distance between two points */
double distance (double f[3], double g[3])
{
	int k;
	double fg[3];

	for (k = 0; k < 3; k++)
		fg[k] = f[k] - g[k];
	return (norm (fg));
}

/* cross product of two vectors */
void cross (double f[3], double g[3], double h[3])
{
	h[0] = f[1] * g[2] - f[2] * g[1];
	h[1] = f[2] * g[0] - f[0] * g[2];
	h[2] = f[0] * g[1] - f[1] * g[0];
}

/* triple product of three vectors */
double triple_product (double f[3], double g[3], double h[3])
{
	double t;
	double fg[3];

	cross (f, g, fg);
	t = h[0] * fg[0] + h[1] * fg[1] + h[2] * fg[2];
	return (t);
}

/* return area of triangle */
double triangle_area (double p[3], double q[3], double r[3])
{
	int k;
	double a;
	double u[3], v[3], uv[3];

	for (k = 0; k < 3; k++) {
		u[k] = q[k] - p[k];
		v[k] = r[k] - p[k];
	}
	cross (u, v, uv);
	a = norm (uv) / 2.0;
	return (a);
}

/* return signed area of triangle */
double signed_triangle_area (double p[3], double q[3], double r[3], double axis[3])
{
	int k;
	double a, dt;
	double u[3], v[3], uv[3];

	for (k = 0; k < 3; k++) {
		u[k] = q[k] - p[k];
		v[k] = r[k] - p[k];
	}
	cross (u, v, uv);
	dt = dot_product (uv, axis);
	a = norm (uv) / 2.0;
	if (dt < 0.0) a = (-a);
	return (a);
}

/* volume of tetrahedron */
double tetrahedron_volume (double p[3], double q[3], double r[3], double s[3])
{
	int k;
	double t;
	double u[3], v[3], w[3];

	for (k = 0; k < 3; k++) {
		u[k] = q[k] - p[k];
		v[k] = r[k] - p[k];
		w[k] = s[k] - p[k];
	}
	t = (-triple_product (u, v, w)) / 6.0;
	return (t);
}

/* set up axis given triangle */
/* solid angle of tetrahedron seen by first vertex */

double tetra_solid (double v1[3], double v2[3], double v3[3], double v4[3])
{
	int k;
	double sum, ang1, ang2, ang3, angle;
	double vect1[3], vect2[3], vect3[3];
	double axis1[3], axis2[3], axis3[3];

	/* compute vectors from first vertex to other three */

	for (k = 0; k < 3; k++) {
		vect1[k] = v2[k] - v1[k];
		vect2[k] = v3[k] - v1[k];
		vect3[k] = v4[k] - v1[k];
	}
	if (!normalize (vect1)) {
		return (0.0);
	}
	if (!normalize (vect2)) {
		return (0.0);
	}
	if (!normalize (vect3)) {
		return (0.0);
	}

	/* compute normals to three faces */
	cross (vect1, vect2, axis1);
	cross (vect2, vect3, axis2);
	cross (vect3, vect1, axis3);
	if (!normalize (axis1)) {
		return (0.0);
	}
	if (!normalize (axis2)) {
		return (0.0);
	}
	if (!normalize (axis3)) {
		return (0.0);
	}

	/* compute interior angles from dihedral angles */

	/* later: check fourth argument */
	ang1 = PI - odd_angle (axis3, axis1, vect1, (double) 1.0);
	ang2 = PI - odd_angle (axis1, axis2, vect2, (double) 1.0);
	ang3 = PI - odd_angle (axis2, axis3, vect3, (double) 1.0);

	sum = ang1 + ang2 + ang3;

	/* subtract the pi to get excess
	   (use 5 pi for negative orientations) */

	if (sum >= 3 * PI) angle = (sum - 5 * PI);
	else angle = (sum - PI);
	return (angle);
}

int clockwise (double p1[2], double p2[2], double p3[2])
{
	int k;
	double v12[2], v13[2];
	double crs;
	
	for (k = 0; k < 2; k++) {
		v12[k] = p2[k] - p1[k];
		v13[k] = p3[k] - p1[k];
	}
	crs = v12[0] * v13[1] - v12[1] * v13[0];
	return (crs <= 0.0);
}

int pit (double pnt[2], int up, double t1[2], double t2[2], double t3[2])
{
	if (up == clockwise (pnt, t1, t2)) return (0);
	if (up == clockwise (pnt, t2, t3)) return (0);
	if (up == clockwise (pnt, t3, t1)) return (0);
	return (1);
}

int piq (double pnt[2], int up, double t1[2], double t2[2], double t3[2], double t4[2])
{
	if (up == clockwise (pnt, t1, t2)) return (0);
	if (up == clockwise (pnt, t2, t3)) return (0);
	if (up == clockwise (pnt, t3, t4)) return (0);
	if (up == clockwise (pnt, t4, t1)) return (0);
	return (1);
}

int lines_intersect (double p1[2], double p2[2], double p3[2], double p4[2])
{
	int k;
	for (k = 0; k < 2; k++) {
		if (p1[k] == p3[k]) return (0);
		if (p1[k] == p4[k]) return (0);
		if (p2[k] == p3[k]) return (0);
		if (p2[k] == p4[k]) return (0);
	}
	if (clockwise (p1, p3, p4) == clockwise (p2, p3, p4)) return (0);
	if (clockwise (p3, p1, p2) == clockwise (p4, p1, p2)) return (0);
	return (1);
}

/* return the fraction of the distance
   from p1 to p2 that the intersection point is */
double line_intersection (double p1[2], double p2[2], double p3[2], double p4[2])
{
	int k;
	double nperp, d13, d12, ratio;
	double v12[2], v13[2], v34[2];
	double perp[2];
	
	for (k = 0; k < 2; k++) {
		v34[k] = p4[k] - p3[k];
		v13[k] = p3[k] - p1[k];
		v12[k] = p2[k] - p1[k];
	}
	/* compute unit vector perpendicular to v34 */
	perp[0] = v34[1]; perp[1] = - v34[0];
	nperp = perp[0] * perp[0] + perp[1] * perp[1];
	if (nperp <= 0.0) return (0.5);
	nperp = sqrt (nperp);
	if (nperp <= 0.0) return (0.5);
	for (k = 0; k < 2; k++)
		perp[k] /= nperp;
	d13 = v13[0] * perp[0] + v13[1] * perp[1];
	d12 = v12[0] * perp[0] + v12[1] * perp[1];
	if (d12 == 0.0) return (0.5);
	ratio = d13/d12;
	return (ratio);
}

/* solid angle of tetrahedron seen by first vertex */

double tetra_solid_angle (double v1[3], double v2[3], double v3[3], double v4[3])
{
	int k, nnegative;
	double sum, omega;
	double ong1, ong2, ong3;
	double ing1, ing2, ing3;
	double vect1[3], vect2[3], vect3[3];
	double axis1[3], axis2[3], axis3[3];

	omega = 0.0;

	/* compute vectors from first vertex to other three */
	for (k = 0; k < 3; k++) {
		vect1[k] = v2[k] - v1[k];
		vect2[k] = v3[k] - v1[k];
		vect3[k] = v4[k] - v1[k];
	}
	if (!normalize (vect1)) {
		return (0.0);
	}
	if (!normalize (vect2)) {
		return (0.0);
	}
	if (!normalize (vect3)) {
		return (0.0);
	}

	/* compute normals to three faces */
	cross (vect1, vect2, axis1);
	cross (vect2, vect3, axis2);
	cross (vect3, vect1, axis3);
	if (!normalize (axis1)) {
		return (0.0);
	}
	if (!normalize (axis2)) {
		return (0.0);
	}
	if (!normalize (axis3)) {
		return (0.0);
	}

	/* compute interior angles from dihedral angles */

	nnegative = 0;
	ong1 = odd_angle (axis3, axis1, vect1, 1.0);
	if (ong1 < 0.0) nnegative++;
	if (ong1 >= 0.0) ing1 = PI - ong1;
	else if (ong1 < 0.0) ing1 = PI + ong1;
	ong2 = odd_angle (axis1, axis2, vect2, 1.0);
	if (ong2 < 0.0) nnegative++;
	if (ong2 >= 0.0) ing2 = PI - ong2;
	else if (ong2 < 0.0) ing2 = PI + ong2;
	ong3 = odd_angle (axis2, axis3, vect3, 1.0);
	if (ong3 < 0.0) nnegative++;
	if (ong3 >= 0.0) ing3 = PI - ong3;
	else if (ong3 < 0.0) ing3 = PI + ong3;
	/*
	pang1 = positive_angle (axis3, axis1, vect1);
	pang2 = positive_angle (axis1, axis2, vect2);
	pang3 = positive_angle (axis2, axis3, vect3);
	*/
	if (fabs (ong1) < 0.01) return (0.0);
	if (fabs (ong2) < 0.01) return (0.0);
	if (fabs (ong3) < 0.01) return (0.0);

	sum = (ing1 + ing2 + ing3);
	if (nnegative == 0) omega = sum - PI;
	else if (nnegative == 3) omega = PI - sum;
	else omega = 0.0;

	return (omega);
}

/* in range from minus pi to pi */
/* axis and v0, v1 are assumed to be of unit length */
double triple_angle (double v0[3], double v1[3], double axis[3])
{
	double t, d, angle;
	if (v1[0] == v0[0] && v1[1] == v0[1] && v1[2] == v0[2])
		return (0.0);
	d = dot_product (v0, v1);
	if (d < -1.0) d = (-1.0);
	else if (d > 1.0) d = 1.0;
	angle = acos (d);
	t = triple_product (v0, v1, axis);
	if (t < 0) angle =  (- angle);
	return (angle);
}


void setup_axis (double p[3], double q[3], double r[3], double a[3])
{
	int k;
	double nw;
	double u[3], v[3], w[3];

	for (k = 0; k < 3; k++) {
		u[k] = q[k] - p[k];
		v[k] = r[k] - p[k];
	}

	cross (u, v, w);
	nw = norm (w);
	if (nw <= 0.0) {
		return;
	}

	for (k = 0; k < 3; k++)
		a[k] = w[k] / nw;
}

double do_axis (double trico[3][3], double axis[3])
{
	int k;
	int i0, i1, i2, three;
	int okays[3];
	double area, areas[3], axes[3][3];
	double vect1[3], vect2[3];
	for (three = 0; three < 3; three++) {
		if (three == 0) {
			i0 = 0; i1 = 1; i2 = 2;
		}
		else if (three == 1) {
			i0 = 1; i1 = 2; i2 = 0;
		}
		else if (three == 2) {
			i0 = 2; i1 = 0; i2 = 1;
		}
		for (k = 0; k < 3; k++) {
			vect1[k] = trico[i1][k] - trico[i0][k];
			vect2[k] = trico[i2][k] - trico[i0][k];
		}
		cross (vect1, vect2, axes[three]);
		area = norm (axes[three]) / 2.0;
		areas[three] = fabs (area);
		okays[three] = normalize (axes[three]);
	}
	for (three = 0; three < 3; three++) {
		area = areas[three];
		if (okays[three] && area > 0.0) {
			for (k = 0; k < 3; k++)
				axis[k] = axes[three][k];
			return (area);
		}
	}
	return (0.0);
}

/* return angle between vectors */
double simple_angle (double vect1[3], double vect2[3])
{
	double dt, norm1, norm2, ang;

	norm1 = norm (vect1);
	if (norm1 <= 0.0) {
		return (0.0);
	}
	norm2 = norm (vect2);
	if (norm2 <= 0.0) {
		return (0.0);
	}

	dt = dot_product (vect1, vect2);
	dt /= (norm1 * norm2);
	if (dt < -1.0) dt = -1.0;
	if (dt > 1.0) dt = 1.0;
	ang = acos (dt);
	return (ang);
}


/* return angle between vectors */
double positive_angle (double vect1[3], double vect2[3], double axis[3])
{
	double dt, norm1, norm2, norm12, trip, ang;

	norm1 = norm (vect1);
	if (norm1 <= 0.0) {
		return (0.0);
	}
	norm2 = norm (vect2);
	if (norm2 <= 0.0) {
		return (0.0);
	}
	norm12 = norm1 * norm2;
	if (norm12 <= 0.0) {
		return (0.0);
	}

	dt = dot_product (vect1, vect2);
	dt /= (norm12);
	if (dt < -1.0) dt = -1.0;
	if (dt > 1.0) dt = 1.0;
	ang = acos (dt);
	trip = triple_product (vect1, vect2, axis);
	if (trip < 0.0) ang = 2 * PI - ang;
	return (ang);
}

/* return angle between vectors */
double odd_angle (double vect1[3], double vect2[3], double axis[3], double sgn)
{
	double ang;
	double axis2[3];
	int k;

	for (k = 0; k < 3; k++)
		axis2[k] = sgn * axis[k];
	ang = positive_angle (vect1, vect2, axis2);
	if (ang > PI) ang -= 2 * PI;
	return (ang);
}

/* arbitrary perpendicular */
void arbprp (double given[3], double perp[3])
{
	int k;
	double dott;
	double try[3];

	/* my own algorithm */
	try[0] = given[1] * given[1] + given[2] * given[2];
	try[1] = given[0] * given[0] + given[2] * given[2];
	try[2] = given[0] * given[0] + given[1] * given[1];

	if (!normalize (try)) return;

	/* subtract parallel component of first try */
	dott = dot_product (given, try);
	for (k = 0; k < 3; k++)
		try[k] = try[k] - dott * given[k];

	if (!normalize (try)) return;

	/* check correctness */
	if (fabs (dot_product (given, try)) > 0.001) {
		return;
	}

	for (k = 0; k < 3; k++)
		perp[k] = try[k];
}

/* create rotation matrix given axis and angle */
int generate_rotation (double direction[3], double angle, double rotation_matrix[3][3])
{
	int j, k;
	double a[3][3], b[3][3], z[3][3], za[3][3];

	if (norm(direction) <= 0.0) return (0);
	for (k = 0; k < 3; k++)
	        a[2][k] = direction[k];
	arbprp (a[2], a[0]);
	cross (a[2], a[0], a[1]);
	if (norm(a[1]) <= 0.0) return (0);
	for (j = 0; j < 3; j++)
		for (k = 0; k < 3; k++)
      	   b[j][k] = a[k][j];
	for (j = 0; j < 3; j++)
		for (k = 0; k < 3; k++)
		 z[j][k] = (j == k);
	z[0][0] = cos (angle);
	z[1][1] = cos (angle);
	z[0][1] = -sin (angle);
	z[1][0] = sin (angle);
	for (j = 0; j < 3; j++)
	for (k = 0; k < 3; k++)
		za[j][k] = z[j][0] * a[0][k] + z[j][1] * a[1][k] + z[j][2] * a[2][k];
    for (j = 0; j < 3; j++)
        for (k = 0; k < 3; k++)
            rotation_matrix[j][k] =
                     b[j][0] * za[0][k] + b[j][1] * za[1][k] + b[j][2] * za[2][k];
	return (1);
}

double quad_angle (v1, v2, v3, v4)
double v1[3], v2[3], v3[3], v4[3];
{
	int k;
	double qa;
	double vect23[3], vect24[3], vect41[3], vect12[3];
	double axis[3];

	for (k = 0; k < 3; k++) {
		vect23[k] = v3[k] - v2[k];
		vect24[k] = v4[k] - v2[k];
		vect41[k] = v1[k] - v4[k];
		vect12[k] = v2[k] - v1[k];
	}
	if (!normalize (vect23)) return (0.0);
	if (!normalize (vect24)) return (0.0);
	if (!normalize (vect41)) return (0.0);
	if (!normalize (vect12)) return (0.0);
	cross (vect23, vect24, axis);
	if (norm (axis) <= 0.0) return (0.0);  /* kludge */
	if (!normalize (axis)) return (0.0);
	qa = odd_angle (vect41, vect12, axis, (double) 1.0);
	return (qa);
}

void make_matrix (double angle, double matrix[3][3])
{
	int j, k;
	double cos_ang, sin_ang, ang;

	ang = angle;

	cos_ang = cos (ang);
	sin_ang = sin (ang);
	for (j = 0; j < 3; j++)
		for (k = 0; k < 3; k++)
			matrix[j][k] = (j == k);
	matrix[0][0] = cos_ang;
	matrix[0][2] = sin_ang;
	matrix[2][0] = -sin_ang;
	matrix[2][2] = cos_ang;
	
}

/* transform vector */

void trvec (double rotation[3][3], double vec[3])
{
	int k;
	double out[3];

	for (k = 0; k < 3; k++)
		out[k] = rotation[k][0] * vec[0]
			   + rotation[k][1] * vec[1]
			   + rotation[k][2] * vec[2];

	for (k = 0; k < 3; k++)
		vec[k] = out[k];
}

/* transform point */

void trpnt (double rotation[3][3], double center[3], double translate[3], double pnt[3])
{
	int k;
	double tmp1[3], tmp2[3];

	for (k = 0; k < 3; k++)
		tmp1[k] = pnt[k] - center[k];

	for (k = 0; k < 3; k++)
		tmp2[k] = rotation[k][0] * tmp1[0]
		        + rotation[k][1] * tmp1[1]
				+ rotation[k][2] * tmp1[2];

	for (k = 0; k < 3; k++)
		pnt[k] = tmp2[k] + center[k] + translate[k];
}

int pclipped (double clip_center[3], double clip_axis[3], double pnt[3])
{
	int k;
	double center_to_pnt[3];
	double dt;

	for (k = 0; k < 3; k++)
		center_to_pnt[k] = pnt[k] - clip_center[k];
	dt = dot_product (center_to_pnt, clip_axis);
	if (dt > 0.0) return (1);

	return (0);
}

/* find the z coordinate of the intersection of the z axis through
 * the point and the clipping plane
 */
double zclipped (double clip_center[3], double clip_axis[3], double pnt[3])
{
	int k;
	double nv, alpha, beta, sin_alpha, sin_beta, d;
	double z[3], v[3], vz[3], w[3], q[3];

	z[0] = 0.0;
	z[1] = 0.0;
	z[2] = 1.0;
	for (k = 0; k < 3; k++)
		v[k] = clip_center[k] - pnt[k];
	nv = norm (v);
	if (!normalize(v)) {
		return (clip_center[2]);
	}
	cross (v, z, vz);
	if (!normalize(vz)) {
		return (clip_center[2]);
	}
	cross (vz, clip_axis, w);
	if (!normalize (w)) {
		return (clip_center[2]);
	}
	alpha = dot_product (v, z);
	beta = dot_product (w, z);
	sin_alpha = sin (alpha);
	sin_beta = sin (beta);
	if (sin_beta == 0.0) {
		return (clip_center[2]);
	}
	d = (sin_alpha / sin_beta) * nv;
	for (k = 0; k < 3; k++)
		q[k] = clip_center[k] + d * w[k];
	return (q[2]);
}

/*
	MSP

	Molecular Surface Package
	Copyright 1986, 1989 by Michael L. Connolly
	All rights reserved

*/

