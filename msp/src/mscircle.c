/*
	Molecular Surface Package
	Copyright 1986, 1989 by Michael L. Connolly
	All rights reserved
	November 28, 2001
*/

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"


/* return angle of arc (in range of 0 to 2 pi) */

double arc_ang (struct arc *a)
{
	int k;
	double ang, angle;
	double vect1[3], vect2[3];
	struct circle *cir;
	struct vertex *v1, *v2;

	cir = a -> cir;
	if (cir -> radius <= 0.0) return (0.0);
	v1 = a -> vtx[0];
	v2 = a -> vtx[1];
	if (a -> shape == STRAIGHT)
		return (0.0);
	/* check for complete circle */
	if (v1 == NULL || v2 == NULL) {
		angle = 2 * PI;
		return (angle);
	}
	/* check for arc beginning and ending at same vertex */
	if (v1 == v2) {
		angle = 2 * PI;
		return (angle);
	}
	/* compute unit vectors from circle center to vertices */
	for (k = 0; k < 3; k++) {
		vect1[k] = v1 -> center[k] - cir -> center[k];
		vect2[k] = v2 -> center[k] - cir -> center[k];
	}
	/* call positive angle function to do the work */
	ang = positive_angle (vect1, vect2, cir -> axis);
	return (ang);
}

/* get tangent to arc */
int get_tangent (struct arc *arc_ptr, int orn, int which_vertex, double vector[3])
{
	struct vertex *vtx;
	struct circle *cir;
	double outward[3];
	int k;
    struct cept *ex;

	if (arc_ptr == NULL) {
        ex = new_cept (PARAMETER_ERROR, NULL_VALUE, FATAL_SEVERITY);
        add_function (ex, "get_tangent");
		add_source (ex, "mscircle.c");
		return(0);
	}
	if (which_vertex > 1 || which_vertex < 0) {
        ex = new_cept (PARAMETER_ERROR, INVALID_VALUE, FATAL_SEVERITY);
        add_function (ex, "get_tangent");
		add_source (ex, "mscircle.c");
        add_long (ex, "which vertex", (long) which_vertex);
		return(0);
	}
	vtx = arc_ptr -> vtx[which_vertex];
	if (vtx == NULL) {
        ex = new_cept (PARAMETER_ERROR, NULL_VALUE, FATAL_SEVERITY);
        add_function (ex, "get_tangent");
		add_source (ex, "mscircle.c");
        add_message (ex, "arc vertex");
		return(0);
	}
	cir = arc_ptr -> cir;
	if (cir == NULL) {
        ex = new_cept (PARAMETER_ERROR, NULL_VALUE, FATAL_SEVERITY);
        add_function (ex, "get_tangent");
        add_message (ex, "circle that the arc lies on");
		return(0);
	}
	for (k = 0; k < 3; k++)
		outward[k] = vtx -> center[k] - cir -> center[k];
	if (orn) cross (outward, cir -> axis, vector);
	else cross (cir -> axis, outward, vector);
	if (!normalize (vector)) {
        ex = new_cept (GEOMETRY_ERROR, NULL_VALUE, FATAL_SEVERITY);
        add_function (ex, "get_tangent");
		add_source (ex, "mscircle.c");
		add_message(ex, "normalize tangent fails");
		return(0);
	}
	return (1);
}

/* theta value of circle (related to geodesic curvature) */
double theta_circle (struct circle *cir, int orn, struct variety *vty, double fsgn)
{
	int k;
	double sgn, displacement, angle, radius;
	double dc_vector[3];

	sgn = 1 - 2 * orn;
	radius = cir -> radius;
	for (k = 0; k < 3; k++)
		dc_vector[k] = cir -> center[k] - vty -> center[k];
	displacement = dot_product (dc_vector, cir -> axis) * sgn * (-1.0) * fsgn;
	angle = atan2 (displacement, radius);
	return (angle);
}


/* return length of arc */

double arc_length (struct arc *a)
{
	double ang, l;
	struct circle *cir;
	struct vertex *v1, *v2;
    struct cept *ex;

	/* move info to local variable */
	cir = a -> cir;
	v1 = a -> vtx[0];
	v2 = a -> vtx[1];
	if (a -> shape == STRAIGHT) {
		if (v1 == NULL && v2 == NULL) {
			ex = new_cept (PARAMETER_ERROR, NULL_VALUE, FATAL_SEVERITY);
			add_function (ex, "arc_length");
			add_source (ex, "mscircle.c");
            add_message (ex, "line segment missing both vertices");
			return (0.0);
		}
		if (v1 == NULL || v2 == NULL) {
			ex = new_cept (PARAMETER_ERROR, NULL_VALUE, FATAL_SEVERITY);
			add_function (ex, "arc_length");
			add_source (ex, "mscircle.c");
            add_message (ex, "line segment missing one vertex");
			return (0.0);
		}
		l = distance (v1 -> center, v2 -> center);
		return (l);
	}
	ang = arc_ang (a);
	l = ang * cir -> radius;
	return (l);
}

/* compute midpoint of arc */
int middle (struct arc *a, double returned_point[3])
{
	int k;
	double dh;
	double half[3], vect[3];
	struct circle *cir;
	struct vertex *vtx1, *vtx2;
    struct cept *ex;

	if (a -> shape == STRAIGHT) {
        ex = new_cept (PARAMETER_ERROR, INVALID_VALUE, FATAL_SEVERITY);
        add_function (ex, "middle");
		add_source (ex, "mscircle.c");
        add_message (ex, "arc shape is straight, should be curved");
		return (0);
	}
	cir = a -> cir;
	vtx1 = a -> vtx[0];
	vtx2 = a -> vtx[1];
	if (vtx1 == NULL && vtx2 == NULL) {
        ex = new_cept (GEOMETRY_ERROR, NULL_VALUE, FATAL_SEVERITY);
        add_function (ex, "middle");
		add_source (ex, "mscircle.c");
        add_message (ex, "both vertices are null");
		return (0);
	}
	if (vtx1 == NULL || vtx2 == NULL) {
        ex = new_cept (GEOMETRY_ERROR, NULL_VALUE, FATAL_SEVERITY);
        add_function (ex, "middle");
		add_source (ex, "mscircle.c");
        add_message (ex, "one vertex is null");
		return (0);
	}
	for (k = 0; k < 3; k++)
		half[k] = (vtx1 -> center[k] + vtx2 -> center[k]) / 2;
	dh = distance (half, cir -> center);
	if (dh <= 0.0) {
        ex = new_cept (GEOMETRY_ERROR, DEGENERACY, FATAL_SEVERITY);
        add_function (ex, "middle");
		add_source (ex, "mscircle.c");
		return (0);
	}
	for (k = 0; k < 3; k++)
		vect[k] = (cir -> radius / dh) * (half[k] - cir -> center[k]);
	for (k = 0; k < 3; k++)
		returned_point[k] = cir -> center[k] + vect[k];
	return (1);
}

/* circle-circle intersection */
int circle_circle (struct circle *cir1, struct circle *cir2, double xpnt[2][3])
{
	int k;
	double l, h, s, d2, r;
	double ccvect[3], lvect[3], svect[3], mco[3];

	/* calculate circle-circle vector */
	for (k = 0; k < 3; k++)
		ccvect[k] = cir2 -> center[k] - cir1 -> center[k];
	h = dot_product (ccvect, cir1 -> axis);
	/* compute line of plane-plane intersection */
	cross (cir2 -> axis, cir1 -> axis, lvect);
	l = norm (lvect);
	if (l <= 0.0)
		return (0);		/* no intersection: planes parallel */
	for (k = 0; k < 3; k++)
		lvect[k] /= l;
	/* compute vector from center of first circle to line */
	cross (cir2 -> axis, lvect, svect);
	s = norm (svect);
	if (s <= 0.0)
		return (0);
	for (k = 0; k < 3; k++)
		svect[k] /= s;
	/* compute coordinates of midpoint (center of s0) */
	for (k = 0; k < 3; k++)
		mco[k] = cir2 -> center[k] + (h / l) * svect[k];
	/* compute radius of s0 */
	d2 = distance_squared (mco, cir1 -> center);
	r = cir1 -> radius * cir1 -> radius - d2;
	if (r <= 0.0)
		return (0);		/* line misses first circle, no intersection */
	/* intersection */
	r = sqrt (r);
	/* compute coordinates of two intersection points */
	for (k = 0; k < 3; k++) {
		xpnt[0][k] = mco[k] - r * lvect[k];
		xpnt[1][k] = mco[k] + r * lvect[k];
	}
	return (1);
}

int arc_arc (struct arc *arc1, struct arc *arc2, double xpnt[2][3]) 	/* arc-arc intersection */
{
	int k, m;
	int x[2];
	struct circle *cir1, *cir2;
	struct vertex *vtx1, *vtx2, *vtx3, *vtx4;

	/* store pointers in local variables */
	cir1 = arc1 -> cir;
	cir2 = arc2 -> cir;
	if (cir1 == cir2) return (0);
	vtx1 = arc1 -> vtx[0];
	vtx2 = arc1 -> vtx[1];
	vtx3 = arc2 -> vtx[0];
	vtx4 = arc2 -> vtx[1];
	if (vtx1 != NULL && vtx1 == vtx3) return (0);
	if (vtx1 != NULL && vtx1 == vtx4) return (0);
	if (vtx2 != NULL && vtx2 == vtx3) return (0);
	if (vtx2 != NULL && vtx2 == vtx4) return (0);
	/* calculate circle-circle intersection */
	if (!circle_circle (cir1, cir2, xpnt))
		return (0);
	/* circles intersect */
	/* initialize to: points in both arcs */
	x[0] = 1;
	x[1] = 1;
	/* check which points lie within arcs */
	for (m = 0; m < 2; m++) {	/* point */
		if (!point_in_arc (xpnt[m], arc1))
			x[m] = 0;		/* failure for 1st arc */
		if (!point_in_arc (xpnt[m], arc2))
			x[m] = 0;		/* failure for 2nd arc */
	}
	/* return # of intersection points */
	if (!x[0] && !x[1])
		return (0);
	if (x[0] && x[1])
		return (2);
	if (x[0] && !x[1])
		return (1);
	/* x0 out and x1 in: shift coordinates */
	for (k = 0; k < 3; k++)
		xpnt[0][k] = xpnt[1][k];
	return (1);
}

int edg_edg (struct edge *edg1, struct edge *edg2, double xpnt[2][3]) 	/* arc-arc intersection */
{
	int k, m, orn1, orn2;
	int x[2];
	double sgn1, sgn2;
	double axis1[3], axis2[3];
	struct arc *arc1, *arc2;
	struct circle *cir1, *cir2;
	struct circle tcir1, tcir2;
	struct vertex *vtx1, *vtx2, *vtx3, *vtx4;

	/* store pointers in local variables */
	arc1 = edg1 -> arcptr;
	arc2 = edg2 -> arcptr;
	orn1 = edg1 -> orn;
	orn2 = edg2 -> orn;
	cir1 = arc1 -> cir;
	cir2 = arc2 -> cir;
	if (cir1 == cir2) return (0);
	sgn1 = 1 - 2 * orn1;
	sgn2 = 1 - 2 * orn2;
	for (k = 0; k < 3; k++) {
		axis1[k] = sgn1 * cir1 -> axis[k];
		axis2[k] = sgn2 * cir2 -> axis[k];
	}
	vtx1 = arc1 -> vtx[orn1];
	vtx2 = arc1 -> vtx[1-orn1];
	vtx3 = arc2 -> vtx[orn2];
	vtx4 = arc2 -> vtx[1-orn2];
	if (vtx1 != NULL && vtx1 == vtx3) return (0);
	if (vtx1 != NULL && vtx1 == vtx4) return (0);
	if (vtx2 != NULL && vtx2 == vtx3) return (0);
	if (vtx2 != NULL && vtx2 == vtx4) return (0);
	/* set up temporary structures */
	store_circle (&tcir1, cir1 -> center, cir1 -> radius, axis1);
	store_circle (&tcir2, cir2 -> center, cir2 -> radius, axis2);
	/* calculate circle-circle intersection */
	if (!circle_circle (&tcir1, &tcir2, xpnt))
		return (0);
	/* circles intersect */
	/* initialize to: points in both arcs */
	x[0] = 1;
	x[1] = 1;
	/* check which points lie within arcs */
	for (m = 0; m < 2; m++) {	/* point */
		if (!point_in_arc (xpnt[m], arc1))
			x[m] = 0;		/* failure for 1st arc */
		if (!point_in_arc (xpnt[m], arc2))
			x[m] = 0;		/* failure for 2nd arc */
	}
	/* return # of intersection points */
	if (!x[0] && !x[1])
		return (0);
	if (x[0] && x[1])
		return (2);
	if (x[0] && !x[1])
		return (1);
	/* x0 out and x1 in: shift coordinates */
	for (k = 0; k < 3; k++)
		xpnt[0][k] = xpnt[1][k];
	return (1);
}

int point_in_arc (double xpnt[3], struct arc *a)			/* is point in arc ? */
{
	int k;
	double ang, angx, r;
	double v0[3], vx[3], v1[3];
	struct circle *cir;
	struct vertex *vtx0, *vtx1;

	/* store pointers in local variables */
	vtx0 = a -> vtx[0];
	vtx1 = a -> vtx[1];
	cir = a -> cir;
	r = cir -> radius;
	/* point is always within boundaryless arc */
	if (vtx0 == NULL || vtx1 == NULL)
		return (1);
	/* calculate vectors from circle center */
	for (k = 0; k < 3; k++) {
		v0[k] = (vtx0 -> center[k] - cir -> center[k]) / r;
		v1[k] = (vtx1 -> center[k] - cir -> center[k]) / r;
		vx[k] = (xpnt[k] - cir -> center[k]) / r;
	}
	/* compute angles */
	ang = positive_angle (v0, v1, cir -> axis);
	angx = positive_angle (v0, vx, cir -> axis);
	/* compare angles */
	return ( (int) (angx <= ang));
}

/* sort points on circle */
void sort_points (double center[3], double radius, double axis[3], int n_points, double circle_points[][3], short orn[], short indices[])
{
	int k, index0, index1, base, n_e, edg_idx;
	int n_start, n_stop, pnt;
	int n_sorted;
	double minimum_angle;
	int point_used[100];
	double point_angle[100];
	double point_vector[100][3];
	char message[MAXLINE];
    struct cept *ex;

	/* should be even number of points */
	if (n_points % 2 != 0) {
        ex = new_cept (GEOMETRY_ERROR, INVALID_VALUE, FATAL_SEVERITY);
        add_function (ex, "sort_points");
        add_message (ex, "odd number of points on circle");
		return;
	}
	/* initialize numbers of starts and stops */
	n_start = 0;
	n_stop = 0;
	for (pnt = 0; pnt < n_points; ++pnt) {
		/* count start and stop vertices */
		if (orn[pnt]) n_stop++; else n_start++;
		/* compute vector from circle center to vertex */
		for (k = 0; k < 3; k++)
			point_vector[pnt][k] =
				(circle_points[pnt][k] - center[k]) / radius;
	}
	/* whatever starts must stop */
	if (n_start != n_stop) {
        ex = new_cept (GEOMETRY_ERROR, INCONSISTENCY, FATAL_SEVERITY);
        add_function (ex, "sort_points");
        add_long (ex, "number of start points", (long) n_start);
        add_long (ex, "number of stop points", (long) n_stop);
		return;
	}
	/* use first start vertex as base vertex for angle sort */
	base = -1;
	for (pnt = 0; pnt < n_points; pnt++)
		if (orn[pnt] == 0) {
			base = pnt;
			break;
		}
	/* check for failure to find base vertex */
	if (base < 0) {
        ex = new_cept (GEOMETRY_ERROR, INCONSISTENCY, FATAL_SEVERITY);
        add_function (ex, "sort_points");
        add_message(ex, "base point not found");
		return;
	}
	/* compute angles of vertex vectors relative to base vector */
	for (pnt = 0; pnt < n_points; pnt++)
		point_angle[pnt] =
			positive_angle (point_vector[base], point_vector[pnt], axis);
	point_angle[base] = 0.0;
	/* initialization */
	for (pnt = 0; pnt < n_points; pnt++)
		point_used[pnt] = 0;
	/* sort vertex list */
	n_e = n_points / 2;			/* number of edges */
	n_sorted = 0;
	for (edg_idx = 0; edg_idx < n_e; edg_idx++) {
		index0 = -1;
		/* find starting vertex */
		for (pnt = 0; pnt < n_points; pnt++) {
			if (point_used[pnt]) continue;		/* skip used vertices */
			if (orn[pnt] == 1) continue;	/* skip stop vertices */
			index0 = pnt;						/* got one */
			break;
		}
		/* check for failure to find starting vertex for arc */
		if (index0 < 0) {
			ex = new_cept (GEOMETRY_ERROR, INCONSISTENCY, FATAL_SEVERITY);
			add_function (ex, "sort_points");
			add_message (ex, "starting vertex not found");
			return;
		}
		/* the first time through we should get the base vertex */
		if (edg_idx == 0 && index0 != base) {
			ex = new_cept (GEOMETRY_ERROR, INCONSISTENCY, FATAL_SEVERITY);
            add_function (ex, "sort_points");
			add_message (ex, "first start point is not the base point");
			return;
		}
		point_used[index0] = 1;		/* mark vertex as used */
		minimum_angle = 2 * PI;		/* initialize to high value */
		index1 = -1;			/* initialize to illegal value */
		/* find ending vertex */
		for (pnt = 0; pnt < n_points; pnt++) {
			if (point_used[pnt]) continue;		/* skip used vertices */
			if (orn[pnt] == 0) continue;	/* skip start vertices */
			/* want minimum > start angle */
			if (point_angle[pnt] <= point_angle[index0]) continue;
			/* store if better than best so far */
			if (point_angle[pnt] < minimum_angle) {
				index1 = pnt;
				minimum_angle = point_angle[pnt];
			}
		}
		/* check for failure to find stop vertex */
		if (index1 < 0) {
			ex = new_cept (GEOMETRY_ERROR, INCONSISTENCY, FATAL_SEVERITY);
            add_function (ex, "sort_points");
			add_message (ex, "ending point not found");
			return;
		}
		point_used[index1] = 1;		/* mark vertex as used */
		/* store info in sorted array */
		if (n_sorted >= n_points - 1) {
			ex = new_cept (LOGIC_ERROR, MSOVERFLOW, FATAL_SEVERITY);
            add_function (ex, "sort_points");
            add_long (ex, "number of points", (long) n_points);
            add_long (ex, "number of sorted points", (long) n_sorted);
			return;
		}
		indices[n_sorted] = (short) index0;
		n_sorted++;
		indices[n_sorted] = (short) index1;
		n_sorted++;
	}
	if (n_sorted != n_points) {
		ex = new_cept (LOGIC_ERROR, MSUNDERFLOW, FATAL_SEVERITY);
        add_function (ex, "sort_points");
        add_source (ex, "mscircle.c");
        add_long (ex, "number of points", (long) n_points);
        add_long (ex, "number of sorted points", (long) n_sorted);
		return;
	}
}

/* compute angle between tangent vector to two edges */

double edge_delta (struct edge *edg1, struct edge *edg2, double normal_vector[3])
{
	int orn1, orn2, result;
	double ang;
	double tang1[3], tang2[3];

	if (edg1 -> arcptr -> cir == edg2 -> arcptr -> cir)
		return (0.0);
	/* transfer to local variables */
	orn1 = edg1 -> orn;
	orn2 = edg2 -> orn;
	result = get_tangent (edg1 -> arcptr, orn1, 1 - orn1, tang1);
	if (!result) return(0.0);
	result = get_tangent (edg2 -> arcptr, orn2, orn2, tang2);
	if (!result) return(0.0);
	ang = odd_angle (tang1, tang2, normal_vector, 1.0);
	return (ang);
}

/* Angle, ARC AND CIRCLE FUNCTIONS */


int arc_plane (struct arc *arc_ptr, double pcen[3], double paxis[3], double xpnts[2][3])
{
	int k, m;
	int x[2];
	int ncx;
	char message[MAXLINE];
	struct circle *cir;
    struct cept *ex;

	/* store pointer in local variables */
	cir = arc_ptr -> cir;

	/* calculate circle-plane intersection */
	ncx = circle_plane (cir, pcen, paxis, xpnts);
	if (error()) return(0);

	if (ncx <= 0) {
		return (0);
	}
	if (ncx != 2) {
        ex = new_cept (RETURN_ERROR, INVALID_VALUE, FATAL_SEVERITY);
        add_function (ex, "arc_plane");
        add_source (ex, "mscircle.c");
        add_long (ex, "number of intersection points", (long) ncx);
        add_message (ex, "should be 2");
		return(0);
	}

	/* circle intersects plane */

	/* initialize to: points in both arcs */
	x[0] = 1;
	x[1] = 1;

	/* check which points lie within arc */
	for (m = 0; m < 2; m++) {	/* point */
		if (!point_in_arc (xpnts[m], arc_ptr))
			x[m] = 0;			/* failure for arc */
	}

	/* return # of intersection points */
	if (!x[0] && !x[1])
		return (0);
	if (x[0] && x[1])
		return (2);
	if (x[0] && !x[1])
		return (1);

	/* x0 out and x1 in: shift coordinates */
	for (k = 0; k < 3; k++)
		xpnts[0][k] = xpnts[1][k];

	return (1);
}

int circle_plane (struct circle *cir, double c2[3], double a2[3], double xpnts[2][3])
{
	int k;
	double dot1;
	double r1, r3;
	double a1[3], c1[3];
	double a3[3], c3[3];
	double a13[3], c12[3], t123, s;
    struct cept *ex;

	for (k = 0; k < 3; k++) {
		c1[k] = cir -> center[k];
		a1[k] = cir -> axis[k];
	}

	r1 = cir -> radius;
	if (r1 <= 0.0) {
        ex = new_cept (GEOMETRY_ERROR, DEGENERACY, FATAL_SEVERITY);
        add_function (ex, "circle_plane");
        add_source (ex, "mscircle.c");
        add_double (ex, "radius", r1);
		return(0);
	}

	cross (a2, a1, a3);
	if (norm (a3) <= 0.0) {
        ex = new_cept (GEOMETRY_ERROR, DEGENERACY, FATAL_SEVERITY);
        add_function (ex, "circle_plane");
        add_source (ex, "mscircle.c");
		add_message (ex, "circle || plane");
		return (0);
	}
	normalize (a3);

	for (k = 0; k < 3; k++)
		c12[k] = c2[k] - c1[k];
	dot1 = dot_product (c12, a2);
	t123 = triple_product (a1, a3, a2);
	if (t123 == 0.0) {
        ex = new_cept (GEOMETRY_ERROR, DEGENERACY, FATAL_SEVERITY);
        add_function (ex, "circle_plane");
        add_source (ex, "mscircle.c");
		add_message (ex, "triple product");
		return(0);
	}
	s = dot1 / t123;
	r3 = r1 * r1 - s * s;
	if (r3 <= 0.0) return (0);

	r3 = sqrt (r3);
	cross (a1, a3, a13);
	for (k = 0; k < 3; k++)
		c3[k] = c1[k] + s * a13[k];

	for (k = 0; k < 3; k++) {
		xpnts[0][k] = c3[k] - r3 * a3[k];
		xpnts[1][k] = c3[k] + r3 * a3[k];
	}
	return (2);
}


/* MEASUREMENT */


/* compute circumference of cycle */

double circum (struct cycle *cyc)
{
	double c;
	struct edge *edg;
	struct arc *arc_ptr;
	struct circle *cir;

	if (cyc == NULL) return (0.0);

	c = 0.0;		/* initialize */

	for (edg = cyc -> first_edge; edg != NULL; edg = edg -> next) {
		arc_ptr = edg -> arcptr;
		cir = arc_ptr -> cir;
		c += arc_ang (arc_ptr) * cir -> radius;
	}
	return (c);
}

int edges_in_cycle (struct cycle *cyc)
{
	int n;
	struct edge *e;

	if (cyc == NULL) return (0);
	n = 0;
	for (e = cyc -> first_edge; e != NULL; e = e -> next)
		n++;
	return (n);
}

/* count the number of edges in the face */
int edges_in_face(struct face *fac)
{
	int n;
	struct cycle *cyc;

	if (fac == NULL) return (0);
	n = 0;
	for (cyc = fac -> first_cycle; cyc != NULL; cyc = cyc -> next)
		n += edges_in_cycle (cyc);
	return (n);
}



/* ARC SUBDIVISION */

/* subdivide arc */
/* this routine is concerned only with pointers, not geometry */
int subdivide_arc (struct surface *this_srf, struct arc *a)
{
	int n_point, i, shape;
	long lfn, ofn, ffn;
	double alpha;
	char message[MAXLINE];
	struct vertex **vertex_list;
	struct circle *cir;
	struct arc **arc_list, *arcptr;
	struct edge *next_edge, *previous_edge, *edg0, *edg1;
	struct edge **forward_list, **backward_list;
	struct face *fac0, *fac1;
    struct cept *ex;

	a -> subdivided = 1;				/* flag not to subdivide again */
	/* do not subdivide short arcs */
	if (a -> small) return(2);

	/* store info in local variables */
	edg0 = a -> edg[0];		/* old direct edge */
	edg1 = a -> edg[1];		/* old reversed edge */
	fac0 = edg0 -> fac;
	fac1 = edg1 -> fac;
	alpha = a -> alpha;
	shape = a -> shape;
	
	/* set up last face number to use in new arc call */
	lfn = a -> lfn;
	/* set up original face number to use in new arc call */
	ofn = a -> ofn;
	/* set up first face number to use in new arc call */
	ffn = a -> ffn;
	
	if (lfn == 0L) {
        ex = new_cept (PARAMETER_ERROR, NULL_VALUE, FATAL_SEVERITY);
        add_function (ex, "subdivide_arc");
        add_source (ex, "mscircle.c");
		add_double (ex, "last face number", lfn);
		return(0);
	}
	if (ffn == 0L) {
        ex = new_cept (PARAMETER_ERROR, NULL_VALUE, FATAL_SEVERITY);
        add_function (ex, "subdivide_arc");
        add_source (ex, "mscircle.c");
		add_double (ex, "first face number", ffn);
		return(0);
	}

	/* geometric subdivision of arc returns point coordinates */
	/* for straight arcs, ask for only a bisection */
	if (a -> shape == STRAIGHT)
		n_point = straight_subdivision (this_srf, a, 2, &vertex_list);
	else n_point = curved_subdivision (this_srf, a, &vertex_list);
	if (error()) return(0);

	/* check for end points only */
	if (n_point <= 2) {
		a -> small = 1;
		return(2);
	}


	/* circle pointer */
	cir = a -> cir;

	/* allocate local arrays */
	forward_list = (struct edge **)
		allocate_pointers (EDGE, n_point);
	if (forward_list == NULL) {
        ex = new_cept (MEMORY_ERROR, ALLOCATION, FATAL_SEVERITY);
        add_function (ex, "subdivide_arc");
        add_source (ex, "mscircle.c");
        add_long (ex, "n_point", n_point);
        add_object (ex, EDGE, "forward list");
		return(0);
	}

	backward_list = (struct edge **)
		allocate_pointers (EDGE, n_point);
	if (backward_list == NULL) {
        ex = new_cept (MEMORY_ERROR, ALLOCATION, FATAL_SEVERITY);
        add_function (ex, "subdivide_arc");
        add_source (ex, "mscircle.c");
        add_long (ex, "n_point", n_point);
        add_object (ex, EDGE, "backward list");
		return(0);
	}

	arc_list = (struct arc **)
		allocate_pointers (ARC, n_point);
	if (arc_list == NULL) {
        ex = new_cept (MEMORY_ERROR, ALLOCATION, FATAL_SEVERITY);
        add_function (ex, "subdivide_arc");
        add_source (ex, "mscircle.c");
        add_long (ex, "n_point", n_point);
        add_object (ex, ARC, "arc list");
		return(0);
	}

	/* save old info */
	next_edge = a -> edg[0] -> next;
	previous_edge = a -> edg[1] -> next;

	/* create new arcs and edges */

	for (i = 1; i < n_point - 1; i++) {
		arcptr = new_arc (cir, *(vertex_list+i), *(vertex_list+i+1), shape, 1, alpha, lfn, ofn, ffn);
		*(arc_list + i) = arcptr;
		link_arc (this_srf, arcptr);
		*(forward_list + i) = new_edge (*(arc_list+i), 0, fac0, NULL);
	}

	/* store old arc pointer at beginning of array */
	*arc_list = a;
	a -> small = 1;

	/* create new reverse edges */
	for (i = 0; i < n_point - 2; i++) {
		*(backward_list+i) = new_edge (*(arc_list+i), 1, fac1,NULL);
	}
	/* store old direct edge pointer at beginning of array */
	*forward_list = edg0;

	/* store old reverse edge pointer at end of array */
	*(backward_list + n_point - 2) = edg1;

	/* old reverse edge structure now points at last arc */
	edg1 -> arcptr = (*(arc_list + n_point - 2));

	/* and vice-versa */
	(*(arc_list+n_point - 2)) -> edg[1] = edg1;

	/* replace old ending vertex of original arc */
	a -> vtx[1] = *(vertex_list + 1);

	/* confirm or set old starting vertex of original arc */
	a -> vtx[0] = *(vertex_list + 0);

	/* set up links from edge to edge in both directions */

	for (i = 0; i < n_point - 2; i++)
		(*(forward_list+i)) -> next = *(forward_list + i + 1);

	for (i = 1; i < n_point - 1; i++)
		(*(backward_list+i)) -> next = *(backward_list + i - 1);

	/* links at ends */
	(*(forward_list + n_point - 2)) -> next = next_edge;
	(*backward_list) -> next = previous_edge;

	/* free temporary memory */
	free_pointers (VERTEX, vertex_list);
	free_pointers (EDGE, forward_list);
	free_pointers (ARC, arc_list);
	free_pointers (EDGE, backward_list);
	return(n_point);
}

/* geometric subdivision */

/* curved arc subdivision */

int curved_subdivision (struct surface *this_srf, struct arc *given_arc, struct vertex ***vertex_list_address)
{
	int i, n_division, n_point, without_ends;
	double alpha;
	double *verts;
	struct vertex *vtx;
	struct vertex **vertex_list;
    struct cept *ex;

	alpha = given_arc -> alpha;
	if (alpha <= 0.0) return (0);
	/* straight and small arcs are not subdivided */
	if (given_arc -> shape == STRAIGHT) return (0);
	if (given_arc -> small) return (0);
	/* if angle is less than input parameter, don't ever subdivide */
	if (!given_arc -> lune && arc_ang(given_arc) <= alpha) {
		given_arc -> small = 1;
		return (0);
	}
	n_point = count_sub_arc (given_arc, alpha);
	if (error()) return (0);
	n_division = n_point - 1;
	verts = geo_sub_arc (given_arc, n_point);
	if (error()) return (0);

	/* no end flag for full circles */
	without_ends = ( given_arc -> vtx[0] == NULL ||
		given_arc -> vtx[1] == NULL);

	/* allocate memory */
	vertex_list = (struct vertex **)
		allocate_pointers (VERTEX, n_point);
	if (vertex_list == NULL) {
        ex = new_cept (MEMORY_ERROR, ALLOCATION, FATAL_SEVERITY);
        add_function (ex, "curved_subdivision");
        add_source (ex, "mscircle.c");
        add_long (ex, "n_point", n_point);
        add_object (ex, VERTEX, "vertex list");
		return(0);
	}

	if (without_ends) {
		for (i = 0; i < n_division; i++) {
			vtx = new_vertex(verts+3*i, NULL, NULL, given_arc, NULL);
			if (vtx == NULL) return (0);
			link_vertex (this_srf, vtx);
			*(vertex_list+i) = vtx;
		}
		*(vertex_list + n_division) = *(vertex_list);
	}
	else {
		*vertex_list = given_arc -> vtx[0];
		*(vertex_list + n_division) = given_arc -> vtx[1];
		for (i = 1; i < n_division; i++) {
			vtx = new_vertex (verts+3*i, NULL, NULL, given_arc, NULL);
			if (vtx == NULL) return (0);
			link_vertex (this_srf, vtx);
			*(vertex_list+i) = vtx;
		}
	}
	/* store vertex array pointer in calling routine */
	*vertex_list_address = vertex_list;
	if (verts != NULL) free_doubles (verts, 0, VERTS);
	return (n_point);
}

/* geometric subdivision */

int count_sub_arc (struct arc *a, double alpha)
{
	double ang;
	int ndiv, npoint;

	/* determine angle subtended by arc */
	ang = arc_ang(a);
	if (error()) return(0);

	/* calculate number of subdivisions and number of points */
	ndiv = ang / alpha + 1;
	if (ndiv < 2) ndiv = 2;
	npoint = ndiv + 1;
	return (npoint);
}

double *geo_sub_arc (struct arc *a, int npoint)
{
	int ndiv, noend, i, k;
	double ang, del, theta, ctheta, stheta;
	double base[3], zenith[3], vcen[3];
	double *verts;
	struct circle *cir;
    struct cept *ex;

	ndiv = npoint - 1;
	/* determine angle subtended by arc */
	ang = arc_ang(a);
	if (error()) return(0);

	/* allocate memory */
	verts = allocate_doubles (3 * npoint, 0, VERTS);
	if (verts == NULL) {
        ex = new_cept (MEMORY_ERROR, ALLOCATION, FATAL_SEVERITY);
        add_function (ex, "geo_sub_arc");
        add_source (ex, "mscircle.c");
        add_long (ex, "npoint", npoint);
        add_variable (ex, VERTS, "verts");
		return(0);
	}

	/* angle subtended by each new, short arc */
	del = ang / ndiv;

	/* pointer to circle */
	cir = a -> cir;
	/* no end flag for full circles */
	noend = ( a -> vtx[0] == NULL || a -> vtx[1] == NULL);

	/* calculate orthonormal frame */

	/* base vector */
	if (noend) arbprp (cir -> axis, base);
	else
	for (k = 0; k < 3; k++)
		base[k] = (a -> vtx[0] -> center[k] - cir -> center[k])
			/ cir -> radius;

	/* zenith at right angles to base, in plane of circle */
	cross (cir -> axis, base, zenith);

	/* two cases: 360 degrees and arc with endpoints */
	if (noend) {
		theta = 0.0;	/* initialize angle */
		for (i = 0; i < ndiv; i++, theta += del) {
			/* trigonometry */
			ctheta = cos (theta);
			stheta = sin (theta);
			for (k = 0; k < 3; k++)
				vcen[k] = cir -> center[k] + cir -> radius *
					(ctheta * base[k] + stheta * zenith[k]);
			for (k = 0; k < 3; k++)
				*(verts + 3 * i + k) = vcen[k];
		}
		for (k = 0; k < 3; k++)
			*(verts + 3 * ndiv + k) = *(verts+k);
	}
	else {
		for (k = 0; k < 3; k++) {
			*(verts + k) = a -> vtx[0] -> center[k];
			*(verts + 3 * ndiv + k) = a -> vtx[1] -> center[k];
		}
		theta = del;		/* initialize angle */
		for (i = 1; i < ndiv; i++, theta += del) {
			/* trigonometry */
			ctheta = cos (theta);
			stheta = sin (theta);
			for (k = 0; k < 3; k++)
				vcen[k] = cir -> center[k] + cir -> radius *
					(ctheta * base[k] + stheta * zenith[k]);
			for (k = 0; k < 3; k++)
				*(verts+3*i+k) = vcen[k];
		}
	}
	return (verts);
}



/* straight arc subdivision - pass number of pieces */

int straight_subdivision (struct surface *this_srf, struct arc *given_arc, int n_pieces, struct vertex ***vertex_list_address)
{
	int i, k, n_division, n_point;
	double len, d, increment;
	double base[3], vertex_center[3];
	struct vertex **vertex_list;
	struct vertex *vtx0, *vtx1, *vtx;
    struct cept *ex;
	
	/* curved and small arcs are not subdivided */
	if (given_arc -> shape != STRAIGHT) return (0);
	if (given_arc -> small) return (0);
	if (n_pieces < 2) return (0);

	vtx0 = given_arc -> vtx[0];
	vtx1 = given_arc -> vtx[1];

	if (vtx0 == NULL && vtx1 == NULL) {
        ex = new_cept (PARAMETER_ERROR, NULL_VALUE, FATAL_SEVERITY);
        add_function (ex, "straight_subdivision");
        add_source (ex, "mscircle.c");
        add_message (ex, "two null vertices in line segment");
		return(0);
	}
	if (vtx0 == NULL || vtx1 == NULL) {
        ex = new_cept (PARAMETER_ERROR, NULL_VALUE, FATAL_SEVERITY);
        add_function (ex, "straight_subdivision");
        add_source (ex, "mscircle.c");
        add_message (ex, "one null vertex in line segment");
		return(0);
	}

	/* determine length of arc */
	len = arc_length (given_arc);
	if (error()) return(0);

	/* calculate number of subdivisions and number of points */

	n_division = n_pieces;
	if (n_division < 2) n_division = 2;
	n_point = n_division + 1;

	/* allocate memory */
	vertex_list = (struct vertex **)
		allocate_pointers (VERTEX, n_point);
	if (vertex_list == NULL) {
        ex = new_cept (MEMORY_ERROR, ALLOCATION, FATAL_SEVERITY);
        add_function (ex, "straight subdivision");
        add_source (ex, "mscircle.c");
        add_long (ex, "n_point", n_point);
        add_object (ex, VERTEX, "vertex_list");
		return(0);
	}

	/* angle subtended by each new, short arc */
	increment = len / n_division;

	/* calculate orthonormal frame */

	/* base vector */
	for (k = 0; k < 3; k++)
		base[k] = (vtx1 -> center[k] - vtx0 -> center[k]) / len;

	*vertex_list = vtx0;
	*(vertex_list + n_division) = vtx1;

	d = increment;
	for (i = 1; i < n_division; i++, d += increment) {
		for (k = 0; k < 3; k++)
			vertex_center[k] = vtx0 -> center[k] + d *  base[k];
		vtx = new_vertex (vertex_center, NULL, NULL, given_arc, (struct face *) NULL);
		if (vtx == NULL) return (0);
		link_vertex (this_srf, vtx);
		*(vertex_list+i) = vtx;
	}

	/* store vertex array pointer in calling routine */
	*vertex_list_address = vertex_list;
	return (n_point);
}


/* make a straight arc from vtx1 to vtx2 */
struct arc *make_straight_arc (struct surface *this_srf, struct vertex *vtx1, struct vertex *vtx2)
{
	int shape, small, k;
	long lfn, ofn, ffn;
	double circle_radius;
	double circle_center[3], circle_axis[3];
	struct circle *cir;
	struct arc *straight_arc;

	shape = STRAIGHT;
	small = 1;
	/* use large (last) face number */
	lfn = 0;
	if (vtx1 -> lfn > lfn) lfn = vtx1 -> lfn;
	if (vtx2 -> lfn > lfn) lfn = vtx2 -> lfn;
	/* use smaller number because it is likely more concave */
	ofn = vtx1 -> ofn;
	if (vtx2 -> ofn < ofn) ofn = vtx2 -> ofn;
	ffn = ofn;
	for (k = 0; k < 3; k++) {
		circle_center[k] = (vtx1 -> center[k] + vtx2 -> center[k]) / 2;
		circle_axis[k] = 0.0;
	}
	/* circle_radius */
	circle_radius = 0.0;
	/* new circle and arc */
	cir = new_circle (circle_center, circle_radius, circle_axis);
	if (error()) return (NULL);
	link_circle (this_srf, cir);
	straight_arc = new_arc (cir, vtx1, vtx2, shape, small, (double) 0.0, lfn, ofn, ffn);
	if (error()) return(NULL);
	link_arc (this_srf, straight_arc);
	return (straight_arc);
}

/* make a straight edge from vtx1 to vtx2 */
struct edge *make_straight_edge (struct surface *this_srf, struct vertex *vtx1, struct vertex *vtx2)
{
	struct edge *straight_edge;
	struct arc *straight_arc;

	straight_arc =  make_straight_arc (this_srf, vtx1, vtx2);
	if (error()) return(NULL);
	straight_edge = new_edge (straight_arc, 0, (struct face *) NULL, NULL);
	if (error()) return(NULL);
	return (straight_edge);
}

void store_arc (struct circle *cir, int shape, struct vertex *v1, struct vertex *v2, struct arc *a)
{
	/* set up fields */
	a -> cir = cir;
	a -> shape = (short) shape;
	a -> vtx[0] = v1;
	a -> vtx[1] = v2;
	/* compute phi angle */
	if (cir -> radius <= 0.0) a -> phi = 0.0;
	else a -> phi = arc_ang (a);
}



void store_circle (struct circle *cir, double center[3], double radius, double axis[3])
{
	int k;

	/* set up fields */
	cir -> radius = radius;
	for (k = 0; k < 3; k++) {
		cir -> center[k] = center[k];
		cir -> axis[k] = axis[k];
	}
}

struct circle *allocate_circle ()
{
	struct circle *cir;
    struct cept *ex;

	/* allocate memory */
	cir = (struct circle *) allocate_object (CIRCLE);
	if (cir == NULL) {
        ex = new_cept (MEMORY_ERROR, ALLOCATION, FATAL_SEVERITY);
        add_function (ex, "allocate_circle");
        add_source (ex, "mscircle.c");
		return(NULL);
	}
	return (cir);
}

struct circle *new_circle (double center[3], double radius, double axis[3])
{
	struct circle *cir;
    struct cept *ex;

	/* allocate memory */
	cir = allocate_circle ();
	if (error()) return (NULL);
	if (cir == NULL) {
        ex = new_cept (MEMORY_ERROR, ALLOCATION, FATAL_SEVERITY);
        add_function (ex, "new_circle");
        add_source (ex, "mscircle.c");
		add_source (ex, "mscircle.c");
		return(NULL);
	}
	store_circle (cir, center, radius, axis);
	return (cir);
}

int link_circle (struct surface *srf, struct circle *cir)
{
	/* link into list */
	if (srf -> head_circle == NULL) srf -> head_circle = cir;
	else srf -> tail_circle -> next = cir;
	srf -> tail_circle = cir;
	/* increment number of circles for molecule */
	srf -> n_circle++;
	cir -> number = srf -> n_circle;
	return (1);
}

void free_circle (struct circle *cir)
{
	free_object (CIRCLE, (short *) cir);
}

