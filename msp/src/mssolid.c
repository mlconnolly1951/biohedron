#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

/* Molecular Surface Package Copyright 1995 by Michael L. Connolly */
/* March 10, 2000 */


/* default evaluation spheres */

void set_omega_radii (struct surface *msphn, double eval_radius, double expansion_radius, double extension_radius)
{
	int k;
	long i;
	char message[MAXLINE];
	struct phnvtx *vtx;
	struct evalpnt *e;

	if (msphn == NULL) return;
	msphn -> n_evaluation = msphn -> n_phnvtx;

	msphn -> evaluation_spheres = (struct evalpnt *)
		allocate_objects (EVALPNT, msphn -> n_evaluation);
	if (msphn -> evaluation_spheres == NULL) {
		set_error1 ("(default_eval): no memory");
		return;
	}

	/* save evaluation radius for polyhedron */
	msphn -> eval_radius = eval_radius;
	/* get the coordinates */
	for (i = 0; i < msphn -> n_evaluation; i++) {
		vtx = num2phnvtx (msphn, i + 1);
		if (vtx == NULL) return;
		e = msphn -> evaluation_spheres + i;
		for (k = 0; k < 3; k++)
			e -> center[k] = vtx -> center[k] + extension_radius * vtx -> outward[k];
		e -> radius = eval_radius;
		e -> expansion = expansion_radius;
		e -> vertex_number = i + 1;
		e -> vtx = vtx;
	}
	sprintf (message,"%8ld evaluation spheres", msphn -> n_evaluation);
	inform(message);
}


void do_evaluation (struct surface *msphn1, struct surface *msphn2, struct surface *msphn3, int buried)
{
	int i, result;
	long n_numerical;
	double omega_min, omega_max;
	char message[MAXLINE];
	struct evalpnt *e;

	n_numerical = 0; omega_min = 4 * PI; omega_max = 0.0;
	for (i = 0; i < msphn1 -> n_evaluation; i++) {
		/* pointer to evaluation point */
		e = msphn1 -> evaluation_spheres + i;
		if (msphn1 != NULL && msphn2 == NULL && msphn3 == NULL) {
			result = one_evaluation (e, msphn1, 0, 0);
			n_numerical += result;
				if (e -> omegas[0] < omega_min) omega_min = e -> omegas[0];
				if (e -> omegas[0] > omega_max) omega_max = e -> omegas[0];
			if (e -> expansion > 0.0) {
				result = one_evaluation (e, msphn1, 1, 1);
				n_numerical += result;
				if (e -> omegas[1] < omega_min) omega_min = e -> omegas[1];
				if (e -> omegas[1] > omega_max) omega_max = e -> omegas[1];
				result = one_evaluation (e, msphn1, 2, 2);
				n_numerical += result;
				if (e -> omegas[2] < omega_min) omega_min = e -> omegas[2];
				if (e -> omegas[2] > omega_max) omega_max = e -> omegas[2];
			}
		}
		else if (msphn1 != NULL && msphn2 != NULL && msphn3 == NULL && !buried) {
			result = one_evaluation (e, msphn1, 0, 0);
			n_numerical += result;
			if (e -> omegas[0] < omega_min) omega_min = e -> omegas[0];
			if (e -> omegas[0] > omega_max) omega_max = e -> omegas[0];
			result = one_evaluation (e, msphn2, 1, 0);
			n_numerical += result;
			if (e -> omegas[1] < omega_min) omega_min = e -> omegas[1];
			if (e -> omegas[1] > omega_max) omega_max = e -> omegas[1];
			e -> omegas[2] = e -> omegas[0] + e -> omegas[1];
			if (e -> vtx != (struct phnvtx *) NULL)
				e -> vtx -> values[2] = e -> omegas[2];
		}
		else if (msphn1 != NULL && msphn2 != NULL && msphn3 == NULL && buried) {
			result = one_evaluation (e, msphn2, 1, 0);
			n_numerical += result;
			if (e -> omegas[1] < omega_min) omega_min = e -> omegas[1];
			if (e -> omegas[1] > omega_max) omega_max = e -> omegas[1];
		}
		else if (msphn1 != NULL && msphn2 != NULL && msphn3 != NULL) {
			result = one_evaluation (e, msphn2, 0, 0);
			n_numerical += result;
			if (e -> omegas[0] < omega_min) omega_min = e -> omegas[0];
			if (e -> omegas[0] > omega_max) omega_max = e -> omegas[0];
			result = one_evaluation (e, msphn3, 1, 0);
			n_numerical += result;
			if (e -> omegas[1] < omega_min) omega_min = e -> omegas[1];
			if (e -> omegas[1] > omega_max) omega_max = e -> omegas[1];
			e -> omegas[2] = e -> omegas[0] + e -> omegas[1];
			if (e -> vtx != (struct phnvtx *) NULL)
				e -> vtx -> values[2] = e -> omegas[2];
		}
	}
	sprintf (message, "%8ld numerical computations of omega", n_numerical);
	inform(message);
	sprintf (message, "%8.4f maximum omega", omega_max);
	inform(message);
	sprintf (message, "%8.4f minimum omega", omega_min);
	inform(message);
}

/* return 1 for analytical and 0 for numerical computation of solid angle */

int one_evaluation (struct evalpnt *e, struct surface *msphn, int which_value, int which_sphere)
{
	int k, numerical, result;
	double radius;
	char aorn[32];
	char message[256];
	struct variety *vty;
	struct solid_angle *sa;

	numerical = 0;
	/* allocate memory */
	vty = allocate_variety ();
	if (vty == NULL) {
		set_error2 ("one_evaluation: ran out of memory");
		return(0);
	}
	/* construct variety */
	for (k = 0; k < 3; k++) {
		vty -> center[k] = e -> center[k];
		vty -> axis[k] = 0.0;
	}
	radius = e -> radius + which_sphere * e -> expansion;
	vty -> radii[0] = radius;
	vty -> type = SPHERE;
	/* get sphere-polyhedron intersection face */
	sa = sphere_polyhedron (msphn, vty);
	if (error()) { set_error2("quit"); return (0);}
	if (sa == NULL) {
		set_error1("null solid_angle pointer");
		return (0);
	}
	result = compute_solid (msphn, sa);
	if (error()) return (0);
	e -> omegas[which_value] = sa -> omega;
	/* possibly store in polyhedron */
	if (e -> vtx != (struct phnvtx *) NULL)
		e -> vtx -> values[which_value] = e -> omegas[which_value];
	/* clean up triangles and edges */
	cleanup_phn (msphn);
	numerical = sa -> numerical;
	if (numerical) strcpy (aorn, "n");
	else strcpy (aorn, "a");
	sprintf (message,
		"%8ld %2ld edg %2ld cir %2ld cyc %2ld fac, sub %7.3f %2s omega %7.3f",
		e -> vertex_number, sa -> n_edge, sa -> n_circle, sa -> n_cycle, sa -> n_face,
		sa -> subtends, aorn, sa -> omega);
	informd(message);
	free_solid_angle (sa);
	return (numerical);
}

/* reset fields to zero */

void cleanup_phn (struct surface *msphn)
{
	int i;
	struct phntri *tri;
	struct phnedg *e;

	for (i = 0; i < msphn -> n_phntri; i++) {
		tri = num2phntri (msphn, i + 1);
		if (tri == NULL) {
			set_error1 ("invalid triangle number");
			return;
		}
		if (tri -> may_intersect) {
			tri -> may_intersect = 0;
			tri -> cir = NULL;
		}
	}

	for (i = 0; i < msphn -> n_phnedg; i++) {
		e = num2phnedg (msphn, i + 1);
		if (e == NULL) {
			set_error1 ("invalid edge number");
			return;
		}
		if (e -> may_intersect) {
			e -> nvtx = 0;
			e -> vtx[0] = NULL;
			e -> vtx[1] = NULL;
			e -> enter[0] = 0;
			e -> enter[1] = 0;
			e -> may_intersect = 0;
		}
	}
}

int compute_solid (struct surface *msphn, struct solid_angle *sa)
{
	double face_omega, radius;
	struct variety *vty;
	struct face *fac;
	char message[MAXLINE];
	
	radius = sa -> radius;
	if (radius < 0.0) return (0);
	vty = sa -> vty;
	if (sa -> numerical) {
		sa -> omega = numerical_omega (msphn, vty);
		if (sa -> omega < 0.0) {
			sprintf (message, "numerical_omega fails");
			set_error1 (message);
			return (0);
		}
		return (1);
	}
	else if (sa -> head_face == NULL) {
		sa -> omega = 0.0;
		return (0);
	}
	else {
		for (fac = sa -> head_face; fac != NULL; fac = fac -> next) {
			face_omega = compute_omega (0.0, fac);	/* zero ignored because convex */
			if (sa -> numerical) break;
			if (face_omega < 0.0 || face_omega > 4 * PI)
				{ sa -> numerical = 1; break;}
			sa -> omega += face_omega;
		}
	}
	if (sa -> omega < 0.0 || sa -> omega > 4 * PI)
			sa -> numerical = 1;
	if (sa -> numerical) {
		sa -> omega = numerical_omega (msphn, vty);
		if (sa -> omega < 0.0) {
			sprintf (message, "numerical_omega fails");
			set_error1 (message);
			return (0);
		}
		return (1);
	}
	return (0);
}

/* sphere contains polyhedron ? */

int polyhedron_in_sphere (struct surface *msphn, struct variety *vty)
{
	int i;
	double ds;
	struct phnvtx *vtx;

	for (i = 0; i < msphn -> n_phnvtx; i++) {
		vtx = num2phnvtx (msphn, i + 1);
		if (vtx == NULL) return (0);
		ds = distance (vty -> center, vtx -> center);
		if (ds > vty -> radii[0]) return (0);
	}
	return (1);
}

/* sphere intersect polyhedron */

struct solid_angle *sphere_polyhedron (struct surface *msphn, struct variety *vty)
{
	int i;
	struct phntri *tri;
	struct phnedg *e;
	struct solid_angle *sa;

	sa = new_solid_angle (vty);
	if (sa == NULL) return (NULL);
	/* mark triangles and edges that maybe intersect sphere */
	mark_maybe (msphn, vty);
	/* compute sphere-edge intersections */
	for (i = 0; i < msphn -> n_phnedg; i++) {
		e = num2phnedg (msphn, i + 1);
		if (e == NULL) return (NULL);
		if (!e -> may_intersect) continue;
		if (!sphere_edge (sa, vty, e)) {
			sa -> numerical = 1;
			return (sa);
		}
	}
	/* compute sphere - triangle intersections */
	for (i = 0; i < msphn -> n_phntri; i++) {
		tri = num2phntri (msphn, i + 1);
		if (tri == NULL) return (NULL);
		if (tri -> may_intersect) {
			if (!sphere_triangle (sa, tri)) {
				sa -> numerical = 1;
				return (sa);
			}
		}
	}
	/* if no edges, and sphere contains polyhedron or
		sphere center is not inside polyhedron, then intersecton is null */
	if (sa -> head_edge == NULL) {
		sa -> pis = polyhedron_in_sphere (msphn, vty);
		sa -> subtends = polyhedron_subtends (msphn, vty -> center);
		sa -> pip = (sa -> subtends > 2 * PI);
		if (sa -> pis) return (sa);
		if (!(sa -> pip)) return (sa);
	}
	/* form cycles from edges of face */
	form_cycles (NULL, sa);
	if (error()) return (NULL);
	if (sa -> numerical) return (sa);
	/* group cycles of face */
	group_cycles (NULL, NULL, sa);
	if (error()) return (NULL);
	if (sa -> numerical) return (sa);
	return (sa);
}

/* INTERSECTION PART */

/* sphere intersect edge */
int sphere_edge (struct solid_angle *sa, struct variety *vty, struct phnedg *e)
{
	int in0, in1, k, nx, j, inside;
	int xin[2];
	struct phnvtx *v0, *v1;
	double x, dot01, s0rad, nvectc1, dis0, dis1;
	double line[3], vectcv[3], vectc1[3], xpnts[2][3];
	double s0cen[3], vect0[3], vect1[3];
	char message[MAXLINE];
	struct vertex *vtx;

	if (sa == NULL) return (0);
	/* get pointers to polyhedron vertices */
	v0 = e -> pvt[0];
	v1 = e -> pvt[1];
	/* compute whether vertex is inside sphere */
	dis0 = distance (v0 -> center, vty -> center);
	dis1 = distance (v1 -> center, vty -> center);
	in0 = (dis0 < vty -> radii[0]);
	in1 = (dis1 < vty -> radii[0]);
	/* if both ends are in, whole line segment is inside,
		no intersection */
	if (in0 && in1) {
		e -> nvtx = 0;
		return (1);
	}
	/* geometry of line-sphere intersection */
	for (k = 0; k < 3; k++)
		line[k] = v1 -> center[k] - v0 -> center[k];
	if (!normalize (line)) {
		informd ("(sphere_edge): normalize fails");
		e -> nvtx = 0;
		return (0);
	}
	for (k = 0; k < 3; k++)
		vectcv[k] = v0 -> center[k] - vty -> center[k];
	x = dot_product (vectcv, line);
	for (k = 0; k < 3; k++)
		vectc1[k] = vectcv[k] - x * line[k];
	nvectc1 = norm (vectc1);
	/* check for line too far away from sphere */
	if (nvectc1 >= vty -> radii[0]) {
		e -> nvtx = 0;
		return (1);
	}
	/* compute intersection points */
	s0rad = vty -> radii[0] * vty -> radii[0] - nvectc1 * nvectc1;
	if (s0rad <= 0.0) {
		e -> nvtx = 0;
		return (1);
	}
	s0rad = sqrt (s0rad);
	for (k = 0; k < 3; k++) {
		s0cen[k] = vty -> center[k] + vectc1[k];
		xpnts[0][k] = s0cen[k] - s0rad * line[k];
		xpnts[1][k] = s0cen[k] + s0rad * line[k];
	}

	/* check which intersection points are in segment */
	nx = 0;
	e -> vtx[0] = NULL;
	e -> vtx[1] = NULL;

	for (j = 0; j < 2; j++) {
		for (k = 0; k < 3; k++) {
			vect0[k] = v0 -> center[k] - xpnts[j][k];
			vect1[k] = v1 -> center[k] - xpnts[j][k];
		}
		/* from inside segment,
			ends are in opposite directions */
		dot01 = dot_product (vect0, vect1);
		inside = (dot01 <= 0.0);
		xin[j] = inside;
		if (xin[j]) {
			vtx = new_vertex (xpnts[j], NULL, NULL, NULL, NULL);
			if (error()) return(0);
			e -> vtx[nx] = vtx;
			e -> enter[nx] = (j == 0);
			nx++;
		}
	}
	/* error checking */
	if (nx % 2 != (in0 != in1)) {
		informd ("(sphere_edge): warning: parity check fails");
		sprintf (message, "ion0 = %2d; in1 = %2d; nx = %2d", in0, in1, nx);
		informd(message);
		/* try saving no intersection points */
		e -> nvtx = 0;
		return (0);
	}
	/* store number of intersections for edge */
	e -> nvtx = nx;
	return (1);
}

/* sphere intersect plane */

int sphere_plane (struct solid_angle *sa, struct variety *vty, struct phntri *tri)
{
	int j, k, orn, somex, retval;
	double dt, rad, maxdis, ds;
	double vxp[3], cen[3];
	double cir_axis[3];
	char message[MAXLINE];
	struct phnedg *e;
	struct phnvtx *vtx;

	for (k = 0; k < 3; k++)
		cir_axis[k] = ( - tri -> axis[k]);
	for (k = 0; k < 3; k++)
		vxp[k] = tri -> center[k] - vty -> center[k];
	dt = dot_product (vxp, tri -> axis);
	for (k = 0; k < 3; k++)
		cen[k] = vty -> center[k] + dt * tri -> axis[k];
	rad = vty -> radii[0] * vty -> radii[0] - dt * dt;
	/* if circle radius imaginary, then no intersection */
	if (rad <= 0.0) return (1);	/* should be zero ? */
	rad = sqrt (rad);
	/* initialize some-intersections variable */
	somex = 0;
	for (j = 0; j < 3; j++)
		if (tri -> edg[j] -> nvtx > 0) somex = 1;
	/* generate circle if it intersects at least one edge */
	if (somex) {
		tri -> cir = new_circle (cen, rad, cir_axis);
		if (error()) return (0);
		link_sa_circle (sa, tri -> cir);
		return (1);
	}
	/* circle does not intersect triangle */
	/* containment checks */
	/* if circle center not in triangle, no intersection is possible */
	retval = point_in_triangle (cen, tri);
	if (retval < 0) {
		informd ("(sphere_plane): point_in_triangle fails");
		sprintf (message, "circle center  = %9.5f %9.5f %9.5f",
			cen[0], cen[1], cen[2]);
		informd(message);
		sprintf (message, "circle radius  = %9.5f", rad);
		informd(message);
		sprintf (message, "variety center = %9.5f %9.5f %9.5f",
			vty -> center[0], vty -> center[1], vty -> center[2]);
		informd(message);
		return (0);
	}
	if (!retval) return (1);
	/* one contains the other, determine which */
	if (rad >= tri -> radius) return (1);
	maxdis = 0.0;
	for (j = 0; j < 3; j++) {
		e = tri -> edg[j];
		orn = tri -> orn[j];
		vtx = e -> pvt[orn];
		ds = distance (vtx -> center, cen);
		if (ds > maxdis) maxdis = ds;
	}
	if (rad >= maxdis) return (1);
	/* circle inside triangle, special kind of intersection */
	tri -> cir = new_circle (cen, rad, cir_axis);
	if (error()) return (0);
	link_sa_circle (sa, tri -> cir);
	return (1);
}


/* sphere intersects triangle */

int sphere_triangle (struct solid_angle *sa, struct phntri *tri)
{
	int j, k, nx, nstart, nstop, orn, ibase, istart, istop;
	int start[6], used[6];
	int ncreate;
	double sangle, best;
	double base[3], angle[6], vect[3];
	char message[MAXLINE];
	struct phnedg *e;
	struct variety *vty;
	struct vertex *these[6];
	struct circle *cir;
	struct edge *edg;
	struct arc *a;

	/* get pointer to sphere variety */
	vty = sa -> vty;
	/* initialize circle pointer to null */
	tri -> cir = NULL;
	/* compute sphere-plane intersection */
	if (!sphere_plane (sa, vty, tri)) {
		informd ("(sphere_triangle): sphere_plane error");
		return (0);
	};
	if (tri -> cir == NULL) return (1);
	/* sphere intersects triangle: create edges */
	/* count intersections */
	nx = 0;
	nstart = 0;
	nstop = 0;
	for (j = 0; j < 3; j++) {
		e = tri -> edg[j];
		orn = tri -> orn[j];
		for (k = 0; k < e -> nvtx; k++) {
			/* mark whether vertex is start of edge */
			start[nx] = (orn != e -> enter[k]);
			if (start[nx]) nstart++;
			else nstop++;
			/* store pointer to vertex */
			these[nx] = e -> vtx[k];
			nx++;
		}
	}
	/* error checking */
	if (nstart != nstop) {
		sprintf (message, "(sphere_triangle): nstart (%d) != nstop (%d)", nstart, nstop);
		informd (message);
		return (0);
	}
	cir = tri -> cir;
	/* arc is entire circle if no intersection points */
	if (nx <= 0) {
		a = new_arc (cir, NULL, NULL, CONVEX, 0, (double) 0.0, 0L, 0L, 0L);
		if (error()) return (0);
		edg = new_edge (a, 0, NULL, NULL);
		if (error()) return(0);
		link_sa_edge (sa, edg);
		return (1);
	}
	/* arcs with edge_handles */
	/* initialize intersection points to not used in arc */
	for (j = 0; j < nx; j++)
		used[j] = 0;
	/* find base vector */
	ibase = -1;
	for (j = 0; j < nx; j++)
		if (start[j]) {
			ibase = j;
			break;
		}
	if (ibase < 0) {
		informd("(sphere_triangle): base not found");
		return (0);
	}
	for (k = 0; k < 3; k++)
		base[k] = these[ibase] -> center[k] - cir -> center[k];
	if (!normalize (base)) {
		informd("(sphere_triangle): normalize fails");
		return (0);
	}
	/* compute angles relative to base vector */
	for (j = 0; j < nx; j++) {
		if (j == ibase) angle[j] = 0.0;
		else {
			for (k = 0; k < 3; k++)
				vect[k] = these[j] -> center[k] - cir -> center[k];
			if (!normalize (vect)) {
				informd("(sphere_triangle): normalize fails");
				return (0);
			}
			angle[j] = positive_angle (base, vect, cir -> axis);
		}
	}
	/* initialize number of arcs created */
	ncreate = 0;
	while (ncreate < nx / 2) {
		/* find start vertex */
		istart = -1;
		for (j = 0; j < nx; j++) {
			if (used[j]) continue;
			if (start[j]) {
				istart = j;
				used[j] = 1;
				break;
			}
		}
		if (istart < 0) {
			informd("(sphere_triangle): start not found");
			return (0);
		}
		/* store start vertex angle */
		sangle = angle[istart];
		/* find stop vertex */
		istop = -1;
		best = 2 * PI;
		for (j = 0; j < nx; j++) {
			if (used[j]) continue;
			if (angle[j] < sangle) continue;
			if (start[j]) continue;
			if (angle[j] < best) {
				best = angle[j];
				istop = j;
			}
		}
		if (istop < 0) {
			informd("(sphere_triangle): warning stop not found");
			return (0);
		}
		used[istop] = 1;	/* mark used */
		if (best >= 2 * PI) {
			/* bad arc, don't actually create */
			informd("(sphere_triangle): warning: skip bad arc");
			return (0);
		}
		/* create new arc and new edge */
		a = new_arc (cir, these[istart], these[istop], CONVEX, 0, (double) 0.0, 0L, 0L, 0L);
		if (error()) return (0);
		edg = new_edge (a, 0, NULL, NULL);
		if (error()) return (0);
		/* link arc into list for solid angle */
		link_sa_edge (sa, edg);
		ncreate++;
	}
	return (1);
}

/* mark triangles and edges that maybe intersect omega sphere */

void mark_maybe (struct surface *msphn, struct variety *vty)
{
	int i, j;
	double d;
	struct phnedg *e;
	struct phntri *tri;

	for (i = 0; i < msphn -> n_phnedg; i++) {
		e = num2phnedg (msphn, i + 1);
		if (e == NULL) return;
		e -> may_intersect = 0;
	}

	for (i = 0; i < msphn -> n_phntri; i++) {
		tri = num2phntri (msphn, i + 1);
		if (tri == NULL) return;
		tri -> may_intersect = 0;
		d = distance (tri -> center, vty -> center);
		if (d <= vty -> radii[0] + tri -> radius) {
			tri -> may_intersect = 1;
			for (j = 0; j < 3; j++)
				tri -> edg[j] -> may_intersect = 1;
		}
	}
}


struct solid_angle *new_solid_angle (struct variety *vty)
{
	struct solid_angle *sa;

	sa = (struct solid_angle *) allocate_object (SOLID_ANGLE);
	sa -> vty = vty;
	sa -> radius = vty -> radii[0];
	return (sa);
}

void free_solid_angle (struct solid_angle *sa)
{
	free_sa_faces (sa);
	free_sa_circles (sa);
	free_variety (sa -> vty);
	free_object (SOLID_ANGLE, (short *) sa);
}

int link_sa_circle (struct solid_angle *sa, struct circle *cir)
{
	/* link into list */
	if (sa -> head_circle == NULL) sa -> head_circle = cir;
	else sa -> tail_circle -> next = cir;
	sa -> tail_circle = cir;
	/* increment number of circles for solid angle */
	sa -> n_circle++;
	cir -> number = sa -> n_circle;
	return (1);
}


int link_sa_edge (struct solid_angle *sa, struct edge *edg)
{
	/* link into list */
	if (sa -> head_edge == NULL) sa -> head_edge = edg;
	else sa -> tail_edge -> next = edg;
	sa -> tail_edge = edg;
	/* increment number of circles for solid angle */
	sa -> n_edge++;
	return (1);
}


void free_sa_circles (struct solid_angle *sa)
{
	struct circle *cir, *next_cir;

	next_cir = NULL;
	for (cir = sa -> head_circle; cir != NULL; cir = next_cir) {
		next_cir = cir -> next;
		free_circle (cir);
	}
	sa -> head_circle = NULL;
	sa -> tail_circle = NULL;
	sa -> n_circle = 0;
}


void free_sa_faces (struct solid_angle *sa)
{
	struct face *fac, *next_fac;

	next_fac = NULL;
	for (fac =  sa -> head_face; fac != NULL; fac = next_fac) {
		next_fac = fac -> next;
		deep_free_face (fac);
	}
	sa -> head_face = NULL;
	sa -> tail_face = NULL;
	sa -> n_face = 0;
} 

