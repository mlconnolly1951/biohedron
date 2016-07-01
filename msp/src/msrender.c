/*
 * MSDraw
 * Copyright 1986 by Michael L. Connolly
 * All Rights Reserved
 * January 8, 2002
 */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

/* RENDERING */

void render_surface (struct msscene *ms, struct surface *current_surface)
{
	int is_vdw, type;
	double fine_pixel;
	char message[MAXLINE];
	struct face *fac;
	struct variety *vty;
    struct cept *ex;

	fine_pixel = ms -> pixel_width / ms -> fineness;
	is_vdw = (current_surface -> type == PQMS_SURFACE && current_surface -> probe_radius <= 0.0);

	for (fac = current_surface -> head_face; fac != NULL; fac = fac -> next) {
		vty = fac -> vty;
		if (fac -> problem) {
			sprintf (message,"skip problem face of atoms %5d %5d %5d %5d",
				vty -> atmnum[0], vty -> atmnum[1], vty -> atmnum[2], vty -> atmnum[3]);
			inform(message);
			continue;
		}
		if (is_vdw && (fac -> shape == CONCAVE || fac -> shape == SADDLE))
			continue;
		type = (int) vty -> type;
		switch (type) {
		case SPHERE:
			slice_sphere (ms, current_surface, fine_pixel, fac);
			break;
		case TORUS:
			slice_torus (ms, current_surface, fine_pixel, current_surface -> probe_radius, fac);
			break;
		case CYLINDER:
			slice_cylinder (ms, current_surface, fine_pixel, fac);
			break;
		default:
			ex = new_cept (ENUM_ERROR,  INVALID_VALUE,  FATAL_SEVERITY);
			add_function (ex, "render_surface");
			add_source (ex, "msrender.c");
            add_long (ex, "variety type", (long) type);
			return;
		}
		if (error()) return;
	}
}

/* FACE SLICING */

/* slice spherical face into leaves */

void slice_sphere (struct msscene *ms, struct surface *current_surface, double fine_pixel, struct face *fac)
{
	int k, j;
	double yinc, y, rad, fsgn;
	double cir_center[3], cir_axis[3];
	char message[MAXLINE];
	struct leaf *lf;
	struct circle *cir;
	struct variety *vty;
    struct cept *ex;

	vty = fac -> vty;
	if (vty -> type != SPHERE) {
		ex = new_cept (ENUM_ERROR,  INVALID_VALUE,  FATAL_SEVERITY);
		add_function (ex, "slice_sphere");
		add_source (ex, "msrender.c");
        add_long (ex, "variety type", (long) vty -> type);
		return;
	}


	/* check versus window */
	for (k = 0; k < 3; k++) {
		if (vty -> center[k] + vty -> radii[0] < ms -> window[0][k]) return;
		if (vty -> center[k] - vty -> radii[0] > ms -> window[1][k]) return;
	}

	/* plus one for convex; minus one for concave */
	fsgn = ((fac -> shape == CONVEX) ? 1.0 : -1.0);


	/* set up circle for leaf */
	for (k = 0; k < 3; k++) {
		cir_center[k] = vty -> center[k];
		cir_axis[k] = ((k == 1) ? fsgn : 0.0);
	}
	cir = new_circle (cir_center, (double) 0.0, cir_axis);
	if (cir == NULL) {
		add_object (tail_cept, CIRCLE, "leaf circle");
		add_function (tail_cept, "slice_sphere");
		return;
	}
	lf = allocate_leaf ();
	if (lf == NULL) {
		add_object (tail_cept, LEAF, "leaf");
		add_function (tail_cept, "slice_sphere");
		return;
	}

	/* copy atom number from variety to leaf */
	for (j = 0; j < MAXPA; j++)
		lf -> atmnum[j] = fac -> vty -> atmnum[j];
	for (k = 0; k < 3; k++)
		lf -> focus[k] = vty -> center[k];

	/* set up leaf fields */
	lf -> cir = cir;
	lf -> shape = fac -> shape;
	lf -> type = fac -> vty -> type;
	lf -> fac = fac;
	lf -> side = OUTSIDE;
	lf -> comp = fac -> comp;
	lf -> input_hue = fac -> input_hue;

	/* y increment for lines of latitude on sphere */
	yinc = fine_pixel;

	/* one leaf per line of latitude */

	for (y = vty -> center[1] - vty -> radii[0] - yinc / 2;
		y < vty -> center[1] + vty -> radii[0]; y += yinc) {
		/* change circle center */
		cir -> center[1] = y;
		/* radius of circle of latitude */
		rad = (vty -> radii[0] * vty -> radii[0]) - (y - vty -> center[1])
			 * (y - vty -> center[1]);
		if (rad <= 0.0) continue;
		rad = sqrt (rad);
		if (rad <= 0.0) continue;
		cir -> radius = rad;
		/* leaf endpoints: west and east */
		for (j = 0; j < 2; j++)
			for (k = 0; k < 3; k++)
				lf -> ends[j][k] = cir -> center[k];
		lf -> ends[0][0] -= rad;
		lf -> ends[1][0] += rad;
		/* determine accessibility of endpoints of leaf */
		for (j = 0; j < 2; j++) {
			lf -> where[j] = point_in_face (lf -> ends[j], fac, 1);
			if (lf -> where[j] < 0) {
				ms -> n_bad_projection++;
				lf -> where[j] = 0;
			}
		}

		/* cut, clip and render (outer) leaf */
		lf -> cep = 0;
		lf -> clip_ep = 0;
		cut_leaf (ms, current_surface, fine_pixel, lf);
		if (error()) return;

		if (current_surface -> clipping) {
			/* cut, clip and render (inner) leaf */
			for (k = 0; k < 3; k++)
				cir -> axis[k] = ((k == 1) ? -fsgn : 0.0);
			lf -> cep = 0;
			lf -> clip_ep = 0;
			lf -> side = INSIDE;
			cut_leaf (ms, current_surface, fine_pixel, lf);
			if (error()) return;
			/* reset what we changed */
			lf -> side = OUTSIDE;
			for (k = 0; k < 3; k++)
				cir -> axis[k] = ((k == 1) ? fsgn : 0.0);
		}
	}
	free_leaf (lf);
	free_circle (cir);
	return;
}

/* slice toroidal face into leaves */

void slice_torus (struct msscene *ms, struct surface *current_surface, double fine_pixel, double probe_radius, struct face *fac)
{
	int k, j, i, nfocus, near1, naif;
	double anginc, bigrad;
	double focus[3], vect1[3], vect2[3], vect[3], qvect[3];
	double dtq, tcv[3];
	double *foci = (double *) NULL;
	char message[MAXLINE];
	struct leaf *lf;
	struct circle *cir1, *cir2, *cir3;
	struct circle *lfcir, *torcir;
	struct variety *vty, *atm1, *atm2;
	struct arc *a, *nxta;
	struct arc *torarc;
	struct vertex *torvtx[2];
	struct vertex *qvtx;
	struct vertex *conevtx;
	struct cycle *cyc;
	struct edge *edg;
    struct cept *ex;

	vty = fac -> vty;
	if (vty -> type != TORUS) {
		ex = new_cept (GEOMETRY_ERROR,  INCONSISTENCY,  FATAL_SEVERITY);
		add_function (ex, "slice_torus");
		add_source (ex, "msrender.c");
        add_long (ex, "variety type", (long) vty -> type);
		return;
	}
	if (vty -> tube) {
		slice_elbow (ms, current_surface, fine_pixel, fac);
		return;
	}

	if (debug >= 2) {
		sprintf (message,"render saddle face for atoms %5d %5d",
			vty -> atmnum[0], vty -> atmnum[1]);
		inform(message);
	}

	/* get pointers to atom varieties */
	atm1 = *(current_surface -> variety_handles + fac -> vty -> atmnum[0] - 1);
	atm2 = *(current_surface -> variety_handles + fac -> vty -> atmnum[1] - 1);

	/* check versus window */
	bigrad = distance (atm1 -> center, atm2 -> center) +
		atm1 -> radii[0] + atm2 -> radii[0];

	for (k = 0; k < 3; k++) {
		if (vty -> center[k] + bigrad < ms -> window[0][k]) return;
		if (vty -> center[k] - bigrad > ms -> window[1][k]) return;
	}
	/* leaf circle */
	lfcir = allocate_circle ();
	if (error()) {
		add_object (tail_cept, CIRCLE, "leaf circle");
		add_function (tail_cept, "slice_torus");
		return;
	}
	/* leaf */
	lf = allocate_leaf ();
	if (error()) {
		add_object (tail_cept, LEAF, "leaf");
		add_function (tail_cept, "slice_sphere");
		return;
	}
	/* torus circle radius, center, axis */
	torcir = new_circle (vty -> center, vty -> radii[0], vty -> axis);
	if (torcir == NULL) {
		add_object (tail_cept, CIRCLE, "torus circle");
		add_function (tail_cept, "slice_circle");
		return;
	}
	/* torus arc */
	torarc = allocate_arc ();
	if (error()) {
		add_object (tail_cept, ARC, "torus arc");
		add_function (tail_cept, "slice_torus");
		add_source (tail_cept, "msrender.c");
		return;
	}
	for (j = 0; j < 2; j++) {
		torvtx[j] = allocate_vertex ();
		if (error()) {
			add_object (tail_cept, VERTEX, "torus vertex");
			add_function (tail_cept, "slice_torus");
			add_source (tail_cept, "msrender.c");
			return;
		}
	}
	torarc -> cir = torcir;
	/* copy atom numbers from variety to leaf */
	for (k = 0; k < MAXPA; k++)
		lf -> atmnum[k] = fac -> vty -> atmnum[k];

	/* set up leaf fields */
	lf -> cir = lfcir;
	lf -> shape = fac -> shape;
	lf -> type = fac -> vty -> type;
	lf -> fac = fac;
	lf -> cep = 0;
	lf -> clip_ep = 0;
	lf -> side = OUTSIDE;
	lf -> comp = fac -> comp;
	lf -> input_hue = fac -> input_hue;

	/* both endpoints of saddle face leaf are always accessible */
	for (j = 0; j < 2; j++)
		lf -> where[j] = ACCESSIBLE;

	/* angular increment for rotation of leaf about torus axis */
	anginc = fine_pixel / (vty -> radii[0]);

	/* next we need endpoints for torus arc */
	/* get them from concave arcs bounding saddle face */

	/* intialization */
	cir1 = NULL;
	cir2 = NULL;
	cir3 = NULL;
	qvtx = NULL;
	conevtx = NULL;
	near1 = 0;

	/* look for concave arcs */
	naif = 0;
	for (cyc = fac -> first_cycle; cyc != NULL; cyc = cyc -> next)
		for (edg = cyc -> first_edge; edg != NULL; edg = edg -> next) {
			naif++;
			a = edg -> arcptr;
			if (a -> shape == CONVEX) {
				cir3 = a -> cir;
				continue;
			}
			if (edg -> next == NULL)
				nxta = cyc -> first_edge -> arcptr;
			else
				nxta = edg -> next -> arcptr;
			if (along (edg, vty -> axis))
				cir2 = a -> cir;
			else
				cir1 = a -> cir;
			/* check for cusp vertex */
			if (a -> shape == CONCAVE && nxta -> shape == CONCAVE) {
				/* cusp point joints two concave arcs */
				qvtx = a -> vtx[1-edg->orn];
			}
		}

	dtq = probe_radius * probe_radius - vty -> radii[0] * vty -> radii[0];

	/* later: note: check PI in bubbles */

	if (naif == 1) {
		if (dtq <= 0.0) {
			ex = new_cept (GEOMETRY_ERROR,  INCONSISTENCY,  FATAL_SEVERITY);
			add_function (ex, "slice_torus");
			add_source (ex, "msrender.c");
			add_message(ex, "toroidal face with only one arc, but not cone");
			return;
		}
		if (cir3 == NULL) {
			ex = new_cept (GEOMETRY_ERROR,  INCONSISTENCY,  FATAL_SEVERITY);
			add_function (ex, "slice_torus");
			add_source (ex, "msrender.c");
			add_message(ex, "toroidal face with only one arc, but no contact circle");
			return;
		}
		/* cone */
		qvtx = allocate_vertex ();
		if (error()) {
			add_object (tail_cept, VERTEX, "CUSP VERTEX");
			add_function (tail_cept, "slice_torus");
			add_source (tail_cept, "msrender.c");
			return;
		}
		conevtx = qvtx;
		dtq = sqrt (dtq);
		for (k = 0; k < 3; k++)
			tcv[k] = cir3 -> center[k] - torcir -> center[k];
		normalize (tcv);
		for (k = 0; k < 3; k++)
			qvtx -> center[k] = torcir -> center[k] + dtq * tcv[k];
		/* hope this is enough */
	}
	if (cir1 == NULL) informd2 ("cir1 null");
	if (cir2 == NULL) informd2 ("cir2 null");
	if (qvtx != NULL) informd2 ("cusp present");

	/* check for cusp vertex */
	if (qvtx != NULL) {
		for (k = 0; k < 3; k++)
			qvect[k] = qvtx -> center[k] - vty -> center[k];
		near1 = (dot_product (qvect, vty -> axis) < 0.0);
	}

	/* check for hoop saddle face */
	if (cir1 == NULL || cir2 == NULL) {
		for (j = 0; j < 2; j++)
			torarc -> vtx[j] = NULL;
		informd2 ("complete toroidal hoop");
	}
	else {
		/* concave arc circle centers are endpoints of sphere rolling */
		for (k = 0; k < 3; k++) {
			torvtx[0] -> center[k] = cir1 -> center[k];
			torvtx[1] -> center[k] = cir2 -> center[k];
		}
		for (j = 0; j < 2; j++)
			torarc -> vtx[j] = torvtx[j];
		sprintf (message, "saddle rendering (from): %8.3f %8.3f %8.3f",
			cir1 -> center[0], cir1 -> center[1], cir1 -> center[2]);
		informd2 (message);
		sprintf (message, "saddle rendering (to)  : %8.3f %8.3f %8.3f",
			cir2 -> center[0], cir2 -> center[1], cir2 -> center[2]);
		informd2 (message);
	}

	/* the probe sphere centers are the foci of the leaves */
	nfocus = render_sub_arc (torarc, &foci, anginc);
	if (nfocus < 2) {
		ex = new_cept (LOGIC_ERROR, MSUNDERFLOW, FATAL_SEVERITY);
        add_function (ex, "slice_torus");
        add_source (ex, "msrender.c");
        add_long (ex, "number of foci", (long) nfocus);
		return;
	}
	sprintf (message, "nfocus = %d", nfocus);
	informd2 (message);

	/* create leaves */
	for (i = 0; i < nfocus; i++) {
		for (k = 0; k < 3; k++) {
			focus[k] = (*(foci + 3 * i + k));
			lfcir -> center[k] = focus[k];
			lf -> focus[k] = focus[k];
		}

		/* unit vectors from focus toward atoms */
		for (k = 0; k < 3; k++) {
			vect1[k] = atm1 -> center[k] - focus[k];
			vect2[k] = atm2 -> center[k] - focus[k];
		}
		/* correct for cusp vertex */
		if (qvtx != NULL) {
			if (near1)
				for (k = 0; k < 3; k++)
					vect2[k] = qvtx -> center[k] - focus[k];
			else
				for (k = 0; k < 3; k++)
					vect1[k] = qvtx -> center[k] - focus[k];
		}
		/* normalize vectors to unit length */
		normalize (vect1);
		normalize (vect2);

		/* leaf circle radius is probe radius */
		lfcir -> radius = probe_radius;
		/* set up endpoints of leaf */
		for (k = 0; k < 3; k++) {
			lf -> ends[0][k] = focus[k] + lfcir -> radius * vect1[k];
			lf -> ends[1][k] = focus[k] + lfcir -> radius * vect2[k];
		}
		/* compute leaf circle axis */
		for (k = 0; k < 3; k++)
			vect[k] = focus[k] - vty -> center[k];
		cross (vty -> axis, vect, lfcir -> axis);
		normalize (lfcir -> axis);

		/* clip and render leaf */
		clip_leaf (ms, current_surface, fine_pixel, lf);
		if (error()) return;
	}

	/* return temporary memory */
	if (!free_doubles (foci, 0, VERTS)) {
		ex = new_cept (MEMORY_ERROR,  FREEING,  FATAL_SEVERITY);
		add_variable (ex, VERTS, "foci");
		add_function (ex, "slice_torus");
		add_source (ex, "msrender.c");
		return;
	}

	free_leaf (lf);
	free_arc (torarc);
	free_circle (torcir);
	free_circle (lfcir);
	for (j = 0; j < 2; j++)
		free_vertex (torvtx[j]);
	if (conevtx != NULL) free_vertex (conevtx);
	return;
}

/* slice elbow (tube) tori */

void slice_elbow (struct msscene *ms, struct surface *current_surface, double fine_pixel, struct face *fac)
{
	char message[MAXLINE];
	struct leaf *lf;
	struct circle *lfcir, *torcir;
	struct arc *torarc;
	struct vertex *torvtx[2];
	int i, j, k, nfocus, atmnum;
	double anginc, bigrad;
	double atmcen[3], ccens[2][3], cradii[2];
	double focus[3], vect[3], z_axis[3], base[3];
	struct variety *vty;
	double *foci;
    struct cept *ex;
	
	vty = fac -> vty;
	atmnum = vty -> atmnum[0];
	if (debug >= 2) {
		sprintf (message,"render elbow face for atom %5d", atmnum);
		inform(message);
	}
	for (k = 0; k < 3; k++)
		atmcen[k] = *(current_surface -> atom_centers + 3 * (atmnum - 1) + k);
	for (j = 0; j < 2; j++)
		cradii[j] = vty -> radii[1];
	for (j = 0; j < 2; j++)
		for (k = 0; k < 3; k++) {
			ccens[j][k] = vty -> ccens[j][k];
		}
	bigrad = distance (atmcen, ccens[0]) + cradii[0];
	for (k = 0; k < 3; k++) {
		if (atmcen[k] + bigrad < ms -> window[0][k]) return;
		if (atmcen[k] - bigrad > ms -> window[1][k]) return;
	}
	/* leaf circle */
	lfcir = allocate_circle ();
	if (lfcir == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
		add_object (ex, CIRCLE, "leaf circle");
		add_function (ex, "slice_elbow");
		add_source (ex, "msrender.c");
		return;
	}
	/* leaf */
	lf = allocate_leaf ();
	if (lf == NULL) {
		add_object (tail_cept, LEAF, "leaf");
		add_function (tail_cept, "slice_sphere");
		return;
	}
	/* torus circle radius, center, axis */
	torcir = new_circle (vty -> center, vty -> radii[0], vty -> axis);
	if (torcir == NULL) {
		add_object (tail_cept, CIRCLE, "torus circle");
		add_function (tail_cept, "slice_elbow");
		return;
	}
	/* torus arc */
	torarc = allocate_arc ();
	if (torarc == NULL) {
		add_object (tail_cept, ARC, "torus arc");
		add_function (tail_cept, "slice_elbow");
		add_source (tail_cept, "msrender.c");
		return;
	}
	for (j = 0; j < 2; j++) {
		torvtx[j] = allocate_vertex ();
		if (error()) {
			add_object (tail_cept, VERTEX, "torus vertex");
			add_function (tail_cept, "slice_elbow");
			add_source (tail_cept, "msrender.c");
			return;
		}
	}
	/* set up leaf fields */
	for (k = 0; k < MAXPA; k++)
		lf -> atmnum[k] = vty -> atmnum[k];
	lf -> cir = lfcir;
	lf -> shape = CONVEX;		/* to avoid reversing normal vector */
	lf -> type = vty -> type;
	lf -> fac = fac;
	lf -> cep = 0;
	lf -> clip_ep = 0;
	lf -> side = OUTSIDE;
	lf -> comp = fac -> comp;
	lf -> input_hue = fac -> input_hue;
	for (j = 0; j < 2; j++)
		lf -> where[j] = ACCESSIBLE;

	/* setup torus central circle for subdivision */
	anginc = fine_pixel / ((vty -> radii[0] + cradii[0]));
	torcir -> radius = vty -> radii[0];
	for (k = 0; k < 3; k++) {
		torcir -> center[k] = vty -> center[k];
		torcir -> axis[k] = vty -> axis[k];
	}
	torarc -> cir = torcir;
	for (j = 0; j < 2; j++) {
		torarc -> vtx[j] = torvtx[j];
		for (k = 0; k < 3; k++)
			torvtx[j] -> center[k] = ccens[j][k];
	}
	foci = (double *) NULL;
	nfocus = render_sub_arc (torarc, &foci, anginc);
	if (nfocus < 2) {
		ex = new_cept (LOGIC_ERROR, MSUNDERFLOW, FATAL_SEVERITY);
        add_function (ex, "slice_elbow");
        add_source (ex, "msrender.c");
        add_long (ex, "number of foci", (long) nfocus);
		return;
	}

	/* create leaves */
	for (i = 0; i < nfocus; i++) {
		for (k = 0; k < 3; k++) {
			focus[k] = (*(foci + 3 * i + k));
			lfcir -> center[k] = focus[k];
			lf -> focus[k] = focus[k];
		}
		lfcir -> radius = cradii[0];
		/* compute tangent to torus central circle */
		for (k = 0; k < 3; k++)
			vect[k] = focus[k] - vty -> center[k];
		cross (vty -> axis, vect, lfcir -> axis);
		normalize (lfcir -> axis);
		for (k = 0; k < 3; k++)
			z_axis[k] = ((k == 2) ? 1.0 : 0.0);
		cross (lfcir -> axis, z_axis, base);
		if (norm (base) <= 0.0) {
			continue;
		}
		normalize (base);
		for (k = 0; k < 3; k++) {
			lf -> ends[0][k] = lfcir -> center[k] - lfcir -> radius * base[k];
			lf -> ends[1][k] = lfcir -> center[k] + lfcir -> radius * base[k];
		}
		/* clip and render leaf */
		clip_leaf (ms, current_surface, fine_pixel, lf);
		if (error()) return;
	}
	free_doubles (foci, 0, VERTS);
	free_leaf (lf);
	free_arc (torarc);
	free_circle (torcir);
	free_circle (lfcir);
	for (j = 0; j < 2; j++)
		free_vertex (torvtx[j]);
}
	
/* slice cylindrical face into leaves */

void slice_cylinder (struct msscene *ms, struct surface *current_surface, double fine_pixel, struct face *fac)
{
	int k, j;
	double axis_inc, a;
	double half_length, big_radius, cyl_radius;
	double cyl_axis[3], z_axis[3], base[3];
	struct leaf *lf;
	struct circle *lfcir;
	struct variety *vty;
    struct cept *ex;

	vty = fac -> vty;
	cyl_radius = vty -> radii[0];
	half_length = vty -> length / 2;
	big_radius = ((half_length > cyl_radius) ? half_length : cyl_radius);


	/* check versus window */
	for (k = 0; k < 3; k++) {
		if (vty -> center[k] + big_radius < ms -> window[0][k]) return;
		if (vty -> center[k] - big_radius > ms -> window[1][k]) return;
	}
	/* leaf circle */
	lfcir = allocate_circle ();
	if (lfcir == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
		add_object (ex, CIRCLE, "leaf circle");
		add_function (ex, "slice_cylinder");
		add_source (ex, "msrender.c");
		return;
	}
	/* leaf */
	lf = allocate_leaf ();
	if (lf == NULL) {
		add_object (tail_cept, LEAF, "leaf");
		add_function (tail_cept, "slice_sphere");
		return;
	}

	/* set up circle for leaf */
	lfcir -> radius = cyl_radius;
	for (k = 0; k < 3; k++) {
		lfcir -> center[k] = vty -> center[k];
		lfcir -> axis[k] = vty -> axis[k];
		cyl_axis[k] = vty -> axis[k];
		lf -> focus[k] = vty -> center[k];
	}
	/* check for end-on */
	if (cyl_axis[0] == 0.0 && cyl_axis[1] == 0.0) {
		free_leaf (lf);
		free_circle (lfcir);
		return;
	}
	for (k = 0; k < 3; k++)
		z_axis[k] = ((k == 2) ? 1.0 : 0.0);
	cross (cyl_axis, z_axis, base);
	if (norm (base) <= 0.0) {
		free_leaf (lf);
		free_circle (lfcir);
		return;
	}
	normalize (base);

	/* copy 3 atom number from variety to leaf */
	for (k = 0; k < MAXPA; k++)
		lf -> atmnum[k] = vty -> atmnum[k];

	/* set up leaf fields */
	lf -> cir = lfcir;
	lf -> shape = fac -> shape;
	lf -> type = fac -> vty -> type;
	lf -> fac = fac;
	lf -> side = OUTSIDE;
	lf -> comp = fac -> comp;
	lf -> input_hue = fac -> input_hue;

	/* increment for lines of latitude on cylinder */
	axis_inc = fine_pixel;

	/* one leaf per line of latitude */

	for (a = (-half_length-axis_inc/2); a < half_length; a += axis_inc) {
		/* change circle center */
		for (k = 0; k < 3; k++) {
			lfcir -> center[k] = vty -> center[k] + a * cyl_axis[k];
			lf -> focus[k] = lfcir -> center[k];
		}
		/* leaf endpoints: west and east */
		for (k = 0; k < 3; k++) {
			lf -> ends[0][k] = lfcir -> center[k] - cyl_radius * base[k];
			lf -> ends[1][k] = lfcir -> center[k] + cyl_radius * base[k];
		}

		/* determine accessibility of endpoints of leaf */
		for (j = 0; j < 2; j++)
			lf -> where[j] = 1;

		/* cut, clip and render (outer) leaf */
		lf -> cep = 0;
		lf -> clip_ep = 0;
		clip_leaf (ms, current_surface, fine_pixel, lf);
		if (error()) return;

		if (current_surface -> clipping) {
			/* cut, clip and render (inner) leaf */
			for (k = 0; k < 3; k++)
				lfcir -> axis[k] = (-cyl_axis[k]);
			lf -> cep = 0;
			lf -> clip_ep = 0;
			lf -> side = INSIDE;
			clip_leaf (ms, current_surface, fine_pixel, lf);
			/* reset what we changed */
			lf -> side = OUTSIDE;
			for (k = 0; k < 3; k++)
				lfcir -> axis[k] = cyl_axis[k];
		}
	}
	free_leaf (lf);
	free_circle (lfcir);
	return;
}

void render_polyhedron (struct msscene *ms, struct surface *phn, double fine_pixel, long interpolate)
{
	int j, k, check, comp, hue, f, atm;
	int orn0, orn1, orn2;
	long t;
	long n_back, n_clipped, n_rendered;
	int vclipped[3], vback;
	double tcen[3];
	double trad;
	double tvertices[3][3];
	double tnormals[4][3];
	double tvalues[3];
	char message[MAXLINE];
	struct phnedg *edg0, *edg1, *edg2;
	struct phntri *tri;
	struct phnvtx *vtx, *vtxs[3];
    struct cept *ex;

	n_back = 0;
	n_clipped = 0;
	n_rendered = 0;
	f = (phn -> function[0]);
	if (phn -> phntri_handles == NULL) return;
	for (t = 0; t < phn -> n_phntri; t++) {
		tri = num2phntri (phn, t+1);
		if (tri == NULL) {
			ex = new_cept (POINTER_ERROR,  NULL_POINTER,  FATAL_SEVERITY);
			add_object (ex, PHNTRI, "tri");
			add_function (ex, "render_polyhedron");
			add_source (ex, "msrender.c");
			return;
		}
		for (k = 0; k < 3; k++)
			tcen[k] = tri -> center[k];
		trad = tri -> radius;
		edg0 = tri -> edg[0];
		edg1 = tri -> edg[1];
		edg2 = tri -> edg[2];
		orn0 = tri -> orn[0];
		orn1 = tri -> orn[1];
		orn2 = tri -> orn[2];
		vtxs[0] = edg0 -> pvt[orn0];
		vtxs[1] = edg1 -> pvt[orn1];
		vtxs[2] = edg2 -> pvt[orn2];
		for (k = 0; k < 3; k++)
			tnormals[3][k] = tri -> axis[k];
		for (j = 0; j < 3; j++) {
			vtx = vtxs[j];
			if      (f == 'x') tvalues[j] = vtx -> center[0];
			else if (f == 'y') tvalues[j] = vtx -> center[1];
			else if (f == 'z') tvalues[j] = vtx -> center[2];
			else if (f == 'u') tvalues[j] = vtx -> values[0];
			else if (f == 'v') tvalues[j] = vtx -> values[1];
			else if (f == 'w') tvalues[j] = vtx -> values[2];
			else tvalues[j] = vtx -> center[2]; /* default */
			for (k = 0; k < 3; k++) {
				tvertices[j][k] = vtx -> center[k];
				tnormals[j][k] = vtx -> outward[k];
			}
			/* check back-facing for triangle vertex normal */
			/* check clipping for triangle vertex */
			vclipped[j] =  ms -> clipping && phn -> clipping && clipped (ms, tvertices[j]);
		}
		vback =  (tnormals[3][2] < 0.0);
		if (ms -> clipping && phn -> clipping) {
			/* if all three vertices are clipped, triangle is clipped */
			check = (vclipped[0] || vclipped[1] || vclipped[2]);
			if (vclipped[0] && vclipped[1] && vclipped[2]) {
				n_clipped++;
				continue;
			}
		}
		else {
			check = 0;	/* no need to check triangle for being partially clipped */
			/* if back-facing, skip triangle */
			if (vback) {
				n_back++;
				continue;
			}
		}
		comp = tri -> comp;
		atm = tri -> atm;
		if (atm <= 0) {
			ex = new_cept (LOGIC_ERROR,  BOUNDS,  FATAL_SEVERITY);
			add_function (ex, "render_polyhedron");
			add_source (ex, "msrender.c");
			add_long (ex, "atm", (long) atm);
			add_message(ex, "invalid atom number for triangle");
			return;
		}
		hue = tri -> hue;
		render_triangle (ms, phn, fine_pixel, interpolate, 
			tcen, trad, tvertices, tnormals, tvalues, check, comp, hue, atm);
		if (error()) return;
		n_rendered++;
	}
	if (ms -> clipping && phn -> clipping)
		sprintf (message,"%8ld triangles rendered, %6ld clipped",
			n_rendered, n_clipped);
	else
		sprintf (message,"%8ld triangles rendered, %6ld back-facing",
			n_rendered, n_back);
	inform(message);
}

void render_triangle (struct msscene *ms, struct surface *phn, double fine_pixel, long interpolate,
double tcen[3], double trad, double tvertices[3][3], double tnormals[4][3], double tvalues[3],
int check, int comp, int trihue, int atm)
{
	int k, point_clipped, x, y, z;
	int hue, shade, inner, shape;
	double input_opacity;
	double xf, yf, zf;
	double xmin, xmax, ymin, ymax, xp, yp;
	double xinc, yinc, zp, value;
	double norvect[3], pnt[3], horvect[3];
	double a0, a1, a2, b0, b1, b2, c0, c1, c2;
	double b1c1, b0c0, a0c0, a1c1, xpc0, ypc1;
	double denom, numa, numb;
	double fa, fb, fc;
	double opacity;
	struct object_scheme *scheme;
	
	scheme = phn -> scheme;
	inner = 0;

	xmin = tcen[0] - trad; xmax = tcen[0] + trad;
	ymin = tcen[1] - trad; ymax = tcen[1] + trad;

	/* extra kludge factor of 0.7 */
	xinc = 0.7 * fine_pixel;
	yinc = 0.7 * fine_pixel;
	if (xinc <= 0.0) xinc = EPSILON;
	if (yinc <= 0.0) yinc = EPSILON;
	a0 = tvertices[0][0]; a1 = tvertices[0][1]; a2 = tvertices[0][2];
	b0 = tvertices[1][0]; b1 = tvertices[1][1]; b2 = tvertices[1][2];
	c0 = tvertices[2][0]; c1 = tvertices[2][1]; c2 = tvertices[2][2];
	denom = (a0 - c0) * (b1 - c1) - (a1 - c1) * (b0 - c0);
	if (denom == 0.0) return;
	b1c1 = b1 - c1;
	b0c0 = b0 - c0;
	a0c0 = a0 - c0;
	a1c1 = a1 - c1;
	for (xp = xmin; xp <= xmax; xp += xinc) {
		xpc0 = xp - c0;
		for (yp = ymin; yp <= ymax; yp += yinc) {
			ypc1 = yp - c1;
			numa = b1c1 * xpc0 - b0c0 * ypc1;
			fa = numa / denom;
			if (fa <= 0.0) continue; if (fa >= 1.0) continue;
			numb = a0c0 * ypc1 - a1c1 * xpc0;
			fb = numb / denom;
			if (fb <= 0.0) continue; if (fb >= 1.0) continue;
			fc = 1.0 - fa - fb;
			if (fc <= 0.0) continue; if (fc >= 1.0) continue;
			zp = fa * a2 + fb * b2 + fc * c2;
			pnt[0] = xp;
			pnt[1] = yp;
			pnt[2] = zp;
			if (check) {
				point_clipped = ms -> clipping && phn -> clipping && clipped (ms, pnt);
				if (point_clipped) continue;
			}
			if (interpolate) {
				for (k = 0; k < 3; k++)
					norvect[k] = fa * tnormals[0][k] +
								 fb * tnormals[1][k] +
								 fc * tnormals[2][k];
			}
			else {
				for (k = 0; k < 3; k++)
					norvect[k] = tnormals[3][k];
			}
			normalize(norvect);
			if (norvect[2] < 0.0) {
				if (!phn -> clipping) continue;
				if (!ms -> clipping) continue;
				/* reverse direction of vector */
				inner = 1;
				for (k = 0; k < 3; k++)
					norvect[k] *= (-1.0);
				for (k = 0; k < 3; k++) 
					horvect[k] = norvect[k];
				horvect[2] = 0.0;
				/* move inner side in one & half pixel widths */
				xf = xp + phn -> surface_thickness * ms -> pixel_width * horvect[0];
				yf = yp + phn -> surface_thickness * ms -> pixel_width * horvect[1];
				zf = zp + phn -> surface_thickness * ms -> pixel_width * horvect[2];
			}
			else {
				inner = 0;
				xf = xp;
				yf = yp;
				zf = zp;
			}
			/* convert from surface to screen coordinates */
			x = ftoi (ms, xf, 0);
			y = ftoi (ms, yf, 1);
			z = ftoi (ms, zf, 2);
			
			value = fa * tvalues[0] + fb * tvalues[1] + fc * tvalues[2];
			shape = 1;
			
			/* hue from color scheme */
			hue = determine_hue (ms -> table, scheme, atm, comp, shape, inner, trihue, value);
			if (error()) return;

			/* get shade (intensity) */
			shade = detsh(ms, norvect, z);
			if (error()) return;
			input_opacity = 1 - trihue % 2;
			opacity = detopac (atm, inner, comp, shape, phn -> scheme, input_opacity);
			
			/* put pixel into depth buffer */
			putpix (ms, phn, inner, shade, hue, opacity, x, y, z);
			if (error()) return;
		}
	}
}


/* CLIPPING and CUTTING */

void cut_leaf (struct msscene *ms, struct surface *current_surface, double fine_pixel, struct leaf *lf)
{
	int j, k, n, nx, nw, orn, sgn, i, used0, used1, nused;
	double t, xsmin, sangle, fsgn;
	double base[3], xpnts[2][3], vect[3];
	char message[MAXLINE];
	struct circle *cir;
	struct arc *al, *a;
	struct variety *vty;
	struct face *fac;
	struct leaf *lfcut;
	struct lax *xsptr, *xsend, *xsp;
	struct cycle *cyc;
	struct edge *edg;
	struct lax laxes[MAX_LAX];
    struct cept *ex;

	fac = lf -> fac;

	/* count number of edges for face */
	n = 2 * edges_in_face (fac);

	/* if toroidal, no cutting */

	if (fac -> shape == SADDLE || n <= 0) {
		/* clip and render leaf */
		clip_leaf (ms, current_surface, fine_pixel, lf);	
		return;
	}

	/* determine sign (plus for convex, minus for concave) */
	fsgn = ((lf -> shape == CONVEX) ? 1.0 : -1.0);

	/* create arc for leaf */
	al = leaf_to_arc (lf);
	vty = fac -> vty;
	cir = lf -> cir;
	for (k = 0; k < 3; k++)
		base[k] = lf -> ends[0][k] - cir -> center[k];
	normalize (base);

	/* determine intersection points between leaf and arcs of face */
	xsptr = &(laxes[0]);
	for (cyc = fac -> first_cycle; cyc != NULL; cyc = cyc -> next)
		for (edg = cyc -> first_edge; edg != NULL; edg = edg -> next) {
			a = edg -> arcptr;
			orn = edg -> orn;
			sgn = 1 - 2 * orn;
			/* call arc-arc intersection function */
			nx = arc_arc (al, a, xpnts);
			if (nx <= 0) continue;	/* no intersections */
			/* fill data into list */
			for (i = 0; i < nx; i++) {
				if (xsptr - &(laxes[0]) >= n) {
					ex = new_cept (ARRAY_ERROR,  MSOVERFLOW,  FATAL_SEVERITY);
					add_function (ex, "cut_leaf");
					add_source (ex, "msrender.c");
                    add_long (ex, "maximum number of leaf-arc intersections", n);
					return;
				}
				for (k = 0; k < 3; k++) {
					xsptr -> co[k] = xpnts[i][k];
					vect[k] = fsgn * (xsptr -> co[k] - vty -> center[k]);
				}
				t = sgn * triple_product (al -> cir -> axis, a -> cir -> axis, vect);
				xsptr -> ent = (t < 0);
				xsptr -> used = 0;
				xsptr++;
			}
		}

	xsend = xsptr;		/* mark end of list of intersection points */
	nx = xsend - &(laxes[0]);	/* compute number of intersection points */

	/* intersection parity error checking */
	if ((lf -> where[0] == lf -> where[1]) == (nx % 2)) {
		/* parity check fails,
		   flag leaf to check accessibility of every pixel of leaf */
		ms -> n_bdy_parity++;
		lf -> cep = 1;
		clip_leaf (ms, current_surface, fine_pixel, lf);
		lf -> cep = 0;
		frearc (al);
		return;
	}
	/* kludge to work around undiscovered bug for quartet cusps */
	if (lf -> atmnum[3] != 0) {
		/* flag leaf to check accessibility of every pixel of leaf */
		lf -> cep = 1;
		clip_leaf (ms, current_surface, fine_pixel, lf);
		lf -> cep = 0;
		frearc (al);
		return;
	}


	/* check for no intersection points */
	if (nx == 0) {
		if (lf -> where[0] == INACCESSIBLE) {
			/* entire leaf inaccessible: nothing to render */
			/* free temporary memory */
			frearc (al);
			return;
		}
		/* entire leaf accessible: render entire leaf */
		clip_leaf (ms, current_surface, fine_pixel, lf);
		if (error()) return;
		/* free temporary memory */
		frearc (al);
		return;
	}

	/* calculate angles for transitions between accessible
		and inaccessible */
	for (xsptr = &(laxes[0]); xsptr < xsend; xsptr++) {
		for (k = 0; k < 3; k++)
			vect[k] = xsptr -> co[k] - cir -> center[k];
		normalize (vect);
		xsptr -> angle = positive_angle (base, vect, cir -> axis);
	}

	/* initialization */
	used0 = 0;
	used1 = 0;
	nused = 0;
	nw = 0;

	/* count number of accessible endpoints of originial leaf */
	for (j = 0; j < 2; j++)
		if (lf -> where[j]) nw++;

	/* create new leaves */

	while (nused < nx + nw) {
		/* get starting point */
		if (!used0 && lf -> where[0] == ACCESSIBLE) {
			/* start at original endpoint */
			used0 = 1;		/* first originial endpoint used */
			nused++;		/* increment number used */
			lfcut = duplicate_leaf (lf);	/* duplicate leaf */
			sangle = 0.0;			/* starting angle */
		}
		else {
			/* look for unused cut vertex */
			xsmin = 4 * PI;
			xsp = NULL;
			for (xsptr = &(laxes[0]); xsptr < xsend; xsptr++) {
				if (xsptr -> used) continue;
				if (xsptr -> angle < xsmin) {
					xsmin = xsptr -> angle;
					xsp = xsptr;
				}
			}

			if (xsp == NULL) {
				if (lf -> where[1] == INACCESSIBLE) break;
				ms -> n_missing_leaves++;
				return;
			}
			if (!xsp -> ent) {
				ms -> n_missing_leaves++;
				return;
			}
			xsp -> used = 1;		/* mark as used */
			nused++;			/* increment number used */
			sangle = xsp -> angle;		/* starting angle */
			lfcut = duplicate_leaf (lf);	/* duplicate leaf */
			for (k = 0; k < 3; k++)
				lfcut -> ends[0][k] = xsp -> co[k];
		}

		/* get ending point */

		/* initialization */
		xsmin = 4 * PI;
		xsp = NULL;

		for (xsptr = &(laxes[0]); xsptr < xsend; xsptr++) {
			if (xsptr -> used) continue;		/* already used */
			if (xsptr -> angle < sangle) continue;	/* end after start */
			if (xsptr -> angle < xsmin) {
				/* best so far, save */
				xsmin = xsptr -> angle;
				xsp = xsptr;
			}
		}

		if (xsp == NULL) {
			if (used1 || lf -> where[1] == INACCESSIBLE) {
				ms -> n_missing_leaves++;
				return;
			}
			/* use East */
			used1 = 1;		/* mark East original endpoint used */
			nused++;		/* increment number used */
			for (k = 0; k < 3; k++)
				lfcut -> ends[1][k] = lf -> ends[1][k];
		}
		else {	/* use cut point */
			for (k = 0; k < 3; k++)
				lfcut -> ends[1][k] = xsp -> co[k];
			xsp -> used = 1;		/* mark as used */
			nused++;				/* increment number used */
		}
		/* we have a cut leaf; clip and render leaf */
		clip_leaf (ms, current_surface, fine_pixel, lfcut);
		free_leaf (lfcut);				/* free temporary memory */
	}

	/* free temporary memory */
	frearc (al);
	return;
}

int clipped (struct msscene *ms, double pnt[3])
{
	int k, result;
	double clip_center[3];
	double clip_axis[3];

	for (k = 0; k < 3; k++) {
		clip_center[k] = ms -> clip_center[k];
		clip_axis[k] = ms -> clip_axis[k];
	}
	result = pclipped (clip_center, clip_axis, pnt);
	return (result);
}

/* clip and render leaf */

void clip_leaf (struct msscene *ms, struct surface *current_surface, double fine_pixel, struct leaf *lf)
{
	int j, k, n, nx, nw, i, used0, used1, nused;
	double t, xsmin, sangle;
	double base[3], xpnts[2][3], vect[3];
	double clip_center[3], clip_axis[3];
	struct circle *cir;
	struct arc *al;
	struct leaf *lfcut;
	struct lax *xsptr, *xsend, *xsp;
	struct lax laxes[MAX_LAX];

	/* count double number of clipping planes */
	n = 0;
	if (ms -> clipping) n += 2;

	if (n <= 0 || !(current_surface -> clipping)) {
		lf -> clip_ep = 0;
		render_leaf (ms, current_surface, fine_pixel, lf);
		return;
	}

	/* determine accessibility of endpoints of leaf */
	for (j = 0; j < 2; j++)
		lf -> when[j] = !(current_surface -> clipping && clipped (ms, lf -> ends[j]));

	/* create arc for leaf */
	al = leaf_to_arc (lf);
	cir = lf -> cir;
	for (k = 0; k < 3; k++)
		base[k] = lf -> ends[0][k] - cir -> center[k];
	normalize (base);

	/* determine intersection points between leaf and clipping planes */
	xsptr = &(laxes[0]);
	if (ms -> clipping) {
		for (k = 0; k < 3; k++) {
			clip_center[k] = ms -> clip_center[k];
			clip_axis[k] = ms -> clip_axis[k];
		}
		/* call arc-plane intersection function */
		nx = arc_plane (al, clip_center, clip_axis, xpnts);
		if (nx < 0) return;
		/* if (nx == 0) continue;	no intersections */
		/* fill data into list */
		/* figure out enter or exit */
		for (i = 0; i < nx; i++) {
			if (xsptr - &(laxes[0]) >= n) {
				ms -> n_missing_leaves++;
				return;
			}
			for (k = 0; k < 3; k++)
				xsptr -> co[k] = xpnts[i][k];
			for (k = 0; k < 3; k++)
				vect[k] = xsptr -> co[k] - al -> cir -> center[k];
			t = triple_product (al -> cir -> axis, clip_axis, vect);
			xsptr -> ent = (t > 0); /* or reverse? */
			xsptr -> used = 0;
			xsptr++;
		}
	}

	xsend = xsptr;		/* mark end of list of intersection points */
	nx = xsend - &(laxes[0]);	/* compute number of intersection points */

	/* error checking */
	if ((lf -> when[0] == lf -> when[1]) == (nx % 2)) {
		/* parity check fails,
		   flag leaf to check clipping of every pixel of leaf */
		ms -> n_clip_parity++;
		lf -> clip_ep = 1;
		render_leaf (ms, current_surface, fine_pixel, lf);
		if (error()) return;
		lf -> clip_ep = 0;
		frearc (al);
		return;
	}

	/* check for no intersection points */
	if (nx == 0) {
		if (lf -> when[0] == CLIPPED) {
			/* entire leaf clipped: nothing to render */
			/* free temporary memory */
			frearc (al);
			return;
		}
		/* entire leaf unclipped: render entire leaf */
		render_leaf (ms, current_surface, fine_pixel, lf);
		if (error()) return;
		/* free temporary memory */
		frearc (al);
		return;
	}

	/* calculate angles for transitions between unclipped and clipped */
	for (xsptr = &(laxes[0]); xsptr < xsend; xsptr++) {
		for (k = 0; k < 3; k++)
			vect[k] = xsptr -> co[k] - cir -> center[k];
		normalize (vect);
		xsptr -> angle = positive_angle (base, vect, cir -> axis);
	}

	/* initialization */
	used0 = 0;
	used1 = 0;
	nused = 0;
	nw = 0;

	/* count number of unclipped endpoints of originial leaf */
	for (j = 0; j < 2; j++)
		if (lf -> when[j]) nw++;

	/* create new leaves */

	while (nused < nx + nw) {
		/* get starting point */
		if (!used0 && lf -> when[0] == UNCLIPPED) {
			/* start at original endpoint */
			used0 = 1;		/* first originial endpoint used */
			nused++;		/* increment number used */
			lfcut = duplicate_leaf (lf);	/* duplicate leaf */
			sangle = 0.0;			/* starting angle */
		}
		else {
			/* look for unused cut vertex */
			xsmin = 4 * PI;
			xsp = NULL;
			for (xsptr = &(laxes[0]); xsptr < xsend; xsptr++) {
				if (xsptr -> used) continue;
				if (xsptr -> angle < xsmin) {
					xsmin = xsptr -> angle;
					xsp = xsptr;
				}
			}

			if (xsp == NULL) {
				if (lf -> when[1] == CLIPPED) break;
				ms -> n_missing_leaves++;
				return;
			}
			if (!xsp -> ent) {
				ms -> n_missing_leaves++;
				return;
			}
			xsp -> used = 1;		/* mark as used */
			nused++;		/* increment number used */
			sangle = xsp -> angle;		/* starting angle */
			lfcut = duplicate_leaf (lf);		/* duplicate leaf */
			for (k = 0; k < 3; k++)
				lfcut -> ends[0][k] = xsp -> co[k];
		}

		/* get ending point */

		/* initialization */
		xsmin = 4 * PI;
		xsp = NULL;

		for (xsptr = &(laxes[0]); xsptr < xsend; xsptr++) {
			if (xsptr -> used) continue;		/* already used */
			if (xsptr -> angle < sangle) continue;	/* end after start */
			if (xsptr -> angle < xsmin) {
				/* best so far, save */
				xsmin = xsptr -> angle;
				xsp = xsptr;
			}
		}

		if (xsp == NULL) {
			if (used1 || lf -> when[1] == CLIPPED) {
				ms -> n_missing_leaves++;
				return;
			}
			/* use East */
			used1 = 1;	/* mark East original endpoint used */
			nused++;	/* increment number used */
			for (k = 0; k < 3; k++)
				lfcut -> ends[1][k] = lf -> ends[1][k];
			}
			else {	/* use cut point */
				for (k = 0; k < 3; k++)
					lfcut -> ends[1][k] = xsp -> co[k];
				xsp -> used = 1;	/* mark as used */
				nused++;		/* increment number used */
			}
			/* we have a cut leaf */
			/* clip and render leaf */
			render_leaf (ms, current_surface, fine_pixel, lfcut);
			if (error()) return;
			free_leaf (lfcut);	/* free temporary memory */
		}

	/* free temporary memory */
	frearc (al);

	return;
}

void render_leaf (struct msscene *ms, struct surface *current_surface, double fine_pixel, struct leaf *lf)
{
	int i, k, npoint, x, y, z, shade, hue, atm;
	int inner, point_clipped, point_if;
	int shape, comp, input_hue;
	double input_opacity;
	double opacity;
	double norvect[3], horvect[3], pnt[3];
	static double *points = (double *) NULL;
	struct circle *cir;

	cir = lf -> cir;
	if (cir -> radius <= 0.0) {
		ms -> n_missing_leaves++;
		return;
	}


	/* check leaf circle versus window */
	for (k = 0; k < 3; k++) {
		if (cir -> center[k] + cir -> radius < ms -> window[0][k]) return;
		if (cir -> center[k] - cir -> radius > ms -> window[1][k]) return;
	}

	/* subdivide leaf into points (pixels) */

	npoint = subdivide_leaf (fine_pixel, lf, &points);
	if (npoint < 2) {
		ms -> n_missing_leaves++;
		return;
	}

	/* for each point, compute normal vector, shade, hue */
	for (i = 0; i < npoint; i++) {
		for (k = 0; k < 3; k++)
			pnt[k] = *(points+3*i+k);
		for (k = 0; k < 3; k++) norvect[k] = (*(points+3*i+k)
			- lf -> focus[k]);
		normalize (norvect);
		if (lf -> shape != CONVEX) reverse (norvect);
		/* check every point ? */
		if (lf -> cep) {
			/* check accessibility for problem leaf */
			point_if = point_in_face (points + 3 * i, lf -> fac, 1);
			if (error()) return;
			if (!point_if) continue;
		}
		if (lf -> clip_ep) {
			/* check clipping for problem leaf */
			point_clipped =  current_surface -> clipping && clipped (ms, points + 3 * i);
			if (point_clipped) continue;
		}
		inner = 0;
		/* if no clipping, do not render backward facing pixels */
		if (norvect[2] <= 0.0) {
			if (!(current_surface -> clipping)) continue;
			/* reverse direction of vector */
			inner = 1;
			for (k = 0; k < 3; k++)
				norvect[k] *= (-1.0);
			for (k = 0; k < 3; k++) 
				horvect[k] = norvect[k];
			horvect[2] = 0.0;
			/* move inner side in one & half pixel widths */
			for (k = 0; k < 3; k++)
				*(points+3*i+k) +=
					 current_surface -> surface_thickness * ms -> pixel_width * horvect[k];
		}

		/* convert from surface to screen coordinates */
		x = ftoi (ms, *(points + 3 * i), 0);
		y = ftoi (ms, *(points + 3 * i + 1), 1);
		z = ftoi (ms, *(points + 3 * i + 2), 2);
		
		/* closest atom of (1, 2, 3) atoms to point */
		atm = closest (current_surface, lf, points + 3 * i);
		if (error()) return;
		if (atm <= 0) continue;

		shape = lf -> shape;
		if (lf -> type == CYLINDER) shape = 2;
		if (shape < 1) shape = 1;
		if (shape > 3) shape = 3;
		comp = lf -> comp;
		input_hue = lf -> input_hue;
		input_opacity = 1 - input_hue % 2;
	
		opacity = detopac (atm, inner, comp, shape, current_surface -> scheme, input_opacity);
		if (error()) return;

		/* hue from color scheme */
		hue = determine_hue (ms -> table, current_surface -> scheme, atm, comp, shape, inner, input_hue, 0.0);
		if (error()) return;

		/* get shade (intensity) */
		shade = detsh(ms, norvect, z);
		if (error()) return;

		/* put pixel into depth buffer */
		putpix (ms, current_surface, inner, shade, hue, opacity, x, y, z);
		if (error()) return;
	}
	/* free points of leaf */
	free_doubles (points, 0, VERTS);
}

/* subdivide leaf into points */

int subdivide_leaf (double fine_pixel, struct leaf *lf, double **pntptr)
{
	int j, k, n;
	double alpha;
	struct arc *lfarc;
	struct vertex *lfvtx[2];
    struct cept *ex;

	if (lf -> cir -> radius <= 0.0) {
		return (0);
	}

	lfarc = allocate_arc ();
	if (lfarc == NULL) {
		add_object (tail_cept, ARC, "leaf arc");
		add_function (tail_cept, "subdivide_leaf");
		add_source (tail_cept, "msrender.c");
		return (0);
	}
	for (j = 0; j < 2; j++) {
		lfvtx[j] = allocate_vertex ();
		if (error()) {
			add_object (tail_cept, VERTEX, "leaf vertex");
			add_function (tail_cept, "subdivide_leaf");
			add_source (tail_cept, "msrender.c");
			return (0);
		}
	}
	/* convert to arc for subdivision */
	lfarc -> cir = lf -> cir;
	for (k = 0; k < 3; k++) {
		lfvtx[0] -> center[k] = lf -> ends[0][k];
		lfvtx[1] -> center[k] = lf -> ends[1][k];
	}

	for (j = 0; j < 2; j++)
		lfarc -> vtx[j] = lfvtx[j];

	/* angular parameter for subdivision */
	alpha = fine_pixel / (lf -> cir -> radius);

	/* call geometric subdivision function */
	n = render_sub_arc (lfarc, pntptr, alpha);
	if (error()) return(0);

	free_arc (lfarc);
	for (j = 0; j < 2; j++)
		free_vertex (lfvtx[j]);
	
	return (n); 	/* return number of subdivision points */

}

int render_sub_arc (struct arc *a, double **vaddr, double alpha)
{
	int ndiv, npoint;
	double *verts;
	char message[MAXLINE];

	npoint = count_sub_arc (a, alpha);
	if (error()) return (0);
	ndiv = npoint - 1;

	verts = geo_sub_arc (a, npoint);
	if (error()) return (0);

	/* store vertex array pointer in calling routine */
	*vaddr = verts;
	return (npoint);
}


/* return an arc corresponding to leaf */

struct arc *leaf_to_arc (struct leaf *lf)
{
	int j, k;
	struct arc *a;
	struct vertex *v;
    struct cept *ex;

	a = allocate_arc ();
	if (a == NULL) {
		add_object (tail_cept, ARC, "leaf arc");
		add_function (tail_cept, "leaf_to_arc");
		add_source (tail_cept, "msrender.c");
		return(NULL);
	}

	a -> cir = lf -> cir;

	for (j = 0; j < 2; j++) {
		v = allocate_vertex ();
		if (error()) {
			add_object (tail_cept, VERTEX, "leaf vertex");
			add_function (tail_cept, "leaf_to_arc");
			add_source (tail_cept, "msrender.c");
			return(NULL);
		}
		for (k = 0; k < 3; k++)
			v -> center[k] = lf -> ends[j][k];
		a -> vtx[j] = v;
	}

	return (a);
}

/* duplicate leaf */

struct leaf *duplicate_leaf (struct leaf *lf)
{
	int j, k;
	struct leaf *lf2;
    struct cept *ex;

	lf2 = allocate_leaf ();
	if (lf2 == NULL) {
		add_object (tail_cept, LEAF, "leaf");
		add_function (tail_cept, "duplicate_leaf");
		return(NULL);
	}

	lf2 -> cir = lf -> cir;
	lf2 -> shape = lf -> shape;
	lf2 -> type = lf -> type;
	lf2 -> fac = lf -> fac;
	lf2 -> side = lf -> side;
	lf2 -> comp = lf -> comp;
	lf2 -> input_hue = lf -> input_hue;
	lf2 -> cep = lf -> cep;
	lf2 -> clip_ep = lf -> clip_ep;

	for (j = 0; j < 2; j++) {
		lf2 -> where[j] = lf -> where[j];
		lf2 -> when[j] = lf -> when[j];
		for (k = 0; k < 3; k++)
			lf2 -> ends[j][k] = lf -> ends[j][k];
	}

	for (k = 0; k < 3; k++)
		lf2 -> focus[k] = lf -> focus[k];

	for (k = 0; k < MAXPA; k++)
		lf2 -> atmnum[k] = lf -> atmnum[k];

	return (lf2);
}



/* determine closes of (1, 2, 3) atoms to point */

int closest (struct surface *current_surface, struct leaf *lf, double pnt[3])
{
	int j, k, jmin, nch;
	double dmin;
	double dta[MAXPA];
	double several_centers[MAXPA][3], several_radii[MAXPA];
    struct cept *ex;

	if (lf -> type == SPHERE && lf -> shape == CONVEX) nch = 1;
	else if (lf -> shape == SADDLE) nch = 2;
	else if (lf -> type == SPHERE && lf -> shape == CONCAVE) nch = 3;
	else if (lf -> type == CYLINDER) nch = 2;
	else if (lf -> type == TORUS && lf -> shape == CONVEX) nch = 1;
	else {
		ex = new_cept (LOGIC_ERROR,  NOT_FOUND,  FATAL_SEVERITY);
		add_function (ex, "closest");
		add_source (ex, "msrender.c");
        add_long (ex, "leaf type", (long) lf -> type);
        add_long (ex, "leaf shape", (long) lf -> shape);
		return(0);
	}

	if (nch <= 1)
		return (lf -> atmnum[0]);

	if (lf -> atmnum[3] > 0)
		nch = 4;
	for (j = 0; j < nch; j++) {
		for (k = 0; k < 3; k++)
			several_centers[j][k] = *(current_surface -> atom_centers + 3 * (lf -> atmnum[j] - 1) + k);
		several_radii[j] = *(current_surface -> atom_radii + (lf -> atmnum[j] - 1));
	}

	/* find closest atom */

	for (j = 0; j < nch; j++)
		dta[j] = distance (pnt, several_centers[j]) - several_radii[j];

	/* initialization */
	dmin = 1000000.0;
	jmin = -1;

	for (j = 0; j < nch; j++)
		if (dta[j] < dmin) {
			dmin = dta[j];
			jmin = j;
		}

	if (jmin < 0 || jmin >= nch) {
		ex = new_cept (LOGIC_ERROR,  NOT_FOUND,  FATAL_SEVERITY);
		add_function (ex, "closest");
		add_source (ex, "msrender.c");
		add_message (ex, "cannot find closest atom");
		return (0);
	}
	return (lf -> atmnum[jmin]);
}

/* determine shade from shading model */

int detsh (struct msscene *ms, double norvect[3], int z)
{
	int sh;
	double cang1, cang2, dt1, f, proj, perp;

	if (norvect[2] <= 0.0) return (1); /* 0 is excluded because it means background */

	/* cosine of angle for diffuse component of shading model */
	cang1 = dot_product (norvect, ms -> lsource);
	if (cang1 <= 0.0) cang1 = 0.0;

	/* cosine of angle for specular component of shading model */
	dt1 = dot_product (ms -> lsource, norvect);
	proj = dt1 * norvect[2];
	perp = ms -> lsource[2] - proj;
	cang2 = proj - perp;
	if (cang2 <= 0.0) cang2 = 0.0;

	/* f lies between 0.0 and 1.0 */
	f = ms -> model[0] * ((double) z / (double) (ms -> zview))
	  + ms -> model[1] + ms -> model[2] * cang1
	  + ms -> model[3] * pow (cang2, ms -> model[4]);

	/* convert to shade number */
	sh = f * 256;

	/* make sure in valid range for shades */
	if (sh < 1) sh = 1;	/* 0 is excluded because it means background */
	if (sh > 255) sh = 255;

	return (sh);
}

/* does edge point along axis ? */

int along (struct edge *edg, double axis[3])
{
	int orn, k, sgn;
	double dt;
	double vect[3];
	struct vertex *vtx1, *vtx2;
	struct arc *a;
    struct cept *ex;

	a = edg -> arcptr;
	vtx1 = a -> vtx[0];
	vtx2 = a -> vtx[1];
	if (vtx1 == NULL || vtx2 == NULL) {
		ex = new_cept (PARAMETER_ERROR,  NULL_VALUE,  FATAL_SEVERITY);
		add_object (ex, VERTEX, "vtx1 or vtx2");
		add_function (ex, "along");
		add_source (ex, "msrender.c");
		return(0);
	}
	orn = edg -> orn;
	sgn = 1 - 2 * orn;

	for (k = 0; k < 3; k++)
		vect[k] = vtx2 -> center[k] - vtx1 -> center[k];
	dt = dot_product (vect, axis) * sgn;
	return (dt > 0.0);
}

struct leaf *allocate_leaf ()
{
	struct leaf *lef;
    struct cept *ex;

	lef = (struct leaf *) allocate_object (LEAF);
	if (lef == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
		add_function (ex, "allocate_leaf");
		add_source (ex, "msrender.c");
		return(NULL);
	}
	return (lef);
}

void free_leaf (struct leaf *lef)
{
	free_object (LEAF, (short *) lef);
	if (error()) return;
}


/*
	MSDraw
	Copyright 1986 by Michael L. Connolly
	All Rights Reserved
*/
