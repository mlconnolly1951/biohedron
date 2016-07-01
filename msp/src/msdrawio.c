/* Molecular Surface Package
 * Copyright 1986 by Michael L. Connolly
 * All Rights Reserved
 * December 26, 2001
 */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"
 

/* INPUT */

/* read piecewise quartic molecular surface */

struct surface *surf_input (FILE *fpsurf)
{
	char molecule_name[64];
	double surface_center[3];
	struct surface *srf;

	srf = read_surface (fpsurf, molecule_name, surface_center);
	if (srf == NULL) return (NULL);
	srf -> clipping = 1;
	srf -> type = PQMS_SURFACE;
	srf -> scheme = define_scheme ((int) UNIFORM_COLORING, 2, -1, 0.0, 0.0);
	return (srf);
}

struct surface *read_surface (FILE *fp_pqms, char molecule_name[64], double surface_center[3])
{
	long i, fcynum, fednum, iatm;
	int k;
	short *n_edge_per_cycle;
	long *first_cycle_number;
	long *first_edge_number;
	struct variety *vty;
	struct face *fac;
	struct cycle *cyc;
	struct surface *srf;

	srf = read_header (fp_pqms, molecule_name);
	if (srf == NULL) {
		return (NULL);
	}
	srf -> variety_handles = read_varieties (fp_pqms, srf, surface_center);
	if (error()) return (NULL);
	srf -> atom_centers = allocate_doubles (srf -> n_atom * 3, 0, ATOM_CENTERS);
	if (srf -> atom_centers == NULL) {
		set_error1 ("(read_surface) memory allocation failure");
		return(NULL);
	}
	srf -> atom_radii = allocate_doubles (srf -> n_atom, 0, ATOM_RADII);
	if (srf -> atom_radii == NULL) {
		set_error1 ("(read_surface) memory allocation failure");
		return(NULL);
	}
	iatm = 0; /* assumes that atoms are first varieties */
	for (i = 0; i < srf -> n_atom; i++) {
		vty = *(srf -> variety_handles + i);
		if (vty -> type == SPHERE) {
			if (iatm >= srf -> n_atom) {
				set_error1 ("(read_surface) memory allocation failure");
				return(NULL);
			}
			/* store atomic coordinates */
			for (k = 0; k < 3; k++)
				*(srf -> atom_centers + 3 * iatm + k) = vty -> center[k];
			*(srf -> atom_radii + iatm) = vty -> radii[0];
			iatm++;
		}
	}
	
	srf -> vertex_handles = read_vertices (fp_pqms, srf);
	if (error()) return (NULL);
	srf -> circle_handles = read_circles (fp_pqms, srf);
	if (error()) return (NULL);
	srf -> arc_handles = read_arcs (fp_pqms, srf, srf -> circle_handles, srf -> vertex_handles);
	if (error()) return (NULL);
	if (srf -> n_face > 0) {
		first_cycle_number = allocate_longs (srf -> n_face, 0, FIRST_CYCLE_NUMBER);
		if (first_cycle_number == NULL) {
			set_error1 ("read_surface: memory allocation failure");
			return (NULL);
		}
	}
	else {
		first_cycle_number = NULL;
	}
	srf -> face_handles = read_faces (fp_pqms, srf, first_cycle_number, srf -> variety_handles);
	if (error()) return (NULL);
	if (srf -> n_cycle > 0) {
		n_edge_per_cycle = allocate_shorts (srf -> n_cycle);
		if (n_edge_per_cycle == NULL) {
			set_error1 ("read_surface: memory allocation failure");
			return (NULL);
		}
		first_edge_number = allocate_longs (srf -> n_cycle, 0, FIRST_EDGE_NUMBER);
		if (first_edge_number == NULL) {
			set_error1 ("read_surface: memory allocation failure");
			return (NULL);
		}
	}
	else {
		n_edge_per_cycle = NULL;
		first_edge_number = NULL;
	}
	srf -> cycle_handles = read_cycles (fp_pqms, srf, n_edge_per_cycle, first_edge_number);
	if (error()) return (NULL);
	for (i = 0; i < srf -> n_face; i++) {
		fac = *(srf -> face_handles + i);
		fcynum = *(first_cycle_number + i);
		fac -> first_cycle = (fcynum == 0) ? NULL : *(srf -> cycle_handles + fcynum - 1);
	}
	srf -> edge_handles = read_edges (fp_pqms, srf, n_edge_per_cycle, srf -> arc_handles);
	if (error()) return (NULL);
	for (i = 0; i < srf -> n_cycle; i++) {
		cyc = *(srf -> cycle_handles + i);
		fednum = *(first_edge_number + i);
		cyc -> first_edge = *(srf -> edge_handles + fednum - 1);
	}
	srf -> component_handles = read_components (fp_pqms, srf);
	if (error()) return (NULL);

	/* free temporary memory */
	free_shorts (n_edge_per_cycle);
	free_longs (first_edge_number, 0, FIRST_EDGE_NUMBER);
	free_longs (first_cycle_number, 0, FIRST_CYCLE_NUMBER);
	return(srf);
}

struct surface *read_header (FILE *fp_pqms, char molecule_name[64])
{
	long i;
	char message[MAXLINE];
	struct surface *srf;
	struct header headerin;

	srf = new_surface ();
	if (srf == NULL) {
		return (NULL);
	}
	/* initialize surface structure */
	srf -> weight = DEFAULT_WEIGHT;
	srf -> large_omega = DEFAULT_LARGE;
	/* read header record */
	fread ((char *) &headerin, sizeof (struct header), 1, fp_pqms);
	if (strcmp (headerin.filetype, "PQMS") != 0) {
		set_error1 ("read_header:  bad format for surface file");
		return (NULL);
	}
	if (headerin.version > VERSION ||
		(headerin.version == VERSION && headerin.subversion > SUBVERSION) ||
		headerin.version < PREVERSION ||
		(headerin.version == PREVERSION && headerin.subversion < PRESUBVERSION)) {
		sprintf (message, "MSP  version = %d.%d",
			(int) VERSION, (int) SUBVERSION);
		set_error1(message);
		sprintf (message, "file version = %d.%d",
			(int) headerin.version, (int) headerin.subversion);
		set_error2(message);
		return (NULL);
	}
	if (headerin.version < 3) shorty = 1;
	else if (headerin.version == 3 && headerin.subversion < 5) shorty = 1;
	else shorty = 0;
	for (i = 0; i < 9; i++)
		if (headerin.counters[i] < 0) {
			set_error1 ("read_surface: negative counter");
			return (NULL);
		}
	srf -> probe_radius = headerin.probe_radius;
	/* transfer to global counters */
	srf -> n_variety = headerin.counters[0];
	sprintf (message, "%8ld varieties to read", srf -> n_variety);
	informd2 (message);
	srf -> n_vertex = headerin.counters[1];
	sprintf (message, "%8ld vertices to read", srf -> n_vertex);
	informd2 (message);
	srf -> n_circle = headerin.counters[2];
	sprintf (message, "%8ld circles to read", srf -> n_circle);
	informd2 (message);
	srf -> n_arc = headerin.counters[3];
	sprintf (message, "%8ld arcs to read", srf -> n_arc);
	informd2 (message);
	srf -> n_face = headerin.counters[4];
	sprintf (message, "%8ld faces to read", srf -> n_face);
	informd2 (message);
	srf -> n_cycle = headerin.counters[5];
	sprintf (message, "%8ld cycles to read", srf -> n_cycle);
	informd2 (message);
	srf -> n_edge = headerin.counters[6];
	sprintf (message, "%8ld edges to read", srf -> n_edge);
	informd2 (message);
	srf -> n_component = headerin.counters[7];
	sprintf (message, "%8ld components to read", srf -> n_component);
	informd2 (message);
	srf -> n_atom = headerin.counters[8];
	sprintf (message, "%8ld atoms to read", srf -> n_atom);
	informd2 (message);
	if (srf -> probe_radius == 0.0) {
		srf ->  van_der_Waals = 1;
	}
	if (srf -> n_edge != 2 * srf -> n_arc) {
		set_error1 ("read_surface: bad surface file");
		set_error2 ("read_surface: #edges != double #arcs");
		return (NULL);
	}
	sprintf (message,"%6ld.%1ld surface file version",
		headerin.version, headerin.subversion);
	inform(message);
	sprintf (message,"         %s", headerin.name);
	inform(message);
	sprintf (message,"%8.3f probe radius", srf -> probe_radius);
	inform(message);
	sprintf (message,"%8ld atoms", srf -> n_atom);
	inform(message);
	strcpy (molecule_name, headerin.name);
	srf -> surface_type = PQMS_SURFACE;
	return (srf);
}

struct variety **read_varieties (FILE *fp_pqms, struct surface *srf, double surface_center[3])
{
	int k, nunpacked;
	long i;
	char message[MAXLINE];
	struct variety **variety_handles;
	struct vtybin vtyin;
	struct vtyshort vtyi;
	struct variety *vty;

	/* allocate memory for handles */
	variety_handles = (struct variety **)
		allocate_pointers (VARIETY, srf -> n_variety);
	if (variety_handles == (struct variety **) NULL) {
		set_error1 ("read_surface: memory allocation failure");
		set_error2 ("for variety handles");
		return (NULL);
	}
	for (i = 0; i < srf -> n_variety; i++) {
		vty = allocate_variety ();
		if (vty == (struct variety *) NULL) {
			set_error1 ("read_surface: variety allocation failure");
			return (NULL);
		}
		*(variety_handles + i) = vty;
	}
	/* read varieties */
	/* center of molecule */
	for (k = 0; k < 3; k++)
		surface_center[k] = 0.0;
	for (i = 0; i < srf -> n_variety; i++) {
		if (shorty) {
			fread ((char *) &vtyi, sizeof (vtyi), 1, fp_pqms);
			vtyin.type = vtyi.type;
			vtyin.radius = vtyi.radius;
			for (k = 0; k < 3; k++) {
				vtyin.center[k] = vtyi.center[k];
				vtyin.axis[k] = vtyi.axis[k];
			}
			for (k = 0; k < 3; k++)
				vtyin.atmnum[k] = vtyi.atmnum[k];
		}
		else fread ((char *) &vtyin, sizeof (vtyin), 1, fp_pqms);
		if (feof (fp_pqms)) {
			set_error1 ("read_varieties: premature eof on surface file");
			return (NULL);
		}
		vty = *(variety_handles + i);
		if (vtyin.type == ATOM_TYPE || vtyin.type == PROBE_TYPE)
			vty -> type = SPHERE;
		else if (vtyin.type == TORUS_TYPE)
			vty -> type = TORUS;
		else {
			set_error1 ("invalid variety type");
			return (NULL);
		}
		for (k = 0; k < 3; k++) {
			vty -> center[k] = vtyin.center[k];
			vty -> axis[k] = vtyin.axis[k];
			/* vty -> atmnum[k] = vtyin.atmnum[k]; */
		}
		nunpacked = unpack_atom_numbers (&vtyin, vty);
		if (nunpacked > 3) {
			sprintf (message, "%8ld variety has %d atoms",
				i + 1, nunpacked);
			informd (message);
		}
		vty -> radii[0] = vtyin.radius;
		if (vtyin.type == TORUS_TYPE)
			vty -> radii[1] = srf -> probe_radius;
		if (vtyin.type == ATOM_TYPE)
			for (k = 0; k < 3; k++)
				surface_center[k] += vty -> center[k];
	}
	if (srf -> n_atom <= 0) {
		sprintf (message,
			"read_varieties: invalid number of atoms: %ld",
				srf -> n_atom);
		set_error1(message);
		return (NULL);
	}
	for (k = 0; k < 3; k++)
		surface_center[k] /= srf -> n_atom;
	sprintf (message,"%8.3f %8.3f %8.3f surface center",
		surface_center[0], surface_center[1], surface_center[2]);
	inform(message);
	sprintf (message,"%8ld varieties read", srf -> n_variety);
	inform(message);
	return (variety_handles);
}

/* fix this code */
int unpack_atom_numbers (struct vtybin *vtyin, struct variety *vty)
{
	int j, n;
	unsigned short an;
	struct packer packers[3];

	/* initialize */
	n = 0; 
	for (j = 0; j < 3; j++)
		packers[j].shortlong.onelong = 0L;
	/* store the numbers in local structure */
	for (j = 0; j < 3; j++) {
		packers[j].shortlong.onelong = vtyin -> atmnum[j];
	}
	/* transfer number from local structure */
	for (j = 0; j < MAXPA; j++) {
		if (j < 3) {
			an = packers[j].shortlong.twoshort[1];
		}
		else {
			an = packers[j-3].shortlong.twoshort[0];
		}
		if (an == 0) continue;
		n++;
		vty -> atmnum[j] = an;
	}
	vty -> natom = (short) n;
	return (n);
}


struct vertex **read_vertices (FILE *fp_pqms, struct surface *srf)
{
	int k;
	long i;
	char message[MAXLINE];
	struct vertex **vertex_handles;
	struct vertex *vtx;
	struct vtxshort vtxi;
	struct vtxbin vtxin;

	if (srf -> n_vertex > 0) {
		vertex_handles = (struct vertex **)
			allocate_pointers (VERTEX, srf -> n_vertex);
		if (vertex_handles == (struct vertex **) NULL) {
			set_error1 ("read_vertices: memory allocation failure");
			set_error2 ("for vertex handles");
			return (NULL);
		}
	}
	else vertex_handles = (struct vertex **) NULL;
	for (i = 0; i < srf -> n_vertex; i++) {
		vtx = allocate_vertex ();
		if (vtx == (struct vertex *) NULL) {
			set_error1 ("read_surface: vertex allocation failure");
			return (NULL);
		}
		*(vertex_handles + i) = vtx;
	}
	/* set up head and tail pointers */
	srf -> head_vertex = NULL;
	srf -> tail_vertex = NULL;
	/* read vertices */
	for (i = 0; i < srf -> n_vertex; i++) {
		if (shorty) {
			fread ((char *) &vtxi, sizeof (vtxi), 1, fp_pqms);
			vtxin.type = vtxi.type;
			for (k = 0; k < 3; k++)
				vtxin.center[k] = vtxi.center[k];
		}
		else fread ((char *) &vtxin, sizeof (vtxin), 1, fp_pqms);
		if (feof (fp_pqms)) {
			set_error1 ("read_vertices: premature eof on surface file");
			return (NULL);
		}
		vtx = *(vertex_handles + i);
		/* set up head and tail pointers for linked list */
		if (i == 0) srf -> head_vertex = vtx;
		srf -> tail_vertex = vtx;
		if (vtxin.type != VERTEX_TYPE) {
			set_error1 ("invalid vertex type");
			return (NULL);
		}
		for (k = 0; k < 3; k++)
			vtx -> center[k] = vtxin.center[k];
		vtx -> next = ((i < srf -> n_vertex - 1) ? *(vertex_handles + i + 1) : NULL);
		vtx -> number = i + 1;
		vtx -> cusp = 1;	/* initialize all to cusp */
	}
	sprintf (message,"%8ld vertices read", srf -> n_vertex);
	inform(message);
	return (vertex_handles);
}

struct circle **read_circles (FILE *fp_pqms, struct surface *srf)
{
	int k;
	long i;
	char message[MAXLINE];
	struct circle **circle_handles;
	struct cirbin cirin;
	struct cirshort ciri;
	struct circle *cir;

	if (srf -> n_circle > 0) {
		circle_handles = (struct circle **)
			allocate_pointers (CIRCLE, srf -> n_circle);
		if (circle_handles == (struct circle **) NULL) {
			set_error1 ("read_surface: memory allocation failure");
			set_error2 ("for circle handles");
			return (NULL);
		}
	}
	else circle_handles = (struct circle **) NULL;
	for (i = 0; i < srf -> n_circle; i++) {
		cir = allocate_circle ();
		if (cir == (struct circle *) NULL) {
			set_error1 ("read_circles: circle allocation failure");
			return (NULL);
		}
		*(circle_handles + i) = cir;
	}
	/* set up head and tail pointers */
	srf -> head_circle = NULL;
	srf -> tail_circle = NULL;
	/* read circles */
	for (i = 0; i < srf -> n_circle; i++) {
		if (shorty) {
			fread ((char *) &ciri, sizeof (ciri), 1, fp_pqms);
			cirin.type = ciri.type;
			cirin.subtype = ciri.subtype;
			for (k = 0; k < 3; k++) {
				cirin.center[k] = ciri.center[k];
				cirin.axis[k] = ciri.axis[k];
			}
			cirin.radius = ciri.radius;
		}
		else fread ((char *) &cirin, sizeof (cirin), 1, fp_pqms);
		if (feof (fp_pqms)) {
			set_error1 ("read_circles: premature eof on surface file");
			return (NULL);
		}
		cir = *(circle_handles + i);
		/* set up head and tail pointers for linked list */
		if (i == 0) srf -> head_circle = cir;
		srf -> tail_circle = cir;
		if (cirin.type != CIRCLE_TYPE) {
			set_error1 ("invalid circle type");
			return (NULL);
		}
		cir -> subtype = cirin.subtype;
		cir -> next = ((i < srf -> n_circle - 1) ? *(circle_handles + i + 1) : NULL);
		for (k = 0; k < 3; k++) {
			cir -> center[k] = cirin.center[k];
			cir -> axis[k] = cirin.axis[k];
		}
		cir -> radius = cirin.radius;
	}
	sprintf (message,"%8ld circles read", srf -> n_circle);
	inform(message);
	return (circle_handles);
}

struct arc **read_arcs (FILE *fp_pqms, struct surface *srf, struct circle **circle_handles, struct vertex **vertex_handles)
{
	int k;
	long i;
	long cirnum, vtxnum;
	double vector1[3], vector2[3];
	char message[MAXLINE];
	struct arc **arc_handles;
	struct arcbin arcin;
	struct arcshort arci;
	struct vertex *arc_vtx;
	struct vertex *vtx1, *vtx2;
	struct circle *cir;
	struct arc *a;

	if (srf -> n_arc > 0) {
		arc_handles = (struct arc **)
			allocate_pointers (ARC, srf -> n_arc);
		if (arc_handles == (struct arc **) NULL) {
			set_error1 ("read_arcs: memory allocation failure");
			set_error2 ("for arc handles");
			return (NULL);
		}
	}
	else arc_handles = (struct arc **) NULL;
	for (i = 0; i < srf -> n_arc; i++) {
		a = allocate_arc ();
		if (a == (struct arc *) NULL) {
			set_error1 ("read_arcs: arc allocation failure");
			return (NULL);
		}
		*(arc_handles + i) = a;
	}
	/* Note: conversion form index to pointer requires decrement */
	/* set up head and tail pointers */
	srf -> head_arc = NULL;
	srf -> tail_arc = NULL;
	/* read arcs */
	for (i = 0; i < srf -> n_arc; i++) {
		if (shorty) {
			fread ((char *) &arci, sizeof (arci), 1, fp_pqms);
			arcin.type = arci.type;
			arcin.subtype = arci.subtype;
			arcin.cirnum = arci.cirnum;
			for (k = 0; k < 2; k++)
				arcin.vtxnum[k] = arci.vtxnum[k];
			arcin.error = arci.error;
		}
		else fread ((char *) &arcin, sizeof (arcin), 1, fp_pqms);
		if (feof (fp_pqms)) {
			set_error1 ("read_arcs: premature eof on surface file");
			return (NULL);
		}
		a = *(arc_handles + i);
		/* set up head and tail pointers for linked list */
		if (i == 0) srf -> head_arc = a;
		srf -> tail_arc = a;
		if (arcin.type == CONVEX_ARC_TYPE) a -> shape = CONVEX;
		else if (arcin.type == CONCAVE_ARC_TYPE) a -> shape = CONCAVE;
		else {
			set_error1 ("invalic arc type");
			return (NULL);
		}
		/* index to pointer */
		cirnum = arcin.cirnum;
		if (cirnum < 1 || cirnum > srf -> n_circle) {
			sprintf (message, "read_surface: invalid arc-circle number %8ld", cirnum);
			set_error1 (message);
			return (NULL);
		}
		a -> cir = *(circle_handles + cirnum - 1);
		for (k = 0; k < 2; k++) {
			vtxnum = arcin.vtxnum[k];
			if (vtxnum < 0 || vtxnum > srf -> n_vertex) {
				sprintf (message, "read_surface: invalid arc-vertex number %8ld", vtxnum);
				set_error1 (message);
				return (NULL);
			}
			a -> vtx[k] = ((vtxnum == 0) ? NULL : *(vertex_handles + vtxnum - 1));
		}

		/* compute phi for omega */
		vtx1 = a -> vtx[0];
		vtx2 = a -> vtx[1];
		cir = a -> cir;
		if (cir -> radius <= 0.0) a -> phi = 0.0;
		else if (vtx1 != NULL && vtx2 != NULL) {
			for (k = 0; k < 3; k++) {
				vector1[k] = (vtx1 -> center[k] - cir -> center[k]) / cir -> radius;
				vector2[k] = (vtx2 -> center[k] - cir -> center[k]) / cir -> radius;
			}
			a -> phi = positive_angle (vector1, vector2, cir -> axis);
		}
		else a -> phi = 2 * PI;
		a -> next = ((i < srf -> n_arc - 1) ? *(arc_handles + i + 1) : NULL);
		a -> number = i+1;
		a -> original = 1;
	}
	/* mark those cusp vertices that show up in no arc */
	/* these vertices are in saddle cones */
	for (a = srf -> head_arc; a != NULL; a = a -> next)
		for (k = 0; k < 2; k++) {
			arc_vtx = a -> vtx[k];
			if (arc_vtx != NULL)
					arc_vtx -> cusp = 0;
		}
	sprintf (message,"%8ld arcs read", srf -> n_arc);
	inform(message);
	return (arc_handles);
}

struct face **read_faces (FILE *fp_pqms, struct surface *srf, long *first_cycle_number, struct variety **variety_handles)
{
	long i;
	char message[MAXLINE];
	struct face **face_handles;
	struct facbin facin;
	struct facshort faci;
	struct face *fac;

	face_handles = (struct face **)
		allocate_pointers (FACE, srf -> n_face);
	if (face_handles == (struct face **) NULL) {
		set_error1 ("read_faces: memory allocation failure");
		set_error2 ("for face handles");
		return (NULL);
	}
	for (i = 0; i < srf -> n_face; i++) {
		fac = allocate_face ();
		if (fac == (struct face *) NULL) {
			set_error1 ("read_faces: face allocation failure");
			return (NULL);
		}
		*(face_handles + i) = fac;
	}
	/* read faces */
	for (i = 0; i < srf -> n_face; i++) {
		if (shorty) {
			fread ((char *) &faci, sizeof (faci), 1, fp_pqms);
			facin.type = faci.type;
			facin.vtynum = faci.vtynum;
			facin.fcynum = faci.fcynum;
			facin.comp = faci.comp;
			facin.error = faci.error;
		}
		else fread ((char *) &facin, sizeof (facin), 1, fp_pqms);
		if (feof (fp_pqms)) {
			set_error1 ("read_faces: premature eof on surface file");
			return (NULL);
		}
		fac = *(face_handles + i);
		fac -> srf = srf;
		/* set up head and tail pointers for linked list */
		if (i == 0) srf -> head_face = fac;
		srf -> tail_face = fac;
		fac -> vty = *(variety_handles + facin.vtynum - 1);
		if (facin.error < 0) {
			sprintf (message,"read_surface warning: face %ld (error flag set)", i+1);
			inform(message);
			fac -> problem = TRUE;
		}
		else if (facin.error > 0) {
			fac -> input_hue = facin.error;
		}
		switch ((int) facin.type) {
		case CONVEX_FACE_TYPE:
			fac -> shape = CONVEX;
			sprintf (message, "convex face variety number %d", facin.vtynum);
			informd2(message);
			break;
		case SADDLE_FACE_TYPE:
			fac -> shape = SADDLE;
			sprintf (message, "saddle face variety number %d", facin.vtynum);
			informd2(message);
			break;
		case CONCAVE_FACE_TYPE:
			fac -> shape = CONCAVE;
			sprintf (message, "concave face variety number %d", facin.vtynum);
			informd2(message);
			break;
		default:
			set_error1 ("invalid face type");
			return (NULL);
		}
		*(first_cycle_number + i) = facin.fcynum;
		fac -> next = ((i < srf -> n_face - 1) ? *(face_handles + i + 1) : NULL);
		fac -> comp = (unsigned char) facin.comp;
		fac -> original = 1;
		fac -> largeok = 1;			/* ok to use add-circle routine for large solid angle */
	}
	sprintf (message,"%8ld faces read", srf -> n_face);
	inform(message);
	return (face_handles);
}

struct cycle **read_cycles (FILE *fp_pqms, struct surface *srf, short *n_edge_per_cycle, long *first_edge_number)
{
	long i;
	char message[MAXLINE];
	struct cycle **cycle_handles;
	struct cycbin cycin;
	struct cycshort cyci;
	struct cycle *cyc;

	if (srf -> n_cycle > 0) {
		cycle_handles = (struct cycle **)
			allocate_pointers (CYCLE, srf -> n_cycle);
		if (cycle_handles == (struct cycle **) NULL) {
			set_error1 ("read_cycles: memory allocation failure");
			set_error2 ("for cycle handles");
			return (NULL);
		}
	}
	else {
		cycle_handles = (struct cycle **) NULL;
	}
	for (i = 0; i < srf -> n_cycle; i++) {
		cyc = allocate_cycle ();
		if (cyc == (struct cycle *) NULL) {
			set_error1 ("read_cycles: cycle allocation failure");
			return (NULL);
		}
		*(cycle_handles + i) = cyc;
	}
	/* read cycles */
	for (i = 0; i < srf -> n_cycle; i++) {
		if (shorty) {
			fread ((char *) &cyci, sizeof (cyci), 1, fp_pqms);
			cycin.type = cyci.type;
			cycin.next = cyci.next;
			cycin.fednum = cyci.fednum;
			cycin.ned = cyci.ned;
		}
		else fread ((char *) &cycin, sizeof (cycin), 1, fp_pqms);
		if (feof (fp_pqms)) {
			set_error1 ("read_cycles: premature eof on surface file");
			return (NULL);
		}
		cyc = *(cycle_handles + i);
		if (cycin.type != CYCLE_TYPE) {
			set_error1 ("read_cycles: invalid cycle type");
			return (NULL);
		}
		cyc -> next = ((cycin.next == 0) ? NULL :
			*(cycle_handles + cycin.next - 1));
		*(n_edge_per_cycle + i) = cycin.ned;
		*(first_edge_number + i) = cycin.fednum;
	}
	sprintf (message,"%8ld cycles read", srf -> n_cycle);
	inform(message);
	return (cycle_handles);
}

struct edge **read_edges (FILE *fp_pqms, struct surface *srf, short *n_edge_per_cycle, struct arc **arc_handles)
{
	int orn;
	long i, j, e, arcnum, ne;
	struct edge *edg;
	struct arc *a;
	char message[MAXLINE];
	struct edge **edge_handles;
	struct edgbin edgin;
	struct edgshort edgi;

	if (srf -> n_edge > 0) {
		edge_handles = (struct edge **)
			allocate_pointers (EDGE, srf -> n_edge);
		if (edge_handles == (struct edge **) NULL) {
			set_error1 ("read_edges: memory allocation failure");
			set_error2 ("for edge handles");
			return (NULL);
		}
	}
	else edge_handles = (struct edge **) NULL;
	for (i = 0; i < srf -> n_edge; i++) {
		edg = allocate_edge ();
		if (edg == (struct edge *) NULL) {
			set_error1 ("read_edges: edge allocation failure");
			return (NULL);
		}
		*(edge_handles + i) = edg;
	}
	/* read edges */
	e = 0;
	for (i = 0; i < srf -> n_cycle; i++) {
		ne = (*(n_edge_per_cycle+i));
		for (j = 0; j < ne; j++) {
			if (e >= srf -> n_edge) {
				set_error1 ("read_edges: inconsistent number of edges");
				return (NULL);
			}
			if (shorty) {
				fread ((char *) &edgi, sizeof (edgi), 1, fp_pqms);
				edgin.type = edgi.type;
				edgin.error = edgi.error;
				edgin.arcnum = edgi.arcnum;
			}
			else fread ((char *) &edgin, sizeof (edgin), 1, fp_pqms);
			if (feof (fp_pqms)) {
				set_error1 ("read_edges: premature eof on surface file");
				return (NULL);
			}
			edg = *(edge_handles + e);
			if (edgin.type != EDGE_TYPE) {
				set_error1 ("invalid edge type");
				return (NULL);
			}
			/* pointer to next edge in cycle */
			edg -> next = ((j < ne - 1) ? *(edge_handles + e + 1) : NULL);
			/* decode */
			orn = ((edgin.arcnum > 0) ? 0 : 1);
			arcnum = abs ((int) edgin.arcnum);
			if (arcnum <= 0 || arcnum > srf -> n_arc) {
				sprintf (message, "arcnum %8ld invalid", arcnum);
				set_error1 (message);
				return (NULL);
			}
			edg -> orn = (short) orn;
			/* pointers: arc <---> edge */
			a = *(arc_handles + arcnum - 1);
			edg -> arcptr = a;
			a -> edg[orn] = edg;
			e++;
		}
	}
	sprintf (message,"%8ld edges read", srf -> n_edge);
	inform(message);
	return (edge_handles);
}

struct component **read_components (FILE *fp_pqms, struct surface *srf)
{
	int k;
	long i, n;
	char message[MAXLINE];
	struct component **component_handles;
	struct cmpbin cmpin;
	struct cmpshort cmpi;
	struct component *cmp_ptr;

	component_handles = (struct component **)
		allocate_pointers (COMPONENT, srf -> n_component);
	if (component_handles == (struct component **) NULL) {
		set_error1 ("read_components: memory allocation failure");
		set_error2 ("for component handles");
		return (NULL);
	}
	for (i = 0; i < srf -> n_component; i++) {
		cmp_ptr = allocate_component ();
		if (cmp_ptr == (struct component *) NULL) {
			set_error1 ("read_components: component allocation failure");
			return (NULL);
		}
		*(component_handles + i) = cmp_ptr;
	}
	/* read components */
	for (i = 0; i < srf -> n_component; i++) {
		if (shorty) {
			n = (long) fread ((void *) &cmpi, (size_t) sizeof (cmpi), (size_t) 1, fp_pqms);
			cmpin.type = cmpi.type;
			cmpin.subtype = cmpi.subtype;
			cmpin.volume = cmpi.volume;
			cmpin.area = cmpi.area;
			cmpin.accessible = cmpi.accessible;
			for (k = 0; k < 3; k++)
				cmpin.center[k] = cmpi.center[k];
		}
		else {
			n = (long) fread ((void *) &cmpin, (size_t) sizeof (cmpin),
			(size_t) 1, fp_pqms);
		}
		if (feof (fp_pqms)) {
			set_error1 ("read_components: premature eof on surface file");
			return (NULL);
		}
		if (ferror (fp_pqms)) {
			set_error1 ("read_components: error reading surface file");
			return (NULL);
		}
		if (cmpin.type != COMPONENT_TYPE) {
			sprintf (message,
				"read_surface: invalid component (%5ld) type: %5ld, n = %5ld",
					i + 1, (long) cmpin.type, n);
			set_error1(message);
			return (NULL);
		}
		cmp_ptr = *(component_handles + i);
		for (k = 0; k < 3; k++)
			cmp_ptr -> center[k] = cmpin.center[k];
		cmp_ptr -> volume = cmpin.volume;
		cmp_ptr -> area = cmpin.area;
		cmp_ptr -> accessible = cmpin.accessible;
		cmp_ptr -> type = cmpin.subtype;
	}
	sprintf (message,"%8ld components read", srf -> n_component);
	inform(message);
	return (component_handles);
}



