/*

	MSRoll

	Copyright 1986, 1989, 1996 by Michael L. Connolly
	All rights reserved

	Written by Michael L. Connolly.
	March 7, 2000

*/

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"



int triangulate (struct msscene *ms)
{
	presub(ms);
	if (error ()) return (0);
	oncethru(ms);
	if (error ()) return (0);
	postsub(ms);
	if (error ()) return (0);
	return (1);
}

int presub (struct msscene *ms)
{
	unsigned long fidx;
	char message[MAXLINE];
	struct face *face_ptr;
	struct arc *arc_ptr;
	struct edge *e;
	struct cycle *cyc;
	int n_lune;

	if (ms -> this_srf == (struct surface *) NULL) {
		inform ("premature triangulate command");
		return(0);
	}
	if (ms -> this_srf -> triangulation_completed) {
		inform ("triangulation already performed for this surface");
		return(0);
	}

	sprintf (message,"%8.3f distance weight parameter",
		ms -> this_srf -> weight);
	informd(message);
	sprintf (message,"%8.3f large omega parameter",
		ms -> this_srf -> large_omega);
	informd(message);

	/* save original number of arcs */
	ms -> this_srf -> n_arc_original = ms -> this_srf -> n_arc;
	if (!setup_ffn (ms -> this_srf)) return(0);
	
	/* allocate head and tail pointers for triangles */
	ms -> this_srf -> heads = (struct phntri **)
		allocate_pointers (PHNTRI, ms -> this_srf -> n_atom);
	if (ms -> this_srf -> heads == NULL) {
		set_error1 ("triangulate: not enough memory for heads");
		return (0);
	}
	ms -> this_srf -> tails = (struct phntri **)
		allocate_pointers (PHNTRI, ms -> this_srf -> n_atom);
	if (ms -> this_srf -> tails == NULL) {
		set_error1 ("triangulate: not enough memory for tails");
		return (0);
	}

	/* mark lunes */
	n_lune = 0;
	for (fidx = 0; fidx < ms -> this_srf -> n_face; fidx++) {
		face_ptr = *(ms -> this_srf -> face_handles+fidx);
		for (cyc = face_ptr -> first_cycle; cyc != NULL; cyc = cyc -> next) {
			if (edges_in_cycle (cyc) == 2) {
				for (e = cyc -> first_edge; e != NULL; e = e -> next) {
					arc_ptr = e -> arcptr;
					arc_ptr -> lune = 1;
				}
				n_lune++;
			}
		}
	}
	if (n_lune > 0) {
		sprintf (message, "%8d lunes", n_lune);
		inform(message);
	}
	ms -> this_srf -> n_face_original = ms -> this_srf -> n_face;
	sprintf (message,
		"%8ld faces to triangulate", ms -> this_srf -> n_face_original);
	informd(message);
	return (1);
}

int oncethru (struct msscene *ms)
{
	unsigned long fidx;

	for (fidx = 0; fidx < ms -> this_srf -> n_face_original; fidx++) {
		onceface (ms, fidx);
		if (error()) return (0);
	}
	ms -> this_srf -> triangulation_completed = 1;
	return (1);
}

/* bisect this face (recursively) */
int onceface (struct msscene *ms, unsigned long fidx)
{
	int shape;
	unsigned long lfn, ofn;
	double percent;
	char message[MAXLINE];
	struct face *face_ptr;

	ofn = fidx + 1;
	oncearc (ms, fidx);
	if (error()) return (0);
	face_ptr = *(ms -> this_srf -> face_handles+fidx);
	shape = face_ptr -> shape;
	/* for VDW polyhedron will not be connected */
	if (ms -> this_srf -> van_der_Waals) {
		if (shape == CONCAVE) concave_degenerate(ms -> this_srf, face_ptr);
		else if (shape == SADDLE) saddle_degenerate(ms -> this_srf, face_ptr);
		if (error()) return(0);
		if (shape != CONVEX) return(1);
	}
	percent = (double) fidx / (double) (ms -> this_srf -> n_face_original);
	percent *= 100.0;
	switch ( shape ) {
	case CONVEX:
		sprintf (message,"%3.0f%% +%ld", percent, fidx+1);
		break;
	case SADDLE:
		sprintf (message,"%3.0f%% $%ld", percent, fidx+1);
		break;
	case CONCAVE:
		sprintf (message,"%3.0f%% -%ld", percent, fidx+1);
		break;
	case CYLINDRICAL:
		sprintf (message,"%3.0f%% |%ld", percent, fidx+1);
		break;
	default:
		set_error1 ("msroll (do_command): invalid face shape");
		return (0);
	}
	informd2(message);
	bisect_face (face_ptr);
	if (error()) return (0);
	lfn = fidx + 1;
	/* later, make the frequency user-setable */
	if (ms -> purge_frequency > 0 &&
		(lfn % (long) (ms -> purge_frequency)) == 0)
		if (!purge_surface (ms -> this_srf, lfn)) return(0);
	return (1);
}

int oncearc (struct msscene *ms, unsigned long fidx)
{
	unsigned long arc_idx;
	unsigned long ofn;
	int shape;
	long n_subd;
	struct arc *arc_ptr;
	int n_point = 0;
	ofn = fidx + 1;
	/* subdivide some of the arcs */
	n_subd = 0;
	for (arc_idx = 0; arc_idx < ms -> this_srf -> n_arc_original; arc_idx++) {
		arc_ptr = *(ms -> this_srf -> arc_handles+arc_idx);
		/* check for released */
		if (arc_ptr == (struct arc *) NULL) continue;
		/* check for whether we need the arc yet */
		if (arc_ptr -> ffn > ofn) continue;
		/* check for whether we are beyond it */
		if (arc_ptr -> lfn < ofn) continue;
		/* check for already subdivided */
		if (arc_ptr -> subdivided) continue;
		shape = arc_ptr -> shape;
		if (ms -> this_srf -> van_der_Waals && shape != CONVEX) continue;
		n_point = subdivide_arc (ms -> this_srf, arc_ptr); n_subd++;
		if (error()) return(0);
	}
	return (1);
}

int postsub (struct msscene *ms)
{
	long removed1, removed2;
	char message[MAXLINE];

	sprintf (message, "%8ld entire sphere faces",
		ms -> this_srf -> n_entire_face);
	inform(message);
	sprintf (message, "%8ld large convex faces",
		ms -> this_srf -> n_large_face);
	inform(message);
	sprintf (message, "%8ld simplified faces",
		ms -> this_srf -> n_simplified_face);
	inform(message);
	sprintf (message,
		"%8ld triangles %8ld edges %8ld vertices",
		ms -> this_srf -> n_face, ms -> this_srf -> n_arc, ms -> this_srf -> n_vertex);
	inform(message);
	convert_vet (ms -> this_srf);
	if (error()) return (0);
	removed1 = purifyPolyhedron (ms -> this_srf, ms -> coalesce);
	if (error()) return (0);
	sprintf (message,
		"%8ld vertices, edges and triangles removed in first pass", removed1);
	inform (message);
	if (removed1 > 0) {
		removed2 = purifyPolyhedron (ms -> this_srf, ms -> coalesce);
		if (error()) return (0);
		sprintf (message,
			"%8ld vertices, edges and triangles removed in second pass", removed2);
		inform (message);
	}
	if (!measure_polyhedron (ms -> this_srf)) return(0);
	return (1);
}



/*
   MSRoll
   Copyright 1996 by Michael L. Connolly
   All Rights Reserved

*/
