/*
	MSRoll
	Copyright 1986, 1989 by Michael L. Connolly
	All rights reserved
	December 26, 2001
*/

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"




/* write piecewise polynomial molecular surface to output file */
int write_surface (struct surface *srf, FILE *fp_surface)
{
	/* writing */
	write_header (srf, srf -> mol -> name, fp_surface); if (error()) return(0);
	write_varieties (srf, fp_surface); if (error()) return(0);
	write_vertices (srf, fp_surface); if (error()) return(0);
	write_circles (srf, fp_surface); if (error()) return(0);
	write_arcs (srf, fp_surface); if (error()) return(0);
	write_faces (srf, fp_surface); if (error()) return(0);
	write_cycles (srf, fp_surface); if (error()) return(0);
	write_edges (srf, fp_surface); if (error()) return(0);
	write_components (srf, fp_surface); if (error()) return(0);

	fclose (fp_surface);
	return (1);
}


/* write header record */
void write_header (struct surface *srf, char *molecule_name, FILE *fp_surface)
{
	long idx;
	struct header headerout;

	for (idx = 0; idx < 8; idx++)
		headerout.filetype[idx] = (char) 0;
	strcpy (headerout.filetype, "PQMS");
	headerout.version = VERSION;
	headerout.subversion = SUBVERSION;
	headerout.probe_radius = srf -> probe_radius;
	headerout.counters[0] = srf -> n_variety;
	headerout.counters[1] = srf -> n_vertex;
	headerout.counters[2] = srf -> n_circle;
	headerout.counters[3] = srf -> n_arc;
	headerout.counters[4] = srf -> n_face;
	headerout.counters[5] = srf -> n_cycle;
	headerout.counters[6] = srf -> n_edge;
	headerout.counters[7] = srf -> n_component;
	headerout.counters[8] = srf -> n_atom;
	headerout.counters[9] = 0;			/* unused */
	for (idx = 0; idx < 64; idx++)
		headerout.name[idx] = (char) 0;
	if (strlen (molecule_name) >= (unsigned) 64) {
		inform ("molecule name too long, omitted");
		strcpy (headerout.name, " ");
	}
	else strcpy (headerout.name, molecule_name);
	for (idx = 0; idx < 384; idx++)
		headerout.filler[idx] = (char) 0;
	fwrite ((char *) &headerout, sizeof (struct header), 1, fp_surface);
}

/* write varieties */
void write_varieties (struct surface *srf, FILE *fp_surface)
{
	int k, npacked;
	long v;
	struct variety *vty;
	struct vtybin vtyout;
	char message[MAXLINE];

	for (v = 0; v < srf -> n_variety; v++) {
		vty = *(srf -> variety_handles + v);
		if (v < srf -> n_atom)
			vtyout.type = ATOM_TYPE;
		else if (v < srf -> n_atom + srf -> n_tori)
			vtyout.type = TORUS_TYPE;
		else vtyout.type = PROBE_TYPE;
		vtyout.radius = vty -> radii[0];
		for (k = 0; k < 3; k++) {
			vtyout.center[k] = vty -> center[k];
			vtyout.axis[k] = vty -> axis[k];
			/* vtyout.atmnum[k] = vty -> atmnum[k]; */
		}
		npacked = pack_atom_numbers (vty, &vtyout);
		if (npacked > 3) {
			sprintf (message, "%8ld variety has %d atoms",
				v + 1, npacked);
			informd (message);
		}
		fwrite ((char *) &vtyout, sizeof (vtyout), 1, fp_surface);
	}
	sprintf (message, "%8ld varieties written", srf -> n_variety);
	informd (message);
}

int pack_atom_numbers (struct variety *vty, struct vtybin *vtyout) 
{
	int j, n;
	unsigned short an;
	struct packer packers[3];

	/* initialize */
	n = 0; 
	for (j = 0; j < 3; j++)
		packers[j].shortlong.onelong = 0L;
	/* store the numbers in local structure */
	for (j = 0; j < MAXPA; j++) {
		an = (unsigned short) vty -> atmnum[j];
		if (an == 0) continue;
		n++;
		if (j < 3) {
			packers[j].shortlong.twoshort[1] = an;
		}
		else {
			packers[j-3].shortlong.twoshort[0] = an;
		}
	}
	/* transfer number from local structure */
	for (j = 0; j < 3; j++) {
		vtyout -> atmnum[j] = packers[j].shortlong.onelong;
	}
	return (n);
}

/* write vertices */
void write_vertices (struct surface *srf, FILE *fp_surface)
{
	int k;
	long n_vertex;
	char message[MAXLINE];
	struct vertex *vtx;
	struct vtxbin vtxout;

	n_vertex = 0;
	vtxout.type = VERTEX_TYPE;
	for (vtx = srf -> head_vertex; vtx != NULL; vtx = vtx -> next) {
		for (k = 0; k < 3; k++)
			vtxout.center[k] = vtx -> center[k];
		fwrite ((char *) &vtxout, sizeof (vtxout), 1, fp_surface);
		n_vertex++;
	}
	sprintf (message, "%8ld vertices written", n_vertex);
	informd (message);
}

/* write circles */
void write_circles (struct surface *srf, FILE *fp_surface)
{
	int k;
	long n_circle;
	char message[MAXLINE];
	struct circle *cir;
	struct cirbin cirout;

	n_circle = 0;
	cirout.type = CIRCLE_TYPE;
	for (cir = srf -> head_circle; cir != NULL; cir = cir -> next) {
		for (k = 0; k < 3; k++) {
			cirout.center[k] = cir -> center[k];
			cirout.axis[k] = cir -> axis[k];
		}
		cirout.radius = cir -> radius;
		cirout.subtype = cir -> subtype;
		fwrite ((char *) &cirout, sizeof (cirout), 1, fp_surface);
		n_circle++;
	}
	sprintf (message, "%8ld circles written", n_circle);
	informd (message);
}

/* write arcs */
void write_arcs (struct surface *srf, FILE *fp_surface)
{
	int m;
	long l;
	char message[MAXLINE];
	struct arc *a;
	struct arcbin arcout;

	arcout.error = FALSE;
	for (l = 0; l < srf -> n_arc; l++) {
		a = *(srf -> arc_handles + l);
		arcout.error = FALSE;
		if (a -> shape == CONVEX)
			arcout.type = CONVEX_ARC_TYPE;
		else if (a -> shape == CONCAVE)
			arcout.type = CONCAVE_ARC_TYPE;
		else {
			set_error1 ("write_arcs: invalid shape");
			return;
		}
		arcout.cirnum = a -> cir -> number;
		for (m = 0; m < 2; m++)
			arcout.vtxnum[m] =
				((a -> vtx[m] == NULL) ? 0 : a -> vtx[m] -> number);
		fwrite ((char *) &arcout, sizeof (arcout), 1, fp_surface);
	}
	sprintf (message, "%8ld arcs written", l);
	informd (message);
}


/* write face */
void write_faces (struct surface *srf, FILE *fp_surface)
{
	int icyc;
	long n_face;
	char message[MAXLINE];
	struct face *fac;
	struct facbin facout;
	struct variety *vty;
    struct cept *ex;

	n_face = 0;
	icyc = 1;
	for (fac =  srf -> head_face; fac != NULL; fac = fac -> next) {
		vty = fac -> vty;
		if (vty == NULL) {
			ex = new_cept (POINTER_ERROR, NULL_VALUE, FATAL_SEVERITY);
			add_function (ex, "write_faces");
			add_source (ex, "msrollio.c");
			return;
		}
		facout.error = 0;
		facout.vtynum = vty -> number;
		sprintf (message, "%8ld face %8ld variety %1d shape",
			n_face, vty -> number, (int) fac -> shape);
		informd (message);
		switch ((int) fac -> shape) {
		case CONVEX:
			facout.type = CONVEX_FACE_TYPE;
			facout.fcynum = ((fac -> n_cycle > 0) ? icyc : 0);
			icyc += fac -> n_cycle;
			break;
		case SADDLE:
			facout.type = SADDLE_FACE_TYPE;
			facout.fcynum =  icyc;
			if (fac -> n_arc == 2) icyc += 2;
			else icyc++;
			break;
		case CONCAVE:
			facout.type =  CONCAVE_FACE_TYPE;
			/* problem concave face */
			if (fac -> problem) facout.error = 1;
			facout.fcynum =  icyc;
			if (fac -> ptr.prb -> low)
				icyc += fac -> n_cycle;
			else icyc++;
			break;
		default:
			set_error1 ("write_faces: invalid face shape");
			return;
		}
		facout.comp = fac -> comp;
		fwrite ((char *) &facout, sizeof (facout), 1, fp_surface);
		n_face++;
	}
	sprintf (message, "%8ld faces written", n_face);
	informd (message);
} 

/* write cycles */
void write_cycles (struct surface *srf, FILE *fp_surface)
{
	long  iedg;
	long y;
	char message[MAXLINE];
	struct cycle *cyc;
	struct cycbin cycout;
	

	iedg = 1;
	for (y = 0; y < srf -> n_cycle; y++) {
		cyc = *(srf -> cycle_handles + y);
		cycout.type = CYCLE_TYPE;
		cycout.next = ((cyc -> next == NULL) ? 0 : y + 2);
		cycout.fednum = iedg;
		cycout.ned = edges_in_cycle (cyc);
		iedg += edges_in_cycle (cyc);
		fwrite ((char *) &cycout, sizeof (cycout), 1, fp_surface);
	}
	if (iedg != srf -> n_edge + 1) {
		set_error1 ("(write_cycles): edge count inconsistency");
		return;
	}
	sprintf (message, "%8ld cycles written", y);
	informd (message);
}

/* write edges */
void write_edges (struct surface *srf, FILE *fp_surface)
{
	int sgn;
	long e;
	long number;
	struct arc *a;
	struct edge *edg;
	struct edgbin edgout;
	char message[MAXLINE];

	for (e = 0; e < srf -> n_edge; e++) {
		edg = *(srf -> edge_handles + e);
		a = edg -> arcptr;
		if (a -> number <= 0 || a -> number > srf -> n_arc) {
			sprintf (message, "arc number %8ld invalid", a -> number);
			set_error1 (message);
			return;
		}
		sgn = (edg -> orn) ? -1 : 1;
		number = sgn * a -> number;
		if (number == 0) {
			set_error1 ("write_edges: invalid arc number");
			return;
		}
		edgout.type = EDGE_TYPE;
		edgout.error = 0;
		edgout.arcnum = number;
		fwrite ((char *) &edgout, sizeof (edgout), 1, fp_surface);
	}
	sprintf (message, "%8ld edges written", srf -> n_edge);
	informd (message);
}

/* write components */
void write_components (struct surface *srf, FILE *fp_surface)
{
	int k, comp;
	char message[MAXLINE];
	struct component *cmp_ptr;
	struct cmpbin cmpout;

	cmpout.type = COMPONENT_TYPE;

	for (comp = 0; comp < srf -> n_component; comp++) {
		cmp_ptr = get_component_ptr (srf, comp+1);
		cmpout.subtype = (cmp_ptr -> volume > 0.0) ? OUTER_SUBTYPE : INNER_SUBTYPE;
		for (k = 0; k < 3; k++)
			cmpout.center[k] = cmp_ptr -> center[k];
		cmpout.volume = cmp_ptr -> volume;
		cmpout.area = cmp_ptr -> area;
		cmpout.accessible = cmp_ptr -> accessible;
		fwrite ((void *) &cmpout, (size_t) sizeof (cmpout), (size_t) 1, fp_surface);
	}
	sprintf (message, "%8ld components written", srf -> n_component);
	informd (message);
}


