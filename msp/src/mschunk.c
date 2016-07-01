#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

/* Copyright 1995 by Michael L. Connolly */
/* December 13, 2001 */


/* CHUNK routines */

/* gather atoms genuinely participating in reentrant face */

int gather (struct face *fac, int atom_numbers[MAXPA])
{
	long i, n, j, atom_number, there;
	struct circle *cir;
	struct cycle *cyc;
	struct edge *edg;
	struct arc *a;
	struct vertex *vtx;
	struct sphere *atm;

	for (i = 0; i < MAXPA; i++)
		atom_numbers[i] = 0;
	n = 0;

	for (cyc = fac -> first_cycle; cyc != NULL; cyc = cyc -> next)
		for (edg = cyc -> first_edge; edg != NULL; edg = edg -> next) {
			a = edg -> arcptr;
			cir = a -> cir;
			atm = cir -> atm;
			if (atm != NULL) {
				atom_number = atm -> number;
				/* check whether already present */
				there = 0;
				for (i = 0; i < n; i++)
					if (atom_number == atom_numbers[i]) {
						there = 1;
						break;
					}
				if (!there) {
					if (n >= MAXPA) {
						set_error1 ("(aliquot): atom overflow");
						return(0);
					}
					/* add to array */
					atom_numbers[n++] = atom_number;
				}
			}
			for (j = 0; j < 2; j++) {
				vtx = a -> vtx[j];
				if (vtx == NULL) continue;
				atm = vtx -> atm;
				if (atm == NULL) continue;
				atom_number = atm -> number;
				/* check whether already present */
				there = 0;
				for (i = 0; i < n; i++)
					if (atom_number == atom_numbers[i]) {
						there = 1;
						break;
					}
				if (there) continue;
				if (n >= MAXPA) {
					set_error1 ("(aliquot): atom overflow");
					return(0);
				}
				/* add to array */
				atom_numbers[n++] = atom_number;
			}	/* end of j loop */
		}	/* end of edg loop */
	return (n);
}

double aliquot (struct face *fac, int atom_number)
{
	int i, n_numbers;
	int atom_numbers[MAXPA];
	double area_quot;

	n_numbers = gather (fac, atom_numbers); 
	if (error()) return ((double) 0.0);
	if (n_numbers <= 0) return ((double) 0.0);
	area_quot = fac -> area / n_numbers;
	for (i = 0; i < n_numbers; i++)
		if (atom_number == atom_numbers[i])
			return (area_quot);
	return ((double) 0.0);
}

void add_chunk (struct atom *atm_ptr, int component_number, double contact_area, double reentrant_area, double accessible_area)
{
	struct chunk *head_chunk, *this_chunk;
	struct chunk *previous_chunk;
	struct chunk *added_chunk;

	head_chunk = atm_ptr -> head_chunk;

	if (head_chunk == NULL) {
		/* start new linked list */
		added_chunk = new_chunk (atm_ptr, component_number,
			contact_area, reentrant_area, accessible_area);
		if (error()) return;
		atm_ptr -> head_chunk = added_chunk;
		return;
	}
	/* hunt down linked list */
	previous_chunk = NULL;
	for (this_chunk = head_chunk; this_chunk != NULL;
		this_chunk = this_chunk -> next[0]) {
		if (this_chunk -> component_number == component_number) {
			this_chunk -> contact_area += contact_area;
			this_chunk -> reentrant_area += reentrant_area;
			this_chunk -> accessible_area += accessible_area;
			return;
		}
		if (this_chunk -> component_number > component_number) break;
		previous_chunk = this_chunk;
	}
	/* link into list */
	added_chunk = new_chunk (atm_ptr, component_number,
		contact_area, reentrant_area, accessible_area);
	if (previous_chunk == NULL)
		atm_ptr -> head_chunk = added_chunk;
	else previous_chunk -> next[0] = added_chunk;
	added_chunk -> next[0] = this_chunk;
	return;
}

struct chunk *new_chunk (struct atom *atm_ptr, int component_number, double contact_area, double reentrant_area, double accessible_area)
{

	struct chunk *n_chunk;
	int atom_number;

	atom_number = atm_ptr -> number;

	n_chunk = (struct chunk *) allocate_chunk ();
	if (n_chunk == NULL) {
		set_error1 ("new_chunk: ran out of memory");
		return(NULL);
	}

	n_chunk -> atom_number = atom_number;
	strcpy (n_chunk -> labels[0], atm_ptr -> group);
	strcpy (n_chunk -> labels[1], atm_ptr -> sequence);
	strcpy (n_chunk -> labels[2], atm_ptr -> name);
	n_chunk -> component_number = (short) component_number;
	n_chunk -> contact_area = contact_area;
	n_chunk -> reentrant_area = reentrant_area;
	n_chunk -> accessible_area = accessible_area;
	return (n_chunk);
}

/* set up second set of links */

void setup_second_links (struct surface *srf)
{
	long component_index, component_number;
	struct chunk *head_chunk;
	struct chunk *this_chunk, *previous_chunk;
	struct component *cmp_ptr;
	struct atom *atm_ptr;

	for (component_index = 0; component_index < srf -> n_component;
		component_index++) {
		cmp_ptr = get_component_ptr (srf, component_index+1);
		component_number = component_index + 1;
		previous_chunk = NULL;
		for (atm_ptr = (struct atom *) (srf -> head_atom); atm_ptr != NULL;
			atm_ptr = atm_ptr -> next) {
			head_chunk = atm_ptr -> head_chunk;
			for (this_chunk = head_chunk; this_chunk != NULL;
				this_chunk = this_chunk -> next[0]) {
				if (this_chunk -> component_number !=
					component_number) continue;
				/* we like this one */
				if (previous_chunk == NULL)
					cmp_ptr -> head_chunk = this_chunk;
				else previous_chunk -> next[1] = this_chunk;
				previous_chunk = this_chunk;
				break;
			}	/* end of chunk loop */
		} /* end of atom loop */
	} /* end of component loop */
}


int write_area (struct surface *this_srf, FILE *fp_area, char *argument)
{
	if (! this_srf -> surface_completed || ! this_srf -> volume_computed) {
		inform ("premature area command");
		return(0);
	}
	if (argument == NULL || strlen(argument) < (unsigned) 1) {
		if (!write_atom_area (this_srf, fp_area)) return (0);
	}
	else {
		if (strcmp (argument, "by_atom_and_component") == 0 ||
			strcmp (argument, "bac") == 0) {
			if (!write_atom_component (this_srf, fp_area)) return (0);
		}
		else if (strcmp (argument, "by_component_and_atom") == 0
			|| strcmp (argument, "bca") == 0) {
			if (!write_component_atom (this_srf, fp_area)) return (0);
		}
		else {
			inform ("write_area: invalid breakdown type for area output");
			return(0);
		}
	}
	if (error()) return (0);
	return (1);
}


int write_atom_area (struct surface *this_srf, FILE *fp_area)
{
	int atom_number, label_present;
	double contact_area, reentrant_area, accessible_area;
	struct chunk *head_chunk, *this_chunk;
	struct atom *atm_ptr;

	for (atm_ptr = (struct atom *) (this_srf -> head_atom); atm_ptr != NULL;
		atm_ptr = atm_ptr -> next) {
		atom_number = atm_ptr -> number;
		contact_area = 0.0;
		accessible_area = 0.0;
		reentrant_area = 0.0;
		label_present = 0;
		head_chunk = atm_ptr -> head_chunk;
		for (this_chunk = head_chunk; this_chunk != NULL;
			this_chunk = this_chunk -> next[0]) {
			contact_area += this_chunk -> contact_area;
			reentrant_area += this_chunk -> reentrant_area;
			accessible_area += this_chunk -> accessible_area;
		}
		label_present = (strcmp (atm_ptr -> group, " ") != 0);
		/* write atom area output file record */
		if (! label_present)
			fprintf (fp_area, "%5d %8.3f %8.3f %8.3f %8.3f\n",
				atom_number, contact_area, reentrant_area,
				contact_area + reentrant_area,
				accessible_area);
		else fprintf (fp_area, "%-5s %-5s %-5s %8.3f %8.3f %8.3f %8.3f\n",
				atm_ptr -> group, atm_ptr -> sequence,
				atm_ptr -> name, contact_area, reentrant_area,
				contact_area + reentrant_area,
				accessible_area);
	}
	return (1);
}

int write_atom_component (struct surface *this_srf, FILE *fp_area)
{
	int atom_number, component_number, label_present;
	double contact_area, reentrant_area, accessible_area;
	struct chunk *head_chunk, *this_chunk;
	struct atom *atm_ptr;

	for (atm_ptr = (struct atom *) (this_srf -> head_atom); atm_ptr != NULL;
		atm_ptr = atm_ptr -> next) {
		atom_number = atm_ptr -> number;
		head_chunk = atm_ptr -> head_chunk;
		if (head_chunk == NULL) continue;
		for (this_chunk = head_chunk; this_chunk != NULL;
			this_chunk = this_chunk -> next[0]) {
			component_number = this_chunk -> component_number;
			contact_area = this_chunk -> contact_area;
			reentrant_area = this_chunk -> reentrant_area;
			accessible_area = this_chunk -> accessible_area;
			label_present = (strcmp (atm_ptr -> group, " ") != 0);
		if (! label_present)
			fprintf (fp_area, "%5d %3d %8.3f %8.3f %8.3f %8.3f\n",
				atom_number, component_number,
				contact_area, reentrant_area,
				contact_area + reentrant_area, accessible_area);
		else fprintf (fp_area,
				 "%-5s %-5s %-5s %3d %8.3f %8.3f %8.3f %8.3f\n",
				atm_ptr -> group, atm_ptr -> sequence,
				atm_ptr -> name, component_number,
				contact_area, reentrant_area,
				contact_area + reentrant_area, accessible_area);
		}
	}
	return (1);
}

int write_component_atom (struct surface *srf, FILE *fp_area)
{
	int atom_number, component_number, label_present;
	double contact_area, reentrant_area, accessible_area;
	struct chunk *head_chunk, *this_chunk;
	struct component *cmp_ptr;

	for (component_number = 1; component_number <= srf -> n_component;
		component_number++) {
		cmp_ptr = get_component_ptr (srf, component_number);
		head_chunk = cmp_ptr -> head_chunk;
		if (head_chunk == NULL) continue;
		for (this_chunk = head_chunk; this_chunk != NULL;
			this_chunk = this_chunk -> next[1]) {
			atom_number = this_chunk -> atom_number;
			contact_area = this_chunk -> contact_area;
			reentrant_area = this_chunk -> reentrant_area;
			accessible_area = this_chunk -> accessible_area;
			label_present = (strcmp (this_chunk -> labels[0], " ") != 0);
			if (!label_present)
				fprintf (fp_area, "%3d %5d %8.3f %8.3f %8.3f %8.3f\n",
					component_number, atom_number,
					contact_area, reentrant_area,
					contact_area + reentrant_area, accessible_area);
			else fprintf (fp_area,
					"%3d %-5s %-5s %-5s %8.3f %8.3f %8.3f %8.3f\n",
					component_number, this_chunk -> labels[0],
					this_chunk -> labels[1], this_chunk -> labels[2],
					contact_area, reentrant_area,
					contact_area + reentrant_area, accessible_area);
			}
	}
	return (1);
}


int write_volume1 (struct surface *srf, FILE *fp_volume)
{
	int comp;
	struct component *cmp_ptr;
	fprintf (fp_volume, "number of atoms         =     %6ld\n", srf -> n_atom);
	fprintf (fp_volume, "probe radius            =     %10.3f\n", srf -> probe_radius);
	fprintf (fp_volume, "contact area            =     %10.3f\n", srf -> total_contact_area);
	fprintf (fp_volume, "reentrant area          =     %10.3f\n", srf -> total_reentrant_area);
	fprintf (fp_volume, "molecular area          =     %10.3f\n", srf -> molecule_area);
	fprintf (fp_volume, "accessible area         =     %10.3f\n", srf -> total_accessible_area);
	fprintf (fp_volume, "solvent-excluded volume =     %10.3f\n", srf -> mol_vol);
	fprintf (fp_volume, "\n");
	fprintf (fp_volume, "component           center              volume    molecular  accessible\n");
	fprintf (fp_volume, " number      x        y        z                    area       area\n");
	fprintf (fp_volume, "\n");
	for (comp = 0; comp < srf -> n_component; comp++) {
		cmp_ptr = get_component_ptr (srf, comp+1);
		fprintf (fp_volume,
			"%5d    %8.3f %8.3f %8.3f %11.3f %11.3f %11.3f\n",
			comp+1,
			cmp_ptr -> center[0],cmp_ptr -> center[1],cmp_ptr -> center[2],
			cmp_ptr -> volume,cmp_ptr -> area,cmp_ptr -> accessible);
	}
	return(1);
}

int write_volume2 (struct surface *srf)
{
	int component_number;
	double taa, tav;
	char message[MAXLINE];
	struct component *cmp_ptr;

	if (srf == (struct surface *) NULL) {
		inform ("premature volume command");
		return(0);
	}
	sprintf (message,"%8.1f convex area  (polyhedral)", srf -> shape_area[0]);
	inform(message);
	sprintf (message,"%8.1f saddle area  (polyhedral)", srf -> shape_area[1]);
	inform(message);
	sprintf (message,"%8.1f concave area (polyhedral)", srf -> shape_area[2]);
	inform(message);
	taa = 0.0; tav = 0.0;
	inform ("                     area                         volume");
	inform ("              analytic  polyhedral          analytic  polyhedral");
	for (component_number = 1; component_number <= srf -> n_component;
		component_number++) {
		cmp_ptr = *(srf -> component_handles + component_number - 1);
		taa += cmp_ptr -> area;
		tav += cmp_ptr -> volume;
		sprintf (message,"%8d    %10.3f  %10.3f        %10.3f  %10.3f",
			component_number, cmp_ptr -> area, cmp_ptr -> parea,
			cmp_ptr -> volume, cmp_ptr -> pvolume);
		inform(message);
	}
	sprintf (message,"            %10.3f  %10.3f        %10.3f  %10.3f",
		taa, srf -> total_area, tav, srf -> total_volume);
	inform(message);
	return(1);
}


struct chunk *allocate_chunk ()
{
	struct chunk *chk;
	
	chk = (struct chunk *) allocate_object (CHUNK);
	return (chk);
}

void free_chunks (struct surface *srf)
{
	int component_index;
	long n_freed;
	char message[MAXLINE];
	struct chunk *head_chunk;
	struct chunk *this_chunk, *next_chunk;
	struct component *cmp_ptr;

	n_freed = 0;
	for (component_index = 0; component_index < srf -> n_component;
		component_index++) {
		cmp_ptr = get_component_ptr (srf, component_index+1);
		if (cmp_ptr == NULL) continue;
		head_chunk = cmp_ptr -> head_chunk;
		next_chunk = NULL;
		for (this_chunk = head_chunk; this_chunk != NULL;
			this_chunk = next_chunk) {
			next_chunk = this_chunk -> next[1];
			free_chunk (this_chunk); n_freed++;
		}	/* end of chunk loop */
	} /* end of component loop */

	sprintf (message,"%8ld chunks freed", n_freed);
	informd(message);
	free_cache (CHUNK);
}



void free_chunk (struct chunk *chk)
{
	free_object (CHUNK, (short *) chk);
}

