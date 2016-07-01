/*
	Molecular Surface Package

	Copyright 1986, 1989 by Michael L. Connolly
	All rights reserved

	January 5, 2002
*/

/* Quartets of atoms that a probe is simultaneously tangent to */


/* Later: consider altitude when 4 atom centers not coplanar */
/* Later: consider probes that are just barely eliminated by collisions */

/* For each list of pholders,
   compare every pair of entries, 
   and create a quartet if they
   (a) share two atom numbers, and
   (b) are close (less than epsilon)
   We are left with a linked list of quartets.
   Then check quartets for redundancy.
   If two quartets match
   (match means same four atoms and close),
   then one is marked deleted,
   and its probe information is transferred
   to the other quartet.
   After all pairs of quartets have been compared,
   remove flagged-as-deleted quartets from list.
   Finally, print linked list of quartets.
   The quartet code should be called after the rolling,
   but before concave face generation.
 */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

void find_quartets (struct surface *srf)
{
	char message[MAXLINE];
	distribute_probes (srf);
	reduce_quartets (srf);
	cleanup_quartets (srf);
	all_quartet_tori (srf);
	count_quartets (srf);
	add_quartet_probes (srf);
	subtract_quartet_probes (srf);
	renumber_probes (srf);
	free_quartets (srf);
	sprintf (message, "%8ld probe positions", srf -> n_probe);
	inform(message);
}

struct quartet *allocate_quartet ()
{
    struct quartet *quart;

	quart = (struct quartet *) allocate_object (QUARTET);
	if (quart == NULL) {
		return(NULL);
	}
	return (quart);
}

struct pholder *allocate_pholder ()
{
    struct pholder *ph;
	char message[MAXLINE];

	ph = (struct pholder *) allocate_object (PHOLDER);
	if (ph == NULL) {
		return(NULL);
	}
	return (ph);
}

int compare_probes (struct probe *prb1, struct probe *prb2)
{
	double epsilon;
	int nsame;
    int ok;
	double dist;

	nsame = 0;
	epsilon = EPSILON * EPSILON;
	nsame = same_atoms(prb1 -> atom, prb2 -> atom);
	dist = distance (prb1 -> center, prb2 -> center);
	ok = (nsame >= 2 && dist < epsilon);
	return (ok);
}

int same_atoms (atomnum a[MAXPA], atomnum b[MAXPA])
{
	int j, k, nsame;

	nsame = 0;
	for (j = 0; j < MAXPA; j++) {
		if (a[j] <= 0) continue;
		for (k = 0; k < MAXPA; k++) {
			if (b[k] <= 0) continue;
			if (a[j] == b[k]) nsame++;
		}
	}
	return (nsame);
}

void distribute_probes (struct surface *srf) {
	int anum, idx, j;
	long nlist;
	struct pholder **heads;
	struct pholder **handle;
	struct probe *prb;
	char message[MAXLINE];

	nlist = srf -> n_atom;
	heads = (struct pholder **)
		allocate_pointers (PHOLDER, nlist);
	/* add probes to array of linked lists */
	for (prb = srf -> head_probe; prb != NULL; prb = prb -> next) {
		sprintf (message, "%8hd %8hd %8hd %8.3f %8.3f %8.3f probe", 
			prb -> atom[0], prb -> atom[1], prb -> atom[2],
			prb -> center[0], prb -> center[1], prb -> center[2]);
		informd2(message);
		for (j = 0; j < MAXPA; j++) {
			anum = prb -> atom[j];
			if (anum <= 0) continue;
			idx = anum - 1;
			handle = heads + idx;
			*handle = add_pholder (*handle, prb);
		}
	}
	/* create quartets */
	for (anum = 1; anum <= nlist; anum++) {
		idx = anum - 1;
		handle = heads + idx;
		if (*handle == NULL) continue;
		some_pholders (srf, *handle);
		if (error()) return;
		free_pholders (*handle);
		if (error()) return;
	}
	free_pointers (PHOLDER, heads);
}

struct pholder *add_pholder (struct pholder *head, struct probe *prb)
{
	struct pholder *ph;

	ph =  allocate_pholder ();
	if (error()) return (NULL);
	ph -> prb = prb;
	ph -> next = head;
	return (ph);
}

void free_pholders (struct pholder *head)
{
	struct pholder *ph, *next;
	char message[MAXLINE];

	for (ph = head; ph != NULL; ph = next) {
		next = ph -> next;
		free_object (PHOLDER, (short *) ph);
		if (error()) return;
	}
	free_cache (PHOLDER);
}

void some_pholders (struct surface *srf, struct pholder *head) {
	int result;
	struct probe *prb1, *prb2;
	struct pholder *ph1, *ph2;

	for (ph1 = head; ph1 != NULL; ph1 = ph1 -> next) {
		prb1 = ph1 -> prb;
		for (ph2 = ph1 -> next; ph2 != NULL; ph2 = ph2 -> next) {
			prb2 = ph2 -> prb;
			result = compare_probes (prb1, prb2);
			if (result) {
				make_quartet (srf, prb1, prb2);
				if (error()) return;
			}
		}
	}
}

struct quartet *make_quartet (struct surface *srf, struct probe *prb1, struct probe *prb2)
{
	int j, k, not_there_already;
	int reverse1, reverse2;
	int nduplicate;
	int position1[MAXPA];
	int position2[MAXPA];
	int direction1, direction2;
	double volume1, volume2;
	double center[3];
	struct quartet *quart;
	atomnum swap_number;
	struct sphere *swap_pointer;
	atomnum numbers[MAXPA];
	atomnum sorted_numbers[MAXPA];
	struct sphere *pointers[MAXPA];
	struct sphere *sorted_pointers[MAXPA];
	int duplicate[MAXPA];
    struct cept *ex;
	char message[MAXLINE];

	nduplicate = 0;
	for (j = 0; j < MAXPA; j++) {
		duplicate[j] = 0;
		position1[j] = -1;
		position2[j] = -1;
	}
	for (k = 0; k < 3; k++) {
		center[k] = (prb1 -> center[k] + prb2 -> center[k]) / 2.0;
	}
	for (j = 0; j < 3; j++) {
		numbers[j] = prb1 -> atom[j];
		pointers[j] = prb1 -> atm[j];
	}
	for (j = 0; j < 3; j++) {
		not_there_already = 1;
		for (k = 0; k < 3; k++) {
			if (prb2 -> atom[j] == prb1 -> atom[k]) {
				not_there_already = 0;
				duplicate[k] = 1;
				position1[nduplicate] = k;
				position2[nduplicate] = j;
				nduplicate++;
			}
		}
		if (not_there_already) {
			numbers[3] = prb2 -> atom[j];
			pointers[3] = prb2 -> atm[j];
		}
	}
	if (nduplicate != 2) {
		ex = new_cept (GEOMETRY_ERROR,  INCONSISTENCY,  FATAL_SEVERITY);
		add_function (ex, "make_quartet");
		add_long (ex, "nduplicate", nduplicate);
		add_long (ex, "prb1 -> number", prb1 -> number);
		add_long (ex, "prb2 -> number", prb2 -> number);
		add_long (ex, "prb1 -> atom[0]", prb1 -> atom[0]);
		add_long (ex, "prb1 -> atom[1]", prb1 -> atom[1]);
		add_long (ex, "prb1 -> atom[2]", prb1 -> atom[2]);
		add_long (ex, "prb2 -> atom[0]", prb2 -> atom[0]);
		add_long (ex, "prb2 -> atom[1]", prb2 -> atom[1]);
		add_long (ex, "prb2 -> atom[2]", prb2 -> atom[2]);
		add_source (ex, "msquartet.c");
		return(NULL);
	}
	direction1 = position1[1] - position1[0];
	if (direction1 == 2) direction1 -= 3;
	if (direction1 == -2) direction1 += 3;
	direction2 = position2[1] - position2[0];
	if (direction2 == 2) direction2 -= 3;
	if (direction2 == -2) direction2 += 3;
	/* both directions are either +1 or -1 */
	/* if they disagree, the other atoms
       are on opposite sides of the duplicate atoms,
	   but if they agree, the other atoms are on
      the same side, and must be sorted */
	if (direction1 == 1 && direction2 == -1) {
		switch (position1[0]) {
			case 0:
				sorted_numbers[0] = numbers[0];
				sorted_numbers[1] = numbers[3];
				sorted_numbers[2] = numbers[1];
				sorted_numbers[3] = numbers[2];
				sorted_pointers[0] = pointers[0];
				sorted_pointers[1] = pointers[3];
				sorted_pointers[2] = pointers[1];
				sorted_pointers[3] = pointers[2];
				break;
			case 1:
				sorted_numbers[0] = numbers[0];
				sorted_numbers[1] = numbers[1];
				sorted_numbers[2] = numbers[3];
				sorted_numbers[3] = numbers[2];
				sorted_pointers[0] = pointers[0];
				sorted_pointers[1] = pointers[1];
				sorted_pointers[2] = pointers[3];
				sorted_pointers[3] = pointers[2];
				break;
			case 2:
				sorted_numbers[0] = numbers[0];
				sorted_numbers[1] = numbers[1];
				sorted_numbers[2] = numbers[2];
				sorted_numbers[3] = numbers[3];
				sorted_pointers[0] = pointers[0];
				sorted_pointers[1] = pointers[1];
				sorted_pointers[2] = pointers[2];
				sorted_pointers[3] = pointers[3];
				break;
			default:
				break;
		}
	}
	else if (direction1 == -1 && direction2 == 1) {
		switch (position1[0]) {
			case 0:
				sorted_numbers[0] = numbers[0];
				sorted_numbers[1] = numbers[1];
				sorted_numbers[2] = numbers[2];
				sorted_numbers[3] = numbers[3];
				sorted_pointers[0] = pointers[0];
				sorted_pointers[1] = pointers[1];
				sorted_pointers[2] = pointers[2];
				sorted_pointers[3] = pointers[3];
				break;
			case 1:
				sorted_numbers[0] = numbers[0];
				sorted_numbers[1] = numbers[3];
				sorted_numbers[2] = numbers[1];
				sorted_numbers[3] = numbers[2];
				sorted_pointers[0] = pointers[0];
				sorted_pointers[1] = pointers[3];
				sorted_pointers[2] = pointers[1];
				sorted_pointers[3] = pointers[2];
				break;
			case 2:
				sorted_numbers[0] = numbers[0];
				sorted_numbers[1] = numbers[1];
				sorted_numbers[2] = numbers[3];
				sorted_numbers[3] = numbers[2];
				sorted_pointers[0] = pointers[0];
				sorted_pointers[1] = pointers[1];
				sorted_pointers[2] = pointers[3];
				sorted_pointers[3] = pointers[2];
				break;
			default:
				break;
		}
	}
	else {
		/* sort */
		if (duplicate[0] && duplicate[1]) {
			sorted_numbers[0] = numbers[1];
			sorted_numbers[1] = numbers[2];
			sorted_numbers[2] = numbers[3];
			sorted_numbers[3] = numbers[0];
			sorted_pointers[0] = pointers[1];
			sorted_pointers[1] = pointers[2];
			sorted_pointers[2] = pointers[3];
			sorted_pointers[3] = pointers[0];
		}
		else if (duplicate[0] && duplicate[2]) {
			sorted_numbers[0] = numbers[0];
			sorted_numbers[1] = numbers[1];
			sorted_numbers[2] = numbers[3];
			sorted_numbers[3] = numbers[2];
			sorted_pointers[0] = pointers[0];
			sorted_pointers[1] = pointers[1];
			sorted_pointers[2] = pointers[3];
			sorted_pointers[3] = pointers[2];
		}
		else if (duplicate[1] && duplicate[2]) {
			sorted_numbers[0] = numbers[2];
			sorted_numbers[1] = numbers[3];
			sorted_numbers[2] = numbers[0];
			sorted_numbers[3] = numbers[1];
			sorted_pointers[0] = pointers[2];
			sorted_pointers[1] = pointers[3];
			sorted_pointers[2] = pointers[0];
			sorted_pointers[3] = pointers[1];
		}
		else {
			ex = new_cept (GEOMETRY_ERROR,  INCONSISTENCY,  FATAL_SEVERITY);
			add_function (ex, "make_quartet");
			add_source (ex, "msquartet.c");
			add_message (ex, "unexpected sort case");
			return(NULL);
		}
		volume1 = tetrahedron_volume (center, sorted_pointers[0] -> center,
			sorted_pointers[1] -> center, sorted_pointers[2] -> center);
		reverse1 = (volume1 < 0.0);
		volume2 = tetrahedron_volume (center, sorted_pointers[1] -> center,
			sorted_pointers[2] -> center, sorted_pointers[3] -> center);
		reverse2 = (volume2 < 0.0);
		if (reverse1 && reverse2) {
			swap_number = sorted_numbers[1];
			sorted_numbers[1] = sorted_numbers[2];
			sorted_numbers[2] = swap_number;
			swap_pointer = sorted_pointers[1];
			sorted_pointers[1] = sorted_pointers[2];
			sorted_pointers[2] = swap_pointer;
		}
	}
	sprintf (message, "%5d %5d %5d : %5d %5d %5d",
		prb1 -> atom[0], prb1 -> atom[1], prb1 -> atom[2],
		prb2 -> atom[0], prb2 -> atom[1], prb2 -> atom[2]);
	informd (message);
	sprintf (message, "unsorted numbers = %5d %5d %5d %5d",
		numbers[0], numbers[1], numbers[2], numbers[3]);
	informd (message);
	sprintf (message, "  sorted numbers = %5d %5d %5d %5d",
		sorted_numbers[0], sorted_numbers[1], sorted_numbers[2], sorted_numbers[3]);
	informd (message);
	quart = allocate_quartet();
	if (error()) return (NULL);
	for (k = 0; k < 3; k++) {
		quart -> center[k] = center[k];
	}
	quart -> prb[0] = prb1;
	quart -> prb[1] = prb2;
	quart -> probe_number[0] = prb1 -> number;
	quart -> probe_number[1] = prb2 -> number;
	for (j = 0; j < MAXPA; j++) {
		quart -> atm[j] = sorted_numbers[j];
		quart -> atmptr[j] = sorted_pointers[j];
	}
	link_quartet (srf, quart);
	return (quart);
}

void link_quartet (struct surface *srf, struct quartet *quart) {
	if (srf -> head_quartet == NULL) {
		srf -> head_quartet = quart;
		srf -> tail_quartet = quart;
		quart -> next = NULL;
		quart -> previous = NULL;
		return;
	}
	srf -> tail_quartet -> next = quart;
	quart -> previous = srf -> tail_quartet;
	srf -> tail_quartet = quart;
}

void reduce_quartets (struct surface *srf)
{
	long original;
	struct quartet *quart1, *quart2;
	char message[MAXLINE];

	original = 0;
	for (quart1 = srf -> head_quartet; quart1 != NULL; quart1 = quart1 -> next)
		original++;
	sprintf (message, "%8ld original quartets", original);
	informd (message);
	for (quart1 = srf -> head_quartet; quart1 != NULL; quart1 = quart1 -> next) {
		if (quart1 -> delete) continue;
		for (quart2 = quart1 -> next; quart2 != NULL; quart2 = quart2 -> next) {
			if (quart2 -> delete) continue;
			if (quartets_match (quart1, quart2)) {
				merge_quartets (srf, quart1, quart2);
			}
		}
	}
}

struct quartet *merge_quartets (struct surface *srf, struct quartet *quart1, struct quartet *quart2)
{
	int j, k, not_there_already;

	for (j = 0; j < MAXPQ; j++) {
		if (quart2 -> prb[j] == NULL) continue;
		not_there_already = 1;
		for (k = 0; k < MAXPQ; k++) {
			if (quart1 -> prb[k] == NULL) continue;
			if (quart1 -> prb[k] == quart2 -> prb[j]) {
				not_there_already = 0;
				break;
			}
		}
		if (!not_there_already) continue;
		/* quart2 -> prb[j] is not in quart1 */
		/* look for empty slot */
		for (k = 0; k < MAXPQ; k++) {
			if (quart1 -> prb[k] == NULL) {
				quart1 -> prb[k] = quart2 -> prb[j];
				quart1 -> probe_number[k] = quart2 -> prb[j] -> number;
				break;
			}
		}
	}
	quart2 -> delete = 1;
	return (quart1);
}

void count_quartets (struct surface *srf) {
	long n_quartet;
	char message[MAXLINE];
	struct quartet *quart;

	n_quartet = 0;
	for (quart = srf -> head_quartet; quart != NULL; quart = quart -> next) {
		print_quartet (quart);
		n_quartet++;
	}
	sprintf (message, "%8ld probes tangent to four atoms", n_quartet);
	if (n_quartet > 0)
		inform(message);
}

void free_quartets (struct surface *srf) {
	struct quartet *quart, *next;

	next = NULL;
	for (quart = srf -> head_quartet; quart != NULL; quart = next) {
		next = quart -> next;
		free_object (QUARTET, (short *) quart);
	}
	free_cache (QUARTET);
}

void print_quartet (struct quartet *quart)
{
	char message[MAXLINE];
	sprintf (message, "%8d %8d %8d %8d quartet %8ld %8ld %8ld %8ld",
		quart -> atm[0], quart -> atm[1],
		quart -> atm[2], quart -> atm[3],
		quart -> probe_number[0],
		quart -> probe_number[1],
		quart -> probe_number[2],
		quart -> probe_number[3]);
	informd (message);
}

void cleanup_quartets (struct surface *srf)
{
	int ndeleted;
	struct quartet *quart;
	ndeleted = 1;
	while (ndeleted > 0) {
		ndeleted = 0;
		for (quart = srf -> head_quartet; quart != NULL; quart = quart -> next) {
			if (quart -> delete) {
				delete_quartet (srf, quart);
				ndeleted++;
				break;
			}
		}
	}
}

void delete_quartet (struct surface *srf, struct quartet *quart)
{
	if (quart -> previous == NULL)
		srf -> head_quartet = quart -> next;
	else quart -> previous -> next = quart -> next;
	if (quart -> next == NULL)
		srf -> tail_quartet = quart -> previous;
	else quart -> next -> previous = quart -> previous;
	free_object (QUARTET, (short *) quart);
}

int quartets_match (struct quartet *quart1, struct quartet *quart2)
{
	int j, k, nmatch;
	double dist;

	dist = distance (quart1 -> center, quart2 -> center);
	if (dist > EPSILON) return (0);
	nmatch = 0;
	for (j = 0; j < MAXPA; j++)
		for (k = 0; k < MAXPA; k++) 
			if (quart1 -> atm[j] == quart2 -> atm[k]) nmatch++;
	return (nmatch >= MAXPA);
}

void all_quartet_tori (struct surface *srf)
{
	struct quartet *quart;

	for (quart = srf -> head_quartet; quart != NULL; quart = quart -> next) {
		quartet_tori (srf, quart);
	}
}

void quartet_tori (struct surface *srf, struct quartet *quart)
{
	int j, j1, j2, nt, nbdt;
	struct sphere *sphere_ptr;
	struct neighbor *nbr, *first_neighbor, *last_neighbor;
	struct pair *torptr;
	char message[MAXLINE];

	nt = 0;
	nbdt = 0;
	for (sphere_ptr = (struct sphere *) (srf -> head_atom); sphere_ptr != NULL;
		sphere_ptr = sphere_ptr -> next) {
		first_neighbor = sphere_ptr -> first_neighbor;
		last_neighbor = sphere_ptr -> last_neighbor;
		for (nbr = first_neighbor; nbr <= last_neighbor; nbr++) {
			torptr = nbr -> torptr;
			if (torptr == NULL) continue;
			/* mark diagonal tori as buried */
			for (j = 0; j < MAXPA; j++) {
				if (torptr -> buried) continue;
				j2 = j + 2;
				if (j2 >= MAXPA) j2 -= MAXPA;
				if (quart -> atmptr[j] == torptr -> sph[0] &&
						quart -> atmptr[j2] == torptr -> sph[1]) {
					torptr -> buried = 1;
					nbdt++;
				}
				else if (quart -> atmptr[j] == torptr -> sph[1] &&
						quart -> atmptr[j2] == torptr -> sph[0]) {
					torptr -> buried = 1;
					nbdt++;
				}
			}
			/* check only once, not twice */
			if (sphere_ptr != torptr -> sph[0]) continue;
			for (j = 0; j < MAXPA; j++) {
				j1 = j + 1;
				if (j1 >= MAXPA) j1 = 0;
				if (quart -> atmptr[j] == torptr -> sph[0] &&
						quart -> atmptr[j1] == torptr -> sph[1]) {
					quart -> torptr[j] = torptr;
					nt++;
				}
				else if (quart -> atmptr[j] == torptr -> sph[1] &&
						quart -> atmptr[j1] == torptr -> sph[0]) {
					quart -> torptr[j] = torptr;
					nt++;
				}
			}
		}
	}
	sprintf (message, "%8d tori (and %3d buried diagonal tori) found for quartet", nt, nbdt);
	informd (message);
}

struct probe *quartet_to_probe (struct surface *srf, struct quartet *quart) 
{
	int k;
	short natom;
	struct probe *prb, *prb1, *prb2, *prb3, *prb4;
    struct cept *ex;

	/* allocate memory */
	prb = (struct probe *) allocate_object (PROBE);
	if (prb == NULL) {
		ex = new_cept (MEMORY_ERROR,  ALLOCATION,  FATAL_SEVERITY);
		add_function (ex, "quartet_to_probe");
		add_source (ex, "msquartet.c");
		return(NULL);
	}

	prb1 = quart -> prb[0];
	prb2 = quart -> prb[1];
	prb3 = quart -> prb[2];
	prb4 = quart -> prb[3];
	/* set up fields */
	for (k = 0; k < 3; k++)
		prb -> center[k] = quart -> center[k];

	prb -> height = prb1 -> height;
	/* store whether a low probe */
	prb -> low = (prb -> height < srf -> probe_radius);
	prb -> unit_altitude_z = prb1 -> unit_altitude_z;
	natom = 0;
	for (k = 0; k < 4; k++) {
		if (quart -> atm[k] != 0) natom++;
		prb -> atom[k] = quart -> atm[k];
		prb -> atm[k] = quart -> atmptr[k];
		prb -> pairs[k] = quart -> torptr[k];
	}
	prb -> natom = natom;
	return (prb);
}

void subtract_quartet_probes (struct surface *srf)
{
	int j;
	struct probe *prb;
	struct quartet *quart;

	for (quart = srf -> head_quartet; quart != NULL; quart = quart -> next) {
		for (j = 0; j < MAXPQ; j++) {
			prb = quart -> prb[j];
			if (prb == NULL) continue;
			delink_probe (srf, prb);
			/* we cannot free probe, because others point at it */
		}
	}
}

void add_quartet_probes (struct surface *srf)
{
	struct quartet *quart;
	struct probe *prb;

	for (quart = srf -> head_quartet; quart != NULL; quart = quart -> next) {
		prb = quartet_to_probe (srf, quart);
		link_probe (srf, prb);
	}
}

void renumber_probes (struct surface *srf)
{
	long n_probe;
	struct probe *prb;

	n_probe = 0;
	for (prb = srf -> head_probe; prb != NULL; prb = prb -> next) {
		prb -> number = ++n_probe;
	}
}

