/* Molecular Surface Package Copyright 1995 Michael L. Connolly */
/* February 3, 2000 */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"


long collect_edges (struct surface *phn)
{
	long p, first, second, n_polygon, edge_number, n_convert, u;
	int n_side, orn, i, j;
	struct polygon *poly, **polygon_handles;
	struct phnedg *edg, **phnedg_handles;

	init_edger (phn);
	n_polygon = phn -> n_polygon;
	if (n_polygon == 0) return (0);
	polygon_handles = phn -> polygon_handles;
	if (polygon_handles == NULL) return (0);
	for (p = 0; p < n_polygon; p++) {
		poly = *(polygon_handles + p);
		if (poly == NULL) return (0L);
		n_side = poly -> n_side;
		for (i = 0; i < n_side; i++) {
			j = (i < n_side - 1) ? i + 1 : 0;
			first = poly -> vertex_index[i] + 1;
			second = poly -> vertex_index[j] + 1;
			edge_number = lookup_edger (phn, first, second);
			if (error()) return (0L);
			if (edge_number == 0) return (0L);
			poly -> edge_number[i] = edge_number;
		}
	}
	n_convert = convert_edger (phn);
	if (error ())
		return (0L);
	if (n_convert == 0)
		return (0L);
	phnedg_handles = phn -> phnedg_handles;
	if (phnedg_handles == NULL)
		return (0L);
	for (p = 0; p < n_polygon; p++) {
		poly = *(polygon_handles + p);
		if (poly == NULL)
			return (0L);
		n_side = poly -> n_side;
		if (n_side == 0)
			return (0L);
		for (i = 0; i < n_side; i++) {
			edge_number = poly -> edge_number[i];
			if (edge_number == 0)
				return (0L);
			u = abs (edge_number);
			orn = (edge_number > 0) ? 0 : 1;
			edg = *(phnedg_handles + u - 1);
			if (edg == NULL)
				return (0L);
			poly -> edg[i] = edg;
			poly -> orn[i] = (short) orn;
		}
	}
	return (n_convert);
}



int init_edger (struct surface *phn)
{
	phn -> n_edger = 0;
	phn -> root_edger = NULL;
	return (1);
}

struct edger *new_edger (struct surface *phn, long first, long second)
{
	static long counter = 0L;
	char message[MAXLINE];
	struct edger *e;
	
	e = (struct edger *) allocate_object (EDGER);
	if (e == NULL) {
		sprintf(message, "(new_edger) fails at edger %ld", counter);
		set_error1 (message);
		return(NULL);
	}
	counter++;
	e -> first = first;
	e -> second = second;
	e -> sum = first + second;
	e -> number = ++(phn -> n_edger);
	return (e);
}

struct edger *insert_edger (struct surface *phn, struct edger *e, long first, long second)
{
	long sum;
	struct edger *n;
	
	sum = first + second;
	if (e == NULL) {
		n = new_edger (phn, first, second);
		if (n == NULL) return (NULL);
		return (n);
	}
	if (sum < e -> sum) {
		n = insert_edger (phn, e -> left, first, second);
		if (n == NULL) return (NULL);
		if (e -> left == NULL) e -> left = n;
	}
	else if (sum == e -> sum) {
		if (e -> first == first && e -> second == second)
			return (e);
		else if (e -> first == second && e -> second == first)
			return (e);
		n = insert_edger (phn, e -> middle, first, second);
		if (n == NULL) return (NULL);
		if (e -> middle == NULL) e -> middle = n;
	}
	else {
		n = insert_edger (phn, e -> right, first, second);
		if (n == NULL) return (NULL);
		if (e -> right == NULL) e -> right = n;
	}
	return (n);
}

long lookup_edger (struct surface *phn, long first, long second)
{
	struct edger *f;

	f = insert_edger (phn, phn -> root_edger, first, second);
	if (error()) return (0L);
	if (f == NULL) {
		set_error1("lookup_edger: null pointer from insert_edger");
		return (0L);
	}
	if (f -> first == first && f -> second == second)
		return (f -> number);
	else if (f -> first == second && f -> second == first)
		return (-(f -> number));
	else return (0L);
}

long convert_edger (struct surface *phn)
{
	int result;
	long ne, ie;
	char message[MAXLINE];
	struct phnedg *pe;
	struct phnedg **phnedg_handles;

	ne = phn -> n_edger;
	if (ne <= 0)
		return (0);
	phnedg_handles = (struct phnedg **)
		allocate_pointers (PHNEDG, ne);
	if (phnedg_handles == NULL)
		return (0);
	for (ie = 0; ie < ne; ie++) {
		pe = allocate_phnedg ();
		if (pe == NULL) {
			sprintf (message, "(convert_edger) edge %ld of %ld", ie+1, ne);
			set_error2 (message);
			return (0);
		}
		*(phnedg_handles + ie) = pe;
	}
	result = recurse_edger (phn, phn -> root_edger);
	if (!result)
		return (0);
	phn -> phnedg_handles = phnedg_handles;
	return (ne);
}
	
int recurse_edger (struct surface *phn, struct edger *e)
{
	long number;
	struct phnedg *pe;

	if (e == NULL) return (1);
	number = e -> number;
	if (number <= 0) return (0);
	pe = *(phn -> phnedg_handles + number - 1);
	if (pe == NULL) return (0);
	pe -> vtxnum[0] = e -> first;
	pe -> vtxnum[1] = e -> second;
	pe -> pvt[0] = num2phnvtx (phn, pe -> vtxnum[0]);
	pe -> pvt[1] = num2phnvtx (phn, pe -> vtxnum[1]);
	if (!recurse_edger (phn, e -> left)) return (0);
	if (!recurse_edger (phn, e -> middle)) return (0);
	if (!recurse_edger (phn, e -> right)) return (0);
	return (1);
}

