/*
 * Molecular Surface Package
 * Copyright 1986 by Michael L. Connolly
 * All Rights Reserved
 * March 7, 2000
 */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"



struct array *new_array (long type, long length)
{
	char message[MAXLINE];
	struct array *ary;
	
	if (length <= 0) {
		set_error1 ("new_array: null length");
		return (NULL);
	}
	ary = (struct array *) allocate_object (ARRAY);
	if (ary == NULL) {
		set_error1 ("new_array: memory failure");
		return (NULL);
	}
	ary -> type = type;
	ary -> length = length;
	switch (type) {
	case INTEGER:
		ary -> integers = allocate_longs (length, 0, INTEGERS);
		if (ary -> integers == NULL) {
			set_error1 ("new_array: memory failure");
			return (NULL);
		}
		break;
	case REAL:
		ary -> reals = allocate_doubles (length, 0, REALS);
		if (ary -> reals == NULL) {
			set_error1 ("new_array: memory failure");
			return (NULL);
		}
		break;
	case SPHERE:
		ary -> radii = allocate_doubles(length, 0, RADII);
		if (ary -> radii == NULL) {
			set_error1 ("new_array: memory failure");
			return (NULL);
		}
		ary -> centers = allocate_doubles (3 * length, 0, CENTERS);
		if (ary -> centers == NULL) {
			set_error1 ("new_array: memory failure");
			return (NULL);
		}
		break;
	case STRING:
		ary -> strings = (char **) allocate_char_pointers (length);
		if (ary -> strings == NULL) {
			set_error1 ("new_array: memory failure");
			return (NULL);
		}
		break;
	case ATOM:
		ary -> atoms = (struct atom **) allocate_pointers (ATOM, length);
		if (ary -> atoms == NULL) {
			set_error1 ("new_array: memory failure");
			return (NULL);
		}
		break;
	case BOND:
		ary -> bonds = (struct bond **) allocate_pointers (BOND, length);
		if (ary -> bonds == NULL) {
			set_error1 ("new_array: memory failure");
			return (NULL);
		}
		break;
	default:
		set_error1 ("new_array: invalid type");
		return (NULL);
	}
	return (ary);
}


int free_array (struct array *ary)
{
	long l;
	long length;
	char *str;
	char message[MAXLINE];
	
	length = ary -> length;
	if (ary -> integers != NULL) free_longs (ary -> integers, 0, INTEGERS);
	if (ary -> reals != NULL) free_doubles (ary -> reals, 0, REALS);
	if (ary -> radii != NULL) free_doubles (ary -> radii, 0, RADII);
	if (ary -> centers != NULL) free_doubles (ary -> centers, 0, CENTERS);
	if (ary -> strings != NULL) {
		/* in this case, also free the strings themselves */
		for (l = 0; l < length; l++) {
			str = *(ary -> strings + l);
			free_chars (str);
		}
		free_char_pointers (ary -> strings);
	}
	if (ary -> atoms != NULL) free_pointers (ATOM, ary -> atoms);
	if (ary -> bonds != NULL) free_pointers (BOND, ary -> bonds);
	free_object (ARRAY, (short *) ary);
	return (1);
}

int store_integer (struct array *ary, long idx, long value)
{
	if (ary == NULL) {
		set_error1 ("store_integer: null array");
		return (0);
	}
	if (ary -> integers == NULL) {
		set_error1 ("store_integer: null integer array");
		return (0);
	}
	if (idx < 0 || idx >= ary -> length) {
		set_error1 ("store_integer: invalid index");
		return (0);
	}
	*(ary -> integers + idx) = value;
	return (1);
}

long fetch_integer (struct array *ary, long idx)
{
	long value;

	if (ary == NULL) {
		set_error1 ("fetch_integer: null array");
		return (0);
	}
	if (ary -> integers == NULL) {
		set_error1 ("fetch_integer: null integer array");
		return (0);
	}
	if (idx < 0 || idx >= ary -> length) {
		set_error1 ("fetch_integer: invalid index");
		return (0);
	}
	value = *(ary -> integers + idx);
	return (value);
}

int store_real (struct array *ary, long idx, double value)
{
	if (ary == NULL) {
		set_error1 ("store_real: null array");
		return (0);
	}
	if (ary -> reals == NULL) {
		set_error1 ("store_real: null real array");
		return (0);
	}
	if (idx < 0 || idx >= ary -> length) {
		set_error1 ("store_real: invalid index");
		return (0);
	}
	*(ary -> reals + idx) = value;
	return (1);
}

double fetch_real (struct array *ary, long idx)
{
	double value;

	if (ary == NULL) {
		set_error1 ("fetch_real: null array");
		return (0.0);
	}
	if (ary -> reals == NULL) {
		set_error1 ("fetch_real: null real array");
		return (0.0);
	}
	if (idx < 0 || idx >= ary -> length) {
		set_error1 ("fetch_real: invalid index");
		return (0.0);
	}
	value = *(ary -> reals + idx);
	return (value);
}

int store_sphere (struct array *ary, long idx, double center[3], double radius)
{
	int k;

	if (ary == NULL) {
		set_error1 ("store_sphere: null array");
		return (0);
	}
	if (idx < 0 || idx >= ary -> length) {
		set_error1 ("store_sphere: invalid index");
		return (0);
	}
	if (ary -> radii == NULL || ary -> centers == NULL) {
		set_error1 ("store_sphere: null radii array");
		return (0);
	}
	*(ary -> radii + idx) = radius;
	for (k = 0; k < 3; k++)
		*(ary -> centers + 3 * idx + k) = center[k];
	return (1);
}

double fetch_sphere (struct array *ary, long idx, double center[3])
{
	int k;
	double radius;

	if (ary == NULL) {
		set_error1 ("fetch_sphere: null array");
		return (0.0);
	}
	if (ary -> radii == NULL || ary -> centers == NULL) {
		set_error1 ("fetch_sphere: null radius array");
		return (0.0);
	}
	if (idx < 0 || idx >= ary -> length) {
		set_error1 ("fetch_real: invalid index");
		return (0.0);
	}
	for (k = 0; k < 3; k++)
		center[k] = *(ary -> centers + 3 * idx + k);
	radius = *(ary -> radii + idx);
	return (radius);
}

int store_atom (struct array *ary, long idx, struct atom *ptr)
{
	if (ary == NULL) {
		set_error1 ("store_atom: null array");
		return (0);
	}
	if (ary -> atoms == NULL) {
		set_error1 ("store_atom: null atom array");
		return (0);
	}
	if (idx < 0 || idx >= ary -> length) {
		set_error1 ("store_atom: invalid index");
		return (0);
	}
	*(ary -> atoms + idx) = ptr;
	return (1);
}

struct atom *fetch_atom (struct array *ary, long idx)
{
	struct atom *ptr;

	if (ary == NULL) {
		set_error1 ("fetch_atom: null array");
		return (0);
	}
	if (ary -> atoms == NULL) {
		set_error1 ("fetch_atom: null atom array");
		return (0);
	}
	if (idx < 0 || idx >= ary -> length) {
		set_error1 ("fetch_atom: invalid index");
		return (0);
	}
	ptr = *(ary -> atoms + idx);
	return (ptr);
}

int store_bond (struct array *ary, long idx, struct bond *ptr)
{
	if (ary == NULL) {
		set_error1 ("store_bond: null array");
		return (0);
	}
	if (ary -> bonds == NULL) {
		set_error1 ("store_bond: null bond array");
		return (0);
	}
	if (idx < 0 || idx >= ary -> length) {
		set_error1 ("store_bond: invalid index");
		return (0);
	}
	*(ary -> bonds + idx) = ptr;
	return (1);
}

struct bond *fetch_bond (struct array *ary, long idx)
{
	struct bond *ptr;

	if (ary == NULL) {
		set_error1 ("fetch_bond: null array");
		return (0);
	}
	if (ary -> bonds == NULL) {
		set_error1 ("fetch_bond: null bond array");
		return (0);
	}
	if (idx < 0 || idx >= ary -> length) {
		set_error1 ("fetch_bond: invalid index");
		return (0);
	}
	ptr = *(ary -> bonds + idx);
	return (ptr);
}

int store_string (struct array *ary, long idx, char *ptr)
{
	long len;
	char *str;

	if (ary == NULL) {
		set_error1 ("store_string: null array");
		return (0);
	}
	if (ary -> strings == NULL) {
		set_error1 ("store_string: null string array");
		return (0);
	}
	if (idx < 0 || idx >= ary -> length) {
		set_error1 ("store_string: invalid index");
		return (0);
	}
	len = strlen (ptr);
	str = allocate_chars (len+1);
	if (str == NULL) {
		set_error1 ("store_string: memory failure");
		return (0);
	}
	strcpy (str, ptr);
	*(ary -> strings + idx) = str;
	return (1);
}

char *fetch_string (struct array *ary, long idx)
{
	char *ptr;

	if (ary == NULL) {
		set_error1 ("fetch_string: null array");
		return (0);
	}
	if (ary -> strings == NULL) {
		set_error1 ("fetch_string: null string array");
		return (0);
	}
	if (idx < 0 || idx >= ary -> length) {
		set_error1 ("fetch_string: invalid index");
		return (0);
	}
	ptr = *(ary -> strings + idx);
	return (ptr);
}

struct array *read_lines (FILE *fp)
{
	long n_line, i_line, len;
	char line[2 * MAXLINE];
	char message[MAXLINE];
	struct array *lines;

	/* count atoms */
	n_line = 0;
	for (;;) {
		fgets (line, MAXLINE, fp);
		if (feof (fp)) break;
		len = strlen (line);
		if (len >= MAXLINE) {
			sprintf (message, "read_lines: %ld characters in line", len);
			set_error1 (message);
			return (NULL);
		}
		n_line++;
	}
	if (n_line == 0) {
		set_error1 ("read_lines: empty file");
		return(NULL);
	}
	/* rewind file using fseek because rewind() is not everywhere */
	fseek (fp, 0L, 0);

	lines = new_array ((long) STRING, n_line);
	if (error ()) return (NULL);

	for (i_line = 0; i_line < n_line; i_line++) {
		fgets (line, MAXLINE, fp);
		if (feof (fp)) break;
		len = strlen (line);
		if (len >= MAXLINE) {
			sprintf (message, "read_lines: %ld characters in line", len);
			set_error1 (message);
			return (NULL);
		}
		store_string (lines, i_line, line);
		if (error ()) return (NULL);
	}
	return (lines);
}


/*
 * Molecular Surface Package
 * Copyright 1986 by Michael L. Connolly
 * All Rights Reserved
 * June 6, 1998
 */

