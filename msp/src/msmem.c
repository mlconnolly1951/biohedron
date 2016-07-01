/*
	Molecular Surface Package
	Copyright 1986, 1989 by Michael L. Connolly
	All rights reserved
	November 19, 2001
*/

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"



/* memory allocation routines */

int init_mem (int mt)
{
	long shortsize, charsize;
	unsigned long size;
	
	shortsize = sizeof (short);
	charsize = sizeof (char);
	bytes_per_short = shortsize/charsize;
	if (shortsize != bytes_per_short*charsize) {      
			set_error1 ("(init_mem) shorts must be even # of chars");
			return(0);                              
	}                                         
	if (max_type != 0) {
		set_error1 ("(init_mem) second call not allowed");
		return(0);
	}
	if (mt < 1) {
		set_error1 ("(init_mem) invalid maximum number of types");
		return(0);
	}
	max_type = mt;
	size = sizeof (struct class);
	classes = (struct class *) calloc (max_type, size);
	if (classes == (struct class *) NULL) {
		set_error1 ("memory failure in init_mem");
		return(0);
	}
	return(1);
}

int define_type (int type, unsigned long size, char *name)
{
	int c, l, count;
	char message[MAXLINE];

	c = type - 1;
	if (type < 1) {
		sprintf (message, "(define_type): invalid type: %5d", type);
		set_error1(message);
		set_error2 ("minimum allowed type: 1");
		return(0);
	}
	if (type > max_type) {
		sprintf (message, "(define_type): invalid type: %5d", type);
		set_error1(message);
		sprintf (message, "maximum allowed type: %5ld", max_type);
		set_error2(message);
		return(0);
	}
	l = strlen (name);
	if (l < 1 || l >= MAX_TYPE_NAME) {
		set_error1 ("(define_type): invalid type name");
		return(0);
	}
	if ((classes + c) -> size != 0L) {
		sprintf (message, "(define_type): type %d %s already defined",
			type, name);
		set_error1(message);
		return(0);
	}
	if (size < 1) {
		set_error1 ("(define_type): invalid block size");
		return(0);
	}
	(classes + c) -> size = size;
	(classes + c) -> type = type;
	strcpy ((classes + c) -> name, name);
    
	return(1);
}

void get_class_name (int object_type, char object_type_name[MAX_NAME])
{
	int c;
	struct object_header *oh;

	c = object_type - 1;
    if (object_type >= 1 && object_type <= max_type) {
		strcpy (object_type_name, (classes + c) -> name);
    }
	else strcpy (object_type_name, "unknown-type");
}


short *allocate_object (int type)
{
	int c, success;
	unsigned long size;
	short *ptr;

	ptr = (short *) NULL;
	if (type < 1 || type > max_type) {
		set_error1 ("allocate_object: invalid type");
		return (ptr);
	}
	c = type - 1;
	size = (classes + c) -> size;
	if (size == 0) {
		set_error1 ("allocate_object: invalid size");
		return (ptr);
	}
	ptr = get_cache (type);
	if (ptr != NULL) return (ptr);
	ptr = allocate_mem (ONEOBJ, type, size, 1);
	(classes + c) -> count++;
	return (ptr);
}

short *allocate_named_object (int type, int routine, int variable)
{
	int c, success;
	unsigned long size;
	short *ptr;

	ptr = (short *) NULL;
	if (type < 1 || type > max_type) {
		set_error1 ("allocate_object: invalid type");
		return (ptr);
	}
	c = type - 1;
	size = (classes + c) -> size;
	if (size == 0) {
		set_error1 ("allocate_object: invalid size");
		return (ptr);
	}
	ptr = get_cache (type);
	if (ptr != NULL) return (ptr);
	ptr = allocate_mem (ONEOBJ, type, size, 1);
	(classes + c) -> count++;
	variable_count[variable-1]++;
	return (ptr);
}

int free_object (int type, short *ptr)
{
	int c, result;
	unsigned long size;
	char type_name[MAX_TYPE_NAME];
	char message[MAXLINE];

	if (type < 1 || type > max_type) return(0);
	if (ptr == NULL) return (0);
	c = type - 1;
	size = (classes + c) -> size;
	if (size == 0) return(0);
	if (put_cache(type, ptr)) return (1);
	result = free_mem (ONEOBJ, type, ptr);
	(classes + c) -> count--;
	return(result);
}

int free_named_object (int type, short *ptr, int routine, int variable)
{
	int c, result;
	unsigned long size;
	char type_name[MAX_TYPE_NAME];
	char message[MAXLINE];

	if (type < 1 || type > max_type) return(0);
	if (ptr == NULL) return (0);
	c = type - 1;
	size = (classes + c) -> size;
	if (size == 0) return(0);
	if (put_cache(type, ptr)) return (1);
	result = free_mem (ONEOBJ, type, ptr);
	(classes + c) -> count--;
	variable_count[variable-1]--;
	return(result);
}
short *get_cache (int type)
{
	int i;
	struct class *cl;
	short *ptr;

	ptr = (short *) NULL;
	if (type < 1 || type > max_type) {
		return (ptr);
	}
	if (type < 1 || type > max_type) return(0);
	cl = classes + (type-1);
	for (i = 0; i < MAX_CACHE; i++) {
		if (cl -> cache[i] != (void *) NULL) {
			ptr = (short *) (cl -> cache[i]);
			cl -> cache[i] = NULL;
			return (ptr);
		}
	}
	return (ptr);
}

int put_cache (int type, short *ptr)
{
	int i;
	unsigned long size;
	struct class *cl;
	char *cptr, *bptr, *eptr;

	if (type < 1 || type > max_type) return(0);
	if (ptr == NULL) return (0);
	cl = classes + (type-1);
	size = cl -> size;
	for (i = 0; i < MAX_CACHE; i++) {
		if (cl -> cache[i] == (void *) NULL) {
			/* zero the memory */
			bptr = (char *) ptr;
			eptr = bptr + size;
			for (cptr = bptr; cptr < eptr; cptr++)
				*cptr = (char) 0;
			cl -> cache[i] = (void *) ptr;
			return (1);
		}
	}
	return(0);
}

int free_cache (int type)
{
	int i, result;
	struct class *cl;
	short *ptr;

	if (type < 1 || type > max_type) return(0);
	cl = classes + (type-1);
	if (cl -> count <= 0) return (1);
	for (i = 0; i < MAX_CACHE; i++) {
		if (cl -> cache[i] != (void *) NULL) {
			ptr = (short *) (cl -> cache[i]);
			cl -> cache[i] = NULL;
			result = free_mem (ONEOBJ, type, ptr);
			cl -> count--;
		}
	}
	return (1);
}

short *allocate_objects (int type, long count)
{
	int c, success;
	unsigned long size;
	short *ptr;
	struct class *cl;

	ptr = (short *) NULL;
	if (type < 1 || type > max_type) {
		set_error1 ("allocate_object: invalid type");
		return (ptr);
	}
	c = type - 1;
	size = (classes + c) -> size;
	if (size == 0) {
		set_error1 ("allocate_object: invalid size");
		return (ptr);
	}
	cl = (classes + c);
	cl -> count += count;
	return (short *) allocate_mem (OBJECTS, type, size, (unsigned long) count);
}

int free_objects (int type, short *ptr)
{
	int c;
	unsigned long size;
	long count = 0;
	struct class *cl;
	struct premem *pm;
	char type_name[MAX_TYPE_NAME];
	char message[MAXLINE];

	if (type < 1 || type > max_type) return(0);
	if (ptr == NULL) return (0);
	c = type - 1;
	size = (classes + c) -> size;
	if (size == 0) return(0);
	cl = (classes + c);
	if (cl -> count <= 0) return (0);
	pm = (struct premem *) ((char *) ptr - sizeof (struct premem));
	count = pm -> count;
	if (count > 0 && count <= cl -> count)
		cl -> count -= count;
	return free_mem (OBJECTS, type, (short *) ptr);
}

int get_type_name (int type, char type_name[MAX_TYPE_NAME])
{
	if (type > 0 && type <= max_type)
		strcpy (type_name, (classes + type - 1) -> name);
	else strcpy (type_name, "unknown");
	return(1);
}


int copy_bytes (char *from,  char *to, unsigned long n)
{
	int i;
	
	for (i = 0; i < n; i++)
		*(to + i) = *(from + i);
	return(1);
}

void *allocate_pointers (int type, long n_pointers)
{
	int c;
	unsigned long size;
	void *ptrs;
	struct class *cl;

	ptrs = (void *) NULL;
	if (type < 1 || type > max_type) {
		set_error1 ("allocate_pointers: invalid type");
		return (ptrs);
	}
	c = type - 1;
	cl = (classes + c);
	size = cl -> size;
	if (size == 0) {
		set_error1 ("allocate_pointers: invalid size");
		return (ptrs);
	}
	cl -> pointer_count += n_pointers;
	return (void *) allocate_mem (OBJPTRS, type, sizeof (void *), (unsigned long) n_pointers);
}

int free_pointers (int type, void *ptrs)
{
	int c;
	unsigned long size;
	long count = 0;
	struct class *cl;
	struct premem *pm;

	if (ptrs == NULL) return (0);
	if (type < 1 || type > max_type) {
		set_error1 ("free_pointers: invalid type");
		return (0);
	}
	c = type - 1;
	cl = (classes + c);
	size = cl -> size;
	if (size == 0) {
		set_error1 ("free_pointers: invalid size");
		return (0);
	}
	pm = (struct premem *) ((char *) ptrs - sizeof (struct premem));
	count = pm -> count;
	if (count > 0 && count <= cl -> pointer_count)
		cl -> pointer_count -= count;
	return free_mem (OBJPTRS, type, (short *) ptrs);
}

double *allocate_doubles (long n_double, int routine, int variable)
{
	variable_count[variable-1]++;
	return (double *) allocate_mem (DOUBLES, 0, sizeof (double), (unsigned long) n_double);
}

int free_doubles (double *doubles, int routine, int variable)
{
	if (doubles == (double *) NULL) return (0);
	variable_count[variable-1]--;
	return free_mem (DOUBLES, 0, (short *) doubles);
}

float *allocate_floats (long n_float)
{
	return (float *) allocate_mem (FLOATS, 0, sizeof (float), (unsigned long) n_float);
}

int free_floats (float *floats)
{
	if (floats == (float *) NULL) return (0);
	return free_mem (FLOATS, 0, (short *) floats);
}

long *allocate_longs (long n_long, int routine, int variable)
{
	variable_count[variable-1]++;
	return (long *) allocate_mem (LONGS, 0, sizeof (long), (unsigned long) n_long);
}

int free_longs (long *longs, int routine, int variable)
{
	if (longs == (long *) NULL) return (0);
	variable_count[variable-1]--;
	return free_mem (LONGS, 0, (short *) longs);
}

short *allocate_shorts (long n_short)
{
	return (short *) allocate_mem (SHORTS, 0, sizeof (short), (unsigned long) n_short);
}

int free_shorts (short *shorts)
{
	if (shorts == (short *) NULL) return (0);
	return free_mem (SHORTS, 0, (short *) shorts);
}

atomnum *allocate_atomnums (long n_atomnum)
{
	return (atomnum *) allocate_mem (SHORTS, 0, sizeof (atomnum), (unsigned long) n_atomnum);
}

atomnum *reallocate_atomnums (long n_atomnum, atomnum *atomnums)
{
	return (atomnum *) reallocate_mem (SHORTS, 0, sizeof (atomnum), (unsigned long) n_atomnum, (short *) atomnums);
}

int free_atomnums (atomnum *atomnums)
{
	if (atomnums == (atomnum *) NULL) return (0);
	return free_mem (SHORTS, 0, (short *) atomnums);
}

unsigned char *allocate_bytes (long n_byte)
{
	return (unsigned char *) allocate_mem (BYTES, 0, sizeof (unsigned char), (unsigned long) n_byte);
}

int free_bytes (unsigned char *bytes)
{
	if (bytes == (unsigned char *) NULL) return (0);
	return free_mem (BYTES, 0, (short *) bytes);
}

char *allocate_chars (long n_char)
{
	return (char *) allocate_mem (CHARS, 0, sizeof (char), (unsigned long) n_char);
}

int free_chars (char *chars)
{
	if (chars == (char *) NULL) return (0);
	return free_mem (CHARS, 0, (short *) chars);
}

char **allocate_char_pointers (long n_pointers)
{
	return (char **) allocate_mem (STRPTRS, 0, sizeof (char *), (unsigned long) n_pointers);
}

int free_char_pointers (char **ptrs)
{
	if (ptrs == (char **) NULL) return (0);
	return free_mem (STRPTRS, 0, (short *) ptrs);
}

short *allocate_mem (int memtype, int objtype, unsigned long objsize, unsigned long count)
{
	unsigned long size;
	unsigned long adjusted;
	char *ptr;
	struct premem *pm;
	short *rp;
	size = sizeof (struct premem) + count * objsize;
	ptr = calloc (1, size);
	if (ptr == NULL) return (NULL);
	rp = (short *) (ptr + sizeof (struct premem));
	pm = (struct premem *) ptr;
	pm -> memtype = memtype;
	pm -> objtype = objtype;
	pm -> objsize = objsize;
	pm -> count = count;
	memory_count[memtype-1] += count;
	adjusted = size - sizeof (struct premem);
	update_memory (adjusted);
	return (rp);
}

short *reallocate_mem (int memtype, int objtype, unsigned long objsize, unsigned long newcount, short *oldrp)
{
	unsigned long oldsize, newsize, oldcount;
	char *oldptr, *newptr;
	struct premem *oldpm, *newpm;
	short *newrp;

	newrp = NULL;
	if (oldrp == NULL) {
		return (newrp);
	}
	oldptr = (char *) oldrp - sizeof (struct premem);
	oldpm = (struct premem *) oldptr;
	oldcount = oldpm -> count;
	if (oldpm -> memtype != memtype) {
		return (newrp);
	}
	if (oldpm -> objtype != objtype) {
		return (newrp);
	}
	if (oldcount == newcount) {
		return (oldrp);
	}
	oldsize = sizeof (struct premem) + oldcount * objsize;
	newsize = sizeof (struct premem) + newcount * objsize;
	newptr = realloc (oldptr, newsize);
	if (newptr == NULL) {
		return (newrp);
	}
	newpm = (struct premem *) newptr;
	if (newpm -> memtype != memtype) {
		return (newrp);
	}
	if (newpm -> objtype != objtype) {
		return (newrp);
	}
	if (newpm -> objsize != objsize) {
		return (newrp);
	}
	newpm -> count = newcount;
	newrp = (short *) (newptr + sizeof (struct premem));
	memory_count[memtype-1] += (newcount - oldcount);
	update_memory (newsize - oldsize);
	return (newrp);
}

int free_mem (int memtype, int objtype, short *rp)
{
	long count = 0;
	long objsize = 0;
	char *ptr;
	struct premem *pm;
	unsigned long size;
	unsigned long adjusted;

	if (rp == NULL) return (0);
	ptr = (char *) rp - sizeof (struct premem);
	pm = (struct premem *) ptr;
	if (pm -> memtype != memtype) return (0);
	if (pm -> objtype != objtype) return (0);
	count = pm -> count;
	objsize = pm -> objsize;
	size = sizeof (struct premem) + count * objsize;
	/* so not freed again */
	pm -> memtype = 0;
	pm -> objtype = 0;
	memory_count[memtype-1] -= count;
	adjusted = size - sizeof (struct premem);
	update_memory (-adjusted);
	free (ptr);
	return (1);
}

void print_counts ()
{
	int type, c;
	long count, size, pointer_count;
	long memsize;
	struct class *cl;
	char *name;
	char message[MAXLINE];

	free_cache (ATOM);
	free_cache (TORUS);
	free_cache (PROBE);
	free_cache (CIRCLE);
	free_cache (ARC);
	free_cache (VERTEX);
	free_cache (VARIETY);
	free_cache (FACE);
	free_cache (EDGE);
	free_cache (CYCLE);
	free_cache (COMPONENT);
	free_cache (EVALPNT);
	free_cache (PHNVTX);
	free_cache (PHNEDG);
	free_cache (PHNTRI);
	free_cache (HEDVTX);
	free_cache (HEDEDG);
	free_cache (HEDTRI);
	free_cache (POLYGON);
	free_cache (SOLID_ANGLE);
	free_cache (CRITLINK);
	free_cache (SURFACE);
	free_cache (MOLECULE);
	free_cache (OBJECT_SCHEME);
	free_cache (COLOR_RAMP);
	free_cache (MATERIAL_TABLE);
	free_cache (MSSCENE);

	/* free globals, too */
	free_cache(ARRAY);
	free_cache(LEAF);
	free_cache(VERTEX_PAIR);
	free_cache(CIRCLE);
	sprintf (message, "%8ld bytes maximum memory", maximum_memory);
	inform(message);
	sprintf (message, "%8ld bytes of memory leaks", current_memory);
	inform(message);
	for (type = 1; type <= MEMORY_TYPES; type++) {
		c = type - 1;
		count = memory_count[c];
		if (count == 0) continue;
		name = "unknown type";
		switch (type) {
		case ONEOBJ:
			memsize = 0;
			continue;
		case OBJECTS:
			memsize = 0;
			continue;
		case OBJPTRS:
			memsize = 0;
			continue;
		case STRPTRS:
			memsize = sizeof (char *);
			name = "string pointers";
			break;
		case CHARS:
			memsize = sizeof (char);
			name = "characters";
			break;
		case BYTES:
			memsize = sizeof (unsigned char);
			name = "bytes";
			break;
		case SHORTS:
			memsize = sizeof (short);
			name = "shorts";
			break;
		case LONGS:
			memsize = sizeof (long);
			name = "longs";
			break;
		case FLOATS:
			memsize = sizeof (float);
			name = "floats";
			break;
		case DOUBLES:
			memsize = sizeof (double);
			name = "doubles";
			break;
		default:
			memsize = 0;
			continue;
		}
		if (memsize == 0) continue;
		size = count * memsize;
		sprintf (message, "%8ld bytes for %8ld %s",
			size, count,  name);
		inform (message);
	}
	for (type = 1; type <= N_OBJECTS; type++) {
		c = type - 1;
		cl = (classes + c);
		count = cl -> count;
		pointer_count = cl -> pointer_count;
		if (count == 0 && pointer_count == 0) continue;
		size = count * cl -> size + pointer_count * sizeof(void *);
		name = (classes + c) -> name;
		sprintf (message, "%8ld bytes for %8ld [%8ld] %s",
			size, count, pointer_count, name);
		inform (message);
	}
	for (type = 1; type < VARIABLE_TYPES; type++) {
		c = type - 1;
		count = variable_count[c];
		if (count == 0) continue;
		sprintf (message, "%8ld instances for variable %8ld", count, type);
		inform (message);
	}
}

void update_memory (long change) {
	current_memory += change;
	if (current_memory > maximum_memory)
		maximum_memory = current_memory;
}
