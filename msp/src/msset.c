/* 
   Molecular Surface Package
   Protein Shape Library (derived from VOID/NYU)
   Copyright 1992 by Michael L. Connolly
   All Rights Reserved
   
   January 5, 2002
   
*/

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"


int psl_init ()
{
	pi_constant = acos (-1.0);
	allocate_headers ();
    define_set();
    if (error()) return (0);
	define_symbol ();
    if (error()) return (0);
	define_simple();
    if (error()) return (0);
	define_region_type ();
    if (error()) return (0);
	define_mol();
    if (error()) return (0);
	return (1);
}

void define_simple ()
{
	unsigned long size;
	
	size = boolean_size ();
	init_object (BOOLEAN, size);
	if (error ()) return;
		
	size = integer_size ();
	init_object (INTEGER, size);
	if (error ()) return;
	
	size = real_size ();
	init_object (REAL, size);
	if (error ()) return;
	
	size = string_size ();
	init_object (STRING, size);
	if (error ()) return;
	
	memory_init (BOOLEAN);
	if (error ()) {
		return;
	}
	memory_init (INTEGER);
	if (error ()) {
		return;
	}
	memory_init (REAL);
	if (error ()) {
		return;
	}
	memory_init (STRING);
	if (error ()) {
		return;
	}
}


void define_region_type ()
{
	unsigned long size;
	
	
	size = region_size ();
	init_object (REGION, size);
	if (error ()) return;
}

void define_symbol ()
{
	unsigned long size;
	
	size = symbol_size ();
	init_object (SYMBOL, size);
	if (error ()) return;
		
	memory_init (SYMBOL);
	if (error ()) {
		return;
	}
}

void define_set ()
{
	unsigned long size;
	
	size = set_size ();
	init_object (SET, size);
	if (error ()) return;
	
	memory_init (SET);
	if (error ()) {
		return;
	}
}

void define_mol ()
{
	unsigned long size;
	
	size = atom_size ();
	init_object (ATOM, size);
	if (error ()) return;
	
	size = bond_size ();
	init_object (BOND, size);
	if (error ()) return;
	
}

void init_object (int object_type, unsigned long object_size)
{
	struct object_header *oh;
	char message[MAXLINE];
	
	strcpy (function_name, "object_init");
	current_object_type = object_type;

	if (object_size == 0L) {
		sprintf (message, "invalid size: %8ld for type: %5d",
				object_size, object_type);
		set_error1 (message);
		return;
	}
	
	/* set defaults */	
	oh = object_table + (object_type - 1);
	if (oh -> size > 0L) return;	/* already done */
	oh -> deletable = FALSE;
	oh -> definable = FALSE;
	oh -> copyable = TRUE;
	oh -> named = FALSE;
	oh -> number = object_type;
	oh -> size = object_size;
	oh -> initial = STANDARD_INITIAL;
	oh -> increment = STANDARD_INCREMENT;
	oh -> n_allocated = 0;
	oh -> first_block = (struct object_block *) NULL;
	oh -> top = 0;
	
	switch (object_type) {
	case  ATOM:
		strcpy (oh -> name, "atom");
		oh -> named = TRUE;
		oh -> definable = TRUE;
		oh -> deep_copy_only = TRUE;
		oh -> shallow_copy_only = FALSE;
		break;
	case  BOND:
		strcpy (oh -> name, "bond");
		oh -> named = TRUE;
		oh -> definable = TRUE;
		oh -> deep_copy_only = TRUE;
		oh -> shallow_copy_only = FALSE;
		break;
	case  BOOLEAN:
		strcpy (oh -> name, "boolean");
		oh -> named = TRUE;
		oh -> deep_copy_only = TRUE;
		oh -> shallow_copy_only = FALSE;
		break;
	case  INTEGER:
		strcpy (oh -> name, "integer");
		oh -> named = TRUE;
		oh -> deletable = TRUE;
		oh -> deep_copy_only = TRUE;
		oh -> shallow_copy_only = FALSE;
		break;
	case  REAL:
		strcpy (oh -> name, "real");
		oh -> named = TRUE;
		oh -> deep_copy_only = TRUE;
		oh -> shallow_copy_only = FALSE;
		break;
	case  REGION:
		strcpy (oh -> name, "region");
		oh -> named = TRUE;
		oh -> deletable = TRUE;
		oh -> copyable = FALSE;
		break;
	case  SET:
		strcpy (oh -> name, "set");
		oh -> named = TRUE;
		oh -> deletable = TRUE;
		oh -> deep_copy_only = FALSE;
		oh -> shallow_copy_only = FALSE;
		break;
	case  STRING:
		strcpy (oh -> name, "string");
		oh -> named = TRUE;
		oh -> deep_copy_only = TRUE;
		oh -> shallow_copy_only = FALSE;
		break;
	case  SYMBOL:
		strcpy (oh -> name, "symbol");
		oh -> deep_copy_only = TRUE;
		oh -> shallow_copy_only = FALSE;
		break;
	default:
		sprintf (message, "unknown object type: %5d", object_type);
		set_error1 (message);
		return;

	}       /* end of switch */
	
	
}


void allocate_headers ()
{
	object_table = (struct object_header *)
		allocate_objects (OBJECT_HEADER, N_OBJECTS);
	if (object_table == (struct object_header *) NULL) {
		set_error2 ("allocating object table");
		set_error1 ("MEMORY_FAILURE");
		return;
	}
}


void memory_init (int object_type)
{
	int t;
	int temporary_set, all_set, delete_set, free_set;
	struct object_header *oh;
	struct object_block *ob;
	struct set *ptr;
	char message[MAXLINE];

	strcpy (function_name, "memory_init");

	t = object_type - 1;

	/* initial allocation */
	current_object_type = object_type;
	oh = object_table + t;
	oh -> top = 0;
	if (oh -> size <= 0L) return;

	ob = allocate_object_block (object_type, oh -> initial);
	if (error()) return;
	oh -> first_block = ob;
	oh -> last_block = ob;
	oh -> n_allocated = oh -> initial;

	strcpy (function_name, "memory_init");

	if (object_type == SET) {
		/* the first set is the set of all sets */
		oh = object_table + SET - 1;
		all_set = allocate_psl (SET);
		if (all_set == 0) {
			set_error1 ("logic_error: set of all sets");
			return;
		}
		set_set_type (all_set, SET);
		oh -> all_set = all_set;
		/* the set of all sets includes itself as a member */
		include_member (all_set, SET, all_set);
		if (error()) {
			strcpy (function_name, "memory_init");
			set_error2 ("from include_member");
			return;
		}
	}

	/* set up all set */
	oh = object_table + t;		/* oh messed up by set code */
	if (object_type != SET) {
		all_set = allocate_psl (SET);
		if (all_set == 0) {
			set_error1 ("another all set: logic_error");
			return;
		}
		set_set_type (all_set, object_type);
		oh -> all_set = all_set;
		ptr = set_ptr (all_set);
		ptr -> special = TRUE;
	}

	/* set up delete set */
	if (oh -> deletable) {
		delete_set = allocate_psl (SET);
		if (delete_set == 0) {
			set_error1 ("logic_error: delete set");
			return;
		}
		set_set_type (delete_set, t + 1);
		oh -> delete_set = delete_set;
		ptr = set_ptr (delete_set);
		ptr -> special = TRUE;
	}

	/* set up free set */
	if (oh -> deletable) {
		free_set = allocate_psl (SET);
		if (free_set == 0) {
			set_error1 ("logic error: free set");
			return;
		}
		set_set_type (free_set, t + 1);
		oh -> free_set = free_set;
		ptr = set_ptr (free_set);
		ptr -> special = TRUE;
	}

	/* set up temporary set */
	temporary_set = allocate_psl (SET);
	if (temporary_set == 0) {
		set_error1 ("logic_error: temporary set");
		return;
	}
	set_set_type (temporary_set, t + 1);
	oh -> temporary_set = temporary_set;
	ptr = set_ptr (temporary_set);
	ptr -> special = TRUE;
}

struct object_block *allocate_object_block(int object_type, int increment)
{
	long long_size;
	unsigned old_size;
	struct object_header *oh;
	struct object_block *b;

	strcpy (function_name, "allocate_object_block");
	oh = object_table + object_type - 1;
	old_size = oh -> n_allocated;
	old_size *= oh -> size;
	long_size = oh -> n_allocated + increment;
	long_size *= oh -> size;

	/* add a block at the end */
	b = (struct object_block *)
		allocate_object(OBJECT_BLOCK);
	if (b == (struct object_block *) NULL) {
		set_error1 ("MEMORY_FAILURE");
		return (NULL);
	}
	b -> start = (char *) allocate_objects (object_type, increment);
	if (b -> start == (char *) NULL) {
		set_error1 ("MEMORY_FAILURE");
		return (NULL);
	}
	b -> first_object = oh -> n_allocated + 1;
	b -> last_object = oh -> n_allocated + increment;
	return (b);
}

int free_object_block (int object_type, struct object_block *b)
{
	free_objects (object_type, (short *) (b -> start));
	if (error()) return(0);
	free_object (OBJECT_BLOCK, (short *) b);
	if (error()) return(0);
	return (1);
}

int free_object_blocks (int object_type)
{
	struct object_block *b, *nb;
	struct object_header *oh;

	oh = object_table + object_type - 1;
	for (b = oh -> first_block; b != NULL; b = nb) {
		nb = b -> next_block;
		free_object_block (object_type, b);
		if (error()) return(0);
	}
	return (1);
}

int free_object_headers ()
{
	free_objects (OBJECT_HEADER, (short *) object_table);
	if (error()) return (0);
	return (1);
}

int free_all_psl ()
{
	struct element_block *blk, *next_blk;

	for (blk = first_element_block; blk != NULL; blk = next_blk) {
		next_blk = blk -> next_block;
		free_objects (ELEMENT, (short *) (blk -> first_element));
		if (error()) return (0);
		free_object (ELEMENT_BLOCK, (short *) blk);
		if (error()) return (0);
	}

	free_cache(ELEMENT_BLOCK);

	free_object_blocks (BOOLEAN);
	free_object_blocks (INTEGER);
	free_object_blocks (REAL);
	free_object_blocks (STRING);
	free_object_blocks (REGION);
	free_object_blocks (SYMBOL);
	free_object_blocks (SET);
	free_cache(OBJECT_BLOCK);
	free_object_headers();
	return (1);
}

int allocate_psl (int object_type)
{
	int idx, object_number, free_set, object_size, all_set;
	int use_free_set;
	long long_size;
	unsigned old_size;
	unsigned char *start_address, *end_address, *uc;
	struct object_block *b;
	struct object_header *oh;
	char message[MAX_STRING];
	char object_type_name[MAX_NAME];

	strcpy (function_name, "allocate_psl");
	current_object_type = object_type;
	oh = object_table + object_type - 1;
	object_size = oh -> size;
	all_set = oh -> all_set;

	free_set = oh -> free_set;
	use_free_set = some_free (object_type);
	if (error ()) {
		get_object_type_name(object_type, object_type_name);
		sprintf (message, "no objects of type %s", object_type_name);
		set_error2 (message);
		return (0);
	}
	if (use_free_set) {
		object_number = get_first_member (free_set);
		if (object_number == 0) {
			strcpy (function_name, "allocate_psl");
			sprintf (message,
				"invalid free set, object type = %d", object_type);
			set_error1 (message);
			return(0);
		}
		current_object_number = object_number;
		if (object_number < 1 || object_number > oh -> top) {
			strcpy (function_name, "allocate_psl");
			sprintf (message,
				"invalid object number, object type = %d", object_type);
 			set_error1 (message);
		    return (0);
		}
		/* must mark as valid object number before marking not free */
		include_member (all_set, object_type, object_number);
		if (error ()) {
			strcpy (function_name, "allocate_psl");
			set_fatal ();
			return (0);
		}
		exclude_member (free_set, object_type, object_number);
		if (error ()) {
			strcpy (function_name, "allocate_psl");
			sprintf (message,
				"exclude from free set fails for object type = %d",
					object_type);
  						set_error2 (message);
			return (0);
		}
		oh -> n_free = count_set (free_set);
		/* ZERO the old object's memory */
		start_address = (unsigned char *) generic_ptr (object_type,
			object_number);
		if (error ()) return (0);
		end_address = start_address + object_size;
		for (uc = start_address; uc < end_address; uc++)
			*uc = 0;
		return (object_number);
	}

	current_object_number = oh -> top;
	if (oh -> top >= oh -> n_allocated) {
		b = allocate_object_block (object_type, oh -> increment);
		if (error()) return(0);
		oh -> last_block -> next_block = b;
		oh -> last_block = b;
		oh -> n_allocated += oh -> increment;
	}

	idx = oh -> top;
	oh -> top++;
	object_number = idx + 1;
	/* don't expect the set of all sets to have an all_set yet */
	if (! (object_type == SET && object_number == 1)) {
		include_member (all_set, object_type, object_number);
		if (error ()) {
			strcpy (function_name, "allocate_psl");
			set_fatal ();
			return (0);
		}
	}
	return (object_number);
}

long boolean_size ()
{
	return (sizeof (struct boolean));
}

long integer_size ()
{
	return (sizeof (struct integer));
}

long real_size ()
{
	return (sizeof (struct real));
}

long string_size ()
{
	return (sizeof (struct string));
}


int some_free (int object_type)
{
	int free_set;
	struct object_header *oh;

		strcpy (function_name, "some_free");
	oh = object_table + object_type - 1;
	free_set = oh -> free_set;
	if (free_set == 0) return (FALSE);
	if (count_set (free_set) > 0) return (TRUE);
	return (FALSE);
}

int allocate_set (int set_type)
{
	int object_number;

		strcpy (function_name, "allocate_set");
	object_number = allocate_psl (SET);
	if (object_number == 0) return (0);
	set_set_type (object_number, set_type);
	return (object_number);
}

int last_of_objects (int object_type)
{
	int last, all_set;
	struct object_header *oh;

	strcpy (function_name, "last_of_objects");
	oh = object_table + object_type - 1;
	all_set = oh -> all_set;
	last = get_last_member (all_set);
	return (last);
}

int number_of_all_objects ()
{
	int tn, n;
	int object_type;

	tn = 0;
	for (object_type = 1; object_type <= N_OBJECTS;
		object_type++) {
		n = number_of_objects (object_type);
		if (error ()) return (0);
		tn += n;
	}
	return (tn);
}

int number_of_object_types ()
{
	return (N_OBJECTS);
}

int number_of_objects (int object_type)
{
	int n, all_set;
	struct object_header *oh;

	strcpy (function_name, "number_of_objects");
	current_object_type = object_type;
	oh = object_table + object_type - 1;
	all_set = oh -> all_set;
	n = count_set (all_set);
	return (n);
}

int object_exists (int object_type, int object_number)
{
	int all_set, answer;
	struct object_header *oh;

	strcpy (function_name, "object_exists");
	if (object_number < 1) return (FALSE);
	oh = object_table + object_type - 1;
	if (object_number > oh -> top) return (FALSE);
	if (oh -> n_free <= 0) return (TRUE);
	if (!oh -> deletable) return (TRUE);

	current_object_type = object_type;
	current_object_number = object_number;
	all_set = oh -> all_set;
	if (all_set <= 0) return (TRUE);
	/* inquiry about set, take a chance, to save time */
	if (object_type == SET) return (oh -> top >= 1);
		answer = member_of (object_number, all_set);
	return (answer);
}

int type_deletable (int object_type)
{
       int answer;
       struct object_header *oh;
       char message[MAXLINE];

       strcpy (function_name, "type_deletable");
       if (object_type < 1 || object_type > N_OBJECTS) {
		sprintf (message,
			"type_deletable: invalid type: %d", object_type);
		set_error1 (message);
		current_object_type = 0;
		current_object_number = 0;
		return (0);
	}
	oh = object_table + object_type - 1;
	answer = oh -> deletable;
	return (answer);
}


/* return the true or false value of the given boolean object */

int is_true (int boolean_number)
{
	struct boolean *b;

	b = boolean_ptr (boolean_number);
	if (b == (struct boolean *) NULL) return (FALSE);
	return (b -> value);
}

/* return the negation of the value of the given boolean object */

int is_false (int boolean_number)
{
	int value;
	struct boolean *b;

	b = boolean_ptr (boolean_number);
	if (b == (struct boolean *) NULL) return (FALSE);
	value = !b -> value;
	return (value);
}

/* return the value of the specified integer object */

int integer_is (int number)
{
	struct integer *int_ptr;

	int_ptr = integer_ptr (number);
	if (int_ptr == (struct integer *) NULL) return (FALSE);
	return (int_ptr -> value);
}

/* return the value of the specified real object */

double real_is (int number)
{
	struct real *r_ptr;

	r_ptr = real_ptr (number);
	if (r_ptr == (struct real *) NULL) return (0.0);
	return (r_ptr -> value);
}

/* return a pointer to the string of the specified string object */

char *string_is (int number)
{
	struct string *s_ptr;

	s_ptr = string_ptr (number);
	if (s_ptr == (struct string *) NULL) {
		set_error1("NULL_STRING");
		return (NULL);
	}
	return ((char *) (s_ptr -> string_value));
}


struct boolean *boolean_ptr (int number)
{
	return ((struct boolean *) generic_ptr (BOOLEAN, number));
}


struct integer *integer_ptr (int number)
{
	return ((struct integer *) generic_ptr (INTEGER, number));
}


struct real *real_ptr (int number)
{
	return ((struct real *) generic_ptr (REAL, number));
}


struct string *string_ptr (int number)
{
	return ((struct string *) generic_ptr (STRING, number));
}



/* element code */

int element_size ()
{
	int size;
	
	size = sizeof (struct element);
	return (size);
}

struct element *new_element ()
{
	int r_size, e_size;
	struct element_block *blk, *last_blk;
	struct element *elem;
	char message[MAXLINE];

	strcpy (function_name, "new_element");
	/* general error_checking */
	r_size = range_size ();
		e_size = element_size ();
	if (r_size != e_size) {
		sprintf (message, "range size = %5d, element size = %5d",
			r_size, e_size);
 		set_error1 (message);
	    return (NULL);
	}

	last_blk = (struct element_block *) NULL;
	for (blk = first_element_block; blk != (struct element_block *) NULL;
		blk = blk -> next_block) {
		last_blk = blk;
		for (elem = blk -> first_element;
			elem <= blk -> last_element; elem++) {
			/* this algorithm depends on some field being set to non-zero */
			if (elem -> object_number != 0) continue;
			if (elem -> coefficient != 0) continue;
			if (elem -> next != NULL) continue;
			elem -> object_number = 1;      /* mark as used */
			return (elem);
		}
	}

	/* allocate new block */
	blk = (struct element_block *) allocate_object (ELEMENT_BLOCK);
	if (blk == (struct element_block *) NULL) {
		set_error1 ("MEMORY_FAILURE");
		set_error2 ("cannot allocate element block");
		return (NULL);
	}

	if (last_blk != (struct element_block *) NULL)
		last_blk -> next_block = blk;
	else first_element_block = blk;

	blk -> n_element = ELEMENT_BLOCK_SIZE;
	elem = (struct element *) allocate_objects (ELEMENT, blk -> n_element);
	if (elem == NULL) {
		set_error1 ("MEMORY_FAILURE");
		set_error2 ("cannot allocate element memory");
		return (NULL);
	}
	blk -> first_element = elem;
	blk -> last_element = blk -> first_element + (blk -> n_element - 1);
	elem -> object_number = 1;      /* mark as used */
	return (elem);
}

void free_element (struct element *returned_elem)
{
	struct element_block *blk;
	struct element *elem;

	strcpy (function_name, "free_element");
	for (blk = first_element_block; blk != (struct element_block *) NULL;
		blk = blk -> next_block) {
		if (returned_elem < blk -> first_element) continue;
		if (returned_elem > blk -> last_element) continue;
		for (elem = blk -> first_element;
			elem <= blk -> last_element; elem++) {
			if (elem == returned_elem) {
				elem -> object_number = 0;
				elem -> coefficient = 0;
				elem -> next = NULL;
				return;
			}
		}
	}
	set_error1 ("non-allocated element returned");
}

struct generic *generic_ptr (int object_type, int number)
{
	char *cptr;
	struct object_block *b;
	struct object_header *oh; struct generic *ptr;
	char message[MAXLINE];

	/* strcpy (function_name, "generic_ptr"); */
	if (error())
		return ((struct generic *) NULL);
	if (object_type < 1 || object_type > N_OBJECTS) {
		sprintf (message,
			"generic_ptr: invalid type: %d", object_type);
		set_error1 (message);
		set_fatal();
		current_object_type = 0;
		current_object_number = 0;
		return ((struct generic *) NULL);
	}
	current_object_type = object_type;
	current_object_number = number;

	oh = object_table + object_type - 1;
		if (oh -> size <= 0L) {
		sprintf (message, "generic_ptr: invalid type: %d", object_type);
		set_error1 (message);
		set_fatal();
		current_object_type = 0;
		current_object_number = 0;
		return ((struct generic *) NULL);
		}

	if (number <= 0 || number > oh -> top) {
		sprintf (message, "non-existent object number = %d", number);
		set_error1 (message);
		current_object_type = object_type;
		current_object_number = number;
		return ((struct generic *) NULL);
	}

	/* more elaborate check requires time */
	if (debug) {
		if (!object_exists (object_type, number)) {
			sprintf (message, "non-existent object number = %d", number);
			set_error1 (message);
			current_object_type = object_type;
			current_object_number = number;
			return ((struct generic *) NULL);
		}
	}
	current_object_type = object_type;
	current_object_number = number;

	for (b = oh -> first_block; b != (struct object_block *) NULL;
		b = b -> next_block) {
		if (number < b -> first_object) {
			strcpy (message, "generic_ptr: object not in block");
			set_error1 (message);
			current_object_type = object_type;
			current_object_number = number;
			return ((struct generic *) NULL);
		}
		if (number > b -> last_object) continue;
		cptr = b -> start + (number - b -> first_object) * oh -> size;
		ptr = (struct generic *) cptr;
		return (ptr);
	}

	current_object_type = object_type;
	current_object_number = number;
	sprintf (message, "generic_ptr: object %d of type %d not found",
		number, object_type);
	set_error1(message);
	return ((struct generic *) NULL);
}

struct generic *named_ptr (int object_type, int number)
{
	struct object_header *oh;
	struct generic *ptr;

	strcpy (function_name, "named_ptr");
	current_object_type = object_type;
	current_object_number = number;
	oh = object_table + object_type - 1;
	if (!oh -> named) {
		set_error1 ("OBJECT_NOT_NAMEABLE");
		return ((struct generic *) NULL);
	}
	ptr = generic_ptr (object_type, number);
	return (ptr);
}

void get_object_type_name (int object_type, char object_type_name[MAX_NAME])
{
	struct object_header *oh;

	if (object_type > 0) oh = object_table + object_type - 1;
	else oh = (struct object_header *) NULL;
	if (oh != (struct object_header *) NULL) {
		strcpy (object_type_name, oh -> name);
    }
	else strcpy (object_type_name, "unknown-type");
}

int number_of_symbols ()
{
	int n_symbol, all_set;
	struct object_header *oh;
	struct set *all_ptr;

	oh = object_table + SYMBOL - 1;
	all_set = oh -> all_set;
	all_ptr = set_ptr (all_set);
	if (all_ptr == (struct set *) NULL) return (0);
	n_symbol = all_ptr -> n_element;
	return (n_symbol);
}

/* SET routines */

struct set *set_ptr (int number)
{
	return ((struct set *) generic_ptr (SET, number));
}

void init_scratch (int set_type)
{
	if (scratch1 == 0) {
		scratch1 = allocate_set (set_type);
		if (error ()) return;
	}
	else set_set_type (scratch1, set_type);
	
	if (scratch2 == 0) {
		scratch2 = allocate_set (set_type);
		if (error ()) return;
	}
	else set_set_type (scratch2, set_type);
}

int range_size ()
{
	int size;
	
	size = sizeof (struct range);
	return (size);
}

int set_size ()
{
	int size;
	
	size = sizeof (struct set);
	return (size);
}

void set_set_type (int set_number, int set_type)
{
	struct set *ptr;
	
	ptr = set_ptr (set_number);
	if (error()) return;
	if (ptr == (struct set *) NULL) {
		strcpy (function_name, "set_set_type");
		set_error1 ("logic error; set_set_type");
		return;
	}
    ptr -> type = set_type;
	if (set_type == SET)
    ptr -> special = TRUE;
}

get_set_type (int set_number)
{
	int type;
	struct set *ptr;
	
	ptr = set_ptr (set_number);
	if (ptr == (struct set *) NULL) return (0);
    type = ptr -> type;
	return (type);
}

void clear_set (int set_number)
{
	struct set *s_ptr;
	struct range *r, *next;

	s_ptr = set_ptr (set_number);
	if (s_ptr == (struct set *) NULL) return;

	next = NULL;
	for (r = s_ptr -> first; r != (struct range *) NULL; r = next) {
		next = r -> next;
		delete_range (r);
	}

	s_ptr -> n_element = 0;
	s_ptr -> first = (struct range *) NULL;
	clear_for (set_number);
}

int get_first_member (int set_number)
{
	int object_number;
	struct set *s_ptr;
	struct range *r;

	s_ptr = set_ptr (set_number);
	if (s_ptr == (struct set *) NULL) return (0);
	r = s_ptr -> first;
	if (r == (struct range *) NULL) return (0);
	object_number = r -> start;
	return (object_number);
}

int get_last_member (int set_number)
{
	int last;
	struct set *s_ptr;
	struct range *r;

	s_ptr = set_ptr (set_number);
	if (s_ptr == (struct set *) NULL) return (0);
	last = 0;
	for (r = s_ptr -> first; r != (struct range *) NULL; r = r -> next)
		last = r -> end;
	return (last);
}

struct range *make_range (int start, int end)
{
	struct range *r;

	if (start < 1 || start > end) {
		set_error1 ("logic error: make_range");
		return ( (struct range *) NULL);
	}
	r = (struct range *) new_element ();
	if (r == (struct range *) NULL) {
		set_error1("MEMORY_FAILURE");
		return ( (struct range *) NULL);
	}
	n_ranges++;
	r -> start = start;
	r -> end = end;
		r -> next = (struct range *) NULL;
	return (r);
}

void delete_range (struct range *r)
{
	free_element ((struct element *) r);
	n_ranges--;
}

void consolidate (int set_number)
{
	int n_before, n_after;
	struct set *s_ptr;
	struct range *r, *rnext;

	n_before = count_set (set_number);
	if (error ()) return;
	s_ptr = set_ptr (set_number);
	if (s_ptr == (struct set *) NULL) return;
	r = s_ptr -> first;
	while (TRUE) {
		if (r == NULL) break;
		rnext = r -> next;
		if (rnext == NULL) break;
		if (r -> end >= rnext -> start - 1) {
			/* consolidate ranges */
			r -> end = rnext -> end;
			r -> next = rnext -> next;
			delete_range (rnext);
		}
		else r = r -> next;
	}
	n_after = count_set (set_number);
	if (n_before != n_after) {
		set_error1 ("logic error: before != after");
		return;
	}
	s_ptr -> n_element = n_after;
	clear_for (set_number);
}

int count_set (int set_number)
{
	int n_element;
	struct range *r;
	struct set *s_ptr;

	s_ptr = set_ptr (set_number);
	if (s_ptr == (struct set *) NULL) return (0);
	n_element = 0;
	for (r = s_ptr -> first; r != (struct range *) NULL; r = r -> next)
		n_element += (r -> end - r -> start + 1);
	return (n_element);
}


void include_member (int set_number, int object_type, int object_number)
{
	struct set *s_ptr;
	struct range *r, *nr, *next, *prev;

	s_ptr = set_ptr (set_number);
	if (error()) return;
	if (s_ptr == (struct set *) NULL) {
		set_error1("invalid set pointer");
		return;
	}
	if (s_ptr -> type != object_type) {
		set_error1("inconsistent_type");
		return;
	}
	prev = NULL;
	next = NULL;
	for (r = s_ptr -> first; r != (struct range *) NULL; r = next) {
		next = r -> next;
		if (object_number < r -> start - 1) {
			/* new range */
			nr = make_range (object_number, object_number);
			if (nr == (struct range *) NULL) return;
			if (prev == NULL) s_ptr -> first = nr;
			else prev -> next = nr;
			nr -> next = r;
			s_ptr -> n_element++;
			return;
		}
		else if (object_number == r -> start - 1) {
			r -> start = object_number;
			s_ptr -> n_element++;
			return;
		}
		else if (object_number <= r -> end) return;
		else if (object_number == r -> end + 1) {
			if (next == NULL || object_number < next -> start - 1) {
				r -> end = object_number;
				s_ptr -> n_element++;
				return;
			}
			else {
				r -> end = next -> end;
				r -> next = next -> next;
				s_ptr -> n_element++;
				delete_range (next);
				next = r -> next;
				return;
			}
		}
		prev = r;
	}  /* end of for loop */

	/* new range beyond last range */
	nr = make_range (object_number, object_number);
	if (nr == (struct range *) NULL) {
		set_error1("make_range fails");
		return;
	}
	if (prev == NULL) s_ptr -> first = nr;
	else prev -> next = nr;
	nr -> next = NULL;
	s_ptr -> n_element++;
	clear_for (set_number);
}


void exclude_member (int set_number, int object_type, int object_number)
{
	struct set *s_ptr;
	struct range *r, *next, *prev, *nr;

	/* it is okay to exclude a non_member (nop) */
	/* it is also to okay to exclude a freed object */

	s_ptr = set_ptr (set_number);
	if (s_ptr == (struct set *) NULL) return;
	if (s_ptr -> type != object_type) {
		set_error1("INCONSISTENT_TYPE");
		return;
	}
	prev = NULL;
	next = NULL;
	for (r = s_ptr -> first; r != (struct range *) NULL; r = next) {
		next = r -> next;
		if (object_number < r -> start) return;
		else if (object_number > r -> end) {
			prev = r;
			continue;
		}
		else if (object_number == r -> start && object_number == r -> end) {
			if (prev == NULL) s_ptr -> first = next;
			else prev -> next = next;
			delete_range (r);
			s_ptr -> n_element--;
			if (s_ptr -> n_element < 0) {
				set_error1("logic_error");
				set_error2("n_element");
				return;
			}
			return;
		}
		else if (object_number == r -> start) {
			r -> start++;
			s_ptr -> n_element--;
			if (s_ptr -> n_element < 1) {
				set_error1("logic_error");
				set_error2("n_element < 1");
				return;
			}
			return;
		}
		else if (object_number == r -> end) {
			r -> end--;
			s_ptr -> n_element--;
			if (s_ptr -> n_element < 1) {
				set_error1("logic_error");
				set_error2("n_element < 1");
				return;
			}
			return;
		}
		else {  /* extract from middle of range */
			nr = make_range (object_number + 1, r -> end);
			if (nr == (struct range *) NULL) return;
			r -> end = object_number - 1;
			nr -> next = r -> next;
			r -> next = nr;
			if (--s_ptr -> n_element < 2) {
				set_error1("logic_error");
				set_error2("n_element < 2");
				return;
			}
			return;

		}
	} /* end of for loop */
	clear_for (set_number);
}

/* the routine does not depend on no deleted objects */
void set_complement (int set1, int set2)
{
	int start, end, object_number, object_type, exists;
	int excluded_object, change, deletable, type1, type2;
	struct set *sp2;
	struct range *r2;
	char message[MAXLINE];

	strcpy (function_name, "set_complement");
	type1 = get_set_type (set1);
	type2 = get_set_type (set2);
	strcpy (function_name, "set_complement");
	if (type1 != type2) {
			set_error1("INCONSISTENT SET TYPES");
			sprintf (message, "type1 = %5d, type2 = %5d", type1, type2);
			set_error2(message);
			return;
	}
	object_type = type1;
	sp2 = set_ptr (set2);
	if (sp2 == (struct set *) NULL) return;
	pseudo_set_complement (set1, set2);
	deletable = type_deletable (object_type);
	if (error()) return;
	if (!deletable) return;
		

	/* consider the all set */

	excluded_object = 0;
	change = TRUE;
	while (change) {
		change = FALSE;
		for (r2 = sp2 -> first; r2 != (struct range *) NULL; r2 = r2 -> next) {
			start = r2 -> start;
			end = r2 -> end;
			if (end <= excluded_object) continue;
			for (object_number = start; object_number <= end; object_number++) {
				if (object_number <= excluded_object) continue;
				exists = object_exists (object_type, object_number);
				if (error ()) return;
				if (!exists) {
					exclude_member (set2, object_type, object_number);
					if (error ()) return;
					excluded_object = object_number;
					change = TRUE;
				}
				if (change) break;
			}
			if (change) break;
		}
	}
	clear_for (set2);
}

/* this function produces a set that may include freed objects */
void pseudo_set_complement (int set1, int set2)
{
	int start, end, last_object, n1, n2;
	struct set *sp1, *sp2;
	struct range *r1, *r2;
	struct range *first, *prev;

	strcpy (function_name, "pseudo_set_complement");
	sp1 = set_ptr (set1);
	if (sp1 == (struct set *) NULL) return;
	sp2 = set_ptr (set2);
	if (sp2 == (struct set *) NULL) return;
	last_object = last_of_objects (sp1 -> type);
	n1 = count_set (set1);
	first = NULL; prev = NULL;
	start = 1;
		for (r1 = sp1 -> first; r1 != (struct range *) NULL; r1 = r1 -> next) {
		end = r1 -> start - 1;
		if (start <= end) {
			r2 = make_range (start, end);
			if (r2 == (struct range *) NULL) return;
			if (prev == NULL) first = r2;
			else prev -> next = r2;
			prev = r2;
		}
		start = r1 -> end + 1;
	}
	if (start <= last_object) {
		r2 = make_range (start, last_object);
		if (r2 == (struct range *) NULL) return;
		if (prev == NULL) first = r2;
		else prev -> next = r2;
	}
	clear_set (set2);       /* this may be set1 */
	sp2 -> first = first;
	n2 = count_set (set2);
	sp2 -> n_element = n2;
	if (n1 + n2 != last_object) {
		set_error1("logic_error: last_object");
		return;
	}
	clear_for (set2);
}

void set_union (int set1, int set2, int set3)
{
	int min_start, max_end, up_to;
	struct set *sp1, *sp2, *sp3;
	struct range *r1, *r2, *r3;
	struct range *first, *prev;

	sp1 = set_ptr (set1);
	if (sp1 == (struct set *) NULL) return;
	sp2 = set_ptr (set2);
	if (sp2 == (struct set *) NULL) return;
	sp3 = set_ptr (set3);
	if (sp3 == (struct set *) NULL) return;

	first = NULL; prev = NULL;

	r1 = sp1 -> first;
	r2 = sp2 -> first;
	up_to = 0;

	while (TRUE) {
		while (TRUE) {
			if (r1 != NULL && r1 -> end <= up_to) r1 = r1 -> next;
			else if (r2 != NULL && r2 -> end <= up_to) r2 = r2 -> next;
			else break;
		}
		if (r1 == NULL || r2 == NULL) break;
		if (r1 -> start > r2 -> end + 1) {
			r3 = make_range (r2 -> start, r2 -> end);
			up_to = r3 -> end;
			if (r3 == (struct range *) NULL) return;
			if (prev == NULL) first = r3;
			else prev -> next = r3;
			prev = r3;
			continue;
		}
		else if (r2 -> start > r1 -> end + 1) {
			r3 = make_range (r1 -> start, r1 -> end);
			up_to = r3 -> end;
			if (r3 == (struct range *) NULL) return;
			if (prev == NULL) first = r3;
			else prev -> next = r3;
			prev = r3;
			continue;
		}
		min_start = ((r1 -> start < r2 -> start) ? r1 -> start : r2 -> start);
		max_end = ((r1 -> end > r2 -> end) ? r1 -> end : r2 -> end);
		r3 = make_range (min_start, max_end);
		up_to = r3 -> end;
		if (r3 == (struct range *) NULL) return;
		if (prev == NULL) first = r3;
		else prev -> next = r3;
		prev = r3;
	}

	/* finish up 1 or 2 */
	for (; r1 != NULL; r1 = r1 -> next) {
		r3 = make_range (r1 -> start, r1 -> end);
		if (r3 == (struct range *) NULL) return;
		if (prev == NULL) first = r3;
		else prev -> next = r3;
		prev = r3;
	}
	for (; r2 != NULL; r2 = r2 -> next) {
		r3 = make_range (r2 -> start, r2 -> end);
		if (r3 == (struct range *) NULL) return;
		if (prev == NULL) first = r3;
		else prev -> next = r3;
		prev = r3;
	}

	clear_set (set3);       /* this may be set1 or set2 */
	sp3 -> first = first;
	consolidate (set3);
	clear_for (set3);
}

void set_intersection (int set1, int set2, int set3)
{
	int max_start, min_end, up_to;
	struct set *sp1, *sp2, *sp3;
	struct range *r1, *r2, *r3;
	struct range *first, *prev;

	sp1 = set_ptr (set1);
	if (sp1 == (struct set *) NULL) return;
	sp2 = set_ptr (set2);
	if (sp2 == (struct set *) NULL) return;
	sp3 = set_ptr (set3);
	if (sp3 == (struct set *) NULL) return;

	first = NULL; prev = NULL;

	r1 = sp1 -> first;
	r2 = sp2 -> first;

	up_to = 0;
	while (TRUE) {
		while (TRUE) {
			if (r1 != NULL && r1 -> end <= up_to) r1 = r1 -> next;
			else if (r2 != NULL && r2 -> end <= up_to) r2 = r2 -> next;
			else break;
		}
		if (r1 == NULL || r2 == NULL) break;
		if (r1 -> start > r2 -> end) {
			up_to = r1 -> start - 1;
			continue;
		}
		else if (r2 -> start > r1 -> end) {
			up_to = r2 -> start - 1;
			continue;
		}
		max_start = ((r1 -> start > r2 -> start) ? r1 -> start : r2 -> start);
		min_end = ((r1 -> end < r2 -> end) ? r1 -> end : r2 -> end);
		r3 = make_range (max_start, min_end);
		up_to = r3 -> end;
		if (r3 == (struct range *) NULL) return;
		if (prev == NULL) first = r3;
		else prev -> next = r3;
		prev = r3;
	}

	clear_set (set3);       /* this may be set1 or set2 */
	sp3 -> first = first;
	consolidate (set3);
	clear_for (set3);
}

int set_comparison (int number1, int number2, int modifier)
{
	int return_value, set_type;
	struct set *s_ptr1, *s_ptr2, *sc_ptr1, *sc_ptr2;

	strcpy (function_name, "set_comparison");

	sc_ptr1 = set_ptr (scratch1);
	if (sc_ptr1 == (struct set *) NULL) return (FALSE);
	sc_ptr2 = set_ptr (scratch2);
	if (sc_ptr2 == (struct set *) NULL) return (FALSE);

	s_ptr2 = set_ptr (number2);
	if (s_ptr2 == (struct set *) NULL) return (FALSE);
	set_type = s_ptr2 -> type;
		init_scratch (set_type);
	sc_ptr1 -> type = set_type;
	sc_ptr2 -> type = set_type;
	clear_set (scratch1);
	clear_set (scratch2);

	if (modifier == MEMBER_OF) {
		return_value = member_of (number1, number2);
		if (error ()) return (FALSE);
		return (return_value);
	}
	else if (modifier == EQUAL_TO) {
		s_ptr1 = set_ptr (number1);
		if (s_ptr1 == (struct set *) NULL) return (FALSE);
		if (s_ptr1 -> n_element != s_ptr2 -> n_element) return (FALSE);
		if (s_ptr1 -> n_element == 0) return (TRUE);
		set_intersection (number1, number2, scratch1);
		return_value = (sc_ptr1 -> n_element == s_ptr1 -> n_element);
		clear_set (scratch1);
		sc_ptr1 -> type = 0;
		return (return_value);
	}
	else if (modifier == NOT_EQUAL) {
		s_ptr1 = set_ptr (number1);
		if (s_ptr1 == (struct set *) NULL) return (FALSE);
		if (s_ptr1 -> n_element != s_ptr2 -> n_element) return (TRUE);
		if (s_ptr1 -> n_element == 0) return (FALSE);
		set_intersection (number1, number2, scratch1);
		return_value = (sc_ptr1 -> n_element != s_ptr1 -> n_element);
		clear_set (scratch1);
		sc_ptr1 -> type = 0;
		return (return_value);
	}
	else if (modifier == SUBSET_OF) {
		s_ptr1 = set_ptr (number1);
		if (s_ptr1 == (struct set *) NULL) return (FALSE);
		if (s_ptr1 -> n_element > s_ptr2 -> n_element) return (FALSE);
		if (s_ptr1 -> n_element <= 0) return (TRUE);
		set_complement (number2, scratch1);
		if (sc_ptr1 -> n_element <= 0) return (TRUE);
		set_intersection (number1, scratch1, scratch2);
		return_value = (sc_ptr2 -> n_element <= 0);
		clear_set (scratch1);
		sc_ptr1 -> type = 0;
		clear_set (scratch2);
		sc_ptr2 -> type = 0;
		return (return_value);
	}
	else {
		set_error1("LOGIC_ERROR: else");
		return (FALSE);
	}
}

int member_of (int object_number, int set_number)
{
	struct set *s_ptr;
	struct range *r;

	s_ptr = set_ptr (set_number);
	if (s_ptr == (struct set *) NULL) return (FALSE);

	for (r = s_ptr -> first; r != (struct range *) NULL;
		r = r -> next) {
		if (object_number < r -> start) return (FALSE);
		else if (object_number <= r -> end) return (TRUE);
	}
	return (FALSE);
}

void set_subtraction (int set1, int set2, int set3)
{
	int set_type;
	struct set *s_ptr1, *sc_ptr1;

	strcpy (function_name, "set_subtraction");

		
	s_ptr1 = set_ptr (set1);
	if (s_ptr1 == (struct set *) NULL) return;
	strcpy (function_name, "set_subtraction");
	set_type = s_ptr1 -> type;

	init_scratch (set_type);
	if (error ()) return;
	strcpy (function_name, "set_subtraction");
	sc_ptr1 = set_ptr (scratch1);
	if (sc_ptr1 == (struct set *) NULL) return;
	sc_ptr1 -> type = set_type;
	clear_set (scratch1);

	/* utilize pre-existing routines (none of which use scratch sets) */
	set_complement (set2, scratch1);
	if (error ()) return;
	set_intersection (set1, scratch1, set3);
	if (error ()) return;

	/* finished with scratch 1 */
	clear_set (scratch1);
	sc_ptr1 -> type = 0;
	clear_for (set3);
}

void set_difference (int set1, int set2, int set3)
{
	int set_type;
	struct set *s_ptr1, *sc_ptr1, *sc_ptr2;

	strcpy (function_name, "set_difference");

	if (get_set_type (set1) != get_set_type (set2)) return;
	s_ptr1 = set_ptr (set1);
	if (s_ptr1 == (struct set *) NULL) return;
	set_type = s_ptr1 -> type;
		init_scratch (set_type);


	sc_ptr1 = set_ptr (scratch1);
	if (sc_ptr1 == (struct set *) NULL) return;
	sc_ptr1 -> type = set_type;
	clear_set (scratch1);

	sc_ptr2 = set_ptr (scratch2);
	if (sc_ptr2 == (struct set *) NULL) return;
	sc_ptr2 -> type = set_type;
	clear_set (scratch2);

	/* utilize pre-existing routines ---
	(none of which use scratch sets) THIS IS IMPORTANT ! */

	set_complement (set2, scratch1);
	if (error ()) return;
	set_intersection (set1, scratch1, scratch1);
	if (error ()) return;
	set_complement (set1, scratch2);
	if (error ()) return;
	set_intersection (set2, scratch2, scratch2);
	if (error ()) return;
	set_union (scratch1, scratch2, set3);
	if (error ()) return;

	/* finished with scratch 1 and 2 */
	clear_set (scratch1);
	sc_ptr1 -> type = 0;
	clear_set (scratch2);
	sc_ptr2 -> type = 0;
	clear_for (set3);
}

/* automated for loop routines */

void clear_for (int set_number)
{
	struct set *s;
	
	s = set_ptr (set_number);
	if (s == (struct set *) NULL) return;
	s -> current_object = 0;
	s -> current_range = (struct range *) NULL;
}

int init_for (int set_number)
{
	int object_number;
	struct set *s; struct range *r;
	
	clear_for (set_number);
	s = set_ptr (set_number);
	if (s == (struct set *) NULL) return (0);
	s -> current_range = s -> first;
	r = s -> current_range;
	if (r == (struct range *) NULL) return (0);
	s -> current_object = r -> start;
	object_number = s -> current_object;
	return (object_number);
}

int next_for (int set_number)
{
	int object_number;
	struct set *s; struct range *r;
	
	s = set_ptr (set_number);
	if (s == (struct set *) NULL) return (0);
	r = s -> current_range;
	if (r == (struct range *) NULL) return (0);
	object_number = s -> current_object;
	object_number++;
	if (object_number > r -> end) {
		r = r -> next;
		if (r == (struct range *) NULL) {
			s -> current_range = NULL;
			s -> current_object = 0;
			return (0);
		}
		s -> current_range = r;
		object_number = r -> start;
	}
	s -> current_object = object_number;
	return (object_number);
}


/* replace end of line character with null byte */

int null_eol (char *s)
{
	char *a;
	char eol;

	eol = '\n';

	for (a = s; *a != (char) 0; a++)
		if (*a == eol) {
			*a = (char) 0;
			return (TRUE);
		}
	return (FALSE);
}

/* replace first occurence of given character with null byte */

int null_char (char *s, char c)
{
	char *a;

	for (a = s; *a != (char) 0; a++)
		if (*a == c) {
			*a = (char) 0;
			return (TRUE);
		}
	return (FALSE);
}

long symbol_size ()
{
	return (sizeof (struct symbol));
}

struct symbol *symbol_ptr (int number)
{
	struct symbol *sym_ptr;

	sym_ptr = (struct symbol *) generic_ptr (SYMBOL, number);
	return (sym_ptr);
}

/* return object number of specified identifier string */
/* it is assumed that the caller knows the object type */

int object_number_of (char *str)
{
	int object_number, symbol_number, nw;
	char junk[MAX_STRING];
	struct symbol *sym;

	if (strlen (str) == 0) return (0);
	nw = sscanf (str, "%s", junk);
	if (nw <= 0) return (0);
	symbol_number = lookup_symbol (str);
	if (symbol_number == 0) return (0);
	sym = symbol_ptr (symbol_number);
	if (sym == (struct symbol *) NULL) return (0);
	object_number = sym -> number;
	return (object_number);
}

int new_symbol (char *str, int token_type, int sub_type)
{
	int okay, symbol_number;
	struct symbol *sym;

	strcpy (function_name, "new_symbol");
	symbol_number = lookup_symbol (str);
	if (error()) return (0);
	strcpy (function_name, "new_symbol");
	okay = FALSE;
	if (symbol_number != 0) {
		sym = symbol_ptr (symbol_number);
		if (sym == (struct symbol *) NULL) return (0);
		if (sym -> global) okay = TRUE;
		if (!sym -> global && !sym -> dummy_arg)
			return (symbol_number);
		if (!okay) {
			set_error1("DUPLICATE_SYMBOL");
			return (0);
		}
	}

	symbol_number = allocate_psl (SYMBOL);
	sym = symbol_ptr (symbol_number);
	if (sym == (struct symbol *) NULL) return (0);

	/* set up fields */

	strcpy (sym -> name, str);
	sym -> token_type = token_type;
	sym -> sub_type = sub_type;
	sym -> number = 0;
	sym -> global = TRUE;
	sym -> owner = FALSE;
	sym -> initialized = FALSE;
	sym -> dummy_arg = FALSE;
	return (symbol_number);
}

int new_object_symbol (char *str, int object_type)
{
	int symbol_number;

	symbol_number = new_symbol (str, 1, object_type);
	return (symbol_number);
}

int lookup_symbol (char *str)
{
	int number, n_symbol;
	struct symbol *sym;

	n_symbol = number_of_symbols ();

	/* look for global symbol */
	for (number = 1; number <= n_symbol; number++) {
		sym = symbol_ptr (number);
		if (sym == (struct symbol *) NULL) return (0);
		if (!sym -> global) continue;
		if (strcmp (sym -> name, str) == 0) return (number);
	}
	return (0);
}

void change_owner (int object_type, int object_number, int new_owner)
{
	int symbol_number, found, same_type;
	struct symbol *sym, *owner_sym; struct generic *go;

	strcpy (function_name, "change_owner");
	same_type = FALSE;
	if (new_owner > 0) {
		owner_sym = symbol_ptr (new_owner);
		if (owner_sym == (struct symbol *) NULL) return;
		same_type = (owner_sym -> sub_type == object_type);
	}
	else owner_sym = (struct symbol *) NULL;

	go = named_ptr (object_type, object_number);
	if (go == (struct generic *) NULL) return;
	strcpy (function_name, "change_owner");
	found = 0;
	sym = NULL;
	for (symbol_number = go -> symbol_number; symbol_number != 0;
		symbol_number = sym -> next_symbol) {
		sym = symbol_ptr (symbol_number);
		if (sym == (struct symbol *) NULL) return;
		if (sym == owner_sym) {
			if (!same_type) {
				set_error1("INCONSISTENT_TYPE");
				return;
			}
			found = 1;
			sym -> owner = TRUE;
		}
		else sym -> owner = FALSE;
	}
	if (owner_sym == (struct symbol *) NULL) return;
	if (!same_type) return;
	if (!found) {
		set_error1("MISSING_OWNER");
		return;
	}
}

void link_and_set (int object_type, int object_number,int new_owner)
{
	int same_type;
	struct symbol *owner_sym; struct generic *go;
	char message[MAXLINE];

	strcpy (function_name, "link_and_set");
	owner_sym = symbol_ptr (new_owner);
	if (owner_sym == (struct symbol *) NULL) return;
	same_type = (owner_sym -> sub_type == object_type);
	if (!same_type) {
		set_error1("INCONSISTENT_TYPE");
		current_object_type = object_type;
		current_object_number = object_number;
		sprintf (message, "types = %d %d",
			owner_sym -> sub_type, object_type);
				set_error2(message);
		return;
	}
	go = named_ptr (object_type, object_number);
	if (go == (struct generic *) NULL) return;
	strcpy (function_name, "link_and_set");
	owner_sym -> next_symbol = go -> symbol_number;
	go -> symbol_number = new_owner;
	owner_sym -> number = object_number;
	owner_sym -> initialized = TRUE;
}


int allocate_named (int object_type, int symbol_number)
{
	int object_number;
	struct symbol *sym_ptr;
	struct set *s_ptr;

	strcpy (function_name, "allocate_named");

	sym_ptr = symbol_ptr (symbol_number);
	if (sym_ptr == (struct symbol *) NULL) return (0);
	if (sym_ptr -> initialized) {
		set_error1("NO_CLOBBER");
		return (0);
	}
	object_number = allocate_psl (object_type);
	if (object_number == 0) return (0);
	link_and_set (object_type, object_number, symbol_number);
	if (error()) return (0);
	change_owner (object_type, object_number, symbol_number);
	if (error()) return (0);
	if (object_type == SET) {
		s_ptr = set_ptr (object_number);
		if (s_ptr == (struct set *) NULL) return(0);
		s_ptr -> type = sym_ptr -> sub_sub_type;
	}
	return (object_number);
}

void define_real (char *str, double val)
{
	int  symbol_number, object_number;
	struct real *rea_ptr;

	symbol_number = new_symbol (str, 1, REAL);
	if (error()) return;
	object_number = allocate_named (REAL, symbol_number);
	if (object_number == 0) return;
	rea_ptr = real_ptr (object_number);
	if (error()) return;
	rea_ptr -> value = val;
}

