#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"

/* Copyright 2001 by Michael L. Connolly */
/* Last revised: January 4, 2002 */


/* mscept.c */

int error ()
{
	return (error_flag);
}

void set_fatal ()
{
	fatal_flag = 1;
}

void set_error1 (char *str)
{
	int len;
	
	len = strlen (str);
	if (len > MAXLINE - 1) {
		strcpy (error_string, "unknown error");
	}
	else strcpy (error_string, str);
	error_flag = 1;
}

void set_error2 (char *str)
{
	int len;
	
	len = strlen (str);
	if (len > MAXLINE - 1) {
		strcpy (error_message, "unknown error");
	}
	else strcpy (error_message, str);
	error_flag = 1;
}

void print_error ()
{
    if (head_cept == NULL)
        add_old();
    print_cepts ();
}

void inform (char *str)
{
	fprintf(fpinform,"%s\n",str);	
}

void informd (char *str)
{
	if (debug) fprintf(fpdebug,"%s\n",str);	
}

void informd2 (char *str)
{
	if (debug > 1) fprintf(fpdebug,"%s\n",str);	
}

struct cept *allocate_cept ()
{
    struct cept *ex;

	ex = (struct cept *) allocate_object (CEPT);
	if (ex == NULL) {
		return(NULL);
	}
	return (ex);
}

struct cept *new_cept (int type, int subtype, int severity)
{
	struct cept *ex;
	
	ex = allocate_cept ();
	if (ex == NULL) {
		fprintf (stderr, "new_cept: no memory, bailing out!\n");
		exit (1);
	}
	ex -> type = type;
	ex -> subtype = subtype;
	ex -> severity = severity;
	ex -> object_type = 0;
    strcpy (ex -> remedy, "email connolly@biohedron.com this entire output");
	ex -> n_function = 0;
	ex -> n_message = 0;
	ex -> n_double = 0;
	ex -> n_long = 0;
	if (head_cept == NULL) {
		head_cept = ex;
	}
	else {
		tail_cept -> next = ex;
	}
	tail_cept = ex;
    if (severity >= FATAL_SEVERITY) error_flag = 1;
	return (ex);
}

void add_object (struct cept *ex, int object_type, char *object_name) {
    if (ex == NULL) return;
    ex -> object_type = object_type;
    strncpy (ex -> object, object_name, MAX_NAME-1);
}

void add_variable (struct cept *ex, int variable_type, char *variable_name) {
    if (ex == NULL) return;
    ex -> variable_type = variable_type;
    strncpy (ex -> variable, variable_name, MAX_NAME-1);
}

void add_source (struct cept *ex, char *source) {
    if (ex == NULL) return;
    strncpy (ex -> source, source, MAX_NAME-1);
}

void add_function (struct cept *ex, char *function) {
    if (ex == NULL) return;
	if (ex -> n_function >= MAX_EXCEPTION - 1) return;
	strncpy (ex -> functions[ex -> n_function], function, MAX_NAME-1);
	ex -> n_function++;
}

void add_message (struct cept *ex, char *message) {
    if (ex == NULL) return;
	if (ex -> n_message >= MAX_EXCEPTION - 1) return;
	strncpy (ex -> messages[ex -> n_message], message, MAX_NAME-1);
	ex -> n_message++;
}

void add_remedy (struct cept *ex, char *remedy) {
    if (ex == NULL) return;
	strncpy (ex -> remedy, remedy, MAX_NAME-1);
}

void add_double (struct cept *ex, char *double_name, double d) {
    if (ex == NULL) return;
	if (ex -> n_double >= MAX_EXCEPTION - 1) return;
	strncpy (ex -> double_names[ex -> n_double], double_name, MAX_NAME-1);
	ex -> doubles[ex -> n_double] = d;
	ex -> n_double++;
}

void add_long (struct cept *ex, char *long_name, long l) {
    if (ex == NULL) return;
	if (ex -> n_long >= MAX_EXCEPTION - 1) return;
	strncpy (ex -> long_names[ex -> n_long], long_name, MAX_NAME-1);
	ex -> longs[ex -> n_long] = l;
	ex -> n_long++;
}

void add_atom (struct cept *ex, struct sphere *a) {
    if (ex == NULL) return;
	if (ex -> n_atom >= MAX_EXCEPTION - 1) return;
	ex -> atoms[ex -> n_atom] = a;
	ex -> n_atom++;
}

void add_old () {
	struct cept *ex;
	ex = new_cept (0, 0, fatal_flag);
	add_message (ex, error_string);
	add_message (ex, error_message);
}

void print_cepts () {
	struct cept *ex;
    for (ex = head_cept; ex != NULL; ex = ex -> next) {
		print_cept (ex);
    }
}

void print_cept (struct cept *ex) {
    int i;
	char *severity_name;
	char *type_name;
	char *subtype_name;
    struct sphere *atm_ptr;
    char object_type_name[MAX_NAME];

    if (ex == NULL) {
        fprintf (stderr, "unknown error\n");
        exit (1);
    }
	if (ex -> type != 0) {
		severity_name = severity_to_name (ex -> severity);
		type_name = type_to_name (ex -> type);
		subtype_name = subtype_to_name (ex -> subtype);
		fprintf (stderr, "%s %s error: %s\n", 
			severity_name, type_name, subtype_name);
	}
	if (strlen (ex -> source) > 0)
		fprintf (stderr, "C source file: %s\n", ex -> source);
	for (i = 0; i < ex -> n_function; i++) {
        if (i == 0)
			fprintf (stderr, "in function: %s\n", ex -> functions[i]);
        else
			fprintf (stderr, "  called by: %s\n", ex -> functions[i]);
	}
    if (ex -> object_type != 0) {
		get_class_name(ex -> object_type, object_type_name);
	    fprintf (stderr, "object type = %s\n", object_type_name);
    }
    if (ex -> variable_type != 0) {
		fprintf (stderr, "variable = %s\n", 
			variable_to_name (ex -> variable_type));
    }
    if (strlen(ex -> object) > 0) {
	    fprintf (stderr, "object name = %s\n", ex -> object);
    }
    if (strlen(ex -> array_name) > 0) {
	    fprintf (stderr, "array name = %s\n", ex -> array_name);
    }
    if (strlen(ex -> input_name) > 0) {
	    fprintf (stderr, "input file = %s\n", ex -> input_name);
    }
    for (i = 0; i < ex -> n_message; i++) {
    	fprintf (stderr, "%s\n", ex -> messages[i]);
    }
    for (i = 0; i < ex -> n_double; i++) {
    	fprintf (stderr, "%s = %12.6f\n", 
    	ex -> double_names[i], ex -> doubles[i]);
    }
    for (i = 0; i < ex -> n_long; i++) {
    	fprintf (stderr, "%s = %12ld\n", 
    	ex -> long_names[i], ex -> longs[i]);
    }
    for (i = 0; i < ex -> n_atom; i++) {
        atm_ptr = ex -> atoms[i];
        if (atm_ptr == NULL) {
            fprintf (stderr, "null pointer for atom\n");
            continue;
        }
        fprintf (stderr, "%5s %5s %-5s %9.3f %9.3f %9.3f %7.3f %5ld\n",
			atm_ptr -> group, atm_ptr -> sequence, atm_ptr -> name,
			atm_ptr -> center[0], atm_ptr -> center[1], atm_ptr -> center[2],
			atm_ptr -> radius, atm_ptr -> number);
    }
    fprintf (stderr, "remedy: %s\n", ex -> remedy);
}

char *type_to_name (int type) {
    char *name;

	name = allocate_chars (MAX_NAME+1);
	if (name == NULL) {
		fprintf (stderr, "type_to_name: memory failure");
		exit (1);
	}

    switch (type) {
    case UNKNOWN_ERROR:
        strcpy (name, "unknown");
        break;
    case GEOMETRY_ERROR:
        strcpy (name, "geometry");
        break;
    case MEMORY_ERROR:
        strcpy (name, "memory");
        break;
    case INPUT_ERROR:
        strcpy (name, "input");
        break;
    case PARAMETER_ERROR:
        strcpy (name, "parameter");
        break;
    case SYNTAX_ERROR:
        strcpy (name, "syntax");
        break;
    case ARRAY_ERROR:
        strcpy (name, "array");
        break;
    case LOGIC_ERROR:
        strcpy (name, "logic");
        break;
    case GRID_ERROR:
        strcpy (name, "grid");
        break;
    case RETURN_ERROR:
        strcpy (name, "grid");
        break;
    default:
        strcpy (name, "unknown");
    }
    return (name);
}

char *subtype_to_name (int type) {
    char *name;

	name = allocate_chars (MAX_NAME+1);
	if (name == NULL) {
		fprintf (stderr, "subtype_to_name: memory failure");
		exit (1);
	}

    switch (type) {
    case UNKNOWN_SUBTYPE:
        strcpy (name, "unknown");
        break;
    case ALLOCATION:
        strcpy (name, "allocation");
        break;
    case BOUNDS:
        strcpy (name, "bounds");
        break;
    case FREEING:
        strcpy (name, "freeing");
        break;
    case INCONSISTENCY:
        strcpy (name, "inconsistency");
        break;
    case FILE_NAME:
        strcpy (name, "file name");
        break;
    case NULL_VALUE:
        strcpy (name, "null value");
        break;
    case INVALID_VALUE:
        strcpy (name, "invalid value");
        break;
    case NEGATIVE_VALUE:
        strcpy (name, "negative value");
        break;
    case CONTENT:
        strcpy (name, "content");
        break;
    case PREMATURE_EOF:
        strcpy (name, "premature eof");
        break;
    case NULL_POINTER:
        strcpy (name, "null pointer");
        break;
    case MSOVERFLOW:
        strcpy (name, "overflow");
        break;
    case DEGENERACY:
        strcpy (name, "degeneracy");
        break;
    case MSUNDERFLOW:
        strcpy (name, "underflow");
        break;
    default:
        strcpy (name, "unknown");
    }
    return (name);
}

char *severity_to_name (int type) {
    char *name;

	name = allocate_chars (MAX_NAME+1);
	if (name == NULL) {
		fprintf (stderr, "severity_to_name: memory failure");
		exit (1);
	}

    switch (type) {
    case WARNING_SEVERITY:
        strcpy (name, "warning");
        break;
    case FATAL_SEVERITY:
        strcpy (name, "fatal");
        break;
    default:
        strcpy (name, "unknown");
    }
    return (name);
}

char *variable_to_name (int variable) {
    char *name;

	name = allocate_chars (MAX_NAME+1);
	if (name == NULL) {
		fprintf (stderr, "variable_to_name: memory failure");
		exit (1);
	}

    switch (variable) {
    case ARC_DIRECTION:
        strcpy (name, "arc_direction");
        break;
    case ATMCO:
        strcpy (name, "atmco");
        break;
    case ATMRAD:
        strcpy (name, "atmrad");
        break;
    case ATOM_ALPHAS:
        strcpy (name, "atom_alphas");
        break;
    case ATOM_COLORS:
        strcpy (name, "atom_colors");
        break;
    case ATOM_OPACITIES:
        strcpy (name, "atom_opacities");
        break;
    case ATOM_CENTERS:
        strcpy (name, "atom_centers");
        break;
    case ATOM_RADII:
        strcpy (name, "atom_radii");
        break;
    case ATMDEN:
        strcpy (name, "atmden");
        break;
    case BAD_VERTEX:
        strcpy (name, "bad_vertex");
        break;
    case BAD_NEIGHBOR:
        strcpy (name, "bad_neighbor");
        break;
    case BASES:
        strcpy (name, "bases");
        break;
    case BIT1_GRID:
        strcpy (name, "bit1_grid");
        break;
    case BIT2_GRID:
        strcpy (name, "bit2_grid");
        break;
    case BOUND_SAME:
        strcpy (name, "bound_same");
        break;
    case CENTERS:
        strcpy (name, "centers");
        break;
    case CIR:
        strcpy (name, "cir");
        break;
    case CIRCLE1:
        strcpy (name, "circle1");
        break;
    case CIRCLE2:
        strcpy (name, "circle2");
        break;
    case CIRCLE3:
        strcpy (name, "circle3");
        break;
    case CIRCLE_PTR:
        strcpy (name, "circle_ptr");
        break;
    case CONE_CIRCLE:
        strcpy (name, "cone_circle");
        break;
    case CONTAINS:
        strcpy (name, "contains");
        break;
    case CNUMS:
        strcpy (name, "cnums");
        break;
    case CYCLE_USED:
        strcpy (name, "cycle_used");
        break;
    case DISTANCE_BACKWARD:
        strcpy (name, "distance_backward");
        break;
    case DISTANCE_FORWARD:
        strcpy (name, "distance_forward");
        break;
    case EDGES:
        strcpy (name, "edges");
        break;
    case EDGVTX:
        strcpy (name, "edgvtx");
        break;
    case ENUMBERS0:
        strcpy (name, "enumbers0");
        break;
    case ENUMBERS1:
        strcpy (name, "enumbers1");
        break;
    case EXTERIOR_ANGLE:
        strcpy (name, "exterior_angle");
        break;
    case F:
        strcpy (name, "f");
        break;
    case F1:
        strcpy (name, "f1");
        break;
    case F2:
        strcpy (name, "f2");
        break;
    case F3:
        strcpy (name, "f3");
        break;
    case FIRST_CYCLE_NUMBER:
        strcpy (name, "first_cycle_number");
        break;
    case FIRST_EDGE_NUMBER:
        strcpy (name, "first_edge_number");
        break;
    case FN1:
        strcpy (name, "fn1");
        break;
    case FN2:
        strcpy (name, "fn2");
        break;
    case FN3:
        strcpy (name, "fn3");
        break;
    case FOCI:
        strcpy (name, "foci");
        break;
    case HEMI_CIRCLE:
        strcpy (name, "hemi_circle");
        break;
    case INTEGERS:
        strcpy (name, "integers");
        break;
    case LFCIR:
        strcpy (name, "lfcir");
        break;
    case MEDIAN:
        strcpy (name, "median");
        break;
    case NORMALSV:
        strcpy (name, "normalsv");
        break;
    case NEIGHBORS:
        strcpy (name, "neighbors");
        break;
    case NFC:
        strcpy (name, "nfc");
        break;
    case NUMBERS:
        strcpy (name, "numbers");
        break;
    case PERIMETER:
        strcpy (name, "perimeter");
        break;
    case POINTSV:
        strcpy (name, "pointsv");
        break;
    case POLYGON_SIDE:
        strcpy (name, "polygon_side");
        break;
    case POLYGON_TANGENT:
        strcpy (name, "polygon_tangent");
        break;
    case PRBDOT:
        strcpy (name, "prbdot");
        break;
    case PROJECTED_VERTEX_LIST:
        strcpy (name, "projected_vertex_list");
        break;
    case PROJECTED_VERTICES:
        strcpy (name, "projected_vertices");
        break;
    case RADII:
        strcpy (name, "radii");
        break;
    case REALS:
        strcpy (name, "reals");
        break;
    case TORCIR:
        strcpy (name, "torcir");
        break;
    case TRIANGLES:
        strcpy (name, "triangles");
        break;
    case TRIEDGVTX:
        strcpy (name, "triedgvtx");
        break;
    case USED:
        strcpy (name, "used");
        break;
    case VERTEX_CENTERS:
        strcpy (name, "vertex_centers");
        break;
    case VERTEX_NORMALS:
        strcpy (name, "vertex_normals");
        break;
    case VERTICES:
        strcpy (name, "vertices");
        break;
    case VERTS:
        strcpy (name, "verts");
        break;
    case VNUMBERS0:
        strcpy (name, "vnumbers0");
        break;
    case VNUMBERS1:
        strcpy (name, "vnumbers1");
        break;
    case XGRID:
        strcpy (name, "xgrid");
        break;
    case YGRID:
        strcpy (name, "ygrid");
        break;
    case ZGRID:
        strcpy (name, "zgrid");
        break;
    default:
        strcpy (name, "unknown");
    }
    return (name);
}
