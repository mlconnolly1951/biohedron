/* Molecular Surface Package Copyright 1996 by Michael L. Connolly */
/* November 22, 2001 */

/* global variables */

extern int shorty;

extern struct material_table *table;

/* GLOBAL VARIABLES */

extern int current_object_type;
extern int current_object_number;
extern int current_symbol;
extern int n_ranges;
extern unsigned long scratch1, scratch2;
extern double pi_constant;
extern struct object_header *object_table;
extern struct element_block *first_element_block;
extern char function_name[MAX_NAME];

/* error handling */
extern int error_flag;
extern int fatal_flag;
extern long debug;
extern char error_string[MAXLINE];
extern char error_message[MAXLINE];
extern FILE *fperror;
extern FILE *fpinform;
extern FILE *fpdebug;
extern struct cept *head_cept;
extern struct cept *tail_cept;


/* memory */
extern long max_type;
extern long bytes_per_short;
extern long lword_size;
extern long mask_size;
extern long block_size;
extern long extension_size;
extern unsigned long current_memory;
extern unsigned long maximum_memory;
extern struct class *classes;
extern long memory_count[MEMORY_TYPES];
extern long variable_count[VARIABLE_TYPES];

extern int n_atom;				/* total number of atoms */
extern int n_file;				/* number of files read */
extern int n_molecule;			/* number of molecules read */
extern int atoms;				/* all atoms - set number */
extern int these;				/* last file read */
extern int tubers;				/* set of atoms using tube_radius */
extern int unknowns;			/* set of atoms with unknown vdw type */
extern int scratch; 			/* temp set */
extern int bonds;				/* set of bonds */
extern int n_radii;
extern int n_pattern;
extern long n_residue_index;
extern struct residueIndex *residue_table;
extern struct patternRecord *pattern_table;
extern struct radiusRecord *radius_table;

/* global variables corresponding to dot surface descriptor record: */

/* input fields */
extern int     maxatom;
extern int     natom;   
extern double  *atmco;    
extern double  *atmrad;    
extern double  *atmden;
extern short  *atmatt;
extern double   pradius;     
extern int     connected;    
extern int     errflg;  
extern char   *errstr;        
extern char    stage[80];        
extern int     atom;         
extern int     verbose;            
extern int     nsrfpnt;     
extern double   cvxarea;    
extern double   renarea;               
extern struct cluster **cvxsp;
extern struct cluster **rensp;
extern struct hidden  *intern;

/* global variables corresponding to hidden record: */

/* counters: */
extern int     ntori;
extern int     nprobe;
extern int     nlow;
extern int     ndead;
extern struct torus   *hedtor;
extern struct torus   *taltor;
extern struct probe   *hedprb;
extern struct probe   *talprb;
/* pointers to arrays of pointers */
extern struct probe  **lowp;
extern struct run **hedrun;
extern struct cluster **hedcav;
/* for connected rolling */
extern struct plist  **hedpl;
extern int     nplist;
extern short *problem;
/* for octree algorithm */
extern double   minwid;
extern struct box *rootbox;

/* other global variables: */

extern double    density;      

extern struct errmsg errtab[];
extern int ntab;

