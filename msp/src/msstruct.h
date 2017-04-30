/* Molecular Surface Package Copyright 1996 by Michael L. Connolly */
/* last revised January 8, 2002 */

typedef unsigned short   atomnum;

/* command-line arguments */

struct argrec {
	short flag;
	short on;
	long position;
	long longarg;
	double doublearg;
	char stringarg[64];
};

/* memory structs */


struct class {
	long type;
	unsigned long size;	/* in bytes */
	long count;
	long pointer_count;
	char name[MAX_TYPE_NAME];
	void *cache[MAX_CACHE];
};

struct premem {
	unsigned char memtype;
	unsigned char objtype;
	unsigned short objsize;
	unsigned long count;
};

struct component {
	double volume;
	double area;
	double accessible;
	double center[3];
	struct chunk *head_chunk;
	long count;
	int type;
	double pvolume;
	double parea;
};

struct chunk {
	struct chunk *next[2];
	atomnum atom_number;
	char labels[3][MAX_ATNAME]; /* 3 strings, up to 7 characters each */
	atomnum component_number;
	double contact_area;
	double reentrant_area;
	double accessible_area;
};

struct evalpnt {
	double center[3];
	double radius;
	double expansion;
	double omegas[3];
	long vertex_number;
	struct phnvtx *vtx;
};


struct polygon {
	struct polygon *next;
	short hue;
	atomnum comp;
	atomnum atm;
	struct edge *first_edge;
	struct phnedg *edg[MAXGON];
	long vertex_index[MAXGON];
	long edge_number[MAXGON];
	short orn[MAXGON];
	short clipped;
	short n_side;
	long material_index;
	double center[3];
	double axis[3];
	double area;
	double radius;
	double vc[MAXGON][3];
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	double zmin;
	double zmax;
};

/* non-axial intersection structure */

struct nai {
	struct face *fac;
	struct edge *edg1;
	struct edge *edg2;
	struct arc *arc1;
	struct arc *arc2;
	struct circle *cir1;
	struct circle *cir2;
	struct cusp *csp1;
	struct cusp *csp2;
	struct probe *prb;
	struct probe *prb1;
	struct probe *prb2;
	short nx;
	short eaten;
	double xpnt[3];
	double closest_distance;
};


struct surface {
	long blank;
	long centering;
	long converted;
	long clipping;
	long solid_shade;
	long do_cusp_intersections;
	long format;
	long n_arc;
	long n_atom;
	long n_bond;
	long n_central;
	long n_circle;
	long n_component;
	long n_cycle;
	long n_edge;
	long n_face;
	long n_surface;
	long n_nai;
	long n_probe;
	long n_problem_atom;
	long n_problem_face;
	long n_phnvtx;
	long o_phnvtx;
	long n_phnedg;
	long o_phnedg;
	long n_phntri;
	long n_phnctr;
	long n_tori;
	long n_variety;
	long n_vertex;
	long n_evaluation;
	long n_prbdot;
	long n_polygon;
	long o_polygon;
	long n_circle_eaten;	/* number of cusp circles eaten */
	long n_cusp_circle;		/* number of cusp circles */
	long n_cusp_arc;		/* number of cusp circles */
	long n_cusp_resort;		/* number of cusp circle resorts */
	long n_eaten_arc;		/* number of arcs eaten by axial X's */
	long n_low_probe;		/* number of low probes */
	long n_non_axial;		/* number of non-axial intersections */
	long n_nontrimmed;		/* number of cusp intersections not trimmed */
	long n_pair;
	long n_point_torus;		/* number of point cusp tori */
	long n_point_cusp_vertices;
	long n_vertex_eaten;	/* number of cusp vertices eaten */
	long old_volume_computed;
	long n_arc_original;				/* input number of arcs */
	long n_face_original;			/* input number of faces */
	long n_large_face;
	long n_entire_face;
	long n_simplified_face;
	long polyhedron_format;
	long surface_completed;
	long surface_type;		/* PQMS or BAS */
	long triangulation_completed;
	long type;
	long use_grid;
	long van_der_Waals;
	long volume_computed;
	long radius_specified;
	long net_setup;
	long polyhedron_read;
	long evaluation_completed;
	char name[64];
	char function[MAXLINE];
	double eval_radius;
	double linewidth[4];
	double alpha;
	double bit_width;
	double cavity_alpha;
	double cube_width;
	double global_alpha;			/* maximum arc angle parameter */
	double large_omega;			/* large omega parameter */
	double max;
	double min;
	double mol_vol;
	double molecule_area;
	double old_volume;
	double polyhedron_volume;
	double probe_radius;
	double surface_volume;
	double surface_thickness;
	double total_accessible_area;
	double total_area;
	double total_contact_area;
	double total_reentrant_area;
	double total_volume;
	double weight;				/* weight for select function */
	double bounds[2][3];
	double center[3];
	double centroid[3];
	double maxvals[3];
	double minvals[3];
	double rotation[3][3];
	double shape_area[3];
	double translate[3];
	double *prbdot;
	double *atom_centers;
	double *atom_radii;
	struct surface *next;
	struct sphere *head_atom;
	struct sphere *tail_atom;
	struct central *head_central;
	struct central *tail_central;
	struct circle *head_circle;
	struct circle *tail_circle;
	struct cusp *head_cusp;
	struct cusp *tail_cusp;
	struct face *head_face;
	struct face *tail_face;
	struct probe *head_probe;
	struct probe *tail_probe;
	struct torus *head_torus;
	struct torus *tail_torus;
	struct vertex *head_vertex;
	struct vertex *tail_vertex;
	struct arc *head_arc;
	struct arc *tail_arc;
	struct phnvtx *head_phnvtx;
	struct phnvtx *tail_phnvtx;
	struct phnedg *head_phnedg;
	struct phnedg *tail_phnedg;
	struct phntri **heads;
	struct phntri **tails;
	struct polygon *head_polygon;
	struct polygon *tail_polygon;
	struct phnctr *head_phnctr;
	struct phnctr *tail_phnctr;
	struct quartet *head_quartet;
	struct quartet *tail_quartet;
	struct probe *head_delinked;
	struct probe *tail_delinked;
	struct object_scheme *scheme;
	struct variety **variety_handles;
	struct arc **arc_handles;
	struct face **face_handles;
	struct component **component_handles;
	struct vertex **vertex_handles;
	struct circle **circle_handles;
	struct cycle **cycle_handles;
	struct edge **edge_handles;
	struct pair *pair_array;
	struct probe **low_probe_hdl;
	struct torus **point_torus_hdl;
	struct spheregrid *sg;
	struct molecule *mol;
	struct linseg *head_linseg;
	struct linseg *tail_linseg;

	double *vertex_centers;
	double *vertex_normals;
	long *edgvtx;
	long *triedgvtx;
	struct phnvtx **phnvtx_handles;
	struct phnvtx **original_phnvtxs;
	struct phnedg **phnedg_handles;
	struct phnedg **original_phnedgs;
	struct phntri **phntri_handles;
	struct polygon **original_polygons;
	struct polygon **polygon_handles;
	struct evalpnt *evaluation_spheres;

	long origin[3];
	long width[3];
	long limit[3];
	long n_cube;
	float *densities;

	long n_edger;
	struct edger *root_edger;
	struct nai nais[MAX_NAI];
	struct depth_buffer *db;
};


struct msscene {
	long error_flag;
	long fatal_flag;
	long purge_frequency;	/* if 0, no purging */
	long plot_format;
	long raster_format;
	long vector_format;
	long clipping;
	long translucency;
	long horizontal;		/* horizonal dimension in pixels */
	long vertical;			/* vertical dimension in pixels */
	long zview;				/* z dimension depth units */
	long max_type;
	long debug;
	long stereo;			/* left or right or mono */
	long screen_size;
	long n_row;
	long n_column;
	long interpolate;
	long vrml;
	long screendoor;
	long overlap_hue;
	long n_square;
	long current_npolygon;
	long viewport[2][3];
	long n_bdy_parity, n_clip_parity;
	long n_bad_projection;
	long n_missing_points;
	long n_missing_leaves;
	unsigned long size;		/* number of pixels in frame */
	double fineness;
	double coalesce;
	double alignment;
	double zrange;
	double clip_fraction;
	double clip_center[3];
	double clip_axis[3];
	double pixel_width;		/* width of pixel in angstroms */
	double plot_width;
	double center[3];
	double rotation[3][3];
	double stereo_angle;
	double stereomat[3][3];
	double window[2][3];
	double lsource[3];		/* light source */
	double model[5];		/* depth, ambient, diffuse, specular exponent */
	double bounds[2][3];
	double psfactor;
	double border;
	double tiltx;
	double tilty;
	double red_amount;
	double green_amount;
	double blue_amount;
	double current_linewidth;
	char error_string[MAXLINE];
	char error_message[MAXLINE];
	char title[80];
	struct depth_buffer *db;
	FILE *fp_plot, *fp_plot2;
	FILE *fp_raster, *fp_raster2;
	FILE *fp_vector;
	FILE *fperror;
	FILE *fpinform;
	FILE *fpdebug;
	struct polygon **current_polygons;
	struct square *squares;
	struct square *current_square;
	struct class *classes;
	struct molecule *current_molecule;
	struct molecule *head_molecule;
	struct molecule *tail_molecule;
	struct molecule *this_mol;
	struct surface *this_srf;
	struct material_table *table;
	struct color_ramp *head_ramp;
	struct color_ramp *tail_ramp;
};


struct msdata {
	long narray;
	long types[MaxMSData];
	long widths[MaxMSData];
	long counts[MaxMSData];
	long *longs[MaxMSData];
	float *floats[MaxMSData];
	unsigned char *bytes[MaxMSData];
	int routine;
	int variable;
};

struct array {
	long length;
	long type;
	long *integers;
	double *reals;
	double *radii;
	double *centers;
	char **strings;
	struct atom **atoms;
	struct bond **bonds;
};

struct vanity {
	float scalar;
	float vector[3];
};

struct denvtx {
	long number;
	double normal[3];
	double average;
};


/* this structure shows up in dsall, dscvx and dscyc */
struct tlist {
    struct tlist   *next;       /* next torus entry in an atom's list */
    struct torus   *tor;        /* pointer to torus */
};

/* hashing table probe lists for connected trajectory mode */
struct plist {
    struct plist   *next;       /* next in linked list for this hash key */
    struct probe   *prb;        /* pointer to probe */
};


/* endpoint of arc on circle */
struct endpnt {
    double   angle;              /* angle of endpoint relative to base */
    struct endpnt  *next;       /* next endpoint going counter-clockwise */
    atomnum atom;               /* atom probe collided with to make endpnt */
    short   begin;              /* true if beginning of arc */
};


/* run length encoding is used for lists of atom numbers,
   because the atoms are generally near each other in space,
   and so are likely to have contiguous atom numbers */

struct run {
    struct run *next;           /* pointer to next contiguous run */
    atomnum atom[2];            /* beginning and ending atom numbers */
};

/* the box is a cube and is used for octree algorithm for finding
   the neighbors of each atom at the beginning of the calculation */

struct box {
    double   bounds[3][2];       /* bounds of box */
    struct box *child[8];       /* eight boxes inside this one (if not leaf) */
    short   leaf;               /* flag saying that there are no subnodes */
    atomnum * atom;             /* pointer to list of atoms in box (if leaf) */
};

/* new structure for determining number of cycles */
struct edg {
    struct endpnt  *ept1;
    struct endpnt  *ept2;
    struct edg *next;
    int     atom;
    double   ept1co[3];
    double   ept2co[3];
};

struct midpln {             /* mid plane between low probes */
    double   center[3];
    double   axis[3];
};

/* internal record - hidden from calling program */

struct hidden {
    struct dsdesc  *dsdptr;     /* pointer back to surface descriptor */
 /* counters: */
    int     ntori;
    int     nprobe;
    int     nlow;
    int     ndead;
 /* pointers to linked lists */
    struct torus   *hedtor;
    struct torus   *taltor;
    struct probe   *hedprb;
    struct probe   *talprb;
 /* pointers to arrays of pointers */
    struct probe  **lowp;
    struct run **hedrun;
    struct cluster **hedcav;
 /* for connected rolling */
    struct plist  **hedpl;
    int     nplist;
    short *problem;
 /* for octree algorithm */
    double   minwid;
    struct box *rootbox;
};


struct errmsg {
	int     number;
	char   *msg;
};
struct leaf {
	long atmnum[MAXPA];
	short comp;
	short shape;
	short type;
	short side;
	short cep;
	short clip_ep;
	long input_hue;
	struct face *fac;
	struct circle *cir;
	double ends[2][3];
	short where[2];
	short when[2];
	double focus[3];
};


struct lax {
	double co[3];
	short ent;
	short used;
	double angle;
};

/* SGI image format: 512 bytes total */

struct sgi_header {
    unsigned short	magic;	 	/* 0732 */
    unsigned char 	storage;	/* 0 (not rle) */
    unsigned char 	bpc;		/* 1 */
    unsigned short 	dimension;	/* 3 */
    unsigned short 	xsize;		/* horizontal */
    unsigned short 	ysize;		/* vertical */
    unsigned short 	zsize;		/* 3 for RGB */
    unsigned long 	min;		/* 0 - intensity */
    unsigned long 	max;		/* 255 - intensity */
    unsigned long 	dummy;
    char		name[80];	/* name of image */
    unsigned long	colormap;	/* 0 */
    unsigned char	unused[404];
};


/* old format:
struct sgi_header {
    unsigned short	imagic;
    unsigned short 	type;
    unsigned short 	dim;
    unsigned short 	xsize;
    unsigned short 	ysize;
    unsigned short 	zsize;
    unsigned long 	min;
    unsigned long 	max;
    unsigned long	wastebytes;
    char			name[80];
    unsigned long	colormap;
};
*/

/* SUN Microsystem's format: */
struct sun_header {
	long magic;	/*  0x59a66a95 */
	long width;
	long height;
	long depth;
	long length;
	long type;
	long maptype;
	long maplength;
};

/* AVS image format */
struct avs_header {
	long width;
	long height;
};


struct pixel_24 {
	unsigned char ucr;
	unsigned char ucg;
	unsigned char ucb;
};

struct pixel_32 {
	unsigned char ucr;
	unsigned char ucg;
	unsigned char ucb;
	unsigned char uct;
};


/* hidden-line elimination structures */

struct linseg {
	struct linseg *next;
	double ends[2][2];
	double height[2];
	double value;
	short hidden;		/* 1 if hidden */
	atomnum comp;
	atomnum atm;
	short hue;
	short inner;
	short filler;
	struct phnedg *edg;
};
	
struct square {
	double xmin;
	double xmax;
	double ymin;
	double ymax;
	double zmin;
	double zmax;
	long npolygon;
	struct polygon **polygons;
};


struct axis {
        int copy_number;
        int symbol_number;
        double center[3];
        double direction[3];
        struct element *first;
};

struct point {
        int copy_number;
        int symbol_number;
        double coordinates[3];
        struct element *first;
};

struct transformation {
        int copy_number;
        int symbol_number;
        double scaling;
        double center[3];
        double rotation[3][3];
        double translation[3];
};

/* for pqms grid */

struct scube {
	short occupancy;
	short reserved;
	struct snumber *first;
	struct snumber *last;
};

struct snumber {
	long number;
	struct snumber *next;
};

struct sNumberBlock {
	struct snumber *snumbers;
	struct sNumberBlock *next;
};

struct spheregrid {
	/* joint variables (both bitmap and cube) */
	long nspheres;
	double *centers;
	double *radii;
	double bounds[2][3];
	double origin[3];
	/* cube grid variables */
	double cubewidth;
	long n_cube[3];
	long ncubexyz;
	long nsnumbers;
	long usnumbers;
	struct scube *scubes;
	struct snumber *freesn;
	struct sNumberBlock *first_block;
	short *union_numbers;
	/* bitmap variables */
	double bitwidth;
	long nbit[3];
	long nbitxyz;
	long bx;
	long by;
	long bz;
	long bxy;
	unsigned long nbitlong;
	unsigned long *bit1_grid;		/* bit on if bit in interior or on surface */
	unsigned long *bit2_grid;		/* bit on if bit in interior */
	long ninterior;
	long nsurface;
	long nexterior;
};


struct material_table {
	long nmaterial;
	double red[MAX_MATERIAL];
	double green[MAX_MATERIAL];
	double blue[MAX_MATERIAL];
	long hue_count[MAX_MATERIAL];
	unsigned char start[MAX_MATERIAL];
	unsigned char nshades[MAX_MATERIAL];
	char name[MAX_MATERIAL][16];
	/* for pseudocolor -> red-green-blue-transparency */
	unsigned char ucrs[256];
	unsigned char ucgs[256];
	unsigned char ucbs[256];
	unsigned char ucts[256];
};


struct color_ramp {
	struct color_ramp *next;
	char name[MAX_NAME];
	long n_ramp;
	long colors[MAX_LEVELS];
	long vertex_colors[MAX_LEVELS][2];
	double vertex_ranges[MAX_LEVELS][2];
};

struct object_scheme {
	long color_type;
	long opacity_type;
	long n_atom;
	long n_component;
	long which_value;
	double minval;
	double maxval;
	long uniform_colors[2];
	long component_colors[4];
	long shape_colors[6];
	struct color_ramp *ramp;
	long *atom_colors;
	double uniform_opacities[2];
	double component_opacities[4];
	double shape_opacities[6];
	double *atom_opacities;
};

/* atom */

struct atom {
	long number;			/* number of atom in input file */
	long hue;
	short type;
	short chemical_element;
	short anumber;
	short srn;
	short mol;
	short col;
	unsigned char buried;		/* atom is completely inaccessible */
	unsigned char skip;			/* skip atom */
	unsigned char problem;		/* problem atom */
	unsigned char filler;
	float covalent;
	float ball;
	float angle;
	float density;
	float occupancy;
	float tfactor;
	float opacity;
	double ers;				/* expanded radius squared */
	double radius;
	double center[3];
	long bonded_to[MAX_VALENCE];
	long bond_number[MAX_VALENCE];
	char name[MAX_ATNAME];
	char group[MAX_ATNAME];
	char sequence[MAX_ATNAME];
	char pdb[MAX_ATNAME];
	char subunit[MAX_ATNAME];
	char kind[MAX_ATNAME];
	struct atom *next;		/* next atom in linked list */
	struct atom **first_neighbor;	/* first neighbor */
	struct atom **last_neighbor;	/* last neighbor */
	struct pair **first_torus;	/* first torus */
	struct pair **last_torus;	/* last torus */
	struct arc *first_arc;		/* first arc */
	struct face *first_face;	/* first face */
	struct chunk *head_chunk;
};

struct sphere {
	long number;			/* number of atom in input file */
	long hue;
	short type;
	short chemical_element;
	short anumber;
	short srn;
	short mol;
	short col;
	unsigned char buried;		/* atom is completely inaccessible */
	unsigned char skip;			/* skip atom */
	unsigned char problem;		/* problem atom */
	unsigned char filler;
	float covalent;
	float ball;
	float angle;
	float density;
	float occupancy;
	float tfactor;
	float opacity;
	double ers;				/* expanded radius squared */
	double radius;
	double center[3];
	long bonded_to[MAX_VALENCE];
	long bond_number[MAX_VALENCE];
	char name[MAX_ATNAME];
	char group[MAX_ATNAME];
	char sequence[MAX_ATNAME];
	char pdb[MAX_ATNAME];
	char subunit[MAX_ATNAME];
	char kind[MAX_ATNAME];
	struct sphere *next;		/* next sphere in linked list */
	struct neighbor *first_neighbor;	/* first neighbor */
	struct neighbor *last_neighbor;	/* last neighbor */
	struct pair **first_torus;	/* first torus */
	struct pair **last_torus;	/* last torus */
	struct arc *first_arc;		/* first arc */
	struct face *first_face;	/* first face */
	struct chunk *head_chunk;
};

struct patternRecord {
	char residue[MAX_ATNAME];
	char atom_name[MAX_ATNAME];
	long type;
	char kind[MAX_ATNAME];
};

struct radiusRecord {
	long type;
	double radius;
	double covalent;
	char kind[MAX_ATNAME];
};

struct residueIndex {
	char residue[MAX_ATNAME];
	long first_pattern;
	long last_pattern;
};

struct bond {
	long atoms[2];
	short hue1;
	short hue2;
	short modified1;
	short modified2;
	float radius;
	float end1[3];
	float end2[3];
};

struct bas {
	unsigned long n_atom;
	unsigned long n_bond;
	struct atom *atoms;
	struct bond *bonds;
};

struct molecule {
	int copy_number;
	int symbol_number;
	int skeleton_number;
	long n_atom;
	long n_bond;
	long n_component;
	long n_surface;
	long type;
	long blank;
	char name[MAXLINE];
	double outer_width;
	double inner_width;
	double cavity_width;
	double bond_width;
	double elbow;
	double bond_radius;
	double ball_radius;
	double tolerance;
	double alpha;
	double cavity_alpha;
	double global_alpha;			/* maximum arc angle parameter */
	double bounds[2][3];
	double center[3];
	double centroid[3];
	long atom_set;
	long bond_set;
	double *atom_alphas;			/* max subdivision angle */
	char alpha_word[64];		/* angle parameter value or filename */
	struct molecule *next;
	struct atom *head_atom;
	struct atom *tail_atom;
	struct surface *head_surface;
	struct surface *tail_surface;
	struct surface *current_surface;
	struct atom_array *atoms;
	struct bond_array *bonds;
	FILE *fp_alpha;
};

/* polyhedron structures */

struct phnvtx {
	short degree;
	short critter;
	atomnum comp;
	atomnum atm;
	unsigned short hue;
	unsigned long number;
	long material_index;
	double center[3];
	double outward[3];
	double normal[3];
	double base[3];
	double zenith[3];
	double values[3];
	struct phnvtx *next;
	struct phnedg *frsedg;
	struct phnedg *on;
	struct critlink *head;
	unsigned char *byte_density;
};

struct phnedg {
	atomnum comp;
	atomnum atm;
	unsigned short hue;
	unsigned short may_intersect;
	short type;
	short clipped;
	short used;
	short upper;
	long nvtx;
	long opposites[2];
	unsigned short enter[2];
	unsigned long number;
	unsigned long vtxnum[2];
	struct phnvtx *pvt[3];
	struct phnedg *next;
	struct phnedg *next_ctredg;
	struct phnedg *split[2];
	struct vertex *vtx[2];
	struct polygon *on[2];
};

struct phntri {
	atomnum comp;
	atomnum atm;
	unsigned short hue;
	unsigned short may_intersect;
	short shape;
	short orn[3];
	short orns[3];
	unsigned long vtxnum[3];
	long edgnum[3];
	long material_index;
	double center[3];
	double axis[3];
	double radius;
	double area;
	struct phnedg *edg[3];
	struct circle *cir;
	struct phntri *next;
};

struct critlink {
	struct critlink *next;
	struct phnedg *edg;
	int orn;
};


struct edger {
	struct edger *right;
	struct edger *middle;
	struct edger *left;
	long first;
	long second;
	long sum;
	long number;
	
};

struct phnctr {
	struct phnctr *next;
	struct phnedg *head_phnedg;
	struct phnedg *tail_phnedg;
	long material_index;
	long n_phnedg;
};

struct record {
	short type;
	short nfield;
	short ntoken;
	struct token *first_token;
	struct token *last_token;
	struct record *next;
	struct record *head[MAX_FIELD];
	struct record *tail[MAX_FIELD];
	short field_type[MAX_FIELD];
	short field_count[MAX_FIELD];
};

struct token {
	char str[MAX_NAME];
	short type;
	short ivalue;
	float fvalue;
	struct token *next;
};



typedef long integer;
typedef float real;

enum PlexType {
	PolyVertex=0,
	PolyEdge=1,
	PolyTriangle=2,
	DenVertex=3,
	DenEdge=4,
	DenSquare=5,
	ProvVertex=6,
	ProvEdge=7,
	ProvTriangle=8,
	ProvTetrahedron=9,
	BorderVertex=10,
	BorderEdge=11,
	BorderTriangle=12
};

enum CubeType {
	EmptyCube=0,
	PartialCube=1,
	FullCube=2,
	BelowCube=3,
	AboveCube=4,
	JustBelow=5,
	JustAbove=6,
	AtCube=7
};

struct PlexVertex {
	enum PlexType type;
	real center[3];
	float vector[3];
	real value;
	struct PlexLink *head;
	struct Glass *glass;
};

struct PlexEdge {
	enum PlexType type;
	integer vns[2];
};

struct PlexTriangle {
	enum PlexType type;
	integer vns[3];
	integer ens[3];
};

struct PlexLink {
	struct PlexLink *next;
	short pN;
	short used;
	integer vns[2];
};

struct PlexCube {
	integer comp;
	float occupancy;
	real vector[3];
	struct PlexCube *next;
	enum CubeType type;
	short provinceNumber;
};

struct Province {
	integer vertexNumber;
	long cubeIndices[3];
	long n_vertex;
	real measure;
	real center[3];
};

struct PlexJoin {
	enum PlexType type;
	enum PlexType subtype;
	integer vertices[4];
	real measure;
	real apex[3];
	int orn;
};

struct Glass {
	enum PlexType type;
	enum PlexType subtype;
	struct Glass *next;
	struct Glass *contour;
	integer number;
	integer idnumbers[4];
	integer vertexnumbers[4];
	unsigned short pN[8];
	integer apex[2];
	/* int orn; */
	short npn;				/* number of distinct pN */
	short coefficient;
	real measure;
	real center[3];
};

struct GlassBlock {
	struct GlassBlock *next;
	struct Glass *glasses;
};

struct Plexi {
	struct Glass *head;
};

struct locube {
	unsigned long pN;
	integer dVN;
	integer pVN;
	struct PlexCube *cube;
	struct Glass *glass;
};

struct subedge {
	integer number;
	integer vns[2];
	struct subedge *parent;
	struct subedge *child[2];
};

struct subpolygon {
	struct subpolygon *parent;
	struct subpolygon *child[2];
	struct subedge *edges[MAXPED];
	struct subedge *se;
	short orn[MAXPED];
	short n_edge;
};

struct Plex {
	/* input */
	long dimension;
	long maxhash;
	long mvertex;
	long medge;
	long mtriangle;
	long n_vertex;
	long n_edge;
	long n_triangle;
	long *edges;
	long *triangles;
	double *vertices;
	double polyvolume;
	long dimensions[3];
	double bounds[2][3];
	double ctrlev;
	/* local */
	integer nglass[MaxPlexType];
	double cubeVolume;
	long glassblock;
	long maxvertex;
	long maxedge;
	long maxtriangle;
	long n_cube;
	long n_component;
	long n_inner;
	long n_surface;
	long n_province;
	long ijoin;
	long n_join;
	long n_cut;
	double *xgrid;
	double *ygrid;
	double *zgrid;
	struct subedge **subedges;
	struct subpolygon **subpolygons;
	struct PlexVertex *plexvertices;
	struct PlexEdge *plexedges;
	struct PlexTriangle *plextriangles;
	struct PlexCube *plexcubes;
	struct Plexi *plexis;
	struct Glass *glasses;
	struct Glass *freeglass;
	struct GlassBlock *head_glassblock;
	struct GlassBlock *tail_glassblock;
	struct Province *provinces;
	struct PlexJoin *joins;
};

/* structure definitions */

/* header structure:  */

struct header {
	char filetype[8];			/*   8 *  1 =   8       */
	long version;				/*   1 *  4 =   4       */
	long subversion;			/*   1 *  4 =   4       */
	double probe_radius;		/*   1 *  8 =   8 bytes */
	long counters[10];			/*  10 *  4 =  40 bytes */
	char name[64];				/*  64 *  1 =  64 bytes */
	char filler[384];			/* 384 *  1 = 384 bytes */
};								/*            512 bytes */


/* structures for binary surface file */
	
struct vtybin {
	long type;
	long atmnum[3];
	double center[3];
	double radius;
	double axis[3];
};

struct vtyshort {
	short type;
	short atmnum[3];
	float center[3];
	float radius;
	float axis[3];
};

struct vtxbin {
	long type;
	double center[3];
};

struct vtxshort {
	short type;
	float center[3];
};

struct cirbin {
	long type;
	long subtype;
	double center[3];
	double radius;
	double axis[3];
};

struct cirshort {
	short type;
	short subtype;
	float center[3];
	float radius;
	float axis[3];
};

struct arcbin {
	long type;
	long subtype;
	long cirnum;
	long vtxnum[2];
	long error;
};

struct arcshort {
	short type;
	short subtype;
	unsigned short cirnum;
	short vtxnum[2];
	short error;
};


struct facbin {
	long type;
	long vtynum;
	long fcynum;
	long comp;
	long error;
};

struct facshort {
	short type;
	short vtynum;
	short fcynum;
	short comp;
	short error;
};

struct cycbin {
	long type;
	long next;
	long fednum;
	long ned;
};

struct cycshort {
	short type;
	short next;
	long fednum;
	short ned;
};

struct edgbin {
	long type;
	long error;
	long arcnum;
};

struct edgshort {
	short type;
	short error;
	long arcnum;
};

struct cmpbin {
	long type;
	long subtype;
	double volume;
	double area;
	double accessible;
	double center[3];
};

struct cmpshort {
	short type;
	short subtype;
	float volume;
	float area;
	float accessible;
	float center[3];
};

struct neighbor {
	struct sphere *sphptr;
	struct pair *torptr;
};

struct mutual {
	struct sphere *sphptr;
	struct pair *tor13;
	struct pair *tor23;
};


/* temporary torus */
struct pair {
	struct arc *first_arc;			/* first concave arc */
	struct sphere *sph[2];		/* spheres defining torus */
	short buried;				/* completely inaccessible */
	short free;					/* completely accessible */
};

/* central */
struct central {
	short buried;			/* competely inaccessible */
	short free;				/* completely accessible */
	double center[3];
	double radius;
	double axis[3];			/* axis */
	struct central *next;	/* next central in linked list */
	struct arc *first_arc;		/* first concave arc */
	struct sphere *sph[2];	/* spheres defining torus */
};

/* torus */
/* torus created by rolling probe around pair of neighboring atoms */

struct torus {
	long number;
	short buried;			/* competely inaccessible */
	short free;				/* completely accessible */
    short   cusp;               /* flag for self-intersecting surface */
    short   reverse;            /* atom[1] < atom[0] */
    atomnum atom[2];            /* the two atoms defining the torus */
	double radius;
	double center[3];
	double axis[3];			/* interatomic axis */
	struct torus *next;		/* next torus in linked list */
	struct arc *first_arc;		/* first concave arc */
	struct sphere *atm[2];	/* atoms defining torus */
    struct circle  *central;        /* the central circle of the torus */
	struct circle *cir[2];	/* contact circles */
	struct cusp_extension *ce;
};

/* probe position tangent to three atoms */
struct probe {
	long number;
    atomnum atom[MAXPA];
	short low;
	short delinked;
	short natom;
	double area;
	double numerical_area;
	double unit_altitude_z;
	double height;
	double center[3];
	struct probe *next;		/* next probe in linked list */
	struct surface *srf;
    struct probe   *link[MAXPA];    /* pointers used for connected rolling */
	struct sphere *atm[MAXPA];	/* three atoms probe tangent to */
	struct pair *pairs[MAXPA];	/* three pairs (tori) probe tangent to */
	struct face *fac;		/* concave face of probe */
	struct cusp_link *first_cusp;	/* first cusp circle */
};


/* there are several kinds of circles: */
/* atom and probe latitudes and longitudes */
/* trajectories of probe rolling along a pair of atoms */

struct circle {
	long number;
	long lfn;
    atomnum atom;               /* atom circle lies on or near (or EOL) */
    short   gone;               /* flag for circle inaccessibility */
	short unused;
	short subtype;
	short filler;
    double width;                /* width of band for area calculations */
	double theta;				/* see 1983 JAC article */
	double radius;
	double center[3];
	double axis[3];
    double base[3];            /* zero direction for measuring angles */
	struct circle *next;
	struct sphere *atm;		/* pointer to atom (or null) */
    struct endpnt  *head;       /* pointer to head of list of arc endpoints */
};

/* vertex */
struct vertex {
	long number;
	long ofn;
	long lfn;
	unsigned char converted;
	unsigned char cusp;
	double center[3];		/* coordinates */
	struct probe *prb;		/* pointer to probe sphere */
	struct sphere *atm;		/* pointer to atom */
	struct vertex *next;
	struct cusp_extension *ce;
	struct phnvtx *pvt;
};

/* circular arc */
struct arc {
	long number;
	long ffn;
	long ofn;
	long lfn;
	short error;			/* error number */
	short subtype;
	unsigned int shape: 3;
	unsigned int small: 1;
	unsigned int converted: 1;
	unsigned int original: 1;
	unsigned int subdivided: 1;
	unsigned int eaten: 2;
	unsigned int eater: 2;
	unsigned int perm: 1;	/* cusp arc intersections permissable */
	unsigned int lune: 1;
	unsigned int filler: 3;
	double alpha;
	double phi;
	struct circle *cir;		/* circle that arc lies on */
	struct vertex *vtx[2];	/* vertices of arc (or null) */
	struct edge *edg[2];	/* edges pointing to arc */
	struct arc *next;		/* next arc in list for some object */
	struct cusp *csp;		/* pointer to cusp circle */
	struct face *fac;		/* concave face that arc lies on */
	struct torus *tor;		/* torus that arc lies on (or null) */
	struct phnedg *ped;		/* pointer to polyhedron edge */
};

/* edge of circular arc */
struct edge {
	short orn;				/* orientation (0:positive; 1: reverse) */
	short error;			/* error number */
	long lfn;
	long number;
	struct edge *next;		/* next edge in list for some object */
	struct arc *arcptr;		/* pointer to arc */
	struct face *fac;		/* pointer to face edge belongs to */
};

struct variety {
	long number;
	long lfn;
	long tube;
	short type;
	short natom;
	short atmnum[MAXPA];
	double length;
	double ccens[2][3];
	double center[3];
    double radii[2];
	double axis[3];
};

struct region {
	short type;
	double center[3];
    double radii[2];
	double axis[3];
};


/* face */
struct face {
	unsigned int shape: 3;
	unsigned int problem: 1;
	unsigned int original: 1;
	unsigned int largeok: 1;
	unsigned int converted: 1;
	unsigned int changed: 1;
	unsigned char notrim;
	short n_arc;			/* number of arcs forming boundary */
	short n_cycle;			/* number of cycles forming boundary */
	short chi;			/* Euler characteristic */
	short simplified;	/* made simply connected */
	short filler;
	long nedg;
	long ofn;
	long lfn;
	long comp;			/* component */
	long input_hue;
	double area;
	double alpha;
	double vol1;
	double vol2;
	union {				/* algebraic variety face lies on */
		struct sphere *atm;
		struct torus *tor;
		struct probe *prb;
	} ptr;
	struct surface *srf;
	struct cusp_extension *ce;
	struct face *next;		/* next face in linked list */
	struct arc *arcsp[MAX_ARC_ARRAY]; /* arcs of saddle, concave faces */
	struct cycle *first_cycle;	/* first cycle of arcs */
	struct edge *first_edge;	/* first edge (before cycle grouping) */
	struct variety *vty;
};

/* cycle of arc on convex face */
struct cycle {
	struct cycle *next;		/* next cycle for face */
	struct edge *first_edge;	/* pointer to first edge */
	long n_edge;
};


/* cusp circle */
struct cusp {
	struct circle *cir;		/* pointer to circle */
	struct probe *prb[2];	/* probes intersection for form cusp */
	struct cusp *next;	/* next cusp circle for molecule */
	struct edge *first_edge;	/* first arc list entry */
	short n_vertex;
	short full;
	short order;
	short reverse;
	short reversed;
	struct vertex *vtx[MAX_VERTEX_PER_CIRCLE];	/* for cusps */
	short orn[MAX_VERTEX_PER_CIRCLE];			/* 0: start of arc */
};

/* cusp circle list entry */
struct cusp_link {
	struct cusp *csp;		/* pointer to cusp circle */
	struct cusp_link *next;	/* next entry */
	short orn;				/* 1: torus axis agrees with concave tri */
};

struct cusp_extension {
	struct circle *nacirs[MAX_PASS][2];	/* non-axial circles for pass */
	struct probe *prbs[3];		/* pointer to three probe spheres */
	struct cusp_link *first_cusp;	/* first cusp circle */
	struct vertex *vtx[2];	/* cusp vertices */
	struct torus *tor;		/* pointer to torus */
	struct vertex *partner;		/* non-axial & volume code */
};


struct solid_angle {
	struct variety *vty;
	struct face *head_face;
	struct face *tail_face;
	struct circle *head_circle;
	struct circle *tail_circle;
	struct edge *head_edge;
	struct edge *tail_edge;
	struct cycle *head_cycle;
	struct cycle *tail_cycle;
	long n_cycle;
	long n_edge;
	long n_circle;
	long n_face;
	long numerical;
	long pis;
	long pip;
	double subtends;
	double radius;
	double omega;
};

/* structure for pair of vertices that may be okay */
struct vertex_pair {
	int i;
	int j;
	int shape;
	double value;
	double dij;
	struct vertex_pair *next;
	struct vertex *vtx1;
	struct vertex *vtx2;
	struct arc a;
	struct circle c;
};


struct element {
        int object_number;
        int coefficient;
        struct element *next;
};

struct generic {
        int copy_number;
        int symbol_number;
};

struct boolean {
        int copy_number;
        int symbol_number;
        int value;
};

struct integer {
        int copy_number;
        int symbol_number;
        int value;
        int option;
        int field;
        int adjective;
};

struct real {
        int copy_number;
        int symbol_number;
        double value;
};


struct string {
        int copy_number;
        int symbol_number;
        int length;
        char string_value[MAX_STRING];
};

struct set {
        int copy_number;
        int symbol_number;
        int n_element;
        int type;
        int special;
		int current_object;
        struct range *first;
        struct range *current_range;
};

struct range {
        int start;
        int end;
        struct range *next;
};

/* non-object structures: */

struct element_block {
        int n_element;
        struct element *first_element;
        struct element *last_element;
        struct element_block *next_block;
};

struct object_block {
        int number;
        int first_object;
        int last_object;
        char *start;
        struct object_block *next_block;
};

struct object_header {
        int number;
        char name[MAX_NAME];
        int named;
        int deletable;
        int copyable;
        int definable;
        int deep_copy_only;
        int shallow_copy_only;
        unsigned size;
        unsigned initial;
        unsigned increment;
        int top;
        int n_allocated;        /* top should be <= n_allocated */
        int n_free;
        int free_set;
        int all_set;
        int delete_set;
        int temporary_set;
        struct object_block *first_block;
        struct object_block *last_block;
};

struct symbol {
        int copy_number;
        char name[MAX_NAME];
        int token_type;
        int sub_type;
        int sub_sub_type;
        int number;
        int global;
        int procedure_number;
        int owner;
        int initialized;
        int dummy_arg;
        int next_symbol;
};


/* cluster of surface points for atom */

struct cluster {
    int     nmem;               /* number of surface points */
    double   area;               /* contact or reentrant area for atom */
    double  *points;             /* pointer to array of point coordinates */
    double  *normals;            /* pointer to array of unit vectors */
    struct cluster *next;       /* pointer to next cluster in linked list */
};

/* surface descriptor */

struct dsdesc {
 /* input fields */
    int     maxatom;            /* maximum number of atoms (for future use) */
    int     natom;              /* (current) number of atoms for molecule */
    double  *atmco;              /* pointer to atomic coordinate array */
    double  *atmrad;             /* pointer to atomic radii array */
    double  *atmden;             /* pointer to array of surface densities */
    short  *atmatt;             /* pointer to array of attention numbers */
    double   pradius;            /* probe radius */
    int     connected;          /* 1 = connected surface; 0 = complete surf */
 /* output fields */
    int     errflg;             /* non-zero value means fatal error */
    char   *errstr;             /* brief description of error */
    char    stage[80];          /* current stage in calculation */
    int     atom;               /* probable atom involved in error */
    int     verbose;            /* true: print out stages as they start */
    int     nsrfpnt;            /* total number of surface points for mol */
    double   cvxarea;            /* total convex area for molecule */
    double   renarea;            /* total reentrant area for molecule */
    struct cluster **cvxsp;     /* array of pointers to contact clusters */
    struct cluster **rensp;     /* array of pointers to reentrant clusters */
    struct hidden  *intern;     /* pointer to record hidden from user */
};

struct sector {
	double angle;
	double width;
	double value;
};

struct term {
	long order;
	double cos_coefficient;
	double sin_coefficient;
	double modulus;
	double phase;
};

struct bar {
	unsigned long shape;
	double fraction;
	struct bar *next;
};

struct oringe {
	struct oringe *next;
	long nsector;
	long nfrequency;
	long nbin1;
	long nsample;
	long nbar;
	double absnorm;
	double radius;
	double concavity;
	double center[3];
	double frame[3][3];
	double sector_density[MAX_SECTOR];
	struct sector *sectors[MAX_SECTOR];
	struct term *terms[DEFAULT_FREQUENCY];
	struct bar *head_bar;
	struct bar *tail_bar;
};

struct hedvtx {
	long number;
	long newber;			/* new number for contiguous */
	long central;			/* vertex number of nearby vertex to coallesce into */
	short used;
	short iscentral;
	short changed;
	short degree;
	short inc;
	short unreferenced;
	long idx;
	struct phnvtx *vtx;
};

struct hededg {
	long number;
	long newber;			/* new number for contiguous */
	long representative;		/* edge number of other with same vertices */
	short changed;
	short degenerate;
	short duplicate;
	short isrepresentative;
	short unreferenced;
	short filler;
	struct phnedg *edg;
	long vtxnum[2];
	struct hededg *next;
};

struct hedtri {
	long number;
	long newber;			/* new number for contiguous */
	long representative;		/* triangle number of other with same vertices */
	short changed;
	short degenerate;
	short duplicate;
	short isrepresentative;
	struct phntri *tri;
	long vtxnum[3];
	long edgnum[3];
	struct hedtri *next;
};

struct hedron {
	long maxvtx;
	long maxedg;
	long maxtri;
	long nvtx;
	long nedg;
	long ntri;
	long nchgvtx;
	long nchgedg;
	long nchgtri;
	long nvtxgrp;
	long van_der_Waals;
	struct hedvtx **vertices;
	struct hedvtx **changedVertices;
	struct hededg **edges;
	struct hededg **changedEdges;
	struct hedtri **triangles;
	struct hedtri **changedTriangles;
	struct vtxgrp *head_vtxgrp;
	struct vtxgrp *tail_vtxgrp;
	long *neighbors;
	struct hededg **edghash;
	struct hedtri **trihash;
};

struct vtxgrp {
	long central;
	long nvtx;
	double radius;
	struct vtxgrp *next;
	struct hedvtx *vtxptrs[MAXVTXGRP];
};

struct cept {
    int type;
    int subtype;
    int severity;
    int object_type;
    int variable_type;
    int n_function;
    int n_message;
    int n_double;
    int n_long;
    int n_atom;
    struct cept *next;
    char object[MAX_NAME];
    char variable[MAX_NAME];
    char remedy[MAX_NAME];
    char array_name[MAX_NAME];
    char input_name[MAX_NAME];
    char source[MAX_NAME];
    char functions[MAX_EXCEPTION][MAX_NAME];
    char messages[MAX_EXCEPTION][MAX_NAME];
    char double_names[MAX_EXCEPTION][MAX_NAME];
    double doubles[MAX_EXCEPTION];
    char long_names[MAX_EXCEPTION][MAX_NAME];
    long longs[MAX_EXCEPTION];
    struct sphere *atoms[MAX_EXCEPTION];
};

struct pholder {
	struct pholder *next;
	struct probe *prb;
};

struct quartet {
	int delete;
	atomnum atm[MAXPA];
	struct probe *prb[MAXPQ];
	long probe_number[MAXPQ];
	struct sphere *atmptr[MAXPA];
	struct pair *torptr[MAXPA];
	double center[3];
	struct quartet *next;
	struct quartet *previous;
};

struct packer {
	union {
		unsigned short twoshort[2];
		unsigned long onelong;
	} shortlong;
};

struct depth_buffer {
	short *heights;
	unsigned char *shades;
	unsigned char *hues;
	unsigned char *alphas;
	unsigned char *inners;
};



