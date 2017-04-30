/* Molecular Surface Package copyright 1996 by Michael L. Connolly */
/* last revised: January 27, 2006 */
/* should be enums */

#define TRUE 1
#define FALSE 0

/* end of atom list - used by ds */
#define EOL 65535

#define   UNKNOWN_ERROR 0
#define   GEOMETRY_ERROR 1
#define   MEMORY_ERROR 2
#define   INPUT_ERROR 3
#define   PARAMETER_ERROR 4
#define   SYNTAX_ERROR 5
#define   ARRAY_ERROR 6
#define   LOGIC_ERROR 7
#define   GRID_ERROR 8
#define   RETURN_ERROR 9
#define   POINTER_ERROR 10
#define   ENUM_ERROR 11

#define   UNKNOWN_SUBTYPE 0
#define   ALLOCATION 1
#define   BOUNDS 2
#define   FREEING 3
#define   INCONSISTENCY 4
#define   FILE_NAME 5
#define   NULL_VALUE 6
#define   INVALID_VALUE 7
#define   NEGATIVE_VALUE 8
#define   CONTENT 9
#define   PREMATURE_EOF 10
#define   NULL_POINTER 11
#define   MSOVERFLOW 12
#define   DEGENERACY 13
#define   MSUNDERFLOW 14
#define   NOT_FOUND 15
#define   INFINITE_LOOP 16
#define   UNDEFINED_VALUE 17

#define   FATAL_SEVERITY 1
#define   WARNING_SEVERITY 0


/* rendering point locations */
#define OUTSIDE             0
#define INSIDE              1
#define INACCESSIBLE        0
#define ACCESSIBLE          1
#define CLIPPED             0
#define UNCLIPPED           1

/* object scheme: kinds of coloring */
#define UNIFORM_COLORING    1
#define INPUT_COLORING      2
#define ATOM_COLORING       3
#define SHAPE_COLORING      4
#define COMPONENT_COLORING  5
#define VERTEX_COLORING     6

/* rendering: kinds of opacity */
#define UNIFORM_OPACITY    1
#define INPUT_OPACITY      2
#define ATOM_OPACITY       3
#define SHAPE_OPACITY      4
#define COMPONENT_OPACITY  5
#define VERTEX_OPACITY     6

/* plotting: types of triangles */
#define FRONTFACING 1
#define BACKFACING 2
#define SILHOUETTE 3
#define BOUNDARY 4

/* data type for msdata routines */
#define LONGDATA 1
#define FLOATDATA 2
#define BYTEDATA 3

/* atom attension (srn) numbers - used by ds */
#define IGNORE 0                /* ignore atom completely */
#define BLOCKER 1               /* blocks probe path, but needs no surface */
#define AREA 2                  /* compute only contact surface area for atom */
#define POINTS 3                /* compute C & R areas and points for atom */
#define NORMALS 4               /* compute areas, points, normals for atom */

/* kinds of surfaces */
#define PQMS_SURFACE   1
#define BAS_SURFACE    2
#define PHN_SURFACE    3
#define DEN_SURFACE    4
#define CTR_SURFACE    5
#define NML_SURFACE    6
#define TAN_SURFACE    7
#define BTN_SURFACE    8
#define ORG_SURFACE    9

/* surface component types */
#define OUTER_SUBTYPE 1
#define INNER_SUBTYPE 2

/* input and output formats */
#define XYZR_FORMAT    1
#define MS_FORMAT      2
#define PDB_FORMAT     3
#define DS_FORMAT      4
#define HPGL           5
#define POSTSCRIPT     6
#define SUNRASTER      7
#define SGIIMAGE       8
#define AVSIMAGE       9
#define AVSFIELD      10
#define AVSDENSITY    11
#define INVENTOR      12
#define WHATIF        13
#define O             14
#define VET           15
#define BMP           16

/* shapes: */
#define CONVEX 1
#define SADDLE 2
#define CONCAVE 3
#define FLAT 4
#define STRAIGHT 4
#define CYLINDRICAL 5

/* variety types */
#define NOWHERE 1
#define SPACE 31
#define LINEV 19
#define POINTV 23
#define S0 29

/* relational: */
#define EQUAL_TO 1
#define NOT_EQUAL 2
#define LESS_THAN 3
#define GREATER_THAN 4
#define LESS_EQUAL 5
#define GREATER_EQUAL 6
#define MEMBER_OF 7
#define SUBSET_OF 8

#define MAX_FIELD 16
#define BEGIN_STATE 1
#define END_STATE 2
#define WHITE_STATE 3
#define ALPHA_STATE 4
#define INT_STATE 5
#define FLOAT_STATE 6
#define PUNCT_STATE 7
#define INVALID_STATE 8

#define INVALID_TOKEN        1
#define ALPHA_TOKEN          2
#define COLON_TOKEN          3
#define OPEN_BRACE           4
#define CLOSE_BRACE          5
#define EOF_TOKEN            6
#define INT_FIELD            7
#define FLOAT_FIELD          8
#define ORINGE_RECORD       10
#define CENTER_RECORD       11
#define RADIUS_RECORD       12
#define ROTATION_RECORD     13
#define NSECTOR_RECORD      14
#define VERTEX_RECORD       15
#define SECTOR_RECORD       16
#define BAR_RECORD          17
#define POLYGON_RECORD      18
#define FREQUENCY_RECORD    19
#define CONCAVITY_RECORD    20

/* polyhedron_and_density actions regarding storing results */

#define GET_LOCAL            1
#define GET_NORMAL           2
#define GET_FOURIER0         3
#define GET_FOURIER1         4
#define GET_FRAME            5
#define GET_ORINGE           6


/* types and subtypes for binary pqms file */
#define ATOM_TYPE 1
#define TORUS_TYPE 2
#define PROBE_TYPE 3
#define VERTEX_TYPE 4
#define CIRCLE_TYPE 5
#define CONVEX_ARC_TYPE 6
#define CONCAVE_ARC_TYPE 7
#define CONVEX_FACE_TYPE 8
#define SADDLE_FACE_TYPE 9
#define CONCAVE_FACE_TYPE 10
#define CYCLE_TYPE 11
#define EDGE_TYPE 12
#define COMPONENT_TYPE 13

/* circle types: */
#define CONTACT_SUBTYPE 1
#define GREAT_SUBTYPE 2
#define CUSP_SUBTYPE 3

/* arc types: */
#define SHORTENED_SUBTYPE 1
#define AXIAL_SUBTYPE 2
#define NONAXIAL_SUBTYPE 3
#define RESORT_SUBTYPE 4
#define REVERSE_SUBTYPE 5

/* memory types */
#define ONEOBJ 1
#define OBJECTS 2
#define OBJPTRS 3
#define STRPTRS 4
#define CHARS 5
#define BYTES 6
#define SHORTS 7
#define LONGS 8
#define FLOATS 9
#define DOUBLES 10
#define MEMORY_TYPES 10

/* variable types */
#define ARC_DIRECTION 1
#define ATMCO 2
#define ATMRAD 3
#define ATOM_ALPHAS 4
#define ATOM_COLORS 5
#define ATOM_OPACITIES 6
#define ATOM_CENTERS 7
#define ATOM_RADII 8
#define ATMDEN 9
#define BAD_VERTEX 10
#define BAD_NEIGHBOR 11
#define BASES 12
#define BIT1_GRID 13
#define BIT2_GRID 14
#define BOUND_SAME 15
#define CENTERS 16
#define CIR 17
#define CIRCLE1 18
#define CIRCLE2 19
#define CIRCLE3 20
#define CIRCLE_PTR 21
#define CONE_CIRCLE 23
#define CONTAINS 24
#define CNUMS 25
#define CYCLE_USED 26
#define DISTANCE_BACKWARD 27
#define DISTANCE_FORWARD 28
#define EDGES 29
#define EDGVTX 30
#define ENUMBERS0 31
#define ENUMBERS1 32
#define EXTERIOR_ANGLE 33
#define F 34
#define F1 35
#define F2 36
#define F3 37
#define FIRST_CYCLE_NUMBER 38
#define FIRST_EDGE_NUMBER 39
#define FN1 40
#define FN2 41
#define FN3 42
#define FOCI 43
#define HEMI_CIRCLE 44
#define INTEGERS 45
#define LFCIR 46
#define MEDIAN 47
#define NORMALSV 48
#define NEIGHBORS 49
#define NFC 50
#define NUMBERS 51
#define PERIMETER 52
#define POINTSV 53
#define POLYGON_SIDE 54
#define POLYGON_TANGENT 55
#define PRBDOT 56
#define PROJECTED_VERTEX_LIST 57
#define PROJECTED_VERTICES 58
#define RADII 59
#define REALS 60
#define TORCIR 61
#define TRIANGLES 62
#define TRIEDGVTX 63
#define USED 64
#define VERTEX_CENTERS 65
#define VERTEX_NORMALS 66
#define VERTICES 67
#define VERTS 68
#define VNUMBERS0 69
#define VNUMBERS1 70
#define XGRID 71
#define YGRID 72
#define ZGRID 73
#define VARIABLE_TYPES 73

/* object types: */
#define SPHERE 1
#define TORUS 2
#define PLANE 3
#define CYLINDER 4
#define CONE 5
#define ATOM 6
#define PROBE 7
#define CENTRAL 8
#define VERTEX 9
#define CIRCLE 10
#define ARC 11
#define FACE 12
#define CYCLE 13
#define EDGE 14
#define CUSP 15
#define CUSP_LINK 16
#define COMPONENT 17
#define CHUNK 18
#define VARIETY 19
#define PVERTEX 20
#define PEDGE 21
#define PTRIANGLE 22
#define EVALPNT 23
#define PHNVTX 24
#define PHNEDG 25
#define PHNTRI 26
#define LEAF 27
#define OBJECT_SCHEME 28
#define SURFACE 29
#define SOLID_ANGLE 30
#define PHNCTR 31
#define MOLECULE 32
#define PLOT 33
#define LINSEG 34
#define POLYGON 35
#define POLYHEDRON 36
#define CONTOURS 37
#define HAIR 38
#define STICKS 39
#define BOND 40
#define COLOR_RAMP 41
#define REGION 42
#define INTEGER 43
#define ELEMENT 44
#define REAL 45
#define SET 46
#define STRING 47
#define SYMBOL 48
#define BOOLEAN 49
#define ORINGE 50
#define PLEXLINK 51
#define SUBEDGE 52
#define SUBPOLYGON 53
#define RECORD 54
#define TOKEN 55
#define SECTOR 56
#define TERM 57
#define BAR 58
#define HEDRON 59
#define HEDVTX 60
#define HEDEDG 61
#define HEDTRI 62
#define VTXGRP 63
#define EDGER 64
#define NEIGHBOR 65
#define MUTUAL 66
#define PAIR 67
#define VANITY 68
#define LAX 69
#define SQUARE 70
#define PLEX 71
#define PLEXI 72
#define PLEXVERTEX 73
#define PLEXEDGE 74
#define PLEXTRIANGLE 75
#define PLEXCUBE 76
#define PLEXJOIN 77
#define PROVINCE 78
#define GLASS 79
#define SCUBE 80
#define SNUMBER 81
#define SNUMBERBLOCK 82
#define RESIDUEINDEX 83
#define PATTERNRECORD 84
#define RADIUSRECORD 85
#define SPHEREGRID 86
#define DENVTX 87
#define MSDATA 88
#define MSSCENE 89
#define VERTEX_PAIR 90
#define CRITLINK 91
#define CUSP_EXTENSION 92
#define MATERIAL_TABLE 93
#define ARRAY 94
#define DSDESC 95
#define HIDDEN 96
#define BOX 97
#define ELEMENT_BLOCK 98
#define RUN 99
#define TLIST 100
#define PLIST 101
#define EDG 102
#define ENDPNT 103
#define MIDPLN 104
#define CLUSTER 105
#define OBJECT_HEADER 106
#define OBJECT_BLOCK 107
#define GLASS_BLOCK 108
#define CEPT 109
#define PHOLDER 110
#define QUARTET 111
#define DEPTH_BUFFER 112
#define N_OBJECTS 112
