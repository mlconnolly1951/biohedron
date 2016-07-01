/* MSForm Copyright 1998 by Michael L. Connolly */
/* Last revised: September 1, 1999 */

#include "ms.h"
#include "msenum.h"
#include "mslimit.h"
#include "msdefault.h"
#include "msstruct.h"
#include "msglobal.h"
#include "msfunction.h"


/* record */

struct token *get_token (FILE *in)
{
	int c, ch, state;
	int gotalpha, gotint, gotfloat, gotpunct;
	char word[MAX_STRING];
	struct token *tok;

	if (in == NULL) return (NULL);
	tok = allocate_token ();
	if (tok == NULL) {
		set_error1 ("get_token: no more memory");
		return (NULL);
	}
	gotalpha = 0;
	gotint = 0;
	gotfloat = 0;
	gotpunct = 0;
	state = BEGIN_STATE;
	c = 0;
	while (!feof (in)) {
		ch = getc (in);
		if (ferror (in)) {
			set_error1 ("get_token: disk input error");
			return (NULL);
		}
		if (feof (in)) break;
		if (c >= MAX_STRING-1) {
			set_error1 ("get_token: word too long");
			return (NULL);
		}
		/* check for whitespace (includes comma) */
		if (ch == ' ' || ch == '\t' || ch == '\n' || ch == ',') {
			if (state == BEGIN_STATE) {
				state = WHITE_STATE;
			}
			else if (state == WHITE_STATE) continue;
			else if (state == ALPHA_STATE) {
				word[c] = (char) 0;
				state = END_STATE;
				gotalpha = 1;
				break;
			}
			else if (state == INT_STATE) {
				word[c] = (char) 0;
				state = END_STATE;
				gotint = 1;
				break;
			}
			else if (state == FLOAT_STATE) {
				word[c] = (char) 0;
				state = END_STATE;
				gotfloat = 1;
				break;
			}
			else {
				state = INVALID_STATE;
				break;
			}
		}
		/* alphabetic */
		else if (isalpha (ch)) {
			if (state == BEGIN_STATE) {
				word[c++] = (char) ch;
				state = ALPHA_STATE;
			}
			else if (state == WHITE_STATE) {
				word[c++] = (char) ch;
				state = ALPHA_STATE;
			}
			else if (state == ALPHA_STATE) {
				word[c++] = (char) ch;
			}
			else if (state == INT_STATE) {
				state = INVALID_STATE;
				break;
			}
			else if (state == FLOAT_STATE) {
				state = INVALID_STATE;
				break;
			}
			else {
				state = INVALID_STATE;
				break;
			}
		}
		else if (isdigit (ch)) {
			if (state == BEGIN_STATE) {
				word[c++] = (char) ch;
				state = INT_STATE;
			}
			else if (state == WHITE_STATE) {
				word[c++] = (char) ch;
				state = INT_STATE;
			}
			else if (state == INT_STATE) {
				word[c] = (char) ch;
				c++;
			}
			else if (state == ALPHA_STATE) {
				state = INVALID_STATE;
				break;
			}
			else if (state == FLOAT_STATE) {
				word[c] = (char) ch;
				c++;
			}
			else {
				state = INVALID_STATE;
				break;
			}
		}
		else if (ispunct (ch)) {
			if (state == BEGIN_STATE) {
				if (ch == '-') {
					word[c++] = (char) ch;
					state = INT_STATE;
				}
				else {
					gotpunct = 1;
					state = END_STATE;
					break;
				}
			}
			else if (state == WHITE_STATE) {
				if (ch == '-') {
					word[c++] = (char) ch;
					state = INT_STATE;
				}
				else {
					gotpunct = 1;
					state = END_STATE;
					break;
				}
			}
			else if (state == INT_STATE) {
				if (ch == '.') {
					word[c++] = (char) ch;
					state = FLOAT_STATE;
				}
				else {
					ungetc (ch, in);
					word[c] = (char) 0;
					gotint = 1;
					state = END_STATE;
					break;
				}
			}
			else if (state == FLOAT_STATE) {
				ungetc (ch, in);
				word[c] = (char) 0;
				gotfloat = 1;
				state = END_STATE;
				break;
			}
			else if (state == ALPHA_STATE) {
				ungetc (ch, in);
				word[c] = (char) 0;
				gotalpha = 1;
				state = END_STATE;
				break;
			}
			else {
				state = INVALID_STATE;
				break;
			}
		}
	}
	if (feof (in)) {
		tok -> type = EOF_TOKEN;
		return (tok);
	}
	if (state == INVALID_STATE) {
		tok -> type = INVALID_TOKEN;
		return (tok);
	}
	if (state != END_STATE) {
		tok -> type = INVALID_TOKEN;
		return (tok);
	}
	if (gotalpha) {
		strcpy (tok -> str, word);
		tok -> type = ALPHA_TOKEN;
	}
	else if (gotint) {
		tok -> ivalue = (short) atoi (word);
		tok -> type = INT_FIELD;
		/* for floats written without decimal point */
		tok -> fvalue = tok -> ivalue;
	}
	else if (gotfloat) {
		tok -> fvalue = (float) atof (word);
		tok -> type = FLOAT_FIELD;
	}
	else if (gotpunct) {
		if (ch == ':')
			tok -> type = COLON_TOKEN;
		else if (ch == '{')
			tok -> type = OPEN_BRACE;
		else if (ch == '}')
			tok -> type = CLOSE_BRACE;
		else tok -> type = INVALID_TOKEN;
	}
	else {
		tok -> type = INVALID_TOKEN;
	}
	lookup_alpha (tok);
	if (error ()) return (NULL);
	return (tok);
}

int lookup_alpha (struct token *tok)
{
	char *str;
	int type;

	if (tok -> type != ALPHA_TOKEN) return(0);
	str = tok -> str; type = 0;
	if (strcmp (str, "int") == 0) type = INT_FIELD;
	else if (strcmp (str, "ints") == 0) type = INT_FIELD;
	else if (strcmp (str, "float") == 0) type = FLOAT_FIELD;
	else if (strcmp (str, "floats") == 0) type = FLOAT_FIELD;
	else if (strcmp (str, "oringe") == 0) type = ORINGE_RECORD;
	else if (strcmp (str, "polygon") == 0) type = POLYGON_RECORD;
	else if (strcmp (str, "center") == 0) type = CENTER_RECORD;
	else if (strcmp (str, "rotation") == 0) type = ROTATION_RECORD;
	else if (strcmp (str, "concavity") == 0) type = CONCAVITY_RECORD;
	else if (strcmp (str, "radius") == 0) type = RADIUS_RECORD;
	else if (strcmp (str, "radii") == 0) type = RADIUS_RECORD;
	else if (strcmp (str, "nsector") == 0) type = NSECTOR_RECORD;
	else if (strcmp (str, "vertex") == 0) type = VERTEX_RECORD;
	else if (strcmp (str, "vertices") == 0) type = VERTEX_RECORD;
	else if (strcmp (str, "sector") == 0) type = SECTOR_RECORD;
	else if (strcmp (str, "sectors") == 0) type = SECTOR_RECORD;
	else if (strcmp (str, "frequency") == 0) type = FREQUENCY_RECORD;
	else if (strcmp (str, "frequencies") == 0) type = FREQUENCY_RECORD;
	else if (strcmp (str, "bar") == 0) type = BAR_RECORD;
	else if (strcmp (str, "bars") == 0) type = BAR_RECORD;
	return (type);
}

int fundamental (int type)
{
	return (type == INT_FIELD || type == FLOAT_FIELD);
}



struct record *get_record (short expected_type, FILE *in)
{
	int t, i, count;
	int nfundamental, nrecursive;
	int rtype, ftype, type, nfield;
	char message[MAXLINE];
	struct token *recordToken, *colonToken;
	struct token *fieldCountToken, *fieldNameToken, *closingBraceToken;
	struct token *tok;
	struct record *rec, *inner;

	/* read record type */
	recordToken = get_token (in);
	if (error ()) return (NULL);
	if (recordToken -> type == INVALID_TOKEN) {
		sprintf (message, "(get_record) first token invalid: %s", recordToken -> str);
		set_error1 (message);
		return (NULL);
	}
	if (recordToken -> type == EOF_TOKEN) {
		return (NULL);
	}
	if (recordToken -> type != ALPHA_TOKEN) {
		sprintf (message, "(get_record) first token not alpha: %s", recordToken -> str);
		set_error1 (message);
		return (NULL);
	}
	rtype = lookup_alpha (recordToken);
	if (error ()) return (NULL);
	if (rtype == 0) {
		sprintf (message, "(get_record) first token invalid record type: %s", recordToken -> str);
		set_error1 (message);
		return (NULL);
	}
	if (rtype != expected_type) {
		sprintf (message, "(get_record) first token wrong record type: %s", recordToken -> str);
		set_error1 (message);
		return (NULL);
	}
	if (recordToken != NULL) free_token (recordToken);
	colonToken = get_token (in);
	if (error ()) return (NULL);
	if (colonToken -> type != COLON_TOKEN) {
		sprintf (message, "(get_record) second token not colon: %s", colonToken -> str);
		set_error1 (message);
		return (NULL);
	}
	if (colonToken != NULL) free_token (colonToken);
	rec = allocate_record ();
	rec -> type = rtype;
	nfield = 0; nfundamental = 0; nrecursive = 0;
	for (t = 0; t <= MAX_FIELD; t++) {
		fieldCountToken = get_token (in);
		if (error ()) return (NULL);
		if (fieldCountToken -> type == OPEN_BRACE) break;
		if (t == MAX_FIELD) break;
		if (fieldCountToken -> type != INT_FIELD) {
			sprintf (message, "(get_record) token not int type: %s", fieldCountToken -> str);
			set_error1 (message);
			return (NULL);
		}
		rec -> field_count[t] = fieldCountToken -> ivalue;
		fieldNameToken = get_token (in);
		if (error ()) return (NULL);
		if (fieldNameToken -> type != ALPHA_TOKEN) {
			sprintf (message, "(get_record) token not record type: %s", fieldNameToken -> str);
			set_error1 (message);
			return (NULL);
		}
		ftype = lookup_alpha (fieldNameToken);
		if (error ()) return (NULL);
		if (ftype == 0) {
			sprintf (message, "(get_record) field token invalid record type: %s", fieldNameToken -> str);
			set_error1 (message);
			return (NULL);
		}
		rec -> field_type[t] = ftype;
		nfield++;
		if (fundamental (ftype)) nfundamental++;
		else nrecursive++;
		if (fieldNameToken != NULL) free_token (fieldNameToken);
	}
	if (fieldCountToken == NULL || !(fieldCountToken -> type == OPEN_BRACE)) {
		set_error1 ("(get_record) token not open brace");
		return (NULL);
	}
	if (fieldCountToken != NULL) free_token (fieldCountToken);
	if (nfield <= 0) {
		set_error1 ("(get_record) no fields for record");
		return (NULL);
	}
	if (nrecursive > 0 && nfundamental > 0) {
		set_error1 ("(get_record) mixed fundamental & non-fundamental fields for record");
		return (NULL);
	}
	rec -> nfield = nfield;
	rec -> ntoken = 0;
	rec -> first_token = NULL;
	rec -> last_token = NULL;
	for (t = 0; t < nfield; t++) {
		count = rec -> field_count[t];
		type = rec -> field_type[t];
		rec -> head[t] = NULL;
		rec -> tail[t] = NULL;
		if (fundamental (type)) {
			for (i = 0; i < count; i++) {
				tok = get_token (in);
				if (error ()) return (NULL);
				if (tok -> type == INVALID_TOKEN) {
					set_error1("(get_record): invalid token");
					return(NULL);
				}
				if (rec -> first_token == NULL)
					rec -> first_token = tok;
				else rec -> last_token -> next = tok;
				rec -> last_token = tok;
				rec -> ntoken++;
			}
		}
		else {
			/* recursive call */
			for (i = 0; i < count; i++) {
				inner = get_record (type, in);
				if (error ()) return (NULL);
				/* link into tree */
				if (rec -> head[t] == NULL)
					rec -> head[t] = inner;
				else rec -> tail[t] -> next = inner;
				rec -> tail[t] = inner;
			}
		}
	}
	closingBraceToken = get_token (in);
	if (error ()) return (NULL);
	if (closingBraceToken == NULL || !(closingBraceToken -> type == CLOSE_BRACE)) {
		set_error1 ("(get_record) token not close brace");
		return (NULL);
	}
	if (closingBraceToken != NULL) free_token (closingBraceToken);
	return (rec);
}

struct record *allocate_record ()
{
	struct record *rec;

	rec = (struct record *) allocate_object (RECORD);
	if (rec == NULL) {
		set_error1 ("(allocate_record): mem alloc fails");
		return(NULL);
	}
	return (rec);
}

struct token *allocate_token ()
{
	struct token *tok;

	tok = (struct token *) allocate_object (TOKEN);
	if (tok == NULL) {
		set_error1 ("(allocate_token): mem alloc fails");
		return(NULL);
	}
	return (tok);
}

void free_token (struct token *tok)
{
	if (tok != NULL) free_object (TOKEN, (short *) tok);
}

void free_record (struct record *rcd)
{
	if (rcd != NULL) free_object (RECORD, (short *) rcd);
}

void deep_free_record (struct record *rcd)
{
	int f;
	struct token *tok, *next_tok;
	struct record *r, *n;

	for (f = 0; f < rcd -> nfield; f++) {
		if (rcd -> field_count[f] > 0 && rcd -> head[f] != NULL) {
			for (r = rcd -> head[f]; r != NULL; r = n) {
				n = r -> next;
				deep_free_record (r);
			}
			rcd -> head[f] = NULL;
		}
	}
	for (tok = rcd -> first_token; tok != NULL; tok = next_tok) {
		next_tok = tok -> next;
		free_token (tok);
	}
	rcd -> first_token = NULL;
	free_record (rcd);
}

/* MSP Copyright 1998 by Michael L. Connolly */
