#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "basic_structures.h"

#define SKIP_COMMENT_LINE(file) { while(getc(file)!='\n'); }
#define SKIP_REMAINING_CHARS(file) SKIP_COMMENT_LINE(file)
#define SKIP_OPTION(file) { SKIP_COMMENT_LINE(file) SKIP_REMAINING_CHARS(file) }

void get_description(FILE* file, char* description);

void parse_option(FILE* file, const char* pattern, const char* option, void* ptr);

void parse_lambda_option(FILE* file, const char* option);

void parse_exclusion_option(FILE* file_const, const char* option, correlator* pC);

void parse_force_zero_option(FILE* file, const char* option);

bool parse_const_file(FILE* file_const, correlator* pC);

void set_default_values();
