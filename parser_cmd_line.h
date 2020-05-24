#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <algorithm>

//true if option exists, false otherwise
bool option_exists (char** start, char** stop, const std::string& option);

//returns option parameter (only 1 peremter is assumed)
char* option_parameter (char** start, char** stop, const std::string& option);

//parser for the whole command line returns true if everything is correct, false otherwise
bool parse_cmd_line(const int& argc, char ** argv);
