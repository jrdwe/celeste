#ifndef PARSING_H
#define PARSING_H

/* func: parses the optional parameters (b & f)
 *       data that has been passed in
 *        
 * return: void
 */
void parseParameters(char** argv, int argc, int* flag, int* bodies, double* upper, double* lower, int* testing, char** fstring); 

/* func: read in data points from file
 *
 * return: void
 */
void readDataFromFile(char* filename, double* bodies);

/* func: read number of lines in the file
 *
 * return: integer
 */
int countLines(char* filename);

#endif
