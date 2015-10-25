
#define bool int
#define false 0
#define true (!false)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <time.h>
#include "clRTOCT.h"


int main(int argc, char **argv);
int readTwoColumnCSVFile(char* sourceFile, unsigned int* len, unsigned int maxRows, float* col1Output, float* col2Output);
//int saveBitmap(char* path, float* sourceArray, unsigned short width, unsigned short sourceHeight, unsigned outputHeight, float minVal, float maxVal);
