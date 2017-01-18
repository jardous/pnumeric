#ifndef __PNUMERIC_H__
#define __PNUMERIC_H__

#define MAX(a,b) a>b? a:b
#define MIN(a,b) a>b? b:a

//#define DEBUG(...) printf(__VA_ARGS__)
#define DEBUG(...)  // __VA_ARGS__

#include "matrix.h"
#include "vector.h"

#include "cgensupport.h"  // for PyArg_GetDoubleArray

#endif
