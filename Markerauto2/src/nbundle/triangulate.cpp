#include "triangulate.h"
#include <cmath>
#include <ctime>
#include <cstring>
#include <iostream>
#include <cstdlib>
#include "util/qsort.h"

extern "C"{
#include <clapack.h>
#include <fblaswr.h>
#include <f2c.h>
#include <math.h>
#include <sys/stat.h>
}
