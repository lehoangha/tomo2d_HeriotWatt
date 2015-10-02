/*
 * error.cc
 *
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#include <iostream.h>
#include <cstdlib>
#include "error.h"

void error(string s)
{
    cerr << s << '\n';
    exit(1);
}

