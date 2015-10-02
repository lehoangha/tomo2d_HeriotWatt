/*
 * error.cc
 *
 * Jun Korenaga, MIT/WHOI
 * January 1999
 */

#include <iostream>
#include <cstdlib>
#include "error.h"

using namespace std;

void error(string s)
{
    cerr << s << '\n';
    exit(1);
}

