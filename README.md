# tomo2d_HeriotWatt
This is a fork project from the tomo2d, created by Jun Korenaga at Yale university, aimed for helping seismic project(s) of the Institute of Petroleum, the Heriot Watt University 

## ACKNOWLEDGEMENT
The project was originally forked from tomo2d repository at the Yale university, website: [http://people.earth.yale.edu/software/jun-korenaga](http://people.earth.yale.edu/software/jun-korenaga). It is created for supporting the work of Ms. Kim Phung Nguyen, a research assistant at the Institute of Petroleum, the Heriot Watt University. All changes and contribution should be reported back to the mentioned 2 parties.

## TECHNICAL REVIEW

There are changes of source code in order to make the code be possibly compiled by the GNU Compiler Collection (gcc) version 4.9.2:
* the iostream.h declared in several "error.cc" were replace by iostream, which is a newer library
* Added the library cstdlib, and declared "using namespace std;" in some source files:
 * util.cc
 * corrlen.cc
 * lsqr.cc
* in "heap_deque.h", any of push\_back actions were replaced by this->push\_back
