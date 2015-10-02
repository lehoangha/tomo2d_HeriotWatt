# Makefile

CCC = /usr/bin/g++
CC  = /usr/bin/gcc
CFLAGS = -ansi -O3 #-pg
HOME = /home/lehoangha/
PWD = /home/lehoangha/Distrib
INCLUDE = -I${PWD}

SRC_GEN_SMESH = gen_smesh.cc zeltform.cc \
                geom.cc error.cc util.cc
OBJ_GEN_SMESH = ${SRC_GEN_SMESH:.cc=.o} 

SRC_EDIT_SMESH = edit_smesh.cc \
                 geom.cc error.cc util.cc smesh.cc interface.cc corrlen.cc
SRC_EDIT_SMESH_C = d_jacobi.c d_choldc.c nrutil.c
OBJ_EDIT_SMESH = ${SRC_EDIT_SMESH:.cc=.o} ${SRC_EDIT_SMESH_C:.c=.o}

SRC_STAT_SMESH = stat_smesh.cc error.cc interface.cc util.cc
OBJ_STAT_SMESH = ${SRC_STAT_SMESH:.cc=.o}

SRC_TT_FORWARD = tt_forward.cc syngen.cc smesh.cc graph.cc \
                 traveltime.cc betaspline.cc \
                 bend.cc bend_mnbrak.cc bend_brent.cc \
                 interface.cc corrlen.cc \
                 geom.cc error.cc util.cc
SRC_TT_FORWARD_C = d_jacobi.c d_choldc.c nrutil.c
OBJ_TT_FORWARD = ${SRC_TT_FORWARD:.cc=.o} ${SRC_TT_FORWARD_C:.c=.o} 

SRC_TT_INVERSE = tt_inverse.cc inverse.cc sparse_rect.cc lsqr.cc \
                 smesh.cc graph.cc \
                 traveltime.cc betaspline.cc \
                 bend.cc bend_mnbrak.cc bend_brent.cc \
                 interface.cc corrlen.cc \
                 geom.cc error.cc util.cc jgrav.cc
SRC_TT_INVERSE_C = d_jacobi.c d_choldc.c nrutil.c d_realft.c d_four1.c
OBJ_TT_INVERSE = ${SRC_TT_INVERSE:.cc=.o} ${SRC_TT_INVERSE_C:.c=.o} 

.c.o: 
	$(CCC) $(CFLAGS) -c  $< $(INCLUDE)
.cc.o: 
	$(CCC) $(CFLAGS) -c  $< $(INCLUDE)

all: gen_smesh edit_smesh stat_smesh tt_forward tt_inverse

gen_smesh: $(OBJ_GEN_SMESH)
	$(CCC) $(CFLAGS) -o $@ $(OBJ_GEN_SMESH) -lm

edit_smesh: $(OBJ_EDIT_SMESH)
	$(CCC) $(CFLAGS) -o $@ $(OBJ_EDIT_SMESH) -lm

stat_smesh: $(OBJ_STAT_SMESH)
	$(CCC) $(CFLAGS) -o $@ $(OBJ_STAT_SMESH) -lm

tt_forward: $(OBJ_TT_FORWARD)
	$(CCC) $(CFLAGS) -o $@ $(OBJ_TT_FORWARD) -lm
tt_forward_depend: 
	makedepend -f Makefile $(SRC_TT_FORWARD)

tt_inverse: $(OBJ_TT_INVERSE)
	$(CCC) $(CFLAGS) -o $@ $(OBJ_TT_INVERSE) -lm
tt_inverse_depend: 
	makedepend -f Makefile $(SRC_TT_INVERSE)

