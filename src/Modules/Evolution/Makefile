NAME := Evolution
OBJS := evolution.o
HSF := evolution.h
LANGUAGE := CPP

OBJS += Nnp2d2NS.o Nnp3d2NS.o  Nnp3d3NS.o Nnp2d3S.o Nnp3d3S.o 

ifdef PN1
OBJS += 1PNnp2d2NS.o  1PNnp3d2NS.o 1PNnp3d3NS.o 1PNnp2d3S.o
endif

ifdef PN2
OBJS += 2PNnp2d2NS.o 2PNnp3d2NS.o  2PNnp3d3NS.o 
endif

ifdef PN2_5
OBJS += 2_5PNnp2d2NS.o  2_5PNnp3d2NS.o 2_5PNnp3d3NS.o 2_5PNnp2d3S.o
endif

ifdef PNSlo
OBJS += hnp2d3S.o hnp3d3S.o Nnp2d3S.o Nnp3d3S.o 1PNnp3d3S.o 2PNnp2d3S.o 2PNnp3d3S.o 2_5PNnp2d3S.o
endif

ifdef PNSnlo
OBJS += SOnloPNnp2d3S.o #SOnloPNnp3d3S.o #SSloPNnp2d3.o S2loPNnp2d3.o SSloPNnp3d3.o S2loPNnp3d3.o
endif

include $(BASE)/Makefile.subdirs
