

#
# Olliptic makefile
#


# ---------------------------------------------------------------
# Set exe name
# ---------------------------------------------------------------

EXE = PN3BP

# ---------------------------------------------------------------
# Set compiler options
# ---------------------------------------------------------------

MPI_COMPILE_FLAGS = $(shell mpicc --showme:compile)
MPI_LINK_FLAGS = $(shell mpicc --showme:link)

GSL_COMPILE_FLAGS = $(shell gsl-config --cflags)
GSL_LINK_FLAGS = $(shell gsl-config --libs)


ifeq ($(shell uname),Darwin)
BOOST_LIB_PATH=/usr/local/lib
BOOST_INC_PATH=/usr/local/include/boost

HDF5_LIB_PATH=/usr/local/lib
HDF5_INC_PATH=/usr/local/include

CPP	= mpic++ 
CC	= mpicc  
CPPFLAGS =  -I$(HDF5_INC_PATH)  -I$(BOOST_INC_PATH)  $(MPI_COMPILE_FLAGS) $(GSL_COMPILE_FLAGS) 
LDFLAGS =   -fvisibility=default\
			$(BOOST_LIB_PATH)/libboost_system.a\
			$(BOOST_LIB_PATH)/libboost_log_setup.a\
			$(BOOST_LIB_PATH)/libboost_log.a\
			$(BOOST_LIB_PATH)/libboost_date_time.a\
			$(BOOST_LIB_PATH)/libboost_thread.a\
			$(BOOST_LIB_PATH)/libboost_filesystem.a \
			$(BOOST_LIB_PATH)/libboost_program_options.a\
			$(HDF5_LIB_PATH)/libhdf5_hl.a\
			$(HDF5_LIB_PATH)/libhdf5.a -lz -ldl -lm\
			$(MPI_LINK_FLAGS) $(GSL_LINK_FLAGS)
endif
ifeq ($(shell uname),Linux)

BOOST_LIB_PATH=/usr/lib/x86_64-linux-gnu
BOOST_INC_PATH=/usr/include/boost

HDF5_LIB_PATH=/usr/lib/x86_64-linux-gnu/hdf5/openmpi/
HDF5_INC_PATH=/usr/include/hdf5/openmpi/

CPP	= mpic++
CPPFLAGS =  -I$(HDF5_INC_PATH)  -I$(BOOST_INC_PATH)  $(MPI_COMPILE_FLAGS)  $(GSL_COMPILE_FLAGS)  
LDFLAGS =   -fvisibility=default\
			$(BOOST_LIB_PATH)/libboost_system.a\
			$(BOOST_LIB_PATH)/libboost_log_setup.a\
			$(BOOST_LIB_PATH)/libboost_log.a\
			$(BOOST_LIB_PATH)/libboost_date_time.a\
			$(BOOST_LIB_PATH)/libboost_thread.a\
			$(BOOST_LIB_PATH)/libboost_filesystem.a \
			$(BOOST_LIB_PATH)/libboost_program_options.a\
			$(HDF5_LIB_PATH)/libhdf5_hl.a\
			$(HDF5_LIB_PATH)/libhdf5.a -lz -lsz -ldl -lm\
			$(MPI_LINK_FLAGS)  $(GSL_LINK_FLAGS)
endif

AR	= ar
MAKE	= make

OPTFLAGS = -O2


debug: OPTFLAGS = -O0 -g -pg -Wall -DDEBUG




# ---------------------------------------------------------------
# Directories
# ---------------------------------------------------------------
BASE = $(shell /bin/pwd)
OBJD = $(BASE)/lib/obj/
LIBD = $(BASE)/lib
EXED = $(BASE)/exe




# ---------------------------------------------------------------
# Modules (add a new module here)
# ---------------------------------------------------------------
libpaths = src/Modules/Evolution
libpaths += src/Modules/Initial_data
libpaths += src/Modules/Output
libpaths += src/Modules/Utils
libpaths += src/Modules/Analysis
# ---------------------------------------------------------------

PN1 = -DodePN1
##PN2 = -DodePN2
#PN2_5 = -DodePN2_5
##PNSlo = -DodePNSlo
##PNSnlo = -DodePNSnlo

CPPFLAGS += -Wno-deprecated -std=c++11 $(OPTFLAGS)



libdirs = $(dir $(libpaths))
libnames = $(notdir $(libpaths))
libnames_a = $(addsuffix .a,$(libnames))
libnames_h = $(addsuffix .h,$(notdir $(libpaths)))
libnames_h := `echo $(libnames_h) | tr A-Z a-z` 

liblist := $(foreach libname,$(libnames_a),$(libname_a))

HEADERDIR := $(addsuffix /includes,$(addprefix -I$(BASE)/,$(libpaths)))


LIBS+=$(notdir $(libpaths))

export 

new: new_module_name = `echo $(MODULE) | cut -c 1 | tr a-z A-Z``echo $(MODULE) | cut -c 2- | tr A-Z a-z`
new: new_module_file = `echo $(MODULE) | tr A-Z a-z`
new: new_module_capitals = `echo $(MODULE) | tr a-z A-Z`

new: new_class_name = `echo $(CLASS) | tr A-Z a-z | sed -r 's/(_)([a-z])/\U\2/g' `
new: new_class_directory = `echo $(CLASS) | cut -c 1 | tr a-z A-Z``echo $(CLASS) | cut -c 2- | tr A-Z a-z`
new: new_class_file = `echo $(CLASS) | tr A-Z a-z`
new: new_class_capitals = `echo $(CLASS) | tr a-z A-Z`



# ---------------------------------------------------------------
# target 1 : all
# ---------------------------------------------------------------

all: intro setdirs makelibs makemain link done


intro:
	@echo ''
	@echo '===> Making Olliptic ...'	
	@echo ''	

setdirs:
	@mkdir -p $(EXED) $(LIBD) $(OBJD) 
	@for X in $(LIBS); do mkdir -p lib/obj/$$X; done
	@rm -rf $(BASE)/src/main.h

makelibs:
	@echo ''
	@echo '===> Building libs'
	@echo ''
	@echo "-----------------------------------------------------------------------------------"
	@for X in $(libpaths); do $(MAKE) --no-print-directory -s -C $$X; done
	@for X in $(libnames_h); do echo "#include \"$$X\"" >> $(BASE)/src/main.h; done
	@/bin/cat src/skeleton/main_natural_docs.txt >> $(BASE)/src/main.h
	@echo "...done"
	@echo "-----------------------------------------------------------------------------------"

makemain:
	@echo ''
	@echo '===> Building main'
	@echo ''
	@cd $(LIBD); $(CPP) $(HEADERDIR) -I$(BASE)/src/ -o main.o -c $(BASE)/src/main.cc $(CPPFLAGS)


link: 
	@echo ''
	@echo '===> Building the executable'
	@echo ''	
	@cd $(LIBD); $(CPP) $(HEADERDIR) $(VTK_INC) -o $(EXED)/$(EXE).x $(LIBD)/main.o  $(libnames_a) $(VTK_LIB) $(CPPFLAGS) $(LDFLAGS) -w




done:
	@echo ''
	@echo '============ All Done! ============'
	@echo ' '


# ---------------------------------------------------------------
# Target 2 : debug
# ---------------------------------------------------------------

debug: intro clean setdirs makelibs makemain link
	@echo ''
	@echo 'DEBUG:'
	@echo ' gdb [exe] [corefile]'
	@echo ''
	@echo 'Remember the settings for the core output !!!'
	@echo '  (sh/bash)    ulimit -c unlimited'
	@echo '  (tcsh/csh)   limit coredumpsize unlimited'
	@echo ' '
	@echo 'PROFILE:'
	@echo ' time [exe]'
	@echo ' gprof -b [exe] gmon.out'
	@echo ' '


# ---------------------------------------------------------------
# Target 3 : clean tmp files
# ---------------------------------------------------------------

clean: intro
	@rm -f $(BASE)/lib/*.a
	@rm -f $(BASE)/lib/*.o
	@rm -rf $(BASE)/lib/obj/*
	@rm -f $(BASE)/*~
	@rm -f $(BASE)/src/*~
	@rm -f $(BASE)/*.dat
	@rm -f $(BASE)/*.xg
	@echo " === clean done ! === "

# ---------------------------------------------------------------
# Target 4 : delete
# ---------------------------------------------------------------

delete: intro
	@rm -f $(BASE)/lib/*.a
	@rm -f $(BASE)/lib/*.o
	@rm -rf $(BASE)/lib/obj/*
	@rm -f $(BASE)/*~
	@rm -f $(BASE)/src/*~
	@rm -f $(BASE)/*\#
	@rm -f $(BASE)/*.x
	@rm -f $(BASE)/*.out
	@rm -f $(BASE)/*.dat
	@rm -f $(BASE)/*.xg
	@echo " === delete done ! === "

# ---------------------------------------------------------------
# target 5 : help
# ---------------------------------------------------------------

help:	intro
	@echo '	'
	@echo '	---------------------------'
	@echo '	Makefile - HELP'
	@echo '	---------------------------'
	@echo '	'
	@echo '	make				: compile with debug options'
	@echo ' '
	@echo '	make run			: compile with optimizations'
	@echo '	'
	@echo '	make docs			: make the documentation' 
	@echo '					  needs http://www.naturaldocs.org'
	@echo '	'
	@echo '	make new MODULE=new_module	: make a directory and files' 
	@echo '					  for  a new module'
	@echo '	'
	@echo '	make new CLASS=new_class	: make a directory and files' 
	@echo '					  for  a new class'
	@echo '	'
	@echo '	make clean			: clean object, ~ and data files'
	@echo ' '
	@echo '	make delete			: delete object, ~, data and exe file'
	@echo ' '


docs:
	@NaturalDocs -i src/ -o HTML doc/Olliptic/ -p doc/NaturalDocsProject/


# ---------------------------------------------------------------
# target 6 : make new module or class
# ---------------------------------------------------------------




new:
ifeq ($(MODULE)$(CLASS),)
	@echo '	'
	@echo 'Stop: Missing name of module or class!!!'
	@$(MAKE)  --no-print-directory -s help
else
ifneq ($(MODULE),)
	@$(MAKE)  --no-print-directory -s intro
	@echo ''
	@echo '===> Making new module :' $(new_module_name)
	@echo ''
	@mkdir $(BASE)/src/Modules/$(new_module_name)
	@mkdir $(BASE)/src/Modules/$(new_module_name)/src
	@mkdir $(BASE)/src/Modules/$(new_module_name)/includes
	@echo 'NAME := '$(new_module_name) > $(BASE)/src/Modules/$(new_module_name)/Makefile
	@echo 'OBJS := '$(new_module_file)'.o' >> $(BASE)/src/Modules/$(new_module_name)/Makefile
	@echo 'HSF := '$(new_module_file)'.h ' >> $(BASE)/src/Modules/$(new_module_name)/Makefile
	@echo 'LANGUAGE := CPP' >> $(BASE)/src/Modules/$(new_module_name)/Makefile
	@echo 'include $$(BASE)/Makefile.subdirs' >> $(BASE)/src/Modules/$(new_module_name)/Makefile
	@echo '\n\n// '$(new_module_file)'.h' >> $(BASE)/src/Modules/$(new_module_name)/includes/$(new_module_file).h
	@/bin/cat src/skeleton/gnu_license.txt >> $(BASE)/src/Modules/$(new_module_name)/includes/$(new_module_file).h
	@echo '\n\n#ifndef '$(new_module_capitals)'_H' >> $(BASE)/src/Modules/$(new_module_name)/includes/$(new_module_file).h
	@echo '\n#define '$(new_module_capitals)'_H' >> $(BASE)/src/Modules/$(new_module_name)/includes/$(new_module_file).h
	@echo '\n\n#endif' >> $(BASE)/src/Modules/$(new_module_name)/includes/$(new_module_file).h
	@echo '\n\n// '$(new_module_file)'.cc' >> $(BASE)/src/Modules/$(new_module_name)/src/$(new_module_file).cc
	@/bin/cat src/skeleton/gnu_license.txt >> $(BASE)/src/Modules/$(new_module_name)/src/$(new_module_file).cc
	@echo '\n\n#include "'$(new_module_file)'.h"' >> $(BASE)/src/Modules/$(new_module_name)/src/$(new_module_file).cc
	@$(MAKE)  --no-print-directory -s newclose

endif
ifneq ($(CLASS),)
	@$(MAKE)  --no-print-directory -s intro
	@echo ''
	@echo '===> Making new class :' $(new_class_name)
	@echo ''
	@mkdir $(BASE)/src/Modules/$(new_class_directory)
	@mkdir $(BASE)/src/Modules/$(new_class_directory)/src
	@mkdir $(BASE)/src/Modules/$(new_class_directory)/includes
	@echo 'NAME := '$(new_class_directory) > $(BASE)/src/Modules/$(new_class_directory)/Makefile
	@echo 'OBJS := '$(new_class_file)'.o' >> $(BASE)/src/Modules/$(new_class_directory)/Makefile
	@echo 'HSF := '$(new_class_file)'.h ' >> $(BASE)/src/Modules/$(new_class_directory)/Makefile
	@echo 'LANGUAGE := CPP' >> $(BASE)/src/Modules/$(new_class_directory)/Makefile
	@echo 'include $$(BASE)/Makefile.subdirs' >> $(BASE)/src/Modules/$(new_class_directory)/Makefile
	@echo '\n\n// '$(new_class_file)'.h' >> $(BASE)/src/Modules/$(new_class_directory)/includes/$(new_class_file).h
	@/bin/cat src/skeleton/gnu_license.txt >> $(BASE)/src/Modules/$(new_class_directory)/includes/$(new_class_file).h
	@echo '\n\n#ifndef '$(new_class_capitals)'_H' >> $(BASE)/src/Modules/$(new_class_directory)/includes/$(new_class_file).h
	@echo '\n#define '$(new_class_capitals)'_H' >> $(BASE)/src/Modules/$(new_class_directory)/includes/$(new_class_file).h
	@echo '\n\n\nclass '$(new_class_name)'\n{' >> $(BASE)/src/Modules/$(new_class_directory)/includes/$(new_class_file).h
	@echo '\n\n public:\n\n '$(new_class_name)'(){};' >> $(BASE)/src/Modules/$(new_class_directory)/includes/$(new_class_file).h
	@echo '\n ~'$(new_class_name)'(){};\n\n\n};' >> $(BASE)/src/Modules/$(new_class_directory)/includes/$(new_class_file).h
	@echo '\n\n#endif' >> $(BASE)/src/Modules/$(new_class_directory)/includes/$(new_class_file).h
	@echo '\n\n// '$(new_class_file)'.cc' >> $(BASE)/src/Modules/$(new_class_directory)/src/$(new_class_file).cc
	@/bin/cat src/skeleton/gnu_license.txt >> $(BASE)/src/Modules/$(new_class_directory)/src/$(new_class_file).cc
	@echo '\n\n#include "'$(new_class_file)'.h"' >> $(BASE)/src/Modules/$(new_class_directory)/src/$(new_class_file).cc
	@$(MAKE)  --no-print-directory -s newclose

endif
endif


newclose:
	@echo '	'
	@echo '	---------------------------------------------------------'
	@echo '	Makefile - WARNING'
	@echo '	---------------------------------------------------------'
	@echo '	'
	@echo '	Remember: Add a new line in this Makefile:'
	@echo ' '
	@echo '	libpaths += src/Modules/'$(new_module_name)$(new_class_directory)
	@echo '	'
	@echo '	and add the directory to git repository'
	@echo '	'
	@echo '	git add src/Modules/'$(new_module_name)$(new_class_directory)
	@echo '	'
	@echo '	---------------------------------------------------------'
	@echo '	'
	@echo '	'




