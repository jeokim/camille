DIR_HOME = ../..
DIR_OBJ = $(DIR_HOME)/obj
DIR_LIB = $(DIR_HOME)/lib
DIR_CUR = io

include $(DIR_HOME)/Makefile.in

FILE_SRC = tecplot.cpp fluent.cpp plot3d.cpp overture.cpp probe.cpp io.cpp
FILE_SRC_C = plot3d_format.c
FILE_OBJ = $(FILE_SRC:.cpp=.o) $(FILE_SRC_C:.c=.o)
FILE_LIB = lib$(DIR_CUR).a

default:	
	@echo ">>> Compiling "$(DIR_CUR)
	@make compile
	@make library
	@make move_obj

compile:	
	$(CC) $(CCFLAGS) -c $(FILE_SRC_C)
	$(CXX) $(CXXFLAGS) -c $(FILE_SRC)

library:	
	$(AR) $(DIR_LIB)/$(FILE_LIB) $(FILE_OBJ)
	$(RL) $(DIR_LIB)/$(FILE_LIB)

move_obj:	
	mv $(FILE_OBJ) $(DIR_OBJ)/
