DIR_HOME = ../..
DIR_OBJ = $(DIR_HOME)/obj
DIR_LIB = $(DIR_HOME)/lib
DIR_CUR = parallel

include $(DIR_HOME)/Makefile.in

FILE_SRC = parallel.cpp decomp.cpp
FILE_OBJ = $(FILE_SRC:.cpp=.o)
FILE_LIB = lib$(DIR_CUR).a

default:	
	@echo ">>> Compiling "$(DIR_CUR)
	@make compile
	@make library
	@make move_obj

compile:	
	$(CXX) $(CXXFLAGS) -c $(FILE_SRC)

library:	
	$(AR) $(DIR_LIB)/$(FILE_LIB) $(FILE_OBJ)
	$(RL) $(DIR_LIB)/$(FILE_LIB)

move_obj:	
	mv $(FILE_OBJ) $(DIR_OBJ)/
