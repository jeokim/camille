DIR_HOME = ../../..
DIR_OBJ = $(DIR_HOME)/obj
DIR_LIB = $(DIR_HOME)/lib
DIR_BIN = $(DIR_HOME)/bin
DIR_CUR = check_grid_plot3d

include $(DIR_HOME)/Makefile.in

FILE_SRC = check_grid_plot3d.cpp
FILE_OBJ = $(FILE_SRC:.cpp=.o)
FILE_LIB = lib$(DIR_CUR).a

default:	
	@echo ">>> Compiling "$(DIR_CUR)
	@make compile
	@make library
	@make move_obj
	@make exe

compile:	
	$(CXX) $(CXXFLAGS) -c $(FILE_SRC)

library:	
	$(AR) $(DIR_LIB)/$(FILE_LIB) $(FILE_OBJ)
	$(RL) $(DIR_LIB)/$(FILE_LIB)

move_obj:	
	mv $(FILE_OBJ) $(DIR_OBJ)/

exe:
	$(CXX) -o $(DIR_BIN)/$(DIR_CUR) $(DIR_OBJ)/$(DIR_CUR).o $(DIR_LIB)/$(FILE_LIB) $(NAME_LIB:%=$(DIR_LIB)/lib%.a)

