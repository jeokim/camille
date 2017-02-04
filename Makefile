DIR_HOME = .
DIR_OBJ = $(DIR_HOME)/obj
DIR_BIN = $(DIR_HOME)/bin
DIR_LIB = $(DIR_HOME)/lib

include $(DIR_HOME)/Makefile.in

all:	
	@echo ">>> Compiling and installing "$(NAME_EXE)
	@mkdir -p $(DIR_OBJ) $(DIR_BIN) $(DIR_LIB)
	@for I in $(NAME_LIB); do make -C src/$$I; done
	@make exe
	@make -C src/tools

exe:	
	$(CXX) -o $(DIR_BIN)/$(NAME_EXE) $(DIR_OBJ)/main.o $(NAME_LIB:%=$(DIR_LIB)/lib%.a)

clean:	
	rm -f $(DIR_OBJ)/*.o
	rm -f $(DIR_BIN)/*
	rm -f $(DIR_LIB)/*.a
