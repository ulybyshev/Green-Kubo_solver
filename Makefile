#this is the comment
#variables defined as NAME:=  (NAME can be any) and addressed as $(NAME) 
#if you want to recompile all:  "make clean"  then "make all"
#if you want to recompile only complex class without recompilation of all other classes: "make complex.o"  then "make example"

#compilator name (for example could be mpicc instead)
CC:=g++

#LIBRARIES (can also include path to these libraries as -L/SOME_PATH)
LIB:=-L/Users/ulybyshev/Calculations/gsl/lib  -lm  -lgsl  -lgslcblas

#Optimization (-O1, -O2 or -O3)
OPT:=-Wall -O3


#directory for .o files
OBJECT_DIR:=./object

#directory for executables
BIN_DIR:=./bin

#directory for .cpp files
SOURCE_DIR:=./source

#directory for .h files
INCLUDE_DIR:=-I./include  -I/Users/ulybyshev/Calculations/gsl/include


#List of source files
SRC_FILES := $(wildcard $(addsuffix /*.cpp, $(SOURCE_DIR)))
OBJ_FILES := $(patsubst %.cpp, $(OBJECT_DIR)/%.o, $(notdir $(SRC_FILES)))
SRC_FILES_NAMES := $(patsubst $(SOURCE_DIR)/%, %, $(SRC_FILES))
OBJ_FILES_NAMES := $(patsubst %.cpp, %.o, $(SRC_FILES_NAMES))
FILES_NAMES := $(patsubst %.o, %, $(OBJ_FILES_NAMES))

all: gk_solver

gk_solver:  $(OBJ_FILES_NAMES)
	$(CC)  $(OPT)  $(OBJ_FILES)  -o $(BIN_DIR)/gk_solver   $(LIB)

#two equivalent ways: semicolumn of hard-tab (two tabs in mcedit) in the next line

define app_compile_template
 $(1)_NAME = $$(patsubst %.o,%, $(1))
 $(1)_CPP_NAME = $$(SOURCE_DIR)/$$($(1)_NAME).cpp
 $(1)_OBJ_NAME = $$(OBJECT_DIR)/$$($(1)_NAME).o
$(1): $$($(1)_CPP_NAME)
	$$(CC) $$(OPT) -c $$($(1)_CPP_NAME)  -o $$($(1)_OBJ_NAME)  $$(INCLUDE_DIR)
endef

$(foreach app, $(OBJ_FILES_NAMES), $(eval $(call app_compile_template,$(app))))

#cleaning
clean:
	rm -rf $(OBJECT_DIR)/*.o $(BIN_DIR)/gk_solver

list_object_files:
	echo $(OBJ_FILES_NAMES)

