# LIFE makefile

# Compiler command
CC=clang++
CFLAGS=-O3 -std=c++0x -fopenmp -Wall -Wextra

# Executable
EXE=LIFE

# Location of source, header and object files
DIR=.
SDIR=$(DIR)/src
HDIR=$(DIR)/inc
ODIR=$(DIR)/obj

# Get the sources and object files
SRCS:=$(wildcard $(SDIR)/*.cpp)
OBJS:=$(addprefix $(ODIR)/,$(notdir $(SRCS:.cpp=.o)))

# Include and library files
INC+=-I/opt/homebrew/include/
LIB=-llapack -lboost_system -lboost_filesystem
LIB+=-L/opt/homebrew/lib
LIB+=-L/opt/homebrew/opt/libomp/lib
LIB+=-L/opt/homebrew/opt/llvm/lib

# Build LIFE
$(EXE): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LIB)

# Build object files
$(OBJS): $(ODIR)/%.o : $(SDIR)/%.cpp
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<

# Clean the project
.PHONY: clean
clean:
	rm -rf $(EXE) $(ODIR) Results *.out && mkdir $(ODIR)
