# Variable for path to PACS Examples
PACS_ROOT = /mnt/c/Users/samue/OneDrive/Desktop/Samuele/Uni/V_Anno/PACS/pacs-examples #/Examples

# Path for reaching libraries and header provided by pacs-Examples repository
EXAMPLES_INCLUDE = ${PACS_ROOT}/include
EXAMPLES_LIB = ${PACS_ROOT}/lib

# Path for personal header file
MY_INCLUDE = include

# Adding run time search path
ADDITIONAL_PATH = 

# Compiler variables
CXXFLAGS = -std=c++20 -O3 -fPIC 
CPPFLAGS = -DNDEBUG -I${EXAMPLES_INCLUDE} -I${MY_INCLUDE}
LDFLAGS = -L. -Wl,-rpath=$(EXAMPLES_LIB) -L${EXAMPLES_LIB}
LDLIBS = -lmpi -lpacs -lmuparser 

# Parallel compiler variables
P_CXX = mpic++
P_CXXFLAGS = -fopenmp $(CXXFLAGS)
P_CPPFLAGS = $(CPPFLAGS)
P_LDLIBS = -lmuparser -lmpi

# Files
# Source file
SRCS = $(wildcard src/*.cpp)
# Object files from sources
OBJS = $(SRCS:.cpp=.o)
# All headers
HEADS = $(wildcard ${MY_INCLUDE}/*.hpp)
# Executable names
EXEC = main
P_EXEC = main_parallel

.PHONY: all clean parallel

# Compiling routine
all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o $(EXEC) $(LDLIBS)

%.o: %.cpp $(HEADS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@ -I$(MY_INCLUDE)

# Parallel compiling routine
parallel: $(P_EXEC)

$(P_EXEC): $(OBJS)
	$(P_CXX) $(LDFLAGS) $(OBJS) -o $(P_EXEC) $(LDLIBS) $(P_LDLIBS)

%.o: %.cpp $(HEADS)
	$(P_CXX) $(P_CPPFLAGS) $(P_CXXFLAGS) -c $< -o $@ -I$(MY_INCLUDE)

clean:
	$(RM) $(OBJS) $(EXEC) $(P_EXEC)