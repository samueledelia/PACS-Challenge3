CXX      = mpicxx
CXXFLAGS = -std=c++20 -fopenmp 
CPPFLAGS = -O0 -Wall -I./include -Wno-conversion-null -Wno-deprecated-declarations 
SOURCE_DIR=./src/ #set the location of source file
VPATH=$(SOURCE_DIR)
SRCS=$(join $(dir $(SOURCE_DIR)),$(notdir *.cpp)) #source files
OBJS= $(SRCS:%.cpp=%.o) #object files

EXEC= main #I want one executable called "main"
.phony= clean 
.DEFAULT_GOAL = all 
all: $(EXEC)

$(OBJS): $(SRCS) 
	$(CXX) -g -c $(CXXFLAGS) $(CPPFLAGS) $(SRCS)
	mv *.o ./src
	mkdir -p results
$(EXEC): $(OBJS)
	$(CXX) -g $(CXXFLAGS) $(OBJS) -o $(EXEC)

clean:
	$(RM) *.o
	$(RM) $(OBJS)
	$(RM) $(EXEC)
	$(RM) -r ./doc/html ./doc/latex