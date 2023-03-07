CC = g++
LIBS	= -lgurobi_c++ -lgurobi81 -lm #change lgurobi## to version you are using
LIBDIRS	= -L /home/stephen/code/gurobi811/linux64/lib/ #change to directory where you have gurobi header file 
G_INC	= -I /home/stephen/code/gurobi811/linux64/include/ #Change directory to fit where you have gurobi and saucy header files
S_INC 	= -I /home/stephen/code/saucy/
FLAGS	= -std=c++11 

PROG1   = $(PWD)/src/essential-DS-saucy.cpp
TARGET1 = $(PWD)/Classifier_C

PROG2	= $(PWD)/src/essential-DS-no-orbits.cpp
TARGET1 = $(PWD)/Classifier_B

PROG3   = $(PWD)/essential-DS-nacher
TARGET3	= $(PWD)/Classifier_A


SOURCES  = $(shell find $(SRCDIR) -xtype f -name "*".$(SRCEXT) ! -name main* )
OBJECTS  = $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

#$(info $(SOURCES))
#$(info $(OBJECTS))

.SUFFIXES: .cpp .o
.cpp.o:
	$(CC) $(CFLAGS) -c 

all: $(TARGET1) $(TARGET2) $(TARGET3)



$(TARGET1): $(OBJECTS) $(PROG1)
	$(CC) -o $@ $^  $(FLAGS) $(INC) $(LIBDIRS) $(LIBS)

$(TARGET2): $(OBJECTS) $(PROG2)
	$(CC) -o $@ $^  $(FLAGS) $(INC) $(LIBDIRS) $(LIBS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)  #igraph
	@mkdir -p $(BUILDDIR)
	$(CC) -c -o $@ $< $(FLAGS) $(INC) $(LIBDIRS) $(LIBS)


clean:
	$(RM) -r $(BUILDDIR) $(TARGET)
	#$(RM) $(LIBDIR)/lib*
	#$(MAKE) --directory=$(IGRAPHDIR) clean





