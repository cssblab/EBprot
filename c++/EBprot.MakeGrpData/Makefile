
UNAME := $(shell uname) ## will return either Darwin or Linux

ifeq ($(UNAME), Darwin )
	INCLUDES = -I/usr/local/include
	LIBPATH = -L/usr/local/lib
endif

ifeq ($(UNAME), Linux)
	INCLUDES = -I/usr/include/gsl
	LIBPATH = -L/usr/lib64
endif

LIBS = -lgsl -lgslcblas
CPPFLAGS = -O3 -std=c++0x
OBJS = main.o global.o
CPP = g++

all: EBprot.MakeGrpData
EBprot.MakeGrpData: $(OBJS)
	$(CPP) $(CPPFLAGS) $(INCLUDES) $(OBJS) $(LIBPATH) $(LIBS) -o ../bin/EBprot.MakeGrpData
	rm -f *.o
	@echo
	@echo Done
	@echo Execute ./bin/EBprot.MakeGrpData to get command line usage statement
	@echo


main.o:
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c main.cpp

global.o:
	$(CPP) $(CPPFLAGS) $(INCLUDES) -c global.cpp

clean:
	rm -f *.o
