#Makefile
#$gls_directory should be set by the user

COMPILER = g++ 
CFLAGS = -Wall -O3 -g
LDFLAGS = -lm
EXECUTABLE = AMD
OBJECT     = AMD.o countBp.o  Correlation.o stdafx.o OtherTree.o

$(EXECUTABLE): $(OBJECT)
	$(COMPILER) $(LIBS)  $(OBJECT) $(LDFLAGS) -o $(EXECUTABLE)

%.o:    %.c
	$(COMPILER) $(INCL) -o $*.o -c $*.cpp

clean:
	rm -f *.o
