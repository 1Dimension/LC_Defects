
OBJS = Basis.o global.o initial.o lbfgs.o Anderson.o double.o matrix.o 
CC = g++
CFLAGS = -g -O3
LIBS = -lfftw3 -lm

double: $(OBJS)
	$(CC) -o $@ $(OBJS)$ $(LIBS)$

double.o: main.cpp
	$(CC) -o $@ -c $< $(CFLAGS)$

Basis.o: Basis.cpp
	$(CC) -o $@ -c $< $(CFLAGS)$ 

initial.o: initial.cpp
	$(CC) -o $@ -c $< $(CFLAGS)$

global.o: global.cpp Anderson.o
	$(CC) -o $@ -c $< $(CFLAGS)$

Anderson.o: Anderson.cpp matrix.o
	$(CC) -o $@ -c $< $(CFLAGS)$

matrix.o: matrix.cpp 
	$(CC) -o $@ -c $< $(CFLAGS)$

doc:
	doxygen -g
	doxygen Doxyfile

clean:
	rm $(OBJS)
	rm double

.PHINY: clean doc
