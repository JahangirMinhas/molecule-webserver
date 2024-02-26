CC = clang
CFLAGS = -Wall -std=c99 -pedantic

all:  _molecule.so

clean:  
	rm -f *.o *.so molecule_wrap.c _molecule.so molecule.py

mol.o: mol.c mol.h
	$(CC) $(CFLAGS) -c mol.c -fPIC -o mol.o

libmol.so: mol.o
	$(CC) mol.o -shared -o libmol.so

molecule_wrap.c molecule.py: mol.o molecule.i 
	swig -python molecule.i

molecule_wrap.o: mol.h molecule_wrap.c
	$(CC) $(CFLAGS) -c molecule_wrap.c -fPIC -o molecule_wrap.o

_molecule.so: libmol.so molecule_wrap.o
	$(CC) -shared molecule_wrap.o -lmol -o _molecule.so
