CC=gcc
CFLAGS=-O3 -fopenmp 
LIBFLAGS=-lm

OBJDIR = obj

_OBJ = args.o data.o setup.o vtk.o boundary.o airfoil.o
OBJ = $(patsubst %,$(OBJDIR)/%,$(_OBJ))

.PHONY: directories

all: directories airfoil

obj/%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS) 

airfoil: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBFLAGS) 

clean:
	rm -Rf $(OBJDIR)
	rm -f airfoil

configure:
	export OMP_NUM_THREADS=8

directories: $(OBJDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)

