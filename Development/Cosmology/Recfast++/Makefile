#===============================================================================
# set to C++ compiler that should be used
#===============================================================================
CC = g++

#===============================================================================
# compiler and linker flags
#===============================================================================
CXXFLAGS = -Wall -O2
CXXFLAGSLOC = $(CXXFLAGS) -D RECFASTPPPATH=\"$(PWD)/\"
LXXFLAGS =

#===============================================================================
# object-files
#===============================================================================
OBJS = ./src/cosmology.Recfast.o \
	   ./src/evalode.Recfast.o \
	   ./src/recombination.Recfast.o \
	   ./src/ODE_solver.Recfast.o \
	   ./src/DM_annihilation.Recfast.o \
	   ./src/Magnetic_fields.Recfast.o \
	   ./src/Variation_constants.Recfast.o \
	   ./src/Rec_corrs_CT.Recfast.o

#===============================================================================
# program	 
#===============================================================================
all: Recfast++

Recfast++:  $(OBJS) main.o 
			 @echo "Linking..."
			 $(CC) $(LXXFLAGS) $(OBJS) main.o -o Recfast++

pub:
	make -f Makefile.pub all

pubtar:
	make -f Makefile.pub tarball

clean:
	rm -f ./src/*.o

tidy:
	rm -f ./src/*.o ./src/*.o~ 
	rm -f *.o *~ Recfast++

wipeDS:
	find . -type f -name \.DS_Store -print | xargs rm

#===============================================================================
# rules
#===============================================================================
INC_PATH = -I./include
		   
.cpp.o:
	@echo "Producing object-file $@"
	$(CC) $(CXXFLAGSLOC) $(INC_PATH) -c $< -o $@
	
#===============================================================================
#===============================================================================
