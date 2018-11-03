OBJECTS =  picard.o AdvancePosition.o AdvanceVelocity.o BCfields.o BCinflow.o BCoutflow.o BCpotentials.o CalcDistance.o CalcEfield.o CalcEpotential.o ParticlesGridify.o ParticlesGridifyAll.o Initial.o ParticlesModify.o ParticlesProduce.o ParticlesTransfer.o Tools.o Dump.o WriteWarnings.o


TARGET  = picard

#FC       = mpiifort
FC       = $(HOME)/mpich2/bin/mpif90
#OFLAGS   = -O3 # -mcmodel=medium -shared-intel
OFLAGS   = -fast
#FCFLAGS  = -I/usr/local/include/openmpi
FCFLAGS  = #-I$(HOME)/mpich2/include
#LDFLAGS  = -L/usr/local/lib/openmpi
LDFLAGS  = -L$(HOME)/mpich2/lib

LD  = $(FC)  # Linker to use. 


#-------------------------------------------------------------
# Make rules
#
$(TARGET): $(OBJECTS)
	$(LD) $(OBJECTS) $(LDFLAGS) -o $(TARGET)

.PHONY: clean
clean_all:
	@rm -rf $(OBJECTS) $(TARGET) core

clean: 
	@rm -rf $(OBJECTS) core

# Rule used to generate object files from source code:
%.o: %.f
	$(FC) $(FCFLAGS) $(OFLAGS) -c $<

%.o: %.f90
	$(FC) $(FCFLAGS) $(OFLAGS) -c $<

%.o: %.f95
	$(FC) $(FCFLAGS) $(OFLAGS) -c $<

%.o: %.f03
	$(FC) $(FCFLAGS) $(OFLAGS) -c $<

%.o: %.F
	$(FC) $(FCFLAGS) $(OFLAGS) -c $<
