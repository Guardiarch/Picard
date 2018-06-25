OBJECTS =  picard.o AdvancePosition.o AdvanceVelocity.o BCfields.o BCinflow.o BCoutflow.o BCpotentials.o CalcDistance.o CalcEfield.o CalcEpotential.o ParticlesGridify.o ParticlesGridifyAll.o ParticlesInitialize.o ParticlesModify.o ParticlesProduce.o ParticlesTransfer.o Tools.o Write.o WriteWarnings.o


TARGET  = picard

FC       = mpiifort
OFLAGS   = -O3 # -mcmodel=medium -shared-intel
FCFLAGS  = -I/usr/local/include/openmpi
LDFLAGS  = -L/usr/local/lib/openmpi

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
