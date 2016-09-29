
TARGET = anal00

OBJS = Globals.o Mathematics.o DumpManipulation.o \
       mod_InertiaMoment.o \
			 Analysis.o ReadOrder.o Main.o
MODS = globals.mod mathematics.mod dumpmanipulation.mod \
       mod_inertiamoment.mod \
			 analysis.mod readorder.mod

FC = mpifort
FCFLAGS = -O3

.SUFFIXES : .o .f90
.f90.o:
	${FC} -c $<

${TARGET} : ${OBJS}
	${FC} -o $@ ${OBJS} ${FCLAGS}

.PHONY: clean
clean:
	rm ${TARGET} ${OBJS} ${MODS}
