
CASE := test-tiny.cfg
#CASE := test-tiny-long.cfg
#CASE := test-big.cfg
#CASE := test-topo.cfg
FRAME := 4
EXE := ../dynearthsol2d
OMP := 1
TFLAG = #-v
DEBUG := 0
ifeq ($(DEBUG), 1)
	DEBUGFLAG = valgrind --tool=helgrind --log-file=valgrind-out.txt -v
endif

# Use single cpu for omp to avoid round-off errors
RUN_TEST = rm -f benchmark.*; OMP_NUM_THREADS=${OMP}; time ${TFLAG} ${DEBUGFLAG} ${EXE} ${CASE}; echo

ORIG = orig-${CASE}

all: cmp

set:
	# run test, store test output to old/
	${RUN_TEST}
	${MAKE} store

cmp:
	${RUN_TEST}
	python3 compare.py ${ORIG} ${FRAME}

store:
	rm -rf ${ORIG}
	mkdir -p ${ORIG}
	mv benchmark.* ${ORIG}
	cp ../snapshot.diff ${ORIG}

vtk:
	python ../2vtk.py -m  benchmark

.PHONY: all set cmp store
