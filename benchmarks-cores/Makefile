
CASE := test-tiny.cfg
#CASE := test-tiny-long.inp
#CASE := test-big.inp
FRAME := 4
EXE := ../dynearthsol2d
OMP := 1

# Use single cpu for omp to avoid round-off errors
RUN_TEST = rm -f benchmark.*; OMP_NUM_THREADS=${OMP} ${EXE} ${CASE}

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

.PHONY: all set cmp store
