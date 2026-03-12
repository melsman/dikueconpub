FUTHARK=futhark

#FUTHARK=$(HOME)/.local/bin/futhark

.PHONY: all
all: bench_bellmanJ.pdf bench_newton_single.pdf bench_linear_solver.pdf

.PHONY: clean
clean:
	rm -rf *~ *.c *.exe *.tm *.tmavg *.out trmodel_benchmark.fut.*

.PHONY: futversion
futversion:
	$(FUTHARK) --version

.PHONY: test
test:
	$(FUTHARK) test trmodel_benchmark.fut

# ---- compile benchmark executable for each backend ----

trmodel_benchmark-c.exe: trmodel_benchmark.fut
	$(FUTHARK) c -o $@ $<

trmodel_benchmark-m.exe: trmodel_benchmark.fut
	$(FUTHARK) multicore -o $@ $<

trmodel_benchmark-o.exe: trmodel_benchmark.fut
	$(FUTHARK) opencl -o $@ $<

# ---- raw timing -> averaged timing ----

%.tmavg: %.tm
	cat $< | python3 -c "import sys; nums = [int(n) for n in sys.stdin]; print(sum(nums)/len(nums)/1000000)" > $@

# =========================================================
# BellmanJ benchmark
# stem format: n-c-Ax-N
# =========================================================

bench_bellmanJ-%-C.tm: trmodel_benchmark-c.exe
	echo "$*" | sed 's/-/ /g' | ./$< -e bench_bellmanJ -r 5 -t $@ > trmodel_benchmark-$*-bellmanJ-c.out

bench_bellmanJ-%-M.tm: trmodel_benchmark-m.exe
	echo "$*" | sed 's/-/ /g' | ./$< -e bench_bellmanJ -r 5 -t $@ > trmodel_benchmark-$*-bellmanJ-m.out

bench_bellmanJ-%-O.tm: trmodel_benchmark-o.exe
	echo "$*" | sed 's/-/ /g' | ./$< -d AMD -e bench_bellmanJ -r 5 -t $@ > trmodel_benchmark-$*-bellmanJ-o.out

BELLMANJ_DATA=2-20-5-100 2-20-10-100 2-20-15-100 2-20-20-100 2-20-25-100 2-20-30-100 2-20-35-100 2-20-40-100

BELLMANJ_C_FILES=$(BELLMANJ_DATA:%=bench_bellmanJ-%-C.tmavg)
BELLMANJ_M_FILES=$(BELLMANJ_DATA:%=bench_bellmanJ-%-M.tmavg)
BELLMANJ_O_FILES=$(BELLMANJ_DATA:%=bench_bellmanJ-%-O.tmavg)

BELLMANJ_FILES=$(BELLMANJ_M_FILES) $(BELLMANJ_C_FILES)
#BELLMANJ_FILES=$(BELLMANJ_O_FILES) $(BELLMANJ_M_FILES) $(BELLMANJ_C_FILES)

bench_bellmanJ.dat: $(BELLMANJ_FILES)
	rm -f $@
	for n in $^ ; do \
             /bin/echo -n $$n | sed 's/[a-z_]*//g' | sed 's/[-.]/ /g' >> $@ ; \
             cat $$n >> $@ ; \
        done

.PHONY: bench_bellmanJ_all
bench_bellmanJ_all: $(BELLMANJ_FILES)

bench_bellmanJ.pdf: plot_benchJ.gp bench_bellmanJ.dat matlab_benchJ.dat
	gnuplot -e "load \"$<\""

.PHONY: bench_bellmanJ_plot
bench_bellmanJ_plot: bench_bellmanJ.pdf

# =========================================================
# Newton single benchmark
# stem format: n-c-Ax-iter-tol
# =========================================================

bench_newton_single-%-C.tm: trmodel_benchmark-c.exe
	echo "$*" | sed 's/-/ /g' | ./$< -e bench_newton_single -r 5 -t $@ > trmodel_benchmark-$*-newton-c.out

bench_newton_single-%-M.tm: trmodel_benchmark-m.exe
	echo "$*" | sed 's/-/ /g' | ./$< -e bench_newton_single -r 5 -t $@ > trmodel_benchmark-$*-newton-m.out

bench_newton_single-%-O.tm: trmodel_benchmark-o.exe
	echo "$*" | sed 's/-/ /g' | ./$< -d AMD -e bench_newton_single -r 5 -t $@ > trmodel_benchmark-$*-newton-o.out

NEWTON_DATA=1-5-20-60-0 1-10-20-60-0 1-15-20-60-0 1-20-20-60-0 1-25-20-60-0 1-30-20-60-0 1-35-20-60-0 1-40-20-60-0 1-45-20-60-0 1-50-20-60-0

NEWTON_C_FILES=$(NEWTON_DATA:%=bench_newton_single-%-C.tmavg)
NEWTON_M_FILES=$(NEWTON_DATA:%=bench_newton_single-%-M.tmavg)
NEWTON_O_FILES=$(NEWTON_DATA:%=bench_newton_single-%-O.tmavg)

NEWTON_FILES=$(NEWTON_M_FILES) $(NEWTON_C_FILES)
#NEWTON_FILES=$(NEWTON_O_FILES) $(NEWTON_M_FILES) $(NEWTON_C_FILES)

bench_newton_single.dat: $(NEWTON_FILES)
	rm -f $@
	for n in $^ ; do \
             /bin/echo -n $$n | sed 's/[a-z_]*//g' | sed 's/[-.]/ /g' >> $@ ; \
             cat $$n >> $@ ; \
        done

.PHONY: bench_newton_single_all
bench_newton_single_all: $(NEWTON_FILES)

bench_newton_single.pdf: plot_newton_single.gp bench_newton_single.dat matlat_newton_bench.dat
	gnuplot -e "load \"$<\""

.PHONY: bench_newton_single_plot
bench_newton_single_plot: bench_newton_single.pdf

# =========================================================
# Linear solver benchmark
# stem format: n-c-Ax
# =========================================================

bench_linear_solver-%-C.tm: trmodel_benchmark-c.exe
	echo "$*" | sed 's/-/ /g' | ./$< -e bench_linear_solver -r 5 -t $@ > trmodel_benchmark-$*-linear-c.out

bench_linear_solver-%-M.tm: trmodel_benchmark-m.exe
	echo "$*" | sed 's/-/ /g' | ./$< -e bench_linear_solver -r 5 -t $@ > trmodel_benchmark-$*-linear-m.out

bench_linear_solver-%-O.tm: trmodel_benchmark-o.exe
	echo "$*" | sed 's/-/ /g' | ./$< -d AMD -e bench_linear_solver -r 5 -t $@ > trmodel_benchmark-$*-linear-o.out

LINEAR_DATA=1-5-25 1-10-25 1-15-25 1-20-25 1-25-25 1-30-25 1-35-25 1-40-25 1-45-25 1-50-25

LINEAR_C_FILES=$(LINEAR_DATA:%=bench_linear_solver-%-C.tmavg)
LINEAR_M_FILES=$(LINEAR_DATA:%=bench_linear_solver-%-M.tmavg)
LINEAR_O_FILES=$(LINEAR_DATA:%=bench_linear_solver-%-O.tmavg)

LINEAR_FILES=$(LINEAR_M_FILES) $(LINEAR_C_FILES)
#LINEAR_FILES=$(LINEAR_O_FILES) $(LINEAR_M_FILES) $(LINEAR_C_FILES)

bench_linear_solver.dat: $(LINEAR_FILES)
	rm -f $@
	for n in $^ ; do \
             /bin/echo -n $$n | sed 's/[a-z_]*//g' | sed 's/[-.]/ /g' >> $@ ; \
             cat $$n >> $@ ; \
        done

.PHONY: bench_linear_solver_all
bench_linear_solver_all: $(LINEAR_FILES)

bench_linear_solver.pdf: plot_linear_solver.gp bench_linear_solver.dat matlab_linear_bench.dat
	gnuplot -e "load \"$<\""

.PHONY: bench_linear_solver_plot
bench_linear_solver_plot: bench_linear_solver.pdf