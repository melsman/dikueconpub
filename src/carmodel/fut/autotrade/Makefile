FUTHARK=futhark

#FUTHARK=$(HOME)/.local/bin/futhark

.PHONY: all
all:
	$(FUTHARK) c trmodel.fut

.PHONY: clean
clean:
	rm -rf *~ doc *.c trmodel trmodel_test trmodel_static trmodel_test trmodel_benchmark.fut.* \
	       trmodel_benchmark-c.exe trmodel_benchmark-m.exe trmodel_benchmark-o.exe \
	       bench_bellmanJ-*.tm bench_bellmanJ-*.tmavg trmodel_benchmark-*.out \
	       bench_bellmanJ.dat

.PHONY: futversion
futversion:
	$(FUTHARK) --version

.PHONY: test
test:
	$(FUTHARK) test trmodel_test.fut

.PHONY: bench_test
bench_test:
	$(FUTHARK) test trmodel_benchmark.fut

doc:
	$(FUTHARK) doc -v -o doc .

# ---- compile benchmark executable for each backend ----

trmodel_benchmark-c.exe: trmodel_benchmark.fut
	$(FUTHARK) c -o $@ $<

trmodel_benchmark-m.exe: trmodel_benchmark.fut
	$(FUTHARK) multicore -o $@ $<

trmodel_benchmark-o.exe: trmodel_benchmark.fut
	$(FUTHARK) opencl -o $@ $<

# ---- run benchmark for bench_bellmanJ ----
# stem format: n-c-Ax-N
# example: 2-2-2-10

bench_bellmanJ-%-C.tm: trmodel_benchmark-c.exe
	echo "$*" | sed 's/-/ /g' | ./$< -e bench_bellmanJ -r 5 -t $@ > trmodel_benchmark-$*-c.out

bench_bellmanJ-%-M.tm: trmodel_benchmark-m.exe
	echo "$*" | sed 's/-/ /g' | ./$< -e bench_bellmanJ -r 5 -t $@ > trmodel_benchmark-$*-m.out

bench_bellmanJ-%-O.tm: trmodel_benchmark-o.exe
	echo "$*" | sed 's/-/ /g' | ./$< -d AMD -e bench_bellmanJ -r 5 -t $@ > trmodel_benchmark-$*-o.out

# ---- average timing file ----

%.tmavg: %.tm
	cat $< | python3 -c "import sys; nums = [int(n) for n in sys.stdin]; print(sum(nums)/len(nums)/1000000)" > $@

# ---- benchmark collections ----

BELLMANJ_DATA=2-20-5-100 2-20-10-100 2-20-15-100 2-20-20-100 2-20-25-100 2-20-30-100 2-20-35-100 2-20-40-100     

BELLMANJ_C_FILES=$(BELLMANJ_DATA:%=bench_bellmanJ-%-C.tmavg)
BELLMANJ_M_FILES=$(BELLMANJ_DATA:%=bench_bellmanJ-%-M.tmavg)
BELLMANJ_O_FILES=$(BELLMANJ_DATA:%=bench_bellmanJ-%-O.tmavg)

# choose which backends to include
BELLMANJ_FILES=$(BELLMANJ_M_FILES) $(BELLMANJ_C_FILES)
#BELLMANJ_FILES=$(BELLMANJ_O_FILES) $(BELLMANJ_M_FILES) $(BELLMANJ_C_FILES)

bench_bellmanJ.dat: $(BELLMANJ_FILES)
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