
FUTHARK=futhark
GPU_DEVICE?=AMD
RUN_EQUI?=run_equilibrium.fut
BENCH_RUNS?=5
VAL_NAME?=default
EQB_DATA?=EQB_DATA_cars
EQB_VARIANT:=$(EQB_DATA:EQB_DATA_%=%)
SAVED_BENCH_DIR=saved_futhark_benchmarks
SAVED_BENCH_PREFIX=$(basename $(RUN_EQUI))-$(VAL_NAME)-bench$(BENCH_RUNS)

.PHONY: all
all: test
	$(MAKE) -C util dpsolve_demo.pdf

run_bellman_ad:
	$(FUTHARK) c run_bellman_ad.fut

run_bellman_j:
	$(FUTHARK) c run_bellman_j.fut

run_bellman:
	$(FUTHARK) c run_bellman.fut

.PHONY: clean
clean:
	rm -rf *~ *.c run_bellman_ad run_bellman_j run_bellman run_bellman.fut.* run_bellman_ad.fut.* run_bellman_j.fut.* run_equilibrium.fut.* *.exe *.tm *.tmavg *.tmstats plot1.pdf plot1j.pdf plot1ad.pdf plot2.pdf plot_eqb_cars.pdf plot_eqb_households.pdf plot_eqb_age.pdf plots.pdf data1.dat data2.dat *.out data1j.dat data1ad.dat bench_solve_cars.dat bench_solve_households.dat bench_solve_age.dat bench_solve_res_cars.dat bench_solve_res_households.dat bench_solve_res_age.dat bench_solve_sp_cars.dat bench_solve_sp_households.dat bench_solve_sp_age.dat bench_solve_sp_res_cars.dat bench_solve_sp_res_households.dat bench_solve_sp_res_age.dat
	$(MAKE) -C util clean
	$(MAKE) -C autotrade clean

.PHONY: test
test:
	$(FUTHARK) test -c run_bellman.fut run_bellman_ad.fut run_bellman_j.fut

.PHONY: futversion
futversion:
	$(FUTHARK) --version

%-c.exe: %.fut
	$(FUTHARK) c -o $@ $<

%-m.exe: %.fut
	$(FUTHARK) multicore -o $@ $<

%-o.exe: %.fut
	$(FUTHARK) opencl -o $@ $<

%-u.exe: %.fut
	$(FUTHARK) cuda -o $@ $<

run_bellman-%-C.tm: run_bellman-c.exe
	echo "$*" | sed 's/-/ /g' | ./$< -e mainn -r 5 -t $@ > run_bellman-$*-c.out

run_bellman-%-M.tm: run_bellman-m.exe
	echo "$*" | sed 's/-/ /g' | ./$< -e mainn -r 5 -t $@ > run_bellman-$*-m.out

run_bellman-%-O.tm: run_bellman-o.exe
	echo "$*" | sed 's/-/ /g' | ./$< -d $(GPU_DEVICE) -e mainn -r 5 -t $@ > run_bellman-$*-o.out

run_bellman_j-%-C.tm: run_bellman_j-c.exe
	echo "$*" | sed 's/-/ /g' | ./$< -e mainn -r 5 -t $@ > run_bellman_j-$*-c.out

run_bellman_j-%-M.tm: run_bellman_j-m.exe
	echo "$*" | sed 's/-/ /g' | ./$< -e mainn -r 5 -t $@ > run_bellman_j-$*-m.out

run_bellman_j-%-O.tm: run_bellman_j-o.exe
	echo "$*" | sed 's/-/ /g' | ./$< -d $(GPU_DEVICE) -e mainn -r 5 -t $@ > run_bellman_j-$*-o.out

run_bellman_ad-%-C.tm: run_bellman_ad-c.exe
	echo "$*" | sed 's/-/ /g' | ./$< -e mainn -r 5 -t $@ > run_bellman_ad-$*-c.out

run_bellman_ad-%-M.tm: run_bellman_ad-m.exe
	echo "$*" | sed 's/-/ /g' | ./$< -e mainn -r 5 -t $@ > run_bellman_ad-$*-m.out

run_bellman_ad-%-O.tm: run_bellman_ad-o.exe
	echo "$*" | sed 's/-/ /g' | ./$< -d $(GPU_DEVICE) -e mainn -r 5 -t $@ > run_bellman_ad-$*-o.out

%.tmavg: %.tm
	cat $< | python3 -c"import sys; nums = [int(n) for n in sys.stdin]; print(sum(nums)/len(nums)/1000000)" > $@

# Compute mean, stdev, and stdev of the mean (standard error)
%.tmstats: %.tm
	cat $< | python3 -c"import sys, statistics; nums=[int(n)/1e6 for n in sys.stdin]; m=statistics.mean(nums); sd=statistics.stdev(nums) if len(nums)>1 else 0.0; se=sd/len(nums)**0.5; print(f'{m} {sd} {se}')" > $@

%.png: %.pdf
	convert -density 300 $< $@

plots.pdf: plot1.pdf plot1j.pdf plot1ad.pdf plot2.pdf
	pdfunite $^ $@

# N-C-Ax
DATA1=2-10-20 4-10-20 6-10-20 8-10-20 10-10-20 \
     4-5-20 4-15-20 4-20-20 4-25-20 4-30-20 4-35-20

DATA1_C_FILES=$(DATA1:%=run_bellman-%-C.tmstats)
DATA1_M_FILES=$(DATA1:%=run_bellman-%-M.tmstats)
DATA1_O_FILES=$(DATA1:%=run_bellman-%-O.tmstats)
#DATA1_FILES=$(DATA1_O_FILES) $(DATA1_M_FILES) $(DATA1_C_FILES)
DATA1_FILES=$(DATA1_M_FILES) $(DATA1_C_FILES)

data1.dat: $(DATA1_FILES)
	for n in $^ ; do \
             /bin/echo -n $$n | sed 's/[a-z_]*//g' | sed 's/[-.]/ /g' >> $@ ; \
             cat $$n >> $@ ; \
        done

plot1.pdf: plot1-noo.gp data1.dat matlab.dat
	gnuplot -e "load \"$<\""

DATA1_J=4-5-25 4-15-25 4-20-25 4-25-25 4-30-25 4-35-25 4-40-25 4-45-25 4-50-25 4-55-25 4-60-25

DATA1_J_C_FILES=$(DATA1_J:%=run_bellman_j-%-C.tmstats)
DATA1_J_M_FILES=$(DATA1_J:%=run_bellman_j-%-M.tmstats)
DATA1_J_O_FILES=$(DATA1_J:%=run_bellman_j-%-O.tmstats)
#DATA1_J_FILES=$(DATA1_J_O_FILES) $(DATA1_J_M_FILES) $(DATA1_J_C_FILES)
DATA1_J_FILES=$(DATA1_J_M_FILES) $(DATA1_J_C_FILES)

data1j.dat: $(DATA1_J_FILES)
	for n in $^ ; do \
             /bin/echo -n $$n | sed 's/[a-z_]*//g' | sed 's/[-.]/ /g' >> $@ ; \
             cat $$n >> $@ ; \
        done

plot1j.pdf: plot1j-noo.gp data1j.dat matlabj.dat
	gnuplot -e "load \"$<\""

DATA1_AD=4-5-25 4-15-25 4-20-25

DATA1_AD_C_FILES=$(DATA1_AD:%=run_bellman_ad-%-C.tmstats)
DATA1_AD_M_FILES=$(DATA1_AD:%=run_bellman_ad-%-M.tmstats)
DATA1_AD_O_FILES=$(DATA1_AD:%=run_bellman_ad-%-O.tmstats)
#DATA1_AD_FILES=$(DATA1_AD_O_FILES) $(DATA1_AD_M_FILES) $(DATA1_AD_C_FILES)
DATA1_AD_FILES=$(DATA1_AD_M_FILES) $(DATA1_AD_C_FILES)

data1ad.dat: $(DATA1_AD_FILES)
	for n in $^ ; do \
             /bin/echo -n $$n | sed 's/[a-z_]*//g' | sed 's/[-.]/ /g' >> $@ ; \
             cat $$n >> $@ ; \
        done

plot1ad.pdf: plot1ad-noo.gp data1ad.dat matlabj.dat
	gnuplot -e "load \"$<\""

DATA2=2-35-25 4-35-25 6-35-25 8-35-25 10-35-25
DATA2_C_FILES=$(DATA2:%=run_bellman-%-C.tmstats)
DATA2_M_FILES=$(DATA2:%=run_bellman-%-M.tmstats)
DATA2_O_FILES=$(DATA2:%=run_bellman-%-O.tmstats)
#DATA2_FILES=$(DATA2_O_FILES) $(DATA2_M_FILES) $(DATA2_C_FILES)
DATA2_FILES=$(DATA2_M_FILES) $(DATA2_C_FILES)
data2.dat: $(DATA2_FILES)
	for n in $^ ; do \
             /bin/echo -n $$n | sed 's/[a-z_]*//g' | sed 's/[-.]/ /g' >> $@ ; \
             cat $$n >> $@ ; \
        done

plot2.pdf: plot2-noo.gp data2.dat matlab2.dat
	gnuplot -e "load \"$<\""


# =========================================================
# Equilibrium benchmark (matches MATLAB bench_run_illustrations_test)
# Format: n-c-Ax-acc0-transcost
# =========================================================

run_equilibrium-c.exe: $(RUN_EQUI)
	$(FUTHARK) c -o $@ $<

run_equilibrium-m.exe: $(RUN_EQUI)
	$(FUTHARK) multicore -o $@ $<

run_equilibrium-o.exe: $(RUN_EQUI)
	$(FUTHARK) opencl -o $@ $<

run_equilibrium-u.exe: $(RUN_EQUI)
	$(FUTHARK) cuda -o $@ $<

bench_solve-%-C.tm: run_equilibrium-c.exe
	echo "$*" | sed 's/-/ /g' | ./$< -e bench_solve -r $(BENCH_RUNS) -t $@ > run_equilibrium-$*-c.out

bench_solve-%-M.tm: run_equilibrium-m.exe
	echo "$*" | sed 's/-/ /g' | ./$< -e bench_solve -r $(BENCH_RUNS) -t $@ > run_equilibrium-$*-m.out

bench_solve-%-O.tm: run_equilibrium-o.exe
	echo "$*" | sed 's/-/ /g' | ./$< -d $(GPU_DEVICE) -e bench_solve -r $(BENCH_RUNS) -t $@ > run_equilibrium-$*-o.out

bench_solve-%-U.tm: run_equilibrium-u.exe
	echo "$*" | sed 's/-/ /g' | ./$< -d $(GPU_DEVICE) -e bench_solve -r $(BENCH_RUNS) -t $@ > run_equilibrium-$*-u.out

EQB_DATA_cars = 2-1-25-5-0 2-3-25-5-0 2-5-25-5-0 2-7-25-5-0 2-9-25-5-0 2-11-25-5-0 2-13-25-5-0
EQB_DATA_households = 2-7-25-5-0 3-7-25-5-0 4-7-25-5-0 5-7-25-5-0 6-7-25-5-0
EQB_DATA_age = 2-7-10-5-0 2-7-15-5-0 2-7-20-5-0 2-7-25-5-0 2-7-30-5-0 2-7-35-5-0

EQB_C_FILES=$($(EQB_DATA):%=bench_solve-%-C.tmstats)
EQB_M_FILES=$($(EQB_DATA):%=bench_solve-%-M.tmstats)
EQB_O_FILES=$($(EQB_DATA):%=bench_solve-%-O.tmstats)
EQB_U_FILES=$($(EQB_DATA):%=bench_solve-%-U.tmstats)

EQB_FILES=$(EQB_U_FILES) $(EQB_M_FILES) $(EQB_C_FILES)
#EQB_FILES=$(EQB_O_FILES) $(EQB_M_FILES) $(EQB_C_FILES)
#EQB_FILES=$(EQB_M_FILES) $(EQB_C_FILES)

EQB_C_OUT=$($(EQB_DATA):%=run_equilibrium-%-c.out)
EQB_M_OUT=$($(EQB_DATA):%=run_equilibrium-%-m.out)
EQB_O_OUT=$($(EQB_DATA):%=run_equilibrium-%-o.out)
EQB_U_OUT=$($(EQB_DATA):%=run_equilibrium-%-u.out)

EQB_OUT=$(EQB_U_OUT) $(EQB_M_OUT) $(EQB_C_OUT)
#EQB_OUT=$(EQB_O_OUT) $(EQB_M_OUT) $(EQB_C_OUT)
#EQB_OUT=$(EQB_M_OUT) $(EQB_C_OUT)

bench_solve_$(EQB_VARIANT).dat: $(EQB_FILES)
	rm -f $@
	for n in $^ ; do \
             /bin/echo -n $$n | sed 's/[a-z_]*//g' | sed 's/[-.]/ /g' >> $@ ; \
             cat $$n >> $@ ; \
        done

bench_solve_res_$(EQB_VARIANT).dat: $(EQB_OUT)
	rm -f $@
	for n in $^ ; do \
             /bin/echo -n $$n | sed 's/[a-z_]*//g' | sed 's/[-.]/ /g' >> $@ ; \
             cat $$n >> $@ ; \
        done

$(SAVED_BENCH_DIR)/$(SAVED_BENCH_PREFIX)-$(EQB_VARIANT).dat: bench_solve_$(EQB_VARIANT).dat
	mkdir -p $(SAVED_BENCH_DIR)
	cp $< $@

.PHONY: bench_solve_all
bench_solve_all: bench_solve_$(EQB_VARIANT).dat \
                 bench_solve_res_$(EQB_VARIANT).dat \
                 $(SAVED_BENCH_DIR)/$(SAVED_BENCH_PREFIX)-$(EQB_VARIANT).dat

plot_eqb_$(EQB_VARIANT).pdf: plot_eqb_$(EQB_VARIANT).gp bench_solve_$(EQB_VARIANT).dat matlab_eqb_$(EQB_VARIANT).dat $(SAVED_BENCH_DIR)/$(SAVED_BENCH_PREFIX)-$(EQB_VARIANT).dat
	gnuplot -e "load \"$<\""

.PHONY: plot_eqb_all
plot_eqb_all:
	$(MAKE) EQB_DATA=EQB_DATA_cars       VAL_NAME=$(VAL_NAME) RUN_EQUI=$(RUN_EQUI) BENCH_RUNS=$(BENCH_RUNS) plot_eqb_cars.pdf
	$(MAKE) EQB_DATA=EQB_DATA_households VAL_NAME=$(VAL_NAME) RUN_EQUI=$(RUN_EQUI) BENCH_RUNS=$(BENCH_RUNS) plot_eqb_households.pdf
	$(MAKE) EQB_DATA=EQB_DATA_age        VAL_NAME=$(VAL_NAME) RUN_EQUI=$(RUN_EQUI) BENCH_RUNS=$(BENCH_RUNS) plot_eqb_age.pdf

# =========================================================
# Equilibrium validation
# Format: n-c-Ax-acc0-transcost
# =========================================================

# Filename prefix encodes which RUN_EQUI was used and the VAL_NAME tag
# (override VAL_NAME to differentiate runs across machines/configurations).
VAL_PREFIX=$(basename $(RUN_EQUI))-$(VAL_NAME)-validate_solve

validation_files/$(VAL_PREFIX)-%-C.val: $(basename $(RUN_EQUI))-c.exe
	echo "$*" | sed 's/-/ /g' | ./$< -e validate_solve -r 1 > $@

validation_files/$(VAL_PREFIX)-%-M.val: $(basename $(RUN_EQUI))-m.exe
	echo "$*" | sed 's/-/ /g' | ./$< -e validate_solve -r 1 > $@

validation_files/$(VAL_PREFIX)-%-O.val: $(basename $(RUN_EQUI))-o.exe
	echo "$*" | sed 's/-/ /g' | ./$< -d $(GPU_DEVICE) -e validate_solve -r 1 > $@

validation_files/$(VAL_PREFIX)-%-U.val: $(basename $(RUN_EQUI))-u.exe
	echo "$*" | sed 's/-/ /g' | ./$< -d $(GPU_DEVICE) -e validate_solve -r 1 > $@

EQB_VAL_C_FILES=$($(EQB_DATA):%=validation_files/$(VAL_PREFIX)-%-C.val)
EQB_VAL_M_FILES=$($(EQB_DATA):%=validation_files/$(VAL_PREFIX)-%-M.val)
EQB_VAL_O_FILES=$($(EQB_DATA):%=validation_files/$(VAL_PREFIX)-%-O.val)
EQB_VAL_U_FILES=$($(EQB_DATA):%=validation_files/$(VAL_PREFIX)-%-U.val)

# Default: c + m. Uncomment to also validate opencl / cuda backends.
#EQB_VAL_FILES=$(EQB_VAL_U_FILES) $(EQB_VAL_M_FILES) $(EQB_VAL_C_FILES)
#EQB_VAL_FILES=$(EQB_VAL_O_FILES) $(EQB_VAL_M_FILES) $(EQB_VAL_C_FILES)
EQB_VAL_FILES=$(EQB_VAL_M_FILES) $(EQB_VAL_C_FILES)

.PHONY: validate_eqb
validate_eqb: $(EQB_VAL_FILES)

.PHONY: validate_eqb_all
validate_eqb_all:
	$(MAKE) EQB_DATA=EQB_DATA_cars       VAL_NAME=$(VAL_NAME) RUN_EQUI=$(RUN_EQUI) validate_eqb
	$(MAKE) EQB_DATA=EQB_DATA_households VAL_NAME=$(VAL_NAME) RUN_EQUI=$(RUN_EQUI) validate_eqb
	$(MAKE) EQB_DATA=EQB_DATA_age        VAL_NAME=$(VAL_NAME) RUN_EQUI=$(RUN_EQUI) validate_eqb

# =========================================================
# Equilibrium benchmark with simple_prices initial guess
# Format: n-c-Ax-acc0-transcost
# =========================================================

bench_solve_sp-%-C.tm: run_equilibrium-c.exe
	echo "$*" | sed 's/-/ /g' | ./$< -e bench_solve_simple_prices -r $(BENCH_RUNS) -t $@ > run_equilibrium-sp-$*-c.out

bench_solve_sp-%-M.tm: run_equilibrium-m.exe
	echo "$*" | sed 's/-/ /g' | ./$< -e bench_solve_simple_prices -r $(BENCH_RUNS) -t $@ > run_equilibrium-sp-$*-m.out

EQB_SP_C_FILES=$($(EQB_DATA):%=bench_solve_sp-%-C.tmstats)
EQB_SP_M_FILES=$($(EQB_DATA):%=bench_solve_sp-%-M.tmstats)

EQB_SP_FILES=$(EQB_SP_M_FILES) $(EQB_SP_C_FILES)

EQB_SP_C_OUT=$($(EQB_DATA):%=run_equilibrium-sp-%-c.out)

bench_solve_sp_$(EQB_VARIANT).dat: $(EQB_SP_FILES)
	rm -f $@
	for n in $^ ; do \
             /bin/echo -n $$n | sed 's/[a-z_]*//g' | sed 's/[-.]/ /g' >> $@ ; \
             cat $$n >> $@ ; \
        done

bench_solve_sp_res_$(EQB_VARIANT).dat: $(EQB_SP_C_OUT)
	rm -f $@
	for n in $^ ; do \
             /bin/echo -n $$n | sed 's/[a-z_]*//g' | sed 's/[-.]/ /g' >> $@ ; \
             cat $$n >> $@ ; \
        done

.PHONY: bench_solve_sp_all
bench_solve_sp_all: bench_solve_sp_$(EQB_VARIANT).dat bench_solve_sp_res_$(EQB_VARIANT).dat

.PHONY: ad-back-c
ad-back-c:
	$(FUTHARK) c run_bellman_ad.fut
	echo 10 25 200 | ./run_bellman_ad -r 5 -t /dev/stderr

.PHONY: back-opencl
back-opencl:
	$(FUTHARK) opencl run_bellman.fut
	echo 10 25 200 | ./run_bellman -d $(GPU_DEVICE) -r 5 -t /dev/stderr

.PHONY: ad-back-opencl
ad-back-opencl:
	$(FUTHARK) opencl run_bellman_ad.fut
	echo 10 25 200 | ./run_bellman_ad -d $(GPU_DEVICE) -r 5 -t /dev/stderr
