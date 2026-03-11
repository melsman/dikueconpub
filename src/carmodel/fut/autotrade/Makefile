FUTHARK=futhark

#FUTHARK=$(HOME)/.local/bin/futhark

.PHONY: all
all:
	$(FUTHARK) c trmodel.fut

.PHONY: clean
clean:
	rm -rf *~ doc *.c trmodel trmodel_test trmodel_static trmodel_test.fut.*

.PHONY: futversion
futversion:
	$(FUTHARK) --version

.PHONY: test
test:
	$(FUTHARK) test trmodel_test.fut

.PHONY: bench_test
bench_test:
	$(FUTHARK) test trmodel_benchmark.fut

.PHONY: static
static:
	$(FUTHARK) c trmodel_static.fut

doc:
	$(FUTHARK) doc -v -o doc .