# Regenerate all paper-02 figures.
#
#   make           # build all PDFs + PNGs
#   make data      # only re-run the Rust extractor
#   make clean     # remove out/
#
# Each fig script writes its PDF and PNG to out/.

REPO_ROOT := $(shell git rev-parse --show-toplevel)
RSCRIPT   := Rscript

DATA_FILES = data/fit_2017.csv data/fit_1998.csv data/rmse_summary.csv data/bed_profile.csv data/n_split_trajectory.csv

R_SCRIPTS = $(wildcard R/fig*.R)
PDFS      = $(patsubst R/fig%.R,out/fig%.pdf,$(R_SCRIPTS))

.PHONY: all data clean

all: data $(PDFS)

data: $(DATA_FILES)

data/fit_2017.csv data/fit_1998.csv data/rmse_summary.csv:
	cd $(REPO_ROOT) && cargo run --release -p hydroflux-autograd --example gen_figure_data

data/bed_profile.csv:
	cp $(REPO_ROOT)/examples/huasco_channel/output/huasco_longitudinal_profile.csv $@

# Split-Manning trajectory (fig06). ~18 min: 40 iters × 2 forward
# passes over the 21-day event. The example writes the CSV itself.
data/n_split_trajectory.csv:
	cd $(REPO_ROOT) && cargo run --release -p hydroflux-autograd --example calibrate_manning_huasco_2017_n_split

out/fig%.pdf: R/fig%.R R/theme_paper.R
	@mkdir -p out
	$(RSCRIPT) $<

clean:
	rm -rf out/
