## hydroflux paper — EMS submission build
##
## Targets:
##   make             — build paper.pdf (3 xelatex passes + bibtex)
##   make body        — regenerate body_clean.tex from ../manuscript.md
##                      via pandoc + the cleanup pipeline
##   make clean       — remove build artifacts (.aux, .log, .bbl, .blg, .pdf,
##                      pandoc intermediates)
##   make distclean   — clean + remove body_clean.tex

LATEX  := xelatex -interaction=nonstopmode -halt-on-error
BIBTEX := bibtex

PAPER       := paper
MD_SOURCE   := ../manuscript.md
PANDOC_INTERMEDIATE := manuscript_pandoc.tex
BODY_RAW    := body_raw.tex
BODY_FROM_INTRO := body_from_intro.tex
BODY_CLEAN  := body_clean.tex

.PHONY: all body clean distclean

all: $(PAPER).pdf

$(PAPER).pdf: $(PAPER).tex $(BODY_CLEAN) references.bib elsarticle.cls elsarticle-harv.bst
	$(LATEX) $(PAPER).tex
	$(BIBTEX) $(PAPER)
	$(LATEX) $(PAPER).tex
	$(LATEX) $(PAPER).tex

body: $(BODY_CLEAN)

$(BODY_CLEAN): $(MD_SOURCE)
	pandoc $(MD_SOURCE) --to=latex --natbib --standalone --output=$(PANDOC_INTERMEDIATE)
	awk '/^\\hypertarget{introduction}/,/^\\end{document}/' $(PANDOC_INTERMEDIATE) > $(BODY_FROM_INTRO)
	sed -E '/^\\hypertarget\{/d; s/^\\section\{[0-9]+\. /\\section{/g; s/^\\subsection\{[0-9]+\.[0-9]+ ?/\\subsection{/g; s/\}\\label\{[^}]+\}\}$$/}/; s/\\textasciitilde\\?[ ]?/$$\\sim$$/g' $(BODY_FROM_INTRO) > $(BODY_CLEAN).tmp
	awk '/^\\section\{Figures/{exit} {print}' $(BODY_CLEAN).tmp > $(BODY_CLEAN)
	python3 fix_superscripts.py $(BODY_CLEAN)
	rm -f $(BODY_CLEAN).tmp

clean:
	rm -f $(PAPER).aux $(PAPER).log $(PAPER).bbl $(PAPER).blg $(PAPER).pdf $(PAPER).out
	rm -f $(PANDOC_INTERMEDIATE) $(BODY_RAW) $(BODY_FROM_INTRO)

distclean: clean
	rm -f $(BODY_CLEAN)
