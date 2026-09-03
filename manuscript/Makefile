# emacs command
fname := psr-vs-gamms
df := ../results/results-power.rds ../results/choose-basis-fn_targets-25.rds
rf := compute-results.R
ecmd := emacs --batch -q -l dotemacs

/output/$(fname).pdf : $(fname).org $(df) $(rf) refs.bib setup.org abstract.tex
	$(ecmd) $< -f org-latex-export-to-pdf && cp $(fname).pdf /output/

abstract.tex : make-abstract.R $(rf) $(df)
	Rscript $<

clean :
	rm -fv $(fname).pdf abstract.tex $(fname).tex $(fnam3).bbl fig/*
	rm -rfv _minted-*
