# emacs command
fname := supplementary-materials
df := ../results/choose-basis-fn_targets-25.rds
ecmd := emacs --batch -q -l dotemacs

/output/$(fname).pdf : $(fname).org $(df) $(rf) refs.bib
	$(ecmd) $< -f org-latex-export-to-pdf && cp $(fname).pdf /output/

clean :
	rm -fv $(fname).pdf $(fname).tex $(fnam3).bbl fig/*
	rm -rfv _minted-*
