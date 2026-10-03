Overleaf package -- "One-variable equations of length three over symmetric groups"
Journal of Symbolic Computation (Elsevier), elsarticle class, preprint layout.

Files
  p3_article.tex        main document (set as main file in Overleaf)
  p3_refs.bib           bibliography (29 entries, all cited)
  elsarticle.cls        Elsevier class, TeX Live 2022 copy (Overleaf also ships it; the copy pins the version)
  elsarticle-num.bst    numbered bibliography style used by \bibliographystyle
  figures/fig_*.pdf     the five figures (vector); fig_*.png are 300 dpi copies for the submission system
                        (\graphicspath also accepts ../figures/, the layout of the order3-conjugacy repository)
  make_figures.py       regenerates every figure from the recorded JSON outputs in data/ (python make_figures.py data figures)
  data/*.json           recorded outputs from the order3-conjugacy and CPWI repositories that the figures read

Compiling
  Overleaf: compiler pdfLaTeX, main document p3_article.tex. Locally: pdflatex; bibtex p3_article; pdflatex; pdflatex.
  Checked with Tectonic (TeX Live 2022 bundle): 0 errors, 0 warnings, 20 pages; only the elsarticle running-header box
  reports a 2.6 pt overfull, which is the class's own and not the text.

Before submission
  * The Reproducibility paragraph points to https://github.com/gautierfilardo-efrei/order3-conjugacy; add the Zenodo DOI once minted.
  * Switch \documentclass[preprint,12pt] to [review,12pt] if the journal asks for line numbers, or to [final,5p] for the
    two-column final layout; nothing else changes.
  * Highlights (3-5 bullets, <= 85 characters each) and a cover letter are separate files in the Elsevier submission system.
