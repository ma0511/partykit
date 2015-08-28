Sweave("treatment_interactions.Rnw")
tools::texi2dvi("treatment_interactions.tex", pdf=TRUE, clean=T)
