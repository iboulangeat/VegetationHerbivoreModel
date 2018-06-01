library(rmarkdown)
library(knitr)

# Analyses
#=============

knit("Calibration.Rmd", "Calibration.md")
knit("Sensitivity_analysis.Rmd", "Sensitivity_analysis.md")
knit("suppGradients.Rmd", "suppGradients.md")
knit("Fisichelli.Rmd", "Fisichelli.md")
knit("MainAnalysis.Rmd", "Readme.md")
knit("MixedFeeder.Rmd", "MixedFeeder.md")
knit("preference_pR.Rmd", "preference_pR.md")



# Graphiques
#=============
rm(list=ls())
source("plots_article.r")
##---

load("calibration.RData")

pdf(file = "../figs_article/supp_model_impacts.pdf", width = 6, height = 3.5)
plot_impacts()
dev.off()


pdf("../figs_article/supp_fitnesses.pdf", width = 9, height=4)
plot_fitnesses()
dev.off()

##----

load("Sensitivity_analysis.RData")

pdf(file = "../figs_article/supp_SA_1.pdf")
plot_sa.a()
dev.off()

pdf(file = "../figs_article/supp_SA_2.pdf")
plot_sa.b()
dev.off()

##---

load("suppGradients.RData")

pdf("../figs_article/supp_addGradients.pdf", width=10, height=15)
plot_suppGradients()
dev.off()

##---

load("Fisichelli.RData")
pdf(file="../figs_article/Fisichelli.pdf", height=6, width=12)
plot_fisichelli()
dev.off()

##---

load("MainAnalysis.RData")

pdf(file="../figs_article/equilibrium.pdf", height=12, width=5)
plot_equilibrium()
dev.off()

pdf(file="../figs_article/allMetrics.pdf", height=15, width=10)
plot_transient()
dev.off()

pdf(file="../figs_article/correlations.pdf", height=9, width=15)
plot_correlations()
dev.off()

pdf(file="../figs_article/trajectories.pdf", height=4, width=12)
plot_trajectories()
dev.off()

##---

load("MixedFeeder.RData")

pdf(file="../figs_article/equilibrium_mf.pdf", height=12, width=5)
plot_equilibrium()
dev.off()

pdf(file="../figs_article/allMetrics_mf.pdf", height=15, width=10)
plot_transient(ylim.H=5)
dev.off()

##---

load("preference_pR.RData")

pdf("../figs_article/allMetrics_pR.pdf", width=10, height=15)
plot_transient()
dev.off()

pdf("../figs_article/correlations_pR.pdf", width = 15, height=9)
plot_correlations_pR()
dev.off()
