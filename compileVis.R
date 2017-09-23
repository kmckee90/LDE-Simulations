
source("lib/plotSim.R")

load("output/results_basic.RData")

pdf("results/Results_basic.pdf", height=6, width=9)
plotSimRandEffect(sim.1.embedD, plotParams, ylim=c(-0.35, 0.35))
compareSE(sim.1.embedD, plotParams, ylim=c(0, 0.2), xlim=c(-1, 0))

plotSim(sim.1.pEvent, plotParams, ylim=c(-0.35, 0.35))
compareSE(sim.1.pEvent, plotParams, ylim=c(0, 0.2), xlim=c(-1, 0))

plotSim(sim.1.SNR, plotParams, ylim=c(-0.35, 0.35))
compareSE(sim.1.SNR, plotParams, ylim=c(0, 0.2), xlim=c(-1, 0))

plotSim(sim.1.N, plotParams, ylim=c(-0.35, 0.35))
compareSE(sim.1.N, plotParams, ylim=c(0, 0.2), xlim=c(-1, 0))
dev.off()

