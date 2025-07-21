devtools::install_github("heiniglab/scPower")
library(scPower)
power.general.withDoublets(
  nSamples = 3,
  nCells = 6000,
  readDepth = 35000,
  ct.freq = 0.5,
  type="de",
  ref.study=scPower::de.ref.study,
  ref.study.name,
  samplesPerLane,
  read.umi.fit,
  gamma.mixed.fits=scPower::gamma.mixed.fits,
  ct="CD8 T",
  disp.fun.param=scPower::disp.fun.param,
)
