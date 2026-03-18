
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)

DiscoverAlleleInPoolSequencing <- function(PopulationSize, PopulationAlleleFreq,
                                           PoolSize, SequenceCoverage, SequenceGammaParameter,
                                           SequenceErrorRate, ThresholdForDiscovery) {
  # PopulationSize - integer, size of population
  # PopulationAlleleFreq - numeric, allele frequency of the alternative allele
  # PoolSize - integer, size of the pool of individuals that will be sequenced
  # SequenceCoverage - numeric, average sequencing coverage
  # SequenceGammaParameter - numeric, sequencability parameter (see Li et al., ????)
  # SequenceErrorRate - numeric, error rate that will "mangle" sequence reads and we might thereore miss alternative allele
  # ThresholdForDiscovery - numeric, how many reads with alternative alleles do we must see to be confident the alternative allele is not an error

  # Simulate population genotypes ----
  # ... assume Hardy-Weinberg equilibrium
  HWGenotypeFreq <- c((1 - PopulationAlleleFreq)^2,
                      2 * (1 - PopulationAlleleFreq) * PopulationAlleleFreq,
                      PopulationAlleleFreq^2)
  PopulationGenotypes <- sample(x=c(0L, 1L, 2L), size=PopulationSize,
                                prob=HWGenotypeFreq, replace=TRUE)

  # Make a random pool ----
  PoolGenotypes <- sample(x=PopulationGenotypes, size=PoolSize, replace=FALSE)

  # Sequence the pool ----
  # ... convert genotypes to alleles
  PoolAlleles <- matrix(data=0L, nrow=2L, ncol=PoolSize)
  for (Individual in 1L:PoolSize) {
    if        (PoolGenotypes[Individual] == 1L) {
      PoolAlleles[1L, Individual] <- 0L
      PoolAlleles[2L, Individual] <- 1L
    } else if (PoolGenotypes[Individual] == 2L) {
      PoolAlleles[1L, Individual] <- 0L
      PoolAlleles[2L, Individual] <- 1L
    }
  }
  PoolAlleles <- c(PoolAlleles)
  # ... assume binomal sampling
  PoolSequenceReadsBinomial <- sample(x=PoolAlleles, size=SequenceCoverage, replace=TRUE)
  # ... assume poisson sampling to account for variation in realized coverage
  RealizedSequenceCoverage <- rpois(n=1L, lambda=SequenceCoverage)
  PoolSequenceReadsPoisson <- sample(x=PoolAlleles, size=RealizedSequenceCoverage, replace=TRUE)
  # ... assume poisson-gamma sampling to account for variation in realized coverage and variation of sequencability along genome
  RealizedSequencability <- rgamma(n=1L, shape=SequenceGammaParameter, scale=1/SequenceGammaParameter)
  RealizedSequenceCoverage <- rpois(n=1L, lambda=SequenceCoverage * RealizedSequencability)
  PoolSequenceReadsPoissonGamma <- sample(x=PoolAlleles, size=RealizedSequenceCoverage, replace=TRUE)
  # ... apply error/mapping/... rate
  MangleRead <- function(x, ErrorRate) {
    n <- length(x)
    Error <- as.logical(rbinom(n=n, size=1, prob=ErrorRate))
    Test <- Error & x > 0
    x[Test] <- 0 # assume that error causes mapping etc to default to the reference allele
    x
  }
  PoolSequenceReadsBinomial     <- MangleRead(x=PoolSequenceReadsBinomial,     ErrorRate=SequenceErrorRate)
  PoolSequenceReadsPoisson      <- MangleRead(x=PoolSequenceReadsPoisson,      ErrorRate=SequenceErrorRate)
  PoolSequenceReadsPoissonGamma <- MangleRead(x=PoolSequenceReadsPoissonGamma, ErrorRate=SequenceErrorRate)

  # Did we captured the alternative allele? ---
  DiscoveryBinomial     <- sum(PoolSequenceReadsBinomial)     > ThresholdForDiscovery
  DiscoveryPoisson      <- sum(PoolSequenceReadsPoisson)      > ThresholdForDiscovery
  DiscoveryPoissonGamma <- sum(PoolSequenceReadsPoissonGamma) > ThresholdForDiscovery

  # Return ---
  c(DiscoveryBinomial=DiscoveryBinomial,
    DiscoveryPoisson=DiscoveryPoisson,
    DiscoveryPoissonGamma=DiscoveryPoissonGamma)
}


nReplicates <- 1000L*1000L
Tmp <- matrix(data=FALSE, nrow=nReplicates, ncol=3L)
colnames(Tmp) <- c("DiscoveryBinomial", "DiscoveryPoisson", "DiscoveryPoissonGamma")
for (Replicate in 1L:nReplicates) {
  Tmp[Replicate, ] <- DiscoverAlleleInPoolSequencing(PopulationSize=as.integer(args[3]),
                                                     PopulationAlleleFreq=as.numeric(args[1]),
                                                     PoolSize=as.integer(args[2]),
                                                     SequenceCoverage=as.numeric(args[4]),
                                                     SequenceGammaParameter=8,
                                                     SequenceErrorRate=0.01,
                                                     ThresholdForDiscovery=1L)
}
Result <- colMeans((Tmp))
print(Result)
