#   Script to calculate the mutation rate from counts of homozygous
#   SNPs in Fast-Neutron irradiated lines of varying levels of selfing.
#   Derivation of this equation will be shown separately.
#   Written by Thomas JY Kono

#   Take arguments
args <- commandArgs(TRUE)

#   Define a function to estimate the expected number of homozygous
#   spontaneous after a certain number of generations.
#   This is a nested sum of sums and is not very clean.
#   The default mutation rate is from Ossowski et al., 2010, and the
#   default genome size is the diploid size of the queried portion of the
#   Glycine max genome, as defined by sequence coverage.
spontaneous_homozygous <- function(
    gens,
    mutrate=5.94982e-9,
    genomesize=1740930036) {
    summation <- 0
    for(i in 1:gens) {
        for(j in 1:i) {
            add <- mutrate * genomesize * (1/2)^(j-1)
            summation <- summation + add
        }
    }
    return(summation)
}

#   Define a function to implement the equation. The defaults for the
#   mutation rate are from Ossowski et al. 2010 in Science, and the
#   genome size is the diploid size of the Glycine max genome assembly (V1)
#   that pass our alignment-filtering criteria.
mutation_rate <- function(x, mutrate=5.94982e-9, genomesize=1740930036) {
    counts <- as.numeric(x[1])
    gens <- as.numeric(x[2])
    spont_hom <- 0.25 * spontaneous_homozygous(
        gens,
        mutrate=mutrate,
        genomesize=genomesize)
    induced_hom <- counts - spont_hom
    scaled_hom <- (1/2) - (1/2)^(gens+1)
    rate <- (induced_hom) / (genomesize * scaled_hom)
    return(rate)
}

#   Read in the data
fndata <- read.table(args[1], header=TRUE)

#   Set some constants. These are values for A. thaliana that were derived from
#   Ossowski et al. 2010 and the TAIR9 genome assembly (best guess at which
#   version they used). These numbers are haploid numbers.
a_size <- 38072592
t_size <- 38027796
c_size <- 21440282
g_size <- 21419561
wg_size <- a_size + t_size + c_size + g_size

at2gc_mutrate <- 0.09 / (a_size + t_size)   #   Transition      A->G, T->C
cg2ta_mutrate <- 0.41 / (c_size + g_size)   #   Transition      C->T, G->A
at2ta_mutrate <- 0.04 / (a_size + t_size)   #   Transversion    A->T, T->A
cg2at_mutrate <- 0.06 / (c_size + g_size)   #   Transversion    C->A, G->T
at2cg_mutrate <- 0.04 / (a_size + t_size)   #   Transversion    A->C, T->G
cg2gc_mutrate <- 0.05 / (c_size + g_size)   #   Transversion    C->G, G->C
transition_rate <- (0.09 + 0.41) / wg_size
transversion_rate <- (0.04 + 0.06 + 0.04 + 0.05) / wg_size

#   These are the counts of the nucleotides for the Gmax assembly V1.
#   We multiply by 2, since we are estimating diploid rates.
gm_a_size <- 285933037
gm_t_size <- 285940059
gm_c_size <- 149269536
gm_g_size <- 149321878

gm_at <- (gm_a_size + gm_t_size) * 2
gm_cg <- (gm_c_size + gm_g_size) * 2
gm_wg <- (gm_a_size + gm_t_size + gm_c_size + gm_g_size) * 2

#   Calculate whole genome mutation rate
wg_mutrates <- as.numeric(
    apply(
        cbind(fndata$Nmut, fndata$GensSelf),
        1,
        mutation_rate
        )
    )

#   AT->GC mutation rate
at_gc <- as.numeric(
    apply(
        cbind(fndata$A2G+fndata$T2C, fndata$GensSelf),
        1,
        mutation_rate,
        mutrate=at2gc_mutrate,
        genomesize=gm_at
        )
    )
#   GC -> AT mutation rate
gc_at <- as.numeric(
    apply(
        cbind(fndata$G2A+fndata$C2T, fndata$GensSelf),
        1,
        mutation_rate,
        mutrate=cg2ta_mutrate,
        genomesize=gm_cg
        )
    )
#   What is the final transition rate?
ti_mutrates <- as.numeric(
    apply(
        cbind(fndata$Ti, fndata$GensSelf),
        1,
        mutation_rate,
        mutrate=transition_rate,
        genomesize=gm_wg
        )
    )

#   Next, we do transversions
#   AT -> CG
at_cg <- as.numeric(
    apply(
        cbind(fndata$A2C+fndata$T2G, fndata$GensSelf),
        1,
        mutation_rate,
        mutrate=at2cg_mutrate,
        genomesize=gm_at
        )
    )
#   AT -> TA
at_ta <- as.numeric(
    apply(
        cbind(fndata$A2T+fndata$T2A, fndata$GensSelf),
        1,
        mutation_rate,
        mutrate=at2ta_mutrate,
        genomesize=gm_at
        )
    )
#   GC -> TA
gc_ta <- as.numeric(
    apply(
        cbind(fndata$G2T+fndata$C2A, fndata$GensSelf),
        1,
        mutation_rate,
        mutrate=cg2at_mutrate,
        genomesize=gm_at
        )
    )
#   CG -> GC
cg_gc <- as.numeric(
    apply(
        cbind(fndata$C2G+fndata$G2C, fndata$GensSelf),
        1,
        mutation_rate,
        mutrate=cg2gc_mutrate,
        genomesize=gm_cg
        )
    )
#   And the general transversion ratop
tv_mutrates <- as.numeric(
    apply(
        cbind(fndata$Tv, fndata$GensSelf),
        1,
        mutation_rate,
        mutrate=transversion_rate,
        genomesize=gm_wg
        )
    )

#   Stick it all in a pretty data frame
final <- data.frame(
    Lines=fndata$Name,
    Selfing_Gens=fndata$GensSelf,
    Dosage=fndata$Dosage,
    Nmut=fndata$Nmut,
    Mutation_Rate=wg_mutrates,
    Ti=fndata$Ti,
    Ti_Rate=ti_mutrates,
    Tv=fndata$Tv,
    Tv_Rate=tv_mutrates,
    TiTv_Ratio=ti_mutrates/tv_mutrates,
    AT2GC_Rate=at_gc,
    GC2AT_Rate=gc_at,
    AT2CG_Rate=at_cg,
    AT2TA_Rate=at_ta,
    GC2TA_Rate=gc_ta,
    CG2GC_Rate=cg_gc
    )

#   And finally, write it to a table!
fname <- gsub(".txt", "_MutRates_Grouped.txt", args[1])
write.table(
    final,
    file=fname,
    quote=FALSE,
    sep="\t",
    col.names=TRUE,
    row.names=FALSE
    )
print(final)
