library(sequenza)


test <- sequenza.extract(snakemake@input[["input"]], verbose = FALSE)


CP <- sequenza.fit(test)


sequenza.results(sequenza.extract = test,
    cp.table = CP, sample.id = snakemake@wildcards[["sample"]],
    out.dir=snakemake@params[["dir"]])
