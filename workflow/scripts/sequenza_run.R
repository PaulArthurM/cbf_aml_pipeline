library(sequenza)


test <- sequenza.extract(snakemake@input[[1]], verbose = FALSE)


CP <- sequenza.fit(test)


sequenza.results(sequenza.extract = test,
    cp.table = CP, sample.id = snakemake@wildcards[[0]],
    out.dir=snakemake@output[[1]])
