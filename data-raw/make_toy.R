# Builds the simulated dataset cepip_toy (data/cepip_toy.rda).
# 1,200 cells from 12 donors. x is a gene of interest on a log1p-like scale.
# Eight partner genes follow x with a slope that increases above x = 1.8
# (an elbow), a donor offset and noise that grows with x. The composite y is
# their row mean. The values are simulated and not derived from any dataset.
set.seed(1)
n_donor <- 12L
n_per <- 100L
donor <- rep(sprintf("D%02d", seq_len(n_donor)), each = n_per)
off <- rnorm(n_donor, 0, 0.06)[as.integer(factor(donor))]
x <- pmax(0.05, rgamma(n_donor * n_per, shape = 5, rate = 2.2))
mu <- 0.6 + 0.12 * x + 0.5 * pmax(0, x - 1.8) + off
noise_sd <- 0.3 * (0.6 + 0.2 * x)
shared <- rnorm(length(x), 0, noise_sd)
partners <- sapply(1:8, function(j) pmax(0, mu + 0.9 * shared + rnorm(length(x), 0, 0.15)))
colnames(partners) <- paste0("P", 1:8)
cepip_toy <- data.frame(id = sprintf("cell%04d", seq_along(x)), donor = donor,
                        region = rep(c("R1", "R2"), length.out = length(x)),
                        x = round(x, 5), round(partners, 5), stringsAsFactors = FALSE)
cepip_toy$y <- rowMeans(as.matrix(cepip_toy[, paste0("P", 1:8)]))
save(cepip_toy, file = file.path("data", "cepip_toy.rda"), compress = "xz")
