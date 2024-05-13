library(tidyverse)
library(RcppAlgos)

# 3 observed trees
no <- 3
y <- c(0, 1, 1)
# 5 available trees
nt <- 5
x <- rnorm(nt)
b <- rnorm(1)
probs <- plogis(b * x)

# unique combinations (order doesn't matter within dead or alive) ---------

dead <- sum(y == 0)
alive <- sum(y)

combs_dead <- combn(1:nt, dead) %>% t
combs_alive <- combn(1:nt, alive) %>% t

cols_merge <- expand.grid(dd = 1:nrow(combs_dead),
                          aa = 1:nrow(combs_alive))

combs_full <- matrix(NA, nrow(cols_merge), length(y))
for(j in 1:nrow(cols_merge)) {
  # j = 1
  combs_full[j, ] <- c(combs_dead[cols_merge$dd[j], ],
                       combs_alive[cols_merge$aa[j], ])
}
unique_len <- apply(combs_full, 1, function(x) length(unique(x)))
combs_uni <- combs_full[unique_len == length(y), ]

n_uni <- nrow(combs_uni)
probs_uni <- numeric(n_uni)
for(i in 1:n_uni) probs_uni[i] <- exp(sum(dbinom(y, prob = probs[combs_uni[i, ]],
                                                 size = 1, log = T)))

# Repeating combinations --------------------------------------------------

combs_all <- permuteGeneral(1:nt, no)
n_all <- nrow(combs_all)
probs_all <- numeric(n_all)
for(i in 1:n_all) probs_all[i] <- exp(sum(dbinom(y, prob = probs[combs_all[i, ]],
                                                 size = 1, log = T)))

probs_all[c(50, 55)]

# compare
mean(probs_all)
mean(probs_uni)
