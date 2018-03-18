## ------------------------------------------------------------------------
koala.trans.mat <- matrix(c(0.000,0.000,0.302,0.302,
                              0.940,0.000,0.000,0.000,
                              0.000,0.884,0.000,0.000,
                              0.000,0.000,0.793,0.793),
                            nrow = 4, ncol = 4, byrow = TRUE)
colnames(koala.trans.mat) <- rownames(koala.trans.mat) <- c('Stage_0_1','Stage_1_2','Stage_2_3','Stage_3_')

koala.trans.mat.es <- matrix(c(0.000,0.000,1,1,
                                 1,0.000,0.000,0.000,
                                 0.000,1,0.000,0.000,
                                 0.000,0.000,1,1),
                               nrow = 4, ncol = 4, byrow = TRUE)
colnames(koala.trans.mat.es) <- rownames(koala.trans.mat.es) <- c('Stage_0_1','Stage_1_2','Stage_2_3','Stage_3_')


## ---- fig.show='hold'----------------------------------------------------
plot(1:10)
plot(10:1)

