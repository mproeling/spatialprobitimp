# Spatialprobitimp
Spatial Probit Imputation functions belonging to Roeling & Nicholls, Imputation of attributes in networked data using Bayesian autocorrelation regression models, Social Networks, 62: https://www.sciencedirect.com/science/article/abs/pii/S0378873320300149

this function is an adaption from sar_combined_mcmc from the spatialprobit package: https://github.com/cran/spatialprobit

##### CUT MODEL = sar_imp_cut
##### Full Bayes MODEL = sar_imp_fbs

- it takes the same objects but performs imputation as explained in Roeling & Nicholls (2020), Social Networks (62)
- the main difference is that most objects are lists so that multiple Ys with missing data can be used, but I guess this requires some more tweaking, with a proper use case this is easily implemented.
- the traces of rho are written to a file
- the first column of X is filled with 1 (for the intercept)
- the imputation output is in $results$y_imp
- starting values for rho come from Dittrich et al. (2017), Social Networks
- if y = 0, 1 make sure it is 0 1 not 1 2
- make sure the attribute data in Y, X, are sorted so that they agree with W (so that Wy calculation is correct)
- recode categorical independent variables to dummy if categories > 2
- method has to be defined and has no default eg method = c("probit")
- model = "edge_conditioned_25_flag0" is a  string value used to write away the Rho traces (to a file), but would be more logical to add this to a vectosrs called rhodraws (such as bdraws)

example :
impmodel.cluster.normal.10 = sar_combined_mcmc_imp_cut(y = yf.cluster, x=X, W, ndraw=10000, burn.in=100, thinning=1,start = start,
                                                       prior=list(a1=1, a2=1, c=rep(0, ncol(X)), T=diag(ncol(X))*1e12, lflag = 1),
                                                       m=10, computeMarginalEffects=TRUE, showProgress=TRUE, method = c("probit"),
                                                       model = "cluster10.normal")

write.table(impmodel.cluster.normal.10$y_imp[[1]], "y_imp_normalcut.csv", sep = ";", quote = F, col.names = F, row.names = F)

at the moment im busy finishing my phd so making this a nice package has low priority but if anybody has the ambition let me know.
