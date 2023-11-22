# devtools::document()
devtools::load_all()



simdata4 <- sim.LogNormal(n=100, p=4, r=1, snr=10, d=3, seed=1, verbose=TRUE)
with(simdata4, quaternary3d(X2, vhat=Vlist, mu=mu))


simdata3 <- sim.LogNormal(n=100, p=3, r=1, snr=10, d=3, seed=1, verbose=TRUE)
ternary(simdata3$X2, vhat=simdata3$Vlist, mu=simdata3$mu)



fit <- CPCA(simdata3$X2, nrank=1 )
reconstruction_error(fit, data = simdata3)


projection(fit, nrank=1)


quaternary3d(X, vhat=fit$vhat, xhat=fit$xhat)

