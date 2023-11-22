# devtools::document()
devtools::load_all()


simdata3 <- sim.LogNormal(n=100, p=3, r=1, snr=10, d=3, seed=1, verbose=TRUE)
ternary(simdata3$X2, vhat=simdata3$Vlist, mu=simdata3$mu)

simdata4 <- sim.LogNormal(n=100, p=4, r=1, snr=10, d=3, seed=1, verbose=TRUE)
with(simdata4, quaternary3d(X2, vhat=Vlist, mu=mu))




fit3 <- CPCA(simdata3$X2, nrank=1)
fit4 <- CPCA(simdata4$X2, nrank=2)

reconstruction_error(fit3, data = simdata3)
reconstruction_error(fit4, data = simdata4)

ternary(simdata3$X2, vhat=fit3$vhat, xhat=fit3$xhat)
quaternary3d(simdata4$X2, vhat=fit4$vhat, xhat=fit4$xhat)





set.seed(1)
simdata.Linear <- sim.Linear(n=100, p=4, r=2, snr=5, d=3, seed=1, verbose=TRUE)
set.seed(1)
simdata.LogNormal <- sim.LogNormal(n=100, p=4, r=2, snr=5, d=3, seed=1, verbose=TRUE)
with(simdata.Linear, quaternary3d(X2, vhat=V, mu=mu))
with(simdata.LogNormal, quaternary3d(X2, vhat=Vlist, mu=mu))

fit.Linear <- fit_all(simdata.Linear$X2, nrank=2)
fit.LogNormal <- fit_all(simdata.LogNormal$X2, nrank=2)

reconstruction_error(fit.Linear, data = simdata.Linear, n.test = 1000)
reconstruction_error(fit.LogNormal, data = simdata.LogNormal, n.test = 1000)



projection(simdata.Linear$X2, fit.Linear$`lrPCA (0)`)$xhat %>% tail
fit.Linear$`lrPCA (0)`$xhat %>% tail

projection(simdata.Linear$X2, fit.Linear$crPCA)$xhat %>% tail
fit.Linear$crPCA$xhat %>% tail

projection(simdata.Linear$X2, fit.Linear$CPCA)$xhat %>% tail
fit.Linear$CPCA$xhat %>% tail

