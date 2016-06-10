library(jaatha)
dm<-dm.createDemographicModel(sample.sizes=c(185,107), loci.num=215,seq.length=130)
dm<-dm.addSpeciationEvent(dm,.1,2,new.time.point.name='T')
dm<-dm.addSymmetricMigration(dm,.01,10,new.par.name='m')
dm<-dm.addMutation(dm,.01,1,new.par.name='theta')
# dm<-dm.addRecombination(dm,fixed=1)


jsfs=read.table('jsfs2D.txt',sep=',')

test1=Jaatha.initialize(dm, jsfs, 1, cores = 0, scaling.factor = 1, use.shm = FALSE, folded = TRUE, smoothing = FALSE)

Jaatha.initialSearch(test1, sim = 200, blocks.per.par = 2, rerun = FALSE)

library(jaatha)
source('script/simulator2.R')

sim_func <- function(x) rpois(10, x)

sum_stats <- list(create_jaatha_stat("id", function(x, opts) x))

par_ranges <- matrix(c(0.1, 0.1, 10, 10), 2, 2)
rownames(par_ranges) <- c("x", "y")
colnames(par_ranges) <- c("min", "max")

jaatha_model <- create_jaatha_model(sim_func, par_ranges, sum_stats)



jaatha_data <- create_jaatha_data(data_obs, jaatha_model)


estimates <- jaatha(jaatha_model, jaatha_data, sim = 100, repetitions = 2, verbose = TRUE,cores=3)








par.ranges <- matrix(c(0.1, 0.1, 0.1, 10, 10, 10), 3, 2)
rownames(par.ranges) <- c('x', 'y', 'z')
colnames(par.ranges) <- c('min', 'max')
create_jaatha_model( function(sim.pars) list(poisson.vector=sampleFromModel(sim.pars[1], sim.pars[2], sim.pars[3])),
par_ranges=par.ranges, sum_stats=sum.stats)
create_jaatha_model( function(sim.pars) sampleFromModel(sim.pars[1], sim.pars[2], sim.pars[3]),
par_ranges=par.ranges, sum_stats=sum.stats)

create_jaatha_model( function(x) c(rpois(10, x[1]), rpois(10, x[2]), rpois(10, x[3])) ,
par_ranges=par.ranges, sum_stats=sum.stats)

create_jaatha_model( function(x) rpois(10, x) ,par_ranges=par.ranges, sum_stats=list("id",identity))

# jaatha <- new('Jaatha', sim.func, par.ranges, sum.stats)
jaatha <- Jaatha.initialSearch(jaatha, 100, 2)
jaatha <- Jaatha.refinedSearch(jaatha, 1, 100)


create_jaatha_model(function(x) rpois(10, x),
par_ranges = matrix(c(0.1, 0.1, 10, 10), 2, 2),
sum_stats = list(create_jaatha_stat("sum", sum)))

