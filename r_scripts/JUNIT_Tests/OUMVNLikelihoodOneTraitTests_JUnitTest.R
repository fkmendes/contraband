# author: Fabio K. Mendes
# This R script gives us the expected values for JUnit tests
# (1) 'OUMVNLikelihoodOneTraitTest'
# (2) 'OUMVNLikelihoodOneTraitTest2'
# (3) 'OUMVNLikelihoodOneTraitTest3'

library(mvMORPH)
library(TreeSim)

## (1) 'OUMVNLikelihoodOneTraitTest'

## tr <- read.tree(text="(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):1.0,sp4:3.0);")
## plotTree(tr, node.numbers=T) # to see node numbers
tr <- paintSubTree(read.tree(text="(((sp1:1.0,sp2:1.0):1.0,sp3:2.0):1.0,sp4:3.0);"), node=6, state=2)
tr <- paintSubTree(tr, node=3, state=3, stem=TRUE)
plotSimmap(tr) # to see

## so we have a 4-taxon tree, with 3 optima
set.seed(123)
dat <- mvSIM(tr, model="OUM", param=list(sigma=0.1, theta=c(0, 0.2, 0.4, 0.6), alpha=1, root=TRUE))
row.names(dat) <- c("sp1","sp2","sp3","sp4")

### lnLk1
res <- mvOU(tr, dat, model="OUM", param=list(vcv="fixedRoot", root=TRUE))
res$sigma # 0.006082604
res$alpha # 7.390366
res$theta # 3.182460e-10 --> this is the root value! (first element in theta vector, despite not being an optimum!!!)
          # 1        2.062229e-01
          # 2        2.663341e-01
          # 3        8.812254e-01 note that they have to be in increasing order in Java

res$LogLik # 9.916106

### lnLk2
res <- mvOU(tr, dat, model="OUM", param=list(vcv="fixedRoot", root=FALSE))
res$sigma # 0.006082499
res$alpha # 7.39037
res$theta # root value is assumed to be = optimum1
          # 1        0.2062229
          # 2        0.2663341
          # 3        0.8812254 note that they have to be in increasing order in Java

res$LogLik # 15.52115

### lnLk3
res <- mvOU(tr, dat, model="OUM", param=list(vcv="randomRoot", root=TRUE))
res$sigma # 0.008287661
res$alpha # 10.07163
res$theta # root value is assumed to be = optimum1
          # 1        2.062229e-01
          # 2        2.663341e-01
          # 3        8.812251e-01 note that they have to be in increasing order in Java

res$LogLik # 9.916107

### lnLk4
res <- mvOU(tr, dat, model="OUM", param=list(vcv="randomRoot", root=FALSE))
res$sigma # 0.008287661
res$alpha # 10.07163
res$theta # root value is assumed to be = optimum1
          # 1        2.062229e-01
          # 2        2.663341e-01
          # 3        8.812251e-01 note that they have to be in increasing order in Java

res$LogLik # 9.916107

## (2) 'OUMVNLikelihoodOneTraitTest2'

tr <- read.tree(text="(((sp1:2.0,sp2:1.0):1.0,sp3:4.0):1.0,sp4:3.0);")

res <- mvOU(tr, dat, model="OU1", param=list(vcv="randomRoot", root=TRUE))
res$sigma # 4.601164
res$alpha # 0.8609833
res$theta # theta_0 -46.965464
          # theta_1   6.201598
res$LogLik # -7.63854

res <- mvOU(tr, dat, model="OU1", param=list(vcv="randomRoot", root=FALSE))
res$sigma # 6.841867
res$alpha # 0.7085376
res$theta # theta_0 3.586504
res$LogLik # -8.817273

res <- mvOU(tr, dat, model="OU1", param=list(vcv="fixedRoot", root=TRUE))
res$sigma # 4.003551
res$alpha #  0.7465763
res$theta # theta_0 -33.591241
          # theta_1   6.449917
res$LogLik # -7.630117

res <- mvOU(tr, dat, model="OU1", param=list(vcv="fixedRoot", root=FALSE))
res$sigma # 1.237864
res$alpha # 1.40338e-08
res$theta # theta_0 2.792045
res$LogLik # -8.457486

## (3) 'OUMVNLikelihoodOneTraitTest3'

set.seed(123)

## simulating birth-death tree
lambda <- rexp(1, rate=80) # lambda from exponential (mean = 1/80)
mu <- rexp(1, rate=100) # mu from exponential (mean = 1/100)
tr <- sim.bd.taxa.age(50, 1, lambda, mu, age=100, mrca=TRUE)[[1]]; write.tree(tr) # mrca=TRUE means the process starts at the root (i.e., no root edge)

# arbitrary values below
sigma <- 0.02914682
alpha <- 4.289784
rv <- 2.219841
th <- 1.442039

dat <- mvSIM(tr, param=list(sigma=sigma, alpha=alpha, theta=c(rv,th), root=TRUE), model="OU1")

# for the unit test
paste(tr$tip.label, collapse=",") # "t37,t42,t24,t19,t12,t5,t18,t25,t10,t47,t9,t20,t14,t43,t38,t13,t41,t50,t8,t35,t32,t34,t6,t48,t39,t17,t15,t40,t46,t29,t22,t16,t28,t31,t45,t23,t7,t4,t36,t11,t49,t3,t27,t26,t33,t21,t2,t44,t30,t1"
paste(dat, collapse=",") # 1.46036403223613,1.46684738870839,1.51293402720201,1.38645983622075,1.52360274878184,1.38604516251543,1.34674375692537,1.57061239621065,1.48962576974379,1.48653875766911,1.47882999523348,1.42987407092074,1.500830609038,1.52573554902885,1.35101467964535,1.3754167761214,1.39992186810313,1.47268752575019,1.30540867869118,1.43286309430355,1.35747497708841,1.42191922014755,1.34837424935406,1.5102597527159,1.4259500689268,1.43069041000211,1.52900698556248,1.39163255825086,1.35189645355176,1.40593882022866,1.40399457029073,1.46245982107446,1.42363339595256,1.35978645751607,1.43131557265059,1.4768092041642,1.36373666377089,1.50038165340804,1.48327796797891,1.52949573694115,1.42292947537042,1.4666710783077,1.39114425864009,1.42708261181083,1.38061320641582,1.50472419662598,1.35762598538961,1.47786657493344,1.40219215464719,1.37165305130316

res <- mvOU(tr, dat, model="OU1", param=list(vcv="fixedRoot", root=TRUE))
res$sigma # 0.01507664
res$alpha # 1.988807
res$theta # theta_0 6.082900e-87
          # theta_1 1.435158e+00
res$LogLik # 69.16449

