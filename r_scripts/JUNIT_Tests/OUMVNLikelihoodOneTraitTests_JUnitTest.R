# author: Fabio K. Mendes
# This R script gives us the expected values for JUnit tests
# (1) 'OUMVNLikelihoodOneTraitTest1'
# (2) 'OUMVNLikelihoodOneTraitTest2'
# (3) 'OUMVNLikelihoodOneTraitTest3'

library(mvMORPH)
library(TreeSim)
library(FossilSim)

## (1) 'OUMVNLikelihoodOneTraitTest1'

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
res$sigma # 0.02879764
res$alpha # 0.4316411
res$theta # theta_0 -1.924925
          # theta_1   1.035041
res$LogLik # 1.171086

res <- mvOU(tr, dat, model="OU1", param=list(vcv="randomRoot", root=FALSE))
res$sigma # 0.09143114
res$alpha # 0.599265
res$theta # theta_0 0.4152632
res$LogLik # -0.5143806

res <- mvOU(tr, dat, model="OU1", param=list(vcv="fixedRoot", root=TRUE))
res$sigma # 0.009103832
res$alpha # 1.818089e-11
res$theta # theta_0 -3.726855e-01
          # theta_1 1.142785e+10
res$LogLik # 1.367972

res <- mvOU(tr, dat, model="OU1", param=list(vcv="fixedRoot", root=FALSE))
res$sigma # 0.01940518
res$alpha # 2.707329e-11
res$theta # theta_0 0.3497826
res$LogLik # -0.1459122

## (3) OUMVNLikelihoodOneTraitTest3

## simulating fossilized birth-death tree
set.seed(123)
lambda <- rexp(1, rate=80)[1]
mu <- rexp(1, rate=100)[1]
psi <- rexp(1, rate=150)[1]

tr <- sim.fbd.taxa(50, 1, lambda, mu, psi, complete=TRUE)[[1]] # 50 species, many fossil/SA tips
## (t72_1:33.95517374,((t45_2:105.7970875,t45_1:0):61.92884103,((t16_2:64.31266032,t16_1:0):33.18670286,((((((((t33_2:17.24690108,t33_1:0):34.29978643,((((t29_5:8.042489968,t29_4:0):48.76232459,t29_3:0):2.870515606,t29_2:0):113.0118958,t29_1:0):78.8636162):2.752286226,t104_1:0):13.2284639,t76_1:19.98684163):35.40421183,t99_1:0):22.37220403,(((((((t25_2:5.545398329,t25_1:0):29.73789406,(((t14_1:19.83955306,(t68_2:7.642029053,t68_1:0):12.197524):4.493086625,t15_1:24.33263968):7.659951961,t32_1:31.99259164):4.589970848):16.34481241,t18_1:35.83487487):68.91853002,((t41_1:26.38487177,(((t39_1:21.69802025,t20_1:21.69802025):22.43410663,t160_2:0):0.5011266949,t160_1:0):24.01552498):43.01583713,t138_1:0):10.18128923):70.28879289,((((t40_3:5.674797941,t40_2:0):85.82961996,t40_1:0):9.144286626,((t22_1:6.416294394,t62_1:6.416294394):52.92804634,t169_1:0):41.3043638):44.11984404,t126_1:0):47.36614923):93.92949444,t113_1:0):24.15812702,(((((t7_1:96.64737931,(t56_1:9.124505504,t60_1:9.124505504):87.5228738):142.4776409,t48_1:53.57593652):13.68925583,((((((t5_1:1.47740974,t24_1:1.47740974):67.36924043,t172_1:0):0.6094754868,((t81_3:23.44109494,t81_2:0):11.8683594,t81_1:0):34.14667132):94.37344028,(((t42_1:24.136168,((t78_2:67.1575288,t78_1:0):6.07403593,t13_1:11.99646428):19.79173068):16.01164487,t8_1:15.88519099):53.96855772,(((((t64_2:57.50834491,t64_1:0):14.48295323,(t44_1:66.77751843,((t19_1:12.37766384,(t3_2:21.04342944,t3_1:0):34.8018801):10.2073654,t146_1:0):0.7248434878):5.213779703):6.124517947,(((((t85_3:30.86302382,t85_2:0):2.746255125,t85_1:0):7.42992738,(t23_2:35.51946648,t23_1:0):5.519739842):5.318122342,(t63_1:29.4370073,t4_1:29.4370073):16.92032137):10.87142916,(((t27_2:0.7705516481,t27_1:0):16.68364597,t65_1:17.45419762):0.8973042275,(t30_1:5.872210534,(t9_1:5.379640887,(t49_1:0.03872308848,t54_1:0.03872308848):5.340917798):0.4925696479):12.47929131):38.87725598):20.88705826):57.37473544,t132_2:0):11.64855352,t132_1:0):15.86439296):0.8260679418):12.34991723,(t67_2:54.95630021,t67_1:0):46.89744564):72.69817908,t115_1:0):3.936613795):37.64004304,((((t73_3:142.7476486,t73_2:0):23.26030434,t73_1:0):33.70759497,(t86_1:29.84643822,((t11_2:52.84278453,t11_1:0):65.0721915,t31_1:54.03135372):52.51395311):30.03845768):0.7587008348,t100_1:0):32.43803262):10.22767547,t6_1:74.01326892):9.540324705):15.08568888):1.038342891,t92_1:0):8.58602375,(((t52_1:16.45509171,(((t46_1:148.9059947,t66_1:169.3974901):0.6996035147,t80_1:71.3880344):1.442748076,t109_1:0):8.282393923):76.16579002,(((t47_1:87.75883092,(((t37_1:19.21090631,t71_1:19.21090631):2.716558918,t162_1:0):52.25715577,(t50_2:14.18543288,t50_1:0):0.7051498905):13.57420992):70.93504021,(t84_1:91.53880307,(((t43_3:45.6907903,t43_2:0):27.11086441,t43_1:0):26.97698465,t82_1:57.43850143):42.53791402):16.37731775):67.82194637,t120_1:0):70.63623304):19.65198808,(((t12_1:85.77104457,((((((t26_1:6.654829005,((t36_1:12.77851341,t70_1:12.77851341):20.88269631,(t53_1:31.23167166,((t74_1:11.18633708,t2_1:11.18633708):13.86162972,(t59_2:8.290822624,t59_1:0):16.75714418):6.183704852):2.429538064):33.88601089):65.52318133,((t17_1:47.74620153,t77_1:47.74620153):26.04800054,((t87_1:24.52842258,t69_1:24.52842258):40.93917936,t79_1:65.46760194):8.326600137):59.27619986):2.038654818,t34_1:33.0666319):2.275848245,t122_1:0):59.7289996,((((t57_2:9.620141517,t57_1:0):45.82168925,((t58_2:5.178405873,t58_1:0):55.47396517,((t21_1:14.00132991,t28_1:14.00132991):36.16842295,t165_1:0):10.48261818):3.289390269):95.35491899,t142_1:0):22.11781524,t35_1:6.324292886):15.69940907):2.23954769,((t83_1:0.07202234145,((t1_1:4.659677209,(((t61_4:25.59732606,t61_3:0):8.674379048,t61_2:0):4.825121609,t61_1:0):4.738801364):14.43560063,t150_1:0):20.9778315):78.47955114,t131_1:0):41.62484094):37.0609):69.68067928,t106_1:0):0.7776192694,((t10_1:9.54139831,(t38_2:26.07965526,t38_1:0):76.6686444):85.31864731,(((t75_1:63.31636704,t55_1:21.6953316):38.78938204,t143_1:0):151.8995115,t51_1:49.73077652):1.996244886):50.8711454):9.931387773):18.12833616):44.30923714):67.66883573):36.06056073):70.8162184;
## tr <- paintSubTree(tr, node=165, state=2, stem=FALSE)

dat <- mvSIM(tr, model="OU1", param=list(theta=c(0, 0.2), sigma=0.1, alpha=1, root=TRUE))

res <- mvOU(tr, dat, model="OU1", param=list(vcv="fixedRoot", root=TRUE))
res$sigma # 0.06733451
res$alpha # 0.7972125
res$theta # theta_0     1.960019e-15
          # theta_1     1.596677e-01
res$LogLik # 25.41701

## (4) OUMVNLikelihoodOneTraitTest4 (small non-ultrametric tree with sampled ancestor and fossil tip)

tr <- read.tree(text="(((sp1:2.0,sp2:1.0):1.0,sp3:0.0):1.0,sp4:3.0);")
tr <- paintSubTree(tr, node=6, state=2, stem=TRUE)

set.seed(123)
dat <- mvSIM(tr, model="OUM", param=list(theta=c(0, 0.2, 0.4), sigma=0.1, alpha=1, root=TRUE))

res <- mvOU(tr, dat, model="OUM", param=list(vcv="fixedRoot", root=TRUE))
res$sigma # 1.11256e-08
res$alpha # 0.5564338
res$theta # theta_0 0.87579783
          # theta_1 0.05779027
          # theta_2 0.19382641
res$LogLik # 31.19386
