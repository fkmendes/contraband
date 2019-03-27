library(stringr)
library(geiger)
library(ggplot2)
library(plotrix)

# --- Simulation functions --- #
make_cmd_str <- function(mspath=NULL, rep.i=NULL, ils.i=NULL, nspp=NULL, spop=NULL, nloci=NULL, dfsplits=NULL, shrinkfactor=NULL, seeds=NULL) {
        
    ## preparing string for splits (for when anc Ne's are made larger)
    splits.str <- character()
    for (i in 1:length(dfsplits)) {
        splits.str = paste(splits.str, "-ej", dfsplits[[i]][1], dfsplits[[i]][2], dfsplits[[i]][3], sep=" ")
    }
    splits.str = str_trim(splits.str, "l")
    # print(splits.str) # checking

    ## preparing string for pop size changes (comment block if not changing pop size)
    shrink.str <- character()
    shrink.str = paste(shrink.str, "-en", dfsplits[[1]][1], "1", shrinkfactor, sep=" ") # it won't go back to the original size because we're simulating just the ingroup
    shrink.str = str_trim(shrink.str, "l")
    ## print(shrink.str) # checking

    output.file.name = paste0("data/ils", ils.i-1, "_", nloci, "loci_", spop, "ind_rep", rep.i, ".txt") # leave uncommented
    cmd.str <- paste(mspath, nspp*spop, nloci, "-I", nspp, paste(rep(spop, nspp), collapse=" "), "-t 4", splits.str, shrink.str, "-seeds", paste(seeds, collapse=" "), "-T >", output.file.name, sep=" ")
    ## print(cmd.str) # checking
    
    return(c(cmd.str, output.file.name))
}
## testing:
## make_cmd_str(mspath=ms.path, rep.i=1, ils.i=1, nspp=n.spp, spop=s.pop, nloci=n.loci, dfsplits=split.5, shrinkfactor=5.2, seeds=seeds.list.100[[1]]) # for changing anc Ne's

batch_ms <- function(msexec=NULL, rep.i=NULL, ils.i=NULL, nspp=NULL, spop=NULL, nloci=NULL, dfsplits=NULL, shrinkfactor=NULL, seeds=NULL, kernel=rnorm, ...) {

    ## initializing return objects
    ss.mtx = vector("list", nspp) # 1st level: species, 2nd level: locus, 3rd level: haplotypes
    for (i in 1:length(ss.mtx)) {
        for (j in 1:nloci) {
            ss.mtx[[i]][[j]] <- vector("list", length=spop)
        }
    }
    ## print(ss.mtx) # checking...

    ## read in ms output
    cmdstr.outfile = make_cmd_str(mspath=msexec, rep.i=rep.i, ils.i=ils.i, nspp=nspp, spop=spop, nloci=nloci, dfsplits=dfsplits, shrinkfactor=shrinkfactor, seeds=seeds)

    ## print(cmdstr.outfile[1]) # printing ms command

    system(cmdstr.outfile[1]) # running ms
    ms.out = readLines(cmdstr.outfile[2])
    write(sprintf("Replicate %d, nl=%d, nspp=%d", rep.i, nloci, nspp), "")

    nrow = length(ms.out)
    ind.idx = 0 # FKM
    sp.idx = 0 # FKM
    locus.idx = 0
    skip.lines = TRUE
    nsegsites = c()

    for (i in 1:nrow) {
        line = ms.out[i]
        first.char = substr(line, 1, 1)

        ## new sim
        if (substr(line,1,2) == "//") {
            skip.lines = FALSE
            locus.idx <- locus.idx + 1
            sp.idx <- 1
        }

        ## skip ms header
        else if (skip.lines == TRUE) {
            next
        }

        ## num segsites (FKM: won't use it)
        else if (first.char == "s") {
            next
        }

        ## segsites
        else if (first.char == "0" || first.char == "1") {
            ind.idx <- ind.idx + 1
            sample_segsites = as.numeric(unlist(strsplit(line,"")))
            ss.mtx[[sp.idx]][[locus.idx]][[ind.idx]] = sample_segsites

            # FKM: prepping for next sp
            if (ind.idx%%spop == 0) {
                sp.idx <- sp.idx + 1
                ind.idx <- 0
            }
        }
    }

    ## print(ss.mtx) # checking
    return(ss.mtx)
}
## testing
## m <- batch_ms(msexec=ms.path, rep.i=1, nspp=n.spp, spop=s.pop, nloci=n.loci, dfsplits=split.5, shrinkfactor=5.2, kernel=rnorm, seeds=seeds.list.100[[1]])    

header <- c("rep", "ils", "sp", "ind", "value")
batch_ms_wrapper <- function(nrep=NULL, nloci=NULL, nspp=NULL, spop=NULL, dfsplits=NULL, shrinkfactors=NULL, seeds=NULL, kernel=NULL, scale=NULL, mspath=NULL, alpha=2.0, alphasn=0.0, xisn=0.0, outfiledir="data/batch_ms_out/") {

    traits.no.epi <- vector("list", 5) # 5 ils conditions
    for (ils in 1:5) {        
        for (r in 1:nrep) {
            reps.this.ils = data.frame(matrix(nrow=spop*nspp,ncol=length(header)))
            names(reps.this.ils) = header
            reps.this.ils$rep = rep(r, nspp*spop)
            reps.this.ils$sp = rep((1:nspp), each=spop)
            reps.this.ils$ind = rep((1:spop), times=nspp)
            reps.this.ils$ils = ils
            reps.this.ils$value = rep(0, times=spop*nspp)

            alleles = batch_ms(msexec=ms.path, rep.i=r, ils.i=ils, nspp=nspp, spop=spop, nloci=nloci, dfsplits=dfsplits, shrinkfactor=shrinkfactors[ils], seeds=seeds[[r]], kernel=kernel) # 1st level = sp, 2nd level = locus, 3rd level = individual
            for (loc in 1:nloci) {
                n.segsites = length(alleles[[1]][[loc]][[1]]) # all sp and ind have same # of seg sites, so can take n.segsites from 1st sp and 1st ind
            
                scale.per.locus = scale/sqrt(nloci) # drawing 1 mut effect PER LOCUS

                mut.effects = kernel(n=n.segsites, sd=scale.per.locus) # getting the effect of the "1" potentially observed at each SEG SITE
                ## print(paste0("Mutation effects for locus ", loc, " are: ")); print(mut.effects) # checking...

                ## if I want to have a look at the effects (writes a bunch of files and requires a folder batch_ms_out
                ## outeffects.path <- paste0(outfile_dir, "effects_locus", loc, "_l", nloci, "_norm_ils", ils, "_", r, ".txt")
                ## cat(mut.effects, file=outeffects.path, sep="\n")

                for (sp in 1:nspp) {
                    ind.trait = rep(0, spop) # will store the trait value for individuals of this sp 
                
                    seg.sites.mtx = matrix(unlist(alleles[[sp]][[loc]]), ncol = n.segsites, byrow=TRUE)
                    ## print(seg.sites.mtx) # checking against ms output
                    
                    for (ind in 1:spop) {
                        ## print(paste0("Alleles from ind ", ind, " are ")); print(seg.sites.mtx[ind,])
                        loc.effect = seg.sites.mtx[ind,] %*% mut.effects
                        ## print(paste0("The effect of locus ", loc, " on individual ", ind, " from sp ", sp, " is ", loc.effect))
                        ind.trait[ind] = ind.trait[ind] + loc.effect # one value
                        reps.this.ils[reps.this.ils$sp==sp,]$value[ind] = reps.this.ils[reps.this.ils$sp==sp,]$value[ind] + ind.trait[ind]
                    }
                    
                    ## print(reps.this.ils) # checking
                    ## print(paste0("The trait values from locus ", loc, " of individual(s) from sp ", sp, " is/are ", ind.trait)) # checking
                }
            }

            ## print(paste0("The final trait values of individual(s) of ILS=", ils, " replicate ", r, " is/are:")); print(reps.this.ils) # checking
            
            traits.no.epi[[ils]][[r]] = reps.this.ils

            ## if I want to check output files (prints a bunch of files)
            outvalue.path = paste0(outfiledir, "l", nloci, "_norm_ils", ils-1, "_", r, ".txt")
            write.table(reps.this.ils, file=outvalue.path, row.names=FALSE, col.names=TRUE, quote=FALSE)
        }
    }

    ## print(traits.no.epi) # checking
    return(traits.no.epi)
}

# --- PCM analyses functions --- #
get.rate <- function(dataset, tree) {
    rates <- vector("list", length=5)
	
    for (ils in 1:5) {
        for (i in 1:length(dataset[[ils]])) {
            values <- dataset[[ils]][[i]]$value
            names(values) <- dataset[[ils]][[i]]$sp
            fit.bm <- fitContinuous(tree, values, model="BM")
            
            rates[[ils]][[i]] <- fit.bm$opt$sigsq
            print(paste0("ILS=", ils, " rep=", i))
        }
    }

    return(rates)
}

make.mean.rate.df <- function(fulldataset, nloci) {
    ggplot.df = data.frame(matrix(ncol=4, nrow=5)) # col: mean, se, ils, nloci; row: each ils

    i = 1
    for (ils in 1:5) {
        ggplot.df[i,] = c(mean(fulldataset[[ils]]), std.error(fulldataset[[ils]]), nloci, ils)
        i <- i+1
    }
    names(ggplot.df) = c("mean.rate", "se", "nloci", "ILS")
    return(ggplot.df)
}

# --- Sanity check functions --- #
mle.analytical <- function(obs.data, sigma.mat) {
    if (dim(sigma.mat)[1] != dim(sigma.mat)[2]) { stop("Var-cov matrix is not a square matrix. Exiting...") }
    if (length(obs.data) != dim(sigma.mat)[1]){ stop("Number of observed data points and var-cov matrix dimensions don't match. Exiting...") }
  
    n = length(obs.data)
    ones = rep(1, n)
    z0 = solve(t(ones)%*%solve(sigma.mat)%*%ones) %*% (t(ones)%*%solve(sigma.mat)%*%obs.data) # add LaTeX formula here
    vecz0 = rep(z0, n)

    sigmasq = t(obs.data - vecz0)%*%solve(sigma.mat)%*%(obs.data - vecz0)/n # add LaTeX formula here
    likelihood = exp(-0.5*t(obs.data - vecz0)%*%solve(sigma.mat)%*%(obs.data - vecz0)/sigmasq)/((2*pi*sigmasq)^(n/2)*det(sigma.mat)^(1/2)) # add LaTeX formula here
  
    res = c(z0, sigmasq, log(likelihood))
    names(res) = c("mean", "sigma.mat", "lnL")
  
    return(res)
}

# --- Functions for obtaining "approximate" vcv matrix from tree --- #
amplify.vcvmat <- function(mat) {
    res = cbind(rbind(mat,0),0) # add a bottom row and right column of 0's
    diag(res) = mat[1,1] # all diags are the same as variance within species 1
    res
}

## "approximates" matrix in the limit (if cov12 != cov21, pick min, max or avg of both)
mat.limit <- function(mat, method = "min") {
    ut = upper.tri(mat, diag=FALSE)
    ut[ut == FALSE] = NA # everything but upper triangle is NA
    to.change = mat*ut # upper triangle-populated matrix
    diagonal = mat[1,1]*diag(dim(mat)[1]) # put variance within species 1 along diagonal

    if (method == "min") {
        index = apply(to.change[,-1], 2, function(x){ unique(min(x, na.rm = T)) }) # 2 indicates apply to column
    }
    
    if (method == "max") {
        index = apply(to.change[,-1], 2, function(x){ unique(max(x, na.rm = T)) })
    }
    
    if (method == "average") {
        index = apply(to.change[,-1], 2, function(x){ unique(mean(x, na.rm = T)) })
    }
    
    for (i in 2:dim(mat)[1]) { ut[,i] = ut[,i]*index[i-1] }
  
    ut[which(is.na(ut))] = 0
    res = ut + t(ut) + diagonal
  
    return(res)
}

## converts var-cov matrix into distance matrix for upgma
vcvmat.to.upgmamat <- function(vcvmat) {
    d = matrix(data=1, nrow=dim(vcvmat)[1], ncol=dim(vcvmat)[2])
    diag(d) = 0
    diagonal = vcvmat[1,1]
    diag(vcvmat) = 0
    res = 2*(diagonal*d - vcvmat)

    return(res)
}
    
## gets distance matrix by calling functions above, then calls upgma to get a tree to then prune with (here we get the "correct" vcv matrix, "approximate" it to what it'd be in the coalescent limit, then convert it into a distance matrix, then call upgma on it to get a tree -- then use the tree in pruning)
get.corrected.rate.pruning <- function(mats, method="min") {

    res = vector("list", length=5)
  
    for (ils in 1:5) {
        for (r in 1:100) {
            m = amplify.vcvmat(mats[[ils]][[r]]) # adds phantom branch entries to matrix
            m = mat.limit(m, method) # "approximates" matrix in the limit (if 12 != 21, pick min, max or avg of both)
            m = vcvmat.to.upgmamat(m) # makes distance matrix out of vcv matrix
      
            tree = upgma(m, method = "average")
            this.data = c(l100.5.rnorm[[ils]][[r]][,5], 0)
            names(this.data) <- c("1", "2", "3", "4", "5", "6")
      
            fit.res = fitContinuous(phy = tree, dat = this.data, model = "BM")
            res[[ils]][r] = fit.res$opt$sigsq
        }
    }
    
    return(res)
}
