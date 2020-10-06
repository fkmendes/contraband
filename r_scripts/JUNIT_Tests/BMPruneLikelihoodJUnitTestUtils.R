calculate.correlations <- function(rate.matrix) {
  trait.nr = ncol(rate.matrix)
  res = c()
  
  idx = 1
  for (i in 1:(trait.nr - 1)) {
    for (j in (i + 1) : trait.nr){
      res[idx] = rate.matrix[i, j] / (sqrt(rate.matrix[i, i]) * sqrt(rate.matrix[j, j]))
      idx = idx + 1
    }
  }
  
  return (res)
}

populate.upper.matrix <- function(trait.nr, sigmasq, correlation){
  res = diag(trait.nr)
  m = 1
  for (i in 1:trait.nr){
    res[i,i] = sigmasq[i]
    if(i < trait.nr) {
      for (j in (i+1):trait.nr) {
        res[i,j] = sqrt(sigmasq[i]) * sqrt(sigmasq[j]) * correlation[m]
        m = m + 1
      }
    }
  }
  return(res)
}

write.script.for.mcmcTree <- function(main.path, file.name, names.list, ages.list, M, shrinkage.correlation.matrix, this.tree, use.tip.date, burn.in, sample.freq, sample.nr, popVar){
  if(is.null(popVar)){
    if(is.null(names.list)){
      write.morpho(M = M, filename = paste0(main.path, "input/", file.name, ".aln"), names = this.tree$tip.label, age = ages.list, R = shrinkage.correlation.matrix, method = "chol")
    } else {
    write.morpho(M = M, filename = paste0(main.path, "input/", file.name, ".aln"), names = names.list, R = shrinkage.correlation.matrix, method = "chol")
    }
  }
  else {
    #if(is.null(names.list)){
      #write.morpho(M = M, filename = paste0(main.path, "input/", file.name, ".aln"), names = this.tree$tip.label, age = ages.list, R = shrinkage.correlation.matrix, method = "chol")
    #} else {
      write.morpho(M = M, filename = paste0(main.path, "input/", file.name, ".aln"), names = ages.list, R = shrinkage.correlation.matrix, c = popVar, method = "chol")

  }

  tree.topolpy <- this.tree
  tree.topolpy$edge.length <- NULL
  write.tree(tree.topolpy, file = paste0(main.path, "input/", file.name, "_topology.trees"))
  treeMCMCtree(tree = paste0(main.path, "input/", file.name, "_topology.trees"), aln = paste0(main.path, "input/", file.name, ".aln"), filename = paste0(main.path, "input/", file.name, "_divtimes.trees"))
  
  
  if (file.exists(paste0(main.path, "input/", file.name, ".ctl"))) {
    file.remove(paste0(main.path, "input/", file.name, ".ctl"))
  }
  
  template.lines = readLines(paste0(main.path,"input/ControlFile_template.ctl"))
  for (line in template.lines) {
    line = gsub("\\[SeqFilePathHere\\]", paste0(main.path, "input/", file.name, ".aln"), line)
    line = gsub("\\[TreeFilePathHere\\]", paste0(main.path, "input/", file.name, "_divtimes.trees"), line)
    line = gsub("\\[OutFilePathHere\\]", paste0(main.path, "output/", file.name, "_out.txt"), line)
    line = gsub("\\[MCMCFilePathHere\\]", paste0(main.path, "output/", file.name, "_mcmc.txt"), line)
    line = gsub("\\[TipDateHere\\]", use.tip.date, line)
    line = gsub("\\[ClockHere\\]", "1", line)
    line = gsub("\\[BurnInHere\\]", burn.in, line)
    line = gsub("\\[SampFreqHere\\]", sample.freq, line)
    line = gsub("\\[NSampleHere\\]", sample.nr, line)
    write(line, file=paste0(main.path, "input/", file.name, ".ctl"), append=TRUE)
  }
}

get.tip.ages.str <- function(tr) {
  depths = node.depth.edgelength(tr)
  tr.height = max(depths)
  fossil.idxs = as.numeric(as.character(depths[1:length(tr$tip.label)])) < as.numeric(as.character(tr.height))
  tip.ages = tr.height - depths[1:length(tr$tip.label)]
  ages.list = c()
  for (i in 1:length(tr$tip.label)) {
    if(fossil.idxs[i]) {
      # it is a fossil
      ages.list[i] = paste0(tr$tip.label[i], "^", depths[i])
    } else {
      # it is an extant species
      ages.list[i] = paste0(tr$tip.label[i], "^", tr.height)
    }
  }
  return(ages.list)
}

get.tip.ages <- function(tr) {
  depths = node.depth.edgelength(tr)
  tr.height = max(depths)
  fossil.idxs = as.numeric(as.character(depths[1:length(tr$tip.label)])) < as.numeric(as.character(tr.height))
  tip.ages = tr.height - depths[1:length(tr$tip.label)]
  return(tip.ages )
}

write.bash.4.mcmcTree <- function(mcmcTree.path, main.path, file.name) {
  header = "#!/bin/sh"
  to.path = paste0("cd ", mcmcTree.path)
  run.command = paste0("./mcmctree ", main.path, "input/", file.name, ".ctl")
  res = paste0(c(header, to.path, run.command), collapse = "\n")
  write(res, file=paste0(main.path, "bashscript/", file.name, ".sh"))
}

print.trait.values <- function(dat){
  res = c()
  for(i in 1:dim(dat)[1]) {
    res = c(res, paste0(dat[i,], collapse = ", "))
  }
  return(paste0(res, collapse = ", "))
}

print.trait.values4.xml <- function(dat){
  res = c()
  for(i in 1:dim(dat)[1]) {
    res = c(res, paste0(dat[i,], collapse = " "))
  }
  return(paste0(res, collapse = " "))
}

get.extant.str <- function(this.tree){
  ages.list = vector("list", length(this.tree$tip.label))
  for (i in 1:length(this.tree$tip.label)) {
    # all tips are extant species
    ages.list[[i]] = 0
  }
  names(ages.list) = this.tree$tip.label
  return(ages.list)
}

