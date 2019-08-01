library(HDInterval)

get.95 <- function(a.vector) {
    res = hdi(a.vector, cresMass=.95)

    return(c(res[[1]], res[[2]], mean(a.vector)))
}

write.shell.script <- function(shell.scripts.path, sim.idx, time.needed, job.prefix, jar.path, xml.file.path, xml.file.prefix) {
    shell.file.name = paste0(shell.scripts.path, xml.file.prefix, sim.idx, ".sh")
    if (file.exists(shell.file.name)) {
        file.remove(shell.file.name)
    }

    xml.file.name = paste0(xml.file.path, xml.file.prefix, sim.idx, ".xml")

    write(file=shell.file.name, paste(
        paste0("#!/bin/bash -e\n#SBATCH -J ", job.prefix, sim.idx),
        "#SBATCH -A nesi00390",
        paste0("#SBATCH --time=", time.needed),
        "#SBATCH --mem-per-cpu=1G",
        "#SBATCH --cpus-per-task=1",
        "#SBATCH --ntasks=1",
        "#SBATCH --hint=nomultithread",
        "#SBATCH -D ./",
        paste0("#SBATCH -o ", job.prefix, sim.idx, "_out.txt"),
        paste0("#SBATCH -e ", job.prefix, sim.idx, "_err.txt"),
        paste0("\nsrun java -jar ", jar.path, " ", xml.file.name),
        sep="\n")
        )
}

get.plot <- function(x.name, y.name, x.min, x.max, y.min, y.max, x.lab, prior.mean, data.df, plot.hdi) {
    lower = data.df[,paste0("lower.",x.name)]
    upper = data.df[,paste0("upper.",x.name)]
    x = data.df[,x.name]; y = data.df[,y.name]
    reg.df = data.frame(cbind(x,y,lower,upper,plot.hdi))


    print((x.max+x.min)/2)

    pl = ggplot() + geom_point(mapping=aes(x=x, y=y), shape=20) +
    coord_cartesian(ylim=c(y.min, y.max)) +
    xlab(x.lab) + ylab("Posterior mean") +
    geom_abline(slope=1, linetype="dotted") +
    geom_abline(slope=0, intercept=prior.mean, color="blue") +
    theme(
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        plot.title = element_text(hjust=0.5),
        axis.line = element_line(),
        axis.ticks = element_line(color="black"),
        axis.text.x = element_text(color="black", size=10),
        axis.text.y = element_text(color="black", size=10),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12)
    ) + scale_x_continuous(labels = function(x) round(as.numeric(x), digits=3), breaks=c(x.min,((x.max+x.min)/2),x.max), limits=c(x.min,x.max)) +
        scale_y_continuous(labels = function(x) round(as.numeric(x), digits=3), breaks=c(y.min,((y.max+y.min)/2),y.max), limits=c(y.min,y.max)) +
        geom_linerange(data=reg.df[plot.hdi,], mapping=aes(x=x, ymax=upper, ymin=lower), color="red", alpha=.4, size=1.5) +
        geom_linerange(data=reg.df[!plot.hdi,], mapping=aes(x=x, ymax=upper, ymin=lower), color="lightgray", alpha=.4)
    return(pl)
}

get.plot.no.hdi <- function(x.name, y.name, x.min, x.max, y.min, y.max, x.lab, prior.mean, data.df) {
    lower = data.df[,paste0("lower.",x.name)]
    upper = data.df[,paste0("upper.",x.name)]
    x = data.df[,x.name]; y = data.df[,y.name]
    reg.df = data.frame(cbind(x,y,lower,upper,plot.hdi))

    pl = ggplot() + geom_point(mapping=aes(x=x, y=y), shape=20) +
    coord_cartesian(ylim=c(y.min, y.max)) +
    xlab(x.lab) + ylab("MLE") +
    geom_abline(slope=1, linetype="dotted") +
    geom_abline(slope=0, intercept=prior.mean, color="blue") +
    theme(
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(),
        plot.title = element_text(hjust=0.5),
        axis.line = element_line(),
        axis.ticks = element_line(color="black"),
        axis.text.x = element_text(color="black", size=10),
        axis.text.y = element_text(color="black", size=10),
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12)
    ) + scale_x_continuous(labels = function(x) round(as.numeric(x), digits=3), breaks=c(x.min,(x.max+x.min)/2,x.max), limits=c(x.min,x.max)) +
    scale_y_continuous(labels = function(x) round(as.numeric(x), digits=3), breaks=c(y.min,(y.max+y.min)/2,y.max), limits=c(y.min,y.max))
    return(pl)
}

flip.trace.var <- function(a.df, cols) {
    a.df.res = a.df
    for (i in 1:nrow(a.df)) {
        for (col in cols) {
            tmp = a.df[i,col]
            ## print(tmp)
            ## print(a.df[i,col])
            if (tmp > a.df[i,col+1]) {
                ## cat("flipped.")
                a.df.res[i,col] = a.df[i,(col+1)]
                a.df.res[i,(col+1)] = tmp
            }
        }
    }
    return(a.df.res)
}

flip.w.bool <- function(a.df, bool.vec, col1, col2) {
    a.df.res = a.df
    for (i in 1:nrow(a.df.res)) {
        if (bool.vec[i]) {
            print(paste0("col 1 before = ", a.df.res[i,col1]))
            print(paste0("col 2 before = ", a.df.res[i,col2]))
            tmp = a.df[i,col2]
            a.df.res[i,col2] = a.df[i,col1]
            a.df.res[i,col1] = tmp
            print(paste0("col 1 after = ", a.df.res[i,col1]))
            print(paste0("col 2 after = ", a.df.res[i,col2]))
        }
    }
    return(a.df.res)
}

