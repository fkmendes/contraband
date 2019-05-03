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

get.plot <- function(x.name, y.name, x.min, x.max, y.min, y.max, x.lab, prior.mean, data.df) {
    x = data.df[,x.name]; y = data.df[,y.name]
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
    ) + scale_x_continuous(labels = function(x) round(as.numeric(x), digits=3), limits=c(x.min,x.max))
    return(pl)
}
