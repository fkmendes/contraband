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
        "#SBATCH --mem-per-cpu=12288",
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
