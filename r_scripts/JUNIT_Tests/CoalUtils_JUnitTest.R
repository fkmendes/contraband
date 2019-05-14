# x is T from Rosenberg
sigma <- function(x, j, i) {
    result = 0
    for (k in j:i) {
        result = result + inside.sigma(x, k, j, i)
    }
    return (result)
}

# geometric progression parenthesis (a_(x))
g.prog.par <- function(base, k) {
    if (k == 0) {
        return (1)
    }
    else {
        z = 0
        prog = 1
        while ((base+z) <= (base+k-1)) {
            prog = prog * (base+z)
            z = z+1
        }
        return (prog)
    }
}

# geometric progression bracket (a_[x])
g.prog.bra <- function(base, k) {
    if (k == 0) {
        return (1)
    }
    else {
        z = 0
        prog = 1
        while ((base-z) >= (base-k+1)) {
            prog = prog * (base-z)
            z = z+1
        }
        return (prog)
    }
}

inside.sigma <- function(x, k, j, i) {
    return (
        exp(-k * (k-1) * x/2) *
            ((2*k-1) * (-1)^(k-j) * g.prog.par(j, (k-1)) * g.prog.bra(i, k))
            /
            (factorial(j) * factorial((k-j)) * g.prog.par(i, k))
    )
}

# finally...
gij <- function(x, j, i) {
    return (sigma(x, j, i))
}

## Now testing!
## i lineages coalesce into j lineages
## i.e., j ancestral lineages, i descendant lineages if forward in time)
# g21, i=2, j=1, x=T=1
gij(1, 1, 2)
1-exp(-1) # 0.6321206

# g22, i=2, j=2, x=T=1
gij(1, 2, 2)
exp(-1) # 0.3678794

# g110, i=10, j=1, x=T=1
gij(1, 1, 10) # 0.2277612

## J. Heled's python code for gij (he calls it pn2k, where n lineages in the present, k in the past)
from numpy import exp
def pn2k(n, k, dt, pop) :
  # n-k rates
  # backwards rates
  rs = [(i*(i-1))/(2) for i in range(k,n+1)]
  print rs
  # work backwards
  # rate of no coal (last)
  cis = [1]
  for i in range(1, n-k+1) :
    r = rs[i]
    for j in range(len(cis)):
      cis[j] *= r/(r-rs[j])
    cis.append(-sum(cis))
  return sum([c*exp((-r/pop)*dt) for r,c in zip(rs,cis)])
