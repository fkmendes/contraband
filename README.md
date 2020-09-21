# contraband

The 'contraband' (**cont**inuous **tra**its **b**rowni**an** mo**d**els) BEAST 2 package provides a suite of models for the study of continuous trait evolution in a Bayesian framework.

## Citation

[[ref]](https://besjournals.onlinelibrary.wiley.com/doi/epdf/10.1111/2041-210X.13458) Théo Gaboriau, Fábio K. Mendes, Simon Joly, Daniele Silvestro and Nicolas Salamin (2020). A multi-platform package for the analysis of intra- and interspecific trait evolution. *Methods Ecol. Evol.*, 1-9.

## Building

In order to build an executable .jar file for the *contraband* package:

+ (1) Clone and build dependencies:
    + (1.1) BEAST2 from [here](https://github.com/CompEvol/beast2);
    + (1.2) BEASTLabs from [here](examples/testing/BMMVNLikelihoodOneTrait_fixedtree.xml);
    + (1.3) sampled-ancestors from [here](https://github.com/CompEvol/sampled-ancestors);

+ (2) Clone the *contraband* repository side-by-side with the dependencies' (above) directories resulting from step (1);

+ (3) From *contraband/*, type:

    ```
    $ ant
    ```

    This command invokes the *ant* tool, which executes the building instructions inside *build.xml*. The new *contraband* .jar file will be put inside *build/dist/*.

## Running an analysis with *bdtree*

From your shell or Terminal, type:

```
$ /path/to/build/dist/contraband.v0.0.1.jar /path/to/contraband/examples/testing/BMMVNLikelihoodOneTrait_fixedtree.xml
```
