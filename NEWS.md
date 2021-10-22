# antigen.garnish 2.3.0

Bugfix:

- update NetMHC tool data directory paths
- fix permissions on configuration script (was non-executable)
- minor documentation improvements

# antigen.garnish 2.2.0

Bugfix:

- remove netMHC binary check from foreignness_score function that does not require netMHC binaries

# antigen.garnish 2.1.0

Added functionality:

- add Dockerfile used to build publically available Docker image (andrewrech/antigen.garnish)
- relax parsing of SnpEff annotations to reduce errors across VCFs and/or SnpEff versions
- simplify internal functions that configure netMHC tools
- use Docker to run tests

# antigen.garnish 2.0.0

Breaking changes:

- new filtering criteria based on TESLA: https://www.parkerici.org/research-project/tumor-neoantigen-selection-alliance-tesla/, see ?garnish_antigens
- update Ensembl transcript CDS database to release 101 (any prior versioned transcript release is acceptable)
- remove support for mitochondrial CDS, xlsx input, JAFFA input, PureCN input to focus on peptide/transcript tables and VCFs
- improve SnpEff annotation parsing

Additional updates:

- allow setting thread use via environment variable AG_THREADS
- include percentile rank filtering at default netMHC thresholds
- simplify installation instructions
- remove 7 dependencies
- additional integration tests
- do not export internal functions

# antigen.garnish 1.1.1

- update to include configuration for netMHC dependencies
- support for paired wild-type and mutant protein level input

# antigen.garnish 1.1.0

- added dissimilarity_score
- updated readme citation, link to manuscript repository
- externalized smith-waterman alignment function
- appropriate test updates
- remove deprecated garnish_fitness function
- update to LICENSE

# antigen.garnish 1.0.0

- prepare for the release of the antigen.garnish manuscript
- multiple metrics of antigen quality are now computed
- dissimilarity_score dummy function present
- pulled most unnecessary mclapply loops to prevent parallelization failures
- functionalized foreignness_score code into external function
- garnish_plot, garnish_summary now incorporate antigen quality metrics
- garnish_antigens now returns ranked antigens by quality metrics
- numerous edge case, efficiency updates and bug fixes

# antigen.garnish 0.0.6

- clonality filter added to garnish_affinity
- RNA expression filter added to garnish_affinity
- fitness model implemented in R vs. Python
- garnish_antigens function added to rank neoantigens
- garnish_score function added to return sample level immune fitness summary
- updated summary and plot functions with better output
- added wiki with installation instructions
- numerous edge case bug fixes
- remove biomaRt dependencies

# antigen.garnish 0.0.5

- add ncbi-blast functionality to determine neoantigen near matches
- add ncbi-blast functionality to determine known IEDB matches
- implement antigen fitness model of Lukza et al. _Nature_ 2017

# antigen.garnish 0.0.4

- add summary plots
- improve test coverage

# antigen.garnish 0.0.3

- inverse match against global normal proteins
- improve testing, test formatting

# antigen.garnish 0.0.2

- add prediction deduplication to garnish_affinity
- add gene fusions as a source of neoantigens using garnish_jaffa
- stable API

# antigen.garnish 0.0.1

- initial version
