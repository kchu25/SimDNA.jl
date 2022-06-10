# SimDNA

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://kchu25.github.io/SimDNA.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://kchu25.github.io/SimDNA.jl/dev/)
[![Build Status](https://github.com/kchu25/SimDNA.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/kchu25/SimDNA.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/kchu25/SimDNA.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/kchu25/SimDNA.jl)


Yes, we can generate synthetic DNA sequence motifs datasets in the following way -- i.i.d. background, and a profile that represents a motif generated from a product multinomial (i.e., PWMs) -- and then plant a realization of that profile at some randomly chosen position for each generated background sequence. But this motif problem is way too easy to tackle. How about we simulate the motif as a mixture of profiles, where each profile may share some identical patterns (i.e., overlaps)? Moreover, what if a motif has a blocked structure such that variable spacings exist between each two adajacent blocks (i.e., gaps)? Maybe let's simulate a mixture of blocked-structured profiles as our ground truth motif? This package creates such patterns.

# Basic examples

Coming soon

# Undetectable patterns 

Coming soon