# Kash.jl
_Kmer minhashing in Julia_

[![Build Status](https://travis-ci.com/jakobnissen/Kash.jl.svg?branch=master)](https://travis-ci.com/jakobnissen/Kash.jl)
[![Codecov](https://codecov.io/gh/jakobnissen/Kash.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jakobnissen/Kash.jl)

__This package is work in progress. I'll work on it when I have time.__

The purpose of this package is to:

* Provide a simple Julia interface
* Provide methods for minhashing DNA sequences and FASTA files
* Provide methods for distance calculation between minhash sketches
* Be highly efficient

The purpose of this package is NOT:

* To hash amino acids or arbitrary alphabets
* To hash reads. This functionality may be added later.

## Usage
_Note: The API is not stable. Version is not yet 1.0_

To be determined.

## Benchmarks
__Case 1: Sketching 1 core, 51k sequences, 330 Mbp, 1 sketch__

```
Hashes     1e3     1e4
Mash 2.2: 15.4s   15.2s
Kash    :  2.6s    2.6s
```

__Case 2: Sketching 1 core, 51k sequences, 330 Mbp, 51k sketches__
```
Hashes     1e2     1e3
Mash 2.2: 20.4s   53.5s
Kash    :  4.2s   12.5s
```

__Case 3: Distance matrix, 1k sketches, 1k hashes__
```
Mash 2.2: 17.87s
Kash    :   .17s
```
