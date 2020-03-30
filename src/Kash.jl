module Kash

using FASTX
using BioSequences
using BSON
using MinHash

# Need to use CodecZlib for Gzip integration!

include("kmersketch.jl")
include("serialize.jl")

#= Keep the commented while developing for each of testing
export Sketcher,
update!,
KmerSketch,
sketch,
overlaps, # one-to-one overlaps of hashes
overlap_matrix,
load_sketches,
save_sketches,
=#

a = randdnaseq(1000)
b = randdnaseq(1000)
c = kmer_minhash(a, 10, Val(10))
d = kmer_minhash(b, 10, Val(10))
e = [c, d]

end # module
