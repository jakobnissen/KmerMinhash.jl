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

end # module
