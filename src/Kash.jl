module Kash

using FASTX
using BioSequences
using BSON
using MinHash

# Need to use CodecZlib for Gzip integration!

include("kmersketch.jl")
include("serialize.jl")

export FastaIterator,
    CanonicalKmerIterator,
    KmerSketcher,
    KmerHashes,
    update!,
    kmer_minhash,
    kmer_minhash_each,
    intersectionlength

end # module
