mutable struct FastaIterator{A}
    reader::FASTA.Reader # UnionAll type, but whatever
    record::FASTA.Record
    seq::LongSequence{A}

    function FastaIterator{A}(reader, record, seq) where {A <: BioSequences.Alphabet, R}
        new(reader, record, seq)
    end
end

function FastaIterator{A}(reader::FASTA.Reader) where {A <: BioSequences.Alphabet}
    return FastaIterator{A}(reader, FASTA.Record(), LongSequence{A}())
end

function FastaIterator{A}(io::IO) where {A <: BioSequences.Alphabet}
    return FastaIterator{A}(FASTA.Reader(io))
end

function FastaIterator{A}(path::AbstractString) where {A <: BioSequences.Alphabet}
    return FastaIterator{A}(FASTA.Reader(open(path)))
end

Base.eltype(::Type{FastaIterator{A}}) where A = LongSequence{A}
Base.IteratorEltype(::Type{FastaIterator{A}}) where A = Base.HasEltype()
Base.IteratorSize(::Type{<:FastaIterator}) = Base.SizeUnknown()

function Base.iterate(it::FastaIterator, ::Nothing=nothing)
    if eof(it.reader)
        close(it.reader)
        return nothing
    end
    read!(it.reader, it.record)
    copy!(it.seq, it.record)
    return it.seq, nothing
end

struct CanonicalKmerIterator{M}
    it::M
end

function CanonicalKmerIterator(it::M) where {M <: BioSequences.AbstractMerIterator}
    return CanonicalKmerIterator{M}(it)
end

Base.IteratorEltype(::Type{CanonicalKmerIterator{M}}) where M = Base.HasEltype()
Base.IteratorSize(::Type{CanonicalKmerIterator{M}}) where M = Base.IteratorSize(M)
Base.length(x::CanonicalKmerIterator{M}) where M = length(x.it)

@inline function Base.iterate(it::CanonicalKmerIterator, s...)
    itval = iterate(it.it, s...)
    itval === nothing && return nothing
    itval = itval::Tuple{Any, Any} # hack from Base.generator to help inference
    return canonical(itval[1]), itval[2]
end

mutable struct KmerSketcher{M,F}
    hasher::MinHasher{F}
    bases::Int

    function KmerSketcher{M,F}(h::MinHasher{F}) where {M <: Mer, F}
        new{M,F}(h, 0)
    end
end

KmerSketcher{M}(h::MinHasher{F}) where {M,F} = KmerSketcher{M,F}(h)
function KmerSketcher{M,F}(s::Integer) where {M <: Mer, F}
    KmerSketcher{M,F}(MinHasher(s))
end

KmerSketcher{M}(s::Integer) where {M <: Mer} = KmerSketcher{M,hash}(s)

function Base.empty!(sk::KmerSketcher)
    empty!(sk.hasher)
    sk.bases = 0
    return sk
end

function update!(sketcher::KmerSketcher{M}, seq::BioSequence) where M
    it = CanonicalKmerIterator(each(M, seq))
    MinHash.update!(sketcher.hasher, it)
    sketcher.bases += length(seq)
    return sketcher
end

function update!(sketcher::KmerSketcher, io::IO)
    it = FastaIterator{DNAAlphabet{4}}(io)
    for seq in it
        update!(sketcher, seq)
    end
    return sketcher
end

# Parameterize to allow different hash functions. Rarely used, so there are few convenience
# functions for this feature.
struct KmerHashes{K,F}
    bases::Int
    sketch::MinHash.MinHashSketch
end

# This also exists in BioSequences, but is an internal method
ksize(::Type{<:AbstractMer{A,K}}) where {A,K} = K

function KmerHashes(sk::KmerSketcher{M,F}) where {M,F}
    return KmerHashes{ksize(M),F}(sk.bases, MinHash.MinHashSketch(sk.hasher))
end

function kmer_minhash(seq::Union{IO, BioSequence}, s::Integer, ::Val{K}) where K
    sketcher = KmerSketcher{DNAMer{K}}(s)
    update!(sketcher, seq)
    return KmerHashes(sketcher)
end

function kmer_minhash_each(it, s::Integer, ::Val{K}) where K
    sketcher = KmerSketcher{DNAMer{K}}(s)
    result = KmerHashes{K,hash}[]
    for i in it
        update!(sketcher, i)
        push!(result, KmerHashes(sketcher))
        empty!(sketcher)
    end
    return result
end

function kmer_minhash_each(io::IO, s::Integer, v::Val)
    it = FastaIterator{DNAAlphabet{4}}(io)
    return kmer_minhash_each(it, s, v)
end

#=
function _add_minhash(result::Vector, path::String, i::Integer, s::Integer, v::Val)
    result[i] = open(path) do file
        kmer_minhash(file, s, v)
    end
end

function kmer_minhash_paths(paths::Vector{String}, s::Integer, v::Val{K}) where K
    result = Vector{KmerHashes{K,hash}}(undef, length(paths))
    processes = Vector{Task}(undef, length(paths))
    for i in eachindex(paths)
        processes[i] = Threads.@spawn _add_minhash(result, paths[i], i, s, v)
    end
    foreach(wait, processes)
    return result
end


function _kmer_minhash(path::String, s::Integer, v::Val)
    return open(path) do file
        kmer_minhash(file, s, v)
    end
end
=#

function kmer_minhash_paths(paths::Vector{String}, s::Integer, v::Val{K}) where K
    processes = Task[]
    for path in paths
        push!(processes, Threads.@spawn open(x -> kmer_minhash(x, s, v), path))
    end
    result = KmerHashes{K,hash}[]
    for process in processes
        push!(result, fetch(process))
    end
    return result
end

function intersectionlength(x::T, y::T) where {T<:KmerHashes}
    return MinHash.intersectionlength(x.sketch, y.sketch)
end

function intersectionlength(v::Vector{KmerHashes{K,F}}) where {K,F}
    v2 = [i.sketch for i in v]
    return MinHash.intersectionlength(v2)
end

function mash_length(h1::KmerHashes{K,F}, h2::KmerHashes{K,F}) where {K,F}
    w = intersectionlength(h1, h2)
    j = w / (length(h1.sketch) + length(h2.sketch) - w)
    return iszero(w) ? 1.0 : -log(2j / (j+1)) / k
end
