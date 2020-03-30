mutable struct FastaIterator{A}
    reader::FASTA.Reader
    record::FASTA.Record
    seq::LongSequence{A}

    function FastaIterator{A}(reader, record, seq) where {A <: BioSequences.Alphabet}
        new(reader, record, seq)
    end
end

function FastaIterator{A}(reader::FASTA.Reader) where {A <: BioSequences.Alphabet}
    return FastaIterator{A}(reader, FASTA.Record(), LongSequence{A}())
end

function FastaIterator{A}(io::IO) where {A <: BioSequences.Alphabet}
    return FastaIterator{A}(FASTA.Reader(io), FASTA.Record(), LongSequence{A}())
end

function FastaIterator{A}(path::AbstractString) where {A <: BioSequences.Alphabet}
    return FastaIterator{A}(FASTA.Reader(open(path)), FASTA.Record(), LongSequence{A}())
end

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

function Base.iterate(it::CanonicalKmerIterator, state...)
    itval = iterate(it.it, state...)
    itval == nothing && return nothing
    val, state = itval
    return canonical(val), state
end


mutable struct KmerSketcher{M,F}
    hasher::MinHasher{F}
    bases::Int

    function KmerSketcher{M,F}(h::MinHasher) where {M <: Mer, F}
        new{M,F}(h, 0)
    end
end

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
    it = FastaIterator{DNAAlphabet{2}}(io)
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
    it = FastaIterator{DNAAlphabet{2}}(io)
    return kmer_minhash_each(it, s, v)
end

function intersectionlength(x::T, y::T) where {T<:KmerHashes}
    return MinHash.intersectionlength(x.sketch, y.sketch)
end

function intersectionlength(v::Vector{KmerHashes{K,F}}) where {K,F}
    v2 = [i.sketch for i in v]
    return MinHash.intersectionlength(v2)
end
#= TODO:

=#

#= How would I want the interface to be?

Sketch all sequences in a file:
    kmer_minhash(io::IO, s, ::Val{K})

Sketch a sequence
    kmer_minhash(seq::BioSequence, s, ::Val{K})

Maybe something to get all sketches of a file re-using the same sketcher?
    kmer_minhash(io::IO, s, ::Val{K})

=#
