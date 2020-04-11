using Kash
using MinHash
using Test
using FASTX
using BioSequences

## FastaIterator
@testset "FastaIterator" begin
    function test_fastaiterator(one::FastaIterator{T}, two) where T
        # FastaIterators mutate their output, so must copy
        A = [copy(i) for i in one]
        B = [FASTA.sequence(LongSequence{T}, i) for i in two]
        @test length(A) == length(B)
        @test all(i == j for (i, j) in zip(A, B))
        close(one.reader)
        close(two)
    end


    path = "/Users/jakni/Documents/scripts/Kash/test/data/test.fna"
    for T in [DNAAlphabet{2}, DNAAlphabet{4}]
        # Instantiate from FASTA Reader
        reader_one = FastaIterator{T}(FASTA.Reader(open(path)))
        reader_two = FASTA.Reader(open(path))
        test_fastaiterator(reader_one, reader_two)

        # Instantiate from IO
        reader_one = FastaIterator{T}(open(path))
        reader_two = FASTA.Reader(open(path))
        test_fastaiterator(reader_one, reader_two)

        # Instantiate from other IO
        io = open(path) do file
            IOBuffer(read(file))
        end
        reader_one = FastaIterator{T}(io)
        reader_two = FASTA.Reader(open(path))
        test_fastaiterator(reader_one, reader_two)

        # Instantiate from IO
        reader_one = FastaIterator{T}(path)
        reader_two = FASTA.Reader(open(path))
        test_fastaiterator(reader_one, reader_two)
    end
end

## CanonicalKmerIterator
@testset "CanonicalKmerIterator" begin
    function test_can_kmerit(it1::CanonicalKmerIterator{M}, it2::M) where M
        A = collect(it1)
        B = [canonical(i) for i in it2]
        @test A == B
    end

    it2 = each(DNAMer{11}, dna"TAGTAGGCGCGCGGCCGTATAGATGATGCTAGAAAGGC")
    it1 = CanonicalKmerIterator(it2)
    test_can_kmerit(it1, it2)

    it2 = each(DNAMer{3}, dna"TAGTAGNGCGCGCGGCNCGTATAGATGATGCTAGAAAGNGC")
    it1 = CanonicalKmerIterator(it2)
    test_can_kmerit(it1, it2)

    it2 = each(DNAMer{4}, dna"TAGTCGAAGTGCTGAGAGATCTCTAAAGAGAGCGCTCTGAAAAA", 3, 5)
    it1 = CanonicalKmerIterator(it2)
    test_can_kmerit(it1, it2)
end

## Kmersketcher
@testset "KmerSketcher" begin
    @testset "Instantiation" begin
        h = MinHasher{identity}(100)
        @test isa(KmerSketcher{DNAMer{4}}(h), KmerSketcher{DNAMer{4}, identity})

        @test_throws MethodError KmerSketcher{DNAMer{3},hash}(h)
        @test_throws MethodError KmerSketcher{4,identity}(h)

        sk = KmerSketcher{DNAMer{5}}(100)
        @test isa(sk.hasher, MinHasher{hash})
    end

    @testset "empty!" begin
        sk = KmerSketcher{DNAMer{3}}(10)
        seq = dna"TAGGCGTAGTGCGTATATAGCGAAAGAGCTCTA"
        Kash.update!(sk, seq)

        @test sk.bases == length(seq)
        @test sk.hasher.filled ≤ length(seq) - 2
        empty!(sk)
        @test sk.bases == 0
        @test sk.hasher.filled == 0
    end

    @testset "update! seq" begin
        function test_update(sk::KmerSketcher{M,F}, seq::BioSequence) where {M,F}
            hashes_before = KmerHashes(sk)
            Kash.update!(sk, seq)
            hashes_after = KmerHashes(sk)

            # Test bases is correct
            @test hashes_before.bases + length(seq) == hashes_after.bases

            # Test content of hashes
            new_hashes = [F(canonical(i)) for i in each(M, seq)]
            should_be = sort!(collect(Set(append!(hashes_before.sketch.hashes, new_hashes))))

            # at most max hashes
            new_length = min(length(should_be), length(sk.hasher.heap))

            should_be = should_be[1:new_length]
            @test hashes_after.sketch.hashes == should_be
        end

        sk = KmerSketcher{DNAMer{3}}(10)
        seq = dna"TAGGCGTAGTNGCGTATATAGCGAAAGAGCTCNTA"
        test_update(sk, seq)

        empty!(sk)
        Kash.update!(sk, randdnaseq(10))
        test_update(sk, randdnaseq(100))
    end
end

## Kmerhashes

## High-level API


@testset "Kash.jl" begin
    # Write your own tests here.
end
