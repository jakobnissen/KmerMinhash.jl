const DEFAULT_VERSION = 0

"Converts between bitstypes, and little endian to host."
function convert_array(T::Type{<:Base.BitInteger64}, src::Vector{<:Base.BitInteger64})
    nbytes = length(src) * sizeof(eltype(src))
    dst = Vector{T}(undef, nbytes ÷ sizeof(T))
    GC.@preserve dst src begin
        srcptr = Ptr{UInt}(pointer(src))
        dstptr = Ptr{UInt}(pointer(dst))
        for i in 1:nbytes ÷ sizeof(UInt)
            unsafe_store!(dstptr, htol(unsafe_load(srcptr, i)), i)
        end
    end
    return dst
end

convert_hashfunction(sk::KmerHashes{K}) where K = KmerHashes{K,hash}(sk.bases, sk.sketch)
convert_hashfunction(sk::Vector{<:KmerHashes}) = [convert_hashfunction(s) for s in sk]

@noinline function validate_hashfunction(F)
    F === hash || error("""Kash currently only supports serializing with hashfunction
    Base.hash, see `convert_hashfunction`""")
end

function bson_dict(v::Val{0}, x::KmerHashes)
    d = Dict{Symbol,Any}()
    d[:r] = Int64(x.sketch.requested)
    d[:b] = Int64(x.bases)
    d[:h] = convert_array(UInt8, x.sketch.hashes)
    return d
end

function bson_dict(v::Val{0}, sketches::Vector{KmerHashes{K,F}}) where {K,F}
    validate_hashfunction(F)
    d = Dict{Symbol,Any}()
    d[:m] = Int32(1752392011) # UTF-8 for "Kash", magic header
    d[:v] = Int32(0)
    d[:k] = Int32(K)

    array = Vector{Any}() # must be a Vector{Any} for BSON
    for sketch in sketches
        push!(array, bson_dict(Val(0), sketch))
    end
    d[:h] = array
    return d
end

# This function should only verify magic header and check version.
# Then it dispatches to the version-specific parser
function load_sketches_bson(dict::Dict{Symbol,Any})
    if !haskey(dict, :m) || (dict[:m] != Int32(1752392011))
        throw(ValueError("Magic number 1752392011 not found in BSON"))
    end
    v = Int(dict[:v])
    if !(v in (0,))
        throw(ValueError("This version of Kash cannot read Kash BSON file v. $v"))
    end
    return load_sketches_bson(Val{v}(), dict)
end

function load_sketches_bson(v::Val{0}, dict::Dict{Symbol,Any})
    K = Int(dict[:k])
    T = KmerHashes{K,hash}
    array = dict[:h]
    sketches = Vector{T}(undef, length(array))
    for (i, sketchdict) in enumerate(array)
        sketch = load_sketch_bson(Val(0), T, sketchdict)
        sketches[i] = sketch
    end
    return sketches
end

function load_sketch_bson(v::Val{0}, T::Type{<:KmerHashes}, dict::Dict{Symbol,Any})
    r = Int(dict[:r])
    b = Int(dict[:b])
    h = convert_array(UInt, dict[:h])
    sketch = MinHash.MinHashSketch(hash, r, h)
    return T(b, sketch)
end

function load_sketches(io::IO)
    d = BSON.parse(io)
    load_sketches_bson(d)
end

load_sketches(path::String) = open(load_sketches, path)

function save_sketches(io::IO, sketches::Vector{KmerHashes{K,F}}) where {K,F}
    d = bson_dict(Val(DEFAULT_VERSION), sketches)
    bson(io, d)
end

function save_sketches(path::String, sketches::Vector{KmerHashes{K,F}}) where {K,F}
    open(io -> save_sketches(io, sketches), path, "w")
end
