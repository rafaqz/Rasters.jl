const ASCII_DEFAULT_CRS = EPSG(4326)

has_layers(::Type{ASCIIfile}) = false

struct ASCIIparams{T, F}
    filename::F
    params::NamedTuple
    write::Bool
end

function ASCIIparams(filename::AbstractString; write=false)
    pars = ASCIIrasters.read_ascii(filename; lazy = true)
    T = typeof(pars[:nodatavalue])
    ASCIIparams{T, typeof(filename)}(filename, pars, write)
end

params(ap::ASCIIparams) = ap.params
filename(ap::ASCIIparams) = ap.filename

# DimensionalData methods
#########################

function DD.dims(ap::ASCIIparams, crs=nothing, mappedcrs=nothing)
    crs = crs isa Nothing ? ASCII_DEFAULT_CRS : crs

    nc, nr = size(ap)

    pars = params(ap)

    xbounds = (pars[:xll], pars[:xll] + pars[:dx] * nc)
    ybounds = (pars[:yll], pars[:yll] + pars[:dy] * nr)

    # Always intervals
    xspan = (xbounds[2] - xbounds[1]) / nc
    yspan = (ybounds[2] - ybounds[1]) / nr

    # Not implemented yet
    xy_metadata = Metadata{ASCIIfile}(Dict())

    xindex = LinRange(xbounds[1], xbounds[2] - xspan, nc)
    yindex = LinRange(ybounds[2] - yspan, ybounds[1], nr)

    xlookup = Projected(xindex;
        order=ForwardOrdered(),
        span=Regular(xspan),
        sampling=Intervals(Start()),
        metadata=xy_metadata,
        crs=crs,
        mappedcrs=mappedcrs,
        dim=X()
    )
    ylookup = Projected(yindex;
        order= ReverseOrdered(),
        span=Regular(yspan),
        sampling=Intervals(Start()),
        metadata=xy_metadata,
        crs=crs,
        mappedcrs=mappedcrs,
        dim=Y()
    )

    x = X(xlookup)
    y = Y(ylookup)

    return x,y

end

missingval(ap::ASCIIparams) = params(ap)[:nodatavalue]

# no metadata for now
DD.metadata(ap::ASCIIparams) = Metadata{ASCIIfile}()
DD.ndims(ap::ASCIIparams) = 2

# Base util methods

Base.eltype(p::ASCIIparams{T}) where T = T

function Base.size(ap::ASCIIparams)
    params(ap)[:ncols], params(ap)[:nrows]
end

Base.Array(ap::ASCIIparams) = _asciigrid(a -> Array(a), ap)

# Array
#######

function FileArray(ap::ASCIIparams, filename = filename(ap); kw...)
    size_ = size(ap)
    eachchunk = DiskArrays.GridChunks(size_, size_)
    haschunks = DiskArrays.Unchunked()
    T = eltype(ap)
    N = ndims(ap)
    FileArray{ASCIIfile, T, N}(filename, size_; eachchunk, haschunks, kw...)
end

# Base i/o methods
###################

# data (ASCIIfile) and metadata (ASCIIparams) objects are separate 
# so data is opened using _asciigrid
function Base.open(f::Function, A::FileArray{ASCIIfile}, key...; write = A.write)
    _asciigrid(dat -> f(RasterDiskArray{ASCIIfile}(dat, A.eachchunk, A.haschunks)), A; write)
end

# also called by _open(f, fa::FileArray{ASCIIfile})
function _open(f, ::Type{ASCIIfile}, filename; write = false, kw...)
    isfile(filename) || _filenotfound_error(filename)
    _open(f, ASCIIfile, ASCIIparams(filename; write))
end

_open(f, ::Type{ASCIIfile}, ap::ASCIIparams; kw...) = f(ap)

function _asciigrid(f::Function, ap::Union{ASCIIparams, FileArray}; kw...)
    _asciigrid(f, filename(ap), eltype(ap), size(ap); kw...)
end

function _asciigrid(f, filename::AbstractString, T::Type, size::Tuple; write = false, kw...)
    dat, pars = open(filename, "r") do io
        dat, pars = ASCIIrasters.read_ascii(filename; lazy = false)
        # dat is a nr x nc matrix, we want a nc x nr matrix for use
        # as Raster.data
        mat = _flip(dat, size, _detect_datatype(pars))
        output = f(mat)
        output, pars
    end
    if write
        mat = _flip(dat, (size[2], size[1]), _detect_datatype(pars))
        ASCIIrasters.write_ascii(filename, mat; pars...)
    end
    dat
end

# Utilities
###########

function _detect_datatype(pars)
    Float64
end

function _flip(mat, size, type)
    new_nr, new_nc = size # size is a raster size so nc x nr
    out = Matrix{type}(undef,new_nr, new_nc)
    for r in 1:new_nr
        for c in 1:new_nc
            out[r,c] = mat[c,r]
        end
    end
    out
end