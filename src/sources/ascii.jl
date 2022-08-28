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

# Base methods