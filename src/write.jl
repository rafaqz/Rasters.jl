"""
    Base.write(filename::AbstractString, A::AbstractRaster; kw...)

Write an [`AbstractRaster`](@ref) to file, guessing the backend from the file extension.

Keyword arguments are passed to the `write` method for the backend.
"""
function Base.write(
    filename::AbstractString, A::AbstractRaster; kw...
)
    write(filename, _sourcetype(filename), A; kw...)
end
Base.write(A::AbstractRaster) = write(filename(A), A)

"""
    Base.write(filename::AbstractString, s::AbstractRasterStack; suffix, kw...)

Write any [`AbstractRasterStack`](@ref) to file, guessing the backend from the file extension.

## Keywords

- `suffix`: suffix to append to file names. By default the layer key is used. 

Other keyword arguments are passed to the `write` method for the backend.

If the source can't be saved as a stack-like object, individual array layers will be saved.
"""
function Base.write(filename::AbstractString, s::AbstractRasterStack; suffix=nothing, kw...)
    base, ext = splitext(filename)
    T = _sourcetype(filename)
    if haslayers(T)
        write(filename, _sourcetype(filename), s; kw...)
    else
        # Otherwise write separate files for each layer
        if suffix === nothing
            suffix = map(k -> string("_", k), keys(s))
        end
        @warn string("Cannot write stacks to \"", ext, "\", writing layers as individual files")
        map(keys(s), suffix) do key, sfx
            fn = string(base, sfx, ext)
            write(fn, _sourcetype(filename), s[key])
        end
    end
end

# Trait for source data that has stack layers
haslayers(T) = false
