function _show_dimname(io, dim::Dim)
    color = _dimcolor(io)
    printstyled(io, "Dim{"; color=color)
    printstyled(io, string(":", name(dim)); color=:yellow)
    printstyled(io, "}"; color=color)
end
function _show_dimname(io, dim::Dimension)
    printstyled(io, DD.dim2key(dim); color = _dimcolor(io))
end

function Base.show(io::IO, mime::MIME"text/plain", stack::AbstractGeoStack)
    nlayers = length(keys(stack))
    layers_str = nlayers == 1 ? "layer" : "layers"
    printstyled(io, nameof(typeof(stack)), color=:blue)
    print(io, " with $nlayers $(childtype(stack)) $layers_str:\n")
    for var in keys(stack)
        printstyled(io, "  :$var", color=:yellow)

        field_dims = dims(stack, var)
        n_dims = length(field_dims)
        dims_str = n_dims == 1 ? "dim" : "dims"
        print(io, " with $dims_str: ")
        if n_dims > 0
            for (d, dim) in enumerate(field_dims)
                _show_dimname(io, dim)
                d != length(field_dims) && print(io, ", ")
            end
            print(io, " (")
            for (d, dim) in enumerate(field_dims)
                print(io, "$(length(dim))")
                d != length(field_dims) && print(io, '×')
            end
            print(io, ')')
        end
        print(io, '\n')
    end

    if stack isa DiskGeoStack 
        if filename(stack) isa AbstractString
            println(io, "\n" * filename(stack))
        else
            for (key, fn) in pairs(filename(stack))
                print(io, "\n ")
                printstyled(io, " :$key", color=:yellow)
                print(io, " = ", fn)
            end
        end
    end

    n_windows = length(window(stack))
    if n_windows > 0
        print(io, "with window:\n")
        for dim in window(stack)
            print(io, ' ')
            show(IOContext(io; :compact=>true), mime, dim)
            print(io, '\n')
        end
    end

    md = metadata(stack)
    if !(md isa NoMetadata)
        n_metadata = length(md)
        if n_metadata > 0
            print(io, "\nwith ")
            show(io, mime, md)
        end
    end
end

# Stack types can be enourmous. Just use nameof(T)
function Base.summary(io::IO, ser::AbstractGeoSeries{T,N}) where {T,N}
    if N == 0  
        print(io, "0-dimensional ")
    elseif N == 1
        print(io, size(ser, 1), "-element ")
    else
        print(io, join(size(ser), "×"), " ")
    end
    printstyled(io, string(nameof(typeof(ser)), "{$(nameof(T)),$N}"); color=:blue)
    return nothing
end

function Base.show(io::IO, mime::MIME"text/plain", A::AbstractGeoSeries{T,N}) where {T,N}
    summary(io, A)
    # TODO: these need a proper interface in DD
    _printname(io, name(A))
    _printdims(io, mime, dims(A))
    !(isempty(dims(A)) || isempty(refdims(A))) && println(io)
    _printrefdims(io, mime, refdims(A))
    return nothing
end

# Move to proper interface in DD
function _printname(io::IO, name)
    if !(name == Symbol("") || name isa NoName)
        printstyled(io, string(" :", name); color=:yellow)
    end
end
function _printdims(io::IO, mime, dims::Tuple)
    if isempty(dims) 
        print(io, ": ")
        return 0
    end
    printstyled(io, " with dimensions: "; color=:light_black)
    return _layout_dims(io, mime, dims)
end
function _printrefdims(io::IO, mime, refdims::Tuple)
    if isempty(refdims) 
        return 0
    end
    printstyled(io, "and reference dimensions: "; color=:light_black)
    ctx = IOContext(io, :is_ref_dim=>true, :show_dim_val=>true)
    lines = _layout_dims(ctx, mime, refdims)
    return lines
end
function _layout_dims(io, mime, dims::Tuple)
    length(dims) > 0 || return 0
    ctx = IOContext(io, :compact=>true)
    if all(m -> m isa NoIndex, mode(dims))
        for d in dims[1:end-1]
            show(ctx, mime, d)
            print(io, ", ")
        end
        show(ctx, mime, dims[end])
        return 0
    else # Dims get a line each
        lines = 3
        print(io, "\n  ")
        for d in dims[1:end-1]
            show(ctx, mime, d)
            print(io, ",")
            lines += 2 # Often they wrap
            print(io, "\n  ")
        end
        show(ctx, mime, dims[end])
        return lines
    end
end

_dimcolor(io) = get(io, :is_ref_dim, false) ? :magenta : :red

function Base.show(io::IO, mime::MIME"text/plain", mode::AbstractProjected)
    DD.show(io, mime, Sampled(mode.order, mode.span, mode.sampling))
    if !(crs(mode) isa Nothing)
        print(io, " crs: ", nameof(typeof(crs(mode))))
    end
    if !(mappedcrs(mode) isa Nothing)
        print(io, " mappedcrs: ", nameof(typeof(mappedcrs(mode))))
    end
    return nothing
end
