
# Add method to avoid printing from disk
function Base.show(io::IO, mime::MIME"text/plain", A::DiskGeoArray{T,N}) where {T,N}
    printstyled(io, string(nameof(typeof(A)), "{$T,$N}"); color=:blue)
    DD._printname(io, name(A))
    DD._printdims(io, mime, dims(A))
    DD._printrefdims(io, mime, refdims(A))
    print(io, "\nFrom file: $(filename(A))")
    println(io)
    if !(metadata(A) isa NoMetadata) 
        print(io, "\nwith ")
        show(io, mime, metadata(A))
    end
end

function Base.show(io::IO, mime::MIME"text/plain", stack::AbstractGeoStack)
    nlayers = length(keys(stack))
    layers_str = nlayers == 1 ? "layer" : "layers"
    printstyled(io, nameof(typeof(stack)), color=:blue)
    if stack isa DiskGeoStack 
        print(io, " with $nlayers $(childtype(stack)) $layers_str:")
        if !(filename(stack) isa String)
            for (key, fn) in pairs(filename(stack))
                print(io, "\n ")
                printstyled(io, " $key", color=:green)
                print(io, " = ", fn)
            end
        end
    else
        print(io, " with $nlayers $layers_str")
    end
    print(io, '\n')

    for var in keys(stack)
        printstyled(io, "  $var", color=:green)

        field_dims = dims(stack, var)
        n_dims = length(field_dims)
        dims_str = n_dims == 1 ? "dimension" : "dimensions"
        print(io, " with $n_dims $dims_str: ")
        if n_dims > 0
            for (d, dim) in enumerate(field_dims)
                printstyled(io, "$(name(dim))", color=:red)
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

function Base.show(io::IO, mime::MIME"text/plain", mode::AbstractProjected)
    DD._printmode(io, mode)
    print(io, " - ")
    DD._printorder(io, mode)
    print(io, " ")
    DD._printspan(io, mode)
    print(io, " ")
    DD._printsampling(io, mode)
    print(io, " crs: ", nameof(typeof(crs(mode))))
    print(io, " mappedcrs: ", nameof(typeof(mappedcrs(mode))))
end
