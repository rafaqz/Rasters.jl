
# Add method to avoid printing from disk
# function Base.show(io::IO, mime::MIME"text/plain", A::DiskGeoArray{T,N}) where {T,N}
#     lines = _print_array_info(io, mime, A)
#     if !(metadata(A) isa NoMetadata) 
#         print(io, "\nwith ")
#         show(io, mime, metadata(A))
#     end
#     println(io)
#     print(io, "\n$(filename(A))")
# end

function _show_dimname(io, dim::Dim)
    color = DD._dimcolor(io)
    printstyled(io, "Dim{"; color=color)
    printstyled(io, string(":", name(dim)); color=:yellow)
    printstyled(io, "}"; color=color)
end
function _show_dimname(io, dim::Dimension)
    printstyled(io, DD.dim2key(dim); color = DD._dimcolor(io))
end

function Base.show(io::IO, mime::MIME"text/plain", stack::AbstractGeoStack)
    nlayers = length(keys(stack))
    layers_str = nlayers == 1 ? "layer" : "layers"
    printstyled(io, nameof(typeof(stack)), color=:blue)
    # print(io, " with $nlayers $(childtype(stack)) $layers_str:\n")
    for var in keys(stack)
        printstyled(io, "  :$var", color=:yellow)

        field_dims = DD.layerdims(stack, var)
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

    if data(stack) isa FileStack 
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
            show(IOContext(io, :compact=>true), mime, dim)
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
function Base.summary(io::IO, A::AbstractGeoSeries{T,N}) where {T,N}
    if N == 0  
        print(io, "0-dimensional ")
    elseif N == 1
        print(io, size(A, 1), "-element ")
    else
        print(io, join(size(A), "×"), " ")
    end
    printstyled(io, string(nameof(typeof(A)), "{$(nameof(T)),$N}"); color=:blue)
end

function Base.show(io::IO, mime::MIME"text/plain", A::AbstractGeoSeries{T,N}) where {T,N}
    lines = _print_array_info(io, mime, A)
    ds = displaysize(io)
    ioctx = IOContext(io, :displaysize => (ds[1] - lines, ds[2]))
    if T <: AbstractString
        DD._show_array(ioctx, mime, parent(A))
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

function _print_array_info(io, mime, A)
    lines = 0
    summary(io, A)
    DD._printname(io, name(A))
    lines += DD._printdims(io, mime, dims(A))
    !(isempty(dims(A)) || isempty(refdims(A))) && println(io)
    lines += DD._printrefdims(io, mime, refdims(A))
    println(io)
    return lines
end
