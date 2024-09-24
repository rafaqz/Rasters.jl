
function DD.show_after(io::IO, mime::MIME"text/plain", A::AbstractRaster; kw...)
    blockwidth = get(io, :blockwidth, 0)
    print_geo(io, mime, A; blockwidth)
    DD.print_block_close(io, blockwidth)
    ndims(A) > 0 && println(io)
    if !(parent(A) isa DiskArrays.AbstractDiskArray)
        DD.print_array(io, mime, A)
    end
end

function DD.show_after(io::IO, mime, stack::AbstractRasterStack; kw...) 
    blockwidth = get(io, :blockwidth, 0)
    print_geo(io, mime, stack; blockwidth)
    DD.print_block_close(io, blockwidth)
end

function print_geo(io, mime, A; blockwidth) 
    DD.print_block_separator(io, "raster", blockwidth)
    printstyled(io, "\n  extent: "; color=:light_black)
    show(io, mime, Extents.extent(A))
    println(io)
    if missingval(A) !== nothing
        printstyled(io, "  missingval: "; color=:light_black)
        show(io, mime, missingval(A))
    end
    if crs(A) !== nothing
        printstyled(io, "\n  crs: "; color=:light_black)
        print(io, convert(String, crs(A)))
    end
    if mappedcrs(A) !== nothing
        printstyled(io, "\n  mappedcrs: "; color=:light_black)
        print(io, convert(String, mappedcrs(A)))
    end
    if isdisk(A)
        fn = filename(A)
        if !(fn == "")
            printstyled(io, "\n  filename: "; color=:light_black)
            print(io, )
        end
    end
    println(io)
end


# Stack types can be enormous. Just use nameof(T)
function Base.summary(io::IO, ser::AbstractRasterSeries{T,N}) where {T,N}
    if N == 0  
        print(io, "0-dimensional ")
    elseif N == 1
        print(io, size(ser, 1), "-element ")
    else
        print(io, join(size(ser), "×"), " ")
    end
    print(io, _series_eltype_string(ser))
    return nothing
end

function DD.show_after(io::IO, mime, ser::AbstractRasterSeries; kw...) 
    if length(ser) > 0
        blockwidth = get(io, :blockwidth, 0)
    end
    if length(ser) > 0
        DD.print_block_separator(io, "holding objects", blockwidth)
        println(io)
        show(io, mime, first(ser))
        println(io)
        println(io)
        DD.print_block_close(io, blockwidth)
        println(io)
        T = _series_eltype_string(ser)
        DD.print_array(io, mime, map(_ -> true, ser))
    else
        DD.print_block_close(io, blockwidth)
    end
end

_series_eltype_string(ser::AbstractRasterSeries{T,N}) where {T,N} =
    string(nameof(typeof(ser)), "{$(nameof(T)),$N}")
