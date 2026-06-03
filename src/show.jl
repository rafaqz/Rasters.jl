
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
    if missingval(A) !== nothing
        printstyled(io, "\n  missingval: "; color=:light_black)
        show(io, mime, missingval(A))
    end
    printstyled(io, "\n  extent: "; color=:light_black)
    str = sprint(show, mime, Extents.extent(A))
    _print_restricted(io, str; blockwidth, used=10)
    if crs(A) !== nothing
        printstyled(io, "\n  crs: "; color=:light_black)
        str = convert(String, crs(A))
        _print_restricted(io, str; blockwidth, used=7)
    end
    if mappedcrs(A) !== nothing
        printstyled(io, "\n  mappedcrs: "; color=:light_black)
        str = convert(String, mappedcrs(A))
        _print_restricted(io, str; blockwidth, used=13)
    end
    if isdisk(A)
        fn = filename(A)
        if !isnothing(fn) && !(fn == "") 
            printstyled(io, "\n  filename: "; color=:light_black)
            _print_restricted(io, fn; blockwidth, used=12)
        end
    end
    println(io)
end

function _print_restricted(io, str; blockwidth, used)
    if length(str) > (blockwidth - used)
        print(io, str[1:min(blockwidth - used - 3, end)] * "...")
    else
        print(io, str)
    end
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
    print(io, string(nameof(typeof(ser)), "{$(nameof(T)),$N}"))
end

function DD.show_after(io::IO, mime, series::AbstractRasterSeries; kw...) 
    if length(series) > 0
        blockwidth = get(io, :blockwidth, 0)
        print_geo(io, mime, first(series); blockwidth)
        DD.print_block_close(io, blockwidth)
    end
    ndims(series) > 0 && println(io)
    DD.print_array(io, mime, series)
end
