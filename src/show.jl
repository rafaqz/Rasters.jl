
function DD.show_after(io::IO, mime::MIME"text/plain", A::AbstractRaster; kw...)
    blockwidth = get(io, :blockwidth, 0)
    print_geo(io, mime, A; blockwidth)
    DD.print_block_close(io, blockwidth)
    ndims(A) > 0 && println(io)
    if parent(A) isa DiskArrays.AbstractDiskArray 
        if parent(A) isa FileArray 
            printstyled(io, "from file:\n"; color=:light_black)
            print(io, filename(parent(A)))
        else
            show(io, mime, parent(A))
        end
    else
        DD.print_array(io, mime, A)
    end
end

function DD.show_after(io, mime, stack::AbstractRasterStack; kw...) 
    blockwidth = get(io, :blockwidth, 0)
    print_geo(io, mime, stack; blockwidth)
    if parent(stack) isa FileStack 
        printstyled(io, "from file:\n"; color=:light_black)
        println(io, filename(stack))
    end
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
    println(io)
end


# Stack types can be enourmous. Just use nameof(T)
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

function Base.show(io::IO, mime::MIME"text/plain", A::AbstractRasterSeries{T,N}) where {T,N}
    lines = 0
    summary(io, A)
    DD.print_name(io, name(A))
    lines += Dimensions.print_dims(io, mime, dims(A))
    !(isempty(dims(A)) || isempty(refdims(A))) && println(io)
    lines += Dimensions.print_refdims(io, mime, refdims(A))
    println(io)
    ds = displaysize(io)
    ioctx = IOContext(io, :displaysize => (ds[1] - lines, ds[2]))
end
