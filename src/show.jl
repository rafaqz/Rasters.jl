
function DD.show_after(io::IO, mime::MIME"text/plain", A::AbstractRaster)

    printstyled(io, "extent: "; color=:light_black)
    show(io, mime, Extents.extent(A))
    println()
    if missingval(A) !== nothing
        printstyled(io, "missingval: "; color=:light_black)
        show(io, mime, missingval(A))
        println()
    end
    if crs(A) !== nothing
        printstyled(io, "crs: "; color=:light_black)
        print(io, convert(String, crs(A)), "\n")
    end
    if mappedcrs(A) !== nothing
        printstyled(io, "mappedcrs: "; color=:light_black)
        print(io, convert(String, mappedcrs(A)), "\n")
    end
    if parent(A) isa DiskArrays.AbstractDiskArray 
        if parent(A) isa FileArray 
            printstyled(io, "from file:\n"; color=:light_black)
            print(io, filename(parent(A)))
        end
    else
        printstyled(io, "parent:\n"; color=:light_black)
        DD.print_array(io, mime, A)
    end
end

function DD.show_after(io, mime, stack::AbstractRasterStack) 
    if parent(stack) isa FileStack 
        printstyled(io, "from file:\n"; color=:light_black)
        println(io, filename(stack))
    end
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

function LA.show_properties(io::IO, lookup::AbstractProjected)
    print(io, " ")
    LA.print_order(io, lookup)
    print(io, " ")
    LA.print_span(io, lookup)
    print(io, " ")
    LA.print_sampling(io, lookup)
    if !(crs(lookup) isa Nothing)
        print(io, " crs: ", nameof(typeof(crs(lookup))))
    end
    if !(mappedcrs(lookup) isa Nothing)
        print(io, " mappedcrs: ", nameof(typeof(mappedcrs(lookup))))
    end
end
