
export ASCIIarray, ASCIImetadata, ASCIIdimMetadata

# Metadata ########################################################################
struct ASCIImetadata{K,V} <: ArrayMetadata{K,V}
    val::Dict{K,V}
end

struct ASCIIdimMetadata{K,V} <: DimMetadata{K,V}
    val::Dict{K,V}
end


# Array ########################################################################

struct ASCIIarray{T,N,A,D<:Tuple,R<:Tuple,Me,Mi,Na,W,SS} <: DiskGeoArray{T,N,D,LazyArray{T,N}}
    filename::A
    dims::D
    refdims::R
    metadata::Me
    missingval::Mi
    name::Na
    window::W
    size::S
    source_size::SS
end

ASCIIarray(x::AbstractString; refdims=(), name="Unnamed", window=()) =
    open(x, "r") do file
        ncols = parse(Int, match(r"NCOLS (.+)", readline(file)).captures[1])
        nrows = parse(Int, match(r"NROWS (.+)", readline(file)).captures[1])
        xllstring = readline(file)
        yllstring = readline(file)
        step = parse(Float64, match(r"CELLSIZE (.+)", readline(file)).captures[1])
        nodata = parse(Float64, match(r"NODATA_value (.+)", readline(file)).captures[1])

        xllmatch = match(r"XLLCORNER (.+)", xllstring).captures[1]
        if xllmatch === nothing 
            xllmatch = match(r"XLLCENTER (.+)", xllstring)
            xlocus = Center()
        else
            xlocus = Start()
        end
        xll = parse(Float64, xllmatch.captures[1])

        yllmatch = match(r"YLLCORNER (.+)", yllstring).captures[1]
        if yllmatch === nothing 
            yllmatch = match(r"YLLCENTER (.+)", yllstring)
            ylocus = Center()
        else
            ylocus = Start()
        end
        yll = parse(Float64, yllmatch.captures[1])

        dims = Lon(xll; mode=RegularIndex(locus=xlocus, step=step)), 
               Lat(yll; mode=RegularIndex(locus=ylocus, step=step))
        size = nrows, ncols
        if window != ()
            window = to_indices(dataset, dims2indices(dims, window))
            println(window)
            size = windowsize(window)
        end
        T = AG.getdatatype(AG.getband(dataset, 1))
        N = length(size)

        ASCIIarray(; filename=filename, dims=dims, refdims=refdims,
                   metadata=metadata, missingval=nodata, window=window, size=size)
    end

filename(A::ASCIIarray) = A.filename
crs(A::ASCIIarray) = nothing

Base.size(A::ASCIIarray) = A.size

Base.parent(A::ASCIIarray) =
    asciiapply(A) do data
        _window = maybewindow2indices(file, dims(A), window(A))
        readwindowed(data, _window)
    end

Base.getindex(A::ASCIIarray, I::Vararg{<:Union{<:Integer,<:AbstractArray}}) =
    asciiapply(A) do data
        _window = maybewindow2indices(file, dims(A), window(A))
        # Slice for both window and indices
        _dims, _refdims = slicedims(slicedims(dims(A), refdims(A), _window)..., I)
        data = readwindowed(data, _window, I...)
        rebuild(A, data, _dims, _refdims)
    end
Base.getindex(A::ASCIIarray, i1::Integer, I::Vararg{<:Integer}) =
    asciiapply(A) do data
        _window = maybewindow2indices(data, dims(A), window(A))
        readwindowed(, _window, i1, I...)
    end

# This simply loads the entire file and runs the function on the array
asciiapply(f, A) = 
    open(filename(A)) do stream
        A = Matrix{Float64}(undef, A.sourcesize...)
        for row in 1:nrow
            A[row, :] .= parse.(Float64, split(readline(stream), " "))
        end
        f(A)
    end

Base.getindex(stream::ASCIIstream, I...) = 

Base.write(filename::String, ::Type{GrdArray}, A::AbstractGeoArray{Any,2}) = begin
    # Standardise dimensions
    A = PermutedDimsArray(A, (Lon(), Lat()))
    ncols, nrows = size(A)
    x = dims(A, Lat())
    y = dims(A, Lon())
    xmin = bounds(lon)[1]
    ymin = bounds(lon)[1]
    xlocus = locus(mode(x))
    xll = if xlocus == Start() 
        "XLLCORNER"
    elseif xlocus == Center()
        "XLLCENTER"
    else
        error("$xlocus is not a valid locus for ASCII files. Use `Start` or `Center`")
    end

    ylocus = locus(mode(y))
    yll = if ylocus == Start() 
        "YLLCORNER"
    elseif ylocus == Center()
        "YLLCENTER"
    else
        error("$ylocus is not a valid locus for ASCII files. Use `Start` or `Center`")
    end

    nodata = missingval(A)

    open(filename, "w") do io
        write(io,
            """
            NCOLS $ncols
            NROWS $nrows
            $xll $xmin
            $yll $ymin
            CELLSIZE $step
            NODATA_VALUE $nodata
            """
        )
        # TODO: deal with reverse order data
        for row in 1:nrows
            write(io, join(parent(A)[row, :]), " ")
        end

    end
    return
end
