
function DD.show_after(io::IO, mime::MIME"text/plain", A::AbstractGeoArray)
    print(io, "\nFrom file:\n$(filename(A))")
end

function DD.show_after(io, mime, stack::AbstractGeoStack) 
    if data(stack) isa FileStack 
        if filename(stack) isa AbstractString
            println(io, "\nFrom file:\n" * filename(stack))
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
