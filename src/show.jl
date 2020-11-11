
# Don't show the data in DiskGeoArray. It 
# defeats the purpose of loading them lazily.
Base.show(io::IO, A::DiskGeoArray) = begin
    l = nameof(typeof(A))
    printstyled(io, nameof(typeof(A)); color=:blue)
    if label(A) != ""
        print(io, " (named ")
        printstyled(io, label(A); color=:blue)
        print(io, ")")
    end
    print(io, " with dimensions:\n")
    for d in dims(A)
        print(io, " ", d, "\n")
    end
    if !isempty(refdims(A))
        print(io, "and referenced dimensions:\n")
        for d in refdims(A)
            print(io, " ", d, "\n")
        end
    end
    print(io, "\n  From file: $(filename(A))")
end

function Base.show(io::IO, stack::AbstractGeoStack)
    printstyled(io, "$(Base.typename(typeof(stack)))", color=:blue)
    print(io, " with $(length(keys(stack))) field(s): $(stack.filename)\n")

    for var in keys(stack)
        printstyled(io, " $var", color=:green)

        field_dims = dims(stack, var)
        if length(field_dims) > 0
            print(io, " with dimension(s) ")
            for (d, dim) in enumerate(field_dims)
                printstyled(io, "$(GeoData.name(dim))($(length(dim)))", color=:red)
                d != length(field_dims) && print(io, ", ")
            end
        end
        print(io, '\n')
    end

    n_windows = length(stack.window)
    if n_windows > 0
        print(io, "and with $n_windows window(s):\n")
        for window in stack.window
            print(io, ' ')
            show(window)
            print(io, '\n')
        end
    end

    n_metadata = length(stack.metadata)
    if n_metadata > 0
        print(io, "and $n_metadata metadata entries:\n")
        show(io, stack.metadata)
    end
end
