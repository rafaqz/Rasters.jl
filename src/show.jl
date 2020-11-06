
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
