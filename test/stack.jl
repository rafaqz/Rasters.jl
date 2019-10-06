stack = GeoStack(ga1, ga2; keys=(:ga1, :ga2))

@testset "stack layers" begin
    @test typeof(parent(stack)) <: NamedTuple
    @test length(parent(stack)) == 2
    @test stack[:ga1] == ga1
    @test stack[:ga2] == ga2
    @test parent(stack[:ga1]) == data1
    @test typeof(parent(stack[:ga1])) <: Array{Float64,2}
    @test keys(stack) == (:ga1, :ga2)
    @test names(stack) == (:ga1, :ga2)
    @test collect(values(stack)) == [ga1, ga2]
end

@testset "stack fields " begin
    @test dims(stack) == formatdims(data1, dims1)
    @test dims(stack, :ga1) == formatdims(data1, dims1)
    @test window(stack) == ()
    @test refdims(stack) == ()
    @test metadata(stack) == nothing
    @test metadata(stack, :ga1) == nothing
end

@testset "indexing" begin
    # Indexing the stack is the same as indexing its child array
    a = stack[:ga1][Lon<|2:4, Lat<|5:6]
    @test a == stack[:ga1, Lon<|2:4, Lat<|5:6]

    @inferred stack[:ga1][Lon<|2:4, Lat<|5:6] 
    # FIXME: This isn't inferred, the constants don't propagate like they 
    # do in the above call. Probably due to the anonymous wrapper function. 
    @test_broken @inferred stack[:ga1, Lon<|2:4, Lat<|5:6] 

    # Getindex for a whole stack of new GeoArrays
    a = stack[Lon<|2:4, Lat<|5:6]
    @test typeof(a) <: GeoStack
    @test typeof(a[:ga1]) <: GeoArray
    @test typeof(parent(a[:ga1])) <: Array
    @test a[:ga1] == data1[2:4, 5:6]
    @test a[:ga2] == data2[2:4, 5:6, 1:1]

    @testset "select new arrays for the whole stack" begin
        s = stack[Lat<|Between(-10, 10.0), Time<|At<|DateTime(2019)]
        stack[Lat<|Between(-10, 10.0), Time<|At<|DateTime(2019)]
        @test typeof(s) <: GeoStack
        @test typeof(s[:ga1]) <: GeoArray
        @test typeof(parent(s[:ga1])) <: Array
        @test s[:ga1] == data1[:, 5:7]
        @test s[:ga2] == data2[:, 5:7, 1]
        @test dims(s[:ga2]) == (Lon<|LinRange(10.0, 100.0, 10), 
                                Lat<|LinRange(-10.0, 10.0, 3))
        @test dims(s, :ga2) == dims(s[:ga2])
        @test refdims(s[:ga2]) == (Time<|DateTime(2019),)
        @test ismissing(missingval(s, :ga2)) && ismissing(missingval(s[:ga2]))
    end

    @testset "select views of arrays for the whole stack" begin
        sv = view(stack, Lat<|Between(-4.0, 27.0), Time<|At<|DateTime(2019))
        @test typeof(sv) <: GeoStack
        @test typeof(sv[:ga1]) <: GeoArray
        @test typeof(parent(sv[:ga1])) <: SubArray
        @test sv[:ga1] == data1[:, 6:8]
        @test sv[:ga2] == data2[:, 6:8, 1]
        @test dims(sv[:ga2]) == (Lon<|LinRange(10.0, 100.0, 10), Lat<|LinRange(0.0, 20.0, 3))
        @test refdims(sv[:ga2]) == (Time<|DateTime(2019),)

        # Stack of view-based GeoArrays
        v = view(stack, Lon(2:4), Lat(5:6))
        # TODO fix type inference
        @test_broken @inferred view(stack, Lon(2:4), Lat(5:6))
        @test typeof(v) <: GeoStack
        @test typeof(v[:ga1]) <: GeoArray
        @test typeof(parent(v[:ga1])) <: SubArray
        @test v[:ga1] == view(data1, 2:4, 5:6)
        @test v[:ga2] == view(data2, 2:4, 5:6, 1:1)
    end
end

@testset "subset stack with specific key(s)" begin
    s1 = GeoStack(stack; keys=(:ga2,))
    @test keys(s1) == (:ga2,)
    @test length(values(s1)) == 1
    s2 = GeoStack(stack; keys=(:ga1, :ga2))
    @test keys(s2) == (:ga1, :ga2)
    @test length(values(s2)) == 2
end
