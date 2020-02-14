using GeoData, Test, Statistics, Dates
using GeoData: Time, formatdims, data, dims2indices, rebuild, window, name, source

data1 = cumsum(cumsum(ones(10, 11); dims=1); dims=2)
data2 = 2cumsum(cumsum(ones(10, 11, 1); dims=1); dims=2)
dims1 = Lon<|(10, 100), Lat<|(-50, 50) 
dims2 = (dims1..., Time<|[DateTime(2019)])
refdimz = ()
nme = "test"
mval = -9999.0
meta = nothing

# Formatting only occurs in shorthand constructors
ga1 = GeoArray(data1, formatdims(data1, dims1), refdimz, nme, meta, mval)
ga2 = GeoArray(data2, dims2)

stack = GeoStack(ga1, ga2; keys=(:ga1, :ga2))

@testset "stack layers" begin
    @test typeof(source(stack)) <: NamedTuple
    @test length(source(stack)) == 2
    @test stack[:ga1] == ga1
    @test stack[:ga2] == ga2
    @test data(stack[:ga1]) == data1
    @test typeof(data(stack[:ga1])) <: Array{Float64,2}
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
    # @test metadata(stack, :ga1) == nothing
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
    @test typeof(data(a[:ga1])) <: Array
    @test a[:ga1] == data1[2:4, 5:6]
    @test a[:ga2] == data2[2:4, 5:6, 1:1]

    @testset "select new arrays for the whole stack" begin
        s = stack[Lat<|Between(-10, 10.0), Time<|At(DateTime(2019))]
        stack[Lat<|Between(-10, 10.0), Time<|At<|DateTime(2019)]
        @test typeof(s) <: GeoStack
        @test typeof(s[:ga1]) <: GeoArray
        @test typeof(data(s[:ga1])) <: Array
        @test s[:ga1] == data1[:, 5:7]
        @test s[:ga2] == data2[:, 5:7, 1]
        @test dims(s[:ga2]) == (Lon(LinRange(10.0, 100.0, 10); grid=RegularGrid(; step=10.0)), 
                                Lat(LinRange(-10.0, 10.0, 3); grid=RegularGrid(; step=10.0)))
        @test dims(s, :ga2) == dims(s[:ga2])
        @test refdims(s[:ga2]) == (Time(DateTime(2019); grid=AlignedGrid()),)
        @test ismissing(missingval(s, :ga2)) && ismissing(missingval(s[:ga2]))
    end

    @testset "select views of arrays for the whole stack" begin
        sv = view(stack, Lat<|Between(-4.0, 27.0), Time<|At<|DateTime(2019))
        @test typeof(sv) <: GeoStack
        @test typeof(sv[:ga1]) <: GeoArray
        @test typeof(data(sv[:ga1])) <: SubArray
        @test sv[:ga1] == data1[:, 6:8]
        @test sv[:ga2] == data2[:, 6:8, 1]
        @test dims(sv[:ga2]) == (Lon(LinRange(10.0, 100.0, 10); grid=RegularGrid(; step=10.0)), 
                                 Lat(LinRange(0.0, 20.0, 3); grid=RegularGrid(; step=10.0)))
        @test refdims(sv[:ga2]) == (Time(DateTime(2019); grid=AlignedGrid()),)
        # Stack of view-based GeoArrays
        v = view(stack, Lon(2:4), Lat(5:6))
        # TODO fix type inference
        @test_broken @inferred view(stack, Lon(2:4), Lat(5:6))
        @test typeof(v) <: GeoStack
        @test typeof(v[:ga1]) <: GeoArray
        @test typeof(data(v[:ga1])) <: SubArray
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

@testset "concatenate stacks" begin
    dims1b = Lon<|(110, 200), Lat<|(-50, 50) 
    dims2b = (dims1b..., Time<|[DateTime(2019)])
    stack_a = GeoStack((ga1=ga1, ga2=ga2))
    stack_b = GeoStack((ga1=GeoArray(data1 .+ 10, dims1b), ga2=GeoArray(data2 .+ 20, dims2b)))
    catstack = cat(stack_a, stack_b; dims=Lon())
    @test size(first(catstack)) == (20, 11)
    @test val(dims(first(catstack), Lon())) ≈ 10.0:10.0:200.0
    @test step(dims(first(catstack), Lon())) == 10.0
    @test bounds(first(catstack), Lon()) == (10.0, 210.0)
    @test catstack[:ga1][Lat(1)] == 1.0:20.0
    @test catstack[:ga2][Lat(1), Time(1)] == 2.0:2.0:40.0
    catstack = cat(stack_a, stack_b; dims=(Lon(), Lat()))
    @test size(first(catstack)) == (20, 22)
end
