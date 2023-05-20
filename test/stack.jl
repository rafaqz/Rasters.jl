using Rasters, DimensionalData, Test, Statistics, Dates
using Rasters.LookupArrays, Rasters.Dimensions

data1 = cumsum(cumsum(ones(10, 11); dims=1); dims=2)
data2 = 2cumsum(cumsum(ones(10, 11, 1); dims=1); dims=2)
dims1 = X(10.0:10.0:100.0), Y(-50.0:10.0:50.0)
dims2 = (dims1..., Ti([DateTime(2019)]))
refdimz = ()
nme = :test
mval = -9999.0
meta = NoMetadata()

# Formatting only occurs in shorthand constructors
ga1 = Raster(data1, dims1; refdims=refdimz, name=nme, metadata=meta, missingval=mval)
ga2 = Raster(data2, dims2)

@testset "stack constructors" begin
    @test_throws ArgumentError RasterStack("notastack")
    # Maybe too many ways to define a stack...
    st1 = RasterStack((ga1, ga2); name=(:ga1, :ga2))
    st2 = RasterStack((ga1, ga2); keys=(:ga1, :ga2))
    st3 = RasterStack((data1, data2), dims2; keys=[:ga1, :ga2])
    st4 = RasterStack(st3)
    st5 = RasterStack(st3, dims2)
    st6 = RasterStack(RasterStack(ga1)..., RasterStack(ga2)...; name=(:ga1, :ga2))
    st7 = RasterStack(RasterStack(ga1; name=:ga1)..., RasterStack(ga2; name=:ga2)...)
    st8 = RasterStack((; ga1, ga2))
    st9 = RasterStack((ga1=data1, ga2=data2), dims2)
    @test st1 == st2 == st3 == st4 == st5 == st6 == st7 == st8 == st9

    # The dimension differences are lost because the table
    # is tidy - every column is the same length
    table_st = RasterStack(DimTable(st3), dims2)
    @test dims(table_st[:ga1]) isa Tuple{<:X,<:Y,<:Ti}
end

st = RasterStack((ga1, ga2); name=(:ga1, :ga2))

@testset "stack layers" begin
    @test length(st) == 2
    @test first(st) == ga1
    @test last(st) == ga2
    @test DimensionalData.layers(st) isa NamedTuple
    @test st.ga1 == ga1
    @test st[:ga2] == ga2
    @test parent(st[:ga1]) == data1
    @test parent(st[:ga1]) isa Array{Float64,2}
    @test keys(st) == (:ga1, :ga2)
    @test haskey(st, :ga1)
    @test names(st) == (:ga1, :ga2)
    @test collect(values(st)) == [ga1, ga2]
    # We can iterate/splat stacks as GeArrays
    @test st == RasterStack(st...)
end

@testset "st fields " begin
    @test DimensionalData.layerdims(st, :ga1) == DimensionalData.format(dims1, data1)
    @test metadata(st) == NoMetadata()
    @test metadata(st, :ga1) == NoMetadata()
end

@testset "indexing" begin
    # Indexing the st is the same as indexing its child array
    a = st[:ga1][X(2:4), Y(5:6)]
    @test a == st[:ga1, X(2:4), Y(5:6)]

    @inferred st[:ga1][X=2:4, Y=5:6]
    # FIXME: This isn't inferred, the constants don't propagate like they
    # do in the above call. Probably due to the anonymous wrapper function.
    if VERSION < v"1.9.0-"
        @test_broken @inferred st[:ga1, X(2:4), Y(5:6)] isa Raster
    else
        @test @inferred st[:ga1, X(2:4), Y(5:6)] isa Raster
    end

    # Getindex for a whole st of new Rasters
    a = st[X=2:4, Y=5:6]
    @test a isa RasterStack
    @test a[:ga1] isa Raster
    @test parent(a[:ga1]) isa Array
    @test a[:ga1] == data1[2:4, 5:6]
    @test a[:ga2] == data2[2:4, 5:6, 1:1]

    @testset "select new arrays for the whole st" begin
        s = st[Y=Between(-10, 10.0), Ti=At(DateTime(2019))]
        @test s isa RasterStack
        @test s.ga1 isa Raster
        @test parent(s[:ga1]) isa Array
        @test s[:ga1] == data1[:, 5:7]
        @test s[:ga2] == data2[:, 5:7, 1]
        @test dims(s[:ga2]) == (X(Sampled(10.0:10.0:100.0, ForwardOrdered(), Regular(10.0), Points(), NoMetadata())),
                                Y(Sampled(-10.0:10.0:10.0, ForwardOrdered(), Regular(10.0), Points(), NoMetadata())))
        @test refdims(s[:ga2]) == (Ti(Sampled([DateTime(2019)], ForwardOrdered(), Irregular(nothing, nothing), Points(), NoMetadata())),)
        @test ismissing(missingval(s, :ga2)) && ismissing(missingval(s[:ga2]))
    end

    @testset "select views of arrays for the whole st" begin
        sv = view(st, Y=Between(-4.0, 27.0), Ti=At(DateTime(2019)))
        @test sv isa RasterStack
        @test sv.ga1 isa Raster
        @test parent(sv.ga1) isa SubArray
        @test sv[:ga1] == data1[:, 6:8]
        @test sv[:ga2] == data2[:, 6:8, 1]
        @test dims(sv.ga2) == (X(Sampled(10.0:10:100.0, ForwardOrdered(), Regular(10.0), Points(), NoMetadata())),
                               Y(Sampled(0.0:10:20.0, ForwardOrdered(), Regular(10.0), Points(), NoMetadata())))
        @test refdims(sv[:ga2])[1] == Ti(Sampled(view([DateTime(2019)], 1:1), ForwardOrdered(), Irregular((nothing, nothing)), Points(), NoMetadata()))
        # Stack of view-based Rasters
        v = view(st, X(2:4), Y(5:6))
        # @inferred view(st, X(2:4), Y(5:6))
        @test v isa RasterStack
        @test v[:ga1] isa Raster
        @test parent(v[:ga1]) isa SubArray
        @test v[:ga1] == view(data1, 2:4, 5:6)
        @test v[:ga2] == view(data2, 2:4, 5:6, 1:1)
    end

end

@testset "subset st with specific key(s)" begin
    s1 = RasterStack(st; keys=(:ga2,))
    @test keys(s1) == (:ga2,)
    @test length(values(s1)) == 1
    s2 = RasterStack(st; keys=(:ga1, :ga2))
    @test keys(s2) == (:ga1, :ga2)
    @test length(values(s2)) == 2
end

@testset "concatenate stacks" begin
    dims1b = X(110:10:200), Y(-50:10:50)
    dims2b = (dims1b..., Ti([DateTime(2019)]))
    stack_a = RasterStack((ga1=ga1, ga2=ga2))
    stack_b = RasterStack((ga1=Raster(data1 .+ 10, dims1b), ga2=Raster(data2 .+ 20, dims2b)))
    catstack = cat(stack_a, stack_b; dims=X())
    @test size(first(catstack)) == (20, 11)
    @test val(dims(catstack, X)) ≈ 10.0:10.0:200.0
    #@test step(dims(first(catstack), X())) == 10.0
    @test DimensionalData.bounds(dims(first(catstack), X)) == (10.0, 200.0)
    @test catstack[:ga1][Y(1)] == 1.0:20.0
    @test catstack[:ga2][Y(1), Ti(1)] == 2.0:2.0:40.0
    catstack = cat(stack_a, stack_b; dims=(X(), Y()))
    @test size(first(catstack)) == (20, 22)
end

@testset "copy" begin
    cp = copy(st)
    @test all(st[:ga1] .=== cp[:ga1])
    @test st[:ga1] !== cp[:ga1]
end

@testset "show" begin
    sh = sprint(show, st)
    # Test but don't lock this down too much
    @test occursin("RasterStack", sh)
    @test occursin("Y", sh)
    @test occursin("X", sh)
end
