using Rasters, DimensionalData, Test, Statistics, Dates, ArchGDAL
using Rasters.Lookups, Rasters.Dimensions

data1 = cumsum(cumsum(ones(10, 11); dims=1); dims=2)
data2 = 2cumsum(cumsum(ones(10, 11, 1); dims=1); dims=2)
dims1 = X(10.0:10.0:100.0), Y(-50.0:10.0:50.0)
dims2 = (dims1..., Ti([DateTime(2019)]))
refdimz = ()
nme = :test
mval = -9999.0
meta = NoMetadata()

# Formatting only occurs in shorthand constructors
raster1 = Raster(data1, dims1; refdims=refdimz, name=nme, metadata=meta, missingval=mval)
raster2 = Raster(data2, dims2)

st = RasterStack((raster1, raster2); name=(:r1, :r2))

@testset "constructors and keywords" begin
    @test_throws ArgumentError RasterStack("notastack")
    md = Dict("a" => 1)
    # Maybe too many ways to define a stack...
    kw = (; missingval=(r1=mval, r2=mval), refdims=(Ti(),), metadata=md, crs=EPSG(4326), mappedcrs=EPSG(3857), layermetadata=(r1=md, r2=md))
    st1 = RasterStack((raster1, raster2); name=(:r1, :r2), kw...)
    st2 = RasterStack((data1, data2), dims2; name=[:r1, :r2], kw...)
    st3 = RasterStack(st2; kw...)
    st4 = RasterStack(st2, dims2; kw...)
    st5 = RasterStack((r1=data1, r2=data2), dims2; kw...)
    st6 = RasterStack((; r1=raster1, r2=raster2); kw...)
    stacks = (st1, st2, st3, st4, st5, st6)
    @test st1 == st2 == st3 == st4 == st5 == st6
    @test all(==((:r1, :r2)), map(name, stacks))
    @test all(==(mval), map(missingval, stacks))
    @test all(==((Ti(),)), map(refdims, stacks))
    @test all(==(md), map(metadata, stacks))
    @test all(==(EPSG(4326)), map(crs, stacks))
    @test all(==(EPSG(3857)), map(mappedcrs, stacks))
    @test all(==((r1=md, r2=md)), map(DimensionalData.layermetadata, stacks))

    # The dimension differences are lost because the table
    # is tidy - every column is the same length
    table_st = RasterStack(DimTable(st3), dims2)
    @test dims(table_st.r1) isa Tuple{<:X,<:Y,<:Ti}
end

@testset "stack layers" begin
    @test length(layers(st)) == 2
    @test first(layers(st)) == raster1
    @test last(layers(st)) == raster2
    @test DimensionalData.layers(st) isa NamedTuple
    @test st.r1 == raster1
    @test st[:r2] == raster2
    @test parent(st[:r1]) == data1
    @test parent(st[:r1]) isa Array{Float64,2}
    @test keys(st) == (:r1, :r2)
    @test haskey(st, :r1)
    @test names(st) == (:r1, :r2)
    @test collect(values(st)) == [raster1, raster2]
end

@testset "st fields " begin
    @test DimensionalData.layerdims(st, :r1) == DimensionalData.format(dims1, data1)
    @test metadata(st) == NoMetadata()
    @test metadata(st, :r1) == NoMetadata()
end

@testset "indexing" begin
    # Indexing the st is the same as indexing its child array
    a = st[:r1][X(2:4), Y(5:6)]
    @inferred st[:r1][X=2:4, Y=5:6]

    # Getindex for a whole st of new Rasters
    a = st[X=2:4, Y=5:6]
    @test a isa RasterStack
    @test a[:r1] isa Raster
    @test parent(a[:r1]) isa Array
    @test a[:r1] == data1[2:4, 5:6]
    @test a[:r2] == data2[2:4, 5:6, 1:1]

    @testset "select new arrays for the whole st" begin
        s = st[Y=Between(-10, 10.0), Ti=At(DateTime(2019))]
        @test s isa RasterStack
        @test s.r1 isa Raster
        @test parent(s[:r1]) isa Array
        @test s[:r1] == data1[:, 5:7]
        @test s[:r2] == data2[:, 5:7, 1]
        @test dims(s[:r2]) == (X(Sampled(10.0:10.0:100.0, ForwardOrdered(), Regular(10.0), Points(), NoMetadata())),
                               Y(Sampled(-10.0:10.0:10.0, ForwardOrdered(), Regular(10.0), Points(), NoMetadata())))
        @test refdims(s[:r2]) == 
            (Ti(Sampled([DateTime(2019)], ForwardOrdered(), Irregular((DateTime(2019), DateTime(2019))), Points(), NoMetadata())),)
        @test isnothing(missingval(s, :r2)) && isnothing(missingval(s[:r2]))
    end

    @testset "select views of arrays for the whole st" begin
        sv = view(st, Y=Between(-4.0, 27.0), Ti=At(DateTime(2019)))
        @test sv isa RasterStack
        @test sv.r1 isa Raster
        @test parent(sv.r1) isa SubArray
        @test sv[:r1] == data1[:, 6:8]
        @test sv[:r2] == data2[:, 6:8, 1]
        @test dims(sv.r2) == (X(Sampled(10.0:10:100.0, ForwardOrdered(), Regular(10.0), Points(), NoMetadata())),
                               Y(Sampled(0.0:10:20.0, ForwardOrdered(), Regular(10.0), Points(), NoMetadata())))
        @test refdims(sv[:r2])[1] == 
            Ti(Sampled(view([DateTime(2019)], 1:1), ForwardOrdered(), Irregular((DateTime(2019), DateTime(2019))), Points(), NoMetadata()))
        # Stack of view-based Rasters
        v = view(st, X(2:4), Y(5:6))
        # @inferred view(st, X(2:4), Y(5:6))
        @test v isa RasterStack
        @test v[:r1] isa Raster
        @test parent(v[:r1]) isa SubArray
        @test v[:r1] == view(data1, 2:4, 5:6)
        @test v[:r2] == view(data2, 2:4, 5:6, 1:1)
    end
end

@testset "subset st with specific key(s)" begin
    s1 = RasterStack(st; name=(:r2,))
    @test keys(s1) == (:r2,)
    @test length(values(s1)) == 1
    s2 = RasterStack(st; name=(:r1, :r2))
    @test keys(s2) == (:r1, :r2)
    @test length(values(s2)) == 2
end

@testset "concatenate stacks" begin
    dims1b = X(110:10:200), Y(-50:10:50)
    dims2b = (dims1b..., Ti([DateTime(2019)]))
    stack_a = RasterStack((l1=raster1, l2=raster2))
    stack_b = RasterStack((l1=Raster(data1 .+ 10, dims1b), l2=Raster(data2 .+ 20, dims2b)))
    catstack = cat(stack_a, stack_b; dims=X)
    @test size(catstack.l1) == (20, 11)
    @test size(catstack.l2) == (20, 11, 1)
    @test val(dims(catstack, X)) ≈ 10.0:10.0:200.0
    #@test step(dims(first(catstack), X())) == 10.0
    @test DimensionalData.bounds(dims(layers(catstack, 1), X)) == (10.0, 200.0)
    @test catstack.l1[Y(1)] == 1.0:20.0
    @test catstack.l2[Y(1), Ti(1)] == 2.0:2.0:40.0
    dims2c = (dims1b..., Ti([DateTime(2019)]))
    stack_c = set(stack_b, X=>110:10:200, Y=>60:10:160)
    catstack = cat(stack_a, stack_c; dims=(X, Y))
    @test size(catstack.l1) == (20, 22)
    @test size(catstack.l2) == (20, 22, 1)
    @test dims(catstack) == 
        (X(Sampled(10:10.0:200, ForwardOrdered(), Regular(10.0), Points(), NoMetadata())),
         Y(Sampled(-50:10.0:160, ForwardOrdered(), Regular(10.0), Points(), NoMetadata())),
         Ti(Sampled([DateTime(2019)], ForwardOrdered(), Irregular(nothing, nothing), Points(), NoMetadata())))
end

@testset "copy" begin
    cp = copy(st)
    @test all(st[:r1] .=== cp[:r1])
    @test st[:r1] !== cp[:r1]
end

@testset "show" begin
    sh = sprint(show, st)
    # Test but don't lock this down too much
    @test occursin("RasterStack", sh)
    @test occursin("Y", sh)
    @test occursin("X", sh)
end
