using Rasters, Test, Dates, DiskArrays, Extents, ArchGDAL, NCDatasets
using Rasters.Lookups, Rasters.Dimensions
using Rasters: isdisk, ismem, filename

@testset "create Raster" begin
    rast = Rasters.create(Int32, Extents.Extent(X=(0, 10), Y=(0, 5));
        size=(1024, 1024),
        crs=EPSG(4326),
        chunks=(X=128, Y=128),
        force=true,
        name=:testname,
        fill=Int32(2),
    ) do A
        A .*= 3
    end
    @test all(rast .=== Int32(6))
    @test crs(rast) == EPSG(4326)
    @test size(rast) == (1024, 1024)
    @test Rasters.name(rast) == :testname
    @test missingval(rast) === nothing
    @test ispoints(rast)

    rast = @test_nowarn Rasters.create(Float64, Extents.Extent(X=(0, 10), Y=(0, 5), Ti=(DateTime(2001), DateTime(2002)));
        res=(X=0.2, Y=0.1, Ti=Month(1)),
        crs=EPSG(4326),
        force=true,
        sampling=Intervals(Start()),
        name=:testname,
        missingval=missing,
        reverse_y=false,
        fill=2.0,
    ) do A
        A .*= 3
    end
    @test all(rast .=== 6.0)
    @test crs(rast) == EPSG(4326)
    @test size(rast) == (50, 50, 12)
    @test Rasters.name(rast) == :testname
    @test missingval(rast) === missing
    @test isintervals(rast)
    @test map(step, lookup(rast)) == (0.2, 0.1, Month(1))

    D = (Ti(DateTime(2000):Month(1):DateTime(2000, 12); sampling=Intervals(Start())), X(0.0:0.01:10.0), Y(0.0:0.01:10))
    rast = Rasters.create(Int32, D; fill=1, missingval=missing, crs=EPSG(4326), name=:testname)
    map(length, Rasters.dims(rast))
    @test crs(rast) == EPSG(4326)
    @test size(rast) == (12, 1001, 1001)
    @test Rasters.name(rast) == :testname
    @test missingval(rast) === missing
    @test isintervals(rast, Ti)
    @test ispoints(rast, (X, Y))
    @test map(step, lookup(rast)) == (Month(1), 0.01, 0.01)
    @test all(x -> x === Int32(1), rast)

    rast1 = Rasters.create(rast)
    @test dims(rast1) == dims(rast)
    @test eltype(rast1) == eltype(rast)
end

@testset "create RasterStack" begin
    st = Rasters.create((a=Int32, b=Float64, c=Bool), Extents.Extent(X=(0, 10), Y=(0, 5));
        size=(X=1024, Y=1024),
        sampling=(X=Points(), Y=Intervals()),
        crs=EPSG(4326),
        force=true,
        verbose=false,
        missingval=(a=Int32(-9999), b=Float64(-9999), c=false),
        fill=(a=Int32(-9999), b=0, c=false),
    ) do st
        st.c .= true
    end
    @test crs(st) == EPSG(4326)
    @test size(st) == (1024, 1024)
    @test Rasters.name(st) == (:a, :b, :c)
    @test eltype(st) === @NamedTuple{a::Int32,b::Float64, c::Bool}
    @test missingval(st) === (a=Int32(-9999), b=-9999.0, c=false)
    @test ispoints(st, X)
    @test isintervals(st, Y)
    @test all(x -> x === Int32(-9999), st.a)
    @test all(x -> x === 0.0, st.b)
    @test all(x -> x === true, st.c)

    st2 = Rasters.create((a=UInt8, b=Float32), st;
        layerdims=(a=(X(), Y()), b=(Y(),)),
        missingval=(a=UInt8(0), b=1.0f0)
    )
    @test basedims(st2.a) == (X(), Y())
    @test basedims(st2.b) == (Y(),)
    @test eltype(st2) === @NamedTuple{a::UInt8, b::Float32}
    @test missingval(st2) === (a=UInt8(0), b=1.0f0)

    @testset "from template with new dims" begin
        st1 = Rasters.create(st;
            layerdims=(a=(X, Y), b=(Y,), c=(X,)),
        )
        @test crs(st1) == EPSG(4326)
        @test size(st1) == (1024, 1024)
        @test Rasters.name(st1) == (:a, :b, :c)
        @test missingval(st1) === (a=Int32(-9999), b=-9999.0, c=false)
        @test ispoints(st1, X)
        @test isintervals(st1, Y)
        @test basedims(st1.a) == (X(), Y())
        @test basedims(st1.b) == (Y(),)
        @test basedims(st1.c) == (X(),)
    end

    @testset "from template with new layers" begin
        st1 = Rasters.create((c=UInt8, d=Int16), st;
            missingval=(c=0x00, d=Int16(1)),
        )
        @test crs(st1) == EPSG(4326)
        @test size(st1) == (1024, 1024)
        @test Rasters.name(st1) == (:c, :d)
        @test eltype(st1) == @NamedTuple{c::UInt8,d::Int16}
        @test missingval(st1) === (c=0x00, d=Int16(1))
    end

    @testset "from template with new dims and layers" begin
        st1 = Rasters.create((c=UInt8, d=Int16), st;
            layerdims=(c=(X, Y), d=(Y,)),
            missingval=(c=UInt8(0), d=Int16(1)),
        )
        @test crs(st1) == EPSG(4326)
        @test size(st1) == (1024, 1024)
        @test Rasters.name(st1) == (:c, :d)
        @test missingval(st1) === (c=0x00, d=Int16(1))
    end
end

ext = ".nc"
for ext in (".nc", ".tif", ".grd")
    @testset "create $ext" begin
        fn = "created$ext"
        created = Rasters.create(fn, UInt8, (X(1:10), Y(1:10));
            missingval=0xff,
            maskingval=nothing,
            fill=0x01,
            force=true
        )
        @test all(Raster(fn; maskingval=nothing) .=== 0x01)
        @test missingval(created) === 0xff

        if ext == ".grd"
            created = Rasters.create(fn, Int16, (X(1:10), Y(1:10));
                missingval=typemax(Int16),
                force=true,
            );
            open(created; write=true) do O
                O .= 2
                nothing
            end
            @test all(Raster(fn) .=== Int16(2))
            @test missingval(Raster(fn; maskingval=nothing)) === typemax(Int16)
        else
            @time created = Rasters.create(fn, Int16, (X(1:10), Y(1:10));
                missingval=typemax(Int16),
                scale=0.1,
                offset=5.0,
                fill=Int16(1),
                force=true,
            ) do C
                C .*= 3
            end
            @test all(Raster(fn) .=== 3.0)
            @test all(Raster(fn; scaled=false) .== Int16(-20))
            @test missingval(Raster(fn; maskingval=nothing, scaled=false)) === typemax(Int16)
        end
    end
end


@testset "create .nc stack" begin
    created = Rasters.create("created.nc", (a=UInt8, b=Float32), (X(1:10), Y(1:10));
        missingval=(a=0xff, b=typemax(Float32)),
        maskingval=nothing,
        fill=(a=0x01, b=1.0f0),
        layerdims=(a=(X,), b=(X, Y)),
        force=true,
    )
    @test missingval(created) == (a=0xff, b=typemax(Float32))
    @test size(created.a) == (10,)
    @test size(created.b) == (10, 10)
    @test all(created.a .=== 0x01)
    @test all(created.b .=== 1.0f0)
    st = RasterStack("created.nc"; maskingval=nothing)
    @test missingval(st) == (a=0xff, b=typemax(Float32))

    created = Rasters.create("created.nc", (a=UInt8, b=Float32), (X(1:10), Y(1:10));
        missingval=(a=0xff, b=typemax(Float32)),
        fill=(a=0x01, b=1.0f0),
        layerdims=(a=(X,), b=(X, Y)),
        force=true,
    )
    @test missingval(created) === missing
    @test size(created.a) == (10,)
    @test size(created.b) == (10, 10)
    @test all(created.a .=== 0x01)
    @test all(created.b .=== 1.0f0)
    st = RasterStack("created.nc"; maskingval=nothing)
    @test missingval(st) == (a=0xff, b=typemax(Float32))

    @testset "with a function" begin
        created = Rasters.create("created.nc", (a=UInt8, b=Float32), (X(1:10), Y(1:10));
            missingval=(a=0xff, b=typemax(Float32)),
            maskingval=nothing,
            fill=(a=0x01, b=1.0f0),
            layerdims=(a=(X,), b=(X, Y)),
            force=true,
        ) do st
             map(layers(st)) do A
                 A .*= 2
             end
        end
        @test all(read(created.a) .=== 0x02)
        @test all(read(created.b) .=== 2.0f0)
    end
end


