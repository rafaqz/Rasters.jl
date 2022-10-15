using Rasters, Test, Statistics, Dates, Plots, DiskArrays, RasterDataSources, CoordinateTransformations
using Rasters.LookupArrays, Rasters.Dimensions
import ArchGDAL, NCDatasets
using Rasters: FileArray, GDALfile, crs, bounds

include(joinpath(dirname(pathof(Rasters)), "../test/test_utils.jl"))
url = "https://download.osgeo.org/geotiff/samples/gdal_eg/cea.tif"
gdalpath = maybedownload(url)

@testset "array" begin
    @test_throws ArgumentError Raster("notafile.tif")

    @time gdalarray = Raster(gdalpath; name=:test)
    @time lazyarray = Raster(gdalpath; lazy=true);
    @time eagerarray = Raster(gdalpath; lazy=false);

    @testset "lazyness" begin
        # Eager is the default
        @test parent(gdalarray) isa Array
        @test parent(lazyarray) isa FileArray
        @test parent(eagerarray) isa Array
    end
    
    @testset "load from url" begin
        A = Raster("/vsicurl/" * url)
        B = Raster(url)
        @test read(A) == read(B) == gdalarray
    end

    @testset "open" begin
        @test open(A -> A[Y=1], gdalarray) == gdalarray[:, 1, :]
        tempfile = tempname() * ".tif"
        cp(gdalpath, tempfile)
        gdalwritearray = Raster(tempfile; lazy=true)
        @test_throws ArchGDAL.GDAL.GDALError gdalwritearray .*= UInt8(2)
        open(gdalwritearray; write=true) do A
            A .*= UInt8(2)
            nothing
        end
        @test Raster(tempfile) == gdalarray .* UInt8(2)
    end

    @testset "read" begin
        @time A = read(gdalarray);
        @test A isa Raster
        @test parent(A) isa Array
        A2 = zero(A)
        @time read!(gdalarray, A2);
        A3 = zero(A)
        @time read!(gdalpath, A3);
        @test A == A2 == A3
    end

    @testset "read and write band names" begin
        A = set(cat(gdalarray, gdalarray; dims=Band), Band=>1:2)
        named = set(A, Band => string.(Ref("layer_"), dims(A, Band)))
        tempfile = tempname() * ".tif"
        write(tempfile, named)
        @test parent(dims(Raster(tempfile), Band)) == ["layer_1", "layer_2"]
        @test keys(RasterStack(tempfile; layersfrom=Band)) == (:layer_1, :layer_2)
    end

    @testset "view of disk array" begin
        A = view(lazyarray, 1:10, 1:10, 1)
        @test A isa Raster
        @test parent(A) isa DiskArrays.SubDiskArray
        @test parent(parent(A)) isa Rasters.FileArray
    end

    @testset "array properties" begin
        @test size(gdalarray) == (514, 515, 1)
        @test gdalarray isa Raster{UInt8,3}
    end

    @testset "dimensions" begin
        @test length(dims(gdalarray, X)) == 514
        @test ndims(gdalarray) == 3
        @test dims(gdalarray) isa Tuple{<:X,<:Y,<:Band}
        @test lookup(gdalarray, Band) isa DimensionalData.Categorical;
        # @test span(gdalarray, (Y, X)) ==
            # (Regular(-60.02213698319351), Regular(60.02213698319374))
        @test sampling(gdalarray, (Y, X)) ==
            (Intervals(Start()), Intervals(Start()))
        @test refdims(gdalarray) == ()
        # Bounds calculated in python using rasterio
        @test all(bounds(gdalarray, Y) .≈ (4224973.143255847, 4255884.5438021915))
        @test all(bounds(gdalarray, X) .≈ (-28493.166784412522, 2358.211624949061))
    end

    @testset "other fields" begin
        # This file has an inorrect missing value
        @test missingval(gdalarray) == nothing
        @test metadata(gdalarray) isa Metadata{GDALfile,Dict{String,Any}} 
        @test basename(metadata(gdalarray)["filepath"]) == "cea.tif"
        metadata(gdalarray)["filepath"]
        @test name(gdalarray) == :test
        @test label(gdalarray) == "test"
        @test units(gdalarray) == nothing
        @test crs(dims(gdalarray, Y)) isa WellKnownText 
        @test crs(dims(gdalarray, X)) isa WellKnownText
        @test crs(gdalarray) isa WellKnownText
        @test crs(gdalarray[Y(1)]) isa WellKnownText
        @test mappedcrs(gdalarray) === nothing
        @test mappedcrs(gdalarray[Y(1), X(1)]) === nothing
        @test crs(gdalarray[Y(1), X(1)]) === nothing
    end

    @testset "indexing" begin
        @test gdalarray[Band(1)] isa Raster{UInt8,2}
        @test gdalarray[Y(1), Band(1)] isa Raster{UInt8,1}
        @test gdalarray[X(1), Band(1)] isa Raster{UInt8,1}
        @test gdalarray[X(1), Y(1), Band(1)] isa UInt8
        @test gdalarray[1, 1, 1] isa UInt8
        # Value indexed in python with rasterio
        @test gdalarray[Y(21), X(21), Band(1)] == 82
    end

    @testset "selectors" begin
        # TODO verify the value with R/gdal etc
        @test gdalarray[X(Contains(-28492)), Y(Contains(4.225e6)), Band(1)] isa UInt8
        @test gdalarray[Y(4.224e6..4.226e6), Band(1)] isa Raster
    end

    @testset "methods" begin
        @testset "mean" begin
            @test all(mean(gdalarray; dims=Y) .=== mean(parent(gdalarray); dims=2))
        end

        @testset "trim, crop, extend" begin
            a = read(replace_missing(gdalarray, zero(eltype(gdalarray))))
            a[X(1:100)] .= missingval(a)
            trimmed = trim(a)
            @test size(trimmed) == (414, 514, 1)
            cropped = Rasters.crop(a; to=trimmed)
            @test size(cropped) == (414, 514, 1)
            @test all(collect(cropped .=== trimmed))
            extended = extend(cropped; to=a)
            @test all(collect(extended .=== a))
        end

        @testset "mask and mask! to disk" begin
            msk = read(replace_missing(gdalarray, missing))
            msk[X(1:100), Y([1, 5, 95])] .= missingval(msk)
            @test !any(gdalarray[X(1:100)] .=== missingval(msk))
            masked = mask(gdalarray; with=msk)
            @test all(masked[X(1:100), Y([1, 5, 95])] .=== missingval(msk))
            tempfile = tempname() * ".tif"
            cp(gdalpath, tempfile)
            @test !all(Raster(tempfile)[X(1:100), Y([1, 5, 95])] .=== 0x00)
            open(Raster(tempfile; lazy=true); write=true) do A
                mask!(A; with=msk, missingval=0x00)
            end
            @test all(Raster(tempfile)[X(1:100), Y([1, 5, 95])] .=== 0x00)
            rm(tempfile)
        end

        @testset "polygon mask" begin
            A = read(Raster(gdalpath; name=:test))
            ds = map(dims(A)) do d
                DimensionalData.maybeshiftlocus(Center(), d)
            end
            A = set(set(A, ds), X => Points, Y => Points)
            polymask = ArchGDAL.createpolygon([
                [-20000, 4.23e6],
                [-20000, 4.245e6],
                [0.0, 4.245e6],
                [0.0, 4.23e6],
                [-20000, 4.23e6]
            ])

            rastermask = replace_missing(copy(A), missing)
            section = X((-20000.0..0.0)), Y((4.23e6..4.245e6)), Band(1)
            rastermask .= missing
            rastermask[section...] .= A[section...]
            pmasked = mask(A; with=polymask);
            revX_pmasked = reverse(mask(reverse(A; dims=X); with=polymask); dims=X);
            revY_pmasked = reverse(mask(reverse(A; dims=Y); with=polymask); dims=Y);
            perm_pmasked1 = permutedims(mask(permutedims(A, (Y, X, Band)); with=polymask), (X, Y, Band));
            perm_pmasked2 = permutedims(mask(permutedims(A, (Band, Y, X)); with=polymask), (X, Y, Band));
            rmasked = mask(A; with=rastermask)
            @test all(rmasked .=== pmasked .=== revX_pmasked .=== revY_pmasked .=== perm_pmasked1 .=== perm_pmasked2)
        end

        @testset "classify! to disk" begin
            tempfile = tempname() * ".tif"
            cp(gdalpath, tempfile)
            open(Raster(tempfile; lazy=true); write=true) do A
                classify!(A, [0x01 0xcf 0x00; 0xd0 0xff 0xff])
            end
            @test count(==(0x00), Raster(tempfile)) + count(==(0xff), Raster(tempfile)) == length(Raster(tempfile))
        end

        @testset "aggregate" begin
            ag = aggregate(mean, gdalarray, 4)
            @test ag == aggregate(mean, gdalarray, (X(4), Y(4), Band(1)))
            ag = set(ag, Band => string.(Ref("layer_"), dims(ag, Band)))
            tempfile = tempname() * ".tif"
            write(tempfile, ag)
            open(Raster(tempfile); write=true) do dst
                aggregate!(mean, dst, gdalarray, 4)
            end
            @test Raster(tempfile) == ag
        end

        @testset "mosaic" begin
            @time gdalarray = Raster(gdalpath; name=:test)
            A1 = gdalarray[X(1:300), Y(1:200)]
            A2 = gdalarray[X(57:500), Y(101:301)]
            tempfile = tempname() * ".tif"
            Afile = mosaic(first, A1, A2; missingval=0x00, atol=1e-8, filename=tempfile)
            Amem = mosaic(first, A1, A2; missingval=0x00, atol=1e-8)
            Atest = gdalarray[X(1:500), Y(1:301)]
            Atest[X(1:56), Y(201:301)] .= 0x00
            Atest[X(301:500), Y(1:100)] .= 0x00
            @test all(Atest .=== Amem .=== Afile)
        end

        @testset "rasterize round trip" begin
            A = rebuild(read(gdalarray); missingval=0x00)
            R = rasterize(A[Band(1)]; to=A, fill=:test)
            @test all(A .===  R .=== gdalarray)
            R = rasterize(A[Band(1)]; to=A, fill=:test)
            @test all(A .=== R .== gdalarray)
            B = rebuild(read(gdalarray) .= 0x00; missingval=0x00)
            rasterize!(B, read(gdalarray[Band(1)]); fill=:test)
            @test all(B .=== gdalarray |> collect)
        end

    end # methods

    @testset "conversion to Raster" begin
        geoA = gdalarray[X(1:50), Y(1:1), Band(1)]
        @test size(geoA) == (50, 1)
        @test eltype(geoA) <: UInt8
        @time geoA isa Raster{UInt8,1}
        @test dims(geoA) isa Tuple{<:X,Y}
        @test refdims(geoA) isa Tuple{<:Band}
        @test metadata(geoA) == metadata(gdalarray)
        @test missingval(geoA) == nothing
        @test name(geoA) == :test
    end

    @testset "write" begin
        gdalarray = Raster(gdalpath; name=:test);

        @testset "2d" begin
            geoA = view(gdalarray, Band(1))
            filename = tempname() * ".asc"
            @time write(filename, geoA)
            saved1 = Raster(filename)[Band(1)];
            eltype(saved1)
            eltype(geoA)
            @test all(saved1 .== geoA)
            # @test typeof(saved1) == typeof(geoA)
            @test val(dims(saved1, X)) ≈ val(dims(geoA, X))
            @test val(dims(saved1, Y)) ≈ val(dims(geoA, Y))
            @test missingval(saved1) === missingval(geoA)
            @test refdims(saved1) == refdims(geoA)
        end

        @testset "3d, with subsetting" begin
            geoA2 = gdalarray[Y(4.224e6..4.226e6), X(-28492..0)]
            filename2 = tempname() * ".tif"
            write(filename2, geoA2)
            saved2 = read(Raster(filename2; name=:test))
            @test size(saved2) == size(geoA2) == length.(dims(saved2)) == length.(dims(geoA2))
            @test refdims(saved2) == refdims(geoA2)
            #TODO test a file with more metadata
            @test val(metadata(saved2))["filepath"] == filename2
            @test missingval(saved2) === missingval(geoA2)
            @test Rasters.name(saved2) == Rasters.name(geoA2)
            @test step(lookup(dims(saved2, Y))) ≈ step(lookup(dims(geoA2, Y)))
            @test step(lookup(dims(saved2, X))) ≈ step(lookup(dims(geoA2, X)))
            @test typeof(dims(saved2)) == typeof(dims(geoA2))
            @test all(val(dims(saved2, Band)) .≈ val(dims(geoA2, Band)))
            @test all(val(dims(saved2, X)) .≈ val(dims(geoA2, X)))
            @test all(val(dims(saved2, Y)) .≈ val(dims(geoA2, Y)))
            @test parent(saved2) == parent(geoA2)
            @test typeof(saved2) == typeof(geoA2)
            filename3 = tempname() * ".tif"
            geoA3 = cat(gdalarray[Band(1)], gdalarray[Band(1)], gdalarray[Band(1)]; dims=Band(1:3))
            write(filename3, geoA3)
            saved3 = read(Raster(filename3))
            @test all(saved3 .== geoA3)
            @test val(dims(saved3, Band)) == 1:3
        end

        @testset "custom gdal options" begin
            gdalarray = Raster(gdalpath; name=:test);
            filename = tempname() * ".tif"
            write(filename, gdalarray; driver="GTiff", options=Dict("TILED"=>"YES","BLOCKXSIZE"=>"128","BLOCKYSIZE"=>"128"))
            @test isfile(filename)
            if isfile(filename)
                r = Raster(filename; lazy=true) # lazy is important here, otherwise it's not chunked
                block_x, block_y = DiskArrays.max_chunksize(DiskArrays.eachchunk(r))
                @test block_x == block_y == 128
                rm(filename)
            end
            @test_throws ArgumentError write(filename, gdalarray; driver="GTiff", options=Dict("COMPRESS"=>"FOOBAR"))
        end

        @testset "resave current" begin
            filename = tempname() * ".rst"
            write(filename, gdalarray)
            gdalarray2 = Raster(filename; lazy=true)
            write(gdalarray2)
            @test read(Raster(filename)) == read(gdalarray2)
        end

        @testset "to grd" begin
            write("testgrd.gri", gdalarray)
            grdarray = Raster("testgrd.gri")
            @test crs(grdarray) == convert(ProjString, crs(gdalarray))
            @test all(map((a, b) -> all(a .≈ b), bounds(grdarray), bounds(gdalarray)))
            @test index(grdarray, Y) ≈ index(gdalarray, Y)
            @test val(dims(grdarray, X)) ≈ val(dims(gdalarray, X))
            @test grdarray == gdalarray
        end

        @testset "from Raster" begin
            filename = tempname() * ".tiff"
            ga = Raster(rand(100, 200), (X, Y))
            write(filename, ga)
            saved = Raster(filename)
            @test all(reverse(saved[Band(1)]; dims=Y) .=== ga)
            @test saved[1, end, 1] == saved[At(1.0), At(1.0), At(1.0)]
            @test saved[100, 1, 1] == saved[At(100), At(200), At(1)]
            filename2 = tempname() * ".tif"
            ga2 = Raster(rand(100, 200), (X(Sampled(101:200)), Y(Sampled(1:200))))
            write(filename2, ga2)
            @test reverse(Raster(filename2)[Band(1)]; dims=Y) == ga2
        end

        @testset "to netcdf" begin
            filename2 = tempname() * ".nc"
            write(filename2, gdalarray[Band(1)])
            saved = Raster(filename2; crs=crs(gdalarray), mappedcrs=crs(gdalarray))
            @test size(saved) == size(gdalarray[Band(1)])
            @test saved ≈ gdalarray[Band(1)]
            clat, clon = DimensionalData.shiftlocus.(Ref(Center()), dims(gdalarray, (Y, X)))
            @test index(clat) ≈ index(saved, Y)
            @test index(clon) ≈ index(saved, X)
            @test all(bounds(saved, X) .≈ bounds(clon))
            @test all(bounds(saved, Y) .≈ bounds(clat))
            @test projectedindex(clon) ≈ projectedindex(saved, X)
            @test all(projectedbounds(clon) .≈ projectedbounds(saved, X))
            # reason lat crs conversion is less accrurate than lon TODO investigate further
            @test all(map((a, b) -> isapprox(a, b; rtol=1e-6),
                projectedindex(gdalarray, Y),
                projectedindex(DimensionalData.shiftlocus(Start(), dims(saved, Y)))
            ))
            @test all(map((a, b) -> isapprox(a, b; rtol=1e-6), projectedbounds(saved, Y),  projectedbounds(gdalarray, Y)))
        end

        @testset "write missing" begin
            A = read(replace_missing(gdalarray, missing))
            filename = tempname() * ".tif"
            write(filename, A)
            @test missingval(Raster(filename)) === typemin(UInt8)
            rm(filename)
        end

    end

    @testset "show" begin
        sh = sprint(show, MIME("text/plain"), gdalarray)
        # Test but don't lock this down too much
        @test occursin("Raster", sh)
        @test occursin("Y", sh)
        @test occursin("X", sh)
        @test occursin("Band", sh)
    end

    @testset "plot" begin # TODO write some tests for this
        gdalarray |> plot
        gdalarray[Y(1)] |> plot
    end

    @testset "nodatavalue type matches the array type" begin
        A = Raster(WorldClim{Climate}, :tavg; res="10m", month=1)
        @test typeof(missingval(A)) === eltype(A)
        @test missingval(A) === -3.4f38
    end

    @testset "rotations" begin
        am = AffineMap([60.0 20; 40 60], [first.(bounds(gdalarray, (X, Y)))...])
        ap = Rasters.AffineProjected(am; crs=crs(gdalarray))
        affine_dims = DimensionalData.format((X(ap), Y(ap), Band(1:1)), gdalarray)
        rotated = rebuild(gdalarray; dims=affine_dims)
        @test rotated[X=At(-1e4; atol=0.5), Y=Near(4.24e6), Band=1] == 0x8c
        plot(rotated)
        write("rotated.tif", rotated)
        newrotated = Raster("rotated.tif")
        plot(newrotated)
        @test rotated == newrotated
        @test lookup(rotated, X).affinemap.linear == lookup(newrotated, X).affinemap.linear
        @test lookup(rotated, X).affinemap.translation == lookup(newrotated, X).affinemap.translation
    end

    @testset "South up/ForwardOrdered Y rasters still work" begin
        # This is not common in the wild, but should work anyway
        ArchGDAL.create(
            "test.tif", driver = ArchGDAL.getdriver("GTiff"), width=240, height=180, nbands=1, dtype=Float32
        ) do dataset
            ArchGDAL.write!(dataset, rand(240, 180), 1)
        end
        rast = Raster("test.tif")
        @test order(dims(rast)) == (ForwardOrdered(), ForwardOrdered(), ForwardOrdered())
        @test span(rast) == (Regular(1.0), Regular(1.0), NoSpan())
        @test sampling(rast) == (Intervals(Start()), Intervals(Start()), NoSampling())
        @test index(rast) == (LinRange(0.0, 239.0, 240), LinRange(0.0, 179.0, 180), 1:1)
    end

end

@testset "stack" begin
    @time gdalstack = RasterStack((a=gdalpath, b=gdalpath))

    @test length(gdalstack) == 2
    @test dims(gdalstack) isa Tuple{<:X,<:Y,<:Band}

    @testset "read" begin
        st = read(gdalstack)
        read!((a=gdalpath, b=gdalpath), st)
        @test st isa RasterStack
        @test st.data isa NamedTuple
        @test first(st.data) isa Array
        @test all(st[:a] .=== gdalstack[:a])
    end

    @testset "child array properties" begin
        @test size(gdalstack[:a]) == (514, 515, 1)
        @test gdalstack[:a] isa Raster{UInt8,3}
    end

    @testset "indexing" begin
        @test gdalstack[:a][Y(2:3), X(1), Band(1)] == [0x00, 0x6b]
        @test gdalstack[:a][Y(1), X(1), Band(1)] == 0x00
        @test gdalstack[:b, Band(1)] == gdalstack[:b][Band(1)]
        @test typeof(gdalstack[:b, Band(1)]) == typeof(gdalstack[:b][Band(1)])
        @test view(gdalstack, Y(2:3), X(1), Band(1))[:a] == [0x00, 0x6b]
    end

    @testset "methods" begin
        @testset "mean" begin
            means = map(A -> mean(parent(A); dims=2), gdalstack)
            @test map((a, b) -> all(a .== b), mean(gdalstack; dims=Y), means) |> all
        end
        @testset "trim, crop, extend" begin
            mv = zero(eltype(gdalstack[:a]))
            st = read(replace_missing(gdalstack, mv))
            st = map(A -> (view(A, X(1:100)) .= mv; A), st)
            trimmed = trim(st)
            @test size(trimmed) == (414, 514, 1)
            cropped = crop(st; to=trimmed)
            @test size(cropped) == (414, 514, 1)
            @test map((c, t) -> all(collect(c .=== t)), cropped, trimmed) |> all
            extended = extend(read(cropped); to=st)
            @test all(map((s, e) -> all(s .=== e), st, extended))
        end
        @testset "mask and mask!" begin
            st = read(gdalstack)
            msk = read(replace_missing(gdalstack, missing))[:a]
            msk[X(1:100), Y([1, 5, 95])] .= missingval(msk)
            @test !any(st[:b][X(1:100)] .=== missingval(msk))
            masked = mask(st; with=msk)
            masked[:b][X(1:100), Y([1, 5, 95])]
            @test all(masked[:b][X(1:100), Y([1, 5, 95])] .=== missing)
            mask!(st; with=msk, missingval=0x00)
            @test all(st[:a][X(1:100), Y([1, 5, 95])] .=== 0x00)
            @test all(st[:b][X(1:100), Y([1, 5, 95])] .=== 0x00)
        end

        @testset "rasterize roud trip" begin
            st = map(A -> rebuild(A; missingval=0x00), gdalstack) |> read
            # We round-trip rasterise the Tables.jl form of st
            r_st = rasterize(read(gdalstack); to=st, fill=keys(gdalstack))
            @test all(map((a, b, c) -> all(a .=== b .=== c), st, r_st, read(gdalstack)))
            r_st = rasterize(read(gdalstack); to=st, fill=(:a, :b))
            @test all(map((a, b, c) -> all(a .=== b .=== c), st, r_st, read(gdalstack)))
            st = map(A -> rebuild(A .* 0x00; missingval=0x00), gdalstack) |> read
            rasterize!(st, read(gdalstack), fill=keys(st))
            @test all(map((a, b) -> all(a .=== b), st, gdalstack))

            # Getting the band column works if we force it
            # name of Symbol gives a Raster, Tuple gives a RasterStack
            b_r = rasterize(read(gdalstack); to=st, fill=:Band)
            @test b_r isa Raster
            b_st = rasterize(read(gdalstack); to=st, fill=(:Band, ))
            @test b_st isa RasterStack
            @test b_r == b_st[:Band]
        end

        @testset "classify" begin
            cstack = classify(gdalstack, 0x01..0xcf => 0x00, 0xd0..0xff => 0xff)
            @test count(==(0x00), cstack[:a]) + count(==(0xff), cstack[:a]) == length(gdalstack[:a])
            @test count(==(0x00), cstack[:b]) + count(==(0xff), cstack[:b]) == length(gdalstack[:b])
            cstack = copy(gdalstack)
            cstack = classify!(cstack, 0x01..0xcf => 0x00, 0xd0..0xff => 0xff)
            @test count(==(0x00), cstack[:a]) + count(==(0xff), cstack[:a]) == length(gdalstack[:a])
            @test count(==(0x00), cstack[:b]) + count(==(0xff), cstack[:b]) == length(gdalstack[:b])
        end
    end

    if VERSION > v"1.1-"
        @testset "copy" begin
            geoA = zero(Raster(gdalstack[:a]))
            copy!(geoA, gdalstack, :a)
            # First wrap with Raster() here or == loads from disk for each cell.
            # we need a general way of avoiding this in all disk-based sources
            @test all(geoA .== Raster(gdalstack[:a]))
        end
    end

    @testset "write" begin
        geoA = read(gdalstack[:a])
        @testset "write multiple files" begin
            filename = tempname() * ".tif"
            write(filename, gdalstack)
            base, ext = splitext(filename)
            filename_b = string(base, "_b", ext)
            saved = read(Raster(filename_b))
            @test all(saved .== geoA)
        end

        @testset "write multiple files with custom suffix" begin
            filename = tempname() * ".tif"
            write(filename, gdalstack; suffix=("_first", "_second"))
            base, ext = splitext(filename)
            filename_b = string(base, "_second", ext)
            saved = read(Raster(filename_b))
            @test all(saved .== geoA)
        end

        @testset "write netcdf" begin
            filename = tempname() * ".nc"
            write(filename, gdalstack);
            saved = RasterStack(filename);
            @test all(read(saved[:a]) .== geoA)
            rm(filename)
        end

    end

    @testset "show" begin
        sh = sprint(show, MIME("text/plain"), gdalstack)
        # Test but don't lock this down too much
        @test occursin("RasterStack", sh)
        @test occursin("Y", sh)
        @test occursin("X", sh)
        @test occursin("Band", sh)
        @test occursin(":a", sh)
        @test occursin(":b", sh)
    end

end

@testset "resample and warp" begin
    gdalarray = read(Raster(gdalpath; name=:test))
    gdalstack = read(RasterStack((a=Raster(gdalpath), b=Raster(gdalpath) .* Int32(2))))
    gdalser = read(RasterSeries([gdalpath, gdalpath], (Ti(),); mappedcrs=EPSG(4326), name=:test))

    output_res = 0.0027
    output_crs = EPSG(4326)
    resample_method = "near"

    ## Resample cea.tif manually with ArchGDAL
    wkt = convert(String, convert(WellKnownText, output_crs))
    AG_output = ArchGDAL.read(gdalpath) do dataset
        ArchGDAL.gdalwarp([dataset], ["-t_srs", "$(wkt)",
                                "-tr", "$(output_res)", "$(output_res)",
                                "-r", "$(resample_method)"]) do warped
            ArchGDAL.read(ArchGDAL.getband(warped, 1))
        end
    end

    ## Resample cea.tif using resample
    raster_output = resample(gdalarray, output_res; crs=output_crs, method=resample_method)
    disk_output = resample(gdalarray, output_res; crs=output_crs, method=resample_method, filename="resample.tif")
    stack_output = resample(gdalstack, output_res; crs=output_crs, method=resample_method)
    written_stack_output = resample(gdalstack, output_res; crs=output_crs, method=resample_method, filename="resample.tif")
    series_output = resample(gdalser, output_res; crs=output_crs, method=resample_method)

    extradim_raster = cat(gdalarray, gdalarray, gdalarray; dims=Z)
    extradim_output = resample(extradim_raster, output_res; crs=output_crs, method=resample_method)

    permuted_raster = permutedims(gdalarray, (Y, X, Band))
    permuted_output = resample(permuted_raster, output_res; crs=output_crs, method=resample_method)

    # Compare ArchGDAL, resample and permuted resample 
    @test AG_output ==
        raster_output[Band(1)] == disk_output[Band(1)] ==
        stack_output[:a][Band(1)] ==
        written_stack_output[:a][Band(1)] ==
        series_output[1][Band(1)] ==
        extradim_output[Z(3), Band(1)] ==
        permutedims(permuted_output, (X, Y, Band))[Band(1)]

    @test stack_output[:b][Band(1)] == written_stack_output[:b][Band(1)] == AG_output .* 2
    @test abs(step(dims(raster_output, Y))) ≈
        abs(step(dims(raster_output, X))) ≈ 
        abs(step(dims(disk_output, X))) ≈ 
        abs(step(dims(disk_output, Y))) ≈ 
        abs(step(dims(stack_output[:a], X))) ≈ 
        abs(step(dims(stack_output[:a], Y))) ≈ 
        abs(step(dims(stack_output[:b], X))) ≈ 
        abs(step(dims(stack_output[:b], Y))) ≈ 
        abs(step(dims(series_output[1], X))) ≈ 
        abs(step(dims(series_output[2], Y))) ≈ 
        abs(step(dims(series_output[1], X))) ≈ 
        abs(step(dims(series_output[2], Y))) ≈ 
        abs(step(dims(extradim_output, X))) ≈ 
        abs(step(dims(extradim_output, Y))) ≈ 
        abs(step(dims(permuted_output, X))) ≈ 
        abs(step(dims(permuted_output, Y))) ≈ output_res

    rm("resample.tif")
    rm("resample_a.tif")
    rm("resample_b.tif")

    @testset "snapped size and dim index match" begin
        snaptarget = aggregate(Center(), read(gdalarray), 2)
        snapped = resample(read(gdalarray); to=snaptarget)
        disk_snapped = resample(gdalarray; to=snaptarget, filename="snap_resample.tif")
        stack_snapped = resample(read(gdalstack); to=snaptarget, filename="snap_resample.tif")
        ser_snapped = resample(read(gdalser); to=snaptarget)
        extradim_snapped = resample(extradim_raster; to=snaptarget)
        @test size(snapped) == size(disk_snapped) == size(snaptarget)
        @test isapprox(index(snaptarget, Y), index(snapped, Y))
        @test isapprox(index(snaptarget, X), index(snapped, X))
        @test isapprox(index(snaptarget, Y), index(stack_snapped, Y))
        @test isapprox(index(snaptarget, X), index(stack_snapped, X))
        @test isapprox(index(snaptarget, Y), index(first(ser_snapped), Y))
        @test isapprox(index(snaptarget, X), index(first(ser_snapped), X))
        @test isapprox(index(snaptarget, Y), index(disk_snapped, Y))
        @test isapprox(index(snaptarget, X), index(disk_snapped, X))
        @test isapprox(index(snaptarget, Y), index(extradim_snapped, Y))
        @test isapprox(index(snaptarget, X), index(extradim_snapped, X))
        rm("snap_resample.tif")
        rm("snap_resample_a.tif")
        rm("snap_resample_b.tif")
    end
end

@testset "series" begin
    gdalser = RasterSeries([gdalpath, gdalpath], (Ti(),); mappedcrs=EPSG(4326), name=:test)
    @test read(gdalser[Ti(1)]) == read(Raster(gdalpath; mappedcrs=EPSG(4326), name=:test))
    @test read(gdalser[Ti(1)]) == read(Raster(gdalpath; mappedcrs=EPSG(4326), name=:test))

    stackser = RasterSeries((a=[gdalpath, gdalpath], b=[gdalpath, gdalpath]), (Ti,))
    # Rebuild the ser by wrapping the disk data in Array.
    # `modify` forces `rebuild` on all containers as in-Memory variants
    modified_ser = modify(Array, stackser)
    @test typeof(modified_ser) <: RasterSeries{<:RasterStack{<:NamedTuple{(:a,:b),<:Tuple{<:Array{UInt8,3},Vararg}}}}

    @testset "read" begin
        ser1 = read(stackser)
        @test ser1 isa RasterSeries{<:RasterStack}
        @test ser1.data isa Vector{<:RasterStack}
        @test first(ser1.data[1].data) isa Array
        ser2 = modify(A -> A .* 0, ser1)
        ser3 = modify(A -> A .* 0, ser1)
        read!([(a=gdalpath, b=gdalpath), (a=gdalpath, b=gdalpath)], ser2)
        read!(ser1, ser3)
        @test map(ser1, ser2, ser3) do st1, st2, st3
            map(st1, st2, st3) do A1, A2, A3
                (A2 .=== A2 .=== A3) |> all
            end |> all
        end |> all
    end

    @testset "methods" begin
        @testset "classify" begin
            cser = classify(gdalser, 0x01..0xcf => 0x00, 0xd0..0xff => 0xff)
            @test count(==(0x00), cser[1]) + count(==(0xff), cser[1]) == length(gdalser[1])
            @test count(==(0x00), cser[2]) + count(==(0xff), cser[2]) == length(gdalser[2])
            cser = copy.(gdalser)
            cser = classify!(cser, 0x01..0xcf => 0x00, 0xd0..0xff => 0xff)
            @test count(==(0x00), cser[1]) + count(==(0xff), cser[1]) == length(gdalser[1])
            @test count(==(0x00), cser[2]) + count(==(0xff), cser[2]) == length(gdalser[2])
        end
        @testset "trim, crop, extend" begin
            mv = zero(eltype(gdalser[1]))
            ser = read(replace_missing(gdalser, mv))
            ser = map(A -> (view(A, X(1:100)) .= mv; A), ser)
            trimmed = trim(ser)
            @test size(trimmed[1]) == (414, 514, 1)
            cropped = crop(ser; to=trimmed[1])
            @test size(cropped[1]) == (414, 514, 1)
            @test map((c, t) -> all(collect(c .=== t)), cropped, trimmed) |> all
            extended = extend(read(cropped); to=ser[1])
            @test all(map((s, e) -> all(s .=== e), ser, extended))
        end
        @testset "replace_missing" begin
            mser = map(x -> rebuild(x; missingval=0x00), read(gdalser))
            repser = replace_missing(mser)
            @test eltype(first(repser)) == Union{Missing,UInt8}
            replace_missing!(repser, 0xff)
            @test count(x -> x == 0xff, first(repser)) == 
                count(x -> x == 0xff, first(mser)) + count(x -> x == 0x00, first(mser))
        end
        @testset "mask and mask!" begin
            ser = read(gdalser)
            msk = first(read(replace_missing(gdalser, missing)))
            msk[X(1:100), Y([1, 5, 95])] .= missingval(msk)
            @test !any(ser[1][X(1:100)] .=== missingval(msk))
            masked = mask(ser; with=msk)
            masked[1][X(1:100), Y([1, 5, 95])]
            @test all(masked[1][X(1:100), Y([1, 5, 95])] .=== missing)
            mask!(ser; with=msk, missingval=0x00)
            @test all(ser[1][X(1:100), Y([1, 5, 95])] .=== 0x00)
            @test all(ser[2][X(1:100), Y([1, 5, 95])] .=== 0x00)
        end
    end

    @testset "show" begin
        sh = sprint(show, MIME("text/plain"), stackser)
        # Test but don't lock this down too much
        @test occursin("RasterSeries", sh)
        @test occursin("RasterStack", sh)
        @test occursin("Ti", sh)
    end

    @testset "crs" begin
        @time gdalarray = Raster(gdalpath; mappedcrs=EPSG(4326), name=:test)
        wkt = WellKnownText(GeoFormatTypes.CRS(), "PROJCS[\"unnamed\",GEOGCS[\"NAD27\",DATUM[\"North_American_Datum_1927\",SPHEROID[\"Clarke 1866\",6378206.4,294.978698213898,AUTHORITY[\"EPSG\",\"7008\"]],AUTHORITY[\"EPSG\",\"6267\"]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4267\"]],PROJECTION[\"Cylindrical_Equal_Area\"],PARAMETER[\"standard_parallel_1\",33.75],PARAMETER[\"central_meridian\",-117.333333333333],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],AXIS[\"Easting\",EAST],AXIS[\"Northing\",NORTH]]")
        @test crs(dims(gdalarray, Y)) == wkt
        @test crs(dims(gdalarray, X)) == wkt
        @test crs(gdalarray) == wkt
        @test crs(gdalarray[Y(1)]) == wkt
    end
end

