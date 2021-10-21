using GeoData, Test, Statistics, Dates, Plots, DiskArrays, RasterDataSources
import ArchGDAL, NCDatasets
using GeoData: mode, span, sampling, name, bounds, FileArray, GDALfile

include(joinpath(dirname(pathof(GeoData)), "../test/test_utils.jl"))
url = "https://download.osgeo.org/geotiff/samples/gdal_eg/cea.tif"
gdalpath = maybedownload(url)

@testset "array" begin

    @time gdalarray = GeoArray(gdalpath; mappedcrs=EPSG(4326), name=:test)
    
    @testset "load from url" begin
        A = GeoArray("/vsicurl/" * url)
        B = GeoArray(url)
        @test read(A) == read(B) == gdalarray
    end

    @testset "open" begin
        @test open(A -> A[Y=1], gdalarray) == gdalarray[:, 1, :]
        tempfile = tempname() * ".tif"
        cp(gdalpath, tempfile)
        gdalwritearray = GeoArray(tempfile)
        @test_throws ArchGDAL.GDAL.GDALError gdalwritearray .*= UInt8(2)
        open(gdalwritearray; write=true) do A
            A .*= UInt8(2)
            nothing
        end
        @test GeoArray(tempfile) == gdalarray .* UInt8(2)
    end

    @testset "read" begin
        @time A = read(gdalarray);
        @test A isa GeoArray
        @test parent(A) isa Array
        A2 = zero(A)
        @time read!(gdalarray, A2);
        A3 = zero(A)
        @time read!(gdalpath, A3);
        @test A == A2 == A3
    end

    @testset "view" begin
        A = view(gdalarray, 1:10, 1:10, 1)
        @test A isa GeoArray
        @test parent(A) isa DiskArrays.SubDiskArray
        @test parent(parent(A)) isa GeoData.FileArray
    end

    @testset "array properties" begin
        @test size(gdalarray) == (514, 515, 1)
        @test gdalarray isa GeoArray{UInt8,3}
    end

    @testset "dimensions" begin
        @test length(dims(gdalarray, X)) == 514
        @test ndims(gdalarray) == 3
        @test dims(gdalarray) isa Tuple{<:X,<:Y,<:Band}
        @test mode(gdalarray, Band) == DimensionalData.Categorical(Ordered())
        @test span(gdalarray, (Y, X)) ==
            (Regular(-60.02213698319351), Regular(60.02213698319374))
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
        @test metadata(gdalarray) isa Metadata{GDALfile}
        @test basename(metadata(gdalarray).val[:filepath]) == "cea.tif"
        @test name(gdalarray) == :test
        @test label(gdalarray) == "test"
        @test units(gdalarray) == nothing
        @test mappedcrs(dims(gdalarray, Y)) == EPSG(4326)
        @test mappedcrs(dims(gdalarray, X)) == EPSG(4326)
        @test mappedcrs(gdalarray) == EPSG(4326)
        @test mappedcrs(gdalarray[Y(1)]) == EPSG(4326)
        @test_throws ErrorException mappedcrs(gdalarray[Y(1), X(1)])
        wkt = WellKnownText(GeoFormatTypes.CRS(),
          "PROJCS[\"unnamed\",GEOGCS[\"NAD27\",DATUM[\"North_American_Datum_1927\",SPHEROID[\"Clarke 1866\",6378206.4,294.978698213898,AUTHORITY[\"EPSG\",\"7008\"]],AUTHORITY[\"EPSG\",\"6267\"]],PRIMEM[\"Greenwich\",0],UNIT[\"degree\",0.0174532925199433,AUTHORITY[\"EPSG\",\"9122\"]],AUTHORITY[\"EPSG\",\"4267\"]],PROJECTION[\"Cylindrical_Equal_Area\"],PARAMETER[\"standard_parallel_1\",33.75],PARAMETER[\"central_meridian\",-117.333333333333],PARAMETER[\"false_easting\",0],PARAMETER[\"false_northing\",0],UNIT[\"metre\",1,AUTHORITY[\"EPSG\",\"9001\"]],AXIS[\"Easting\",EAST],AXIS[\"Northing\",NORTH]]")
        @test crs(dims(gdalarray, Y)) == wkt
        @test crs(dims(gdalarray, X)) == wkt
        @test crs(gdalarray) == wkt
        @test crs(gdalarray[Y(1)]) == wkt
        @test_throws ErrorException crs(gdalarray[Y(1), X(1)])
    end

    @testset "indexing" begin
        @test gdalarray[Band(1)] isa GeoArray{UInt8,2}
        @test gdalarray[Y(1), Band(1)] isa GeoArray{UInt8,1}
        @test gdalarray[X(1), Band(1)] isa GeoArray{UInt8,1}
        @test gdalarray[X(1), Y(1), Band(1)] isa UInt8
        @test gdalarray[1, 1, 1] isa UInt8
        # Value indexed in python with rasterio
        @test gdalarray[Y(21), X(21), Band(1)] == 82
    end

    @testset "selectors" begin
        # TODO verify the value with R/gdal etc
        @test gdalarray[Y(Contains(33.8)), X(Contains(-117.5)), Band(1)] isa UInt8
        @test gdalarray[Y(Between(33.7, 33.9)), Band(1)] isa GeoArray
    end

    @testset "methods" begin
        @testset "mean" begin
            @test all(mean(gdalarray; dims=Y) .=== mean(parent(gdalarray); dims=2))
        end

        @testset "trim, crop, extend" begin
            a = replace_missing(gdalarray, zero(eltype(gdalarray)))
            a[X(1:100)] .= missingval(a)
            trimmed = trim(a)
            @test size(trimmed) == (414, 514, 1)
            cropped = crop(a; to=trimmed)
            @test size(cropped) == (414, 514, 1)
            @test all(collect(cropped .=== trimmed))
            extended = extend(cropped; to=a)
            @test all(collect(extended .=== a))
        end

        @testset "mask and mask! to disk" begin
            msk = replace_missing(gdalarray, missing)
            msk[X(1:100), Y([1, 5, 95])] .= missingval(msk)
            @test !any(gdalarray[X(1:100)] .=== missingval(msk))
            masked = mask(gdalarray; to=msk)
            @test all(masked[X(1:100), Y([1, 5, 95])] .=== missingval(msk))
            tempfile = tempname() * ".tif"
            cp(gdalpath, tempfile)
            @test !all(GeoArray(tempfile)[X(1:100), Y([1, 5, 95])] .=== 0x00)
            open(GeoArray(tempfile); write=true) do A
                mask!(A; to=msk, missingval=0x00)
            end
            @test all(GeoArray(tempfile)[X(1:100), Y([1, 5, 95])] .=== 0x00)
            rm(tempfile)
        end

        @testset "polygon mask" begin
            A = read(GeoArray(gdalpath; name=:test))
            ds = map(dims(A)) do d
                DimensionalData.maybeshiftlocus(Center(), d)
            end
            A = set(A, ds)
            A = set(A, X => Points, Y=> Points)
            polymask = [[-20000, 4.23e6],
                        [-20000, 4.245e6],
                        [0.0, 4.245e6],
                        [0.0, 4.23e6],
                        [-20000, 4.23e6]]
            rastermask = replace_missing(copy(A), missing)
            # mask or `Between` is not exactly accurate yet...
            section = X(Between(-20000.0, 0.0)), Y(Between(4.23005e6, 4.24495e6)), Band(1)
            rastermask .= missing
            rastermask[section...] .= A[section...]
            pmasked = mask(A; to=polymask) |> read
            rmasked = mask(A; to=rastermask) |> read
            @test all(rmasked .=== pmasked)
        end

        @testset "classify! to disk" begin
            tempfile = tempname() * ".tif"
            cp(gdalpath, tempfile)
            open(GeoArray(tempfile); write=true) do A
                classify!(A, [0x01 0xcf 0x00; 0xd0 0xff 0xff])
            end
            @test count(==(0x00), GeoArray(tempfile)) + count(==(0xff), GeoArray(tempfile)) == length(GeoArray(tempfile))
        end

        @testset "aggregate" begin
            ag = aggregate(mean, gdalarray, 4)
            @test ag == aggregate(mean, gdalarray, (X(4), Y(4), Band(1)))
            tempfile = tempname() * ".tif"
            write(tempfile, ag)
            GeoArray(tempfile)
            open(GeoArray(tempfile); write=true) do dst
                aggregate!(mean, dst, gdalarray, 4)
            end
            @test GeoArray(tempfile) == ag
        end

        @testset "mosaic" begin
            @time gdalarray = GeoArray(gdalpath; name=:test)
            A1 = gdalarray[X(1:300), Y(1:200)]
            A2 = gdalarray[X(57:500), Y(101:301)]
            tempfile = tempname() * ".tif"
            Afile = mosaic(first, A1, A2; missingval=0x00, atol=1e-7, filename=tempfile)
            Amem = mosaic(first, A1, A2; missingval=0x00, atol=1e-7)
            Atest = gdalarray[X(1:500), Y(1:301)]
            Atest[X(1:56), Y(201:301)] .= 0x00
            Atest[X(301:500), Y(1:100)] .= 0x00
            @test all(Atest .=== Amem .=== Afile)
        end

        @testset "rasterize" begin
            A = rebuild(read(gdalarray); missingval=0x00)
            # We round-trip rasterise the Tables.jl form of A
            R = rasterize(A; to=A)
            @test all(A .=== R .== gdalarray)
            R = rasterize(A; to=A, name=:test)
            @test all(A .=== R .== gdalarray)
            B = rebuild(read(gdalarray) .= 0x00; missingval=0x00)
            rasterize!(B, read(gdalarray))
            @test all(B .=== gdalarray |> collect)
        end

        @testset "chunk_series" begin
            @test GeoData.chunk_series(gdalarray) isa GeoSeries
            @test size(GeoData.chunk_series(gdalarray)) == (1, 1, 1)
        end
    end

    @testset "conversion to GeoArray" begin
        geoA = gdalarray[X(1:50), Y(1:1), Band(1)]
        @test size(geoA) == (50, 1)
        @test eltype(geoA) <: UInt8
        @time geoA isa GeoArray{UInt8,1}
        @test dims(geoA) isa Tuple{<:X,Y}
        @test refdims(geoA) isa Tuple{<:Band}
        @test metadata(geoA) == metadata(gdalarray)
        @test missingval(geoA) == nothing
        @test name(geoA) == :test
    end

    @testset "write" begin
        gdalarray = GeoArray(gdalpath; mappedcrs=EPSG(4326), name=:test);

        @testset "2d" begin
            geoA = view(gdalarray, Band(1))
            filename = tempname() * ".asc"
            @time write(filename, geoA)
            saved1 = GeoArray(filename; mappedcrs=EPSG(4326))[Band(1)];
            eltype(saved1)
            eltype(geoA)
            @test all(saved1 .== geoA)
            # @test typeof(saved1) == typeof(geoA)
            @test val(dims(saved1, X)) ≈ val(dims(geoA, X))
            @test val(dims(saved1, Y)) ≈ val(dims(geoA, Y))
            @test all(metadata.(dims(saved1)) .== metadata.(dims(geoA)))
            @test metadata(dims(saved1)[1]) == metadata(dims(geoA)[1])
            @test missingval(saved1) === missingval(geoA)
            @test refdims(saved1) == refdims(geoA)
        end

        @testset "3d, with subsetting" begin
            geoA2 = gdalarray[Y(Between(33.7, 33.9)), X(Between(-117.6, -117.4))]
            filename2 = tempname() * ".tif"
            write(filename2, geoA2)
            saved2 = read(GeoArray(filename2; name=:test, mappedcrs=EPSG(4326)))
            @test size(saved2) == size(geoA2) == length.(dims(saved2)) == length.(dims(geoA2))
            @test refdims(saved2) == refdims(geoA2)
            #TODO test a file with more metadata
            @test val(metadata(saved2))[:filepath] == filename2
            @test missingval(saved2) === missingval(geoA2)
            @test GeoData.name(saved2) == GeoData.name(geoA2)
            @test step(mode(dims(saved2, Y))) ≈ step(mode(dims(geoA2, Y)))
            @test step(mode(dims(saved2, X))) ≈ step(mode(dims(geoA2, X)))
            @test typeof(dims(saved2)) == typeof(dims(geoA2))
            @test all(val(dims(saved2, Band)) .≈ val(dims(geoA2, Band)))
            @test all(val(dims(saved2, X)) .≈ val(dims(geoA2, X)))
            @test all(val(dims(saved2, Y)) .≈ val(dims(geoA2, Y)))
            @test all(metadata.(dims(saved2)) .== metadata.(dims(geoA2)))
            @test parent(saved2) == parent(geoA2)
            @test typeof(saved2) == typeof(geoA2)
            filename3 = tempname() * ".tif"
            geoA3 = cat(gdalarray[Band(1)], gdalarray[Band(1)], gdalarray[Band(1)]; dims=Band(1:3))
            write(filename3, geoA3)
            saved3 = read(GeoArray(filename3; mappedcrs=EPSG(4326)))
            @test all(saved3 .== geoA3)
            @test val(dims(saved3, Band)) == 1:3
        end

        @testset "resave current" begin
            filename = tempname() * ".rst"
            write(filename, gdalarray)
            gdalarray2 = GeoArray(filename)
            write(gdalarray2)
            @test read(GeoArray(filename)) == read(gdalarray2)
        end

        @testset "to grd" begin
            write("testgrd.gri", gdalarray)
            grdarray = GeoArray("testgrd.gri")
            @test crs(grdarray) == convert(ProjString, crs(gdalarray))
            @test bounds(grdarray) == (bounds(gdalarray))
            @test val(dims(grdarray, Y)) == reverse(val(dims(gdalarray, Y)))
            @test val(dims(grdarray, X)) ≈ val(dims(gdalarray, X))
            @test grdarray == gdalarray
        end

        @testset "from GeoArray" begin
            filename = tempname() * ".tiff"
            ga = GeoArray(rand(100, 200), (X, Y))
            write(filename, ga)
            saved = GeoArray(filename)
            @test all(reverse(saved[Band(1)]; dims=Y) .=== ga)
            @test saved[1, end, 1] == saved[At(1.0), At(1.0), At(1.0)]
            @test saved[100, 1, 1] == saved[At(100), At(200), At(1)]
            filename2 = tempname() * ".tif"
            ga2 = GeoArray(rand(100, 200), (X(101:200; mode=Sampled()), Y(1:200; mode=Sampled())))
            write(filename2, ga2)
            @test reverse(GeoArray(filename2)[Band(1)]; dims=Y) == ga2
        end

        @testset "to netcdf" begin
            filename2 = tempname() * ".nc"
            write(filename2, gdalarray[Band(1)])
            saved = GeoArray(filename2; crs=crs(gdalarray))
            @test size(saved) == size(gdalarray[Band(1)])
            @test saved ≈ gdalarray[Band(1)]
            clat, clon = DimensionalData.shiftlocus.(Ref(Center()), dims(gdalarray, (Y, X)))
            @test mappedindex(clat) ≈ mappedindex(saved, Y)
            @test mappedindex(clon) ≈ mappedindex(saved, X)
            @test all(mappedbounds(saved, X) .≈ mappedbounds(clon))
            @test all(mappedbounds(saved, Y) .≈ mappedbounds(clat))
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
            A = replace_missing(gdalarray, missing)
            filename = tempname() * ".tif"
            write(filename, A)
            @test missingval(GeoArray(filename)) === typemin(UInt8)
            rm(filename)
        end

    end

    @testset "show" begin
        sh = sprint(show, MIME("text/plain"), gdalarray)
        # Test but don't lock this down too much
        @test occursin("GeoArray", sh)
        @test occursin("Y", sh)
        @test occursin("X", sh)
        @test occursin("Band", sh)
    end

    @testset "plot" begin # TODO write some tests for this
        gdalarray |> plot
        gdalarray[Y(1)] |> plot
    end

    @testset "nodatavalue type matches the array type" begin
        A = GeoArray(WorldClim{Climate}, :tavg; res="10m", month=1)
        @test typeof(missingval(A)) === eltype(A)
        @test missingval(A) === -3.4f38
    end

end

@testset "stack" begin
    @time gdalstack = GeoStack((a=gdalpath, b=gdalpath))

    @test length(gdalstack) == 2
    @test dims(gdalstack) isa Tuple{<:X,<:Y,<:Band}

    @testset "read" begin
        st = read(gdalstack)
        read!((a=gdalpath, b=gdalpath), st)
        @test st isa GeoStack
        @test st.data isa NamedTuple
        @test first(st.data) isa Array
        @test all(st[:a] .=== gdalstack[:a])
    end

    @testset "child array properties" begin
        @test size(gdalstack[:a]) == (514, 515, 1)
        @test gdalstack[:a] isa GeoArray{UInt8,3}
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
            st = replace_missing(gdalstack, mv)
            st = map(A -> (view(A, X(1:100)) .= mv; A), st)
            trimmed = trim(st)
            @test size(trimmed) == (414, 514, 1)
            cropped = crop(st; to=trimmed)
            @test size(cropped) == (414, 514, 1)
            @test map((c, t) -> all(collect(c .=== t)), cropped, trimmed) |> all
            extended = extend(cropped; to=st)
            @test all(collect(extended .== st))
        end
        @testset "mask and mask!" begin
            st = read(gdalstack)
            msk = replace_missing(gdalstack[:a], missing)
            msk[X(1:100), Y([1, 5, 95])] .= missingval(msk)
            @test !any(st[:b][X(1:100)] .=== missingval(msk))
            masked = mask(st; to=msk)
            masked[:b][X(1:100), Y([1, 5, 95])]
            @test all(masked[:b][X(1:100), Y([1, 5, 95])] .=== missing)
            mask!(st; to=msk, missingval=0x00)
            @test all(st[:a][X(1:100), Y([1, 5, 95])] .=== 0x00)
            @test all(st[:b][X(1:100), Y([1, 5, 95])] .=== 0x00)
        end

        @testset "rasterize" begin
            st = map(A -> rebuild(A; missingval=0x00), gdalstack) |> read
            # We round-trip rasterise the Tables.jl form of A
            r_st = rasterize(read(gdalstack); to=st)
            @test all(map((a, b, c) -> all(a .=== b .=== c), st, r_st, gdalstack))
            r_st = rasterize(read(gdalstack); to=st, name=(:a, :b))
            @test all(map((a, b, c) -> all(a .=== b .=== c), st, r_st, gdalstack))
            st = map(A -> rebuild(A .* 0x00; missingval=0x00), gdalstack) |> read
            rasterize!(st, read(gdalstack))
            @test all(map((a, b) -> all(a .=== b), st, gdalstack))
        end
    end

    if VERSION > v"1.1-"
        @testset "copy" begin
            geoA = zero(GeoArray(gdalstack[:a]))
            copy!(geoA, gdalstack, :a)
            # First wrap with GeoArray() here or == loads from disk for each cell.
            # we need a general way of avoiding this in all disk-based sources
            @test all(geoA .== GeoArray(gdalstack[:a]))
        end
    end

    @testset "write" begin
        geoA = read(gdalstack[:a])
        @testset "write multiple files" begin
            filename = tempname() * ".tif"
            write(filename, gdalstack)
            base, ext = splitext(filename)
            filename_b = string(base, "_b", ext)
            saved = read(GeoArray(filename_b))
            @test all(saved .== geoA)
        end

        @testset "write multiple files with custom suffix" begin
            filename = tempname() * ".tif"
            write(filename, gdalstack; suffix=("_first", "_second"))
            base, ext = splitext(filename)
            filename_b = string(base, "_second", ext)
            saved = read(GeoArray(filename_b))
            @test all(saved .== geoA)
        end

        @testset "write netcdf" begin
            filename = tempname() * ".nc"
            write(filename, gdalstack);
            saved = GeoStack(filename);
            @test all(read(saved[:a]) .== geoA)
            rm(filename)
        end
    end

    @testset "show" begin
        sh = sprint(show, MIME("text/plain"), gdalstack)
        # Test but don't lock this down too much
        @test occursin("GeoStack", sh)
        @test occursin("Y", sh)
        @test occursin("X", sh)
        @test occursin("Band", sh)
        @test occursin(":a", sh)
        @test occursin(":b", sh)
    end

end

@testset "series" begin
    gdalser = GeoSeries([gdalpath, gdalpath], (Ti(),); mappedcrs=EPSG(4326), name=:test)
    @test read(gdalser[Ti(1)]) == read(GeoArray(gdalpath; mappedcrs=EPSG(4326), name=:test))
    @test read(gdalser[Ti(1)]) == read(GeoArray(gdalpath; mappedcrs=EPSG(4326), name=:test))

    gdalser = GeoSeries((a=[gdalpath, gdalpath], b=[gdalpath, gdalpath]), (Ti,))
    # Rebuild the ser by wrapping the disk data in Array.
    # `modify` forces `rebuild` on all containers as in-Memory variants
    modified_ser = modify(Array, gdalser)
    @test typeof(modified_ser) <: GeoSeries{<:GeoStack{<:NamedTuple{(:a,:b),<:Tuple{<:Array{UInt8,3},Vararg}}}}

    @testset "read" begin
        ser1 = read(gdalser)
        @test ser1 isa GeoSeries{<:GeoStack}
        @test ser1.data isa Vector{<:GeoStack}
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

    @testset "show" begin
        sh = sprint(show, MIME("text/plain"), gdalser)
        # Test but don't lock this down too much
        @test occursin("GeoSeries", sh)
        @test occursin("GeoStack", sh)
        @test occursin("Ti", sh)
    end
end
