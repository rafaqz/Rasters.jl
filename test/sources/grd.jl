using Rasters, Test, Statistics, Dates, Plots
using Rasters.Lookups, Rasters.Dimensions
using DiskArrays
import NCDatasets, ArchGDAL
using Rasters: FileArray, GRDsource, GDALsource, metadata, trim

testpath = joinpath(dirname(pathof(Rasters)), "../test/")
include(joinpath(testpath, "test_utils.jl"))
const DD = DimensionalData

maybedownload("https://raw.githubusercontent.com/rspatial/raster/master/inst/external/rlogo.grd", "rlogo.grd")
maybedownload("https://github.com/rspatial/raster/raw/master/inst/external/rlogo.gri", "rlogo.gri")
stem = joinpath(testpath, "data/rlogo")
@test isfile(stem * ".grd")
@test isfile(stem * ".gri")
grdpath = stem * ".gri"

@testset "Grd array" begin
    @time grdarray = Raster(grdpath)
    @time lazyarray = Raster(grdpath; lazy=true);
    @time eagerarray = Raster(grdpath; lazy=false);
    @test_throws ArgumentError Raster("notafile.grd")
    @test_throws ArgumentError Raster("notafile.gri")

    @testset "lazyness" begin
        # Eager is the default
        @test parent(grdarray) isa Array
        @test parent(lazyarray) isa FileArray
        @test parent(eagerarray) isa Array
    end

    @testset "maskingval keyword" begin
        @time missingarray = Raster(grdpath)
        @test missingval(missingarray) === missing
        @test eltype(missingarray) === Union{Missing,Float32}
        @time missingarray = Raster(grdpath; maskingval=nothing)
        @test missingval(missingarray) === -3.4f38
        @test eltype(missingarray) === Float32
    end

    @testset "open" begin
        @test all(open(A -> A[Y=1], grdarray) .=== grdarray[:, 1, :])
        tempfile = tempname()
        cp(stem * ".grd", tempfile * ".grd")
        cp(stem * ".gri", tempfile * ".gri")
        grdwritearray = Raster(tempfile * ".gri"; lazy=true)
        open(grdwritearray; write=true) do A
            A .*= 2
        end
        @test Raster(tempfile * ".gri") == grdarray .* 2
    end

    @testset "read" begin
        A = read(grdarray)
        @test A isa Raster
        @test parent(A) isa Array
        A2 = zero(A)
        @time read!(grdarray, A2);
        A3 = zero(A)
        @time read!(grdpath, A3);
        @test A == A2 == A3
    end

    @testset "array properties" begin
        @test grdarray isa Raster{Union{Missing,Float32},3}
    end

    @testset "dimensions" begin
        @test length(val(dims(dims(grdarray), X))) == 101
        @test ndims(grdarray) == 3
        @test dims(grdarray) isa Tuple{<:X,<:Y,<:Band}
        @test refdims(grdarray) == ()
        @test bounds(grdarray) == ((0.0, 101.0), (0.0, 77.0), (1, 3))
    end

    @testset "other fields" begin
        proj = ProjString("+proj=merc +datum=WGS84")
        @test name(grdarray) == Symbol("red:green:blue")
        @test missingval(grdarray) === missing
        @test metadata(grdarray) isa Metadata{GRDsource,Dict{String,Any}}
        @test label(grdarray) == "red:green:blue"
        @test units(grdarray) == nothing
        @test crs(grdarray) == proj
        @test mappedcrs(grdarray) == nothing
        @test crs(dims(grdarray, Y)) == proj
        @test crs(dims(grdarray, X)) == proj
        @test mappedcrs(dims(grdarray, Y)) == nothing
        @test mappedcrs(dims(grdarray, X)) == nothing
    end

    @testset "custom keywords" begin
        customgrdarray = Raster(grdpath; 
            name=:test, crs=EPSG(1000), mappedcrs=EPSG(4326), refdims=(Ti(),),
            write=true, lazy=true, dropband=false, replace_missing=true,
        )
        @test name(customgrdarray) == :test
        @test refdims(customgrdarray) == (Ti(),)
        @test label(customgrdarray) == "test"
        @test crs(customgrdarray) == EPSG(1000)
        @test crs(dims(customgrdarray, Y)) == EPSG(1000)
        @test crs(dims(customgrdarray, X)) == EPSG(1000)
        @test mappedcrs(customgrdarray) == EPSG(4326)
        @test mappedcrs(dims(customgrdarray, Y)) == EPSG(4326)
        @test mappedcrs(dims(customgrdarray, X)) == EPSG(4326)
        @test parent(customgrdarray) isa Rasters.FileArray
        @test eltype(customgrdarray) == Union{Float32,Missing}
        # Needs to be separate as it overrides crs/mappedcrs 
        dimsgrdarray = Raster(grdpath; 
            dims=(Z(), X(), Y()),
        )
        @test dims(dimsgrdarray) isa Tuple{<:Z,X,Y}
    end

    @testset "getindex" begin
        @test grdarray[Band(1)] isa Raster{Union{Missing,Float32},2}
        @test grdarray[Y(1), Band(1)] isa Raster{Union{Missing,Float32},1}
        @test grdarray[X(1), Band(1)] isa Raster{Union{Missing,Float32},1}
        @test grdarray[X(50), Y(30), Band(1)] == 115.0f0
        @test grdarray[1, 1, 1] == 255.0f0
        @test grdarray[Y(At(20.0; atol=1e10)), X(At(20; atol=1e10)), Band(3)] == 255.0f0
        @test grdarray[Y(Contains(60)), X(Contains(20)), Band(1)] == 255.0f0
    end

    @testset "methods" begin 
        @test mean(grdarray; dims=Y) == mean(parent(grdarray); dims=2)
        @testset "trim, crop, extend" begin
            a = read(grdarray)
            a[X(1:20)] .= missingval(a)
            trimmed = trim(a)
            @test size(trimmed) == (81, 77, 3)
            cropped = crop(a; to=trimmed)
            @test size(cropped) == (81, 77, 3)
            kwcropped = crop(a; to=trimmed, dims=(X,))
            @test size(kwcropped) == (81, size(a,Y), 3)
            @test all(collect(cropped .=== trimmed))
            extended = extend(cropped; to=a);
            @test all(collect(extended .=== a))
        end

        @testset "mask and mask! to disk" begin
            msk = read(replace_missing(grdarray, missing))
            msk[X(1:73), Y([1, 5, 77])] .= missingval(msk)
            @test !any(grdarray[X(1:73)] .=== missingval(msk))
            masked = mask(grdarray; with=msk)
            @test all(masked[X(1:73), Y([1, 5, 77])] .=== missingval(masked))
            tn = tempname()
            tempgrd = tn * ".grd"
            tempgri = tn * ".gri"
            cp(stem * ".grd", tempgrd)
            cp(stem * ".gri", tempgri)
            @test !all(Raster(tempgrd)[X(1:73), Y([1, 5, 77])] .=== missingval(grdarray))
            open(Raster(tempgrd; lazy=true); write=true) do A
                mask!(A; with=msk, missingval=missingval(A))
            end
            @test all(Raster(tempgri)[X(1:73), Y([1, 5, 77])] .=== missingval(grdarray))
            rm(tempgrd)
            rm(tempgri)
        end

        @testset "classify! to disk" begin
            tn = tempname()
            tempgrd = tn * ".grd"
            tempgri = tn * ".gri"
            cp(stem * ".grd", tempgrd)
            cp(stem * ".gri", tempgri)
            extrema(Raster(tempgrd))
            open(Raster(tempgrd; lazy=true); write=true) do A
                classify!(A, [0.0f0 100.0f0 100.0f0; 100 300 255.0f0])
            end
            A = Raster(tempgrd)
            @test count(==(100.0f0), A) + count(==(255.0f0), A) == length(A)
        end

        @testset "mosaic" begin
            @time grdarray = Raster(grdpath)
            A1 = grdarray[X(1:40), Y(1:30)]
            A2 = grdarray[X(27:80), Y(25:60)]
            tn = tempname()
            tempgrd = tn * ".grd"
            tempgri = tn * ".gri"
            Afile = mosaic(first, A1, A2; missingval=0.0f0, atol=1e-1, filename=tempgrd, maskingval=nothing)
            Amem = mosaic(first, A1, A2; missingval=0.0f0, atol=1e-1)
            Atest = grdarray[X(1:80), Y(1:60)]
            Atest[X(1:26), Y(31:60)] .= 0.0f0
            Atest[X(41:80), Y(1:24)] .= 0.0f0
            @test size(Atest) == size(Afile) == size(Amem)
            @test all(Atest .=== Amem .== Afile)
            read(Atest .- Afile)
        end

        @testset "rasterize" begin
            # A = read(grdarray)
            # R = rasterize(A; to=A)
            # @test all(A .=== R .== grdarray)
            # B = rebuild(read(grdarray) .= 0x00; missingval=0x00)
            # rasterize!(B, read(grdarray))
            # @test all(B .=== grdarray |> collect)
        end

    end

    @testset "selectors" begin
        geoA = grdarray[Y(Contains(3)), X(:), Band(1)]
        @test geoA isa Raster{Union{Missing,Float32},1}
        @test grdarray[X(Contains(20)), Y(Contains(10)), Band(1)] isa Float32
    end

    @testset "conversion to Raster" begin
        geoA = grdarray[X(1:50), Y(1:1), Band(1)]
        @test size(geoA) == (50, 1)
        @test eltype(geoA) <: Union{Missing,Float32}
        @time geoA isa Raster{Union{Missing,Float32},1}
        @test dims(geoA) isa Tuple{<:X,Y}
        @test refdims(geoA) isa Tuple{<:Band}
        @test metadata(geoA) == metadata(grdarray)
        @test missingval(geoA) === missing
        @test name(geoA) == Symbol("red:green:blue")
    end

    @testset "write" begin
        @testset "2d" begin
            filename2 = tempname() * ".gri"
            write(filename2, grdarray[Band(1)]; force=true)
            saved = Raster(filename2)
            # 1 band is added again on save
            @test size(saved) == size(grdarray[Band(1)])
            @test parent(saved) == parent(grdarray[Band(1)])
            write(filename2, grdarray[Band(1)]; force=true, verbose=false)
            @test (@allocations write(filename2, grdarray[Band(1)]; force=true, verbose=false)) < 1e3
        end

        @testset "3d with subset" begin
            geoA = grdarray[1:100, 1:50, 1:2]
            filename = tempname() * ".grd"
            write(filename, GRDsource(), geoA; force=true)
            saved = Raster(filename)
            @test size(saved) == size(geoA)
            @test refdims(saved) == ()
            @test bounds(saved) == bounds(geoA)
            @test size(saved) == size(geoA)
            @test missingval(saved) === missingval(geoA)
            @test metadata(saved) != metadata(geoA)
            @test metadata(saved)["creator"] == "Rasters.jl"
            @test all(metadata.(dims(saved)) .== metadata.(dims(geoA)))
            @test name(saved) == name(geoA)
            @test all(lookup.(dims(saved)) .== lookup.(dims(geoA)))
            @test dims(saved) isa typeof(dims(geoA))
            @test all(val.(dims(saved)) .== val.(dims(geoA)))
            @test all(lookup.(dims(saved)) .== lookup.(dims(geoA)))
            @test all(metadata.(dims(saved)) .== metadata.(dims(geoA)))
            @test dims(saved) == dims(geoA)
            @test all(parent(saved) .=== parent(geoA))
            @test saved isa typeof(geoA)
            @test parent(saved) == parent(geoA)
            @test (@allocations write(filename, GRDsource(), geoA; force = true)) < 1e3
        end

        @testset "to netcdf" begin
            filename2 = tempname() * ".nc"
            span(grdarray[Band(1)])
            write(filename2, grdarray[Band(1)]; force = true)
            saved = Raster(filename2; crs=crs(grdarray))
            @test size(saved) == size(grdarray[Band(1)])
            @test all(replace_missing(saved, missingval(grdarray)) .≈ grdarray[Band(1)])
            @test index(saved, X) ≈ index(grdarray, X) .+ 0.5
            @test index(saved, Y) ≈ index(grdarray, Y) .+ 0.5
            @test bounds(saved, Y) == bounds(grdarray, Y)
            @test bounds(saved, X) == bounds(grdarray, X)
            @test (@allocations write(filename2, grdarray[Band(1)]; force = true)) < 1e3
        end

        @testset "to gdal" begin
            # No Band
            gdalfilename = tempname() * ".tif"
            write(gdalfilename, GDALsource(), grdarray[Band(1)]; force = true)
            @test (@allocations write(gdalfilename, GDALsource(), grdarray[Band(1)]; force = true)) < 1e4
            gdalarray = Raster(gdalfilename; maskingval=nothing)
            # @test convert(ProjString, crs(gdalarray)) == convert(ProjString, EPSG(4326))
            @test val(dims(gdalarray, X)) ≈ val(dims(grdarray, X))
            @test val(dims(gdalarray, Y)) ≈ val(dims(grdarray, Y))
            @test gdalarray ≈ replace_missing(permutedims(grdarray[Band(1)], [X(), Y()]), typemin(Int32))
            # 3 Bands
            gdalfilename2 = tempname() * ".tif"
            write(gdalfilename2, grdarray)
            gdalarray2 = Raster(gdalfilename2)
            @test all(Raster(gdalarray2) .== Raster(grdarray))
            @test val(dims(gdalarray2, Band)) == 1:3
        end

        @testset "write missing" begin
            A = replace_missing(grdarray, missing)
            filename = tempname() * ".grd"
            write(filename, A)
            @test missingval(Raster(filename)) === missing
            filename = tempname() * ".grd"
            write(filename, A)
            @test missingval(Raster(filename; maskingval=nothing)) === typemin(Float32)
        end

    end

    @testset "show" begin
        sh = sprint(show, MIME("text/plain"), grdarray)
        # Test but don't lock this down too much
        @test occursin("Raster", sh)
        @test occursin("Y", sh)
        @test occursin("X", sh)
        @test occursin("Band", sh)
    end

    @testset "plot" begin
        grdarray |> plot
        grdarray[Band(1)] |> plot
    end

end

@testset "Grd stack" begin
    grdstack = RasterStack((a=grdpath, b=grdpath))

    @testset "lazyness" begin
        # Eager is the default
        @test parent(grdstack[:a]) isa Array
        @time lazystack = RasterStack((a=grdpath, b=grdpath); lazy=true);
        @time eagerstack = RasterStack((a=grdpath, b=grdpath); lazy=false);
        @test parent(lazystack[:a]) isa FileArray
        @test parent(eagerstack[:a]) isa Array
    end

    @testset "replace_missing keyword" begin
        st = RasterStack((a=grdpath, b=grdpath); replace_missing=true)
        @test eltype(st) == @NamedTuple{a::Union{Missing,Float32},b::Union{Missing,Float32}}
        @test missingval(st) === missing
    end

    @test length(layers(grdstack)) == 2
    @test dims(grdstack) isa Tuple{<:X,<:Y,<:Band}

    @testset "read" begin
        st = read(grdstack)
        @test st isa RasterStack
        @test parent(st) isa NamedTuple
        @test first(parent(st)) isa Array
    end

    @testset "indexing" begin
        @test grdstack[:a][Y(20), X(20), Band(3)] == 70.0f0
        @test grdstack[:a][Y([2,3]), X(40), Band(2)] == [240.0f0, 246.0f0]
    end

    @testset "child array properties" begin
        @test size(grdstack[:a]) == size(Raster(grdstack[:a])) == (101, 77, 3)
        @test grdstack[:a] isa Raster{Union{Missing,Float32},3}
    end

    # Stack Constructors
    @testset "conversion to RasterStack" begin
        geostack = RasterStack(grdstack)
        @test Symbol.(Tuple(keys(grdstack))) == keys(geostack)
        smallstack = RasterStack(grdstack; name=(:a,))
        @test keys(smallstack) == (:a,)
    end

    if VERSION > v"1.1-"
        @testset "copy" begin
            geoA = zero(Raster(grdstack[:a]))
            copy!(geoA, grdstack, :a)
            # First wrap with Raster() here or == loads from disk for each cell.
            # we need a general way of avoiding this in all disk-based sources
            @test geoA == Raster(grdstack[:a])
        end
    end

    @testset "write" begin
        geoA = Raster(grdstack[:b])
        filename = tempname() * ".grd"
        write(filename, grdstack)
        base, ext = splitext(filename)
        filename_b = string(base, "_b", ext)
        saved = read(Raster(filename_b))
        @test typeof(read(geoA)) == typeof(saved)
        @test parent(saved) == parent(geoA)
    end

    @testset "show" begin
        sh = sprint(show, MIME("text/plain"), grdstack)
        # Test but don't lock this down too much
        @test occursin("RasterStack", sh)
        @test occursin("Y", sh)
        @test occursin("X", sh)
        @test occursin("Band", sh)
        @test occursin(":a", sh)
        @test occursin(":b", sh)
    end

end


@testset "Grd Band stack" begin
    @test keys(RasterStack(grdpath)) == (Symbol("red:green:blue"),)
    grdstack = RasterStack(grdpath; layersfrom=Band)

    @test length(layers(grdstack)) == 3
    @test dims(grdstack) isa Tuple{<:X,<:Y}

    @testset "read" begin
        st = read(grdstack)
        @test st isa RasterStack
        @test parent(st) isa NamedTuple
        @test first(parent(st)) isa Array
    end

    @testset "indexing" begin
        @test grdstack[:Band_3][Y(20), X(20)] == 70.0f0
        @test grdstack[:Band_2][Y([2,3]), X(40)] == [240.0f0, 246.0f0]
    end

    @testset "child array properties" begin
        @test size(grdstack[:Band_3]) == size(Raster(grdstack[:Band_3])) == (101, 77)
        @test grdstack[:Band_1] isa Raster{Union{Missing,Float32},2}
    end

    # Stack Constructors
    @testset "conversion to RasterStack" begin
        geostack = RasterStack(grdstack)
        @test Symbol.(Tuple(keys(grdstack))) == keys(geostack)
        smallstack = RasterStack(grdstack; name=(:Band_2,))
        @test keys(smallstack) == (:Band_2,)
    end

    if VERSION > v"1.1-"
        @testset "copy" begin
            geoA = zero(Raster(grdstack[:Band_3]))
            copy!(geoA, grdstack, :Band_3)
            # First wrap with Raster() here or == loads from disk for each cell.
            # we need a general way of avoiding this in all disk-based sources
            @test geoA == Raster(grdstack[:Band_3])
        end
    end

    @testset "save" begin
        geoA = Raster(grdstack[:Band_3])
        filename = tempname() * ".grd"
        write(filename, grdstack)
        base, ext = splitext(filename)
        filename_3 = string(base, "_Band_3", ext)
        saved = read(Raster(filename_3))
        @test typeof(rebuild(saved[Band(1)], refdims=())) == typeof(read(geoA))
        @test parent(saved[Band(1)]) == parent(geoA)
    end

    @testset "show" begin
        sh = sprint(show, MIME("text/plain"), grdstack)
        # Test but don't lock this down too much
        @test occursin("RasterStack", sh)
        @test occursin("Y", sh)
        @test occursin("X", sh)
        @test occursin(":Band_1", sh)
        @test occursin(":Band_2", sh)
        @test occursin(":Band_3", sh)
    end

end

@testset "Grd series" begin
    grdpath2 = stem * "2" * ".gri"
    write(grdpath2, 2 .* Raster(grdpath); force=true)
    Raster(grdpath) .* 2 == Raster(grdpath2)
    eager_grdseries = RasterSeries([grdpath, grdpath2], (Ti,); mappedcrs=EPSG(4326))
    lazy_grdseries = RasterSeries([grdpath, grdpath2], (Ti,); mappedcrs=EPSG(4326), lazy=true)
    duplicate_first_grdseries = RasterSeries([grdpath, grdpath2], (Ti,); mappedcrs=EPSG(4326), lazy=true, duplicate_first=true)
    @test eager_grdseries[Ti(1)] == lazy_grdseries[Ti(1)] == duplicate_first_grdseries[Ti(1)] == Raster(grdpath)
    @test eager_grdseries[Ti(2)] == lazy_grdseries[Ti(2)] == duplicate_first_grdseries[Ti(2)] == Raster(grdpath2)
    stacks = [RasterStack((a=grdpath, b=grdpath); mappedcrs=EPSG(4326)), RasterStack((a=grdpath2, b=grdpath2); mappedcrs=EPSG(4326))]

    grdseries2 = RasterSeries(stacks, (Ti,))
    @test all(grdseries2[Ti(1)][:a] .== Raster(grdpath; mappedcrs=EPSG(4326), name=:test))
    modified_ser = modify(x -> Array(1.0f0x), grdseries2)
    @test typeof(modified_ser) <: RasterSeries{<:RasterStack{(:a, :b),@NamedTuple{a::Float32,b::Float32},3,@NamedTuple{a::Array{Float32,3},b::Array{Float32,3}}},1}

    @testset "read" begin
        geoseries = read(grdseries2)
        @test geoseries isa RasterSeries{<:RasterStack}
        @test parent(geoseries) isa Vector{<:RasterStack}
        @test parent(geoseries) isa Vector{<:RasterStack}
        @test first(parent(first(parent(geoseries)))) isa Array 
    end
end


nothing
