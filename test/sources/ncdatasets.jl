using NCDatasets, ArchGDAL, GeoData, Test, Statistics, Dates, CFTime, Plots, GeoFormatTypes
using GeoData: name, window, mode, span, sampling
include(joinpath(dirname(pathof(GeoData)), "../test/test_utils.jl"))

ncexamples = "https://www.unidata.ucar.edu/software/netcdf/examples/"
ncsingle = geturl(joinpath(ncexamples, "tos_O1_2001-2002.nc"))
ncmulti = geturl(joinpath(ncexamples, "test_echam_spectral.nc"))

@testset "NCDarray" begin
    ncarray = NCDarray(ncsingle)

    @testset "array properties" begin
        @test size(ncarray) == (180, 170, 24)
        @test ncarray isa NCDarray
    end

    @testset "dimensions" begin
        @test ndims(ncarray) == 3
        @test length.(val.(dims(ncarray))) == (180, 170, 24)
        @test dims(ncarray) isa Tuple{<:Lon,<:Lat,<:Ti}
        @test refdims(ncarray) == ()
        # TODO detect the time span, and make it Regular
        @test mode.(dims(ncarray)) == 
            (Projected(Ordered(), Regular(2.0), Intervals(Center()), EPSG(4326), nothing),
             Projected(Ordered(), Regular(1.0), Intervals(Center()), EPSG(4326), nothing),
             Sampled(Ordered(), Irregular((DateTime360Day(2001, 1, 16), DateTime360Day(2003, 01, 16))), Intervals(Start())))
        @test bounds(ncarray) == ((0.0, 360.0), (-80.0, 90.0), (DateTime360Day(2001, 1, 16), DateTime360Day(2003, 1, 16)))
    end

    @testset "other fields" begin
        @test ismissing(missingval(ncarray))
        @test metadata(ncarray) isa NCDarrayMetadata # TODO make this a namedtuple
        @test name(ncarray) == "tos"
    end

    @testset "indexing" begin
        @test ncarray[Ti(1)] isa GeoArray{<:Any,2}
        @test ncarray[Lat(1), Ti(1)] isa GeoArray{<:Any,1}
        @test ncarray[Lon(1), Ti(1)] isa GeoArray{<:Any,1}
        @test ncarray[Lon(1), Lat(1), Ti(1)] isa Missing
        @test ncarray[Lon(30), Lat(30), Ti(1)] isa Float32
        @test ncarray[30, 30, 2] === 278.47168f0
    end

    @testset "selectors" begin
        a = ncarray[Lon(At(21.0)), Lat(Between(50, 52)), Ti(Contains(DateTime360Day(2002, 12)))];
        @test bounds(a) == ((50.0, 52.0),)
        x = ncarray[Lon(Contains(150)), Lat(Contains(30)), Ti(1)]
        @test x isa Float32
        # TODO make sure we are getting the right cell.
    end

    @testset "conversion to GeoArray" begin
        geoarray = ncarray[Lon(1:50), Lat(20:20), Ti(1)]
        @test size(geoarray) == (50, 1)
        # @test eltype(geoarray) <: Union{Missing,Float32}
        @time geoarray isa GeoArray{Float32,1}
        @test dims(geoarray) isa Tuple{<:Lon,<:Lat}
        @test refdims(geoarray) isa Tuple{<:Ti}
        @test metadata(geoarray) == metadata(ncarray)
        @test ismissing(missingval(geoarray))
        @test name(geoarray) == "tos"
    end

    @testset "save" begin
        @testset "to netcdf" begin
            # TODO save and load subset
            geoarray = GeoArray(ncarray)
            metadata(geoarray)
            @test size(geoarray) == size(ncarray)
            filename = tempname()
            write(filename, NCDarray, geoarray)
            saved = GeoArray(NCDarray(filename))
            @test size(saved) == size(geoarray)
            @test refdims(saved) == refdims(geoarray)
            @test missingval(saved) === missingval(geoarray)
            @test_broken metadata(saved) == metadata(geoarray)
            @test GeoData.name(saved) == GeoData.name(geoarray)
            @test_broken all(metadata.(dims(saved)) .== metadata.(dims(geoarray)))
            @test all(mode.(dims(saved)) .== mode.(dims(geoarray)))
            @test dims(saved) isa typeof(dims(geoarray))
            @test val(dims(saved)[3]) == val(dims(geoarray)[3])
            @test all(val.(dims(saved)) .== val.(dims(geoarray)))
            @test all(data(saved) .=== data(geoarray))
            @test saved isa typeof(geoarray)
            # TODO test crs
        end
        @testset "to gdal" begin
            nccleaned = replace_missing(ncarray[Ti(1)], -9999.0)
            write("testgdal.tif", GDALarray, nccleaned)
            gdalarray = GDALarray("testgdal.tif")
            # gdalarray WKT is missing one AUTHORITY
            @test_broken crs(gdalarray) == convert(WellKnownText, EPSG(4326))
            # But the Proj representation is the same
            @test convert(ProjString, crs(gdalarray)) == convert(ProjString, EPSG(4326))
            @test val(dims(gdalarray, Lat)) ≈ val(dims(nccleaned, Lat))
            @test val(dims(gdalarray, Lon)) ≈ val(dims(nccleaned, Lon))
            @test GeoArray(gdalarray) ≈ nccleaned
        end
        @testset "to grd" begin
            nccleaned = replace_missing(ncarray[Ti(1)], -9999.0)
            write("testgrd", GrdArray, nccleaned)
            grdarray = GrdArray("testgrd");
            @test crs(grdarray) == convert(ProjString, EPSG(4326))
            @test bounds(grdarray) == (bounds(nccleaned)..., (1, 1))
            @test val(dims(grdarray, Lat)) ≈ val(dims(nccleaned, Lat)) .- 0.5
            @test val(dims(grdarray, Lon)) ≈ val(dims(nccleaned, Lon)) .- 1.0
            @test GeoArray(grdarray) ≈ reverse(nccleaned; dims=Lat)
        end
    end

    @testset "plot" begin
        ncarray |> plot
        ncarray[Ti(1)] |> plot
    end

end

@testset "NCDstack" begin
    ncstack = NCDstack(ncmulti)

    @testset "load ncstack" begin
        @test ncstack isa NCDstack{String}
        @test ismissing(missingval(ncstack))
        @test metadata(ncstack) isa NCDstackMetadata
        @test refdims(ncstack) == ()
        # Loads child as a regular GeoArray
        @test ncstack[:albedo] isa GeoArray{<:Any,3}
        @test ncstack[:albedo, 2, 3, 1] isa Float32
        @test ncstack[:albedo, :, 3, 1] isa GeoArray{<:Any,1}
        @test dims(ncstack, :albedo) isa Tuple{<:Lon,<:Lat,<:Ti}
        @test keys(ncstack) isa NTuple{131,Symbol}
        @test first(keys(ncstack)) == :abso4
        @test metadata(ncstack) isa NCDstackMetadata
        @test metadata(ncstack)["institution"] == "Max-Planck-Institute for Meteorology"
        @test metadata(ncstack, :albedo) isa NCDarrayMetadata
        @test metadata(ncstack, :albedo)["long_name"] == "surface albedo"
        # Test some DimensionalData.jl tools work
        # Time dim should be reduced to length 1 by mean
        @test axes(mean(ncstack[:albedo, Lat(1:20)] , dims=Ti)) ==
              (Base.OneTo(192), Base.OneTo(20), Base.OneTo(1))
        geoarray = ncstack[:albedo][Ti(4:6), Lon(1), Lat(2)]
        @test geoarray == ncstack[:albedo, Ti(4:6), Lon(1), Lat(2)]
        @test size(geoarray) == (3,)
    end

    if VERSION > v"1.1-"
        @testset "copy" begin
            geoarray = ncstack[:albedo]
            copy!(geoarray, ncstack, :albedo);
            # First wrap with GeoArray() here or == loads from disk for each cell.
            # we need a general way of avoiding this in all disk-based sources
            @test geoarray == GeoArray(ncstack[:albedo])
        end
    end

    @testset "indexing" begin
        ncmultistack = NCDstack(ncsingle)
        @test dims(ncmultistack, :tos) isa Tuple{<:Lon,<:Lat,<:Ti}
        @test ncmultistack[:tos] isa GeoArray{<:Any,3}
        @test ncmultistack[:tos, Ti(1)] isa GeoArray{<:Any,2}
        @test ncmultistack[:tos, Lat(1), Ti(1)] isa GeoArray{<:Any,1}
        @test ncmultistack[:tos, 8, 30, 10] isa Float32
    end

    @testset "window" begin
        windowedstack = NCDstack(ncmulti; window=(Lat(1:5), Lon(1:5), Ti(1)))
        @test window(windowedstack) == (Lat(1:5), Lon(1:5), Ti(1))
        windowedarray = windowedstack[:albedo]
        @test size(windowedarray) == (5, 5)
        @test windowedarray[1:3, 2:2] == reshape([0.84936917f0, 0.8776228f0, 0.87498736f0], 3, 1)
        @test windowedarray[1:3, 2] == [0.84936917f0, 0.8776228f0, 0.87498736f0]
        @test windowedarray[1, 2] == 0.84936917f0
        windowedstack = NCDstack(ncmulti; window=(Lat(1:5), Lon(1:5), Ti(1:1)))
        windowedarray = windowedstack[:albedo]
        @test windowedarray[1:3, 2:2, 1:1] == reshape([0.84936917f0, 0.8776228f0, 0.87498736f0], 3, 1, 1)
        @test windowedarray[1:3, 2:2, 1] == reshape([0.84936917f0, 0.8776228f0, 0.87498736f0], 3, 1)
        @test windowedarray[1:3, 2, 1] == [0.84936917f0, 0.8776228f0, 0.87498736f0]
        @test windowedarray[1, 2, 1] == 0.84936917f0
        windowedstack = NCDstack(ncmulti; window=(Ti(1),))
        windowedarray = windowedstack[:albedo]
        @test windowedarray[1:3, 2:2] == reshape([0.84936917f0, 0.8776228f0, 0.87498736f0], 3, 1)
        @test windowedarray[1:3, 2] == [0.84936917f0, 0.8776228f0, 0.87498736f0]
        @test windowedarray[1, 2] ==  0.84936917f0
    end

    @testset "conversion to GeoStack" begin
        geostack = GeoStack(ncstack)
        @test Symbol.(Tuple(keys(geostack))) == keys(ncstack)
        smallstack = GeoStack(ncstack; keys=(:albedo, :evap, :runoff))
        @test keys(smallstack) == (:albedo, :evap, :runoff)
    end

    @testset "save" begin
        geostack = GeoStack(ncstack);
        length(dims(geostack[:aclcac]))
        ndims(geostack[:aclcac])
        filename = tempname()
        write(filename, NCDstack, geostack);
        saved = GeoStack(NCDstack(filename))
        @test keys(saved) == keys(geostack)
        @test metadata(saved)["advection"] == "Lin & Rood"
        @test metadata(saved) == metadata(geostack)
        @test first(values(saved)) == first(values(geostack))
    end

end

@testset "NCD series" begin
    series = GeoSeries([ncmulti, ncmulti], (Ti,);
                       childtype=NCDstack, name="test")
    geoarray = GeoArray(NCDarray(ncmulti, :albedo; name="test"))
    @test series[Ti(1)][:albedo] == geoarray
    @test typeof(series[Ti(1)][:albedo]) == typeof(geoarray)
end

nothing
