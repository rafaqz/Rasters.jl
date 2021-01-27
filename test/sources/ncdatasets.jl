using GeoData, Test, Statistics, Dates, CFTime, Plots, GeoFormatTypes
import ArchGDAL, NCDatasets
using GeoData: name, window, mode, span, sampling, val, Ordered
include(joinpath(dirname(pathof(GeoData)), "../test/test_utils.jl"))

ncexamples = "https://www.unidata.ucar.edu/software/netcdf/examples/"
ncsingle = maybedownload(joinpath(ncexamples, "tos_O1_2001-2002.nc"))
ncmulti = maybedownload(joinpath(ncexamples, "test_echam_spectral.nc"))

stackkeys = (
    :abso4, :aclcac, :aclcov, :ahfcon, :ahfice, :ahfl, :ahfliac, :ahfllac,
    :ahflwac, :ahfres, :ahfs, :ahfsiac, :ahfslac, :ahfswac, :albedo, :albedo_nir,
    :albedo_nir_dif, :albedo_nir_dir, :albedo_vis, :albedo_vis_dif, :albedo_vis_dir,
    :alsobs, :alsoi, :alsol, :alsom, :alsow, :ameltdepth, :ameltfrac, :amlcorac,
    :ao3, :apmeb, :apmegl, :aprc, :aprl, :aprs, :aps, :az0i, :az0l, :az0w,
    :barefrac, :dew2, :drain, :evap, :evapiac, :evaplac, :evapwac, :fage, :friac,
    :geosp, :glac, :gld, :hyai, :hyam, :hybi, :hybm, :lsp, :q, :qres, :qvi, :relhum,
    :runoff, :sd, :seaice, :siced, :sicepdi, :sicepdw, :sicepres, :slm, :sn, :snacl,
    :snc, :sni, :snifrac, :snmel, :sofliac, :sofllac, :soflwac, :srad0, :srad0d,
    :srad0u, :sradl, :srads, :sradsu, :sraf0, :srafl, :srafs, :st, :svo, :t2max,
    :t2min, :temp2, :thvsig, :topmax, :tpot, :trad0, :tradl, :trads, :tradsu,
    :traf0, :trafl, :trafs, :trfliac, :trfllac, :trflwac, :tropo, :tsi, :tsicepdi,
    :tslm1, :tsurf, :tsw, :u10, :ustr, :ustri, :ustrl, :ustrw, :v10, :vdis, :vdisgw,
    :vstr, :vstri, :vstrl, :vstrw, :wimax, :wind10, :wl, :ws, :wsmx, :xi, :xivi,
    :xl, :xlvi
)

@testset "NCDarray" begin
    ncarray = NCDarray(ncsingle)

    @testset "open" begin
        @test all(open(A -> A[Lat=1], ncarray) .=== ncarray[:, 1, :])
    end

    @testset "array properties" begin
        @test size(ncarray) == (180, 170, 24)
        @test ncarray isa NCDarray
        @test val(dims(ncarray, Ti())) == DateTime360Day(2001, 1, 16):Month(1):DateTime360Day(2002, 12, 16)
        @test val(dims(ncarray, Lat())) == -79.5:89.5
        @test val(dims(ncarray, Lon())) == 1.0:2:359.0
        @test bounds(ncarray) == (
            (0.0, 360.0), 
            (-80.0, 90.0), 
            (DateTime360Day(2001, 1, 16), DateTime360Day(2002, 12, 16)),
        )
    end

    @testset "dimensions" begin
        @test ndims(ncarray) == 3
        @test length.(dims(ncarray)) == (180, 170, 24)
        @test dims(ncarray) isa Tuple{<:Lon,<:Lat,<:Ti}
        @test refdims(ncarray) == ()
        # TODO detect the time span, and make it Regular
        @test mode(ncarray) == 
            (Mapped(Ordered(), Regular(2.0), Intervals(Center()), EPSG(4326), EPSG(4326)),
             Mapped(Ordered(), Regular(1.0), Intervals(Center()), EPSG(4326), EPSG(4326)),
             Sampled(Ordered(), Irregular(), Points()))
        @test bounds(ncarray) == ((0.0, 360.0), (-80.0, 90.0), (DateTime360Day(2001, 1, 16), DateTime360Day(2002, 12, 16)))
    end

    @testset "other fields" begin
        @test ismissing(missingval(ncarray))
        @test metadata(ncarray) isa NCDarrayMetadata
        @test name(ncarray) == :tos
    end

    @testset "indexing" begin
        @test ncarray[Ti(1)] isa GeoArray{<:Any,2}
        @test ncarray[Lat(1), Ti(1)] isa GeoArray{<:Any,1}
        @test ncarray[Lon(1), Ti(1)] isa GeoArray{<:Any,1}
        @test ncarray[Lon(1), Lat(1), Ti(1)] isa Missing
        @test ncarray[Lon(30), Lat(30), Ti(1)] isa Float32
        # Russia
        @test ncarray[Lon(50), Lat(100), Ti(1)] isa Missing
        # Alaska
        @test ncarray[Lat(Near(64.2008)), Lon(Near(149.4937)), Ti(1)] isa Missing
        @test ncarray[Ti(2), Lon(At(59.0)), Lat(At(-50.5))] == ncarray[30, 30, 2] === 278.47168f0
    end

    @testset "indexing with reverse lat" begin
        if !haskey(ENV, "CI") # CI downloads fail. But run locally
            ncrevlat = maybedownload("ftp://ftp.cdc.noaa.gov/Datasets/noaa.ersst.v5/sst.mon.ltm.1981-2010.nc")
            ncrevlatarray = NCDstack(ncrevlat; childkwargs=(missingval=-9.96921f36,))[:sst]
            @test order(dims(ncrevlatarray, Lat)) == Ordered(ReverseIndex(), ReverseArray(), ForwardRelation())
            @test ncrevlatarray[Lat(At(40)), Lon(At(100)), Ti(1)] == missingval(ncrevlatarray)
            @test ncrevlatarray[Lat(At(-40)), Lon(At(100)), Ti(1)] == ncrevlatarray[51, 65, 1] == 14.5916605f0
            @test val(span(ncrevlatarray, Ti)) == Month(1)
        end
    end

    @testset "selectors" begin
        a = ncarray[Lon(At(21.0)), Lat(Between(50, 52)), Ti(Near(DateTime360Day(2002, 12)))];
        @test bounds(a) == ((50.0, 52.0),)
        x = ncarray[Lon(Near(150)), Lat(Near(30)), Ti(1)]
        @test x isa Float32
        # TODO make sure we are getting the right cell.
        @test size(ncarray[Lat(Between(-80, 90)), Lon(Between(0, 360)),
            Ti(Between(DateTime360Day(2001, 1, 16), DateTime360Day(2003, 01, 16)))
        ]) == (180, 170, 24)
        nca = ncarray[Lat(Between(-80, -25)), Lon(Between(0, 180)), 
                      Ti(Near(DateTime360Day(2002, 02, 20)))]
        @test size(nca) == (90, 55)
    end

    @testset "conversion to GeoArray" begin
        geoarray = ncarray[Lon(1:50), Lat(20:20), Ti(1)]
        @test size(geoarray) == (50, 1)
        @test eltype(geoarray) <: Union{Missing,Float32}
        @time geoarray isa GeoArray{Float32,1}
        @test dims(geoarray) isa Tuple{<:Lon,<:Lat}
        @test refdims(geoarray) isa Tuple{<:Ti}
        @test metadata(geoarray) == metadata(ncarray)
        @test ismissing(missingval(geoarray))
        @test name(geoarray) == :tos
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
            @test_broken all(metadata.(dims(saved)) .== metadata.(dims(geoarray)))
            @test GeoData.name(saved) == GeoData.name(geoarray)
            @test all(mode.(dims(saved)) .== mode.(dims(geoarray)))
            @test typeof(dims(saved)) <: typeof(dims(geoarray))
            @test val(dims(saved)[3]) == val(dims(geoarray)[3])
            @test all(val.(dims(saved)) .== val.(dims(geoarray)))
            @test all(data(saved) .=== data(geoarray))
            @test saved isa typeof(geoarray)
            # TODO test crs
        end
        @testset "to gdal" begin
            gdalfilename = tempname() * ".tif"
            nccleaned = replace_missing(ncarray[Ti(1)], -9999.0)
            write(gdalfilename, GDALarray, nccleaned)
            gdalarray = GDALarray(gdalfilename)
            # gdalarray WKT is missing one AUTHORITY
            # @test_broken crs(gdalarray) == convert(WellKnownText, EPSG(4326))
            # But the Proj representation is the same
            @test convert(ProjString, crs(gdalarray)) == convert(ProjString, EPSG(4326))
            @test bounds(gdalarray) == (bounds(nccleaned)..., (1, 1))
            # Tiff locus = Start, Netcdf locus = Center
            @test reverse(val(dims(gdalarray, Lat))) .+ 0.5 ≈ val(dims(nccleaned, Lat))
            @test val(dims(gdalarray, Lon)) .+ 1.0  ≈ val(dims(nccleaned, Lon))
            @test reverse(GeoArray(gdalarray); dims=Lat()) ≈ nccleaned
        end
        @testset "to grd" begin
            nccleaned = replace_missing(ncarray[Ti(1)], -9999.0)
            write("testgrd", GRDarray, nccleaned)
            grdarray = GRDarray("testgrd");
            @test crs(grdarray) == convert(ProjString, EPSG(4326))
            @test bounds(grdarray) == (bounds(nccleaned)..., (1, 1))
            @test val(dims(grdarray, Lat)) ≈ val(dims(nccleaned, Lat)) .- 0.5
            @test val(dims(grdarray, Lon)) ≈ val(dims(nccleaned, Lon)) .- 1.0
            @test GeoArray(grdarray) ≈ reverse(nccleaned; dims=Lat)
        end
    end

    @testset "show" begin
        sh = sprint(show, ncarray)
        # Test but don't lock this down too much
        @test occursin("NCDarray", sh)
        @test occursin("Latitude", sh)
        @test occursin("Longitude", sh)
        @test occursin("Time", sh)
    end

    @testset "plot" begin
        ncarray[Ti(1:3:12)] |> plot
        ncarray[Ti(1)] |> plot
        ncarray[Lat(100), Ti(1)] |> plot
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
        @test_throws NCDatasets.NetCDFError ncstack[:not_a_key]
        @test ncstack[:albedo] isa GeoArray{<:Any,3}
        @test ncstack[:albedo, 2, 3, 1] isa Float32
        @test ncstack[:albedo, :, 3, 1] isa GeoArray{<:Any,1}
        @test dims(ncstack, :albedo) isa Tuple{<:Lon,<:Lat,<:Ti}
        @test keys(ncstack) isa NTuple{131,Symbol}
        @test keys(ncstack) == stackkeys
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
        metadata(ncstack)
        geostack = GeoStack(ncstack);
        metadata(geostack)
        length(dims(geostack[:aclcac]))
        ndims(geostack[:aclcac])
        filename = tempname()
        write(filename, NCDstack, geostack);
        saved = GeoStack(NCDstack(filename))
        @test keys(saved) == keys(geostack)
        @test metadata(saved)["advection"] == "Lin & Rood"
        @test metadata(saved) == metadata(geostack) == metadata(ncstack)
        @test all(first(values(saved)) .== first(values(geostack)))
    end

end

@testset "NCD series" begin
    series = GeoSeries([ncmulti, ncmulti], (Ti,); childtype=NCDstack)
    geoarray = GeoArray(NCDarray(ncmulti, :albedo; name=:test))
    @test series[Ti(1)][:albedo] == geoarray
    @test typeof(series[Ti(1)][:albedo]) == typeof(geoarray)
    modified_series = modify(Array, series)
    @test typeof(modified_series) <: GeoSeries{<:GeoStack{<:NamedTuple{stackkeys,<:Tuple{<:GeoArray{Float32,3,<:Tuple,<:Tuple,<:Array{Float32,3}},Vararg}}}}
end

nothing
