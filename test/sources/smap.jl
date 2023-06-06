using Rasters, Test, Statistics, Dates, Plots
using Rasters.LookupArrays, Rasters.Dimensions
import ArchGDAL, NCDatasets, HDF5, CFTime
using Rasters: layerkeys, SMAPsource, FileArray

Ext = Base.get_extension(Rasters, :RastersHDF5Ext)

testpath = joinpath(dirname(pathof(Rasters)), "../test/")
include(joinpath(testpath, "test_utils.jl"))

# TODO example files without a login requirement
path1 = joinpath(testpath, "data/SMAP_L4_SM_gph_20160101T223000_Vv4011_001.h5")
path2 = joinpath(testpath, "data/SMAP_L4_SM_gph_20160102T223000_Vv4011_001.h5")

smapkeys = (
    :baseflow_flux, :heat_flux_ground, :heat_flux_latent, :heat_flux_sensible, 
    :height_lowatmmodlay, :land_evapotranspiration_flux, :land_fraction_saturated, 
    :land_fraction_snow_covered, :land_fraction_unsaturated, :land_fraction_wilting, 
    :leaf_area_index, :net_downward_longwave_flux, :net_downward_shortwave_flux, 
    :overland_runoff_flux, :precipitation_total_surface_flux, :radiation_longwave_absorbed_flux, 
    :radiation_shortwave_downward_flux, :sm_profile, :sm_profile_pctl, :sm_profile_wetness, 
    :sm_rootzone, :sm_rootzone_pctl, :sm_rootzone_wetness, :sm_surface, :sm_surface_wetness, 
    :snow_depth, :snow_mass, :snow_melt_flux, :snowfall_surface_flux, :soil_temp_layer1, 
    :soil_temp_layer2, :soil_temp_layer3, :soil_temp_layer4, :soil_temp_layer5, :soil_temp_layer6, 
    :soil_water_infiltration_flux, :specific_humidity_lowatmmodlay, :surface_pressure, :surface_temp, 
    :temp_lowatmmodlay, :vegetation_greenness_fraction, :windspeed_lowatmmodlay
)

if isfile(path1) && isfile(path2)
    # We need to wrap HDF5 for SMAP, as h5 file may not be SMAP files
    @testset "SMAPhdf5 wrapper" begin
        HDF5.h5open(path1) do f
            ds= Ext.SMAPhdf5(f)
            @test keys(ds) == layerkeys(ds) == smapkeys
            @test dims(ds) isa Tuple{<:X,<:Y}
        end
    end

    @testset "Raster" begin
        @time smaparray = Raster(path1)
        @test_throws ArgumentError Raster("notafile.h5")

        @testset "lazyness" begin
            @time read(Raster(path1));
            @time lazyarray = Raster(path1; lazy=true);
            @time eagerarray = Raster(path1; lazy=false);
            # Eager is the default
            @test parent(smaparray) isa Array
            @test parent(lazyarray) isa FileArray
            @test parent(eagerarray) isa Array
        end

        @testset "open" begin
            @test all(open(A -> A[Y=1], smaparray) .=== smaparray[:, 1])
        end

        @testset "read" begin
            @time A = read(smaparray);
            @test A isa Raster
            @test parent(A) isa Array
            A2 = zero(A)
            read!(smaparray, A2)
            A3 = zero(A)
            @time read!(path1, A3);
            @test A == A2 == A3
        end

        @testset "array properties" begin
            @test smaparray isa Raster
        end

        @testset "dimensions" begin
            @test ndims(smaparray) == 2
            HDF5.h5open(path1) do ds
                @test size(smaparray) == length.(dims(smaparray)) == size(ds["Geophysical_Data/baseflow_flux"])
            end
            @test dims(smaparray) isa Tuple{<:X,<:Y}
            @test refdims(smaparray) == ()
            @test sampling(smaparray) == (Intervals(Center()), Intervals(Center()))
            @test span(smaparray) == (Irregular((-180.0f0, 180.0f0)), Irregular((-85.04456f0, 85.04456f0)))
            @test bounds(smaparray) == ((-180.0f0, 180.0f0), (-85.04456f0, 85.04456f0))
        end

        @testset "other fields" begin @test missingval(smaparray) == -9999.0
            @test metadata(smaparray) isa Metadata{SMAPsource,Dict{String,Any}}
            @test name(smaparray) == :baseflow_flux
        end

        @testset "indexing" begin
            @test smaparray[Ti(1)] isa Raster{<:Any,2}
            @test smaparray[Y(1), Ti(1)] isa Raster{<:Any,1}
            @test smaparray[X(1), Ti(1)] isa Raster{<:Any,1}
            @test smaparray[X(1), Y(1)] == -9999.0
            @test smaparray[X(30), Y(30)] isa Float32
            # Russia
            @test smaparray[X(50), Y(100)] == -9999.0 
            # Alaska
            @test smaparray[Y(Near(64.2008)), X(Near(149.4937))] == 0.0f0
        end

        @testset "selectors" begin
            a = smaparray[X(Near(21.0)), Y(Between(50, 52))]
            @test bounds(a) == ((50.08451f0, 51.977905f0),)
            x = smaparray[X(Near(150)), Y(Near(30))]
            @test x isa Float32
            dimz = X(Between(-180.0, 180)), Y(Between(-90, 90)) 
            @test size(smaparray[dimz...]) == (3856, 1624)
            @test index(smaparray[dimz...]) == index(smaparray)
            nca = smaparray[Y(Between(-80, -25)), X(Between(0, 180))]
        end

        @testset "methods" begin 
            @test all(mean(smaparray; dims=Y) .=== mean(parent(smaparray); dims=2))
            @testset "trim, crop, extend" begin
                a = read(smaparray)
                a[X(1:20)] .= missingval(a)
                trimmed = trim(a)
                @test size(trimmed) == (3836, 1512)
                cropped = crop(a; to=trimmed)
                dims(trimmed, X)
                dims(cropped, X)
                @test size(cropped) == (3836, 1512)
                @test all(collect(cropped .=== trimmed))
                dims(cropped)
                extended = extend(cropped; to=a)
                @test all(collect(extended .=== a))
            end
            @testset "slice" begin
                ser = Rasters.slice(smaparray, X) 
                @test ser isa RasterSeries
                @test size(ser) == (3856,)
                @test index(ser, X) == index(smaparray, X)
                @test bounds(ser, X) == (-180.0f0, 180.0f0)
                A = ser[1]
                @test bounds(A, Y) == (-85.04456f0, 85.04456f0)
            end
        end

        @testset "save" begin
            @testset "to grd" begin
                # TODO save and load subset
                geoA = read(smaparray)
                @test size(geoA) == size(smaparray)
                filename = tempname() * ".grd"
                write(filename, geoA)
                saved = read(Raster(filename; mappedcrs=EPSG(4326))[Band(1)])
                dims(saved)
                @test size(saved) == size(geoA)
                @test missingval(saved) === missingval(geoA)
                @test map(metadata.(dims(saved)), metadata.(dims(Raster))) do s, g
                    all(s .== g)
                end |> all
                @test Rasters.name(saved) == Rasters.name(geoA)
                @test all(lookup.(dims(saved)) .!= lookup.(dims(geoA)))
                @test all(
                          mappedbounds(saved)[1]
                          .≈ 
                          bounds(geoA)[1])
                @test all(mappedbounds(saved)[2] .≈ mappedbounds(geoA)[2])
                bounds(saved)
                @test all(parent(saved) .=== parent(geoA))
            end
            @testset "to gdal" begin
                gdalfilename = tempname() * ".tif"
                @time write(gdalfilename, read(smaparray))
                gdalarray = Raster(gdalfilename; mappedcrs=EPSG(4326))
                # These come out with slightly different format
                # @test convert(ProjString, crs(gdalarray)) == crs(smaparray)
                @test all(map((a, b) -> all(a .≈ b), mappedbounds(dims(gdalarray)) , bounds(smaparray)))
                # Tiff locus = Start, SMAP locus = Center
                @test mappedindex(DimensionalData.shiftlocus(Center(), dims(gdalarray, Y))) ≈ index(smaparray, Y)
                @test mappedindex(DimensionalData.shiftlocus(Center(), dims(gdalarray, X))) ≈ index(smaparray, X)
                @test all(gdalarray .== read(smaparray))
            end
            @testset "to netcdf" begin
                ncdfilename = tempname() * ".nc"
                @time write(ncdfilename, smaparray)
                reorder(smaparray, ForwardOrdered())
                saved = Raster(ncdfilename)
                @test_broken bounds(saved) == bounds(smaparray)
                @test index(saved, Y) == index(smaparray, Y)
                @test index(saved, X) == index(smaparray, X)
            end
        end

        @testset "show" begin
            sh = sprint(show, MIME("text/plain"), smaparray)
            @test occursin("Raster", sh)
            @test occursin("Y", sh)
            @test occursin("X", sh)
        end

        @testset "plot" begin
            smaparray[Ti(1:3:12)] |> plot # FIXME plot is weird
            smaparray[Ti(1)] |> plot
            smaparray[Y(100), Ti(1)] |> plot
        end

    end

    @testset "stack" begin
        @time smapstack = RasterStack(path1)

        @testset "lazyness" begin
            @time read(RasterStack(path1));
            @time lazystack = RasterStack(path1; lazy=true)
            @time eagerstack = RasterStack(path1; lazy=false);
            # Lazy is the default
            @test parent(smapstack[:snow_mass]) isa Array
            @test parent(lazystack[:snow_mass]) isa FileArray
            @test parent(eagerstack[:snow_mass]) isa Array
        end

        @testset "read" begin
            @time st = read(smapstack);
            @test st isa RasterStack
            @test parent(st) isa NamedTuple
            @test first(parent(st)) isa Array
            st2 = map(a -> a .* 0, st)
            @time read!(smapstack, st2);
            st3 = map(a -> a .* 0, st)
            @time st = read!(path1, st3);
            @test all(map((a, b, c) -> all(a .== b .== c), st, st2, st3))
        end

        @testset "conversion to Raster" begin
            smaparray = smapstack[:soil_temp_layer1];
            @test smaparray isa Raster{Float32,2}
            @test dims(smaparray) isa Tuple{<:X{<:Mapped{Float32}}, <:Y{<:Mapped{Float32}}}
            @test span(smaparray) isa Tuple{Irregular{Tuple{Float32,Float32}},Irregular{Tuple{Float32,Float32}}}
            @test span(smaparray) == (Irregular((-180.0f0, 180.0f0)), Irregular((-85.04456f0, 85.04456f0)))
            @test bounds(smaparray) == ((-180.0f0, 180.0f0), (-85.04456f0, 85.04456f0))
            @test index(smaparray) isa Tuple{Vector{Float32},Vector{Float32}}
            @test missingval(smaparray) == -9999.0
            @test smaparray[1, 1] == -9999.0
            @test name(smaparray) == :soil_temp_layer1
            dt = DateTime(2016, 1, 1, 22, 30)
            step_ = Hour(3)
            # Currently empty
            @test metadata(smaparray) isa Metadata{SMAPsource,Dict{String,Any}}
        end

        @testset "conversion to regular RasterStack" begin
            # Stack Constructors
            geostack = RasterStack(smapstack)
            @test parent(geostack) isa NamedTuple;
            @test Symbol.(Tuple(keys(smapstack))) == keys(geostack)
            geostack = RasterStack(smapstack; keys=(:baseflow_flux, :snow_mass, :soil_temp_layer1))
            @test keys(geostack) == (:baseflow_flux, :snow_mass, :soil_temp_layer1)
        end

        if VERSION > v"1.1-"
            @testset "copy" begin
                geoA = zero(smapstack[:soil_temp_layer1])
                @test geoA isa Raster
                @test geoA != smapstack[:soil_temp_layer1]
                copy!(geoA, smapstack, :soil_temp_layer1)
                @test geoA == smapstack[:soil_temp_layer1]
            end
        end

        @testset "show" begin
            sh1 = sprint(show, MIME("text/plain"), smapstack[:soil_temp_layer1])
            @test occursin("Raster", sh1)
            @test occursin("Y", sh1)
            @test occursin("X", sh1)
        end

    end

    @testset "series" begin
        ser = smapseries([path1, path2])
        val.(dims(ser))
        @test ser[1] isa RasterStack
        @test first(bounds(ser, Ti)) == DateTime(2016, 1, 1, 22, 30)
        @test last(bounds(ser, Ti)) == DateTime(2016, 1, 3, 1, 30)
        @time modified_series = modify(Array, ser)
        stackkeys = keys(modified_series[1])
        @test typeof(modified_series) <: RasterSeries{<:RasterStack{<:NamedTuple{stackkeys,<:Tuple{<:Array{Float32,2,},Vararg}}}}
        @testset "read" begin
            # FIXME: uses too much memory
            # Seems to be a HDF5.jl memory leak?
            # @time geoseries = read(ser)
            # @test geoseries isa RasterSeries{<:RasterStack}
            # @test geoseries.data isa Vector{<:RasterStack}
            # @test first(geoseries.data[1].data) isa Array 
        end

        @testset "show" begin
            sh = sprint(show, MIME("text/plain"), ser)
            # Test but don't lock this down too much
            @test occursin("RasterSeries", sh)
        end
    end
end
