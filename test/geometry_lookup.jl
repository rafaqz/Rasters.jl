using NaturalEarth
using Test

using Rasters, DimensionalData
using Rasters.Lookups
import Proj
import GeometryOps as GO, GeoInterface as GI
using Extents

@testset "construction" begin
    # fetch land polygons from Natural Earth
    country_polygons = NaturalEarth.naturalearth("admin_0_countries", 110).geometry
    # create a DimVector from the land polygons
    gl = GeometryLookup(country_polygons)
    @test crs(gl) == EPSG(4326)
    @test all(GO.equals, zip(val(gl), country_polygons))
end

@testset "reprojecting a GeometryLookup" begin
    # fetch land polygons from Natural Earth
    country_polygons = NaturalEarth.naturalearth("admin_0_countries", 110)
    # create a DimVector from the land polygons
    dv = Raster(rand(Dim{:Geometry}(GeometryLookup(country_polygons))))
    # reproject the full vector data cube (vector data vector, in this case :D)
    target_crs = ProjString("+proj=wintri +type=crs")
    reprojected_via_rasters = val(dims(reproject(target_crs, dv), Dim{:Geometry}))
    reprojected_via_geometryops = GO.reproject(country_polygons; source_crs = EPSG(4326), target_crs = target_crs) 
    # test that the reprojected geometries are the same
    @test all(splat(GO.equals), zip(
        reprojected_via_rasters,  # reproject the vdc, get the geometries from it
        reprojected_via_geometryops # reproject the geometries directly
        )
    )
end

# @testset "Tree types" begin
ras2d = Raster(
    rand(
        X(Sampled(-180:179, sampling = Intervals(Start()))),
        Y(Sampled(-90:89, sampling = Intervals(Start()))),
    );
    crs = EPSG(4326)
)

@testset "Different tree types" begin
    country_polygons = NaturalEarth.naturalearth("admin_0_countries", 110)
    country_lookup = GeometryLookup(country_polygons)
    country_lookup_notree = GeometryLookup(country_polygons; tree = nothing)
    country_lookup_flat_notree = GeometryLookup(country_polygons; tree = GO.SpatialTreeInterface.FlatNoTree)

    results = zonal(sum, ras2d; of = country_polygons, threaded = false)

    results_regular = zonal(sum, ras2d; of = country_lookup)
    results_notree = zonal(sum, ras2d; of = country_lookup_notree)
    results_flat_notree = zonal(sum, ras2d; of = country_lookup_flat_notree)

    @testset "Zonal" begin
        @test results_regular isa Raster
        @test hasdim(results_regular, Geometry)
        @test lookup(results_regular, Geometry) isa GeometryLookup
        @test val(lookup(results_regular, Geometry)) == country_lookup
        @test val(lookup(results_regular, Geometry)) == country_polygons.geometry

        @test results == results_regular
        @test results_regular == results_notree
        @test results_regular == results_flat_notree
    end

    @testset "Selectors" begin
        @testset "Touches that is contained in a geometry" begin
            usa_extent = Extent(X = (103, 104), Y = (44, 45)) # squarely within the USA
            # broken because of GeometryOps?????
            @test_broken only(DimensionalData.dims2indices(results_regular, Geometry(Touches(usa_extent)))) == findfirst(==("USA"), country_polygons.ADM0_A3)
        end
    end
end