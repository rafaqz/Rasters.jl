using NaturalEarth
using Test

using Rasters, DimensionalData
using Rasters.Lookups
import Proj
import GeometryOps as GO, GeoInterface as GI
using Extents

import NCDatasets

@testset "construction" begin
    # fetch land polygons from Natural Earth
    country_polygons = NaturalEarth.naturalearth("admin_0_countries", 110).geometry
    # create a DimVector from the land polygons
    gl = GeometryLookup(country_polygons)
    @test crs(gl) == EPSG(4326)
    @test all(GO.equals, zip(val(gl), country_polygons))
end

