_geomindices(geoms) = _geomindices(GI.trait(geoms), geoms) 
_geomindices(::Nothing, geoms) = eachindex(geoms)
_geomindices(::GI.FeatureCollectionTrait, geoms) = 1:GI.nfeature(geoms)
_geomindices(::GI.FeatureTrait, geoms) = _geomindices(GI.geometry(geoms))
_geomindices(::GI.AbstractGeometryTrait, geoms) = 1:GI.ngeom(geoms)

_getgeom(geoms, i::Integer) = __getgeom(GI.trait(geoms), geoms, i)
__getgeom(::GI.FeatureCollectionTrait, geoms, i::Integer) = GI.geometry(GI.getfeature(geoms, i))
__getgeom(::GI.FeatureTrait, geoms, i::Integer) = GI.getgeom(GI.geometry(geoms), i)
__getgeom(::GI.AbstractGeometryTrait, geoms, i::Integer) = GI.getgeom(geoms, i)
__getgeom(::GI.PointTrait, geom, i::Integer) = error("PointTrait should not be reached")
__getgeom(::Nothing, geoms, i::Integer) = geoms[i] # Otherwise we can probably just index?

_getgeom(geoms) = __getgeom(GI.trait(geoms), geoms)
__getgeom(::GI.FeatureCollectionTrait, geoms) = (GI.geometry(f) for f in GI.getfeature(geoms))
__getgeom(::GI.FeatureTrait, geoms) = GI.getgeom(GI.geometry(geoms))
__getgeom(::GI.AbstractGeometryTrait, geoms) = GI.getgeom(geoms)
__getgeom(::GI.PointTrait, geom) = error("PointTrait should not be reached")
__getgeom(::Nothing, geoms) = geoms
