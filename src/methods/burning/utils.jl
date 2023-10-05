_geomindices(geoms) = _geomindices(GI.trait(geoms), geoms) 
_geomindices(::Nothing, geoms) = eachindex(geoms)
_geomindices(::GI.FeatureCollectionTrait, geoms) = 1:GI.nfeature(geoms)
_geomindices(::GI.FeatureTrait, geoms) = _geomindices(GI.geometry(geoms))
_geomindices(::GI.AbstractGeometryTrait, geoms) = 1:GI.ngeom(geoms)

_getgeom(geoms, i::Integer) = _getgeom(GI.trait(geoms), geoms, i)
_getgeom(::GI.FeatureCollectionTrait, geoms, i::Integer) = GI.geometry(GI.getfeature(geoms, i))
_getgeom(::GI.FeatureTrait, geoms, i::Integer) = GI.getgeom(GI.geometry(geoms), i)
_getgeom(::GI.AbstractGeometryTrait, geoms, i::Integer) = GI.getgeom(geom, i)
_getgeom(::GI.PointTrait, geom, i::Integer) = error("PointTrait should not be reached")
_getgeom(::Nothing, geoms, i::Integer) = geoms[i] # Otherwise we can probably just index?

_getgeom(geoms) = _getgeom(GI.trait(geoms), geoms)
_getgeom(::GI.FeatureCollectionTrait, geoms) = (GI.geometry(f) for f in GI.getfeature(geoms))
_getgeom(::GI.FeatureTrait, geoms) = GI.getgeom(GI.geometry(geoms))
_getgeom(::GI.AbstractGeometryTrait, geoms) = GI.getgeom(geoms)
_getgeom(::GI.PointTrait, geom) = error("PointTrait should not be reached")
_getgeom(::Nothing, geoms) = geoms
