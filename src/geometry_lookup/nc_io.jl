using NCDatasets
import NCDatasets.CommonDataModel as CDM

import Rasters
import GeoInterface as GI


import NCDatasets: NetCDF_jll

# generate dataset
output_path = tempname() * ".nc"
testfile = "multipolygons.ncgen"
run(`$(NetCDF_jll.ncgen) -k nc4 -b -o $output_path $(joinpath(dirname(dirname(Base.pathof(Rasters))), "test", "data", testfile))`)
ds = NCDataset(output_path)



var = ds["someData"]

knowndims = Rasters._dims(var)  

unknowndims_idxs = findall(Rasters.isnolookup ∘ Rasters.lookup, knowndims)

if length(unknowndims_idxs) > 1
    throw(ArgumentError("Only one unknown dimension is supported"))
elseif length(unknowndims_idxs) == 0
    return knowndims
end

u_idx = only(unknowndims_idxs)

u_dim_name = CDM.dimnames(var)[u_idx]

has_geometry = haskey(CDM.attribs(var), "geometry")
if !has_geometry
    throw(ArgumentError("Variable $u_dim_name does not have a geometry attribute"))
end
geometry_container_varname = CDM.attribs(var)["geometry"]
geometry_container_var = ds[geometry_container_varname]

geometry_container_attribs = CDM.attribs(geometry_container_var)

haskey(geometry_container_attribs, "geometry_type") && 
geometry_container_attribs["geometry_type"] == "polygon" ||
throw(ArgumentError("We only support polygon geometry types at this time, got $geometry_type"))

geoms = Rasters._geometry_cf_decode(GI.PolygonTrait(), ds, geometry_container_attribs)



# lines now

output_path = tempname() * ".nc"
testfile = "lines.ncgen"
run(`$(NetCDF_jll.ncgen) -k nc4 -b -o $output_path $(joinpath(dirname(dirname(Base.pathof(Rasters))), "test", "data", testfile))`)
ds = NCDataset(output_path)

var = ds["someData"]
geometry_container_varname = CDM.attribs(var)["geometry"]
geometry_container_var = ds[geometry_container_varname]

geometry_container_attribs = CDM.attribs(geometry_container_var)

geoms = Rasters._geometry_cf_decode(GI.LineStringTrait(), ds, geometry_container_attribs)


output_path = tempname() * ".nc"
testfile = "multilines.ncgen"
run(`$(NetCDF_jll.ncgen) -k nc4 -b -o $output_path $(joinpath(dirname(dirname(Base.pathof(Rasters))), "test", "data", testfile))`)
ds = NCDataset(output_path)

var = ds["someData"]
geometry_container_varname = CDM.attribs(var)["geometry"]
geometry_container_var = ds[geometry_container_varname]

geometry_container_attribs = CDM.attribs(geometry_container_var)

geoms = Rasters._geometry_cf_decode(GI.MultiLineStringTrait(), ds, geometry_container_attribs)


output_path = tempname() * ".nc"
testfile = "multipoints.ncgen"
run(`$(NetCDF_jll.ncgen) -k nc4 -b -o $output_path $(joinpath(dirname(dirname(Base.pathof(Rasters))), "test", "data", testfile))`)
ds = NCDataset(output_path)

var = ds["someData"]
geometry_container_varname = CDM.attribs(var)["geometry"]
geometry_container_var = ds[geometry_container_varname]

geometry_container_attribs = CDM.attribs(geometry_container_var)

geoms = Rasters._geometry_cf_decode(GI.MultiPointTrait(), ds, geometry_container_attribs)



output_path = tempname() * ".nc"
testfile = "points.ncgen"
run(`$(NetCDF_jll.ncgen) -k nc4 -b -o $output_path $(joinpath(dirname(dirname(Base.pathof(Rasters))), "test", "data", testfile))`)
ds = NCDataset(output_path)

var = ds["someData"]
geometry_container_varname = CDM.attribs(var)["geometry"]
geometry_container_var = ds[geometry_container_varname]

geometry_container_attribs = CDM.attribs(geometry_container_var)

geoms = Rasters._geometry_cf_decode(GI.PointTrait(), ds, geometry_container_attribs)