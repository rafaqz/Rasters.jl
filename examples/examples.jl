# First time: 
# ] dev http://github.com/rafaqz/CoordinateReferenceSystemsBase.jl
# ] dev http://github.com/rafaqz/DimensionalData.jl
# ] dev http://github.com/rafaqz/GeoData.jl

using GeoData, Test, CoordinateReferenceSystemsBase, 
      NCDatasets, Plots, Statistics, BenchmarkTools

# Load data #######################################################################

# Download an example file (monthly global sea surface temperature)
filename = "tos_O1_2001-2002.nc"
# download("https://www.unidata.ucar.edu/software/netcdf/examples/$filename", filename)
ds = Dataset(filename)

# Convert to a GeoArray 
# This is code that could be added to NCDatasets and similar packages
# so this happens automatically.
ncfields(ds, index) = ds[index][:], ds[index].attrib["units"]
dimz = (Lon(ncfields(ds, "lon")...,), Lat(ncfields(ds, "lat")...,), Time(ncfields(ds, "time")...,))
attrib = ds["tos"].attrib
# The type has some fixed names it want in metadata.
# Have to think about how to inluce all the other arbitrary metadata
meta = Dict(:name=>attrib["long_name"],
            :shortname=>attrib["standard_name"],
            :attrib=>Dict(attrib))  
g = GeoArray(ds["tos"][:,:,:], dimz; units=attrib["units"], metadata=meta); 
close(ds)



# Manipulate data and plot #################################################
pyplot()

# Plot 3d data in a grid
g[Time(1:4)] |> plot

# Or single slices
view(g, Time(4)) |> plot

# Or reduce over some dimension
mean(g; dims=Time) |> plot
std(g; dims=Time) |> plot

# Other things work too
minimum(x->ismissing(x) ? NaN : x, g; dims=Time) |> plot
maximum(x->ismissing(x) ? NaN : x, g; dims=Time) |> plot
reduce(+, g; dims=Time) |> plot

# Or something more complicated: 
# plot the mean sea surface temperature around Australia from August to December 2002
dimranges = Time((DateTime360Day(2002, 08, 1), DateTime360Day(2002, 012, 30))), Lat(-45:0.5), Lon(110:160)
select(g, dimranges) |> x->mean(x; dims=Time) |> plot

# Reorganise the dimensions in the underlying data
# It stil plots the right way up
permutedims(g, [Lat, Lon, Time]) |> plot

# Line plots have useful labels
replace(g[Lat(20)], missing=>NaN) |> surface
g[Lon(170), Time(10)] |> plot
g[Lat(1:80), Lon(170), Time(10)] |> plot


# Save data ################################################################
# needs work
g = replace_missing(g, NaN)
newds = Dataset("out.nc","c")

# Need to add the actual lat and long vars too, but how?
defDim(newds, "lon", size(getdim(g, Lon), 1))
defDim(newds, "lat", size(getdim(g, Lat), 1))

# Define the sea surface temperature variable
v = defVar(newds, shortname(g) , eltype(g) , ("lon", "lat"))
attrib = metadata(g)[:attrib]
attrib["_FillValue"] = Float64(attrib["_FillValue"])
attrib["_missingvalue"] = NaN 
for (key, val) in attrib
    v.attrib[key] = val
end
v[:,:] = parent(g)
close(newds)

# Benchmark #####################################################

dimindex(g) = @inbounds g[Time(20), Lon(17), Lat(80)]
normalindex(g) = @inbounds parent(g)[80, 17, 20]

# gitindex() for a single value has no performance penalty using Lat(80) etc
@btime dimindex($g)
@btime normalindex($g)

# So dimensional indexing is fine to use in spatial simulations.
# Other operations have some penalty, which is unavoidable as we are subsetting
# the dimensions as well as retreiving the data, but can be imroved.

# There are more benchmarks with the DimensionalData.jl tests
