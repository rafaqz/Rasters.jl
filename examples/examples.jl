# First time:
# ] dev http://github.com/rafaqz/CoordinateReferenceSystemsBase.jl
# ] dev http://github.com/rafaqz/DimensionalData.jl
# ] dev http://github.com/rafaqz/GeoData.jl
# Load data #######################################################################

# Download an example file (monthly global sea surface temperature)
using GeoData, NCDatasets, Statistics, Plots

load(fname) = begin
    ncexamples = "https://www.unidata.ucar.edu/software/netcdf/examples/"
    isfile(fname) || download(joinpath(ncexamples, fname), fname)
    Dataset(filename)
end

# Sea surface temperatures collected by PCMDI for use by the IPCC.
filename = "tos_O1_2001-2002.nc"
dimranges = Time<|Between(DateTime360Day(2002, 08, 1), DateTime360Day(2002, 012, 30)), 
            Lat<|Between(-45, 0.5), 
            Lon<|Between(110, 160)
select(NCstack(filename)["tos"], dimranges...) |> x->mean(x; dims=Time) |> plot
stack = NCstack(filename)
keys(stack)

ds = (Time<|[DateTime360Day(2001, 01, 1), DateTime360Day(2001, 02, 1), DateTime360Day(2001, 03, 1)],)
filenames = ["tos_O1_2001-2002.nc","tos_O1_2001-2002.nc","tos_O1_2001-2002.nc"]
series = GeoSeries(filenames, ds; childtype=NCstack)
series[Time(1)]
array = select(series, Time<|Near<|DateTime360Day(2001, 01, 1))["tos"][Lon(1:100), Time(1)]
array = NCarray(filename)
dims(array)

using CFTime
filename = "TA1cm_0pctShade_1990.nc"
d = Dataset(filename)
data = Float64.(d["time"].var[:])
units = d["time"].attrib["units"]
units = "seconds since 1-1-1970"

tunit_mixedcase,starttime = strip.(split(units," since "))
tunit = lowercase(tunit_mixedcase)
replace(starttime, " " => "-")
t0 = CFTime.parseDT(DT,starttime)
DT = DateTimeStandard
t0, plength = CFTime.timeunits(DT,units)
CFTime.timedecode(DT, data, units)

stack = NCstack(filenames)
stack["TA1cm"]

# Example model output from the ECHAM general circulation model. 
# Almost CF, but not quite. Has a spectral coordinate for variables 
# such as temperature (st) and vorticity (svo). 
filename = "test_echam_spectral.nc"
ds = load(filename)
v = ds["sn"][:,:,:]
reverse(v[Time(1)]; dims=Lat) |> plot # FIXME: Upside down
v |> plot # FIXME: Upside down
ds["xl"][Time(2), Dim{:mlev}(25:33)] |> plot
show(v)
reverse(v; dims=Lon)
show(dims(v))
size(v, Lon)

# Surface data for July 2002 from the ECMWF 40 Years Re-Analysis, daily fields. 
# (ECMWF ERA-40 data used in this study/project have been obtained from 
# the ECMWF data server which has these conditions of use.)
filename = "ECMWF_ERA-40_subset.nc"
ds = load(filename)
# Show every second (daytime) time
ds["msl"][Time(1:2:23)] |> plot

# From the Community Climate System Model (CCSM), one time step of 
# precipitation flux, air temperature, and eastward wind. More details.
filename = "sresa1b_ncar_ccsm3-example.nc"
ds = load(filename)
ds["time_bnds"]
ds["pr"][Time(1)] |> plot

# Initialization data for a CAM 3.0 model run
filename = "cami_0000-09-01_64x128_L26_c030918.nc" 
ds = load(filename)
ds["LCWAT"][Time(1), Dim{:lev}(15)] |> plot

# Smith & Sandwell v. 8.2: 1/30-degree topography and bathymetry
filename = "smith_sandwell_topo_v8_2.nc" 
ds = load(filename)
# ds["ROSE"] |> plot # VERY high res and slow

# NCEP/NCAR Reanalysis 1 data for relative humidity at multiple levels, daily averages for 2003
filename = "rhum.2003.nc" 
ds = load(filename)
v = ds["rhum"]
reverse(v[Time(1)]; dims=Lat) |> plot # FIXME: Upside down. use actual_range
 
# Other conventions

# Sample files following other conventions CDL (metadata) 	netCDF file 	Description
filename = "madis-hydro.nc" # Meteorological Assimilation Data Ingest System (MADIS) hydro data for a single hour within a flood control district. The MADIS data emphasizes quality control and uses AWIPS conventions.
filename = "madis-maritime.nc" # MADIS maritime data for a single hour.
filename = "madis-mesonet.nc" # MADIS mesonet data for a single hour.
filename = "madis-metar.nc" # MADIS METAR data for a single hour.
filename = "madis-profiler.nc" # MADIS profiler data for a single hour.
filename = "madis-raob.nc" # MADIS RAOB data for a single hour.
filename = "madis-sao.nc" # MADIS SAO data for a single hour.
filename = "GLASS.nc" # NCAR Mobile GLASS atmospheric sounding generated from the ASPEN system. This file is an upsonde (ground release). Dropsondes have a similar format.
filename = "WMI_Lear.nc" # Aircraft track files used by Zebra and IDV from BAMEX field project. A flight is broken down into separate hourly files (based on clock time not flight time).
filename = "04091217_ruc.nc" # Output from Rapid Update Cycle model run, uses NUWG conventions.
filename = "GOTEX.C130-example.nc" # An example of sensor data from instruments on an aircraft, stored using NCAR-RAF/nimbus conventions.
filename = "ncswp_SPOL_RHI_.nc"
filename = "ncswp_SPOL_PPI_.nc" # NCAR S-band dual-polarized (S-Pol) radar netCDF sweep file format. Range Height Indicator (RHI) files are vertical scans, Plan Position Indicator (PPI) files are horizontal scans. Each file contains one sweep (azimuth (up/down) for RHI or elevation (horizontal) for PPI).
filename = "HRDL_iop12-example.nc" # HRDL Lidar Data collected during the CASES-99 Experiment near Leon, Kansas. Uses unspecified conventions: more details.
filename = "wrfout_v2_Lambert.nc" # WRF output with staggered axes. Grids that are offset from each other, such as in the WRF model output, are currently difficult to deal with because no commonly used convention expresses the relations between staggered grids.
filename = "slim_100897_198.nc" # SLIMCAT Reference Atmosphere for UTLS-Ozone, see catalog record
filename = "sgpsondewnpnC1.nc" # Atmospheric Radiation Measurement (ARM) sounding file. ARM uses netCDF extensively.
filename = "19981111_0045.nc" # Single-banded satellite image example from the NWS AWIPS system.
filename = "IMAGE0002.nc" # Multi-banded satellite image generated by McIDAS program.

 
# Example netCDF-4 files
filename = "test_hgroups.nc" # A couple of aircraft flight data sets stored as individual groups. For the sake of experimenting, the first group has its dimension (recNum) and the time variable defined at root level. All other groups have all dimensions stored inside the group.
filename = "OMI-Aura_L2-example.nc" # Example of a satellite data file with averaging kernels. Has variables with two dimensions of the same name in it (short O3.COLUMN.PARTIAL_AVK(DATETIME, PRESSURE, PRESSURE)).
filename = "test_echam_spectral-deflated.nc" # Example model output from the ECHAM general circulation model. Same as test_echam_spectral.nc netCDF-classic format file above, but this uses netCDF-4 chunking and compression to store the same data in a smaller file, only 43% as large. Subsets of the data can be accessed efficiently without uncompressing the whole file. 
ds = load(filename)


# Manipulate data and plot #################################################
pyplot()
plotlyjs()

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
select(v, dimranges) |> x->mean(x; dims=Time) |> plot

# Reorganise the dimensions in the underlying data
# It stil plots the right way up
permutedims(g, (Lat, Lon, Time) |> plot

# Line plots have useful labels
replace(g[Lat(20)], missing=>NaN) |> surface
g[Lon(170), Time(10)] |> plot
g[Lat(1:80), Lon(170), Time(10)] |> plot


# Save data ################################################################
# needs work
g = replace_missing(v, NaN)
newds = Dataset("out.nc", "c")

# Need to add the actual lat and long vars too, but how?
defDim(newds, "lon", size(val(getdim(v, Lon())), 1))
defDim(newds, "lat", size(val(getdim(v, Lat())), 1))

# Define the sea surface temperature variable
newvar = defVar(newds, , eltype(v) , ("lon", "lat"))
for (key, val) in newvar.attrib
    newvar.attrib[key] = val
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

using HDF5, GeoData, Dates, Plots
using GeoData: Time

# start_date = Date("2016-01-01")
# end_date = Date("2016-12-31")
# end_date = Date("2016-02-28")
filepath = "/home/raf/CESAR/SMAP/SMAP_L4_SM_gph_20160101T223000_Vv4011_001.h5" # raw climatic data
path = "/home/raf/CESAR/SMAP/" # raw climatic data
series = smapseries(path)
metadata(series)
x = series[Time<|1] ["surface_temp", Lon<|1:1, Lat<|1:5]
dropdims(x; dims=Lon())
@code_warntype series[Time<|1]["surface_temp"][Lon<|1:12, Lat<|1:5]
stack = series[Time(1)]
keys(stack)

dims(series)
refdims(stack)
@time array[Lon(1), Lat(1)]

using GeoData, GDAL, ArchGDAL#, Plots
const AG = ArchGDAL

filepath = "/home/raf/CESAR/Raster/limited_growth/limited_growth_2016_01.tif"
filepath = "/home/raf/CESAR/Raster/human_footprint_cea_project.tif"
array = GDALarray(filepath; missingval=-3.4f38)
plot(array)
dims(array)
# heatmap(4size(array, 2):-4:1, val(dims(array)[1]), parent(array))

nt = (a=filepath, b=filepath1, c=filepath)
stack = GDALstack(nt)
stack[:b] |> plot
split(metadata(stack), "\n")
x = CoordinateReferenceSystemsBase.crs(stack)

ArchGDAL.registerdrivers() do
   ArchGDAL.read(filepath) do ds
       println(ArchGDAL.getband(ds, 1))
   end
end
GeoData.gdalrun(filepath) do ds
    ArchGDAL.read(ds, 1, 1:1, 1:2)
end

#        # println(fieldnames(dataset))
#        # AG.read(dataset, 1)
#        # println(names(dataset))
#        # println(ArchGDAL.height(dataset))
#        # println(ArchGDAL.width(dataset))
#        # drv = ArchGDAL.getdriver(dataset)
#        # println(ArchGDAL.shortname(drv))
#        # println(ArchGDAL.longname(drv))
#        # println(ArchGDAL.nraster(dataset))
#        # geotransform = ArchGDAL.getgeotransform(dataset)
#        # WellKnownText <| string(ArchGDAL.importWKT(AG.getproj(dataset)))
#        # GDAL.getmetadatadomainlist(dataset.ptr)
#hjj
gdal_url = "https://download.osgeo.org/geotiff/samples/GeogToWGS84GeoKey/GeogToWGS84GeoKey5.tif"
download(gdal_url, "GeogToWGS84GeoKey5.tif") 
filepath = geturl(gdal_url)

GeoData.gdalrun("GeogToWGS84GeoKey5.tif") do dataset
    band = ArchGDAL.getband(dataset, 1)
    # GDAL.datasetgetlayercount(ds.ptr)
    # drv = ArchGDAL.getdriver(dataset)
end

using REPL
?→ # Primes the symbols_latex cache
assoc = :right
minprec = 0
maxprec = 20

for x in 1:20000
    op = Symbol(Char(x))
    Base.operator_associativity(op) == assoc || continue
    prec = Base.operator_precedence(op)
    prec > minprec && prec < maxprec || continue 
    print((op, prec))
    cstr = string(op)
    if haskey(REPL.symbols_latex, cstr)
        print(" ", REPL.symbols_latex[cstr])
    end
    println()
end

