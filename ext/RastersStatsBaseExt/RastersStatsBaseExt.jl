module RastersStatsBaseExt

@static if isdefined(Base, :get_extension) # julia < 1.9
    using Rasters, StatsBase
else    
    using ..Rasters, ..StatsBase
end
using StatsBase.Random

const RA = Rasters

import Rasters: _True, _False, _booltype
import Rasters.DimensionalData as DD

include("sample.jl")

end # Module
