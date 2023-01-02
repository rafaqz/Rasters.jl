"""
    function convert_z_to_p(raster::Raster)
Convert the depth dimension (`Z`) to pressure using `gsw_p_from_z`. Note that pressure
depends on depth and _latitude_ so the returned pressure is a `Matrix` when a `Raster` also
has the latitude dimension. If there is no latitude dimension a reference latitude may be
passed in as the keyword argument `ref_lat` (which has default value of zero i.e. the equator).
"""
function convert_z_to_p(raster::Raster; ref_lat = 0)

    z = dims(raster, Z)
    p = Array{Float64}(undef, length(z))
    if isnothing(dims(raster, Y))

        p = GibbsSeaWater.gsw_p_from_z.(z, ref_lat)

    else

        lats = dims(raster, Y)
        p = Array{Float64}(undef, length(z), length(lats))
        for (i, lat) ∈ enumerate(lats)
            p[:, i] = GibbsSeaWater.gsw_p_from_z.(z, lat)
        end

    end

    return p

end

"""
    function convert_Sₚ_to_Sₐ(raster::Raster)
Convert a `Raster` of practical salinity (`Sₚ`) to absolute salinity (`Sₐ`) using
`gsw_sa_from_sp`.
"""
function convert_Sₚ_to_Sₐ(Sₚ::Raster)

    p = convert_z_to_p(Sₚ)
    lons, lats, z, time = dims(Sₚ)#dims(Sₚ, X), dims(Sₚ, Y), dims(Sₚ, Z), dims(Sₚ, Ti)
    Sₐ = similar(Array(Sₚ), Union{Float64, Missing})
    for t ∈ time
        for (i, lon) ∈ enumerate(lons), (j, lat) ∈ enumerate(lats)

            Sₚ_profile = Sₚ[X(At(lon)), Y(At(lat)), Ti(t)]
            find_nm = findall(.!ismissing.(Sₚ_profile))
            Sₐ[i, j, find_nm, t] = GibbsSeaWater.gsw_sa_from_sp.(Sₚ_profile, p[find_nm, j], lon, lat)

        end
    end

    return Raster(Sₐ, (lons, lats, z, Ti))

end

"""
    function convert_θ_to_Θ(raster::Raster)
Convert a `Raster` of potential temperature (`θ`) to conservative temperature (`Θ`) using
`gsw_ct_from_pt`. This conversion depends on absolute salinity. If the `Raster` has the
variable absolute salinity pass keyword argument `Sₐ_in_raster = true`, by default this set
to `false` assuming that the salinity is practical salinity. The name of the
variables for potential temperature and salinity (either practical or absolute) must be
passed in as a named tuple of the form `(pt = :potential_temp_name, s = :salt_name)` where
`:potential_temp_name` and `:salt_name` are the name of the potential temperature and
salinity in the `Raster`.
"""
function convert_θ_to_Θ(raster::RasterStack, var_names::NamedTuple; Sₐ_in_raster = false)

    θ = raster[var_names.pt]
    S = raster[var_names.s]
    lons, lats, z, time = dims(θ, X), dims(θ, Y), dims(θ, Z), dims(θ, Ti)
    Θ = similar(Array(θ), Union{Float64, Missing})
    if Sₐ_in_raster

        for t ∈ time
            for (i, lon) ∈ enumerate(lons), (j, lat) ∈ enumerate(lats)

                S_profile = S[X(At(lon)), Y(At(lat)), Ti(t)]
                find_nm = findall(.!ismissing.(S_profile))
                θ_profile = θ[X(At(lon)), Y(At(lat)), Z(find_nm), Ti(t)]
                Θ[i, j, find_nm, t] = GibbsSeaWater.gsw_ct_from_pt.(S_profile[find_nm], θ_profile)

            end
        end

    else

        Sₐ = convert_Sₚ_to_Sₐ(S)

        for t ∈ time
            for (i, lon) ∈ enumerate(lons), (j, lat) ∈ enumerate(lats)

                Sₐ_profile = Sₐ[X(At(lon)), Y(At(lat)), Ti(t)]
                find_nm = findall(.!ismissing.(Sₐ_profile))
                θ_profile = θ[X(At(lon)), Y(At(lat)), Z(find_nm), Ti(t)]
                Θ[i, j, find_nm, t] = GibbsSeaWater.gsw_ct_from_pt.(Sₐ_profile[find_nm], θ_profile)

            end
        end

    end

    return Raster(Θ, (lons, lats, z, Ti))

end

"""
    function convert_ocean_vars(raster::RasterStack, var_names::NamedTuple)
    function convert_ocean_vars(raster::Rasterseries, var_names::NamedTuple)
Convert ocean variables depth, practical salinity and potential temperature to pressure,
absolute salinity and conservative temperature. All conversions are done using the julia
implementation of TEOS-10 via [GibbsSeaWater](https://github.com/TEOS-10/GibbsSeaWater.jl).
The potential temperature and practical salinity variables in the original are replaced by
default but can be kept by setting the keyword argument `replace = false`. As pressure
on latitude and depth, it is added as a new variable --- that is, each longitude, latitude,
depth has a variable for pressure.
"""
function convert_ocean_vars(raster::RasterStack, var_names::NamedTuple; replace = true)

    ## call convert_z_to_p and form this into variable that matches dimensions of
    # Sₐ and Θ
    ## call convert_Sₚ_to_Sₐ this returns the `Raster` you want
    ## call convert_θ_to_Θ this returns the `Raster` you want

    ## rebuild the `RasterStack` with the new variables

end
function convert_ocean_vars(raster::Rasterseries, var_names::NamedTuple; replace = true)

    ## call `convert_ocean_vars(::RasterStack)` for each member of the series

end
