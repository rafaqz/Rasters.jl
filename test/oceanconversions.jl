using Rasters, Test
import GibbsSeaWater

lons, lats, z, time = -180:1.0:180, -90:2.0:90, -0:-10.0:-4000, 1:2
rs_dims = (lons, lats, z, time)
N = 500
Sₚ_vals = range(33, 38, length = N)
θ_vals = range(-2, 20, length = N)
Sₚ = rand(Sₚ_vals, X(lons), Y(lats), Z(z), Ti(time))
θ = rand(θ_vals, X(lons), Y(lats), Z(z), Ti(time))
p_ref = 1000.0
test_vars = (Sₚ = Sₚ, θ = θ)
rs_stack = RasterStack(test_vars, (X(lons), Y(lats), Z(z), Ti(time)))
rs_series = RasterSeries([rs_stack[Ti(t)] for t ∈ time], Ti)

## Output to test
rs_stack_res_in_situ = convert_ocean_vars(rs_stack, (sp = :Sₚ, pt = :θ))
rs_stack_res_pd = convert_ocean_vars(rs_stack, (sp = :Sₚ, pt = :θ); p_ref)
rs_series_res_in_situ = convert_ocean_vars(rs_series, (sp = :Sₚ, pt = :θ))
rs_series_res_pd = convert_ocean_vars(rs_series, (sp = :Sₚ, pt = :θ); p_ref)

test_vars = keys(rs_stack_res_in_situ)

## Transform p, Sₚ and θ then find ρ (in-situ and potential) to test against
p = similar(Array(Sₚ))
Sₐ = similar(Array(Sₚ))
Θ = similar(Array(Sₚ))
for t ∈ time
    for (i, lon) ∈ enumerate(lons), (j, lat) ∈ enumerate(lats)
        p[i, j, :, t] .= GibbsSeaWater.gsw_p_from_z.(z, lat)
        Sₐ[i, j, :, t] .= GibbsSeaWater.gsw_sa_from_sp.(Sₚ[i, j, :, t], p[i, j, :, t], lon, lat)
        Θ[i, j, :, t] .= GibbsSeaWater.gsw_ct_from_pt.(Sₐ[i, j, :, t], θ[i, j, :, t])
    end
end

ρ = GibbsSeaWater.gsw_rho.(Sₐ, Θ, p)
ρ_ref = GibbsSeaWater.gsw_rho.(Sₐ, Θ, p_ref)

vars_in_situ = (p, Sₐ, Θ, ρ)
vars_pd = (p, Sₐ, Θ, ρ_ref)

@testset "conversions" begin
    ## In situ density stack tests
    for (i, var) ∈ enumerate(test_vars)
        @test rs_stack_res_in_situ[var] == vars_in_situ[i]
    end

    ## Potential density stack tests
    for (i, var) ∈ enumerate(test_vars)
        @test rs_stack_res_pd[var] == vars_pd[i]
    end

    ## In situ density series tests
    for t ∈ eachindex(rs_series)
        for (i, var) ∈ enumerate(test_vars)
            @test rs_series_res_in_situ[Ti(t)][var] == vars_in_situ[i][:, :, :, t]
        end
    end

    ## Potenrial density series tests
    for t ∈ eachindex(rs_series)
        for (i, var) ∈ enumerate(test_vars)
            @test rs_series_res_pd[Ti(t)][var] == vars_pd[i][:, :, :, t]
        end
    end

end
