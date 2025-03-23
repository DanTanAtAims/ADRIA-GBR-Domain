using NetCDF, YAXArrays
using CSV, DataFrames
import GeoDataFrames as GDF

using ADRIA
using ProgressMeter

include("common.jl")

if !isdir(CYCLONES_DIR)
    throw(ArgumentError("Could not open cyclones directory. $(CYCLONES_DIR)"))
end

cyclones_files = joinpath.(Ref(CYCLONES_DIR), readdir(CYCLONES_DIR))

first_file = CSV.read(cyclones_files[1], DataFrame, comment="#", header=false)

canonical_gpkg = GDF.read(CANONICAL_GPKG)
canonical_gpkg[!, :RME_GBRMPA_ID] .= replace.(canonical_gpkg.RME_GBRMPA_ID, "-" => "")

n_locations = size(first_file, 1)
n_years     = size(first_file, 2) - 1
n_scenarios = length(cyclones_files)

axlist = (
    Dim{:timesteps}(2001:(2001 + n_years - 1)),
    Dim{:locations}(canonical_gpkg.UNIQUE_ID),
    Dim{:scenarios}(1:n_scenarios)
)

cyclone_categories = YAXArray(
    axlist, zeros(Int64, n_years, n_locations, n_scenarios)
)

function get_cyclone_to_canonical(canonical_gpkg, cyclone_ids)::Vector{Int64}
    cyclone_ids = replace.(cyclone_ids, "-" => "")
    return [
        findfirst(cyclone_ids .== Ref(canon_gbrmpa))
        for canon_gbrmpa in canonical_gpkg.RME_GBRMPA_ID
    ]
end

function get_unique_ids(canonical_gpkg, cyclone_ids)::Vector{String}
    cyclone_ids = replace.(cyclone_ids, "-" => "")
    order = get_cyclone_to_canonical(canonical_gpkg, cyclone_ids)
    return canonical_gpkg.UNIQUE_ID[order]
end

@showprogress desc="Reading cyclone categories." for (scen_idx, cyc_file) in enumerate(cyclones_files)
    cyc_cat = CSV.read(cyc_file, DataFrame, comment="#", header=false)
    cyc_loc_ord = get_cyclone_to_canonical(canonical_gpkg, cyc_cat.Column1)
    cyclone_categories[:, :, scen_idx] .= Matrix(cyc_cat[cyc_loc_ord, 2:end])'
end

function _cyclone_mortality_scens(
    cyclone_scens,
    spatial_data::DataFrame,
    location_ids::Vector{String},
    timeframe::Tuple{Int64,Int64}
)::YAXArray{Float64}
    # Add 1 to every scenarios so they represent indexes in cyclone_mr vectors
    cyclone_scens::YAXArray =
        cyclone_scens[timestep=At(timeframe[1]:timeframe[2])] .+ 1

    species::Vector{String} = String.(ADRIA.functional_group_names())
    cyclone_mortality_scens::YAXArray{Float64} = ADRIA.ZeroDataCube(;
        T=Float64,
        timesteps=timeframe[1]:timeframe[2],
        locations=location_ids,
        species=String.(species),
        scenarios=1:length(cyclone_scens.scenarios)
    )

    # Mortality rate was generated using rrap_dg
    cyclone_mr::NamedTuple = ADRIA._cyclone_mortalities()

    # To filter massives/branchings
    massives::BitVector = contains.(String.(species), ["massives"])
    branchings::BitVector = .!massives

    # Set massives mortality rates
    mr_massives::Vector{Float64} = cyclone_mr[:massives]
    cm_scens_massives = mr_massives[cyclone_scens].data
    for m in species[massives]
        cyclone_mortality_scens[species=At(m)] .= cm_scens_massives
    end

    # Set mortality rates for branching corals at <= 5m depth
    below_5::BitVector = spatial_data.depth_med .<= -5
    if sum(below_5) > 0
        mr_bd5::Vector{Float64} = cyclone_mr[:branching_deeper_than_5]
        cm_scens_bd5::Array{Float64} = mr_bd5[cyclone_scens[location=below_5]].data
        for b in species[branchings]
            cyclone_mortality_scens[locations=below_5, species=At(b)] .= cm_scens_bd5
        end
    end

    # Set mortality rates for branching corals at > 5m depth
    above_5::BitVector = spatial_data.depth_med .> -5
    if sum(above_5) > 0
        mr_bs5::Vector{Float64} = cyclone_mr[:branching_shallower_than_5]
        cm_scens_bs5::Array{Float64} = mr_bs5[cyclone_scens[location=above_5]].data
        for b in species[branchings]
            cyclone_mortality_scens[locations=above_5, species=At(b)] .= cm_scens_bs5
        end
    end

    return cyclone_mortality_scens
end

cyclone_morts = _cyclone_mortality_scens(cyclone_categories, canonical_gpkg, canonical_gpkg.UNIQUE_ID, (2001, 2100))

out_cyclone_dir = joinpath(@__DIR__, "..", "Outputs", "cyclones")
mkpath(out_cyclone_dir)

ds = Dataset(; properties = Dict(), cyclone_mortality = cyclone_morts)

savedataset(ds, path=joinpath(out_cyclone_dir, "cyclone_mortality.nc"), overwrite=true)
