using CSV, DataFrames, NetCDF, YAXArrays

import GeoDataFrames as GDF

include("common.jl")

DHW_FILES = String.(readdir(DHW_DIR))

canonical_gpkg = GDF.read(CANONICAL_GPKG)

RCPs = [
    "19",
    "26",
    "45",
    "70",
    "85"
]

output_names = [
    "dhwRCP19.nc",
    "dhwRCP26.nc",
    "dhwRCP45.nc",
    "dhwRCP70.nc",
    "dhwRCP85.nc",
]

RCP_files = [filter(dhw_fn -> contains(dhw_fn, rcp), DHW_FILES) for rcp in RCPs]

"""
    location_permutation(dhw_order::Vector{String}; canonical_gpkg=canonical_gpkg)::Vector{Int64}

Create the permutation required to reorder the dhw locations to match the canonical gpkg.
"""
function location_permutation(
    dhw_order::Vector{String};
    canonical_gpkg=canonical_gpkg
)::Vector{Int64}
    dhw_order[dhw_order .== "20-198"] .= "20198"

    return [findfirst(x -> x == id, dhw_order) for id in canonical_gpkg.RME_GBRMPA_ID]
end

"""
    write_dhw_into_yax!(filepath::String, destination::YAXArray)::Nothing
"""
function write_dhw_into_yax!(
    filepath::String,
    destination::YAXArray,
    scen_idx::Int64;
    canonical_gpkg=canonical_gpkg
)::YAXArray
    source::DataFrame = CSV.read(filepath, DataFrame, header=true)
    loc_perm::Vector{Int} = location_permutation(String.(source.id); canonical_gpkg=canonical_gpkg)
    destination[member=scen_idx] .= Matrix(source[loc_perm, 2:end])'

    return destination
end

"""
    write_UNIQUE_ID_array(fn, canonical_gpkg)::Nothing

Write the UNIQUE_ID netcdf variable directly into the netcdf file. The YAXArray interface
does not allow for YAXArrays with variable element sizes.

Does not currently work.
"""
function write_UNIQUE_ID_array(
    fn, canonical_gpkg
)::Nothing
    properties::Dict{String, Any} = Dict(
        "units"         => "",
        "coordinates"   => "sites",
        "name"          => "UNIQUE_ID",
        "long_name"     => "unique id",
        "standard_name" => "unique_id"
    )
    dims = (
        Dim{:sites}(1:size(canonical_gpkg, 1)),
    )
    NetCDF.ncwrite(canonical_gpkg.UNIQUE_ID, fn, "UNIQUE_ID")
    for (key, val) in properties
        NetCDF.ncputatt(fn, key, val)
    end

    return nothing
end

"""
    lat_lon_array(canonical_gpkg; lat::Bool=true)::YAXArray

Create YAXArray of latitude or longitude values for the dhw netcdf.
"""
function lat_lon_array(canonical_gpkg; lat::Bool=true)::YAXArray
    var_name::String = lat ? "latitude" : "longitude"
    units::String  = lat ? "degrees_north" : "degrees_east"
    properties::Dict{String, Any} = Dict(
        "units"         => units,
        "name"          => var_name,
        "long_name"     => var_name,
        "standard_name" => var_name,
        "projection"    => "EPSG:7844"
    )
    dims = (
        Dim{:sites}(1:size(canonical_gpkg, 1)),
    )

    return YAXArray(dims, lat ? canonical_gpkg.LAT : canonical_gpkg.LON, properties)
end

"""
    create_dhw_dataset(dir::String, fns::Vector{String}; canonical_gpkg=canonical_gpkg)::Dataset
"""
function create_dhw_dataset(
    dir::String,
    fns::Vector{String};
    canonical_gpkg=canonical_gpkg
)::Dataset
    variable_properties::Dict{String, Any}  = Dict(
        "units"         => "DegC-week",
        "name"          => "dhw",
        "long_name"     =>"degree heating week",
        "standard_name" =>"DHW"
    )
    global_properties::Dict{String, Any} = Dict(
        "source_path" => "rme_ml_2024_06_13/data_files/dhw_csv",
        "source_location" => "https://data.mds.gbrrestoration.org/dataset/102.100.100/637903?view=overview",
        "source_desc" => "ReefModEngine v1.0.33 DHW projects"
    )

    timesteps = 2000:2100

    n_timesteps = length(timesteps)
    n_locations = size(canonical_gpkg, 1)
    n_scenarios = length(fns)

    dims = (
        Dim{:timesteps}(timesteps),
        Dim{:sites}(1:n_locations),
        Dim{:member}(1:n_scenarios)
    )

    dhws = YAXArray(dims, zeros(Float64, n_timesteps, n_locations, n_scenarios), variable_properties)
    for (idx, fn) in enumerate(fns)
        dhws = write_dhw_into_yax!(joinpath(dir, fn), dhws, idx; canonical_gpkg=canonical_gpkg)
    end

    lats = lat_lon_array(canonical_gpkg, lat=true)
    lons = lat_lon_array(canonical_gpkg, lat=false)

    return Dataset(;
        properties=global_properties,
        :dhw=>dhws,
        :longitude=>lons,
        :latitude=>lats
    )
end

write_dir = joinpath(@__DIR__, "..", "Outputs", "DHWs")
mkpath(write_dir)

@showprogress desc="Creating DHW NetCDFs" for (write_loc, read_locs) in zip(output_names, RCP_files)
    dhw_dset = create_dhw_dataset(DHW_DIR, read_locs; canonical_gpkg=canonical_gpkg)
    savedataset(dhw_dset, path=joinpath(write_dir, write_loc), driver=:netcdf, overwrite=true)
end
