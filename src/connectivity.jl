using CSV, DataFrames, ProgressMeter

import GeoDataFrames as GDF

include("common.jl")

canonical_gpkg = GDF.read(CANONICAL_GPKG)

CONNECTIVITY_FILES = readdir(CONNECTIVITY_DIR)

tst_file = CSV.read(
    joinpath(CONNECTIVITY_DIR, CONNECTIVITY_FILES[1]),
    DataFrame,
    header=false,
    comment="#"
)

"""
    name_gbrmpa_id_to_unique_id(name_gbrmpa::String; canonical_gpkg=canonical_gpkg)::String

Convert the id format "Name (GBRMPA ID)" to UNIQUE_ID.
"""
function name_gbrmpa_id_to_unique_id(
    name_gbrmpa::AbstractString;
    canonical_gpkg=canonical_gpkg
)::String
    row_idx = findfirst(x -> x == String(name_gbrmpa), canonical_gpkg.reef_name)
    return canonical_gpkg.UNIQUE_ID[row_idx]
end

function convert_ids(connectivity::DataFrame)::DataFrame
    n_locs::Int64 = size(connectivity, 1) - 1

    for idx in 2:(n_locs + 1)
        connectivity[1, idx] = name_gbrmpa_id_to_unique_id(connectivity[1, idx])
        connectivity[idx, 1] = name_gbrmpa_id_to_unique_id(connectivity[idx, 1])
    end

    return connectivity
end

output_connectivity_dir = joinpath(@__DIR__, "..", "Outputs", "connectivity")
mkpath(output_connectivity_dir)

@showprogress desc="Transforming..." for fn in CONNECTIVITY_FILES
    conn_csv = CSV.read(
        joinpath(CONNECTIVITY_DIR, fn),
        DataFrame,
        header=false,
        comment="#"
    )

    conn_csv = convert_ids(conn_csv)
    CSV.write(
        joinpath(output_connectivity_dir, fn),
        conn_csv,
        writeheader=false
    )
end
