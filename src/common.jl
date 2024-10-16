using TOML

CONFIG = TOML.parsefile(joinpath(@__DIR__, "..", "config.toml"))["dependencies"]

CANONICAL_GPKG = CONFIG["canonical_gpkg"]
CONNECTIVITY_DIR = CONFIG["connectivity_dir"]
