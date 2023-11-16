using XGPaint
using Healpix, Plots, Pixell
using HDF5
using Random
Random.seed!(3)

cat_dir = "/global/cfs/projectdirs/sobs/www/users/Radio_WebSky/matched_catalogs/"
map_dir = "/pscratch/sd/x/xzackli/delta_maps/"

shape, wcs = fullsky_geometry(deg2rad(0.5 / 60))
m = Enmap(zeros(shape), wcs)

flux_cut_143::Float32 = 7 * 1e-3  # in Jy
flux_cut_ref_freq = "143.0"       # for computing the mask
freq_freq_filename = joinpath(cat_dir, "catalog_$(flux_cut_ref_freq).h5")
ref_sources = read(h5open(freq_freq_filename, "r"))
flux_cut_mask = ref_sources["flux"] .> flux_cut_143


pixsizes = pixareamap(m)
for (root, dirs, files) in walkdir(cat_dir)
    println("Files in $root")
    for file in files
        freq = replace(file, ".h5" => "", "catalog_" => "")
        filename = joinpath(root, file)

        sources = convert(Dict{String, Vector{Float64}}, read(h5open(filename, "r")))
        fluxes = sources["flux"]
        fluxes[flux_cut_mask] .= 0
        XGPaint.catalog2map!(m, sources["flux"], sources["theta"], sources["phi"], pixsizes)
        map_path = joinpath(map_dir, "flux_cut_7mJy", "map_radio_0.5arcmin_f$(freq).fits")
        mkpath(joinpath(map_dir, "flux_cut_7mJy"))
        write_map("!$(map_path)", m)
    end
end
