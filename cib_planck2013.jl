using XGPaint
using Healpix, Plots, Pixell
using HDF5

using Random
Random.seed!(3)

# ## Load halos from HDF5 files, establish a CIB model and cosmology
@time halo_pos, halo_mass = read_halo_catalog_hdf5(
    "/global/cfs/cdirs/sobs/www/users/Radio_WebSky/websky_halos-light.hdf5");

# example , ra and dec in radians, halo mass in M200c (Msun)
# ra, dec, redshift, halo_mass = XGPaint.load_example_halos()
# cosmo = get_cosmology(h=0.6774f0, OmegaM=0.3075f0)
# x, y, z = XGPaint.ra_dec_redshift_to_xyz(ra, dec, redshift, cosmo)
# halo_pos = [x[1:400]'; y[1:400]'; z[1:400]';]
# halo_mass = halo_mass[1:400]

cosmo = get_cosmology(h=0.7f0, OmegaM=0.25f0)
model = CIB_Planck2013{Float32}(nside=4096)

## Write one chunk to disk
function write_freq(
    output_dir, model, cosmo,
    pos, mass, freqs)
    # Allocate some arrays and fill them up for centrals and satellites
    @time sources = generate_sources(model, cosmo, pos, mass);

    # Deposit the sources into maps
    fluxes_cen = Array{Float32, 1}(undef, sources.N_cen)
    fluxes_sat = Array{Float32, 1}(undef, sources.N_sat)
    m_hp = HealpixMap{Float64,RingOrder}(model.nside)

    shape, wcs = fullsky_geometry(deg2rad(0.5 / 60))
    m_car = Enmap(zeros(Float64, shape), wcs)


    h5open(joinpath(output_dir, "sources/centrals.h5"), "w") do file
        write(file, "redshift", sources.redshift_cen)
        write(file, "theta", sources.theta_cen)
        write(file, "phi", sources.phi_cen)
    end
    h5open(joinpath(output_dir, "sources/satellites.h5"), "w") do file
        write(file, "redshift", sources.redshift_sat)
        write(file, "theta", sources.theta_sat)
        write(file, "phi", sources.phi_sat)
    end

    # loop over all frequencies and paint sources to appropriate freq map
    @time begin
        for freq in freqs

            XGPaint.paint!(m_hp, parse(Float32, freq) * 1.0f9, model, sources,
                fluxes_cen, fluxes_sat)
            XGPaint.paint!(m_car, parse(Float32, freq) * 1.0f9, model, sources,
                fluxes_cen, fluxes_sat)

            # save sources with mass, redshift, angles
            h5open(joinpath(output_dir, "sources/centrals_flux_$(freq).h5"), "w") do file
                write(file, "flux", fluxes_cen)
            end
            h5open(joinpath(output_dir, "sources/satellites_flux_$(freq).h5"), "w") do file
                write(file, "flux", fluxes_sat)
            end

            filename_hp = joinpath(output_dir, "cib_$(freq)_hp.fits")
            filename_car = joinpath(output_dir, "cib_$(freq)_car.fits")
            Healpix.saveToFITS(m_hp, "!$(filename_hp)", typechar="D")
            write_map(filename_car, m_car)
        end
    end

end

freqs = [
    "18.7", "21.6", "24.5", "27.3", "30.0", "35.9", "41.7", "44.0", "47.4",
    "63.9", "67.8", "70.0", "73.7", "79.6", "90.2", "100", "111", "129", "143",
    "153", "164", "189", "210", "217", "232", "256", "275", "294", "306", "314",
    "340", "353", "375", "409", "467", "525", "545", "584", "643", "729", "817",
    "857", "906", "994", "1080"
]
# freqs = ["143"]
scratch_dir = joinpath(ENV["SCRATCH"], "cib3")
println("SCRATCH: ", scratch_dir)
mkpath(scratch_dir)
mkpath(joinpath(scratch_dir, "sources"))
write_freq(scratch_dir, model, cosmo, halo_pos, halo_mass, freqs)
