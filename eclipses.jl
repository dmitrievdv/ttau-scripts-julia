using TTauUtils
using Plots
using Printf

macro name(var)
    string(var)
end

function savephotometry(filename, eclipse_photometry, filters, mag_pars, anomaly_pars; dir = "eclipses")
    mkpath("$dir/$filename")
    open("$dir/$filename/$filename.pars", "w") do io
        for (n, v) in zip(keys(mag_pars), values(mag_pars))
            @printf(io, "%s = %.2f\n", n, v)
        end
        for (n, v) in zip(keys(anomaly_pars), values(anomaly_pars))
            @printf(io, "%s = %.2f\n", n, v)
        end
    end

    open("$dir/$filename/$filename.eclipse", "w") do io
        print(io, "filters:")
        for filter in filters
            print(io, " $filter")
        end
        println(io, "")
        n_phot = length(eclipse_photometry[filters[1]])
        for i = 1:n_phot
            for filter in filters
                @printf(io, "%10.5f", eclipse_photometry[filter][i])
            end
            println(io, "")
        end
    end
end

function loadphotometry(filename; dir = "eclipses")
    lines = readlines("$dir/$filename/$filename.eclipse")
    n_phot = length(lines) - 1
    line = lines[1]
    filters = join(split(line)[2:end])
    eclipse_photometry = Dict(filter => zeros(n_phot) for filter in filters)
    for i = 1:n_phot
        line = lines[i+1]
        phot = Dict(p for p in zip(filters, parse.(Float64, split(line))))
        for filter in filters
            eclipse_photometry[filter][i] = phot[filter]
        end
    end

    lines = readlines("$dir/$filename/$filename.pars")
    pairs = split.(lines, '=') 
    mag_pars = Dict(strip(p[1]) => parse(Float64, p[2]) for p in pairs[1:4])
    anomaly_pars = Dict(strip(p[1]) => parse(Float64, p[2]) for p in pairs[5:end])
    return mag_pars, anomaly_pars, eclipse_photometry
end

star = TTauUtils.Star("RYLupi", 2, 1.5, 5300, 12)
TTauUtils.Stars.savestar(star)
star = TTauUtils.Star("RYLupi")

# T_spot = TTauUtils.Stars.calcmagspottemperature(star, 1e-7, 4, 5)
R_in = 4.0
W = 1.0
Ṁ = 1e-7
incl_angle = 75

mag_pars = Dict("R_in" => R_in,
                "W" => W,
                "lg Ṁ" => log10(Ṁ),
                "i" => incl_angle)

x_speed_rel_toy_speed = 20

anomaly_screen_position = 1.0
anomaly_τ = 0.0

anomaly_x_size = 0.1
anomaly_obliq = 2
anomaly_y_size = anomaly_x_size/anomaly_obliq 
anomaly_smoothness = 0.1

anomaly_pars = Dict(
    "x_speed_rel_toy_speed" => 20,
    "anomaly_screen_position" => 1.0,
    "anomaly_τ" => 0.0,
    "anomaly_x_size" => 0.1,
    "anomaly_obliq" => 2,
    "anomaly_smoothness" => 0.1
)

name = "test"

used_filters = "BV"
scatter_intencity = 0.06
scatter_filter = TTauUtils.Eclipses.flat_scatter_filter

star = TTauUtils.Stars.MagnetosphereSpotStarFromMdot(star, Ṁ, R_in, R_in + W)
spot_angle_rad = asin(1/√(R_in + W/2))

orientation = TTauUtils.GeometryAndOrientations.Orientation(incl_angle)
incl_angle_rad = 75/180*π

# screens = TTauUtils.Eclipses.GaussianScreen.(10, [0.1, 0.5,1.0,4,10])
screen_no_anomalies = TTauUtils.Eclipses.GaussianScreen(10, 1)
screen = screen_no_anomalies
spot_y = sin(incl_angle_rad - spot_angle_rad)

# n_anomalies = 3

# for i = 1:n_anomalies
#     h_y = randn()*0.05 + 0.1
#     while h_y < 0.0
#         h_y = randn()*0.05 + 0.1
#     end
#     h_x = 3*h_y
#     x = randn()*0.5
#     y = randn()*0.3 + 1.5
#     τ_v = 0.1*TTauUtils.Eclipses.opticaldepthvis(screen_no_anomalies, 0, x, y)
#     anomaly = TTauUtils.Eclipses.ScreenAnomaly(x, y, h_x, h_y, 0)
#     screen = TTauUtils.Eclipses.addanomaly(screen, anomaly)
# end

begin
anomaly = TTauUtils.Eclipses.SmoothBorderScreenAnomaly(anomaly_screen_position*x_speed_rel_toy_speed,
                                                            anomaly_screen_position + spot_y, 
                                                            anomaly_x_size, anomaly_y_size, 
                                                            anomaly_τ, anomaly_smoothness)
                           
screen = TTauUtils.Eclipses.addanomaly(screen_no_anomalies, anomaly)
# anomaly = TTauUtils.Eclipses.ScreenAnomaly(0, 2, 0.15, 0.05, 10)
# screen = TTauUtils.Eclipses.addanomaly(screen, anomaly)



yscreen = [-20:0.5:-5;-5:0.25:-3;-3:0.01:0;0:0.1:1]
xscreen = x_speed_rel_toy_speed*yscreen
eclipse_photometry = TTauUtils.Eclipses.eclipsephotometry(star, orientation, screen, yscreen, xscreen, h = 0.01, 
                                            scatter_filter = scatter_filter, filters = used_filters, scatter = scatter_intencity)

savephotometry(name, eclipse_photometry, used_filters, mag_pars, anomaly_pars)
mag_pars, anomaly_pars, eclipse_photometry = loadphotometry("test")                                        

function createanim(filename, eclipse_photometry, yscreen, xscreen; dir = "eclipses")
    xs = [-2:0.01:2;]
    ys = [-2:0.01:2;]
    n_x = length(xs)
    n_y = length(ys)
    τs = zeros(n_x, n_y)
    phot = zeros(n_x, n_y)

    for i = 1:n_x, j = 1:n_y
        τs[i,j] = TTauUtils.Eclipses.opticaldepthvis(screen, -2, xs[i], ys[j])
        phot[i,j] = if xs[i]^2 + ys[j]^2 < 1.0
            TTauUtils.Stars.pictureplanecontinuum(star, 3e8/550e-9, xs[i], ys[j], orientation.star_axis)
        else
            0.0
        end
    end

    screen_map = heatmap(xs, ys, τs')

    anim = @animate for is = 1:length(yscreen)
        println(is)
        plt_color = plot(eclipse_photometry['B'] .- eclipse_photometry['V'], eclipse_photometry['V'], yflip = true, legend = false)
        scatter!(plt_color, [eclipse_photometry['B'][is] - eclipse_photometry['V'][is]], [eclipse_photometry['V'][is]])
        for i = 1:n_x, j = 1:n_y
            τs[i,j] = TTauUtils.Eclipses.opticaldepthvis(screen, yscreen[is], xs[i] - xscreen[is], ys[j])
            # phot[i,j] = TTauUtils.Stars.pictureplanecontinuum(star, 3e8/550e-9, xs[i], ys[j], orientation.star_axis)
        end
        screen_map = heatmap(xs, ys, τs', clims = (0, 10), aspect_ratio = :equal)
        phot_map = heatmap(xs, ys, (exp.(-τs) .* phot)', aspect_ratio = :equal, clims = (0, TTauUtils.Stars.spotcontinuum(star, 3e8/550e-9)))
        plot!(screen_map, cos.([-2π:0.01:2π;]), sin.([-2π:0.01:2π;]))
        plot(plt_color, screen_map, phot_map, layout = @layout([A B C]), size = (1500, 500))
    end

    gif(anim, "$dir/$filename/$filename.gif")
end
    # TTauUtils.Eclipses.manyscreens(star, orientation, screens, ys, outfile = "V695Per_scat.eclipse",filters = "RI", scatter = 0.1, scatter_filter = TTauUtils.Eclipses.vosh_scatter_filter)

plot(eclipse_photometry['B'] .- eclipse_photometry['V'], eclipse_photometry['V'], yflip = true, legend = false)
createanim(name, eclipse_photometry, yscreen, xscreen)
end