using TTauUtils
using Plots
using Printf

macro name(var)
    string(var)
end

function savepurescreenphotometry(mag_name, screen_name, eclipse_photometry, mag_pars, screen_pars; dir = "eclipses")
    filters = join(keys(eclipse_photometry))

    star_name = mag_pars["star"]
    mkpath("$dir/$star_name/$mag_name/$screen_name")

    mag_path = "$dir/$star_name/$mag_name"
    screen_path = "$mag_path/$screen_name"

    open("$mag_path/$mag_name.pars", "w") do io
        for (n, v) in zip(keys(mag_pars), values(mag_pars))
            if v isa Real
                @printf(io, "%s = %.2f\n", n, v)
            else
                @printf(io, "%s = %s\n", n, v)
            end
        end
    end

    open("$screen_path/$screen_name.pars", "w") do io
        for (n, v) in zip(keys(screen_pars), values(screen_pars))
            @printf(io, "%s = %.2f\n", n, v)
        end
    end

    open("$screen_path/$screen_name.eclipse", "w") do io
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

function saveanomalyphotometry(mag_name, screen_name, anomaly_name, eclipse_photometry, mag_pars, screen_pars, anomaly_pars; dir = "eclipses")
    filters = join(keys(eclipse_photometry))
    
    star_name = mag_pars["star"]
    mkpath("$dir/$star_name/$mag_name/$screen_name")

    mag_path = "$dir/$star_name/$mag_name"
    screen_path = "$mag_path/$screen_name"
    anomaly_path = "$screen_path/$anomaly_name"

    open("$anomaly_path/$anomaly_name.pars", "w") do io
        for (n, v) in zip(keys(anomaly_pars), values(anomaly_pars))
            @printf(io, "%s = %.2f\n", n, v)
        end
    end

    open("$anomaly_path/$anomaly_name.eclipse", "w") do io
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

function loadphotometry(star :: TTauUtils.AbstractStar, args...; kwargs...)
    loadphotometry(star.name, args...; kwargs...)
end

function loadphotometry(star_name, mag_name, screen_name; dir = "eclipses")
    mag_path = "$dir/$star_name/$mag_name"
    screen_path = "$mag_path/$screen_name"
    
    lines = readlines("$mag_path/$mag_name.pars")
    pairs = split.(lines, '=') 
    mag_pars = Dict(strip(p[1]) => p[2] isa Real ? parse(Float64, p[2]) : p[2] for p in pairs)

    lines = readlines("$screen_path/$screen_name.pars")
    pairs = split.(lines, '=') 
    screen_pars = Dict(strip(p[1]) => parse(Float64, p[2]) for p in pairs)

    lines = readlines("$screen_path/$screen_name.eclipse")
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

    return mag_pars, screen_pars, eclipse_photometry
end

star = TTauUtils.Star("RYLupi", 2, 1.5, 5300, 12)
TTauUtils.Stars.savestar(star)
star = TTauUtils.Star("RYLupi")

# T_spot = TTauUtils.Stars.calcmagspottemperature(star, 1e-7, 4, 5)
R_in = 4.0
W = 1.0
Ṁ = 1e-7
incl_angle = 75

star = TTauUtils.Stars.MagnetosphereSpotStarFromMdot(star, Ṁ, R_in, R_in + W)
spot_angle_rad = asin(1/√(R_in + W/2))

orientation = TTauUtils.GeometryAndOrientations.Orientation(incl_angle)
incl_angle_rad = 75/180*π

spot_y = sin(incl_angle_rad - spot_angle_rad)

mag_pars = Dict("star" => star.name,
                "R_in" => R_in,
                "W" => W,
                "lg Ṁ" => log10(Ṁ),
                "i" => incl_angle)

mag_name = "4_1_7_75"

used_filters = "BV"
scatter_intencity = 0.06
scatter_filter = TTauUtils.Eclipses.flat_scatter_filter
RY_Lupi_dispersion = Dict('U' => 1.6, 'B' => 1.3, 'V' => 1.0, 'R' => 0.85, 'I' => 0.7)
# dispersion_filter = TTauUtils.Eclipses.interstellar_dispersion
dispersion_filter = RY_Lupi_dispersion

screen_h = 0.5
screen_τ = 10

screen_pars = Dict(
    "screen_h" => screen_h,
    "screen_τ" => screen_τ,
)

screen_name = "1_10"

screen_no_anomalies = TTauUtils.Eclipses.GaussianScreen(screen_τ, screen_h)

screen = screen_no_anomalies
yscreen = [-20:0.5:-5;-5:0.25:-3;-3:0.01:0;0:0.1:1]
screen_eclipse_photometry = TTauUtils.Eclipses.eclipsephotometry(star, orientation, screen, yscreen, h = 0.01, 
                                            scatter_filter = scatter_filter, dispersion_filter = dispersion_filter,
                                            filters = used_filters, scatter = scatter_intencity)

savepurescreenphotometry(mag_name, screen_name, screen_eclipse_photometry, mag_pars, screen_pars)
mag_pars, screen_pars, screeen_eclipse_photometry = loadphotometry(star, mag_name, screen_name)

plot(screeen_eclipse_photometry['B'] .- screeen_eclipse_photometry['V'], screeen_eclipse_photometry['V'], yflip = true, legend = false)

x_speed_rel_to_y_speed = 0.1
anomaly_screen_position = 1.8*screen_h
anomaly_τ = 9.0
anomaly_x_size = 0.3
anomaly_obliq = 3
anomaly_smoothness = 0.2

anomaly_pars = Dict(
    "x_speed_rel_to_y_speed" => x_speed_rel_to_y_speed,
    "anomaly_screen_position" => anomaly_screen_position,
    "anomaly_τ" => anomaly_τ,
    "anomaly_x_size" => anomaly_x_size,
    "anomaly_obliq" => anomaly_obliq,
    "anomaly_smoothness" => anomaly_smoothness
)

anomaly_name = "test_cloud"

anomaly_y_size = anomaly_x_size/anomaly_obliq 

anomaly = TTauUtils.Eclipses.SmoothBorderScreenAnomaly(0.0, anomaly_screen_position, 
                                                                anomaly_x_size, anomaly_y_size, 
                                                                anomaly_τ, anomaly_smoothness)

anomaly = TTauUtils.Eclipses.XSlitAnomaly(anomaly_screen_position, 0.05, anomaly_τ, 0.5)                                                                

screen_x_start = 0.0
screen_x_end = 0.0

screen_y_start = 0.0
screen_y_end = 0.0

if x_speed_rel_to_y_speed > 1
    screen_x_start = -1 - anomaly_x_size*(1+anomaly_smoothness)
    screen_x_end = -screen_x_start

    screen_y_start = screen_x_start/x_speed_rel_to_y_speed - anomaly_screen_position + spot_y
    screen_y_end = screen_x_end/x_speed_rel_to_y_speed - anomaly_screen_position + spot_y
else
    screen_y_start = -1 - anomaly_y_size*(1+anomaly_smoothness) - anomaly_screen_position
    screen_y_end = -screen_y_start - 2*anomaly_screen_position

    screen_x_start = (screen_y_start + spot_y)*x_speed_rel_to_y_speed
    screen_x_end = (screen_y_end + spot_y)*x_speed_rel_to_y_speed
end


n_phot = 100
xscreen = collect(range(screen_x_start, screen_x_end, n_phot))
yscreen = collect(range(screen_y_start, screen_y_end, n_phot))

screen = TTauUtils.Eclipses.addanomaly(screen_no_anomalies, anomaly)

eclipse_photometry = TTauUtils.Eclipses.eclipsephotometry(star, orientation, screen, yscreen, xscreen, h = 0.01, 
                                            scatter_filter = scatter_filter, dispersion_filter = dispersion_filter,
                                            filters = used_filters, scatter = scatter_intencity)

plot!(eclipse_photometry['B'] .- eclipse_photometry['V'], eclipse_photometry['V'], yflip = true, legend = false)
# # anomaly_τs = [0.0, 10.0]
# # anomaly_screen_position = [1.5]

# begin
# anomaly = TTauUtils.Eclipses.SmoothBorderScreenAnomaly(anomaly_screen_position*x_speed_rel_toy_speed,
#                                                             anomaly_screen_position + spot_y, 
#                                                             anomaly_x_size, anomaly_y_size, 
#                                                             anomaly_τ, anomaly_smoothness)
                           
# screen = TTauUtils.Eclipses.addanomaly(screen_no_anomalies, anomaly)
# # anomaly = TTauUtils.Eclipses.ScreenAnomaly(0, 2, 0.15, 0.05, 10)
# # screen = TTauUtils.Eclipses.addanomaly(screen, anomaly)



# yscreen = [-20:0.5:-5;-5:0.25:-3;-3:0.01:0;0:0.1:1]
# xscreen = x_speed_rel_toy_speed*yscreen
# eclipse_photometry = TTauUtils.Eclipses.eclipsephotometry(star, orientation, screen, yscreen, xscreen, h = 0.01, 
#                                             scatter_filter = scatter_filter, filters = used_filters, scatter = scatter_intencity)

# savephotometry(mag_name, anomaly_name, eclipse_photometry, used_filters, mag_pars, screen_pars, anomaly_pars)
# mag_pars, anomaly_pars, eclipse_photometry = loadphotometry(mag_name, anomaly_name)                                        

function createanim(eclipse_photometry, yscreen, xscreen; dir = "eclipses")
    xs = [-2:0.01:2;]
    ys = [-2:0.01:2;]
    n_x = length(xs)
    n_y = length(ys)
    τs = zeros(n_x, n_y)
    phot = zeros(n_x, n_y)

    for i = 1:n_x, j = 1:n_y
        τs[i,j] = TTauUtils.Eclipses.opticaldepthvis(screen, 0, xs[i], ys[j])
        phot[i,j] = if xs[i]^2 + ys[j]^2 < 1.0
            TTauUtils.Stars.pictureplanecontinuum(star, 3e8/550e-9, xs[i], ys[j], orientation.star_axis)
        else
            TTauUtils.Stars.starcontinuum(star, 3e8/550e-9)*exp(-screen_τ)
        end
    end

    screen_map = heatmap(xs, ys, τs')

    anim = @animate for is = 1:length(yscreen)
        println(is)
        plt_color = plot(screen_eclipse_photometry['B'] .- screen_eclipse_photometry['V'], screen_eclipse_photometry['V'], yflip = true, legend = false)
        plot!(plt_color, eclipse_photometry['B'] .- eclipse_photometry['V'], eclipse_photometry['V'], yflip = true, legend = false)
        scatter!(plt_color, [eclipse_photometry['B'][is] - eclipse_photometry['V'][is]], [eclipse_photometry['V'][is]])
        for i = 1:n_x, j = 1:n_y
            τs[i,j] = TTauUtils.Eclipses.opticaldepthvis(screen, yscreen[is], xs[i] - xscreen[is], ys[j])
            # phot[i,j] = TTauUtils.Stars.pictureplanecontinuum(star, 3e8/550e-9, xs[i], ys[j], orientation.star_axis)
        end
        screen_map = heatmap(xs, ys, τs', clims = (0, 10), aspect_ratio = :equal)
        phot_map = heatmap(xs, ys, (log10.(phot) - log10(ℯ) .* τs)', aspect_ratio = :equal, 
                    clims = (log10(TTauUtils.Stars.starcontinuum(star, 3e8/550e-9)) + log10(exp(-screen_τ)), log10(TTauUtils.Stars.spotcontinuum(star, 3e8/550e-9))))
        plot!(screen_map, cos.([-2π:0.01:2π;]), sin.([-2π:0.01:2π;]))
        plot(plt_color, screen_map, phot_map, layout = @layout([A B C]), size = (1500, 500))
    end

    gif(anim, "eclipse.gif", fps = 5)
end
#     # TTauUtils.Eclipses.manyscreens(star, orientation, screens, ys, outfile = "V695Per_scat.eclipse",filters = "RI", scatter = 0.1, scatter_filter = TTauUtils.Eclipses.vosh_scatter_filter)

# plot(eclipse_photometry['B'] .- eclipse_photometry['V'], eclipse_photometry['V'], yflip = true, legend = false)
createanim(eclipse_photometry, yscreen, xscreen)
# end