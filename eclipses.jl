using TTauUtils

using Plots

star = TTauUtils.Star("RYLupi", 2, 1.5, 5300, 12)
TTauUtils.Stars.savestar(star)
star = TTauUtils.Star("RYLupi")

# T_spot = TTauUtils.Stars.calcmagspottemperature(star, 1e-7, 4, 5)
# star = TTauUtils.Stars.MagnetosphereSpotStarFromMdot(star, 1e-7, 4, 5)
R_in = 4.0
W = 1.0
θ_1 = asin(√(1/(R_in + W)))
θ_2 = asin(√(1/(R_in)))


star = TTauUtils.Stars.MagnetosphereSpotStar(star, 10000, θ_1, θ_2)

angs = [60:5:90;]
spot_position = 

orientation = TTauUtils.GeometryAndOrientations.Orientation(75)

# screens = TTauUtils.Eclipses.GaussianScreen.(10, [0.1, 0.5,1.0,4,10])
screen_no_anomalies = TTauUtils.Eclipses.GaussianScreen(10, 1)
screen = screen_no_anomalies

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
anomaly_dark = TTauUtils.Eclipses.ScreenAnomaly(6, 1.4, 0.3, 0.1, 10)
anomaly_hole = TTauUtils.Eclipses.ScreenAnomaly(0, 0.7, 0.2, 0.07, 0)
screen = TTauUtils.Eclipses.addanomaly(screen_no_anomalies, anomaly_dark)
screen = TTauUtils.Eclipses.addanomaly(screen, anomaly_hole)
# anomaly = TTauUtils.Eclipses.ScreenAnomaly(0, 2, 0.15, 0.05, 10)
# screen = TTauUtils.Eclipses.addanomaly(screen, anomaly)



yscreen = [-20:0.5:-5;-5:0.25:-3;-3:0.01:0;0:0.1:1]
xscreen = 10*yscreen
eclipse_photometry = TTauUtils.Eclipses.eclipsephotometry(star, orientation, screen, yscreen, xscreen, h = 0.01, scatter_filter = TTauUtils.Eclipses.flat_scatter_filter, filters = "BV", scatter = 0.06)

function createanim(eclipse_photometry, yscreen, xscreen)
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

    gif(anim, "eclipse.gif")
end
    # TTauUtils.Eclipses.manyscreens(star, orientation, screens, ys, outfile = "V695Per_scat.eclipse",filters = "RI", scatter = 0.1, scatter_filter = TTauUtils.Eclipses.vosh_scatter_filter)

plot(eclipse_photometry['B'] .- eclipse_photometry['V'], eclipse_photometry['V'], yflip = true, legend = false)
createanim(eclipse_photometry, yscreen, xscreen)
end