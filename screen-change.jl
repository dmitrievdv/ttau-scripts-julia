using TTauUtils
using Plots
using LaTeXStrings

star = Star("RYLupi")

T_star = star.T

R_in = 4.0
W = 2
Ṁ = 1e-7
incl_angle = 80

star = TTauUtils.Stars.MagnetosphereSpotStarFromMdot(star, Ṁ, R_in, R_in + W)

orientation = TTauUtils.GeometryAndOrientations.Orientation(incl_angle)

used_filters = "BV"
scatter_intencity = 0.06
scatter_filter = TTauUtils.Eclipses.flat_scatter_filter
dispersion_filter = TTauUtils.Eclipses.interstellar_dispersion

screen_τ = 10

lg_screen_hs = [-1:0.1:1;]
screen_y_hs = [-5:0.5:-3;-3:0.05:0;]
star_ys = [-1:0.01:1;]

n_star_ys = length(star_ys)

n_hs = length(lg_screen_hs)
n_ys = length(screen_y_hs)

Vs = zeros(n_hs, n_ys + n_star_ys)
Bs = zeros(n_hs, n_ys + n_star_ys)

for i_h = 1:n_hs
    h = 10^lg_screen_hs[i_h]
    screen = TTauUtils.Eclipses.GaussianScreen(screen_τ, h)
    yscreen = screen_y_hs * h .- 1
    yscreen = [yscreen; star_ys]
    println(h)
    eclipse_photometry = TTauUtils.Eclipses.eclipsephotometry(star, orientation, screen, yscreen, h = 0.01, 
                                            scatter_filter = scatter_filter, dispersion_filter = dispersion_filter,
                                            filters = used_filters, scatter = scatter_intencity)
    Vs[i_h,:] .= eclipse_photometry['V']
    Bs[i_h,:] .= eclipse_photometry['B']
end

plt = plot()
for i_h = 1:n_hs
    plot!(plt, Bs[i_h,:] .- Vs[i_h,:], Bs[i_h,:], yflip = true, legend = false)
end
plt

