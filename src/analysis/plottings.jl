@recipe function f(sol::SciMLBase.AbstractODESolution)
    layout := (3, 3)
    size --> (1000,600)
    thickness_scaling --> 1.0
    legend := false
    @series begin
        subplot := 1
        xguide := "[ms]"
        yguide := "[m]"
        title := raw"$x$"
        sol[5,:].*1e3, sol[1, :]
    end
    @series begin
        subplot := 2
        xguide := "[ms]"
        yguide := "[m]"
        title := raw"$y$"
        sol[5, :].*1e3, sol[3, :]
    end

    @series begin
        subplot := 3
        xguide := "[ms]"
        title := raw"$\delta = \Delta\mathrm{KE}/\mathrm{KE}$"
        sol[5,:].*1e3, sol[6, :]
    end

    @series begin
        subplot := 4
        xguide := "[ms]"
        yguide := "[rad]"
        title := raw"$x'$"
        sol[5,:].*1e3, sol[2, :]
    end
    @series begin
        subplot := 5
        xguide := "[ms]"
        yguide := "[rad]"
        title := raw"$y'$"
        sol[5, :].*1e3, sol[4, :]
    end

    @series begin
        subplot := 6
        xguide := "[ms]"
        yguide := "[m]"
        title := raw"$s$"
        sol[5,:].*1e3, sol.t
    end


    @series begin
        subplot := 7
        xguide := "[ms]"
        yguide := "[rad]"
        title := raw"$S_x$"
        sol[5, :].*1e3, sol[7, :]
    end
    @series begin
        subplot := 8
        xguide := "[ms]"
        yguide := "[rad]"
        title := raw"$S_y$"
        sol[5, :].*1e3, sol[8, :]
    end
    @series begin
        subplot := 9
        xguide := "[ms]"
        yguide := "[rad]"
        title := raw"$S_s$"
        sol[5, :].*1e3, sol[9, :]
    end

end

@userplot VerticalPosition
@recipe function f(input::VerticalPosition; bins=48, cw=true, errorbars=true)
    sol = input.args[1];
    p = sol.prob.p;
    detector_bins =  range(0, length = bins, stop = p.ring.total_length);
    average_positions_per_bin = zeros(bins);
    std_per_bin = zeros(bins);
    if cw == false
        for bin_index in 1:bins
            average_positions_per_bin[bin_index] = mean([sol(detector_bins[bin_index] + p.ring.total_length * (turn - 1))[3]
                    for turn in 1:p.turns])
            std_per_bin[bin_index] = std([sol(detector_bins[bin_index] + p.ring.total_length * (turn - 1))[3]
                    for turn in 1:p.turns])
        end
    else
        for bin_index in 1:bins
            average_positions_per_bin[bin_index] = mean([sol(-detector_bins[bin_index] + p.ring.total_length * (turn))[3]
                    for turn in 1:p.turns])
            std_per_bin[bin_index] = std([sol(-detector_bins[bin_index] + p.ring.total_length * (turn))[3]
                    for turn in 1:p.turns])
        end
    end

    if errorbars
        ribbon := std_per_bin .* 1e6
    end
    label :=  cw ? "CW" : "CCW"
    xguide -->  "Azimuthal position [m]"
    yguide -->  "y [μm]"
    xdata = @view detector_bins[:]

    ydata = @view average_positions_per_bin[:]

    @series begin
        xdata, ydata .* 1e6
    end
end


@userplot HorizontalPosition
@recipe function f(input::HorizontalPosition; bins=48, cw=true, errorbars=true)
    sol = input.args[1];
    p = sol.prob.p;
    detector_bins =  range(0, length = bins, stop = p.ring.total_length);
    average_positions_per_bin = zeros(bins);
    std_per_bin = zeros(bins);
    if cw == false
        for bin_index in 1:bins
            average_positions_per_bin[bin_index] = mean([sol(detector_bins[bin_index] + p.ring.total_length * (turn - 1))[1]
                    for turn in 1:p.turns])
            std_per_bin[bin_index] = std([sol(detector_bins[bin_index] + p.ring.total_length * (turn - 1))[1]
                    for turn in 1:p.turns])
        end
    else
        for bin_index in 1:bins
            average_positions_per_bin[bin_index] = mean([sol(-detector_bins[bin_index] + p.ring.total_length * (turn))[1]
                    for turn in 1:p.turns])
            std_per_bin[bin_index] = std([sol(-detector_bins[bin_index] + p.ring.total_length * (turn))[1]
                    for turn in 1:p.turns])
        end
    end
    label :=  cw ? "CW" : "CCW"
    xguide -->  "Azimuthal position [m]"
    yguide -->  "x [μm]"

    if errorbars
        ribbon := std_per_bin .* 1e6
    end
    label :=  cw ? "CW" : "CCW"
    xguide -->  "Azimuthal position [m]"
    yguide -->  "y [μm]"
    xdata = @view detector_bins[:]

    ydata = @view average_positions_per_bin[:]

    @series begin
        xdata, ydata .* 1e6
    end
end


# function plot_precession(sols; index=8, windows=103, errorbars = false, fit=true)
#     w = windows
#     x  = getMovingAverage(sols[1][time_index,:], w)
#     # errorbars = false
#     y1 = getMovingAverage(sols[1][index,:], w, errorbars=errorbars) # 1 CW

#     y2 = getMovingAverage(sols[2][index,:], w, errorbars=errorbars) # 2 CW

#     y3 = getMovingAverage(sols[3][index,:], w, errorbars=errorbars) # 1 CC

#     y4 = getMovingAverage(sols[4][index,:], w, errorbars=errorbars) # 2 CC

#     y = @. 1/4  * (y1 + y2 - y3 - y4)

#     ylabel = index==8 ? raw"$S_y$ [μrad]" : raw"$S_x$ [μrad]"
#     xlabel = raw"Time [ms]"
#     x .= 1e3 .* x

#     function get_fit_data(x, y)
#         fit_params = curve_fit(lin_model, x, y, zeros(2)).param
#         return lin_model(x, fit_params), fit_params[2]
#     end

#     main = Plots.plot(layout=(3, 3), size = (1000, 600), legend=:bottomleft, legendfontsize=8, titlefontsize=12, dpi=100, link=:all);

#     plot!(main[1, 1], x, 1e6.*y1, title=L"CW_1", ylabel=ylabel, label=false, xformatter=_->"")
#     plot!(main[1, 2], x, 1e6.*y2, title=L"CW_2", label=false, xformatter=_->"", yformatter=_->"")
#     plot!(main[1, 3], x, 1e6.* (y1 .+ y2) ./2, title=L"CW_1 + CW_2", label=false, xformatter=_->"", yformatter=_->"")

#     plot!(main[2, 1], x, 1e6.*y3, title=L"CCW_1", ylabel=ylabel, label=false, xformatter=_->"")
#     plot!(main[2, 2], x, 1e6.*y4, title=L"CCW_2", label=false, xformatter=_->"", yformatter=_->"")
#     plot!(main[2, 3], x, 1e6.* (y3 .+ y4) ./ 2, title=L"CCW_1 + CCW_2", label=false, xformatter=_->"", yformatter=_->"")

#     plot!(main[3, 1], x, 1e6.* (y1 .- y3) ./ 2, title=L"CW_1 - CCW_1", ylabel=ylabel, xlabel=xlabel, label=false)
#     plot!(main[3, 2], x, 1e6.* (y2 .- y4) ./ 2, title=L"CW_2 - CCW_2", xlabel=xlabel, label=false, yformatter=_->"")
#     plot!(main[3, 3], x, 1e6.* y, title="Total", xlabel=xlabel, label=false, yformatter=_->"")

#     if fit
#         for i in eachindex(main.layout.grid)
#             y, slope = get_fit_data(x, main.subplots[i].series_list[1].plotattributes[:y])
#             plot!(main[i], x, y, label=latexify(1e6*slope,fmt=FancyNumberFormatter(3)) *" [nrad/s]" )
#         end
#     end

#     main
# end


@userplot BunchPlot
@recipe function f(input::BunchPlot)
    init_states = input.args[1];

    legend := false
    layout := (1, 3)
    size --> (1000,300)
    thickness_scaling --> 1.0

    @series begin
        seriestype := :scatter
        title := raw"$x$ phase space"
        markersize = 1
        subplot := 1
        xguide := "[mm]"
        yguide := "[mrad]"

        init_states[:,1].*1e3, init_states[:,2].*1e3
    end
    @series begin
        seriestype := :scatter
        subplot := 2
        title := raw"$y$ phase space"
        markersize = 1
        xguide := "[mm]"
        yguide := "[mrad]"

        init_states[:,3].*1e3, init_states[:,4].*1e3
    end
    @series begin
        seriestype := :histogram
        subplot := 3
        xguide := xlab=raw"$10^{-4} \times dp/p$"
        title :="Dispersion"

        1e4 .* get_dp_over_p.(init_states[:,6])
    end

end
