include("findpeaks.jl")
# include("emittanceutils.jl")


function almostzeros(i, total)
    r = zeros(total)
    r[i] = 1.0
    return r
end


function get_running_average(x)
    running = zeros(length(x))
    for i in 1:length(x)
        if i == 1
            running[1] = x[1]
        else
            running[i] = (running[i-1] * (i-1) + x[i])/i
        end
    end
    return running
end


function get_moving_average(x, n::Int64; errorbars=false, shift=:none)
    N = length(x)
    window_size = floor(Int, N/n)
    shortened_N = window_size  * n
    shortened_x = @view x[1:shortened_N]

    if errorbars
        output = Vector{Measurement{Float64}}(undef, n)
    else
        output = zeros(n)
    end

    for i in 1:n
        index_range = (i-1)*window_size + 1 : i * window_size
        output_i = mean(shortened_x[index_range])
        error_i = std(shortened_x[index_range]) / sqrt(window_size)
        if errorbars
            output[i] = measurement(output_i, error_i)
        else
            output[i] = output_i
        end

    end

    if shift == :mean
        return output .- mean(output)
    end

    if shift == :first
        return output .- output[1]
    end
    return output
end

function get_moving_average(x, n::Bool; errorbars=false)
    return x
end

function get_moving_average(t, y, n; errorbars=false)
    mt = get_moving_average(t, n, errorbars=false)
    my = get_moving_average(y, n, errorbars=errorbars)

    return mt, my
end

function get_max_min_excursions(s, positions, bin_size, total_length)
    bin_bounds = 0:bin_size:total_length
    number_of_bins = length(bin_bounds) - 1
    maxXs = zeros(number_of_bins)
    minXs = zeros(number_of_bins)
    for i=1:length(positions)
        current_pos = positions[i]
        current_s = s[i] % total_length
        histogram_index = ceil(Int, current_s/bin_size)
        histogram_index == 0 ? histogram_index=1 : histogram_index
        histogram_index == number_of_bins+1 ? histogram_index-=1 : histogram_index
        maxXs[histogram_index] = max(current_pos, maxXs[histogram_index])
        minXs[histogram_index] = min(current_pos, minXs[histogram_index])
    end
    return maxXs, minXs, bin_bounds
end

mean_diff(t) = median(t[2:end] .- t[1:end-1])
function get_fft(t, positions, Δt=mean_diff(t), windowing_function=blackman)
    Interpolations.deduplicate_knots!(t)
    interp = LinearInterpolation(t, positions, extrapolation_bc = Line())
    t = t[1]:Δt:t[end]
    positions = [interp(t_i) for t_i in t]
    n = length(t)
    windowing = windowing_function(n)
    @. positions = positions * windowing
    unique = ceil(Int, (n+1)/2)
    sampling_freq = 1/(t[end]/n)

    freq_amplitudes = abs.(fft(positions)[1:unique])/ sum(windowing) *2

    freq_list = (0:(unique-1)) * sampling_freq/n
    return freq_list/1e3, freq_amplitudes
end
export get_fft


function quinn_fernandes(signal, seed)
    eps = 1.0e-6
    n = length(signal)
    a = 2.0*cos(seed)
    b = 0.1
    signal = vcat([0,0], signal)
    xi = zeros(size(signal))

    it = 0
    while (abs(a-b) > eps)
        it = it+1
        if (it > 1000)
            b = NaN
            break
        end
        a = 2*b - a
        for i=3:(n-1)
            xi[i] = signal[i] + a*xi[i-1] - xi[i-2]
        end
        b = sum([(xi[i]+xi[i-2])*xi[i-1] for i=3:(n-1)])/sum([xi[i-1]^2 for i=3:(n-1)])
    end

    acos(0.5*b)
end
function get_freq_qf(t, positions; Δt=mean_diff(t))
    dcoffset = mean(positions)
    positions = @. positions - dcoffset

    t = copy(t)
    Interpolations.deduplicate_knots!(t)
    interp = LinearInterpolation(t, positions, extrapolation_bc = Line())
    t = t[1]:Δt:t[end]
    positions = [interp(t_i) for t_i in t]
    n = length(t)
    sampling_freq = 1/Δt

    freq = quinn_fernandes(positions, 0.5)

    return freq * sampling_freq / 2pi
end
export get_freq_qf

function get_fft_welch(t, positions, Δt=t[2]-t[1])
    interp = LinearInterpolation(t, positions, extrapolation_bc = Line())
    t = t[1]:Δt:t[end]
    positions = [interp(t_i) for t_i in t]
    n = length(t)
    pgram = welch_pgram(positions, div(length(positions), 8), 0, window=hanning)
    return freq(pgram)/Δt/1e3, power(pgram)
end

function get_freq_fit(t, positions)
    dcoffset = mean(positions)
    positions = @. positions - dcoffset
    x,y = get_fft(t, positions)
    omega=2*pi*x[findmax(y)[2]]*1e3
    p0 = zeros(3)

    p0[2] = omega
    p0[1] = findmax(y)[1]

    oscillation_model(t, p) = @. p[1] * cos(p[2]*t + p[3])

    fitres = curve_fit(oscillation_model, t, positions, p0)
    return fitres.param[2] / 2pi
end
export get_freq_fit

function get_tune(t, positions)
    dcoffset = mean(positions)
    positions = @. positions - dcoffset
    x,y = get_fft(t, positions)
    omega=2*pi*x[findmax(y)[2]]*1e3
    p0 = zeros(3)

    p0[2] = omega
    p0[1] = findmax(y)[1]

    oscillation_model(t, p) = @. p[1] * cos(p[2]*t + p[3])

    fitres = curve_fit(oscillation_model, t, positions, p0)
    return fitres
end

function get_freq_fft(t, positions; maxfreq=1000)
    dcoffset = mean(positions)
    positions = @. positions - dcoffset
    x,y = get_fft(t, positions)
    peaks = findpeaks(y[x .< maxfreq])
    return x[x .< maxfreq][peaks[1]]*1e3
end
export get_freq_fft

function get_tune_fft(t, positions; maxfreq=1000, cyclotron_freq=p.cyclotron_freq)
    return get_freq_fft(t, positions; maxfreq=maxfreq) / cyclotron_freq
end
export get_tune_fft


function get_vertical_tune(sol::ODESolution, p=p; maxfreq=1000)
    t = sol[time_index,:]
    ys = sol[3,:]
    return get_tune_fft(t, ys, maxfreq=maxfreq,cyclotron_freq=p.cyclotron_freq)
end


function get_horizontal_tune(sol::ODESolution, p=p; maxfreq=1000)
    t = sol[time_index,:]
    ys = sol[1,:]
    return get_tune_fft(t, ys, maxfreq=maxfreq,cyclotron_freq=p.cyclotron_freq)
end

function get_2_case_precession_data(sols)
    tdata = sols[1][time_index,:]
    ydata = (sols[1][8,:] .- sols[3][8,:]) ./ 2

    return tdata, ydata
end

function get_4_case_precession_data(sols)
    tdata = sols[1][time_index,:]
    ydata = (sols[1][8,:] .+ sols[2][8,:] .- sols[3][8,:] .- sols[4][8,:]) ./ 4

    return tdata, ydata
end


lin_model(t, p) = @. p[1] + p[2]*t
quad_model(t, p) = @. p[1] + p[2]*t + p[3]*t*t
cubic_model(t ,p) = @. p[1] + p[2]*t + p[3]*t*t + p[4] * t*t*t
power_model(t, p) = @. p[1] + p[2]*t + p[3]/t
sin_model(t, p) = @. p[1] + p[2]*t + p[3] * sin(p[4]*t + p[5])

function model_fit(x, y, model; verbose=true)
    if model==lin_model
        dims=2
    elseif model==quad_model
        dims=3
    else
        dims=4
    end

    p0 = zeros(dims)
    fitresult = curve_fit(model, x, y, p0, maxIter=5000, x_tol=1e-11, g_tol=1e-14)
    as = fitresult.param
    aerrs = stderror(fitresult)

    if verbose
        for i in 1:dims
            printfmtln("a$(i) ± σ = {:.3e} ± {:.3e}", as[i], aerrs[i])
        end   
    end
    return fitresult;
end

function linear_fit(x, y, verbose=false)
    #    x_tol::Real=1e-8`: search tolerance in x
    #   `g_tol::Real=1e-12`: search tolerance in gradient
    #   `maxIter::Integer=1000`: maximum number of iterations
    p0 = zeros(2)
    fitresult = curve_fit(lin_model, x, y, p0, maxIter=5000, x_tol=1e-11, g_tol=1e-14)
    a0, a1 = fitresult.param
    a0err, a1err = stderror(fitresult)
    if verbose
        printfmtln("a0 ± σ = {:.3e} ± {:.3e}", a0, a0err)
        printfmtln("a1 ± σ = {:.3e} ± {:.3e}", a1, a1err)
    end
    return fitresult;
end


function linear_fit(x, y::AbstractVector{Measurement{T}}, verbose=false) where T<:Real
    #    x_tol::Real=1e-8`: search tolerance in x
    #   `g_tol::Real=1e-12`: search tolerance in gradient
    #   `maxIter::Integer=1000`: maximum number of iterations
    y_data = Measurements.value.(y)
    y_err = Measurements.uncertainty.(y)
    p0 = zeros(2)
    fitresult = curve_fit(lin_model, x, y_data, 1 ./y_err.^2, p0, maxIter=5000, x_tol=1e-11, g_tol=1e-14)
    a0, a1 = fitresult.param
    a0err, a1err = stderror(fitresult)
    if verbose
        printfmtln("a0 ± σ = {:.3e} ± {:.3e}", a0, a0err)
        printfmtln("a1 ± σ = {:.3e} ± {:.3e}", a1, a1err)
    end
    return fitresult;
end
export linear_fit

function quadratic_fit(x, y, verbose=true)
    #    x_tol::Real=1e-8`: search tolerance in x
    #   `g_tol::Real=1e-12`: search tolerance in gradient
    #   `maxIter::Integer=1000`: maximum number of iterations
    a0,a1 = linear_fit(x, y, false).param
    p0 = zeros(3)
    p0[1] = a0
    p0[2] = a1
    fitresult = curve_fit(quad_model, x, y, p0, maxIter=5000, x_tol=1e-11, g_tol=1e-14)
    a0, a1, a2 = fitresult.param
    a0err, a1err, a2err = stderror(fitresult)
    if verbose
        printfmtln("a0 ± σ = {:.3e} ± {:.3e}", a0, a0err)
        printfmtln("a1 ± σ = {:.3e} ± {:.3e}", a1, a1err)
        printfmtln("a2 ± σ = {:.3e} ± {:.3e}", a2, a2err)
    end
    return fitresult;
end

include("rate_analysis.jl")

function get_average_per_meter(x, s, bin_size, p)
    positions = s .% p.ring.total_length

    bin_bounds = 0:bin_size:p.ring.total_length
    number_of_bins = length(bin_bounds) - 1
    avgs = zeros(number_of_bins)
    bin_belonging = zeros(length(positions))
    for i in 1:length(positions)
        current_s  =  positions[i]

        histogram_index = ceil(Int, current_s/bin_size)
        # histogram_index == 0 ? histogram_index=1 : histogram_index
        histogram_index == number_of_bins+1 ? histogram_index-=1 : histogram_index

        bin_belonging[i] = histogram_index
    end

    for i in 1:number_of_bins
        avgs[i] = mean(x[bin_belonging .== i] )
    end

    return avgs, bin_bounds
end

function get_dp_over_p(δ, particle)
    eta0 = particle.eta0
    p_magic = particle.momentum
    eta = eta0 * (1 + δ)
    p = p_magic * sqrt(eta*(2+eta)/eta0/(2+eta0))

    return (p - p_magic) / p_magic
end

function get_delta(dp_over_p, particle)
    p_magic = particle.momentum
    m_p = particle.m
    K0 = particle.K0
    p = dp_over_p * p_magic + p_magic
    K = sqrt((p*c)^2 + (m_p*c^2)^2) - (m_p*c^2)
    return (K - K0) / K0
end
export get_delta, get_dp_over_p

get_total_spin(sol) = sqrt.(sol[9,:].^2 + sol[7,:].^2 + sol[8,:].^2)
