using Makie
using GLMakie
using DifferentialEquations
# using BSON
# import REPL

include("ampsolvers.jl")



function calc_T_l(s, T_l, ΔT_l, k)
    return T_l + ΔT_l/2*sin(2π*s/k)
end

function calc_s_arr(k, n_s)
    return range(0, 2k, n_s)
end

function calc_stat_data(n_H, T_l, ΔT_l, k, W, T_s, z, v, n_s = 100; n_levels = 15, simple_ne = false, plane_parallel = true)
    s_arr = range(0, 2k, n_s)
    stat_data_result = zeros(n_s, n_levels+1)

    solve_kwargs = (:n_levels => n_levels, :simple_ne => simple_ne, :plane_parallel => plane_parallel)

    for i_s = 1:n_s
        s = s_arr[i_s]
        popul = solvestationary(n_H, calc_T_l(s, T_l, ΔT_l, k), W, T_s, z, v; solve_kwargs...)
        n_e = n_H - sum(popul)
        stat_data_result[i_s, 1:n_levels] .= popul
        stat_data_result[i_s, n_levels+1] = n_e
    end
    return stat_data_result
end

function calc_nonstat_data(n_H, T_l, ΔT_l, k, W, T_s, z, v, n_s = 100; n_levels = 15, simple_ne = false, plane_parallel = true)
    s_arr = range(0, 2k, n_s)
    nonstat_data_result = zeros(n_s, n_levels+1)

    solve_kwargs = (:n_levels => n_levels, :simple_ne => simple_ne, :plane_parallel => plane_parallel)

    popul = solvestationary(n_H, calc_T_l(0, T_l, ΔT_l, k), W, T_s, z, v; solve_kwargs...)
    
    n_e = n_H - sum(popul)
    n_1 = popul[1]
    f_0 = n_1 / n_H
    l_f_0 = log(f_0/(1 - f_0))

    function diff_f(l_f, p, s)
        # println(s)
        n_H, T_s, W, z, v, T_l, ΔT_l, k = p
        f = 1/(1 + exp(-l_f))
        n_1 = f*n_H
        return calc_logarithmic_derivative(n_1, n_H, calc_T_l(s, T_l, ΔT_l, k), W, T_s, z, v; solve_kwargs...)
    end

    prob = ODEProblem(diff_f, l_f_0*(1 + 1e-3), (0.0, 2k), p = [n_H, T_s, W, z, v, T_l, ΔT_l, k])
    min_s_step = 1e-3*k
    sol = solve(prob, dtmax = min_s_step, abstol = 1e-3)

    

    l_f_arr = sol.(s_arr)
    f_arr = 1 ./ (1 .+ exp.(-l_f_arr))
    nonstat_data_result[:, 1] .= f_arr * n_H

    for i_s = 1:n_s
        s = s_arr[i_s]
        n_1 = nonstat_data_result[i_s, 1]
        popul = solve_stationary_upper_levels(n_1, n_H, calc_T_l(s, T_l, ΔT_l, k), W, T_s, z, v; solve_kwargs...)
        n_e = if simple_ne 
            n_H - n_1
        else
            n_H - sum(popul)
        end
        nonstat_data_result[i_s, 2:n_levels] .= popul[2:n_levels]
        nonstat_data_result[i_s, n_levels+1] = n_e
    end
    return nonstat_data_result
end

function calc_main(n_H, T_l, W, T_s, z, v, k, ΔT_l; n_s = 100, n_levels = 15)

    simple_ne_stat_result = calc_stat_data(n_H, T_l, ΔT_l, k, W, T_s, z, v, 100; n_levels = 15, simple_ne = true, plane_parallel = true)
    stat_result = calc_stat_data(n_H, T_l, ΔT_l, k, W, T_s, z, v, 100; n_levels = 15, simple_ne = false, plane_parallel = true)

    simple_ne_nonstat_result = calc_nonstat_data(n_H, T_l, ΔT_l, k, W, T_s, z, v, 100; n_levels = 15, simple_ne = true, plane_parallel = true)
    nonstat_result = calc_nonstat_data(n_H, T_l, ΔT_l, k, W, T_s, z, v, 100; n_levels = 15, simple_ne = false, plane_parallel = true)

    s_arr = range(0, 2k, n_s)

    return s_arr, simple_ne_stat_result, simple_ne_nonstat_result, stat_result, nonstat_result
end

function plot_main(s_arr, simple_ne_stat_result, simple_ne_nonstat_result, stat_result, nonstat_result, T_l, ΔT_l, k)
    fig = Figure()
    ax_T_l = Axis(fig[1,1])
    ax = Axis(fig[2,1])
    ax_diff = Axis(fig[3,1])

    n_H = sum(stat_result[1, :])

    simple_ne_stat_f_arr = simple_ne_stat_result[:, 1] / n_H
    stat_f_arr = stat_result[:, 1] / n_H

    simple_ne_nonstat_f_arr = simple_ne_nonstat_result[:, 1] / n_H
    nonstat_f_arr = nonstat_result[:, 1] / n_H

    lines!(ax_T_l, s_arr, calc_T_l.(s_arr, T_l, ΔT_l, k), color = :black)

    lines!(ax, s_arr, ( 1 .- stat_f_arr ), color = :black)
    lines!(ax, s_arr, ( 1 .-  simple_ne_stat_f_arr ), color = :black, linestyle = :dash)

    lines!(ax, s_arr, ( 1 .-  nonstat_f_arr ), color = :red)
    lines!(ax, s_arr, ( 1 .-  simple_ne_nonstat_f_arr ), color = :red, linestyle = :dash)

    lines!(ax_diff, s_arr, simple_ne_stat_f_arr - stat_f_arr, color = :black)
    # lines!(ax, s_arr, ( 1 .-  simple_ne_stat_f_arr ), color = :black, linestyle = :dash)

    lines!(ax_diff, s_arr, simple_ne_nonstat_f_arr -  nonstat_f_arr, color = :red)
    # lines!(ax, s_arr, ( 1 .-  simple_ne_nonstat_f_arr ), color = :red, linestyle = :dash)
    # vlines!(ax, k_sl.value)
    # lines!(ax2, deriv_points)
    fig
end

function save_result(name, s_arr, result)
    open("1d-nonstat/$name.bin", "w") do io
        n_s = length(s_arr)
        write(io, Int64(n_s))
        for s in s_arr
            write(io, Float64(s))
        end
        n_levels = length(result[1,:]) - 1
        write(io, Int64(n_levels))
        for i_s = 1:n_s
            for i_lev = 1:n_levels+1
                write(io, Float64(result[i_s, i_lev]))
            end
        end
    end
end

function load_result(name)
    result = open("1d-nonstat/$name.bin", "r") do io
        n_s = read(io, Int64)
        s_arr = zeros(Float64, n_s)
        for i_s = 1:n_s
            s_arr[i_s] = read(io, Float64)
        end
        n_levels = read(io, Int64)
        result = zeros(Float64, (n_s, n_levels+1))
        for i_s = 1:n_s
            for i_lev = 1:n_levels+1
                result[i_s, i_lev] = read(io, Float64)
            end
        end
        result
        end
end

function load_s_arr(name)
    s_arr = open("1d-nonstat/$name.bin", "r") do io
        n_s = read(io, Int64)
        s_arr = zeros(Float64, n_s)
        for i_s = 1:n_s
            s_arr[i_s] = read(io, Float64)
        end
        s_arr
    end
end


begin 
W = 0.01
T_s = 4000
n_H = 1e12

z = 1e11
v = 2e7
T_l = 7000
ΔT_l = 1000
k = 1e11
n_s = 100
n_levels = 15


s_arr, simple_ne_stat_result, simple_ne_nonstat_result, stat_result, nonstat_result = calc_main(n_H, T_l, W, T_s, z, v, k, ΔT_l; n_s = n_s, n_levels = n_levels)


    model_name = "$(round(Int, log10(n_H)))_$(round(Int, T_l))"
    save_result("$(model_name)_stat_simple_ne", s_arr, simple_ne_stat_result)
    save_result("$(model_name)_nonstat_simple_ne", s_arr, simple_ne_nonstat_result)
    save_result("$(model_name)_stat", s_arr, stat_result)
    save_result("$(model_name)_nonstat", s_arr, nonstat_result)

    plot_main(s_arr, simple_ne_stat_result, simple_ne_nonstat_result, stat_result, nonstat_result, T_l, ΔT_l, k)
end

begin
    n_H = 1e10
    T_l = 10000
    ΔT_l = 1000
    k = 1e11
    model_name = "$(round(Int, log10(n_H)))_$(round(Int, T_l))"
s_arr = load_s_arr("$(model_name)_stat_simple_ne")
simple_ne_stat_result = load_result("$(model_name)_stat_simple_ne")
simple_ne_nonstat_result = load_result("$(model_name)_nonstat_simple_ne")
stat_result = load_result("$(model_name)_stat")
nonstat_result = load_result("$(model_name)_nonstat")
plot_main(s_arr, simple_ne_stat_result, simple_ne_nonstat_result, stat_result, nonstat_result, T_l, ΔT_l, k)
end


