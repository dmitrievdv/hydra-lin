using Makie
using GLMakie
using DifferentialEquations
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

function calc_main(n_H, T_l)
W = 0.01
T_s = 4000
# n_H = 1e11

z = 1e11
v = 1e7
# T_l = 8000
ΔT_l = 1000
k = 1e11
n_s = 100
n_levels = 15

simple_ne_stat_result = calc_stat_data(n_H, T_l, ΔT_l, k, W, T_s, z, v, 100; n_levels = 15, simple_ne = true, plane_parallel = true)
stat_result = calc_stat_data(n_H, T_l, ΔT_l, k, W, T_s, z, v, 100; n_levels = 15, simple_ne = false, plane_parallel = true)

simple_ne_nonstat_result = calc_nonstat_data(n_H, T_l, ΔT_l, k, W, T_s, z, v, 100; n_levels = 15, simple_ne = true, plane_parallel = true)
nonstat_result = calc_nonstat_data(n_H, T_l, ΔT_l, k, W, T_s, z, v, 100; n_levels = 15, simple_ne = false, plane_parallel = true)

s_arr = range(0, 2k, n_s)

return s_arr, simple_ne_stat_result, simple_ne_nonstat_result, stat_result, nonstat_result
end

function plot_main(s_arr, simple_ne_stat_result, simple_ne_nonstat_result, stat_result, nonstat_result)
fig = Figure()
ax = Axis(fig[1,1])

n_H = sum(stat_result[1, :])
# min_f = minimum(log.( f_arr_stat ./ (1 .- f_arr_stat)))
# max_f = maximum(log.( f_arr_stat ./ (1 .- f_arr_stat)))
# f_range = max_f - min_f

# f_span = range(min_f - f_range/2, max_f + f_range/2, n_s)

# ax2 = Axis(fig[2,1])

# k_sl = Slider(fig[3,1], range = -0:k/100:k, startvalue = k/2)

# deriv_points = lift(k_sl.value) do k_val
#     deriv_points = [Point2f(0.0,0.0) for i_f = 1:length(f_span)]
#     # dl_f_arr = log10.(abs.(diff_f.(f_span, Ref([n_H, T_s, W, z, v, T_l, ΔT_l, k]), k_val)))

#     for i_f = 1:length(f_span)
#         f = f_span[i_f]
#         dl_f = log10(abs(diff_f(f, [n_H, T_s, W, z, v, T_l, ΔT_l, k], k_val)))
#         deriv_points[i_f] = Point2f(f, dl_f)
#     end
#     println("done")
#     deriv_points
# end

simple_ne_stat_f_arr = simple_ne_stat_result[:, 1] / n_H
stat_f_arr = stat_result[:, 1] / n_H

simple_ne_nonstat_f_arr = simple_ne_nonstat_result[:, 1] / n_H
nonstat_f_arr = nonstat_result[:, 1] / n_H

lines!(ax, s_arr, ( 1 .- stat_f_arr ), color = :black)
lines!(ax, s_arr, ( 1 .-  simple_ne_stat_f_arr ), color = :black, linestyle = :dash)

lines!(ax, s_arr, ( 1 .-  nonstat_f_arr ), color = :red)
lines!(ax, s_arr, ( 1 .-  simple_ne_nonstat_f_arr ), color = :red, linestyle = :dash)

# vlines!(ax, k_sl.value)
# lines!(ax2, deriv_points)
fig
end

s_arr, simple_ne_stat_result, simple_ne_nonstat_result, stat_result, nonstat_result = calc_main(1e9, 11000)
plot_main(s_arr, simple_ne_stat_result, simple_ne_nonstat_result, stat_result, nonstat_result)