import RichardsDAE as RD
import CSV
using DataFrames
using Plots

const COLORS = RD.okabe_ito_colors()

ptotal = plot(layout = (1, 2), size = (900, 400), margin = 5Plots.mm)
p1 = ptotal[1, 1]
p2 = ptotal[1, 2]

df = CSV.read("cases/output/haverkamp-final-head.csv", DataFrame)

# Add the mixed formulations as (sampled) scatter
labels = ("Mixed-BDF1 (Δt=120.0)", "Mixed-BDF2 (Δt=120.0)")
colors = (COLORS[:light_blue], COLORS[:green])
zmid = collect(0.5:1.0:39.5)
for (label, color) in zip(labels, colors)
    scatter!(
        p1,
        zmid,
        df[!, Symbol(label)],
        label = label,
        color = color,
        markersize = 3,
        markerstrokewidth = 0,
    )
end

labels = (
    "Reduced-BDF1 (Δt=1.0)",
    "ψ-based (Δt=120.0)",
    "Reduced-BDF1 (Δt=120.0)",
    "Reduced-BDF2 (Δt=120.0)",
)
colors = (COLORS[:black], COLORS[:dark_orange], COLORS[:light_blue], COLORS[:green])
plot!(p1, xlabel = "Elevation (cm)", ylabel = "Pressure head (cm)", dpi = 300)
for (label, color) in zip(labels, colors)
    plot!(p1, zmid, df[!, Symbol(label)], label = label, color = color, lw = 2)
end


df = CSV.read("cases/output/newmexico-final-head.csv", DataFrame)
columns = [
    "Reduced-BDF1 (dt=0.000001)",
    "Reduced-BDF1 (dt=0.1)",
    "Reduced-BDF2 (dt=0.1)",
    "Reduced-BDF3 (dt=0.1)",
]
labels = [
    "Reduced-BDF1 (Δt=1.0E-6)",
    "Reduced-BDF1 (Δt=0.1)",
    "Reduced-BDF2 (Δt=0.1)",
    "Reduced-BDF3 (Δt=0.1)",
]
plot!(p2, xlabel = "Elevation (m)", ylabel = "Pressure head (m)", dpi = 300)
zmid = collect(1.25e-3:2.5e-3:0.3-1.25e-3)
for (name, label, color) in zip(columns, labels, colors)
    plot!(p2, zmid, df[!, Symbol(name)], label = label, color = color, lw = 2)
end
savefig(ptotal, "cases/output/head.pdf")
