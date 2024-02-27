using SEMTor 

p = StaticParameters()
ens = simulate_ensemble(p)

plot_state(ens[1][1], p)

st = compute_statistics(ens, p)

# example analysis
using CairoMakie
df = st.stats_mean["tissue"]

T = df.time
X = df.apical_height_mean
X_std = df.apical_height_std

lines(T, X, axis = (xlabel = "time [hr]", ylabel = "tissue height [5µm]"))
band!(T, X - X_std, X + X_std, color=(:gray, 0.3))
current_figure()
