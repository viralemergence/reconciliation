import CSV
using DataFrames
using EcologicalNetworks
using EcologicalNetworksPlots
using ProgressMeter
using TSne
using StatsPlots

# Load data
clover = DataFrame(CSV.File("../CLOVER_0.1_MammalViruses_AssociationsFlatFile.csv"))

viruses = unique([v.Virus for v in eachrow(clover)])
hosts = unique([v.Host for v in eachrow(clover)])
A = zeros(Bool, (length(viruses), length(hosts)))
M = BipartiteNetwork(A, viruses, hosts)

networks = Dict([u => copy(M) for u in unique(clover.Database)])

for r in eachrow(clover)
    virus = r.Virus
    host = r.Host
    networks[r.Database][virus,host] = true
    M[virus,host] = true
end

for (k,v) in networks
    simplify!(v)
end

# Color palette
cbfp = [colorant"#e69f00", colorant"#56b4e9", colorant"#009e73", colorant"#f0e442", colorant"#0072b2", colorant"#d55e00", colorant"#cc79a7"]

# t-SNE on the unipartite projection
UM = EcologicalNetworks.mirror(convert(UnipartiteNetwork, M))
x = tsne(convert.(Float64, Array(UM.edges)), 2, 0, 2000, 6)

I = initial(RandomInitialLayout, UM)
for (i,s) in enumerate(species(UM))
    I[s].x = x[i,1]
    I[s].y = x[i,2]
end

scatter(I, M, bipartite=true, nodesize=degree(M), msc=:black, size=(1000, 1000), aspectratio=1, dpi=600, msw=0.5)
savefig("tSNE.pdf")
savefig("tSNE.png")

ps = scatter(I, M, bipartite=true, nodesize=degree(M), msc=:black, size=(1000, 1000), dpi=600, msw=0.3)
scatter!(ps, I, networks["Shaw"], mc=cbfp[1], bipartite=true, nodesize=degree(M), msc=:black, msw=0.3)
savefig(ps, "shaw.pdf")

ph = scatter(I, M, bipartite=true, nodesize=degree(M), msc=:black, size=(1000, 1000), aspectratio=1, dpi=600, msw=0.5)
scatter!(ph, I, networks["HP3"], mc=cbfp[5], bipartite=true, nodesize=degree(M), msc=:black, msw=0.5)
savefig(ph, "hp3.pdf")

pe = scatter(I, M, bipartite=true, nodesize=degree(M), msc=:black, size=(1000, 1000), aspectratio=1, dpi=600, msw=0.5)
scatter!(pe, I, networks["EID2"], mc=cbfp[3], bipartite=true, nodesize=degree(M), msc=:black, msw=0.5)
savefig(pe, "eid2.pdf")

pg = scatter(I, M, bipartite=true, nodesize=degree(M), msc=:black, size=(1000, 1000), aspectratio=1, dpi=600, msw=0.5)
scatter!(pg, I, networks["GMPD2"], mc=cbfp[6], bipartite=true, nodesize=degree(M), msc=:black, msw=0.5)
savefig(pg, "gmpd2.pdf")

title!(ps, "Shaw")
title!(ph, "HP3")
title!(pe, "EID2")
title!(pg, "GMPD2")

plot(ps, ph, pe, pg)
savefig("clover.pdf")
savefig("clover.png")
savefig("clover.jpg")
savefig("clover.tif")
