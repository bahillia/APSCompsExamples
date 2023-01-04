using Plots
function getΨMatch(ΨBinned,levels=[-0.4*i for i=0:11])
    #match matplotlib default:
    #"fills intervals that are closed at the top;
    #that is, for regions z1 and z2 the filled region is z1 < Z <= z2"

    logΨ = log10.(ΨBinned)
    res = zeros(size(logΨ))
    mask = (logΨ .<= levels[2]) .& (logΨ .>= levels[1])
    res[mask] .= (levels[1]+levels[2])/2
    for i=2:length(levels)-1
        mask = (logΨ .< levels[i]) .& (logΨ .>= levels[i+1])
        res[mask] .= (levels[i]+levels[i+1])/2
    end
    mask = (logΨ .< levels[end])
    res[mask] = logΨ[mask]
    return res
end

function makePlot(yBinned,tBinned,ΨBinned,rs,c,stype=:heatmap;levels=[-0.4*i for i=0:11],size=(720,540),clims=(-5,0),A0=1,tlims=(0,40),vlims=(-12,12),
    plotLPExact=false,cList=[:forestgreen,:crimson]) #every previously passed param now a list of length 2 to plot 2 cases side by side

    ΨDiscrete = [getΨMatch(ΨBinned[i],levels) for i=1:2]
    clrticks = ["$(round(level,sigdigits=2))" for level in reverse(levels)]
    n = length(levels)
    yt = range(0,1,n)[1:n] #.+ 0.5/n add the half if levels are centered

    l = @layout [
        [a{0.4w} b{0.4w} c]
        [d{0.4h,0.8w} e]
        ]

    colors = [palette([:white,cList[1],:black],length(levels)-1),palette([:white,cList[2],:black],length(levels)-1)]

    p1 = plot(yBinned[1].*(c/1e6),tBinned[1][(tBinned[1].<=tlims[2]) .& (tBinned[1].>=tlims[1])],ΨDiscrete[1][:,(tBinned[1].<=tlims[2]) .& (tBinned[1].>=tlims[1])]',
        color=colors[1],cbar=false,tickfont="Computer Modern",guidefont="Computer Modern",
        xlims=vlims,ylims=tlims,seriestype=stype,fill=true,levels=n,bottom_margin=0*Plots.Measures.mm,
        xlabel="",ylabel="t [days]",minorticks=true,tickdirection=:out,minorgrid=true,clims=clims,
        framestyle=:box,right_margin=0*Plots.Measures.mm)
    #p2 = plot([NaN], lims=(0,1), framestyle=:none, legendDorodnitsyn=false) -- if you want a colorbar title

    xx = range(0,1,100)
    zz = zero(xx)' .+ xx
    p1 = plot!(title="",inset=(1,bbox(1/40,1/10,0.1,0.5)),titlefont="Computer Modern")
    p1 = plot!(p1[2],xx, xx, zz, ticks=false, ratio=10, legend=false, fc=colors[1], lims=(0,1),title="logΨ",
             framestyle=:box, right_margin=20*Plots.Measures.mm,seriestype=:heatmap,cbar=false,titlefontsize=10)

    for (yi,ti) in zip(yt,clrticks)
        p1=plot!(p1[2],annotations=(1.5,yi,text(ti, 7, "Computer Modern")))
    end

    p2 = plot(yBinned[2].*(c/1e6),tBinned[2][(tBinned[2].<=tlims[2]) .& (tBinned[2].>=tlims[1])],ΨDiscrete[2][:,(tBinned[2].<=tlims[2]) .& (tBinned[2].>=tlims[1])]',
        color=colors[2],cbar=false,tickfont="Computer Modern",guidefont="Computer Modern",
        xlims=vlims,ylims=tlims,seriestype=stype,fill=true,levels=n,bottom_margin=0*Plots.Measures.mm,
        xlabel="",ylabel="",minorticks=true,tickdirection=:out,minorgrid=true,clims=clims,yticks=false,
        framestyle=:box,right_margin=0*Plots.Measures.mm) #xticks=([-10,-5,0,5,10],["-10","-5","0","5","10"])
    #p2 = plot([NaN], lims=(0,1), framestyle=:none, legendDorodnitsyn=false) -- if you want a colorbar title

    xx = range(0,1,100)
    zz = zero(xx)' .+ xx
    p2 = plot!(title="",inset=(1,bbox(1/40,1/10,0.1,0.5)),titlefont="Computer Modern")
    p2 = plot!(p2[2],xx, xx, zz, ticks=false, ratio=10, legend=false, fc=colors[2], lims=(0,1),title="logΨ",
             framestyle=:box, right_margin=20*Plots.Measures.mm,seriestype=:heatmap,cbar=false,titlefontsize=10)

    for (yi,ti) in zip(yt,clrticks)
        p2=plot!(p2[2],annotations=(1.5,yi,text(ti, 7, "Computer Modern")))
    end
    #[annotate!(1.5, yi, text(ti, 7, "Computer Modern")) for (yi,ti) in zip(yt,clrticks)]
    #p1 = plot!(p1[2],annotations=[])
    #annotate!(2.2,0.5,text("logΨ",10,"Computer Modern",rotation=90))

    #now make Ψ(τ)
    #dt1 = [tBinned[2]-tBinned[1] for i=1:nBint]; dt2 = [tBinned[end]-tBinned[end-1] for i=1:64]
    #dt = vcat(dt1,dt2)./(3600*24)

    Ψτ1 = [(sum(ΨBinned[1][:,i])) for i=1:length(tBinned[1])]./(rs[1]/c) #for normalization make unitless again
    τ_mean1 = sum(tBinned[1].*Ψτ1)/sum(Ψτ1)
    p3=plot(Ψτ1,tBinned[1],label="",lw=4,color=cList[1],xlabel="Ψ(t)",framestyle=:box,minorticks=true,minorgrid=true,ymirror=true,
    xflip=true,guidefont="Computer Modern",tickfont="Computer Modern",xlims=(0.,0.01),ylims=tlims,
    xrotation=90,bottom_margin=0*Plots.Measures.mm,left_margin=0*Plots.Measures.mm,tickdirection=:out,titlefont="Computer Modern",titlefontsize=10)

    Ψτ2 = [(sum(ΨBinned[2][:,i])) for i=1:length(tBinned[2])]./(rs[2]/c) #for normalization make unitless again
    τ_mean2 = sum(tBinned[2].*Ψτ2)/sum(Ψτ2)
    p3=plot!(Ψτ2,tBinned[2],label="",lw=4,color=cList[2],title="") #mean delay = ($(round(τ_mean1,sigdigits=3)),$(round(τ_mean2,sigdigits=3))) days

    #now make LP
    LP1 = [sum(ΨBinned[1][i,:]) for i=1:length(yBinned[1])]
    p4=plot(yBinned[1].*(c/1e6),LP1./maximum(LP1),label="",color=cList[1],lw=4,minorticks=true,minorgrid=true,framestyle=:box,guidefont="Computer Modern",
        tickfont="Computer Modern",xlims=vlims,xlabel="Δv [Mm/s]",widen=false,ylims=(0,1.1),ylabel="Normalized flux Ψ(ν)",
        top_margin=0*Plots.Measures.mm,right_margin=0*Plots.Measures.mm,tickdirection=:out,yticks=[0.2*i for i=0:5])

    LP2 = [sum(ΨBinned[2][i,:]) for i=1:length(yBinned[2])]
    p4=plot!(yBinned[2].*(c/1e6),LP2./maximum(LP2),label="",color=cList[2],lw=4)
    if plotLPExact!=false
        p4=plot!(plotLPExact[2].*(c/1e6),plotLPExact[3]./maximum(plotLPExact[3]),label="exact",color=:dodgerblue)
    end
    #p4=plot(left_margin=0*Plots.Measures.mm,right_margin=0*Plots.Measures.mm,top_margin=0*Plots.Measures.mm,bottom_margin=0*Plots.Measures.mm)
    P=plot(p1, p2, p3, p4, layout=l, margins=0*Plots.Measures.mm,size=size,link=:both)
    return P
end
