#!/usr/bin/env julia
"""
    This file contains the code used to generate figures 2, and 6 in the paper.
    For the other figures see `figures.py`
"""

include("DiskWind/src/functions.jl")
using JLD

##################### SHOWING MODEL LINE/PHASE PROFILES (FIGURE 2) #####################

function modelViz(params::Array{Float64,},data::Array{Array{Float64,N} where N, 1},
    bins::Int=200,nr::Int=1024,nϕ::Int=2048,coordsType::Symbol=:polar,scale_type::Symbol=:log,νMin = 0.975, νMax = 1.025;genComparison=false,centered=true,noAbs=false,norm=false,old=false)

    i,r̄,Mfac,rFac,f1,f2,f3,f4,pa,scale,cenShift = params; γ = 1.; A0 = 1.; τ = 10. #some parameters fixed for now
    λCen = 2.172 + cenShift #microns, to compare with data
    #ν = (data[1].-2.172)./λCen #code units, cenShift should be small this is just for calculating min and max
    BLRAng = Mfac*3e8*2e33*6.67e-8/9e20/548/3.09e24 #solar masses * G / c^2 / Mpc -> end units = rad
    α,β,r,ν,ϕ,sini,cosi,dA,rMin,rMax = setup(i,nr,nϕ,r̄,rFac,γ,coordsType,scale_type)
    function getReturn(α,β,r,ν,ϕ,sini,cosi,dA,rMin,rMax,bins,νMin,νMax;f1,f2,f3,f4,total=false,old=false)
        I,γ,A0,τ = old == false ? getIntensity(r,ϕ,sini,cosi,rMin,rMax,γ,A0,τ,f1=f1,f2=f2,f3=f3,f4=f4,noAbs=noAbs) : getIntensityOld(r,ϕ,1.,sini,cosi,rMin,rMax,1.,1.,10.,f1=f1,f2=f2,f3=f3,noAbs=noAbs)
        νEdges,νCenters,flux = histSum(ν,I.*dA,bins=bins,νMin=νMin,νMax=νMax,centered=centered)

        #phases at resolution of binning
        UData = data[2]; VData = data[3]; psf=4e-3/2.35
        X = α.*BLRAng; Y = β.*BLRAng
        dϕList = []
        for i=1:length(UData)
            for ii in [I]
                dϕAvgRaw = phase(ν,ii,dA,X,Y,r,UData[i],VData[i],pa,νMin,νMax,bins) 
                dϕAvg = G1D(dϕAvgRaw,psf/3e5/(νCenters[2]-νCenters[1]))
                push!(dϕList,dϕAvg)
            end
        end
        #line profile
        if total == false
            fline = norm == true ? flux./maximum(flux)*maximum(data[4])*scale : flux
            lineAvg = G1D(fline,psf/3e5/(νCenters[2]-νCenters[1]))
            λ = λCen ./ νCenters
            return [reverse(λ),reverse(lineAvg),[dϕ.*reverse(lineAvg) for dϕ in dϕList]],maximum(flux)
        else
            fline = norm == true ? flux./maximum(flux)*maximum(data[4])*scale : flux
            lineAvg = G1D(fline,psf/3e5/(νCenters[2]-νCenters[1]))
            λ = λCen ./ νCenters
            return [reverse(λ),reverse(lineAvg),[dϕ.*reverse(lineAvg) for dϕ in dϕList]],genComparison
        end
    end

    ret1,overallMax = getReturn(α,β,r,ν,ϕ,sini,cosi,dA,rMin,rMax,bins,νMin,νMax,f1=f1,f2=f2,f3=f3,f4=f4)
    returnList = [ret1]
    if genComparison == true
        push!(returnList,getReturn(α,β,r,ν,ϕ,sini,cosi,dA,rMin,rMax,bins,νMin,νMax,f1=f1,f2=0.,f3=0.,f4=0.,total=overallMax)[1])
        push!(returnList,getReturn(α,β,r,ν,ϕ,sini,cosi,dA,rMin,rMax,bins,νMin,νMax,f1=0.,f2=f2,f3=0.,f4=0.,total=overallMax)[1])
        push!(returnList,getReturn(α,β,r,ν,ϕ,sini,cosi,dA,rMin,rMax,bins,νMin,νMax,f1=0.,f2=0.,f3=f3,f4=0.,total=overallMax)[1])
        push!(returnList,getReturn(α,β,r,ν,ϕ,sini,cosi,dA,rMin,rMax,bins,νMin,νMax,f1=0.,f2=0.,f3=0.,f4=f4,total=overallMax)[1])
    end
    return returnList
end

#paper plots
#1 f1 only
#making them in one plot with twin x
pLine = plot(size=(740,650),minorgrid=false,fontfamily="Computer Modern",tickdirection=:out,title="Line (top) and phase (bottom) profiles",xlabel="Δv [Mm/s]",ylabel="Normalized flux        ",yguidefontvalign = :top,
    left_margin=5*Plots.Measures.mm,ylims=(-1.1,1.1),box=:on,right_margin=20*Plots.Measures.mm,bottom_margin=5*Plots.Measures.mm,tickfontsize=10,labelfontsize=12,widen=:false,legendfontsize=12,
    titlefontsize=16,gridalpha=0.4,minorgridalpha=0.2,xlims=(-0.04,0.04),legend_position=(0.1,0.45),legend_background_color=:transparent,foreground_color_legend=:transparent,grid=false,
    xticks=([-0.0362,-0.0181,0,0.0181,0.0362],["-5.0","-2.5","0","2.5","5.0"]))

data = readPickle("3c273_juljanmarmay_append_gilles_specirf_wide_v6.p")
pPhase = twinx()
θ = [45.,6e3,1.,45.,1.,1.,1.,1.,320.,1.,0.]
retList = modelViz(θ,data,200,1024,2048,genComparison=true,centered=true,noAbs=false)
using LaTeXStrings
labels=["All terms equal",L"$f_1$ only",L"$f_2$ only",L"$f_3$ only", L"$f_4$ only"]
indx=[0,1,2,6,7,8,12,13,14,18,19,20].+1; oindx=[3,4,5,9,10,11,15,16,17,21,22,23].+1
maxf = maximum(retList[1][2])
#gridlines and ticks
p=plot!(pLine,[-0.04,0.04],[0.,0.],color=:grey,alpha=0.5,label="",minorticks=false)
p=plot!(pPhase,[-0.04,0.04],[0.,0.],color=:grey,alpha=0.5,label="")
p=plot!(pLine,[0.,0.],[-1.1,1.1],color=:grey,alpha=0.5,label="")
p=plot!(pLine,yticks=([0.0,0.25,0.5,0.75,1.0],["0.0","0.25","0.5","0.75","1.0"]))
p=plot!(pPhase,yticks=([-0.2,-0.1,0.0,0.1,0.2,0.3,0.4],["-0.2","-0.1","0.0","0.1","0.2","0.3","0.4"]))
for i=1:length(retList)
    λModel,lineModel,phaseListModel = retList[i]
    λModel .-= 2.172
    if i>1
        p=plot!(pLine,λModel,(lineModel)./maxf,label=labels[i],lw=3,α=0.8,c=i)
        phaseOn = mean(phaseListModel[indx],dims=1); phaseOff = mean(phaseListModel[oindx],dims=1)
        p=plot!(pPhase,λModel,phaseOn.*200,label="",lw=3,α=0.8,c=i)
        #p=plot!(pPhase,λModel,phaseOff,label="",lw=2,α=0.8)
    else
        p=plot!(pLine,λModel,(lineModel)./maxf,label=labels[i],lw=3,α=1.,c=i)#,marker=:circle,markerstrokewidth=0.,ms=4)
        phaseOn = mean(phaseListModel[indx],dims=1); phaseOff = mean(phaseListModel[oindx],dims=1)
        p=plot!(pPhase,λModel,phaseOn.*200,label="",lw=3,α=1.,c=i)#,marker=:circle,markerstrokewidth=0.,ms=4)
        #p=plot!(pPhase,λModel,phaseOff,label="",lw=3,α=1.,marker=:circle,markercolor=:black,ms=2)
    end
end
P1 = plot!(pPhase,xticks=:none,ylabel="        ΔΦ [deg]",ylims=(-0.3,1.2),box=:on,minorticks=false,tickdirection=:out,tickfontsize=10,labelfontsize=12,widen=:false,xlims=(-0.04,0.04),yguidefontvalign = :bottom)
savefig(P1,"allTerms.pdf")

#making them in one plot with twin x
pLine = plot(size=(740,650),minorgrid=false,fontfamily="Computer Modern",tickdirection=:out,title="Line (top) and phase (bottom) profiles",xlabel="Δv [Mm/s]",ylabel="Normalized flux        ",
    left_margin=5*Plots.Measures.mm,ylims=(-1.1,1.1),box=:on,right_margin=20*Plots.Measures.mm,bottom_margin=5*Plots.Measures.mm,tickfontsize=10,labelfontsize=12,widen=:false,legendfontsize=12,
    titlefontsize=16,gridalpha=0.4,minorgridalpha=0.2,xlims=(-0.04,0.04),legend_position=(0.1,0.45),legend_background_color=:transparent,foreground_color_legend=:transparent,grid=false,
    xticks=([-0.0362,-0.0181,0,0.0181,0.0362],["-5.0","-2.5","0","2.5","5.0"]),yguidefontvalign = :top)

pPhase = twinx()
#i,r̄,Mfac,rFac,f1,f2,f3,f4,pa,scale,cenShift = params
θ1 = [45.,6e3,1.,45.,1.,1.,1.,1.,320.,1.,0.]
θ2 = [45.,6e3,1.,45.,.25,1.,.25,.25,320.,1.,0.]
θ3 = [45.,6e3,1.,45.,.25,.25,1.,.25,320.,1.,0.]
θ4 = [45.,6e3,1.,45.,.25,1.,1.,.25,320.,1.,0.]
θList = [θ1,θ2,θ3,θ4]
retList = [modelViz(θ,data,200,1024,2048,genComparison=false,centered=true,noAbs=false) for θ in θList]
labels=["All terms equal",L"$f_2 = 4f_{1,2,3}$",L"$f_3 = 4f_{1,2,4}$", L"$f_{2,3} = 4f_{1,4}$"]
indx=[0,1,2,6,7,8,12,13,14,18,19,20].+1; oindx=[3,4,5,9,10,11,15,16,17,21,22,23].+1
maxf = maximum(retList[1][1][2])
#gridlines
p=plot!(pLine,[-0.04,0.04],[0.,0.],color=:grey,alpha=0.5,label="",minorticks=false,yticks=([0.0,0.25,0.5,0.75,1.0],["0.0","0.25","0.5","0.75","1.0"]))
p=plot!(pPhase,[-0.04,0.04],[0.,0.],color=:grey,alpha=0.5,label="",yticks=([-0.2,-0.1,0.0,0.1,0.2,0.3,0.4],["-0.2","-0.1","0.0","0.1","0.2","0.3","0.4"]))
p=plot!(pLine,[0.,0.],[-1.1,1.1],color=:grey,alpha=0.5,label="")
cList = [1,:crimson,:grey,:cyan]
for i=1:length(retList)
    λModel,lineModel,phaseListModel = retList[i][1]
    λModel .-= 2.172
    if i>1
        c = cList[i]
        p=plot!(pLine,λModel,(lineModel)./maxf,label=labels[i],lw=3,α=0.8,c=c)
        phaseOn = mean(phaseListModel[indx],dims=1); phaseOff = mean(phaseListModel[oindx],dims=1)
        p=plot!(pPhase,λModel,phaseOn.*200,label="",lw=3,α=0.8,c=c)
    else
        p=plot!(pLine,λModel,(lineModel)./maxf,label=labels[i],lw=3,α=1.,c=i)#,marker=:circle,markerstrokewidth=0.,ms=4)
        phaseOn = mean(phaseListModel[indx],dims=1); phaseOff = mean(phaseListModel[oindx],dims=1)
        p=plot!(pPhase,λModel,phaseOn.*200,label="",lw=3,α=1.,c=i)#,marker=:circle,markerstrokewidth=0.,ms=4)
    end
end
P2 = plot!(pPhase,xticks=:none,ylabel="        ΔΦ [deg]",ylims=(-0.3,1.2),box=:on,minorticks=false,tickdirection=:out,tickfontsize=10,labelfontsize=12,widen=:false,xlims=(-0.04,0.04),yguidefontvalign = :bottom)
savefig(P2,"f2f3.pdf")

P=plot(P1,P2,size=(1244,550),fontfamily="Computer Modern")
savefig(P,"combinedTerms.pdf")

####################### SHOWING MODEL ECHO IMAGES (FIGURE 6) #######################

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
    plotLPExact=false,cList=[:forestgreen,:crimson],Ψmax=0.01) #every previously passed param now a list of length 2 to plot 2 cases side by side

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
        xlabel="Δv [Mm/s]",ylabel="t [days]",minorticks=true,tickdirection=:out,minorgrid=true,clims=clims,
        framestyle=:box,right_margin=0*Plots.Measures.mm)
    #p2 = plot([NaN], lims=(0,1), framestyle=:none, legendDorodnitsyn=false) -- if you want a colorbar title

    xx = range(0,1,100)
    zz = zero(xx)' .+ xx
    p1 = plot!(title="",inset=(1,bbox(1/20,1/10,0.1,0.5)),titlefont="Computer Modern")
    p1 = plot!(p1[2],xx, xx, zz, ticks=false, ratio=10, legend=false, fc=colors[1], lims=(0,1),title="logΨ",
             framestyle=:box, right_margin=20*Plots.Measures.mm,seriestype=:heatmap,cbar=false,titlefontsize=10)

    for (yi,ti) in zip(yt,clrticks)
        p1=plot!(p1[2],annotations=(1.5,yi,text(ti, 7, "Computer Modern")))
    end

    p2 = plot(yBinned[2].*(c/1e6),tBinned[2][(tBinned[2].<=tlims[2]) .& (tBinned[2].>=tlims[1])],ΨDiscrete[2][:,(tBinned[2].<=tlims[2]) .& (tBinned[2].>=tlims[1])]',
        color=colors[2],cbar=false,tickfont="Computer Modern",guidefont="Computer Modern",
        xlims=vlims,ylims=tlims,seriestype=stype,fill=true,levels=n,bottom_margin=0*Plots.Measures.mm,
        xlabel="Δv [Mm/s]",ylabel="",minorticks=true,tickdirection=:out,minorgrid=true,clims=clims,yticks=false,
        framestyle=:box,right_margin=0*Plots.Measures.mm) #xticks=([-10,-5,0,5,10],["-10","-5","0","5","10"])
    #p2 = plot([NaN], lims=(0,1), framestyle=:none, legendDorodnitsyn=false) -- if you want a colorbar title

    xx = range(0,1,100)
    zz = zero(xx)' .+ xx
    p2 = plot!(title="",inset=(1,bbox(1/20,1/10,0.1,0.5)),titlefont="Computer Modern")
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
    xflip=true,guidefont="Computer Modern",tickfont="Computer Modern",xlims=(0.,Ψmax),ylims=tlims,
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

res_low = load("RMVars_low.jld")
res_high = load("RMVars_high.jld")
#these JLD files are generated from `RMDelaySummit.jl`

yBinned_low,tBinned_low,ΨBinned_low = res_low["result"]
yBinned_high,tBinned_high,ΨBinned_high = res_high["result"]
G = 6.67e-11; c = 3e8;
M_high = 8.34e7*2e30
M_low = 2.82e8*2e30
rs_low = 2*G*M_low/c^2; rs_high = 2*G*M_high/c^2
P = makePlot([yBinned_low,yBinned_high],[tBinned_low,tBinned_high],[ΨBinned_low,ΨBinned_high],[rs_low,rs_high],c,vlims=(-6,6),tlims=(0,100),size=(650,600))
savefig(P,"echoImage.pdf")