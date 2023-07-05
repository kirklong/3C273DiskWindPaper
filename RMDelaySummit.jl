#!/usr/bin/env julia
using JLD
include("DiskWind/src/functions.jl")
function setup2(i::Float64,n1::Int64,n2::Int64,r̄::Float64,rFac::Float64,rMin::Float64,rMax::Float64,coordsType::Symbol=:polar,scale::Symbol=:log)
    #rMin,rMax = get_rMinMax(r̄,rFac,γ)
    i = i/180*π; cosi = cos(i); sini = sin(i) #inclination angle in rad
    α = nothing; β = nothing; r = nothing; ν = nothing; ϕ = nothing; dA = nothing
    if coordsType == :cartesian
        nx = n1; ny = n2; rlim = rMax
        a = nothing; b = nothing

        if scale == :linear
            a = range(-rlim,stop=rlim,length=nx); b = range(-rlim,stop=rlim,length=ny)
        elseif scale == :log
            a = vcat(-reverse(exp.(range(log(rMin*cosi),stop=log(rMax),length=Int(nx/2)))),exp.(range(log(rMin*cosi),stop=log(rMax),length=Int(nx/2))))
            b = vcat(-reverse(exp.(range(log(rMin*cosi),stop=log(rMax),length=Int(ny/2)))),exp.(range(log(rMin*cosi),stop=log(rMax),length=Int(ny/2))))
        else
            println("invalid scale symbol -- should be :linear or :log")
            exit()
        end

        α,β = meshgrid(a,b)
        α = reshape(α,nx,ny); β = reshape(β,nx,ny)

        dA = zeros(size(α))
        for i=1:size(dA)[1]
            for j=1:size(dA)[2]
                Δα = i<size(dA)[1] ? abs(α[i+1,j]-α[i,j]) : abs(α[end,j]-α[end-1,j]) #kinda bs but want n things and linear spacing so fine
                Δβ = j<size(dA)[2] ? abs(β[i,j+1]-β[i,j]) : abs(β[i,end]-β[i,end-1])
                dA[i,j] = Δα*Δβ
            end
        end

        r = reshape(sqrt.(β.^2 ./cosi^2 .+ α.^2),nx,ny); ϕ = reshape(atan.(β./cosi,α),nx,ny)
        ν = 1 .+ sqrt.(1 ./(2 .*r)).*sini.*cos.(ϕ)

    elseif coordsType == :polar
        nr = n1; nϕ = n2
        offset = 0.
        ϕ = range(0+offset,stop=2π+offset,length=nϕ+1)[1:end-1] 
        r = nothing; rGhost = nothing; Δr = nothing; Δlogr = nothing

        if scale == :linear
            r = range(rMin*cosi,stop=rMax,length=nr)
            Δr = r[2]-r[1]
            rGhost = [rMin*cosi-Δr,rMax*cosi+Δr]
        elseif scale == :log
            logr = range(log(rMin*cosi),stop=log(rMax),length=nr)
            Δlogr = logr[2]-logr[1]
            rGhost = exp.([log(rMin*cosi)-Δlogr,log(rMax)+Δlogr])
            r = exp.(logr)
        else
            println("invalid scale symbol -- should be :linear or :log")
            exit()
        end

        rMesh, ϕMesh = meshgrid(r,ϕ)
        rMesh = reshape(rMesh,nr,nϕ); ϕMesh = reshape(ϕMesh,nr,nϕ)
        α,β = rMesh.*cos.(ϕMesh), rMesh.*sin.(ϕMesh)
        Δϕ = ϕ[2]-ϕ[1]
        dA = zeros(size(rMesh))
        for i=1:size(dA)[1]
            for j=1:size(dA)[2]
                if scale == :log
                    Δr = rMesh[i,j]*Δlogr
                end
                dA[i,j] = rMesh[i,j]*Δϕ*Δr
            end
        end

        r = reshape(sqrt.(β.^2/cosi^2 .+ α.^2),nr,nϕ); ϕ = reshape(atan.(β./cosi,α),nr,nϕ)
        ν = 1 .+ sqrt.(1 ./(2 .* r)).*sini.*cos.(ϕ)

    else
        println("invalid coords system -- should be :cartesian or :polar")
        exit()
    end
    return reshape(α,n1,n2),reshape(β,n1,n2),r,ν,ϕ,sini,cosi,dA,rMin,rMax
end

function mainCode(nR,nϕ,nBinv,nBint;tMaxR=20.,directCompare = false,binCenter = true,rMin=0.5e3,rMax=1e4,
    i=78.,rBar=6.4e3,M=8.4e7*2e30,rFac=47.,f1=0.69,f2=0.8,f3=0.65,f4=0.44,coords=:polar,scale=:log,gamma=1.)
    
    println("Running with $(Threads.nthreads()) threads")
    start = time()
    G = 6.67e-11; c = 2.99e8; days =3.6e3*24.
    rs = 2*G*M/c^2
    α,β,rArr,nuArr,ϕArr,sini,cosi,dA,rMin,rMax = directCompare == true ? setup(i,nR,nϕ,rBar,rFac,gamma,coords,scale) : setup2(i,nR,nϕ,rBar,rFac,rMin,rMax,coords,scale)

    ϕ′ = ϕArr .+ π/2 #waters ϕ
    tArr = rArr .* (rs/c/days) .* (1 .- cos.(ϕ′).*sini)
    yArr = -sqrt.(1 ./ (2 .*rArr)).*sin.(ϕ′).*sini

    ICode,γ,A0,τ = getIntensity(rArr,ϕArr,sini,cosi,rMin,rMax,gamma,1.,10.,f1=f1,f2=f2,f3=f3,f4=f4,test=false,noAbs=false) #this does the Waters conversion within

    ΨCode = ICode.*dA

    yMax = maximum(yArr); yMin = minimum(yArr)
    tMax = maximum(tArr); tMin = minimum(tArr)
    yBinned = range(yMin,stop=yMax,length=nBinv);
    tBinned=range(tMin,stop=tMaxR,length=nBint)
    Δy = yBinned[2]-yBinned[1]; Δt = tBinned[2]-tBinned[1]
    Δt = tBinned[2]-tBinned[1]
    ΨBinned = zeros(length(yBinned),length(tBinned))

    if binCenter == true
        Threads.@threads for i=1:length(yBinned)
            #print("$(round(i/length(yBinned)*100,sigdigits=3)) % complete \r")
            yMinTmp = yBinned[i]-Δy/2; yMaxTmp = yBinned[i] + Δy/2
            for j=1:length(tBinned)
                tMinTmp = tBinned[j]-Δt/2; tMaxTmp = tBinned[j] + Δt/2
                mask = (yArr .>= yMinTmp) .& (yArr .< yMaxTmp) .& (tArr .>= tMinTmp) .& (tArr .< tMaxTmp) .& (rArr .> rMin)
                s = sum(ΨCode[mask])
                if s > 0
                    ΨBinned[i,j] = s
                else
                    ΨBinned[i,j] = 1e-30
                end
            end
        end
    else
        ΨBinned = zeros(length(yBinned)-1,length(tBinned)-1)
        Threads.@threads for i=1:length(yBinned)-1
            #print("$(round(i/length(yBinned)*100,sigdigits=3)) % complete \r")
            yMinTmp = yBinned[i]; yMaxTmp = yBinned[i+1]
            for j=1:length(tBinned)-1
                tMinTmp = tBinned[j]; tMaxTmp = tBinned[j+1]
                mask = (yArr .>= yMinTmp) .& (yArr .< yMaxTmp) .& (tArr .>= tMinTmp) .& (tArr .< tMaxTmp) .& (rArr .> rMin)
                s = sum(ΨCode[mask])
                if s > 0
                    ΨBinned[i,j] = s
                else
                    ΨBinned[i,j] = 1e-30
                end
            end
        end
        yBinned = yBinned[1:end-1]; tBinned = tBinned[1:end-1] #only return left edges
    end
    LP = histSum(yArr,ICode.*dA,bins=nBinv,νMin=minimum(yArr),νMax=maximum(yArr))
    A0=maximum(LP[3])
    ΨBinned[ΨBinned .== 0] .= 1e-30 #get rid of zero values so there are no log errors
    ΨBinned = ΨBinned ./ A0 #A0 = max(Ψ)
    finish = time()
    println("took $(round((finish-start)/60,sigdigits=3)) min for (nR,nϕ,nBinv,nBint) = ($nR,$nϕ,$nBinv,$nBint)")
    result = [yBinned,tBinned,ΨBinned]
    save("RMVars.jld","result",result)
    return 0
end

println("starting now")
#for low inclination
mainCode(8192*8,2048*4,512,2048,tMaxR=3e3,directCompare=true,binCenter=false,gamma=1.,i=29.,rBar=1636.,M=0.94*3e8*2e30,rFac=42.,f1=0.54,f2=0.91,f3=0.04,f4=0.)
#for high inclination
mainCode(8192*8,2048*4,512,2048,tMaxR=6e3,directCompare=true,binCenter=false,gamma=1.,i=78.,rBar=6400.,M=8.4e7*2e30,rFac=47.,f1=0.69,f2=0.80,f3=0.65,f4=0.44)
#parameters come from fit results (avgs)