module Giggle

#PACKAGES USED
using Combinatorics
using DataFrames
using Distances
using Gurobi
using JuMP
using MathOptInterface
using Random
using UUIDs
using XLSX
using Query

#STRUCTS USED
struct vtx
    #IDENTIFIERS
    name::String
    type::String

    #PARAMS
    x::Float64 #x coor
    y::Float64 #y coor
    MAX::Float64 #max inventory level
    MIN::Float64 #min inventory level
    START::Float64 #starting inventory level

    #COSTS
    h::Float64
end

struct veh
    #IDENTIFIERS
    name::String
    type::String

    #PARAMS
    cover::Vector{Int64}
    Q::Int64
    start::Int64
    freq::Int64

    #COSTS
    varq::Float64
    vardq::Float64
    vard::Float64
    fix::Float64
end

struct col
    #QUANTITY RELATED
    q::JuMP.Containers.DenseAxisArray
    u::JuMP.Containers.DenseAxisArray
    v::JuMP.Containers.DenseAxisArray

    #0-1 VARS
    p::JuMP.Containers.DenseAxisArray
    y::JuMP.Containers.DenseAxisArray
    z::JuMP.Containers.DenseAxisArray
end

struct dval
    λ::JuMP.Containers.DenseAxisArray
    δ::JuMP.Containers.DenseAxisArray
end

struct stabilizer
    slackCoeff::Float64
    slackLim::JuMP.Containers.DenseAxisArray
    surpCoeff::Float64
    surpLim::JuMP.Containers.DenseAxisArray
end

struct bound
    idx::NamedTuple
    val::Int64
end

struct node
    #IDENTIFIERS
    parent::UUID
    self::UUID

    #BASE
    base::NamedTuple
    bounds::Vector{bound}

    #DYNAMIC SETS
    columns::Vector{col}

    #SUPPORT STRUCTURE
    stblzr::stabilizer
    status::Vector{String}
end

#CONSTANTS USED THROUGHOUT
const M = 9999999
const rng = MersenneTwister(1234)

#function to generate base data
function base(path::String)
    #READ WORKSHEET
    xf = XLSX.readxlsx(path)

    #DATAFRAME DICT
    data = Dict{Symbol,DataFrame}()

    #TURN SHEETS INTO DATAFRAME
    for sheets in XLSX.sheetnames(xf)
        #TRANSFORM TABLE INTO DATAFRAME
        df = DataFrame(XLSX.gettable(xf[sheets])...)

        #DEFINE THE NAME FROM THE WORKSHEET
        data[Symbol(sheets)] = df
    end

    #PROCESS VERTICES
    V = Dict{Int64,vtx}()

    #iterate over  each row of vertices
    for v in eachrow(data[:vertices])
        #direct input intu struct vtx
        V[v.id] = vtx(
            v.name,v.type,
            v.x,v.y,
            v.MAX,v.MIN,v.START,
            v.h
        )
    end

    #PROCESS VEHICLES
    K = Dict{Int64,veh}()

    #INITIATE INDEXING FOR VEHICLES
    idx = 0

    #iterate over each row of vehicles
    for v in eachrow(data[:vehicles])
        #convert string into vector{Int64}
        v.cover = parse.(Int64,split(v.cover))
        v.loadp = parse.(Int64,split(v.loadp))

        for f in v.loadp
            idx += 1
            K[idx] = veh(
                v.name,v.type,
                v.cover,v.Q,f,v.freq,
                v.varq,v.vardq,v.vard,v.fix
            )
        end
    end

    #EXTRACT T SET
    t_start = data[:periods].start[end] #STARTING MONTH
    t_end = data[:periods].start[end] + data[:periods].T[end] - 1 #FINAL MONTH (CALCULATE BASED ON DURATION)
    T = [i for i in t_start:t_end]

    #EXTRACT DEMANDS
    d = JuMP.Containers.DenseAxisArray{Float64}(undef,collect(keys(V)),T)
    for i in keys(V)
        row = @from x in data[:demands] begin
            @where x.point == i
            @select x
            @collect DataFrame
        end

        for t in T
            d[i,t] = row[2:end-3][t][end]
        end
    end

    #GENERATE DISTANCE MATRIX (FROM V)
    earth = 6378.137
    dist = JuMP.Containers.DenseAxisArray{Float64}(undef,collect(keys(V)),collect(keys(V)))
    for i in keys(V), j in keys(V)
        dist[i,j] = haversine([V[i].x,V[i].y],[V[j].x,V[j].y],earth)
    end

    #GENERATE G MATRIX
    G = JuMP.Containers.DenseAxisArray{Float64}(undef,collect(keys(V)),collect(keys(K)))
    G .= M
    for k in keys(K)
        seed = K[k].cover[findmin([dist[K[k].start,j] for j in K[k].cover if j != K[k].start])[2]] #findmin means the closest to dep point
        for x in K[k].cover
            if x != seed
                G[x,k] = min(dist[K[k].start,x]+dist[x,seed]+dist[seed,K[k].start] , dist[K[k].start,seed]+dist[seed,x]+dist[x,K[k].start]) - (dist[K[k].start,seed] + dist[seed,K[k].start])
            else
                G[x,k] = 0
            end
        end
    end

    #GENERATE DELIVERY COST MATRIX
    deli = JuMP.Containers.DenseAxisArray{Float64}(undef,collect(keys(V)),collect(keys(K)))
    deli .= M
    for k in keys(K),i in K[k].cover
        deli[i,k] = K[k].varq + K[k].vardq * dist[K[k].start,i]
    end

    #DATA GENERATION STATUS
    println("there are $(length(K)) vehicle data points with the composition:")
    println("$(length(collect(keys(filter(p -> last(p).type == "truck",K))))) truck(s)")
    println("$(length(collect(keys(filter(p -> last(p).type == "ship",K))))) ship(s)")
    println("$(length(collect(keys(filter(p -> last(p).type == "train",K))))) train(s)")
    println()
    println("there are $(length(V)) vertex data points with the composition:")
    println("$(length(collect(keys(filter(p -> last(p).type == "source",V))))) source(s)")
    println("$(length(collect(keys(filter(p -> last(p).type == "point",V))))) point(s)")

    return (K=K,V=V,T=T,d=d,dist=dist,G=G,deli=deli)
end

#function to initialize stabilizer
function initStab(dt::NamedTuple,vSlCoeff::Float64,vSurpCoeff::Float64)
    slackCoeff = vSlCoeff
    surpCoeff = vSurpCoeff
    slackLim = JuMP.Containers.DenseAxisArray{Float64}(undef,keys(dt.V),dt.T)
    surpLim = JuMP.Containers.DenseAxisArray{Float64}(undef,keys(dt.V),dt.T)
    for i in keys(dt.V),t in dt.T
        slackLim[i,t] = abs(dt.d[i,t])
        surpLim[i,t] = abs(dt.d[i,t])
    end

    return stabilizer(slackCoeff,slackLim,surpCoeff,surpLim)
end

#function to create root node
function root(dt::NamedTuple;slCoeff::Float64,suCoeff::Float64)
    id = uuid1(rng)
    root = node(
        id,id,
        dt,Vector{bound}(),
        Vector{col}(),
        initStab(dt,slCoeff,suCoeff),
        ["UNVISITED"]
    )

    return root
end

#function to update stabilizer
function updateStab(stabilizers::stabilizer,param::Float64)
    for i in first(stabilizers.slackLim.axes),t in last(stabilizers.slackLim.axes)
        stabilizers.slackLim[i,t] = param * stabilizers.slackLim[i,t]
        if stabilizers.slackLim[i,t] < 1
            stabilizers.slackLim[i,t] = 0
        end
    end

    for i in first(stabilizers.surpLim.axes),t in last(stabilizers.surpLim.axes)
        stabilizers.surpLim[i,t] = param * stabilizers.surpLim[i,t]
        if stabilizers.surpLim[i,t] < 1
            stabilizers.surpLim[i,t] = 0
        end
    end

    return stabilizers
end

function buildMaster(n::node;silent::Bool)
    #COLUMN LABELING
    R = Dict(1:length(n.columns) .=> n.columns)

    #MODEL DECLARATION
    mp = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(Gurobi.Env())))
    if silent
        set_silent(mp)
    end

    #VARIABLE DECLARATION
    θ = @variable(mp, θ[keys(R),keys(n.base.K),n.base.T] >= 0) #LABEL 20
    I = @variable(mp, I[keys(n.base.V),vcat(first(n.base.T)-1,n.base.T)]) #NO BOUNDS (bounded below)

    #SLACK SURPLUS DECLARATION
    @variable(mp, 0 <= slack_I[i=keys(n.base.V),t=n.base.T] <= n.stblzr.slackLim[i,t])
    @variable(mp, 0 <= surp_I[i=keys(n.base.V),t=n.base.T] <= n.stblzr.surpLim[i,t])

    #CONSTRAINT 17 & 18 + STARTING INVENTORY
    @constraint(mp, λ[i=keys(n.base.V),t=n.base.T], I[i,t-1] + sum(R[r].q[i,k,t] * θ[r,k,t] for r in keys(R),k in keys(n.base.K)) + slack_I[i,t] - surp_I[i,t] == n.base.d[i,t] + I[i,t])
    @constraint(mp, [i=keys(n.base.V),t=n.base.T], V[i].MIN <= I[i,t] <= V[i].MAX)
    @constraint(mp, [i=keys(n.base.V)], I[i,first(n.base.T)-1] == V[i].START)

    #CONVEXITY CONSTRAINT LABEL 19
    @constraint(mp, δ[k=keys(n.base.K),t=n.base.T], sum(θ[r,k,t] for r in keys(R)) <= n.base.K[k].freq)

    #OBJECTIVE FUNCTION
    begin
        @objective(mp, Min,

        sum(( n.base.K[k].vard * g(R[r].p[:,k,t]) + sum(deli[i,k] * R[r].v[i,k,t] for i in n.base.K[k].cover) + sum(K[k].fix * R[r].z[i,k,t] for i in n.base.K[k].cover) ) * θ[r,k,t] for r in keys(R),k in keys(n.base.K),t in n.base.T) +
        sum(V[i].h * I[i,t] for i in keys(n.base.V),t in n.base.T) +
        sum(n.stblzr.slackCoeff * slack_I[i,t] for i in keys(n.base.V),t in n.base.T) -
        sum(n.stblzr.surpCoeff * surp_I[i,t] for i in keys(n.base.V),t in n.base.T)
        )
    end
end

export base,root,buildMaster

end
