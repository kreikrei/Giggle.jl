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
const GUROBI_ENV = Gurobi.Env()
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

function initialize(base::NamedTuple,valueSlackCoeff::Float64,valueSurplusCoeff::Float64)
    slackCoeff = valueSlackCoeff
    surpCoeff = valueSurplusCoeff
    slackLim = JuMP.Containers.DenseAxisArray{Float64}(undef,keys(base.V),base.T)
    surpLim = JuMP.Containers.DenseAxisArray{Float64}(undef,keys(base.V),base.T)
    for i in keys(base.V),t in base.T
        slackLim[i,t] = abs(base.d[i,t])
        surpLim[i,t] = abs(base.d[i,t])
    end

    return stabilizer(slackCoeff,slackLim,surpCoeff,surpLim)
end

#function to create root node
function root(dt::NamedTuple;slCoeff,suCoeff)
    id = uuid1(rng)
    root = node(
        id,id,
        dt,Vector{bound}(),
        Vector{col}(),
        initialize(dt,slCoeff,suCoeff),
        ["UNVISITED"]
    )

    return root
end

export base,root

end
