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
        if i != j
            dist[i,j] = haversine([V[i].x,V[i].y],[V[j].x,V[j].y],earth)
        else
            dist[i,j] = M
        end
    end

    #GENERATE G MATRIX
    G = JuMP.Containers.DenseAxisArray{Float64}(undef,collect(keys(V)),collect(keys(K)))
    G .= M
    for k in keys(K)
        seed = K[k].cover[findmin([dist[K[k].start,j] for j in K[k].cover])[2]]
        for x in K[k].cover
            if x != seed
                G[x,k] = (min(dist[K[k].start,x]+dist[x,seed]+dist[seed,K[k].start] , dist[K[k].start,seed]+dist[seed,x]+dist[x,K[k].start]) - (dist[K[k].start,seed] + dist[seed,K[k].start]))
            else
                G[x,k] = 0
            end
        end
    end

    #GENERATE DELIVERY COST MATRIX
    deli = JuMP.Containers.DenseAxisArray{Float64}(undef,collect(keys(V)),collect(keys(K)))
    deli .= M
    for k in keys(K)
        for i in K[k].cover
            deli[i,k] = K[k].varq + K[k].vardq * dist[K[k].start,i]
        end
    end

    #DATA GENERATION STATUS
    trucks = sort(collect(keys(filter(p -> last(p).type == "truck",K))))
    ships = sort(collect(keys(filter(p -> last(p).type == "ship",K))))
    trains = sort(collect(keys(filter(p -> last(p).type == "train",K))))
    println("there are $(length(K)) vehicle data points with the composition:")
    println("$(length(trucks)) truck(s) for index $(first(trucks)) to $(last(trucks))")
    println("$(length(ships)) ship(s) for index $(first(ships)) to $(last(ships))")
    println("$(length(trains)) train(s) for index $(first(trains)) to $(last(trains))")
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

#===============MASTER FUNCTION SET=====#
function master(n::node;silent::Bool,env::Gurobi.Env)
    #build and bound
    mp = buildMaster(n;silent=silent,env=env)

    #solve the mp
    optimize!(mp)

    return mp
end

function buildMaster(n::node;silent::Bool,env::Gurobi.Env)
    #COLUMN LABELING
    R = Dict(1:length(n.columns) .=> n.columns)

    #MODEL DECLARATION
    mp = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(env)))
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
    @constraint(
        mp, λ[i=keys(n.base.V),t=n.base.T],
        I[i,t-1] + sum(R[r].q[i,k,t] * θ[r,k,t] for r in keys(R),k in keys(n.base.K)) + slack_I[i,t] - surp_I[i,t] == n.base.d[i,t] + I[i,t]
    )
    @constraint(mp, [i=keys(n.base.V),t=n.base.T], n.base.V[i].MIN <= I[i,t] <= n.base.V[i].MAX)
    @constraint(mp, [i=keys(n.base.V)], I[i,first(n.base.T)-1] == n.base.V[i].START)

    #CONVEXITY CONSTRAINT LABEL 19
    @constraint(mp, δ[k=keys(n.base.K),t=n.base.T], sum(θ[r,k,t] for r in keys(R)) <= n.base.K[k].freq)

    #BOUND GENERATOR
    #empty

    #OBJECTIVE FUNCTION
    begin
        @objective(mp, Min,
            sum(
                (
                    n.base.K[k].vard * g(R[r].y[:,k,t];n=n,k=k) +
                    sum(n.base.deli[i,k] * R[r].u[i,k,t] for i in n.base.K[k].cover) +
                    n.base.K[k].fix * R[r].z[n.base.K[k].start,k,t]
                ) * θ[r,k,t] for r in keys(R),k in keys(n.base.K),t in n.base.T
            ) +
            sum(n.base.V[i].h * I[i,t] for i in keys(n.base.V),t in n.base.T) +
            sum(n.stblzr.slackCoeff * slack_I[i,t] for i in keys(n.base.V),t in n.base.T) -
            sum(n.stblzr.surpCoeff * surp_I[i,t] for i in keys(n.base.V),t in n.base.T)
        )
    end

    return mp
end

function getDual(mp::Model)
    return duals = dval(dual.(mp.obj_dict[:λ]),dual.(mp.obj_dict[:δ]))
end
#===============MASTER FUNCTION SET=====#

#===============SUB FUNCTION SET=====#
function sub(n::node,duals::dval;silent::Bool,env::Gurobi.Env)
    #build and bound
    sp = buildSub(n,duals;silent=silent,env=env)

    #solve the sp
    optimize!(sp)

    return sp
end

function buildSub(n::node,duals::dval;silent::Bool,env::Gurobi.Env)
    #MODEL DECLARATION
    sp = Model(optimizer_with_attributes(() -> Gurobi.Optimizer(env)))
    if silent
        set_silent(sp)
    end

    #VARIABLE DECLARATION
    q = @variable(sp, q[keys(n.base.V),keys(n.base.K),n.base.T],Int)
    u = @variable(sp, u[keys(n.base.V),keys(n.base.K),n.base.T] >= 0,Int)
    v = @variable(sp, v[keys(n.base.V),keys(n.base.K),n.base.T] >= 0,Int)
    y = @variable(sp, y[keys(n.base.V),keys(n.base.K),n.base.T], Bin)
    z = @variable(sp, z[keys(n.base.V),keys(n.base.K),n.base.T], Bin)
    p = @variable(sp, p[keys(n.base.V),keys(n.base.K),n.base.T], Bin)

    #CONSTRAINTS
    @constraint(sp, [k=keys(n.base.K),t=n.base.T], sum(q[i,k,t] for i in n.base.K[k].cover) == 0)
    @constraint(sp, [k=keys(n.base.K),i=n.base.K[k].cover,t=n.base.T], q[i,k,t] == u[i,k,t] - v[i,k,t])
    @constraint(sp, [k=keys(n.base.K),i=n.base.K[k].cover,t=n.base.T], u[i,k,t] <= n.base.K[k].Q * y[i,k,t])
    @constraint(sp, [k=keys(n.base.K),i=n.base.K[k].cover,t=n.base.T], v[i,k,t] <= n.base.K[k].Q * z[i,k,t])
    @constraint( #everything outside the possible starting point is zero
        sp, [k=keys(n.base.K),t=n.base.T],
        sum(z[i,k,t] for i in filter(x -> x != n.base.K[k].start,n.base.K[k].cover)) == 0
    )
    @constraint( #the starting point can be zero or one
        sp, [k=keys(n.base.K),t=n.base.T],
        z[n.base.K[k].start,k,t] <= 1
    )
    @constraint(sp, [k=keys(n.base.K),i=n.base.K[k].cover,t=n.base.T], p[i,k,t] == y[i,k,t] + z[i,k,t])

    #BOUND GENERATOR SUBPROBLEM
    #empty

    #OBJECTIVE FUNCTION
    begin
        @objective(sp,Min,
            sum(
                sum(n.base.K[k].vard * n.base.G[i,k] * y[i,k,t] for i in n.base.K[k].cover) +
                sum(n.base.deli[i,k] * u[i,k,t] for i in n.base.K[k].cover) +
                n.base.K[k].fix * z[n.base.K[k].start,k,t] for k in keys(n.base.K),t in n.base.T
            ) -
            sum(sum(q[i,k,t] * duals.λ[i,t] for i in n.base.K[k].cover) for k in keys(n.base.K),t in n.base.T) -
            sum(duals.δ[k,t] for k in keys(n.base.K),t in n.base.T)
        )
    end

    return sp
end

function getCol(sp::Model)
    #EXTRACT VARIABLES
    q = value.(sp.obj_dict[:q])
    u = value.(sp.obj_dict[:u])
    v = value.(sp.obj_dict[:v])
    p = value.(sp.obj_dict[:p])
    y = value.(sp.obj_dict[:y])
    z = value.(sp.obj_dict[:z])

    return col(q,u,v,p,y,z)
end

function realPrice(n::node,duals::dval;column::col)
    #CALCULATE REAL PRICING
    begin
        price = (
            sum(
                n.base.K[k].vard * g(column.y[:,k,t];n=n,k=k) +
                sum(n.base.deli[i,k] * column.u[i,k,t] for i in n.base.K[k].cover) +
                n.base.K[k].fix * column.z[n.base.K[k].start,k,t] for k in keys(n.base.K),t in n.base.T
            ) -
            sum(sum(column.q[i,k,t] * duals.λ[i,t] for i in n.base.K[k].cover) for k in keys(n.base.K),t in n.base.T) -
            sum(duals.δ[k,t] for k in keys(n.base.K),t in n.base.T)
        )
    end

    return price
end
#===============SUB FUNCTION SET=====#

#===============SUPPORT ESTIMATOR FOR TSP===#
function g(y::JuMP.Containers.DenseAxisArray;n::node,k::Int64)
    x = twoOpt([i for i in first(y.axes) if y[i] == 1.0];n=n,k=k)
    return x
end #CLEARED

function twoOpt(p::Vector{Int64};n::node,k::Int64) #INPUT IS A COLLECTION OF POINTS [15 22 23 24] ORDER IS ARBITRARY
    if length(p) > 0
        #DETERMINE START AND NODE
        start = n.base.K[k].start #passing k stop di sini
        nodes = vcat(start,p)

        #DO NEAREST NEIGHBORHOOD construction
        nnPath = nn(start,nodes;n=n)

        #RETURN ITS IMPROVEMENT
        return trav(improve(nnPath;n=n);n=n)
    else
        return 0
    end
end #CLEARED

function nn(start::Int64,nodes::Vector{Int64};n::node)
    tour = [start]

    unvisited = deepcopy(nodes)
    current_index = 1
    terminate = false

    while !terminate
        current_position = tour[current_index]

        #REMOVE CURRENT NODES AFTER FIRST MOVE
        to_remove = findfirst(x -> x==current_position,unvisited)
        splice!(unvisited,to_remove)

        #if no more unvisited, terminate
        if length(unvisited) == 0
            push!(tour,first(tour))
            terminate = true
        else
            terpilih = selectNn(current_position,unvisited;n=n)
            push!(tour,terpilih)
        end

        #UPDATE INDEX
        current_index += 1
    end

    return tour
end #CLEARED

function selectNn(current::Int64,unvisited::Vector{Int64};n::node)
    chosen = unvisited[1]

    for j in unvisited
        if n.base.dist[current,chosen] > n.base.dist[current,j]
            chosen = j #tuker kalo jaraknya lebih kecil
        end
    end

    return chosen
end #CLEARED

function improve(p::Vector{Int64};n::node) #INPUT: NNPATH i.e. [15 27 28 29 15]
    a = view(p,2:length(p)-1,1)

    for i in 1:length(a) - 1
        for k in i+1:length(a) - 1
            best = trav(p;n=n)
            reverse!(view(a,i:k,1))
            if trav(p;n=n) >= best
                reverse!(view(a,i:k,1))
            end
        end
    end

    return p
end #CLEARED

function trav(p::Vector{Int64};n::node) #INPUT: NNPATH i.e. [15 27 28 29 15]
    return sum(n.base.dist[p[i],p[i+1]] for i in 1:length(p)-1)
end
#===============SUPPORT ESTIMATOR FOR TSP===#

function check(model::Model)
    return value(sum(model.obj_dict[:slack_I])) + value(sum(model.obj_dict[:surp_I]))
end #CLEARED

function colGen(n::node,maxCG::Float64;silent::Bool,env::Gurobi.Env,track::Bool)
    #INITIALIZE
    terminate = false
    iter = 0

    while !terminate
        #GENERATE MASTER PROBLEM
        mp = master(n;silent=silent,env=env)

        #NODE PROCESSING STATUS CHECK
        if has_values(mp) && has_duals(mp)
            if track #tracking status
                println("obj: $(objective_value(mp))")
            end

            #EXTRACT DUAL
            duals = getDual(mp)

            #GENERATE SUBPROBLEM
            sp = sub(n,duals;silent=silent,env=env)

            #EXTRACT COLS & PRICE
            cols = getCol(sp)
            price = realPrice(n,duals;column=cols)
            slacksurp = check(mp)

            if track #tracking status
                println(price)
                println("nilai slack surp $slacksurp")
            end

            #NEGATIVE COLUMN CHECK
            if isapprox(price, 0, atol = 1e-8) || price > 0
                if isapprox(slacksurp, 0, atol = 1e-8)
                    terminate = true
                    push!(n.status,"EVALUATED")
                    if track
                        println("EVALUATED")
                    end
                else
                    updateStab(n.stblzr,0.4)
                    push!(n.status,"STABILIZED")
                    if track
                        println("STABILIZED")
                    end
                end
            else
                push!(n.columns,cols)
                push!(n.status,"ADD_COLUMN")
                if track
                    println("ADD_COLUMN")
                end
            end

            #UPDATE ITER
            iter += 1

            #MAXIMUM ITER
            if iter >= maxCG
                terminate = true
                pop!(n.columns)
                push!(n.status,"EVALUATED-TIME OUT")
                if track
                    println("EVALUATED-TIME OUT")
                end
            end
        else
            terminate = true
            push!(n.status,"NO_SOLUTION")
            if track
                println("NO_SOLUTION")
            end
        end
    end

    if n.status[end] == "EVALUATED" || n.status[end] == "EVALUATED-TIME OUT"
        println("NODE $(n.self) FINISHED.")
    else
        println("NODE $(n.self) FAILED.")
    end
end

export base
export root

export master
export buildMaster
export getDual

export sub
export buildSub
export getCol

export realPrice

export colGen

end
