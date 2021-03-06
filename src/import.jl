

"""
setgeometry!(ap::T) where T<:AbstractElement
"""
function setgeometry!(ap::T) where T<:AbstractElement
    π = ap.π
    for x in π
        π = getπ(ap,x)
        π€ = getπ€(ap,x)
        x.x = π[1]
        x.y = π[2]
        x.z = π[3]
        x.π€ = π€
    end
end

"""
set_memory_π­!(ap::T,ss::Symbol...) where T<:AbstractElement
"""
function set_memory_π­!(aps::Vector{T},ss::Symbol...) where T<:AbstractElement
    n = sum(length(ap.π)*length(ap.π) for ap in aps)
    for s in ss
        push!(getfield(aps[1].π[1],:data),s=>(3,zeros(n)))
    end
end

"""
set_memory_π !(aps::Vector{T},ss::Symbol... = keys(aps[1].π )...) where T<:ReproducingKernel
"""
function set_memory_π !(aps::Vector{T},ss::Symbol... = keys(aps[1].π )...) where T<:ReproducingKernel
    set_memory_π !(aps[1],ss...)
end

function set_memory_π !(ap::T,ss::Symbol... = keys(ap[1].π )...) where T<:ReproducingKernel
    n = getππ(ap)
    nβ = getππβ(ap)
    nβ = getππβ(ap)
    empty!(ap.π )
    for s in ss
        if s == :βΜ
            ap.π [s] = SymMat(nβ)
        elseif s β (:βΜΒ²,:ββΜΒ²βΞΎ,:ββΜΒ²βΞ·)
            ap.π [s] = SymMat(nβ)
        else
            ap.π [s] = SymMat(n)
        end
    end
end

## ---------------- msh ---------------
function importmsh(filename::String)
    fid = open(filename,"r")
    readline(fid)
    line = readline(fid)
    v_,f_,d_ = split(line," ")
    version = parse(Float64,v_)
    filetype = parse(Int,f_)
    datasize = parse(Int,d_)
    readline(fid)
    if version == 4.1
        elements,nodes = import_msh_4(fid)
    elseif version == 2.2
        elements,nodes = import_msh_2(fid)
    else
        println("Version does not match!")
    end
    return elements, nodes
end

function import_msh_4(fid::IO) end

function import_msh_2(fid::IO)
    etype = Dict(1=>:Seg2,2=>:Tri3,3=>:Quad,15=>:Poi1)
    nodes = Dict{Symbol,Vector{Float64}}()
    elements = Dict{String,Vector{Tuple{Symbol,Vector{Int}}}}()
    physicalnames = Dict{Int,String}()
    for line in eachline(fid)
        if line == "\$PhysicalNames"
            numPhysicalNames = parse(Int,readline(fid))
            physicalnames=>Dict{Int,String}()
            for i in 1:numPhysicalNames
                line = readline(fid)
                d_,p_,n_ = split(line," ")
                dimension = parse(Int,d_)
                physicalTag = parse(Int,p_)
                name = strip(n_,'\"')
                physicalnames[physicalTag] = name
                elements[name] = Vector{Tuple{Symbol,Vector{Int}}}()
            end
            readline(fid)
        elseif line == "\$Nodes"
            line = readline(fid)
            nβ = parse(Int,line)
            x = zeros(nβ)
            y = zeros(nβ)
            z = zeros(nβ)
            for i in 1:nβ
                line = readline(fid)
                t_,x_,y_,z_ = split(line," ")
                tag = parse(Int,t_)
                x[i] = parse(Float64,x_)
                y[i] = parse(Float64,y_)
                z[i] = parse(Float64,z_)
            end
            nodes[:x] = x
            nodes[:y] = y
            nodes[:z] = z
            readline(fid)
        elseif line == "\$Elements"
            line = readline(fid)
            nβ = parse(Int,line)
            for i in 1:nβ
                line = readline(fid)
                elmN_,elmT_,numT_,phyT_,elmE_,l_... = split(line," ")
                elmNumber = parse(Int,elmN_)
                elmType = parse(Int,elmT_)
                numTag = parse(Int,numT_)
                phyTag = parse(Int,phyT_)
                elmEntary = parse(Int,elmE_)
                nodeList = parse.(Int,l_)
                name = physicalnames[phyTag]
                type = etype[elmType]
                push!(elements[name],(type,nodeList))
            end
            return elements, nodes
        end
    end
end

function importmsh(filename::String,config::Dict{Any,Any})
    elms, nodes = importmsh(filename)
    elements = Dict{String,Any}()
    if haskey(config,"RegularGrid")
        cfg = config["RegularGrid"]
        sp = RegularGrid(nodes[:x],nodes[:y],nodes[:z];n=cfg["n"],Ξ³=cfg["Ξ³"])
        delete!(config,"RegularGrid")
    else
        sp = nothing
    end
    if haskey(config,"IndependentDofs")
        for (k,v) in config["IndependentDofs"]
            dofs = Set{Int}()
            for (type,nodeList) in elms[v]
                union!(dofs,Set(nodeList))
            end
            elms[k] = [(:Poi1,[dof]) for dof in dofs]
        end
        delete!(config,"IndependentDofs")
    end

    nodes = Node(nodes...)
    for (name,cfg) in config
        Type = eval(Meta.parse(cfg["type"]))
        if Type <: ReproducingKernel
            π  = Dict{Symbol,SymMat}()
            elements[name] = [Type([nodes[i] for i in s[2]],π ) for s in elms[cfg["π"]["tag"]]]
        else
            elements[name] = [Type([nodes[i] for i in s[2]]) for s in elms[cfg["π"]["tag"]]]
        end
        sp β  nothing ? sp(elements[name]) : nothing
        if haskey(cfg,"π")
            QType = Meta.parse(cfg["π"]["type"])
            if haskey(cfg["π"],"tag")
                elms_π = [Element{s[1]}([nodes[i] for i in s[2]]) for s in elms[cfg["π"]["tag"]]]
                elements[name] = elements[name]β©elms_π
                setπ!(elms_π,QType)
                setπ!(elements[name],elms_π)
            else
                setπ!(elements[name],QType)
            end
            if haskey(cfg["π"],"π­")
                ss = Meta.parse.(cfg["π"]["π­"])
                Type<:ReproducingKernel ? set_memory_π !(elements[name],ss...) : nothing
                set_memory_π­!(elements[name],ss...)
            end
        end
    end
    return elements, nodes
end
