

"""
setgeometry!(ap::T) where T<:AbstractElement
"""
function setgeometry!(ap::T) where T<:AbstractElement
    𝓖 = ap.𝓖
    for x in 𝓖
        𝒙 = get𝒙(ap,x)
        𝑤 = get𝑤(ap,x)
        x.x = 𝒙[1]
        x.y = 𝒙[2]
        x.z = 𝒙[3]
        x.𝑤 = 𝑤
    end
    if T<:AbstractElement{:Seg2}
        𝐿 = get𝐿(ap)
        for x in 𝓖
            x.𝐿 = 𝐿
        end
    elseif T<:AbstractElement{:Tri3}
        𝐴 = get𝐴(ap)
        for x in 𝓖
            x.𝐴 = 𝐴
        end
    elseif T<:AbstractElement{:Tet4}
        𝑉 = get𝑉(ap)
        for x in 𝓖
            x.𝑉 = 𝑉
        end
    end
end

"""
set_memory_𝝭!(ap::T,ss::Symbol...) where T<:AbstractElement
"""
function set_memory_𝝭!(aps::Vector{T},ss::Symbol...) where T<:AbstractElement
    n = sum(length(ap.𝓒)*length(ap.𝓖) for ap in aps)
    for s in ss
        push!(getfield(aps[1].𝓖[1],:data),s=>(3,zeros(n)))
    end
end

"""
set_memory_𝗠!(aps::Vector{T},ss::Symbol... = keys(aps[1].𝗠)...) where T<:ReproducingKernel
"""
function set_memory_𝗠!(aps::Vector{T},ss::Symbol... = keys(aps[1].𝗠)...) where T<:ReproducingKernel
    set_memory_𝗠!(aps[1],ss...)
end

function set_memory_𝗠!(ap::T,ss::Symbol... = keys(ap[1].𝗠)...) where T<:ReproducingKernel
    n = get𝑛𝒑(ap)
    empty!(ap.𝗠)
    for s in ss
        if s == :∇̃
            n₁ = get𝑛𝒑₁(ap)
            ap.𝗠[s] = SymMat(n₁)
        elseif s ∈ (:∇̃²,:∂∇̃²∂ξ,:∂∇̃²∂η)
            n₂ = get𝑛𝒑₂(ap)
            ap.𝗠[s] = SymMat(n₂)
        else
            ap.𝗠[s] = SymMat(n)
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
    elements = Dict{String,Vector{Int}}()
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
            nₚ = parse(Int,line)
            x = zeros(nₚ)
            y = zeros(nₚ)
            z = zeros(nₚ)
            for i in 1:nₚ
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
            nodes = Node(nodes...)
            readline(fid)
        elseif line == "\$Elements"
            line = readline(fid)
            nₑ = parse(Int,line)
            for i in 1:nₑ
                line = readline(fid)
                entries = split(line," ")
                elmN_ = entries[1]
                elmT_ = entries[2]
                numT_ = entries[3]
                phyT_ = entries[4]
                elmE_ = entries[5]
                l_ = entries[6:end]
                elmNumber = parse(Int,elmN_)
                elmType = parse(Int,elmT_)
                numTag = parse(Int,numT_)
                phyTag = parse(Int,phyT_)
                elmEntary = parse(Int,elmE_)
                nodeList = parse.(Int,l_)
                name = physicalnames[phyTag]
                type = etype[elmType]
                push!(elements[name],Element{type}([nodes[i] for i in nodeList]))
            end
        end
    end
    return elements, nodes
end

function importmsh(filename::String,config::Dict{Any,Any})
    elms, nodes = importmsh(filename)
    if haskey(config,"RegularGrid")
        x = getfield(nodes[1],:data)[:x][2]
        y = getfield(nodes[1],:data)[:y][2]
        z = getfield(nodes[1],:data)[:z][2]
        n = config["RegularGrid"]["n"]
        γ = config["RegularGrid"]["γ"]
        sp = RegularGrid(x,y,z,n=n,γ=γ)
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

    elements = Dict{String,Any}()
    for (name,cfg) in config["elements"]
         # set𝓖
        element_tag = Meta.parse(cfg["𝓒"]["type"])
        integration_type = Meta.parse(cfg["𝓖"]["type"])
        integration_tag = haskey(cfg["𝓖"],"tag") ? elms[cfg["𝓖"]["tag"]] : element_tag
        set𝓖!(elms[integration_tag],integration_type)
        setgeometry!(elms[integration_tag])
        if integration_tag ≠ element_tag
            elms[element_tag*"∩"*integration_tag] = unique!(elms[element_tag]∩elms[integration_tag])
            element_tag = element_tag*"∩"*integration_tag
            set𝓖!(elms[element_tag],elms[integration_tag])
        end

        # set 𝓒
        type = eval(Meta.parse(cfg["type"]))
        if element_type <: Element
            elements[name] = [type([elm.𝓒],elm.𝓖) for elm in elms[element_tag]]
        elseif element_type <: ReproducingKernel
            if haskey(cfg["𝓒"],"type")
                elements[name] = [type(Node[],elm.𝓖) for elm in elms[element_tag]]
                position_type= Meta.parse(cfg["𝓒"]["type"])
                set𝓖!(elms[element_tag],position_type)
            else
                empty!.(elm.𝓖 for elms[element_tag])
            end
            elements[name] = [type(sp(elm,nodes)) elm in elms[element_tag]]
        end

        # set𝓖
        integration_type = Meta.parse(cfg["𝓖"]["type"])
        integration_tag = haskey(cfg["𝓖"],"tag") ? elms[cfg["𝓖"]["tag"]] : element_tag
        set𝓖!(elms[integration_tag],integration_type)
        setgeometry!(elms[integration_tag])
        if integration_tag ≠ element_tag
            set𝓖!(elms[element_tag],elms[integration_tag])
            set𝓖!(elements[name],elms[element_tag])
        else
            set𝓖!(elements[name],elms[integration_tag])
        end

        # set memory
    end
end

function importmsh(filename::String,config::Dict{Any,Any})
    elms, nds = importmsh(filename)
    elements = Dict{String,Any}()
    if haskey(config,"RegularGrid")
        cfg = config["RegularGrid"]
        sp = RegularGrid(nds[:x],nds[:y],nds[:z];n=cfg["n"],γ=cfg["γ"])
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

    nodes = Node(nds...)
    for (name,cfg) in config
        Type = eval(Meta.parse(cfg["type"]))
        elements[name] = [Type([nodes[i] for i in s[2]]) for s in elms[cfg["𝓒"]["tag"]]]
        sp ≠ nothing ? sp(elements[name]) : nothing
        if haskey(cfg,"𝓖")
            QType = Meta.parse(cfg["𝓖"]["type"])
            if haskey(cfg["𝓖"],"tag")
                elms_𝓖 = [Element{s[1]}([nodes[i] for i in s[2]]) for s in elms[cfg["𝓖"]["tag"]]]
                elements[name] = elements[name]∩elms_𝓖
                set𝓖!(elms_𝓖,QType)
                set𝓖!(elements[name],elms_𝓖)
            else
                set𝓖!(elements[name],QType)
            end
        end
    end
    return elements, nodes
end

function importmsh(filename1::String,filename2::String,config::Dict{Any,Any})
    elms, nds_ = importmsh(filename1)
    elms_, nds = importmsh(filename2)
    elements = Dict{String,Any}()
    cfg = config["RegularGrid"]
    sp = RegularGrid(nds[:x],nds[:y],nds[:z];n=cfg["n"],γ=cfg["γ"])
    delete!(config,"RegularGrid")
    nodes = Node(nds...)
    nodes_ = Node(nds_...)
    if haskey(config,"Ωᴳ")
        cfg = config["Ωᴳ"]
        Type = eval(Meta.parse(cfg["type"]))
        elements["Ωᴳ"] = [Type([nodes[i] for i in s[2]]) for s in elms_[cfg["𝓒"]["tag"]]]
        sp(elements["Ωᴳ"])
        QType = Meta.parse(cfg["𝓖"]["type"])
        set𝓖!(elements["Ωᴳ"],QType)
        delete!(config,"Ωᴳ")
    end

    for (name,cfg) in config
        Type = eval(Meta.parse(cfg["type"]))
        elms_𝓖 = [Element{s[1]}([nodes_[i] for i in s[2]]) for s in elms[cfg["𝓖"]["tag"]]]
        if haskey(cfg,"𝓒")
            Type_ = eval(Meta.parse(config[cfg["𝓒"]["tag"]]["type"]))
            if supertype(Type_) ≠ supertype(typeof(elms_𝓖[1]))
                QType_ = Meta.parse(config[cfg["𝓒"]["tag"]]["𝓖"]["type"])
                elms_𝓖_ = [Type_([nodes_[i] for i in s[2]]) for s in elms[cfg["𝓒"]["tag"]]]
                elms_𝓖_ = elms_𝓖_∩elms_𝓖
                set𝓖!(elms_𝓖_,QType_)
                unique!(elms_𝓖_)
                set𝑛ᵢⱼ!(elms_𝓖_)
                elements[name] = [Type(sp(elm,nodes)) for elm in elms_𝓖_]
                name_ = cfg["𝓒"]["tag"]*"∩"*name
                elements[name_] = [Type_(sp(elm,nodes)) for elm in elms_𝓖_] 
                QType = Meta.parse(cfg["𝓖"]["type"])
                set𝓖!(elements[name_],elms_𝓖_)
                set𝓖!(elms_𝓖,QType)
                set𝓖!(elms_𝓖_,elms_𝓖)
                set𝓖!(elements[name],elms_𝓖_)
            else
                QType = Meta.parse(config[cfg["𝓒"]["tag"]]["𝓖"]["type"])
                set𝓖!(elms_𝓖,QType)
                elements[name] = [Type(sp(elm,nodes)) for elm in elms_𝓖]
                QType = Meta.parse(cfg["𝓖"]["type"])
                set𝓖!(elms_𝓖,QType)
                set𝓖!(elements[name],elms_𝓖)
            end
        else
            QType = Meta.parse(cfg["𝓖"]["type"])
            set𝓖!(elms_𝓖,QType)
            set𝑛ᵢⱼ!(elms_𝓖)
            elements[name] = [Type(sp(elm,nodes)) for elm in elms_𝓖]
            set𝓖!(elements[name],elms_𝓖)
        end

        if haskey(cfg["𝓖"],"𝝭")
            ss = Meta.parse.(cfg["𝓖"]["𝝭"])
            Type<:ReproducingKernel ? set_memory_𝗠!(elements[name],ss...) : nothing
            set_memory_𝝭!(elements[name],ss...)
        end
    end
    return elements, nodes
end

function importmsh(filename::String,::Val{:test})
    elems,nodes = importmsh(filename)
    data = Dict([s=>(2,v) for (s,v) in nodes])
    dofs = getboundarydofs2D(elems["Ω"])
    elements = Dict{String,Any}()
    nodes = Node(nodes...)
    gnodes = GNode[]
    elements["∂Ω"] = Vector{Element{:Seg2}}(undef,length(dofs))
    for (dof,i) in dofs
        elements["∂Ω"][i] = Element{:Seg2}([nodes[j] for j in dof])
    end
    set𝓖!(elements["∂Ω"],:SegGI2)
    elements["Γ_λ"] = Element{:Tri3}[]
    elements["Ω"] = DBelement{:Tri3}[]
    elements["Γ"] = DBelement{:Tri3}[]
    elements["Γᵍ"] = DBelement{:Tri3}[]
    haskey(elems,"Γᵗ") ? elements["Γᵗ"] = DBelement{:Tri3}[] : nothing
    for (type,nodeList) in elems["Ω"]
        𝓒 = [GNode((dofs[Set(setdiff(nodeList,i))],i),data) for i in nodeList]
        union!(gnodes,𝓒)
        push!(elements["Ω"],DBelement{:Tri3}(𝓒))
        push!(elements["Γ"],DBelement{:Tri3}(𝓒))
        push!(elements["Γ_λ"],Element{:Tri3}([nodes[i] for i in nodeList]))
        push!(elements["Γᵍ"],DBelement{:Tri3}(𝓒))
        haskey(elems,"Γᵗ") ? push!(elements["Γᵗ"],DBelement{:Tri3}(𝓒)) : nothing
    end
    set𝓖!(elements["Ω"],:TriGI13)
    set𝓖_DB!(elements["Γ"],:SegGI2)
    set𝓖_DB!(elements["Γ_λ"],:SegGI2)
    if haskey(elems,"Γᵗ")
        elms_𝓖 = [Element{type}([nodes[i] for i in nodeList])     for (type,nodeList) in elems["Γᵗ"]]
        elements["Γᵗ"] = elements["Γᵗ"]∩elms_𝓖
        set𝓖!(elms_𝓖,:SegGI2)
        set𝓖!(elements["Γᵗ"],elms_𝓖)
    end

    elms_𝓖 = [Element{type}([nodes[i] for i in nodeList])     for (type,nodeList) in elems["Γᵍ"]]
    elements["Γᵍ"] = elements["Γᵍ"]∩elms_𝓖
    set𝓖!(elms_𝓖,:SegGI2)
    set𝓖!(elements["Γᵍ"],elms_𝓖)

    # elements["Γᵍ"] = DBelement{:Tri3}[]
    # for (type,nodeList) in elems["Γᵍ"]
    #     𝐼 = dofs[Set(nodeList)]
    #     𝓒 = [GNode((0,i),data) for i in nodeList]
    #     push!(𝓒,GNode((𝐼,0),data))
    #     push!(elements["Γᵍ"],DBelement{:Seg2}(𝓒))
    # end
    # set𝓖!(elements["Γᵍ"],:SegGI2)

    set_memory_𝝭!(elements["Ω"],:𝝭,:∂𝝭∂x,:∂𝝭∂y,:∂𝝭∂z)
    haskey(elems,"Γᵗ") ? set_memory_𝝭!(elements["Γᵗ"],:𝝭) : nothing
    set_memory_𝝭!(elements["Γᵍ"],:𝝭)
    set_memory_𝝭!(elements["Γ"],:𝝭)
    return elements, gnodes
end

function getboundarydofs2D(elements::Vector{Tuple{Symbol,Vector{Int}}})
    dofs = Dict{Set{Int},Int}()
    idBoundaries = (Tri3=((1,2),(2,3),(3,1)),)
    n = 0
    for (type,nodeList) in elements
        for bc in idBoundaries[type]
            dof = Set(nodeList[i] for i in bc)
            if ~haskey(dofs,dof)
                n += 1
                dofs[dof] = n
            end
        end
    end
    return dofs
end
