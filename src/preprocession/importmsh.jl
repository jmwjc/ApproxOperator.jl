
"""
importmsh
"""
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
        elements,entities = import_msh_4(fid)
    elseif version == 2.2
        elements,entities = import_msh_2(fid)
    else
        println("Version does not match!")
    end
    return elements, entities
end

function import_msh_4(fid::IO) end

function import_msh_2(fid::IO)
    etype = Dict(1=>:Seg2,2=>:Tri3,3=>:Quad4,8=>:Seg3,9=>:Tri6,15=>:Point)
    elements = Dict{String,Any}()
    entities = Dict{String,Any}()
    physicalnames = Dict{Int,String}()
    x = Float64[]
    y = Float64[]
    z = Float64[]
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
            end
            readline(fid)
        elseif line == "\$Nodes"
            line = readline(fid)
            nₚ = parse(Int,line)
            resize!(x,nₚ)
            resize!(y,nₚ)
            resize!(z,nₚ)
            for i in 1:nₚ
                line = readline(fid)
                i_,x_,y_,z_ = split(line," ")
                i_ = parse(Int,i_)
                x[i_] = parse(Float64,x_)
                y[i_] = parse(Float64,y_)
                z[i_] = parse(Float64,z_)
            end
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
                type = eval(etype[elmType])
                if ~haskey(elements,name)
                    elements[name] = type[]
                    entities[name] = Int[]
                end
                push!(elements[name],type(Tuple(nodeList),x,y,z))
                push!(entities[name],elmEntary)
            end
        end
    end
    return elements, entities
end
