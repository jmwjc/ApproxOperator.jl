
function importcomsol(filename::String)
    fid = open(filename,"r")
    points = Point[]
    elements = Dict(["Ω"=>Tri3[],"Γ"=>Seg2[]])
    entities = Dict(["Ω"=>Int[],"Γ"=>Int[]])
    nₚ = 0
    for line in eachline(fid)
        if occursin("# number of mesh vertices",line)
            nₚ,s = split(line," ")
            nₚ = parse(Int,nₚ)
        elseif line == "# Mesh vertex coordinates"
            for i in 1:nₚ
                line = readline(fid)
                x,y = split(line," ")
                x = parse(Float64,x)
                y = parse(Float64,y)
                push!(points,Point(i,x,y,0.0))
            end
        elseif line == "# Type #1"
            for i in 1:5
                readline(fid)
            end
            line = readline(fid)
            n,s = split(line," ")
            n = parse(Int,n)
            line = readline(fid)
            for i in 1:n
                line = readline(fid)
                i1,i2 = split(line," ")
                i1 = parse(Int,i1)+1
                i2 = parse(Int,i2)+1
                push!(elements["Γ"],Seg2((points[i1],points[i2])))
            end
            for i in 1:3
                readline(fid)
            end
            for i in 1:n
                line = readline(fid)
                j = parse(Int,line)
                push!(entities["Γ"],j)
            end
        elseif line == "# Type #2"
            for i in 1:5
                readline(fid)
            end
            line = readline(fid)
            n,s = split(line," ")
            n = parse(Int,n)
            line = readline(fid)
            for i in 1:n
                line = readline(fid)
                i1,i2,i3 = split(line," ")
                i1 = parse(Int,i1)+1
                i2 = parse(Int,i2)+1
                i3 = parse(Int,i3)+1
                push!(elements["Ω"],Tri3((points[i1],points[i2],points[i3]),(Seg2((points[i2],points[i3])),Seg2((points[i3],points[i1])),Seg2((points[i1],points[i2])))))
            end
        end
    end
    return elements,points,entities
end

function importcomsol_fem(filename::String)
    elms,nds = importcomsol(filename)
    nₚ = length(nds)
    nodes = Node{(:𝐼,),1}[]
    data = Dict([:x=>(1,zeros(nₚ)),:y=>(1,zeros(nₚ)),:z=>(1,zeros(nₚ))])
    for (i,p) in enumerate(nds)
        node = Node{(:𝐼,),1}((i,),data)
        node.x = p.x
        node.y = p.y
        node.z = p.z
        push!(nodes,node)
    end

    elements = Dict(["Ω"=>Element{:Tri3}[],"Γ"=>Element{:Seg2}[]])

    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    g = 0
    G = 0
    s = 0
    ng = 3
    gauss_scheme = :TriGI3
    nₑ = length(elms["Ω"])

    scheme = quadraturerule(gauss_scheme)
    data_𝓖 = Dict([
        :ξ=>(1,scheme[:ξ]),
        :η=>(1,scheme[:η]),
        :w=>(1,scheme[:w]),
        :x=>(2,zeros(ng*nₑ)),
        :y=>(2,zeros(ng*nₑ)),
        :z=>(2,zeros(ng*nₑ)),
        :𝑤=>(2,zeros(ng*nₑ)),
        :𝝭=>(4,zeros(ng*nₑ*3)),
        :∂𝝭∂x=>(4,zeros(ng*nₑ*3)),
        :∂𝝭∂y=>(4,zeros(ng*nₑ*3)),
    ])
    for (C,a) in enumerate(elms["Ω"])
        element = Element{:Tri3}((c,3,𝓒),(g,ng,𝓖))
        for v in a.vertices
            i = v.i
            push!(𝓒,nodes[i])
        end
        c += 3

        𝐴 = get𝐴(a)
        for i in 1:ng
            G += 1
            x = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((i,G,C,s),data_𝓖)
            ξ = x.ξ
            η = x.η
            x_,y_,z_ = a(ξ,η)
            x.x = x_
            x.y = y_
            x.z = z_
            x.𝑤 = 𝐴*x.w
            push!(𝓖,x)
            s += 3
        end
        g += ng
        push!(elements["Ω"],element)
        
    end
    
    
    𝓒 = Node{(:𝐼,),1}[]
    𝓖 = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}[]
    c = 0
    g = 0
    G = 0
    s = 0
    ng = 2 
    gauss_scheme = :SegGI2
    nₑ = length(elms["Γ"])
    scheme = quadraturerule(gauss_scheme)

    data_𝓖 = Dict([
        :ξ=>(1,scheme[:ξ]),
        :w=>(1,scheme[:w]),
        :x=>(2,zeros(ng*nₑ)),
        :y=>(2,zeros(ng*nₑ)),
        :z=>(2,zeros(ng*nₑ)),
        :𝑤=>(2,zeros(ng*nₑ)),
        :𝝭=>(4,zeros(ng*nₑ*2)),
        :∂𝝭∂x=>(4,zeros(ng*nₑ*2)),
        :∂𝝭∂y=>(4,zeros(ng*nₑ*2)),
    ])
    for (C,a) in enumerate(elms["Γ"])
        element = Element{:Seg2}((c,2,𝓒),(g,ng,𝓖))
        for v in a.vertices
            i = v.i
            push!(𝓒,nodes[i])
        end
        c += 2
       
        𝐿 = get𝐿(a)
        for i in 1:ng
            G += 1
            x = Node{(:𝑔,:𝐺,:𝐶,:𝑠),4}((i,G,C,s),data_𝓖)
            ξ = x.ξ
            x_,y_,z_ = a(ξ)
            x.x = x_
            x.y = y_
            x.z = z_
            x.𝑤 = 𝐿*x.w
            push!(𝓖,x)
            s += 2
        end
        g += ng
        push!(elements["Γ"],element)
    end
 
    return elements,nodes
end
