function Operator(::Val{:DataPool})
    data = Dict{Symbol,Any}()
    :physicalnodes
        Symbol=>Vector{Float64}
        :x
        :y
        :z
        :s₁
        :s₂
        :s₃
    :elements
        String=>Vector{Approximator}
    :parametricnodes
        String=>(Symbol=>Vector{Float64})
    :config

    return elements, physicalnodes, parametricnodes


function Operator(::Val{:msh};meshfree::Bool=false)
    ntype = Dict(1=>Node,2=>Node,3=>Node,15=>Node)
    if meshfree
        etype = Dict(1=>SegN,2=>TriN,3=>QuadN,15=>PoiN)
        qtype = Dict(1=>:SegGI5,2=>:TriGI13,3=>:QuadGI5,15=>:PoiGI1)
        push!(data,:basisfunction=>:Linear1D,:kerneltype=>:□,:kernelfunction=>:CubicSpline)
        push!(data,:stype=>[:𝝭])
        push!(data,:spatialpartition=>:RegularGrid,:nᵣ=>1,:γᵣ=>1)
    else
        etype = Dict(1=>Seg2,2=>Tri3,3=>Quad,15=>Poi1)
        qtype = Dict(1=>:SegGI2,2=>:TriGI3,3=>:QuadGI2,15=>:PoiGI1)
    end
    push!(data,:meshfree=>meshfree)
    push!(data,:etype=>etype)
    push!(data,:ntype=>ntype)
    push!(data,:qtype=>qtype)
    return Operator(Val(:msh),data)
end

function (op::Operator{:msh})(filename::String)
    fid = open(filename,"r")
    readline(fid)
    line = readline(fid)
    v_,f_,d_ = split(line," ")
    version = parse(Float64,v_)
    filetype = parse(Int,f_)
    datasize = parse(Int,d_)
    readline(fid)
    if version == 4.1
        import_msh_4(fid,op)
    elseif version == 2.2
        import_msh_2(fid,op)
    else
        println("Version does not match!")
    end
    return op.elements
end

function import_msh_4(fid::IO,op::Operator{:msh})
    aps = Dict{String,Vector{Approximator}}()
    phy = Dict{Int,String}()
    for line in eachline(fid)
        if line == "\$PhysicalNames"
            numPhysicalNames = parse(Int,readline(fid))
            push!(op.data,)
            for i in 1:numPhysicalNames
                line = readline(fid)
                d_,p_,name = split(line," ")
                dimension = parse(Int,d_)
                physicalTag = parse(Int,p_)

                phy[physicalTag] = name
            end
            readline(fid)
        elseif line == "\$Nodes"
            line = readline(fid)
            numE_,numN_,minN_,maxN_ = split(line," ")
            numEntityBlocks = parse(Int,numE_)
            numNodes = parse(Int,numN_)
            minNodeTag = parse(Int,minN_)
            maxNodeTag = parse(Int,maxN_)
            x = Vector{PhysicalNode}()
            for i in 1:numEntityBlocks
                line = readline(fid)
                entityD_,entityE_,para_,numN_ = split(line," ")
                # entityD_
            end
        end
    end
    return aps
end

function import_msh_2(fid::IO,op::Operator{:msh})
    for line in eachline(fid)
        if line == "\$PhysicalNames"
            numPhysicalNames = parse(Int,readline(fid))
            push!(op,:physicaldata=>Dict{Int,Any}())
            push!(op,:elements=>Dict{String,Any}())
            for i in 1:numPhysicalNames
                line = readline(fid)
                d_,p_,n_ = split(line," ")
                dimension = parse(Int,d_)
                physicalTag = parse(Int,p_)
                name = strip(n_,'\"')

                op.physicaldata[physicalTag] = Dict(:name=>name,:nₑ=>0,:nᵢ=>0,:data=>Dict{Symbol,Vector{Float64}}())
            end
            readline(fid)
        elseif line == "\$Nodes"
            line = readline(fid)
            nₚ = parse(Int,line)
            push!(op,:nodes=>Dict(:x=>zeros(nₚ),:y=>zeros(nₚ),:z=>zeros(nₚ)))
            push!(op,:nₚ=>nₚ)
            for i in 1:nₚ
                line = readline(fid)
                t_,x_,y_,z_ = split(line," ")
                tag = parse(Int,t_)
                op.nodes[:x][i] = parse(Float64,x_)
                op.nodes[:y][i] = parse(Float64,y_)
                op.nodes[:z][i] = parse(Float64,z_)
            end
            readline(fid)
        elseif line == "\$Elements"
            line = readline(fid)
            nₑ = parse(Int,line)
            push!(op,:nₑ=>nₑ)
            if op.meshfree
                nₘ = 0
                if op.spatialpartition == :RegularGrid
                    sp = RegularGrid(op.nodes,n=op.nᵣ,γ=op.γᵣ)
                    for c in sp.cells
                        nₘ = max(length(c),nₘ)
                    end
                end
                n = length(get𝒑(Val(op.basisfunction),(0.0,0.0,0.0)))
                push!(op,:temp𝗠=>Dict{Symbol,SymMat}(),:temp𝝭=>Dict{Symbol,Vector{Float64}}())
                for s in op.stype
                    push!(op.data[:temp𝗠],s=>SymMat(n))
                    push!(op.data[:temp𝝭],s=>zeros(nₘ))
                end
                𝗠 = op.temp𝗠
                𝝭 = op.temp𝝭
                𝒑 = op.basisfunction
                𝑠 = op.kerneltype
                𝜙 = op.kernelfunction
                for i in 1:nₑ
                    line = readline(fid)
                    elmN_,elmT_,numT_,phyT_,elmE_,l_... = split(line," ")
                    elmNumber = parse(Int,elmN_)
                    elmType = parse(Int,elmT_)
                    numTag = parse(Int,numT_)
                    phyTag = parse(Int,phyT_)
                    elmEntary = parse(Int,elmE_)
                    nodeList = parse.(Int,l_)
                    name = op.physicaldata[phyTag][:name]
                    op.physicaldata[phyTag][:nₑ] += 1
                    for i in nodeList
                        x = (op.nodes[:x][i],op.nodes[:y][i],op.nodes[:z][i])
                        union!(nodeList,sp(x))
                    end
                    𝓒 = [Node(i,op.nodes) for i in nodeList]
                    quadraturepoints = QuadratureRule[op.qtype[elmType]]
                    data = op.physicaldata[phyTag][:data]
                    𝓖 = op.ntype[elmType][]
                    for ξ in quadraturepoints
                        op.physicaldata[phyTag][:nᵢ] += 1
                        haskey(data,:w) ? push!(data[:w],ξ[1]) : push!(data,:w=>[ξ[1]])
                        haskey(data,:ξ) ? push!(data[:ξ],ξ[2]) : push!(data,:ξ=>[ξ[2]])
                        if length(ξ)≥3
                            haskey(data,:η) ? push!(data[:η],ξ[3]) : push!(data,:η=>[ξ[3]])
                        end
                        if length(ξ)≥4
                            haskey(data,:γ) ? push!(data[:γ],ξ[4]) : push!(data,:γ=>[ξ[4]])
                        end
                        n = length(data[:w])
                        push!(𝓖,op.ntype[elmType](n,data))
                    end
                    haskey(op.elements,name) ? push!(op.elements[name],op.etype[elmType](𝓒,𝓖,𝗠,𝝭,𝒑,𝑠,𝜙)) : push!(op.elements,name=>op.etype[elmType][op.etype[elmType](𝓒,𝓖,𝗠,𝝭,𝒑,𝑠,𝜙)])
                end
            else
                for i in 1:nₑ
                    line = readline(fid)
                    elmN_,elmT_,numT_,phyT_,elmE_,l_... = split(line," ")
                    elmNumber = parse(Int,elmN_)
                    elmType = parse(Int,elmT_)
                    numTag = parse(Int,numT_)
                    phyTag = parse(Int,phyT_)
                    elmEntary = parse(Int,elmE_)
                    nodeList = parse.(Int,l_)
                    name = op.physicaldata[phyTag][:name]
                    op.physicaldata[phyTag][:nₑ] += 1
                    𝓒 = [Node(i,op.nodes) for i in nodeList]
                    quadraturepoints = QuadratureRule[op.qtype[elmType]]
                    data = op.physicaldata[phyTag][:data]
                    𝓖 = op.ntype[elmType][]
                    for ξ in quadraturepoints
                        op.physicaldata[phyTag][:nᵢ] += 1
                        haskey(data,:w) ? push!(data[:w],ξ[1]) : push!(data,:w=>[ξ[1]])
                        haskey(data,:ξ) ? push!(data[:ξ],ξ[2]) : push!(data,:ξ=>[ξ[2]])
                        if length(ξ)≥3
                            haskey(data,:η) ? push!(data[:η],ξ[3]) : push!(data,:η=>[ξ[3]])
                        end
                        if length(ξ)≥4
                            haskey(data,:γ) ? push!(data[:γ],ξ[4]) : push!(data,:γ=>[ξ[4]])
                        end
                        n = length(data[:w])
                        push!(𝓖,op.ntype[elmType](n,data))
                    end
                    haskey(op.elements,name) ? push!(op.elements[name],op.etype[elmType](𝓒,𝓖)) : push!(op.elements,name=>op.etype[elmType][op.etype[elmType](𝓒,𝓖)])
                end
            end
        end
    end
end

## Quadrature Points
const QuadratureRule = Dict(
:PoiGI1 => ((1.0,-1.0),),
:SegGI1 => ((2.0,0.0),),
:SegGI2 =>
(
    (1.0,-0.5773502691896257645091487805),
    (1.0, 0.5773502691896257645091487805)
),
:SegGI3 =>
(
    (0.555555555555555555555555555556,-0.774596669241483377035853079957),
    (0.88888888888888888888888888889,0.0),
    (0.555555555555555555555555555556, 0.774596669241483377035853079957)
),
:SegGI4 =>
(
    (0.347854845137453857373063949222,-0.861136311594052575223946488893),
    (0.652145154862546142626936050778,-0.339981043584856264802665759103),
    (0.652145154862546142626936050778, 0.339981043584856264802665759103),
    (0.347854845137453857373063949222, 0.861136311594052575223946488893)
),
:SegGI5 =>
(
    (0.23692688505618908751426404072,-0.906179845938663992797626878299),
    (0.47862867049936646804129151484,-0.5384693101056830910363144207),
    (0.568888888888888888888888888889, 0.0),
    (0.47862867049936646804129151484, 0.5384693101056830910363144207),
    (0.23692688505618908751426404072, 0.906179845938663992797626878299)
),
:SegGI6 =>
(
    (0.171324492379170345040296142173,-0.932469514203152027812301554494),
    (0.360761573048138607569833513838,-0.6612093864662645136613995950),
    (0.46791393457269104738987034399,-0.238619186083196908630501721681),
    (0.46791393457269104738987034399, 0.238619186083196908630501721681),
    (0.360761573048138607569833513838, 0.66120938646626451366139959502),
    (0.171324492379170345040296142173, 0.932469514203152027812301554494)
),
:SegGI7 =>
(
    (0.129484966168869693270611432679,-0.949107912342758524526189684048),
    (0.27970539148927666790146777142,-0.741531185599394439863864773281),
    (0.38183005050511894495036977549,-0.405845151377397166906606412077),
    (0.417959183673469387755102040816, 0.0),
    (0.381830050505118944950369775489, 0.405845151377397166906606412077),
    (0.279705391489276667901467771424, 0.741531185599394439863864773281),
    (0.129484966168869693270611432679, 0.949107912342758524526189684048)
),
:SegGI8 =>
(
    (0.10122853629037625915253135431,-0.96028985649753623168356086857),
    (0.22238103445337447054435599443,-0.796666477413626739591553936476),
    (0.313706645877887287337962201987,-0.525532409916328985817739049189),
    (0.36268378337836198296515044928,-0.18343464249564980493947614236),
    (0.362683783378361982965150449277, 0.18343464249564980493947614236),
    (0.31370664587788728733796220199, 0.525532409916328985817739049189),
    (0.222381034453374470544355994426, 0.796666477413626739591553936476),
    (0.10122853629037625915253135431, 0.96028985649753623168356086857)
),
:SegGI9 =>
(
    (0.0812743883615744119718921581105,-0.968160239507626089835576202904),
    (0.180648160694857404058472031243,-0.83603110732663579429942978807),
    (0.260610696402935462318742869419,-0.613371432700590397308702039342),
    (0.31234707704000284006863040658,-0.32425342340380892903853801464),
    (0.330239355001259763164525069287, 0.0),
    (0.31234707704000284006863040658, 0.32425342340380892903853801464),
    (0.260610696402935462318742869419, 0.613371432700590397308702039342),
    (0.180648160694857404058472031243, 0.83603110732663579429942978807),
    (0.081274388361574411971892158111, 0.968160239507626089835576202904)
),
:SegGI10 =>
(
    (0.066671344308688137593568809893,-0.973906528517171720077964012085),
    (0.149451349150580593145776339658,-0.865063366688984510732096688424),
    (0.219086362515982043995534934228,-0.679409568299024406234327365115),
    (0.26926671930999635509122692157,-0.433395394129247190799265943166),
    (0.295524224714752870173892994651,-0.14887433898163121088482600113),
    (0.295524224714752870173892994651, 0.14887433898163121088482600113),
    (0.26926671930999635509122692157, 0.433395394129247190799265943166),
    (0.219086362515982043995534934228, 0.679409568299024406234327365115),
    (0.149451349150580593145776339658, 0.865063366688984510732096688424),
    (0.066671344308688137593568809893, 0.973906528517171720077964012085)
),
:TriGI3 =>
(
    (1/3,2/3,1/6),
    (1/3,1/6,2/3),
    (1/3,1/6,1/6)
),
:TriGI4 =>
(
    (-0.562500000000000,0.333333333333333,0.333333333333333),
    (0.520833333333333,0.600000000000000,0.200000000000000),
    (0.520833333333333,0.200000000000000,0.600000000000000),
    (0.520833333333333,0.200000000000000,0.200000000000000)
),
:TriGI6 =>
(
    (0.223381589678011,0.108103018168070,0.445948490915965),
    (0.223381589678011,0.445948490915965,0.108103018168070),
    (0.223381589678011,0.445948490915965,0.445948490915965),
    (0.109951743655322,0.816847572980459,0.091576213509771),
    (0.109951743655322,0.091576213509771,0.816847572980459),
    (0.109951743655322,0.091576213509771,0.091576213509771)
),
:TriGI7 =>
(
    (0.125939180544800,0.101286507323500,0.101286507323500),
    (0.125939180544800,0.797426985353100,0.101286507323500),
    (0.125939180544800,0.101286507323500,0.797426985353100),
    (0.132394152788500,0.470142064105100,0.059715871789800),
    (0.132394152788500,0.470142064105100,0.470142064105100),
    (0.132394152788500,0.059715871789800,0.470142064105100),
    (0.225000000000000,0.333333333333300,0.333333333333300)
),
:TriGI12 =>
(
    (0.116786275726379,0.501426509658179,0.249286745170910),
    (0.116786275726379,0.249286745170910,0.501426509658179),
    (0.116786275726379,0.249286745170910,0.249286745170910),
    (0.050844906370207,0.873821971016996,0.063089014491502),
    (0.050844906370207,0.063089014491502,0.873821971016996),
    (0.050844906370207,0.063089014491502,0.063089014491502),
    (0.082851075618374,0.053145049844817,0.310352451033784),
    (0.082851075618374,0.053145049844817,0.636502499121399),
    (0.082851075618374,0.310352451033784,0.053145049844817),
    (0.082851075618374,0.310352451033784,0.636502499121399),
    (0.082851075618374,0.636502499121399,0.053145049844817),
    (0.082851075618374,0.636502499121399,0.310352451033784)
),
:TriGI13 =>
(
    (0.053347235608800,0.065130102902200,0.065130102902200),
    (0.053347235608800,0.869739794195600,0.065130102902200),
    (0.053347235608800,0.065130102902200,0.869739794195600),
    (0.077113760890300,0.312865496004900,0.048690315425300),
    (0.077113760890300,0.638444188569800,0.312865496004900),
    (0.077113760890300,0.048690315425300,0.638444188569800),
    (0.077113760890300,0.638444188569800,0.048690315425300),
    (0.077113760890300,0.312865496004900,0.638444188569800),
    (0.077113760890300,0.048690315425300,0.312865496004900),
    (0.175615257433200,0.260345966079000,0.260345966079000),
    (0.175615257433200,0.479308067841900,0.260345966079000),
    (0.175615257433200,0.260345966079000,0.479308067841900),
    (-0.14957004446770,0.333333333333300,0.333333333333300)
),
:TriGI16 =>
(
    (0.144315607677787,0.333333333333333,0.333333333333333),
    (0.095091634267285,0.081414823414554,0.459292588292723),
    (0.095091634267285,0.459292588292723,0.081414823414554),
    (0.095091634267285,0.459292588292723,0.459292588292723),
    (0.103217370534718,0.658861384496480,0.170569307751760),
    (0.103217370534718,0.170569307751760,0.658861384496480),
    (0.103217370534718,0.170569307751760,0.170569307751760),
    (0.032458497623198,0.898905543365938,0.050547228317031),
    (0.032458497623198,0.050547228317031,0.898905543365938),
    (0.032458497623198,0.050547228317031,0.050547228317031),
    (0.027230314174435,0.008394777409958,0.263112829634638),
    (0.027230314174435,0.008394777409958,0.728492392955404),
    (0.027230314174435,0.263112829634638,0.008394777409958),
    (0.027230314174435,0.263112829634638,0.728492392955404),
    (0.027230314174435,0.728492392955404,0.008394777409958),
    (0.027230314174435,0.728492392955404,0.263112829634638)
),
:QuadGI1 =>
(
    (2.0,0.0,0.0),
),
:QuadGI2 =>
(
    (1.0,-0.5773502691896258,-0.5773502691896258),
    (1.0, 0.5773502691896258,-0.5773502691896258),
    (1.0, 0.5773502691896258, 0.5773502691896258),
    (1.0,-0.5773502691896258, 0.5773502691896258)
)
)
