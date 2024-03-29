
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
            x = zeros(nₚ)
            y = zeros(nₚ)
            z = zeros(nₚ)
            for i in 1:nₚ
                line = readline(fid)
                x,y = split(line," ")
                x[i] = parse(Float64,x)
                y[i] = parse(Float64,y)
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
                push!(elements["Γ"],Seg2((i1,i2),x,y,z))
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
                push!(elements["Ω"],Tri3((i1,i2,i3),x,y,z))
            end
        end
    end
    return elements, entities
end
