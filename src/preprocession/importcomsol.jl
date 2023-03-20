
function importcomsol(filename::String)
    fid = open(filename,"r")
    points = Point[]
    elements = Dict(["Ω"=>Tri3[],"Γ"=>Seg2[]])
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
                push!(points,Point(x,y,0.0))
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
        end
    end
    return elements,points
end
