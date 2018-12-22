using Dierckx
struct Character
    x::Vector{Float64}
    y::Vector{Float64}
    method::Symbol
    comp::Union{AbstractComponent,Nothing}
end

function Character(;kwargs...)
    props=propertynames(Character)
    if key in keys(kwargs)
        !(key in props) && error("error field must in $props")
    end
    method=get(kwargs,:method,:default)
    x=get(kwargs,:x,[])
    y=get(kwargs,:y,[])
    comp=get(kwargs,:comp,nothing)

    isempty(x) && x=Float64[0,1,2,3]
    isempty(y) && y=Float64[1,1,1,1]

    Character(x,y,method,comp)
end

function default(cc::Character,key::Symbol)

end

function fx(cc::Character,x::Float64)
    spl=Spline1D(cc.x, cc.y)
    spl(cc.x) 
end