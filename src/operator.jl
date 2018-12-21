import Base:+

+(a::Array{AbstractFloat,3},b::Array{AbstractFloat,3})=vcat(a,b)

function +(a::Vector{T} where T<:Any,b...) 
    a=convert(Vector{Float64},a)
    for x in b[1]
        push!(a,convert(Float64),x))
    end
    return a
end