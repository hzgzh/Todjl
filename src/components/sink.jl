mutable struct Sink <: Component
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,Float64}
end


export Sink