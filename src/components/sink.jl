@with_kw mutable struct Sink <: AbstractComponent
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,Float64}
    mode=:design
    method=:no
end

function (cp::Sink)(;kwargs)
    
end

export Sink