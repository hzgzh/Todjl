@with_kw mutable struct Sink <: AbstractComponent
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,Float64}
    mode=:design
    method=:no
end

inlets(cp::Source)=[cp.conns[:in]]

outlets(cp::Source)=[]

portnames(cp::Source)=[:in]

export Sink