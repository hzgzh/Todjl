@with_kw mutable struct Source <: AbstractComponent
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,VarProp}
    Source(s::String)=new(s,Dict{Symbol,Connection}(),Dict{Symbol,Float64}())
end

inlets(cp::Source)=[]

outlets(cp::Source)=[cp.conns[:out]]

portnames(cp::Source)=[:out]

export Source,equations,jacobi,setattrs,addconnection