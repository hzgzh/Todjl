@with_kw mutable struct Source <: AbstractComponent
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,CVar}
    mode=:design
    method=:no
    Source(s::String)=new(s,Dict{Symbol,Connection}(),Dict{Symbol,Float64}())
end

function (cp::Source)(;kwargs)
    
end

inlets(cp::Source)=[]

outlets(cp::Source)=[cp.conns[:out]]

portnames(cp::Source)=[:out]

export Source,equations,jacobi,setattrs,addconnection