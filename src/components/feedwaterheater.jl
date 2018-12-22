

@with_kw mutable struct FeedwaterHeater<:AbstractHeater
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,Var}
    FeedwaterHeater(s::String)=new(s,Dict{Symbol,Connection}(),Dict{Symbol,Float64}())
end

inlets(cp::FeedwaterHeater)=
    haskey(cp.conns,:prevdrain) ? [cp.conns[:hi] cp.conns[:drain] cp.conns[:ci]] : [cp.conns[:hi] cp.conns[:ci] ]

outlets(cp::FeedwaterHeater)=[cp.conns[:ho] cp.conns[:co]]
ports(cp::FeedwaterHeater)=[:hi,:drain,:ci,:ho,:co]
attrs(cp::FeedwaterHeater)=[:tu,:tl,:pr1,:pr2]
designmode(cp::FeedwaterHeater)=[:tu,:tl,:pr1,:pr2]
offdesignmode(cp::FeedwaterHeater)=[]

function aux_equations(cp::FeedWaterHeater)

end

function aux_derivatives(cp::FeedWaterHeater)

end

export FeedwaterHeater,equations,jacobi,setattrs,addconnection