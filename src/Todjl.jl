module Todjl
import CoolProp
import Parameters:@with_kw
import Base:show,+
export AbstractComponent,Connection,Network,Var,Bus
export addcomponents,addconnections,addbuses,connect,equations,jacobi,setattrs
export addconnection,equations,jacobi,setattr

+(a::Array{Float64,3},b::Array{Float64,3})=vcat(a,b)
include("properties.jl")
include("types.jl")
include("network.jl")
include("components/condensor.jl")
include("components/dearator.jl")
include("components/generalheater.jl")
include("components/source.jl")
include("components/sink.jl")
include("components/pipe.jl")
include("components/valve.jl")
include("components/split.jl")
include("components/merge.jl")
include("components/pump.jl")
include("components/turbine.jl")
include("helper.jl")

end # module
