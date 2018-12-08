module Todjl
import CoolProp

export addconnection,equations,jacobi,setattr
export pt,pth,pts,ptv,pht,phs,phv,phq,psh,pst,psv,pqh,pqs,tp,thp,tsp,tqp,ths

include("struct.jl")
include("properties.jl")
include("components/condensor.jl")
include("components/dearator.jl")
include("components/generalheater.jl")
include("components/source.jl")
include("components/sink.jl")
include("components/pipevalve.jl")
include("components/split.jl")
include("components/merge.jl")
include("components/turbine.jl")


end # module
