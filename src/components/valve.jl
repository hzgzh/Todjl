mutable struct Valve<:Component
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,Float64}
    Valve(s::String)=new(s,Dict{Symbol,Connection}(),Dict{Symbol,Float64}())
end

function setattrs(comp::Valve;kwargs...)

end

function addconnection(comp::Valve,port::Symbol,c::Connection)
    ports=[:in,:out]
    if port in port
        comp.c[port]=c
    else
        println("wrong port name")
    end
end

function jacobi(comp::Valve,c::Connection)
    i=comp.conn["in"];o=comp.conn["out"]
    if c==i
        jac[1,1]=1.0;jac[2,1]=i.h.val;jac[2,2]=i.m.val;jac[3,3]=1-comp.dp
    end
    if c==o
        jac[1,1]=-1.0;jac[2,1]=-o.h.val;jac[2,2]=-o.m.val;jac[3,3]=-1.0
    end
end
numberofequations(comp::Valve)=4
function equations(comp::Valve)
    i=comp.conn["in"];o=comp.conn["out"]
    res=[]
    vars=[i.m,i.p,i.h,o.m,o.p,o.h]    
    push!(res,i.m.val-o.m.val)
    push!(res,i.m.val*i.h.val-o.m.val*o.h.val)
    push!(res,i.p.val*(1-dp)-o.p.val)
    res
end


export Valve