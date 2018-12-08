mutable struct Pipe<:Component
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,Float64}
    Pipe(s::String)=new(s,Dict{Symbol,Connection}(),Dict{Symbol,Float64}())
end

function setattrs(comp::Pipe;kwargs...)

end

function addconnection(comp::Pipe,c::Connection,id::Symbol)
    port=[:in,:out]
    if id in port
        comp.c[s]=c
    else
        println("wrong port name")
    end
end

function jacobi(comp::Pipe,c::Connection)
    i=comp.conn["in"];o=comp.conn["out"]
    if c==i
        jac[1,1]=1.0;jac[2,1]=i.h.val;jac[2,2]=i.m.val;jac[3,3]=1-comp.dp
    end
    if c==o
        jac[1,1]=-1.0;jac[2,1]=-o.h.val;jac[2,2]=-o.m.val;jac[3,3]=-1.0
    end
end
numberofequations(comp::Pipe)=4
function equations(comp::Pipe)
    i=comp.conn["in"];o=comp.conn["out"]
    res=[]
    vars=[i.m,i.p,i.h,o.m,o.p,o.h]    
    push!(res,i.m.val-o.m.val)
    push!(res,i.m.val*i.h.val-o.m.val*o.h.val)
    push!(res,i.p.val*(1-dp)-o.p.val)
    res
end


export Pipe