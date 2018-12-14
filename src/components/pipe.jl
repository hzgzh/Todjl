@with_kw mutable struct Pipe<:AbstractComponent
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,Float64}
    mode=:design
    method=:no
    Pipe(s::String)=new(s,Dict{Symbol,Connection}(),Dict{Symbol,Float64}())
end

function (cp::Pipe)(;kwargs)
    
end

function setattrs(comp::Pipe;kwargs...)

end

inlets(comp::Pipe)=[comp.conns[:in]]
outlets(comp::Pipe)=[comp.conns[:out]]

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
    jac=zeros(3,3)
    if c==i
        jac[1,1]=1.0;jac[2,1]=i.h.val;jac[2,2]=i.m.val;jac[3,3]=1-comp.dp
    end
    if c==o
        jac[1,1]=-1.0;jac[2,1]=-o.h.val;jac[2,2]=-o.m.val;jac[3,3]=-1.0
    end
    return jac
end
numberofequations(comp::Pipe)=4
function equations(comp::Pipe)
    i=comp.conn[:in];o=comp.conn[:out]
    res=zeros(0)
     
    push!(res,i.m.val-o.m.val)
    push!(res,i.m.val*i.h.val-o.m.val*o.h.val)
    push!(res,i.p.val*(1-dp)-o.p.val)
    res
end

export Pipe,equations,jacobi,setattrs,addconnection