@with_kw mutable struct Merge <: AbstractComponent
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,Float64}
    mode=:design
    method=:no
    Merge(s::String)=new(s,Dict{Symbol,Connection}(),Dict{Symbol,Float64}())
end

function (cp::Merge)(;kwargs)
    
end

inlets(comp::Merge)=[comp.conns[:in1] comp.conns[:in2]]
outlets(comp::Merge)=[comp.conns[:out1]]

function addconnection(comp::Merge,s::Symbol,c::Connection)
    port=[:in1,:in2,:out1]
    if s in port
        comp.c[s]=c
    else
        println("wrong input")
    end
end

function jacobi(comp::Merge,c::Connection)
    in1=comp.conns[:in1];in2=comp.conn[:in2];o=comp.conn[:out1]
    jac=zeros(4,3)
    if in1==c
        jac[1,1]=1;jac[2,1]=in1.h.val;jac[2,2]=in1.m.val
        jac[3,3]=1.0
    end
    if in2==c
        jac[1,1]=1.0;jac[2,1]=in2.h.val;jac[2,2]=in2.m.val
        jac[4,3]=1.0
    end
    if o==c
        jac[1,1]=-1.0;jac[2,1]=-o.h.val;jac[2,2]=-o.m.val
        jac[3,3]-1.0;jac[4,3]=-1.0
    end
    return jac
end
function equations(comp::Merge)
    in1=comp.conns[:in1];in2=comp.conns[:in2];o=comp.conns[:out1]
    res=zeros(0)
    push!(res,in1.m.val+in2.m.val-o.m.val)
    push!(res,in1.m.val*in1.h.val+in2.m.val*in2.h.val-o.m.val*o.h.val)
    push!(res,in1.p.val-o.p.val)
    push!(res,in2.p.val-o.p.val)
    res
end

initsource(comp::Merge,c::Connection)=[1.0,500,1.0]
inittarget(comp::Merge,c::Connection)=[1.0,500,1.0]


export Merge,equations,jacobi,setattrs,addconnection






