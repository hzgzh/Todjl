mutable struct Split <: Component
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,Float64}
end

function addconnection(comp::Split,port::Symbol,c::Connection)
    ports=[:in,:out1,:out2]
    if port in ports
        comp.c[port]=c
    else
        println("wrong input")
    end
end
function jacobi(comp::Split,c::Connection)
    in=comp.conns[:in];o1=comp.conns[:out1];o2=comp.conns[:out2]
    jac=zeros(4,3)
    if c==in
        jac[1,1]=1.0;jac[2,1]=in.h.val;jac[2,2]=in.m.val
        jac[3,3]=1.0;jac[4,3]=1.0
    end
    if c==o1
        jac[1,1]=-1.0;jac[2,1]=-o1.h.val;jac[2,2]=-o1.m.val
        jac[3,3]=-1.0
    end
    if c==o2
        jac[1,1]=-1.0;jac[2,1]=-o2.h.val;jac[2,2]=-o2.m.val
        jac[4,3]=-1.0
    end
end
function equations(comp::Split)
    in=comp.conns[:in];o1=comp.conns[:out1];o2=comp.conns[:out2]
    res=[]
    push!(res,in.m.val-o1.m.val-o2.m.val)
    push!(res,in.m*in.h.val-o1.m.val*o1.h.val-o2.m.val*o2.h.val)
    push!(res,in.p.val-o1.p.val)
    push!(res,in.p.val-o2.p.val)
end

export Split