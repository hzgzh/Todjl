export GeneralHeater,setattr,addconnection,euqations,jacobi

mutable struct GeneralHeater
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,Float64}
    calcmode::Symbol
    mode::Symbol
end

function setattr(comp::GeneralHeater,sym::Symbol,val::Float64)
    if sym==:temp
        attrs[sym]=val
    end
end
function addconnection(c::Connection,id::String)
    port=[:in,:out]
    if s in port
        comp.conns[s]=c
    else
        println("wrong input")
    end
end
function equations(comp::GeneralHeater)
    i=comp.conns[:in];o=comp.conns[:out];a=comp.attrs
    res=[]
    push!(res,i.m.val-o.m.val)
    push!(res,o.h.val-pth(o.p.val,a[:temp]))
    push!(res,i.p.val-o.p.val)
end

function jacobi(comp::GeneralHeater,c::Connection)
    i=comp.conns[:in];o=comp.conns[:out]
    jac=zeros(3,6)
    if c==i
        jac[1,1]=1;jac[3,3]=1;
    end

    if c==j
        jac[1,4]=-1;
        jac[2,5]=1;jac[2][6]=-dhdp(o.p.val,a[:temp])
        jac[3,6]=-1
    end
    jac
end

