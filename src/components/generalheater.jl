

mutable struct GeneralHeater <: Component
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,Float64}
    GeneralHeater(s::String)=new(s,Dict{Symbol,Connection}(),Dict{Symbol,Float64}())
end

function setattr(comp::GeneralHeater;kwargs...)
    opt=[:temp]
    for key in keys(kwargs)
        if key in opts
            attrs[key]=kwargs[key]
        end
    end
end
function addconnection(c::Connection,id::Symbol)
    port=[:in,:out]
    if id in port
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

export GeneralHeater,setattr,addconnection,euqations,jacobi