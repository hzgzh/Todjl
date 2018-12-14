
@with_kw mutable struct GeneralHeater <: AbstractComponent
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,Float64}
    mode=:design
    method=:no
    GeneralHeater(s::String)=new(s,Dict{Symbol,Connection}(),Dict{Symbol,Float64}())
end

function (cp::GeneralHeater)(;kwargs)
    
end

function setattrs(comp::GeneralHeater;kwargs...)
    opts=[:temp];attrs=comp.attrs
    for key in keys(kwargs)
        if key in opts
            attrs[key]=kwargs[key]
        end
    end
end

function addconnection(comp::GeneralHeater,port::Symbol,c::Connection)
    ports=[:in,:out]
    if port in ports
        comp.conns[port]=c
    else
        println("wrong input")
    end
end

inlets(comp::GeneralHeater)=[comp.conns[:in]]
outlets(comp::GeneralHeater)=[comp.conns[:out]]

function equations(comp::GeneralHeater)
    i=comp.conns[:in];o=comp.conns[:out];attrs=comp.attrs
    res=zeros(0)
    push!(res,i.m.val-o.m.val)
    push!(res,o.h.val-pth(o.p.val,attrs[:temp]))
    push!(res,i.p.val-o.p.val)
    res
end

function jacobi(comp::GeneralHeater,c::Connection)
    i=comp.conns[:in];o=comp.conns[:out];attrs=comp.attrs
    jac=zeros(3,3)
    if c==i
        jac[1,1]=1;jac[3,3]=1;
    end

    if c==o
        jac[1,1]=-1;
        jac[2,2]=1;jac[2,3]=-dhdp(o.p.val,attrs[:temp])
        jac[3,3]=-1
    end
    return jac
end

export GeneralHeater,equations,jacobi,setattrs,addconnection