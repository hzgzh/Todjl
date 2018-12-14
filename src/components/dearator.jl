

@with_kw mutable struct Dearator <: AbstractComponent
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,Float64}
    Dearator(s::String)=new(s,Dict{Symbol,Connection}(),Dict{Symbol,Float64}())
    mode=:design
    method=:steampress
end

function (cp::Dearator)(;kwargs)
    
end

function setattr(comp::Dearator,sym::Symbol,val::Float64)
    attrs=comp.attrs
    if sym==:temp
        attrs[sym]=val
    end
end

inlets(comp::Dearator)=[comp.conns[:waterin] comp.conns[:steam] comp.conns[:drain]]
outlets(comp::Dearator)=[comp.conns[:waterout]]

function addconnection(c::Connection,id::String)
    port=[:steam,:waterin,:drain,:waterout]
    if s in port
        comp.conns[s]=c
    else
        println("wrong input")
    end
end

function equations(comp::Dearator)
    wi=comp.conns[:waterin];si=comp[:steam];di=comp[:drain];wo=comp.conns[:waterout];a=comp.attrs
    res=zeros(0)
    push!(res,wi.m.val+si.m.val+di.m.val-wo.m.val)
    push!(res,wi.m.val*wi.h.val+si.m.val*si.h.val+di.m.val*di.h.val-wo.m.val*wo.h.val)
    push!(res,wo.h.val-pqh(si.p.val,0))
    push!(res,si.p.val-wo.p.val)
    return res
end

conns(comp::Dearator)=[:waterin,:steam,:drain,:waterout]

function jacobi(comp::Dearator,c::Connection)
    wi=comp.conns[:waterin];si=comp[:steam];di=comp[:drain];wo=comp.conns[:waterout];a=comp.attrs
    jac=zeros(4,3)
    if c==wi
        jac[1,1]=1;jac[2,1]=wi.h.val;;jac[2,2]=wi.m.val
    end

    if c==si
        jac[1,1]=1;
        jac[2,1]=si.h.val;jac[2,2]=si.m.val
        jac[3,3]=-sat_dhdp(si.p.val,0)
        jac[4,3]=1
    end

    if c==di
        jac[1,1]=1;jac[2,1]=di.h.val;jac[2,2]=di.m.val;
    end

    if c==wo
        jac[1,1]=-1;jac[2,1]=-wo.h.val;jac[2,2]=-wo.m.val
        jac[3,2]=1;jac[4,3]=-1
    end
    return jac
end

export Dearator,equations,jacobi,setattrs,addconnection