

mutable struct Dearator <: Component
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,Float64}
    Dearator(s::String)=new(s,Dict{Symbol,Connection}(),Dict{Symbol,Float64}())
end

function setattr(comp::Dearator,sym::Symbol,val::Float64)
    if sym==:temp
        attrs[sym]=val
    end
end

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
    res=[]
    push!(res,wi.m.val+si.m.val+di.m.val-wo.m.val)
    push!(res,wi.m.val*wi.h.val+si.m.val*si.h.val+di.m.val*di.h.val-wo.m.val*wo.h.val)
    push!(res,wo.h.val-pqh(si.p.val,0))
    push!(res,si.p.val-wo.p.val)
end

conns(comp::Dearator)=[:waterin,:steam,:drain,:waterout]

function jacobi(comp::Dearator,c::Connection)
    wi=comp.conns[:waterin];si=comp[:steam];di=comp[:drain];wo=comp.conns[:waterout];a=comp.attrs
    jac=zeros(4,12)
    if c==wi
        jac[1,1]=1;jac[2,1]=wi.h.val;;jac[2,2]=wi.m.val
    end

    if c==si
        jac[1,4]=1;
        jac[2,4]=si.h.val;jac[2,5]=si.m.val
        jac[3,6]=-sat_dhdp(si.p.val,0)
        jac[4,6]=1
    end

    if c==di
        jac[1,7]=1;jac[2,7]=di.h.val;jac[2,8]=di.m.val;
    end

    if c==wo
        jac[1,10]=-1;jac[2,10]=-wo.h.val;jac[2,11]=-wo.m.val
        jac[3,11]=1;jac[4,12]=-1
    end
    jac
end

export Dearator