@with_kw mutable struct Dearator <: AbstractComponent
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,Float64}
    Dearator(s::String)=new(s,Dict{Symbol,Connection}(),Dict{Symbol,Float64}())
    mode=:design
    method=:steampress
end

inlets(cp::Dearator)=[cp.conns[:hi] cp.conns[:ci] cp.conns[:drain]]

outlets(cp::Dearator)=[cp.conns[:co]]

portnames(cp::Dearator)=[:hi,:ci,:drain,:co]

function equations(cp::Dearator)
    ci=cp.conns[:ci];hi=cp[:hi];di=cp[:drain];co=cp.conns[:co];a=cp.attrs
    vec=[]
    vec+=mat_res(cp)

    res=0
    for (idx,i) in enumerate(inlets(cp))
        res+=i.m.val*i.h.val
    end
    for (idx,o) in enumerate(outlets(cp))
        res-=o.m.val*o.h.val
    end
    vec+=res
    
    vec+=co.h.val-pqh(hi.p.val,0)
    vec=hi.p.val-co.p.val
    return vec
end

function derivatives(cp::Dearator)
    wi=cp.conns[:waterin];si=cp[:steam];di=cp[:drain];wo=cp.conns[:waterout];a=cp.attrs
    der=mat_deriv(cp)

    e_der=zeros(1,4,3)
    for (idx,i) in enumerate(inlets(cp))
        e_der[1,idx,1]=i.h.val;e_der[1,idx,2]=i.m.val
    end
    for (idx,i) in enumerate(outlets(cp))
        e_der(1,idx+length(inlets(cp)),1)=-i.h.val;e_der(1,idx+length(inlets(cp)),2)=-i.m.val
    end
    der+=e_der
    e_der=zeros(1,4,3)
    e_der[1,4,2]=1.0;e_der[1,1,3]=-dhdp(hi.p.val,0)
    der+=e_der
    p_der=zeros(1,4,3)
    p_der[1,1,3]=1.0;p_der[1,4,3]=-1.0
    der+=e_der
    
    return der
end

export Dearator,equations,jacobi,setattrs,addconnection