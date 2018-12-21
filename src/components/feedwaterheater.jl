

@with_kw mutable struct FeedwaterHeater<:AbstractComponent
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,Var}
    FeedwaterHeater(s::String)=new(s,Dict{Symbol,Connection}(),Dict{Symbol,Float64}())
    mode=:design
    method=:td
end

inlets(cp::FeedwaterHeater)=
    haskey(cp.conns,:prevdrain) ? [cp.conns[:hi] cp.conns[:drain] cp.conns[:ci]] : [cp.conns[:hi] cp.conns[:ci] ]

outlets(cp::FeedwaterHeater)=[cp.conns[:ho] cp.conns[:co]]

portnames(cp::FeedwaterHeater)=[:hi,:drain,:ci,:ho,:co]

function equations(cp::FeedwaterHeater)
    conns=cp.conns;attrs=cp.attrs
    ci=conns[:ci];co=conns[:co];hi=conns[:hi];
    co=conns[:ho];di=haskey(conns,:drain) ? conns[:drain] : nothing;tu=attrs[:tu];td=attrs[:td]
    vec=[]
    vec+=mass_res(cp)
    
    res=0.
    for i in inlets(cp)
        res+=i.m.val*i.h.val
    end
    for o in inlets(cp)
        res-=o.m.val*o.h.val
    end
    vec+=res

    if cp.attrs[:tu].val_set
        vec+=ttd_u_func(cp)
    end

    if cp.attrs[:tl].val_set
        vec+=ttd_l_func(cp)
    end

    vec+=hi.p.val-ho.p.val
    vec+=ci.p.val-co.p.val
 
    return res
end

function derivatives(cp::FeedwaterHeater)
    conns=cp.conns;attrs=cp.attrs
    ci=conns[:ci];co=conns[:co];hi=conns[:hi];
    ho=conns[:ho];di=haskey(conns,:drain) ? conns[:drain] : nothing;tu=attrs[:tu];td=attrs[:td]
    der=mass_deriv(cp)
        
    e_derive=zeros(1,length(inlets(cp)),3)
    for (idx,i) in enumerate(inlets(cp))
        e_der[1,idx,1]=i.h.val;e_derive[1,idx,2]=i.m.val
    end

    for (idx,o) in enumerate(outlets(cp))
        e_der[1,length(inlets(cp))+idx,1]=-o.h.val;e_der[1,length(inlets(cp))+idx,2]=-o.m.val
    end
    der+=e_der

    if cp.attrs[:td].val_set
        der+=ttd_u_deriv(cp)
    end

    if cp.attrs[:tl].val_set
        der+=ttd_l_deriv(cp)
    end

    p_der=zeros(2,length(inlets(cp)),3)
    if haskey(cp.conns,:drain)
        p_der[1,1,3]=1;p_der[1,2,3]=1;p_dere[1,4,3]=-1
        p_der[1,3,3]=1;p_der[1,5,3]=-1
    else
        p_der[1,1,3]=1;p_der[1,3,3]=-1
        p_der[1,2,3]=1;p_der[1,4,3]=-1
    end
    
    der+=p_deriv

    return der
end

function ttd_u_func(cp::FeedwaterHeater)
    co=cp.conns[:co];hi=cp.conns[:hi];tu=cp.attrs[:tu].val
    return co.h.val-pth(co.p.val,pt(hi.p.val)+tu)
end

function ttd_u_deriv(cp::FeedwaterHeater)
    deriv=zeros(1,length(inlets(cp)),3)
    deriv[1,1,3]=derive(cp,ttd_u_func,1,:p)
    if haskey(cp.conns,:drain)
        deriv[1,5,2]=1;deriv[1,5,3]=derive(cp,ttd_u_func,5,:p)
    else
        deriv[1,4,2]=1;deriv[1,4,3]=derive(cp,ttd_u_func,4,:p)
    end
    return deriv
end

function ttd_l_func(cp::FeedwaterHeater)
    ho=cp.conns[:ho];ci=cp.conns[:ci];td=cp.attrs[:tl].val
    return ho.h.val-pth(ho.p.val,pht(ci.p.val,ci.h.val)+td)
end

function ttd_l_deriv(cp::FeedwaterHeater)
    deriv=zeros(1,length(inlets(cp)),3)
    if haskey(cp.conns,:prevdrain)
        deriv[1,4,2]=1;deriv[1,4,3]=derive(cp,ttd_l_func,4,:p)
        deriv[1,3,2]=derive(cp,ttd_l_func,3,:h);deriv[1,3,3]=deriv(cp,ttd_l_func,3,:p)
    else
        deriv[1,3,2]=1;deriv[1,3,3]=derive(cp.ttd_l_func,3,:p)
        deriv[1,2,2]=derive(cp,ttd_l_func,2,:h);deriv[1,2,3]=derive(cp,ttd_l_func,2,:p)
    end
    return deriv
end

function mass_res(cp::FeedwaterHeater)
    conns=cp.conns;attrs=cp.attrs
    ci=conns[:ci];co=conns[:co];hi=conns[:hi];ho=conns[:ho];
    res=[]
    if haskey(cp.conns,:drain)
        di=conns[:drain]
        res+=hi.m.val+di.m.va-ho.m.val
        res+=ci.m.val-co.m.val
    else
        res+=hi.m.val-ho.m.val
        res+=ci.m.val-co.m.val
    end
    return res
end

function mass_deriv(cp::FeedwaterHeater)
    deriv=nothing
    if haskey(cp.conns,:prevdrain)
        deriv=zeros(2,5,3)
        deriv[1,1,1]=1;deriv[1,2,1]=1;derive[1,3,1]=1
        deriv[2,4,1]=-1;derive[2,5,1]=-1
    else
        deriv=zeros(2,4,3)
        deriv[1,1,1]=1;deriv[1,2,1]=1
        deriv[2,3,1]=-1;derive[2,4,1]=-1
    end
    return deriv
end

export FeedwaterHeater,equations,jacobi,setattrs,addconnection