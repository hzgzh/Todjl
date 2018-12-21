@with_kw mutable struct Condensor <: AbstractComponent
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,CVar}
    mode=:design
    method=:fixpress

    Condensor(s::String)=new(s,Dict{Symbol,Connection}(),Dict{Symbol,Float64}())
end

function inlets(cp::Condensor)
    if haskey(conns,:drain)
        return [cp.conns[:hi] cp.conns[:drain] cp.conns[:ci]]
    else
        return [cp.conns[:hi] cp.conns[:ci]]
    end
end

outlets(cp::Condensor)=[cp.conns[:ho] cp.conns[:co]]

portnames(cp::Condensor)=[:hi,:drain,:ci,:ho,:co]



function equations(cp::Condensor)
    hi=cp.conns[:hi];ho=cp.conns[:ho]
    vec=[]
    vec+==mass_res(cp)
    res=0
    for (idx,i) in enumerate(inlets(cp))
        res+=i.m.val*i.h.val
    end
    for (idx,o) in enumerate(outlets(cp))
        res-=o.m.val+o.h.val
    end
    vec+=res
    vec+=pqh(hi.p.val,0)-ho.h.val
    vec+=hi.p.val-ho.p.val
    vec+=ci.p.val-co.p.val
    return res
end



function derivatives(cp::Condensor)
    hi=cp.conns[:hi];ho=cp.conns[:ho];ci=cp.conns[:ci];co=cp.conns[:co]
    di=haskey(cp.conns,:drain) ? cp.conns[:drain] : nothing
    der=mass_deriv(cp)

    e_der=zeros(1,length(inlets(cp)),3)
    for (idx,i) in enumerate(inlets(cp))
        e_der[1,idx,1]=i.h.val;e_derive[1,idx,2]=i.m.val
    end

    for (idx,o) in enumerate(outlets(cp))
        e_der[1,length(inlets(cp))+idx,1]=-o.h.val;e_deriv[1,length(inlets(cp))+idx,2]=-o.m.val
    end
    der+=e_der
    e_der=zeros(1,length(inlet(cp)),3)
    if haskey(cp.conns,:drain)
        e_der[1,1,3]=dhdp(hi.p.val,0);e_der[1,4,2]=-1
    else
        e_der[1,1,3]=dhdp(hi,p.val,0);e_der[1,3,2]=-1
    end
    der+=e_der

    p_der=zeros(2,length(inlets(cp)),3)
    if haskey(cp.conns,:drain)
        p_der[1,1,3]=1;p_der[1,4,3]=-1;p_der[1,3,3]=1;p_der[1,5,3]=-1
    else
        p_der[1,1,3]=1;p_der[1,3,3]=-1;p_der[1,2,3]=1;p_der[1,4,3]=-1
    end
    deriv+=p_der

    return der

end

function mass_res(cp::Condenor)
    si=cp.conns[:steam];wo=cp.conns[:waterout];ci=conns[:coolin];co=conns[:coolout]
    res=[]
    if haskey(cp.conns,:drain)
        di=conns[:drain]
        res+=si.m.val+di.m.val-wo.m.val
    else
        res+=si.m.val-wo.m.val
    end
    res+=ci.m.val-co.m.val
    return res
end

function mass_deriv(cp::Condensor)
    mat_deriv=nothing
    if haskey(cp.conns,:drain)
        mat=zeros(2,5,3)
        ma[1,1,1]=1;mat[1,2,1]=1;mat[1,3,1]=-1
        mat[2,4,1]=1;mat[2,5,1]=-1
    else
        mat=zeros(2,4,3)
        ma[1,1,1]=1;mat[1,2,1]=-1
        mat[2,3,1]=1;mat[2,4,1]=-1
    end
    return mat_deriv
end


export Condensor,equations,jacobi,setattrs,addconnection