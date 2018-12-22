abstract type AbstractHeater <: AbstractComponent end

inlets(cp::AbstractHeater)=
    haskey(cp.conns,:drain) ? [cp..hi,cp.ci,cp..drain] : [cp.hi,cp.ci]

outlets(cp::AbstractHeater)=[cp.ho cp.co]
ports(cp::AbstractHeater)=[:hi,:ci,:drain,:ho,:co]
attrs(cp::AbstractHeater)=[:tu,:tl,:pr1,:pr2,:ka1,:ka2,:zeta1,:zeta2]
design_mode(cp::AbstractHeater)=[:tu,:tl,:pr1,:pr2]
offdesign_mode(cp::AbstractHeater)=[]

function equations(cp::AbstractHeater)
    ci,cp,hi,ho=cp.ci.ci,cp.co,cp.hi,cp.ho
    di=haskey(cp.conns,:drain) ? cp.drain : nothing;
   
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

    cp.tu.val_set && vec+=ttd_u_func(cp)
    cp.tl.val_set &&  vec+=ttd_l_func(cp)
    cp.pr1.val_set && vec+=hi.p.val-ho.p.val
    cp.pr2.val_set && vec+=ci.p.val-co.p.val
    
    res+=aux_equations(cp)
    return res
end

function derivatives(cp::AbstractHeater)
    ci,cp,hi,ho=cp.ci.ci,cp.co,cp.hi,cp.ho
    di=haskey(cp.conns,:drain) ? cp.drain : nothing

    der=mass_deriv(cp)
        
    e_derive=zeros(1,length(inlets(cp))+length(outlets(cp)),3)
    for (idx,i) in enumerate(inlets(cp))
        e_der[1,idx,1]=i.h.val;e_derive[1,idx,2]=i.m.val
    end

    for (idx,o) in enumerate(outlets(cp))
        e_der[1,length(inlets(cp))+idx,1]=-o.h.val;e_der[1,length(inlets(cp))+idx,2]=-o.m.val
    end
    der+=e_der

    cp.td.val_set && der+=ttd_u_deriv(cp)
    cp.tl.val_set && der+=ttd_l_deriv(cp)
    

    
    if cp.pr1.val_set
        p_der=zeros(1,length(inlets(cp))+length(outlets(cp)),3)
        if haskey(cp.conns,:drain)
            p_der[1,1,3]=1;p_dere[1,4,3]=-1
        else
            p_der[1,1,3]=1;p_der[1,3,3]=-1
        end
    end
    der+=p_deriv

    if cp.pr2.val_set
        p_der=zeros(1,length(inlets(cp))+length(outlets(cp)),3)
        if haskey(cp.conns,:drain)
            p_der[1,3,3]=1;p_dere[1,5,3]=-1
        else
            p_der[1,2,3]=1;p_der[1,4,3]=-1
        end
    end
    der+=p_deriv
    der+=aux_derivatives(cp)
    return der
end

function ttd_u_func(cp::AbstractHeater)
    co,hi=cp.co,cp.hi
    return co.h.val-pth(co.p.val,pt(hi.p.val)+cp.tu.val)
end

function ttd_u_deriv(cp::AbstractHeater)
    deriv=zeros(1,length(inlets(cp))+length(outlets(cp)),3)
    deriv[1,1,3]=derive(cp,ttd_u_func,1,:p)
    if haskey(cp.conns,:drain)
        deriv[1,5,2]=1;deriv[1,5,3]=derive(cp,ttd_u_func,5,:p)
    else
        deriv[1,4,2]=1;deriv[1,4,3]=derive(cp,ttd_u_func,4,:p)
    end
    return deriv
end

function ttd_l_func(cp::AbstractHeater)
    ho,ci=cp.ho,cp.ci
    return ho.h.val-pth(ho.p.val,pht(ci.p.val,ci.h.val)+cp.tl.val)
end

function ttd_l_deriv(cp::AbstractHeater)
    deriv=zeros(1,length(inlets(cp))+length(outlets(cp)),3)
    if haskey(cp.conns,:prevdrain)
        deriv[1,4,2]=1;deriv[1,4,3]=derive(cp,ttd_l_func,4,:p)
        deriv[1,3,2]=derive(cp,ttd_l_func,3,:h);deriv[1,3,3]=deriv(cp,ttd_l_func,3,:p)
    else
        deriv[1,3,2]=1;deriv[1,3,3]=derive(cp.ttd_l_func,3,:p)
        deriv[1,2,2]=derive(cp,ttd_l_func,2,:h);deriv[1,2,3]=derive(cp,ttd_l_func,2,:p)
    end
    return deriv
end

function mass_res(cp::AbstractHeater)
    ci,co,hi,ho=cp.ci,cp.co,cp.hi,cp.ho
    res=[]
    if haskey(cp.conns,:drain)
        di=conns.drain
        res+=hi.m.val+di.m.va-ho.m.val
        res+=ci.m.val-co.m.val
    else
        res+=hi.m.val-ho.m.val
        res+=ci.m.val-co.m.val
    end
    return res
end

function mass_deriv(cp::AbstractHeater)
    
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

export AbstractHeater,equations,jacobi,setattrs,addconnection