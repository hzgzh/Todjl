@with_kw mutable struct Turbine<:AbstractComponent
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,Float64}
    mode=:design
    method=:nomethod
    Turbine(s::String)=new(s,Dict{Symbol,Connection}(),Dict{Symbol,Float64}(),:design)
end

function callparamers(cp::Turbine)
    in=cp.conns[:in];out=cp.conns[:out]
    ratio=out.p.val/in.p.val
    criticalratio=0.7ratio
    cq=in.m.val/sqrt(in.p.val,ptv(in.p.val,in.h.val))/sqrt((1-(ratio-criticalration)^2/(1-criticalratio)^2))
    eff=(in.h.val-out.h.val)/(in.h.val-psh(in.p.val,phs(out.p.val,in.h.val)))
    (ratio,criticalratio,cq,eta)
end

inlets(cp::Turbine)=[cp.conns[:in]];
outlets(cp::Turbine)=[cp.conns[:out]]
portnames(cp::Turbine)=[:in,:out]

function equations(cp::Turbine)
    in=cp.conns[:in];out=cp.conns[:out];attrs=cp.attrs
    res=[]
    res+=mass_res(cp)
    if cp.eta.val_set
       res+=eta_func(cp) 
    end

    if cp.pr.val_set
        res+=cp.conns[:pr]*in.p.val-out.p.val
    end

    if cp.cone.val_set
        res+=cone_func(cp)
    end

    return res
end

function derivatives(cp::Turbine)
    in=cp.conns[:in];out=cp.conns[:out];attrs=cp.attrs
    der=mass_der(cp)
    if cp.eta.val_set
        der+=eta_deriv(cp)
    end

    if cp.pr.val_set
        der+=cp.conns[:pr]*in.p.val-out.p.val
    end

    if cp.cone.val_set
        der+=cone_func(cp)
    end

    return der
end

function cone_func(cp::Turbine)
    in=cp.conns[:in];out=cp.conns[:out];cq=cp.attrs[:cq].val;cpr=cp.attrs[:cpr].val
    in.m.val-cq*sqrt(in.p.val/phv(in.p.val,in.h.val))*sqrt(1.0-((out.p.val/in.p.val-cpr)/(1.0-cpr))^2.0)
end

function cone_der(cp::Turbine)
    der=zeros(1,2,3)
    der[1,1,1]=1.0;der[1,1,2]=derive(cp,cone_func,1,:h);der[1,1,3]=derive(cp,cone_func,1,:p)
    der[1,2,2]=derive(cp,cone_func,2,:h);der[1,2,3]=derive(cp,cone_func,2,:p)
    return der
end

function eta_func(cp::Turbine)
    in=cp.conns[:in];out=cp.conns[:out];eta=cp.attrs[:eta].val
    return in.h.val-out.h.val-eta*(in.h.val-psh(out.p.val,phs(in.p.val,in.h.val)))
end

function eta_der(cp::Turbine)
    der=zeros(1,2,3)
    der[1,1,2]=derive(cp,eta_func,1,:h);der[1,1,3]=derive(cp,eta_func,1,:p)
    der[1,2,2]=derive(cp,eta_func,2,:h);der[1,2,3]=derive(cp,eta_func,2,:p)
    return der
end

function bus_func(cp::Turbine)
    i=cp.conns[:in];o=cp.conns[:out]
    return i.m.val*(i.h.val-o.h.val)
end

function bus_der(cp::Turbine)
    i=cp.conns[:in];o=cp.conns[:out]
    der=zeros(1,1,2)
    der[1,1,1]=i.h.val-o.h.val;der[1,1,2]=i.m.val
    der[1,2,2]=-i.m.val
    return der
end

function checkconverge(cp::Turbine)
    in,out=cp[:in],cp[:out]
    if !in.p.val_set && in.p.val<1
        i.p.val=1
    end
    if in.p.val<=out.p.val && !out.p.val_set
        out.p.val=in.p.val/2
    end
    if in.h.val<100 && !in.h.val_set
        in.h.val=100
    end
    if in.h.val<out.h.val && out.h.val_set
        out.h.val=in.h.val*0.75
    end
end

initsource(cp::Turbine,c::Connection)=[1.0,3000.,0.5]

inittarget(cp::Turbine,c::Connection)=[1.0,3200.,25.]


export Turbine,equations,jacobi,setattrs,addconnection