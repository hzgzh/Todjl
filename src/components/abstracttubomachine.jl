abstract type AbstractTurboMachine <: AbstractComponent

function callparamers(cp::AbstractTurboMachine)
    in=cp.conns[:in];out=cp.conns[:out]
    ratio=out.p.val/in.p.val
    criticalratio=0.7ratio
    cq=in.m.val/sqrt(in.p.val,ptv(in.p.val,in.h.val))/sqrt((1-(ratio-criticalration)^2/(1-criticalratio)^2))
    eff=(in.h.val-out.h.val)/(in.h.val-psh(in.p.val,phs(out.p.val,in.h.val)))
    (ratio,criticalratio,cq,eta)
end

inlets(cp::AbstractTurboMachine)=[cp.in]
outlets(cp::AbstractTurboMachine)=[cp.out]
ports(cp::AbstractTurboMachine)=[:in,:out]
attrs(cp::AbstractTurboMachine)=[:pr,:eta,:eta_char]

design_mode=[:pr,:eta]
offdesign_mode=[:eta_char]

function equations(cp::AbstractTurboMachine)
    in=cp.conns[:in];out=cp.conns[:out];attrs=cp.attrs
    res=[]
    res+=mass_res(cp)
    if cp.eta.val_set
       res+=eta_func(cp) 
    end

    if cp.pr.val_set
        res+=cp.conns[:pr]*in.p.val-out.p.val
    end
    rec+=additional_equations(cp)
    return res
end

function derivatives(cp::AbstractTurboMachine)
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
    der+=additional_derivatives(cp)
    return der
end

function eta_func(cp::AbstractTurboMachine)
    in=cp.in;out=cp.out;eta=cp.eta.val
    return in.h.val-out.h.val-eta*(in.h.val-psh(out.p.val,phs(in.p.val,in.h.val)))
end

function eta_der(cp::AbstractTurboMachine)
    der=zeros(1,2,3)
    der[1,1,2]=derive(cp,eta_func,1,:h);der[1,1,3]=derive(cp,eta_func,1,:p)
    der[1,2,2]=derive(cp,eta_func,2,:h);der[1,2,3]=derive(cp,eta_func,2,:p)
    return der
end

function bus_func(cp::AbstractTurboMachine)
    i=cp.in;o=cp.out
    return i.m.val*(i.h.val-o.h.val)
end

function bus_der(cp::AbstractTurboMachine)
    i=cp.in;o=cp.out
    der=zeros(1,2,2)
    der[1,1,1]=i.h.val-o.h.val;der[1,1,2]=i.m.val
    der[1,2,2]=-i.m.val
    return der
end

export AbstractTurboMachine,equations,jacobi,setattrs,addconnection