@with_kw mutable struct Pump<:AbstractTurboMachine
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,Float64}
    Pump(s::String)=new(s,Dict{Symbol,Connection}(),Dict{Symbol,Float64}())
end

inlets(cp::Pump)=[cp.conns[:in]]
outlets(cp::Pump)=[cp.conns[:out]]
ports(cp::Pump)=[:in,:out]
attrs(cp::Pump)=[:pr,:eta]

function additional_equations(cp::Pump)
    in=cp.in;out=cp.out;
    vec=[]
        
    if cp.eta.val_set
        eta_res=eta_func(cp)
        vec+=eta_res
    end

    return vec
end

function additional_derivatives(cp::Pump,c::Connection)
    in=cp.in;out=cp.out
    if cp.eta.val_set
        e_der=eta_deriv(cp)
        der+=e_der
    end

    return der
end

function eta_func(cp::Pump)
    i=cp.in;o=cp.out
    return o.h.val-i.h.val-cp.eta.val*(o.h.val-psh(o.p.val,phs(i.p.val,i.h.val)))
end

function eta_der(cp::Pump)
    der=zeros(1,2,3)
    der[1,1,2]=derive(cp,eta_func,1,:h)
    der[1,1,3]=derive(cp,eta_func,1,:p)
    der[1,2,2]=derive(cp,eta_func,2,:h)
    der[1,2,3]=derive(cp,eta_func,2,:p)
    return der
end

function checkconverge(cp::Pump)
    in,out=cp.in,cp.out
    if !out.p.val_set && out.p.val<in.p.val
        out.p.val=2.0in.p.val
    end
    if !in.p.val_set && out.p.val<in.p.val
        in.p.val=0.5out.p.val
    end
    if !out.h.val_set && out.h.val<in.h.val
        out.h.val=1.1*in.h.val
    end
    if !in.h.val_set && out.h.val<in.h.val
        in.h.val=0.9*out.h.val
    end
end

function bus_func(cp::Pump)
    i=cp.in;o=cp.out
    return i.m.val*(o.h.val-i.h.val)
end

function bus_deriv(cp::Pump)
    i=cp.in;o=cp.out
    b_der=zeros(1,2,3)
    b_der[1,1,1]=o.h.val-i.h.val
    b_der[1,1,2]=-i.m.val
    b_der[1,2,2]=i.m.val
    return b_der
end

initsource(cp::Pump,c::Connection)=[1.0,300.,10.]

inittarget(cp::Pump,c::Connection)=[1.0,290.,1.]

export Pump,equations,jacobi,setattrs,addconnection