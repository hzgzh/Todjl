@with_kw mutable struct Turbine<:AbstractComponent
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,Float64}
    mode=:design
    method=:nomethod
    Turbine(s::String)=new(s,Dict{Symbol,Connection}(),Dict{Symbol,Float64}(),:design)
end

function (cp::Turbine)(;kwargs)
    
end

function setattrs(comp::Turbine;kwargs...)
    #cq-流量系数 pr-压比 cpr-临界压比 eta-设计效率
    opts=[:cq,:pr,:cpr,:eta];attrs=comp.attrs
    for key in keys(kwargs)
        if key in opts
            attrs[key]=kwargs[key]
        end
    end
end

function callparamers(comp::Turbine)
    in=comp.conns[:in];out=comp.conns[:out]
    ratio=out.p.val/in.p.val
    criticalratio=0.7ratio
    cq=in.m.val/sqrt(in.p.val,ptv(in.p.val,in.h.val))/sqrt((1-(ratio-criticalration)^2/(1-criticalratio)^2))
    eff=(in.h.val-out.h.val)/(in.h.val-psh(in.p.val,phs(out.p.val,in.h.val)))
    (ratio,criticalratio,cq,eta)
end

inlets(comp::Turbine)=[comp.conns[:in]];outlets(comp::Turbine)=[comp.conns[:out]]

function addconnection(comp::Turbine,port::Symbol,c::Connection)
    ports=[:in,:out]
    if port in ports
        comp.conns[port]=c
    else
        print("wrong name")
    end
end
function equations(comp::Turbine)
    in=comp.conns[:in];out=comp.conns[:out];attrs=comp.attrs
    mode=comp.mode
   
    res=zeros(0)
    if mode == :design
        pr=attrs[:pr];eta=attrs[:eta]
        push!(res,in.m.val-out.m.val)
        push!(res,etasfunc(comp))
        push!(res,pr*in.p.val-out.p.val)
        attrs[:cpr]=pr*0.7
        attrs[:cq]=in.m.val/sqrt(in.p.val/phv(in.p.val,in.h.val))*sqrt(1-((out.p.val/in.p.val-0.7pr)/(1.0-0.7pr))^2)
    end
    if mode == :offdesign
        cq=attrs[:cq];;cpr=attrs[:cpr]
        push!(res,in.m.val-out.m.val)
        push!(res,conefun(comp))
        push!(res,etasfunc(comp))
    end
    return res
end

function jacobi(comp::Turbine,c::Connection)
    in=comp.conns[:in];out=comp.conns[:out];attrs=comp.attrs
    mode=comp.mode;pr=attrs[:pr]
    
    jac=zeros(3,3)
    if in==c
        if mode==:design
            gradeta=etasderive(comp,c)
            jac[1,1]=1;jac[2,2]=gradeta[1];jac[2,3]=gradeta[2]
            jac[3,1]=pr;
        end
        if mode==:offdesign
            gradcone=conederive(comp,c);gradetas=etasderive(comp,c)
            jac[1,1]=1
            jac[2,1]=gradcone[1];jac[2,2]=gradcone[2];jac[2,3]=gradcone[3]
            jac[3,2]=gradetas[1];jac[3,3]=gradetas[2]
        end
    end

    if out==c
        if mode==:design
            gradeta=etasderive(comp,c)
            jac[1,1]=-1.0;jac[2,2]=gradeta[1];jac[2,3]=gradeta[2]
            jac[3,3]=-1.0
        end
        if mode==:offdesign
            gradcone=conederive(comp,c);gradetas=etasderive(comp,c)
            jac[1,1]=-1.0;
            jac[2,1]=gradcone[1];jac[2,2]=gradcone[2];jac[2,3]=gradcone[3]
            jac[3,2]=gradetas[1];jac[3,3]=gradetas[2]
        end
    end
    return jac
end

function conefunc(cp::Turbine)
    in=cp.conns[:in];out=cp.conns[:out];cq=cp.attrs[:cq];cpr=cp.attrs[:cpr]
    in.m.val-cq*sqrt(in.p.val/phv(in.p.val,in.h.val))*sqrt(1.0-((out.p.val/in.p.val-cpr)/(1.0-cpr))^2.0)
end

function conederive(cp::Turbine)
    dm,dh,dp=0.,0.,0.
    if c==cp.conns[:in]
        dm=derive(cp,conefunc,1,:m)
        dh=derive(cp,conefunc,1,:h)
        dp=derive(cp,conefunc,1,:p)
    end
    
    if c==cp.conns[:out]
        dm=derive(cp,conefunc,1,:m)
        dh=derive(cp,conefunc,1,:h)
        dp=derive(cp,conefunc,1,:p)
    end
    (dm,dh,dp)
end

function etasfunc(cp::Turbine)
    in=cp.conns[:in];out=cp.conns[:out];eta=cp.attrs[:eta]
    in.h.val-out.h.val-eta*(in.h.val-psh(out.p.val,phs(in.p.val,in.h.val)))
end

function etasderive(cp::Turbine,c::Connection)
    dh,dp=0.,0.
    if c==cp.conns[:in]
        dh=derive(cp,etasfunc,1,:h)
        dp=derive(cp,etasfunc,1,:p)    
    end

    if c==cp.conns[:out]
        dh=derive(cp,etasfunc,2,:h)
        dp=derive(cp,etasfunc,2,:p)
    end
    (dh,dp)
end

function busfunc(comp::Turbine)
    in=comp.conns[:in];out=comp.conns[:out];attrs=comp.attrs
    cq=attrs[:cq];eta=attrs[:eta];cpr=attrs[:cpr]

    res=[]
    push!(in.m.val*(in.h.val-out.h.val))
    res
end

function busjacaobi(comp::Turbine,c::Connection)
    in=comp.conns[:in];out=comp.conns[:out];attrs=comp.attrs
    cq=attrs[:cq];eta=attrs[:eta];cpr=attrs[:cpr]
    jac=zeros(1,3)
    if c==in
        jac[1,1]=in.h_val-out.h.val;jac[1,2]=in.m.val
    end
    if c==out
        jac[1,2]=-in.m.val
    end
end

function checkconverge(comp::Turbine)
    in,out=comp[:in],comp[:out]
    if !in.p.isset && in.p.val<1
        i.p.val=1
    end
    if in.p.val<=out.p.val && !out.p.isset
        out.p.val=in.p.val/2
    end
    if in.h.val<100 && !in.h.isset
        in.h.val=100
    end
    if in.h.val<out.h.val && out.h.isset
        out.h.val=in.h.val*0.75
    end
end

initsource(comp::Turbine,c::Connection)=[1.0,3000.,0.5]

inittarget(comp::Turbine,c::Connection)=[1.0,3200.,25.]


export Turbine,equations,jacobi,setattrs,addconnection