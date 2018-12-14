

@with_kw mutable struct FeedwaterHeater<:AbstractComponent
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,Float64}
    FeedwaterHeater(s::String)=new(s,Dict{Symbol,Connection}(),Dict{Symbol,Float64}())
    mode=:design
    method=:td
end

function (cp::FeedwaterHeater)(;kwargs)
    
end

function setattrs(comp::FeedwaterHeater;kwargs...)
    opts=[:tu,:td];attrs=comp.attrs
    for key in keys(kwargs)
        if key in opts
            attrs[key]=kwargs[key]
        end
    end
end

inlets(comp::FeedwaterHeater)=[comp.conns[:waterin] comp.conns[:steam] comp.conns[:prevdrain]]
outlets(comp::FeedwaterHeater)=[comp.conns[:waterout],comp.conns[:drain]]

function addconnection(comp::FeedwaterHeater,s::Symbol,c::Connection)
    ports=[:waterin,:waterout,:steam,:drain,:prevdrain]
    if s in ports
        comp.conns[s]=c
    else
        println("wrong port name")
    end
end
function jacobi(cp::FeedwaterHeater,c::Connection)
    conns=cp.conns;attrs=cp.attrs
    wi=conns[:waterin];wo=conns[:waterout];si=conns[:steam];
    dwo=conns[:drain];dpo=haskey(conns,:prevdrain) ? conns[:prevdrain] : nothing;tu=attrs[:tu];td=attrs[:td]
    jac=zeros[7,3]
    if wi==c
        jac[1,1]=1.0;jac[3,1]=wi.h.val;jac[3,2]=wi.m.val;jac[4,3]=1.0;
        jac[7,2]=-dhdt(dwo.p,pht(wi.p.val,wi.h.val))*dtdh(wi.p.val,wi.h.val)
        jac[7,3]=-dhdt(dwo.p,pht(wi.p.val,wi.h.val)+fwh.dd)*dtdp(wi.p.val,wi.h.val)
    end
    if wo==c
        jac[1,1]=-1.0;jac[3,1]=-wo.h.val;jac[3,2]=-wo.m.val
        jac[4,3]=-1.0;
        jac[5,2]=1.0
        jac[5,3]=-dhdp(wo.p,pt(si.p.val)+fwh.td)
    end
    if si==c
        jac[2,1]=1.0;jac[3,1]=si.h.val;jac[3,2]=si.m.val
        jac[5,3]=1.0
        jac[6,3]=-dhdt(wo.p.val,pt(si.p.val)+fwh.td)*sat_dtdp(si.p.val)
    end
    if dwo=c
        jac[2,1]=-1.0;jac[3,1]=-dwo.h.val;jac[3,2]=-dwo.m.val
        jac[5,3]=-1.0
        jac[7,2]=1.0;jac[7,3]=-dhdp(dwo.p.val,pht(wi.p.val,wi.h.val)+fwh.dd)
    end
    if dpo==c
        jac[2,1]=1.0;jac[3,1]=dpo.h.val;jac[3,2]=dpo.m.val
    end
    return jac
end
function equations(cp::FeedwaterHeater,c::Connection)
    conns=cp.conns;attrs=cp.attrs
    wi=conns[:waterin];wo=conns[:waterout];si=conns[:steam];
    dwo=conns[:drain];dpo=haskey(conns,:prevdrain) ? conns[:prevdrain] : nothing;tu=attrs[:tu];td=attrs[:td]
    res=zeros(0)
    push!(res,wi.m.val-wo.m.val)
    if haskey(conns,:prevdrain)
        dpo=conns[:prevdrain]
        push!(res,si.m.val+dpo.m.val-dwo.m.val)
        push!(res,wi.m.val*wi.h.val-wo.m.val*wo.h.val+si.m.val*si.h.val-dwo.m.val*dwo.h.val+dpo.m.val*dpo.h.val)
    else
        push!(res,si.m.val-dwo.m.val)
        push!(res,wi.m.val*wi.h.val-wo.m.val*wo.h.val+si.m.val*si.h.val-dwo.m.val*dwo.h.val)
    end
    push!(res,wi.p.val-wo.p.val)
    push!(res,si.p.val-dwo.p.val)
    push!(res,wo.h.val-pth(wo.p.val,pt(si.p.val)+tu))
    push!(res,dwo.h.val-pth(dwo.p.val,pht(wi.p.val,wi.h.val)+td))
    return res
end
numberofequations(comp::FeedwaterHeater)=7

export FeedwaterHeater,equations,jacobi,setattrs,addconnection