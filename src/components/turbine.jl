export Turbine,setattr,addconnection,euqations,jacobi

mutable struct Turbine<:Component
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,Float64}
    
    flow_ratio::Float64
    design_ratio::Float64
    critical_ratio::Float64
    design_eta::Float64

end

function calldesignparamers(comp::Turbine)
    in=comp.conns["in"];out=comp.conns["out"]
    ratio=out.p.val/in.p.val
    criticalratio=0.7ratio
    cq=in.m.val/sqrt(in.p.val,ptv(in.p.val,in.h.val))/sqrt((1-(ratio-criticalration)^2/(1-criticalratio)^2))
    eff=(in.h.val-out.h.val)/(in.h.val-psh(in.p.val,phs(out.p.val,in.h.val)))
    (ratio,criticalratio,cq,eta)
end

function addconnection(comp::Turbine,s::String,c:Connection)
    ports=[:in,:out]
    if s in ports
        comp.conns[s]=c
    else
        print("wrong name")
    end
end
function equations(comp::Turbine)
    in=comp.conns[:in];out=comp.conns[:out]
    flow_ratio=comp.flow_ratio;design_eta=comp.design_eta;
    critical_ratio=comp.critical_ratio;

    res=[]
    push!(res,in.m.val-out.m.val)
    push!(res,in.m.val-flow_ratio*sqrt(in.p.val/phv(in.p.val,in.h.val))*sqrt(1-((out.p.val/in.p.val-critical_ratio)/(1-critical_ratio))^2))
    push!(res,in.h.val-out.h.val-design_eta*(in.h.val-psh(out.p.val,phs(in.p.val,in.h.val))))
end



function jacobi(comp::Turbine,c::Connection)
    in=comp.conns["in"];out=comp.conns["out"]
    flow_ratio=comp.flow_ratio;design_eta=comp.design_eta;
    critical_ratio=comp.critical_ratio

    f1(p0,h0,p1)=-flow_ratio*sqrt(p0/phv(p0,h0))*sqrt(1-((p1/p1-critical_ratio)/(1-critical_ratio))^2)
    f2(p0,h0,p1,h1)=h0-h1-design_ratio*(h0-psh(p1,phs(p0,h0)))

    grad1=ForwardDiff.grad(f1,[in.p.val,in.h.val,out.p.val])
    grad2=ForwardDiff.grad(f2,[in.p.val,in.h.val,out.p.val,out.h.val])
    jac=zeros(6,3)
    if in==c
        jac[1,1]=1
        jac[2,1]=1;jac[2,2]=grad1[2];jac[2,3]=grad1[1]
        jac[3,2]=grad2[2];jac[3,3]=grad2[1]
    end

    if out==c
        jac[1,1]=-1.0;
        jac[2,3]=grad1[3]
        jac[3,2]=grad2[4];jac[3,3]=grad2[3]
    end
end