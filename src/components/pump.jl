@with_kw mutable struct Pump<:AbstractComponent
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,Float64}
    mode=:design
    method=:no
    Pump(s::String)=new(s,Dict{Symbol,Connection}(),Dict{Symbol,Float64}())
end

function (cp::Pump)(;kwargs)
    
end

function setattrs(comp::Pump;kwargs...)
    
    opts=[:prise,:eta];attrs=comp.attrs
    for key in keys(kwargs)
        if key in opts
            attrs[key]=kwargs[key]
        end
    end
end

inlets(comp::Pump)=[comp.conns[:in]];outlets(comp::Pump)=[comp.conns[:out]]

function addconnection(comp::Pump,port::Symbol,c::Connection)
    ports=[:in,:out]
    if port in ports
        comp.conns[port]=c
    else
        print("wrong name")
    end
end

function equations(comp::Pump)
    in=comp.conns[:in];out=comp.conns[:out];attrs=comp.attrs
    prise=attrs[:prise];eta=attrs[:eta]
    res=zeros(0)
    push!(res,in.m.val-out.m.val)
    push!(res,out.h.val-in.h.val-9.81*prise/1000.)
    push!(res,out.p.val-in.p.val-9.8*prise/100)
    return res
end



function jacobi(comp::Pump,c::Connection)
    in=comp.conns[:in];out=comp.conns[:out];attrs=comp.attrs
    prise=attrs[:prise];eta=attrs[:eta]
    jac=zeros(3,3)
    if in==c
        jac[1,1]=1;jac[2,2]=-1.0;jac[3,3]=-1.0
    end

    if out==c
        jac[1,1]=-1.0;jac[2,2]=1.0;jac[3,3]=1.0
    end
    return jac
end

function checkconverge(comp::Pump)
    in,out=comp.conns[:in],comp.conns[:out]
    if !out.p.isset && out.p.val<in.p.val
        out.p.val=2.0in.p.val
    end
    if !in.p.isset && out.p.val<in.p.val
        in.p.val=0.5out.p.val
    end
    if !out.h.isset && out.h.val<in.h.val
        out.h.val=1.1*in.h.val
    end
    if !in.h.isset && out.h.val<in.h.val
        in.h.val=0.9*out.h.val
    end
end
function busfunc(comp::AbstractComponent)

end

initsource(comp::Pump,c::Connection)=[1.0,295.,1.0]
inittarget(comp::Pump,c::Connection)=[1.0,300.,10.]

export Pump,equations,jacobi,setattrs,addconnection