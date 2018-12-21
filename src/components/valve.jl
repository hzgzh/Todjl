@with_kw mutable struct Valve<:AbstractComponent
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,AbstractComponent}
    Valve(s::String)=new(s,Dict{Symbol,Connection}(),Dict{Symbol,AbstractComponent}())
end


inlets(cp::Valve)=[cp.conns[:in]]
outlets(cp::Valve)=[cp.conns[:out]]
portnames(cp::Valve)=[:out]


function equations(cp::Valve)
    i=cp.conn[:in];o=cp.conn[:out]
    res=[]
    res+=mass_res(cp)  
    res+=i.h.val-o.h.val
    if cp.conns[:pr].isset
        res+=i.p.val*cp.conns[:pr]-o.p.val
    end
    res
end

function derivativescp::Valve)
    i=cp.conn[:in];o=cp.conn[:out]
    der=mass_deriv(cp)
    e_der=zeros(1,2,3)
    e_der[1,1,2]=1;e_der[1,2,2]=-1
    der+=e_der

    if cp.conns[:pr].isset
        p_der=zeros(1,2,3)
        p_der[1,1,3]=cp.conns[:pr];p_der[1,2,3]=-1
        der+=p_der
    end
    der
end

export Valve,equations,jacobi,setattrs,addconnection