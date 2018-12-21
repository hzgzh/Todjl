@with_kw mutable struct SimpleHeater <: AbstractComponent
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,CVar}
    
    SimpleHeater(s::String)=new(s,Dict{Symbol,Connection}(),Dict{Symbol,Float64}())
end


inlets(cp::SimpleHeater)=[cp.conns[:in]]
outlets(cp::SimpleHeater)=[cp.conns[:out]]
portnames(cp::SimpleHeater)=[:in,:out]
function equations(cp::SimpleHeater)
    vec=[]
    vec+=mass_res(cp)
    if cp.attrs[:pr].isset
        p_res=cp.attrs[:pr].val*inlets(cp)[0].p.val-outlets(cp).p.val
        vec+=p_res
    end
    return vec
end

function derivatives(cp::SimpleHeater)
    i=cp.conns[:in];o=cp.conns[:out];attrs=cp.attrs
    der=mass_deriv(cp)
    if cp.attrs[:pr].isset
        p_der=zeros(1,2,3)
        p_der[1,1,3]=cp.attrs[:pr];p_der[1,2,3]=-1
        der+=p_der
    end
    return der
end

export SimpleHeater,equations,jacobi,setattrs,addconnection