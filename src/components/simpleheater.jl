@with_kw mutable struct SimpleHeater <: AbstractComponent
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,CVar}
    
    SimpleHeater(s::String)=new(s,Dict{Symbol,Connection}(),Dict{Symbol,Float64}())
end


inlets(cp::SimpleHeater)=[cp.conns[:in]]
outlets(cp::SimpleHeater)=[cp.conns[:out]]
ports(cp::SimpleHeater)=[:in,:out]
attrs(cp::SimpleHeater)=[:pr]

function equations(cp::SimpleHeater)
    vec=[]
    vec+=mass_res(cp)
    if cp.pr.val_set
        p_res=cp.pr.val*inlets(cp)[0].p.val-outlets(cp)[0].p.val
        vec+=p_res
    end
    return vec
end

function derivatives(cp::SimpleHeater)
    i=cp.in;o=cp.out
    der=mass_deriv(cp)
    if cp.pr.val_set
        p_der=zeros(1,2,3)
        p_der[1,1,3]=cp.pr.val;p_der[1,2,3]=-1
        der+=p_der
    end
    return der
end

export SimpleHeater,equations,jacobi,setattrs,addconnection