@with_kw mutable struct Pipe<:AbstractComponent
    label::String
    conns::Dict{Symbol,Connection}
    attrs::Dict{Symbol,CVar}
    mode=:design
    method=:no
    Pipe(s::String)=new(s,Dict{Symbol,Connection}(),Dict{Symbol,Float64}())
end


inlets(cp::Pipe)=[cp.conns[:in]]
outlets(cp::Pipe)=[cp.conns[:out]]
ports(cp::Pipe)=[:in,:out]
attrs(cp::Pipe)=[:pr]
function equations(cp::Pipe)
    i=cp.conn[:in];o=cp.conn[:out]
    vec=[]
    vec+=mass_res(cp)
    vec+=i.h.val-o.h.val
    vec+=cp.pr.val*i.p.val-o.p.val
    vec
end

function derivatives(cp::Pipe)
    i,o=cp.in,cp.out
    der=mass_deriv(cp)
    e_der=zeros(1,2,3)
    e_der[1,1,2]=1;e_der[1,2,2]=-1
    der+=e_der
    p_der=zeros(1,2,3)
    p_der[1,1,3]=cp.pr.val;p_der[1,2,3]=-1
    der+=p_der
    return der
end

export Pipe,equations,jacobi,setattrs,addconnection