using ForwardDiff
using Parameters



abstract type AbstractComponent end

mutable struct Var{T<:Real}
    val::T
    val0::T
    isset::Bool
end
Var()=Var(0.0,0.0,false)
Var(set)=Var(set,set,true)

mutable struct Connection 
    source::AbstractComponent
    sid::Symbol
    target::AbstractComponent
    tid::Symbol
    m::Var
    h::Var
    p::Var
    function Connection(source::AbstractComponent,sid::Symbol,target::AbstractComponent,tid::Symbol)
        new(source,sid,target,tid,Var(),Var(),Var())
    end
end

function setattrs(c::Connection;kwargs...)
    opts=[:m,:h,:p]
    for kw in keys(kwargs)
        if kw==:m
            print("m")
            c.m=Var(kwargs[kw])
        end
        if kw==:h
            c.h=Var(kwargs[kw])
        end
        if kw==:p
            c.p=Var(kwargs[kw])
        end
    end
end

mutable struct Bus
    label::String
    P::Var
    comps::Vector{AbstractComponent}
end

mutable struct Network
    conns::Vector{Connection}
    comps::Vector{AbstractComponent}
    buses::Vector{Bus}
    function Network()
        conns=Vector{Connection}()
        comps=Vector{AbstractComponent}()
        buses=Vector{Bus}()
        new(conns,comps,buses)
    end
end

function Bus(label::String)
    Bus(label,Var(undef,undef,false,false),Vector{Component}())
end

addcomponent(bus::Bus,cp::AbstractComponent)=push!(bus.comps,cp)

function addcomponent(bus::Bus,comps...)
    for comp in comps
        addcomponent(bus,comp)
    end
end


function connect(s::AbstractComponent,sid::Symbol,t::AbstractComponent,tid::Symbol)
    c=Connection(s,sid,t,tid)
    c.m.val=float(1);c.h.val=float(500);c.p.val=float(1)
    addconnection(s,sid,c)
    addconnection(t,tid,c)
    return c
end

equations(comp::AbstractComponent)=[]

jacobi(comp::AbstractComponent,c::Connection)=[]

initcomponent(comp::AbstractComponent)=nothing

initsource(comp::AbstractComponent,c::Connection)=zeros(3)

inittarget(comp::AbstractComponent,c::Connection)=zeros(3)

function derive(cp::AbstractComponent,func::Function,pos::Int64,dx::Symbol)

    dm,dh,dp=0.,0.,0.
    if dx == :m
        dm=1.e-4
    elseif dx==:h
        dh=1.e-3
    elseif dx==:p
        dp=1.e-5
    else
        dx=1.e-5
    end

    conns=vcat(inlets(cp),outlets(cp))
    exp=0.0
    conns[pos].m.val+=dm
    conns[pos].h.val+=dh
    conns[pos].p.val+=dp
    
    exp+=func(cp)
    
    conns[pos].m.val-=2.0*dm
    conns[pos].h.val-=2.0*dh
    conns[pos].p.val-=-2.0*dp
    exp-=func(cp)
    
    derive=exp/(dm+dh+dp)/2.0
    
    return derive
end

