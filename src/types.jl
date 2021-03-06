using ForwardDiff
using Parameters



abstract type AbstractComponent end



abstract type AbstractVariable end

mutable struct Ref
    conn::Connection
    factor::Float64
    delta::Float64
end

Ref(c,f,d)=Ref(c,convert(Float64,f),convert(Float64,d))

@with_kw mutable struct VarProp <:AbstractVariable
    val::Float64=0
    val_set::Bool=false
    ref::Union{Nothing,Ref}=nothing
    ref_set::Bool=false
end


@with_kw mutable struct VarComp <: AbstractVariable
    val::Float64=0
    val_set::Bool=false
    isvar::Bool=false
end

@with_kw mutable struct VarChar <: AbstractVariable
    func::Function=nothing
    x::Vector{Float64}=Vector{Float64}()
    y::Vector{Float64}=Vector{Float64}()
    isset::Bool=false
end

mutable struct Connection 
    source::AbstractComponent
    sid::Symbol
    target::AbstractComponent
    tid::Symbol
    m::VarProp
    h::VarProp
    p::VarProp
    T::VarProp
    
    function Connection(source::AbstractComponent,sid::Symbol,target::AbstractComponent,tid::Symbol)
        
        !isa(source,AbstractComponent) && error("$source is not a type of AbstractComponent")
        !isa(target,AbstractComponent) && error("$target is not a type of AbstractComponent")
        !in(sid,portnams(source)) && error("$sid is not a port of $source")
        !in(tid,portnames(target)) && error("$tid is not a port of $target")

        new(source,sid,target,tid,VarProp(),VarProp(),VarProp(),VarProp())
    end
end

function setattrs(comp::AbstractComponent;kwarg...)
    attrs=comp.attrs
    for (key,value) in kwargs
        if key in opts
            attrs[key]=VarComp(convert(Float64,value))
        end
    end
end

function setattrs(c::Connection;kwargs...)
    for (key,value) in kwargs
        if key==:m
            if isa(value,Ref)
                c.m.ref=value
                c.m.ref_set=true
            else
                c.m.val=convert(Float64,value)
                c.m.val_set=true
            end
        end
        if key==:h
            if isa(value,Ref)
                c.h.ref=value
                c.h.ref_set=true
            else
                c.h.val=convert(Float64,value)
                c.h.val_set=true
            end
        end
        if key==:p
            if isa(value,Ref)
                c.p.ref=value
                c.p.ref_set=true
            else
                c.p.val=convert(Float64,value)
                c.p.val_set=true
            end
        end
        if key==:T
            if isa(value,Ref)
                c.T.ref=value
                c.T.ref_set=true
            else
                c.T.value=convert(Float64,value)
                c.T.val_set=true
            end
        end
    end
end

to_flow(c::Connection)=[c.m.val,c.h.val,c.p.val]

Base.show(c::Connection)=println("source=$(c.source.label) target=$(c.target.label) m=$(c.m.val) h=$(c.h.val) p=$(c.p.val)")

mutable struct Bus
    label::String
    P::VarComp
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

function mass_res(cp::AbstractComponent)
    if isa(cp,Source) || isa(cp,Sink)
        return []
    end

    if isa(cp,Split) || isa(cp,Merge) || isa(cp,Dearator) || isa(cp,Turbine) 
        || isa(cp,Pipe) || isa(cp,Valve)
        res=0
        for i in inlets(cp)
            res+=i.m.val
        end
        for o in outlets(cp)
            res+=o.m.val
        end
        return res
    end
end

function mass_deriv(comp::AbstractComponent)
    if isa(cp,Split) || isa(cp,Merge) || isa(cp,Dearator) || isa(cp,Turbine) 
        || isa(cp,Pipe) || isa(cp,Valve)
        mat_deriv=zeros(1,length(inlets(cp))+length(outles(cp)),3)
        j=1
        for i in inlets(cp)
            mat_deriv[1,j,1]=1.0
            j+=1
        end
        
        for o in outlets(cp)
            mat_deriv[1,j,1]=-1
            j+=1
        end
        return mat_deriv
    end
end

function addconnection(cp::AbstractComponent,port::Symbol,c::Connection)
    ports=ports(cp)
    if port in ports
        comp.conns[port]=c
    else
        println("wrong port name")
    end

end

equations(comp::AbstractComponent)=[]

jacobi(comp::AbstractComponent,c::Connection)=[]

initcomponent(comp::AbstractComponent)=nothing

initsource(comp::AbstractComponent,c::Connection)=zeros(3)

inittarget(comp::AbstractComponent,c::Connection)=zeros(3)

function derive(cp::AbstractComponent,func::Function,pos::Int64,dx::Symbol)

    dm,dh,dp=0.,0.,0.
    if dx == :m
        dm=cbrt(eps())
    elseif dx==:h
        dh=cbrt(eps())
    elseif dx==:p
        dp=cbrt(eps())
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
    conns[pos].p.val-=2.0*dp
    exp-=func(cp)
    
    derive=exp/(dm+dh+dp)/2.0
    
    return derive
end

function Base.getproperty(a::AbstractComponent,x::Symbol)
    props=propertynames(a)
    if x in props
        return getfield(a,x)
    elseif x in ports(a)
        return getfield(a,conns)[x]
    elseif x in attrs(cp)
        return getfield(a,attrs)[x]
    else
        error("no such filed")
    end
end

function Base.setproperty!(a::AbstractComponent,name::Symbol,x)
    props=propertynames(a)
    if name in props
        setfie        
    elseif name in attrs
        getfield(a,:attrs)[name]=x
    elseif name in ports
        getfield(a,:conns)[name]=x
    end
end

len(cp::AbstractComponent)=length(inlets(cp))+length(outlets(cp))

function zeta_func(cp::AbstractComponent)
        r"""
        calculates pressure drop from zeta (zeta1 for heat exchangers)

        :returns: residual value for the pressure drop

        .. math::

            \zeta = \frac{\Delta p \cdot v \cdot 2}{c^2}\\
            c = \frac{\dot{m} \cdot v}{A}

        As the cross sectional area A will not change from design to offdesign
        calculation, it is possible to handle this the following way:

        .. math::
            0 = \zeta - \frac{(p_{in} - p_{out}) \cdot \pi^2}{8 \cdot
            \dot{m}_{in}^2 \cdot \frac{v_{in} + v_{out}}{2}}
        """
        i = inlets(cp)[1]
        o = outlets(cp)[1]
        if haskey(cp.attrs,:zeta):
            val = self.zeta.val
        else:
            val = self.zeta1.val
        end
        return (val - (i.p.val - o.p.val) * pi^2 /
                (8 * i.m.val^2 * (phv(i.p.val,i.h.val) + phv(o.p.val,o.h.val)) / 2))
end

function zeta2_func(cp::AbstractComponent)
    r"""
    calculates pressure drop from zeta2
    :returns: residual value for the pressure drop
    .. math::
        \zeta_2 = \frac{\Delta p_2 \cdot v_2 \cdot 2}{c_2^2}\\
        c_2 = \frac{\dot{m}_2 \cdot v_2}{A_2}
    As the cross sectional area A will not change from design to offdesign
    calculation, it is possible to handle this the following way:
    .. math::
        0 = \zeta_2 - \frac{(p_{2,in} - p_{2,out}) \cdot \pi^2}{8 \cdot
        \dot{m}_{2,in}^2 \cdot \frac{v_{2,in} + v_{2,out}}{2}}
    """
    i = inlets(cp)[2]
    o = outlets(cp)[2]
    return (cp.zeta2.val - (i.p.val - o.p.val) * pi^2 /
            (8 * (i.m.val)^2 * (phv(i.p.val,i.h.val) + phv(o.p.val,o.h.val)) / 2))
end
