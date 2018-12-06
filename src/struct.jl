using ForwardDiff
export Component,Connection,Network,Var,Bus
export addcomponents,addconnections

abstract type Component end

mutable struct Var
    #变量值
    val::Float64
    #初始值
    val0::Float64
    #是否强制设置
    val_set::Float64
    is_set::Bool
end

mutable struct Bus

end


mutable struct Connection 
    s::Component
    sid::String
    t::Component
    tid::String
    m::Var
    h::Var
    p::Var
    function Connection(source,sid,target,tid)
        new(source,sid,target,tid)
    end
end

function equations(conn::Connection)
    m.val=ifelse(m.is_est,m.val_set,1.0)
    h.val=ifelse(h.is_set,h.val_set,500)
    p.val=ifelse(p.is_set,p.val_set,1.0)
    []
end

function jacobi(conn::Connection)
    []
end


mutable struct Network
    conns::Vector{Connection}
    comps::Vector{Component}
    buses::Vector{Bus}
    function Network()
        conns=Vector{Connection}()
        comps=Vector{Components}()
        buses=Vector{Bus}
    end
end

function addcomponents(nw::Network,comps...)
    for comp in comps
        push!(nw.comps,comp)
    end
end
function addconnections(nw::Network,conns...)
    for conn in conns
        push!(nw.conns,conn)
    end
end

function connect(s::Component,sid::String,t::Component,tid::String)
    c=Connection(s,t)
    addconnection(s,sid,c)
    addconnection(t,tid,c)
end

function assemble(nw::Network)
    conns=nw.conns;comps=nw.comps
    numequations=sum(numberofequations(comp) for comp in comps)
    numvars=3*length(conns)
    mat=zeros(numequations,numvars)
    row=0;col=0
    vec=[]
    for conn in conns
        jac=jacobi(c.source,conn)
        push!(vec,c.source)
        r,c=size(jac)
        mat[row:row+r,col:col+c]=jac
        row+=r;col+=c
        jac=jacobi(c.target,conn)
        push!(vec,c.target)
        r,c=size(jac)
        mat[row:row+r,col:col+c]=jac
        row=row+r;col=col+c
    end

    x=similar(vec)
    num=0
    for (index,conn) in enumerate(conns)
        conn.m.val=x[index]
        conn.h.val=x[index+1]
        conn.p.val=x[index+2]
    end
end

