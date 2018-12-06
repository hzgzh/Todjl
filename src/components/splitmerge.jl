mutable struct Merge
    conn
    function Merge()
        c=Dict{String,Connection}()
        new(c)
    end
end
function addconnection(comp::Split,s::String,c::Connection)
    port=["in1","in2","out1"]
    if s in port
        comp.c[s]=c
    else
        println("wrong input")
    end
end

function jacobi(comp::Merge,c::Connection)
    in1=comp.conn["in1"];in2=comp.conn["in2"];o=comp.conn["out1"]
    if in1==c
        jac[1,1]=1;jac[2,1]=in1.h.val;jac[2,2]=in1.m.val
        jac[3,3]=1.0
    end
    if in2==c
        jac[1,1]=1.0;jac[2,1]=in2.h.val;jac[2,2]=in2.m.val
        jac[4,3]=1.0
    end
    if o==c
        jac[1,1]=-1.0;jac[2,1]=-o.h.val;jac[2,2]=-o.m.val
        jac[3,3]-1.0;jac[4,3]=-1.0
    end
end
function equations(comp::Merge)
    in1=comp.conn["in1"];in2=comp.conn["in2"];o=comp.conn["out1"]
    res=[]
    push!(res,in1.m.val+in2.m.val-o.m.val)
    push!(res,in1.m.val*in1.h.val+in2.m.val*in2.h.val-o.m.val*o.h.val)
    push!(res,in1.p.val-o.p.val)
    push!(res,in2.p.val-o.p.val)
end

function numofequations(comp::Split,mode)
    4
end

mutable struct Split
    conn
    function Split()
        c=Dict{String,Connection}()
        new(c)
    end
end
function addconnection(comp::Split,s::String,c::Connection)
    port=["in","out1","out2"]
    if s in port
        comp.c[s]=c
    else
        println("wrong input")
    end
end
function jacobi(comp::Split,c::Connection)
    in=comp.conn["in"];o1=comp.conn["out1"];o2=comp.conn["out2"]
    jac=zeros(4,3)
    if c==in
        jac[1,1]=1.0;jac[2,1]=in.h.val;jac[2,2]=in.m.val
        jac[3,3]=1.0;jac[4,3]=1.0
    end
    if c==o1
        jac[1,1]=-1.0;jac[2,1]=-o1.h.val;jac[2,2]=-o1.m.val
        jac[3,3]=-1.0
    end
    if c==o2
        jac[1,1]=-1.0;jac[2,1]=-o2.h.val;jac[2,2]=-o2.m.val
        jac[4,3]=-1.0
    end
end
function equations(comp::Split)
    in=comp.conn["in"];o1=comp.conn["out1"];o2=comp.conn["out2"]
    res=[]
    push!(res,in.m.val-o1.m.val-o2.m.val)
    push!(res,in.m*in.h.val-o1.m.val*o1.h.val-o2.m.val*o2.h.val)
    push!(res,in.p.val-o1.p.val)
    push!(res,in.p.val-o2.p.val)
end
numberofequations(comp::Split)=4

