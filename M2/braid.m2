loadPackage("Graphs",Reload => true);

-- why is this not predefined
isNull = method(Dispatch=>Thing);
isNull(Nothing) := x -> true;
isNull(Thing) := x -> false;

-- return a new list with the i-th entry replaced by x
update = method();
update(List, ZZ, Thing) := (l, i, x) -> insert(i, x, drop(l, {i,i}));

-- braid constructor
Braid = new Type of HashTable;
braid = method();
braid(ZZ, List) := (n, word) -> (
    return new Braid from hashTable{(global n)=>n, (global word)=>word};
);

-- vertex object
Vertex = new Type of HashTable;
vertex = method();
vertex(ZZ,ZZ) := (row, column) -> (
    return new Vertex from hashTable{(global row) => row, (global column) => column};
);

toString(Vertex) := v -> concatenate("(", toString(v.row), ",", toString(v.column), ")");

Vertex == Vertex := (v1,v2) -> (
    return if v1.row == v2.row and v1.column == v2.column then true else false;
);

Vertex ? Vertex := (v1,v2) -> (
    if (v1.row < v2.row) or (v1.row == v2.row and v1.column < v2.column) then
    	return symbol <
    else if v1.row == v2.row and v1.column == v2.column then
    	return symbol ==
    else
    	return symbol >
);

-- edge object
Edge = new Type of HashTable;
edge = method();
-- Thing should really be some parent of RingElement
edge(Vertex, Vertex, Thing) := (s, t, var) -> (
    return new Edge from hashTable{
	(global s) => s,
	(global t) => t,
	(global var) => var
    };
);

Edge == Edge := (e1, e2) -> (e1.s == e2.s and e1.t == e2.t and e1.var == e2.var);

BraidRes = new Type of MutableHashTable;

emptyRes = method();
emptyRes(Ring) := (r) -> (
    return new BraidRes from hashTable{
	-- Vertex -> {Edge}
	-- Edges are ordered {upper-left, upper-right, lower-right, lower-left}
        (global adjacent) => new MutableHashTable,
	-- base ring
	(global r) => r,
	-- RingElement -> Edge
	(global edges) => new MutableHashTable
    };
);

copyRes = method();
copyRes(BraidRes) := (br) -> (
    return new BraidRes from hashTable{
	(global adjacent) => copy(br.adjacent),
	(global r) => br.r,
	(global edges) => copy(br.edges)
    };
);

addBrVertex = method();
addBrVertex(BraidRes, Vertex) := (br, v) -> (
    br.adjacent#v = {null, null, null, null};
);

addBrEdge = method();
addBrEdge(BraidRes, Edge, ZZ, ZZ) := (br, e, ps, pt) -> (
    br.adjacent#(e.s) = update(br.adjacent#(e.s), ps, e);
    br.adjacent#(e.t) = update(br.adjacent#(e.t), pt, e);
    br.edges#(e.var) = e;
);

removeBrVertex = method();
removeBrVertex(BraidRes, Vertex) := (br, v) -> (
    remove(br.adjacent, v);
);

removeBrEdge = method();
removeBrEdge(BraidRes, Edge) := (br, e) -> (
    ps := position(br.adjacent#(e.s), a -> (not isNull(a)) and (a.var == e.var));
    br.adjacent#(e.s) = update(br.adjacent#(e.s), ps, null);
    pt := position(br.adjacent#(e.t), a -> (not isNull(a)) and (a.var == e.var));
    br.adjacent#(e.t) = update(br.adjacent#(e.t), pt, null);
    remove(br.edges, e.var);
);

-- returns the underlying undirected graph of the given braid resolution
ugraph = method();
ugraph(BraidRes) := (br) -> (
    g := graph({});
    vertexIndices := new MutableHashTable;
        
    vertices := sort(keys(br.adjacent));
    for i from 0 to #vertices-1 do (
	v := vertices#i;
	g = addVertex(g, i);
	vertexIndices#v = i;
    );
    
    edgeVars = keys(br.edges);
    for i from 0 to #edgeVars-1 do (
	e := br.edges#(edgeVars#i);
	sIndex := vertexIndices#(e.s);
	tIndex := vertexIndices#(e.t);
	g = addEdge(g, set {sIndex, tIndex});
    );

    return (g,vertexIndices);
);

-- build the fully singular resolution of the given braid
-- Braid -> BraidRes
singularResolution = method();
singularResolution(Braid) := (b) -> (
    n := b.n;
    word := b.word;
    R := QQ[for i from 1 to 6*n-4+2*#word list "u"|i];

    br := emptyRes(R);
    
    -- add vertices
    for i from 1 to n do(
	addBrVertex(br, vertex(-3,i))
    );
    strands := new MutableList;
    for i from 0 to 2*n-1 do(
	strands#i = (vertex(-3,i//2+1), (i+1)%2+2);
    );
    nextEdgeIndex := 2*n;
    for i from 1 to n-1 do(
	v := vertex(-2,i);
	(v1, p1) := strands#(2*i-1);
	(v2, p2) := strands#(2*i);
        addBrVertex(br, v);
	e1 := edge(v, v1, R_nextEdgeIndex);
	addBrEdge(br, e1, 0, p1);
	nextEdgeIndex = nextEdgeIndex + 1;
	e2 := edge(v, v2, R_nextEdgeIndex);
	addBrEdge(br, e2, 1, p2);
	nextEdgeIndex = nextEdgeIndex + 1;
	strands#(2*i-1) = (v, 3);
	strands#(2*i) = (v, 2);
    );
    for i from 1 to n-1 do(
	v := vertex(-1,i);
	(v1, p1) := strands#(2*i);
	(v2, p2) := strands#(2*i+1);
	addBrVertex(br, v);
	e1 := edge(v, v1, R_nextEdgeIndex);
	addBrEdge(br, e1, 0, p1);
	nextEdgeIndex = nextEdgeIndex + 1;
	e2 := edge(v, v2, R_nextEdgeIndex);
	addBrEdge(br, e2, 1, p2);
	nextEdgeIndex = nextEdgeIndex + 1;
	strands#(2*i) = (v, 3);
	strands#(2*i+1) = (v, 2);
    );
    for i from 0 to #word-1 do(
	v := vertex(i,0);
	(v1, p1) := strands#(abs(word#i)-1);
	(v2, p2) := strands#(abs(word#i));
	addBrVertex(br, v);
	e1 := edge(v, v1, R_nextEdgeIndex);
	addBrEdge(br, e1, 0, p1);
	nextEdgeIndex = nextEdgeIndex + 1;
	e2 := edge(v, v2, R_nextEdgeIndex);
	addBrEdge(br, e2, 1, p2);
	nextEdgeIndex = nextEdgeIndex + 1;
	strands#(abs(word#i)-1) = (v, 3);
	strands#(abs(word#i)) = (v, 2);
    );
    for i from 0 to #strands-1 do(
	v := vertex(-3, i//2+1);
	(v1, p1) := strands#i;
	e1 := edge(v, v1, R_i);
	addBrEdge(br, e1, i%2, p1);
    );
    
    return br;
);

-- splits the crossing x in place
splitCrossing = method();
splitCrossing(BraidRes, ZZ) := (br, x) -> (
    v := vertex(x, 0);
    vl := vertex(x, -1);
    vr := vertex(x, 1);
    
    e0 := br.adjacent#v#0;
    p0 := position((br.adjacent)#(e0.t), a -> (not isNull(a)) and (a.var == e0.var));
    e1 := br.adjacent#v#1;
    p1 := position((br.adjacent)#(e1.t), a -> (not isNull(a)) and (a.var == e1.var));
    e2 := br.adjacent#v#2;
    p2 := position((br.adjacent)#(e2.s), a -> (not isNull(a)) and (a.var == e2.var));
    e3 := br.adjacent#v#3;
    p3 := position((br.adjacent)#(e3.s), a -> (not isNull(a)) and (a.var == e3.var));
    
    addBrVertex(br, vl);
    addBrVertex(br, vr);
    
    addBrEdge(br, edge(vl, e0.t, e0.var), 0, p0);
    addBrEdge(br, edge(vr, e1.t, e1.var), 1, p1);
    addBrEdge(br, edge(e2.s, vr, e2.var), p2, 2);
    addBrEdge(br, edge(e3.s, vl, e3.var), p3, 3);
    
    removeBrVertex(br, v);
);

-- joins the crossing x in place
joinCrossing = method();
joinCrossing(BraidRes, ZZ) := (br, x) -> (
    v := vertex(x, 0);
    vl := vertex(x, -1);
    vr := vertex(x, 1);
        
    e0 := br.adjacent#vl#0;
    p0 := position((br.adjacent)#(e0.t), a -> (not isNull(a)) and (a.var == e0.var));
    e1 := br.adjacent#vr#1;
    p1 := position((br.adjacent)#(e1.t), a -> (not isNull(a)) and (a.var == e1.var));
    e2 := br.adjacent#vr#2;
    p2 := position((br.adjacent)#(e2.s), a -> (not isNull(a)) and (a.var == e2.var));
    e3 := br.adjacent#vl#3;
    p3 := position((br.adjacent)#(e3.s), a -> (not isNull(a)) and (a.var == e3.var));
    
    addBrVertex(br, v);
    
    addBrEdge(br, edge(v, e0.t, e0.var), 0, p0);
    addBrEdge(br, edge(v, e1.t, e1.var), 1, p1);
    addBrEdge(br, edge(e2.s, v, e2.var), p2, 2);
    addBrEdge(br, edge(e3.s, v, e3.var), p3, 3);
    
    removeBrVertex(br, vl);
    removeBrVertex(br, vr);
);