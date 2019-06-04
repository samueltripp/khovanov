isNull := method(Dispatch=>Thing);
isNull(Nothing) := x -> true;
isNull(Thing) := x -> false;
indexOf := (pred, l) -> (
    for i from 0 to #l-1 do if (not isNull(l#i) and pred(l#i)) then return i;
);
addEdge := (ht, s, t, ps, pt, v) -> (
    if not ht#?s then ht#s = new MutableList from {,,,};
    if not ht#?t then ht#t = new MutableList from {,,,};
    (ht#s)#ps = (t, v);
    (ht#t)#pt = (s, v);
);

Braid = new Type of HashTable;
braid = method();
braid(ZZ, List) := (n, word) -> (
    return new Braid from hashTable{(global n)=>n, (global word)=>word};
);

BraidRes = new Type of HashTable;
braidRes = method();
braidRes(Digraph, HashTable) := (g, edges) -> (
    return new BraidRes from hashTable{(global graph) => g, (global edges) => edges};
);

-- build the fully singular resolution of the given braid
-- returns the base ring and the resolution
-- Braid -> Ring x Resolution
singularResolution = method();
singularResolution(Braid) := (b) -> (
    n := b.n;
    word := b.word;
    R := QQ[for i from 1 to 6*n-4+2*#word list "u"|i];

    -- Vertex -> {{Vertex,RingElement}}
    edgeVars := new MutableHashTable;
     
    sing := digraph({});

    -- add vertices
    for i from 1 to n do(
	sing = addVertex(sing, (-3,i));
    );
    strands = new MutableList;
    for i from 0 to 2*n-1 do(
	strands#i = ((-3,i//2+1), i%2+2);
    );
    nextEdgeIndex = 2*n;
    for i from 1 to n-1 do(
	v = (-2,i);
	(v1, p1) = strands#(2*i-1);
	(v2, p2) = strands#(2*i);
	sing = addVertex(sing, v);
	sing = addEdges'(sing, {{v, v1}, {v, v2}});
	addEdge(edgeVars, v, v1, 0, p1, R_nextEdgeIndex);
	nextEdgeIndex = nextEdgeIndex + 1;
	addEdge(edgeVars, v, v2, 1, p2, R_nextEdgeIndex);
	nextEdgeIndex = nextEdgeIndex + 1;
	strands#(2*i-1) = (v, 2);
	strands#(2*i) = (v, 3);
    );
    for i from 1 to n-1 do(
	v = (-1,i);
	(v1, p1) = strands#(2*i);
	(v2, p2) = strands#(2*i+1);
	sing = addVertex(sing, v);
	sing = addEdges'(sing, {{v, v1}, {v, v2}});
	addEdge(edgeVars, v, v1, 0, p1, R_nextEdgeIndex);
	nextEdgeIndex = nextEdgeIndex + 1;
	addEdge(edgeVars, v, v2, 1, p2, R_nextEdgeIndex);
	nextEdgeIndex = nextEdgeIndex + 1;
	strands#(2*i) = (v, 2);
	strands#(2*i+1) = (v, 3);
    );
    for i from 0 to #word-1 do(
	v = (i,0);
	(v1, p1) = strands#(abs(word#i)-1);
	(v2, p2) = strands#(abs(word#i));
	sing = addVertex(sing, v);
	sing = addEdges'(sing, {{v, v1}, {v, v2}});
	addEdge(edgeVars, v, v1, 0, p1, R_nextEdgeIndex);
	nextEdgeIndex = nextEdgeIndex + 1;
	addEdge(edgeVars, v, v2, 1, p2, R_nextEdgeIndex);
	nextEdgeIndex = nextEdgeIndex + 1;
	strands#(abs(word#i)-1) = (v, 2);
	strands#(abs(word#i)) = (v, 3);
    );
    for i from 0 to #strands-1 do(
	v = (-3, i//2+1);
	(v1, p1) = strands#i;
	sing = addEdges'(sing, {{v, v1}});
	addEdge(edgeVars, v, v1, i%2, p1, R_i);
    );
    
    return (R, braidRes(sing, edgeVars));
);

splitCrossing = method();
splitCrossing(BraidRes, ZZ) := (br, x) -> (    
    g := br.graph;
    edges := br.edges;
    v = (x, 0);
    vl = (x, -1);
    vr = (x, 1);
    
    (w0, var0) = edges#v#0;
    p0 = indexOf((a,b) -> b == var0, edges#w0);
    (w1, var1) = edges#v#1;
    p1 = indexOf((a,b) -> b == var1, edges#w1);
    (w2, var2) = edges#v#2;
    p2 = indexOf((a,b) -> b == var2, edges#w2);
    (w3, var3) = edges#v#3;
    p3 = indexOf((a,b) -> b == var3, edges#w3);
        
    addEdge(edges, vl, w0, 0, p0, var0);
    addEdge(edges, vr, w1, 1, p1, var1);
    addEdge(edges, vl, w2, 2, p2, var2);
    addEdge(edges, vr, w3, 3, p3, var3);
    
    g = addVertex(g, vl);
    g = addVertex(g, vr);
    g = addEdges'(g, {{vl, w0}, {vr, w1}, {w2, vl}, {w3, vr}});
    
    g = deleteVertex(g, v);
    remove(edges, v);
    
    return braidRes(g, edges);
);

joinCrossing = method();
joinCrossing(BraidRes, ZZ) := (br, x) -> (    
    g := br.graph;
    edges := br.edges;
    v = (x, 0);
    vl = (x, -1);
    vr = (x, 1);
    
    (w0, var0) = edges#vl#0;
    p0 = indexOf((a,b) -> b == var0, edges#w0);
    (w1, var1) = edges#vr#1;
    p1 = indexOf((a,b) -> b == var1, edges#w1);
    (w2, var2) = edges#vl#2;
    p2 = indexOf((a,b) -> b == var2, edges#w2);
    (w3, var3) = edges#vr#3;
    p3 = indexOf((a,b) -> b == var3, edges#w3);
        
    addEdge(edges, v, w0, 0, p0, var0);
    addEdge(edges, v, w1, 1, p1, var1);
    addEdge(edges, v, w2, 2, p2, var2);
    addEdge(edges, v, w3, 3, p3, var3);
    
    g = addVertex(g, v);
    g = addEdges'(g, {{v, w0}, {v, w1}, {w2, v}, {w3, v}});
    
    g = deleteVertex(g, vl);
    g = deleteVertex(g, vr);
    remove(edges, vl);
    remove(edges, vr);
    
    return braidRes(g, edges);
);

(R, sing) = singularResolution(braid(2,{2,2,2}));
sing = splitCrossing(sing, 1);
sing = joinCrossing(sing, 1);
