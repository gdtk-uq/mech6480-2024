-- Create the nodes that define key points for our geometry.
A = Vector3:new{x=-0.07, y=0.0};
B = Vector3:new{x=-0.05, y=0.0}
C = Vector3:new{x=0.0, y=0.0};
D = Vector3:new{x=0.005, y=0.012}
E = Vector3:new{x=0.1, y=0.03};
F = Vector3:new{x=0.202, y=0.03}
G = Vector3:new{x=0.207, y=0.0};
H = Vector3:new{x=0.3, y=0.0}
I = Vector3:new{x=-0.07, y=0.1};
J = Vector3:new{x=-0.05, y=0.1}
K = Vector3:new{x=0.1, y=0.1};
L = Vector3:new{x=0.202, y=0.1}
M = Vector3:new{x=0.3, y=0.1};
N = Vector3:new{x=0.3, y=0.03}

-- Some interior Bezier control points
CD_b1 = Vector3:new{x=0.0, y=0.006}
GF_b1 = Vector3:new{x=0.207, y=0.027}
DE_b1 = Vector3:new{x=0.0064, y=0.012}
DE_b2 = Vector3:new{x=0.0658, y=0.0164}
DE_b3 = Vector3:new{x=0.0727, y=0.0173}


-- Paths defining the beer bottle shape
CD = Bezier:new{points={C, CD_b1, D}} -- top of bottle
DE = Bezier:new{points={D, DE_b1, DE_b2, DE_b3, E}} -- neck of bottle
EF = Line:new{p0=E, p1=F} -- side of bottle
GF = ArcLengthParameterizedPath:new{
   underlying_path=Bezier:new{points={G, GF_b1, F}}}
GH = Line:new{p0=G, p1=H}

-- Upper boundary of domain
KL = Line:new{p0=K, p1=L}
LM = Line:new{p0=L, p1=M} 

-- Dividing Paths (interior block joins)
EK = Line:new{p0=E, p1=K}
FL = Line:new{p0=F, p1=L} 
FN = Line:new{p0=F, p1=N}

-- Paths at outflow
HN = Line:new{p0=H, p1=N}
NM = Line:new{p0=N, p1=M}

-- Define the blocks, boundary conditions and set the discretisation.
n0 = 10; n1 = 4; n2 = 20; n3 = 20; n4 = 20; n5 = 12; n6 = 8

patch = {}
-- patch[0] = 
-- patch[1] = 
-- patch[2] = 
patch[3] = CoonsPatch:new{north=KL, east=FL, south=EF, west=EK}
patch[4] = CoonsPatch:new{north=LM, east=NM, south=FN, west=FL}
patch[5] = CoonsPatch:new{north=FN, east=HN, south=GH, west=GF}

grid = {}
-- grid[0] = 
-- grid[1] = 
-- grid[2] = 
grid[3] = StructuredGrid:new{psurface=patch[3], niv=n4+1, njv=n2+1}
grid[4] = StructuredGrid:new{psurface=patch[4], niv=n5+1, njv=n2+1}
grid[5] = StructuredGrid:new{psurface=patch[5], niv=n5+1, njv=n6+1}

for ib = 3, 5 do
   fileName = string.format("beer-bottle-blk-%d.vtk", ib)
   grid[ib]:write_to_vtk_file(fileName)
end

