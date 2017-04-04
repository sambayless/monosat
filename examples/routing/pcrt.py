#Read RUC-style PCRT constrained routing files, based on those from Alexander Nadel's  "Routing Under Constraints" (FMCAD 2016).
#From README:
#		; 2D Grid definition
#		G X Y; Comments
# 		; Nets
# 		N vid1 vid2 ... vidn; Comments
# 		; Constraints
# 		C vid1 vid2 ... vidn; Comments
#		; Disabled vertices
#		D vid1; Comments

#Returns a tuple, ((GridWidth,GridHeight),(Tuples of vertices to route), (Tuples of vertices that cannot be used simultaneously))
def read(filename_or_stream):

    if isinstance(filename_or_stream, str):
        file = open(filename_or_stream,"r", encoding="utf-8")
    else:
        file = filename_or_stream
    X=None
    Y=None
    all_netids=[]
    all_constraints=[]
    all_disabled=[]
    allow_diagonals=False

    def toVertex(vID):
        # A vertex ID for (x,y) is y*X+x.
        assert(X is not None)
        assert (Y is not None)
        assert(vID>=0)
        y = vID//X
        x = vID - (y*X)
        assert(x<=X)
        assert(x>=0)
        assert(y>=0)
        assert(y<=Y)
        return (x,y)

    try:
        for line in file:
            line = line.rstrip()
            if line.startswith(";"):
                continue
            line = line.partition(';')[0].strip()
            if (len(line)==0):
                continue
            parts = line.split()
            if parts[0] == "G":
                #this is a 2D grid definition
                #A 2D grid (X,Y) where a vertex is associated with every grid point and an edge is created between every two (neighbouring) vertices
                assert(len(parts)==3 or len(parts)==4)
                X = int(parts[1])
                Y = int(parts[2])
                if len(parts)>=4 and parts[3]=='45':
                    #this file supports 45 degree routing
                    allow_diagonals=True
            elif parts[0] == "N":
                #this is a net definition - a list of vertex IDs to be connected together.
                netids =  list(map(toVertex, map(int,parts[1:])))
                all_netids.append(netids)
            elif parts[0] == "C":
                #constraints definition, which is a list of one or more vertex IDs that cannot be used simultaneously.
                #A vertex ID for (x,y) is y*X+x.
                cids = list(map(toVertex, map(int,parts[1:])))
                all_constraints.append(cids)
            elif parts[0] == "D":
                #a disabled vertex
                assert(len(parts)==2)
                all_disabled.append(toVertex(int(parts[1])))

    finally:
        file.close()
    return ((X,Y),allow_diagonals,all_netids,all_constraints,all_disabled)

if __name__ == '__main__':
    import sys
    print(read(sys.argv[1]))
