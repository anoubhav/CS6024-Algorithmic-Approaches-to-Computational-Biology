import sys

""" 
Pseudocode
 NeighborJoining(D,n)
 if n = 2
  T ← tree consisting of a single edge of length D1,2
  return T
 D' ← neighbor-joining matrix constructed from the distance matrix D
 find elements i and j such that D'i,j is a minimum non-diagonal element of D'
 Δ ← (TotalDistanceD(i) - TotalDistanceD(j)) /(n - 2)
 limbLengthi ← (1/2)(Di,j + Δ)
 limbLengthj ← (1/2)(Di,j - Δ)
 add a new row/column m to D so that Dk,m = Dm,k = (1/2)(Dk,i + Dk,j - Di,j) for any k
 remove rows i and j from D
 remove columns i and j from D
 T ← NeighborJoining(D, n - 1)
 add two new limbs (connecting node m with leaves i and j) to the tree T
 assign length limbLengthi to Limb(i)
 assign length limbLengthj to Limb(j)
 return T 
 """

class Tree():
    """ Tree class
    N: Number of nodes,
    bidrectional: flag variable indicates undirected/directed tree 
    """
    def __init__(self,N=0,bidirectional=True):
        self.nodes=list(range(N))
        self.edges={}
        self.bidirectional=bidirectional
        self.N = N

    def link(self,start,end,weight=1): 
        """ Link two nodes to form an edge """
        self.half_link(start,end,weight)
        if self.bidirectional:
            self.half_link(end,start,weight)

    def half_link(self,a,b,weight=1):
        if a not in self.nodes:
            self.nodes.append(a)        
        if a in self.edges:               
            self.edges[a]=[(b0,w0) for (b0,w0) in self.edges[a] if b0 != b] + [(b,weight)]
        else:
            self.edges[a]=[(b,weight)]

    def print_adjacency(self):
        """ Print the adjacency list; start node sorted then end node """
        self.nodes.sort()
        for node in self.nodes:
            if node in self.edges:
                self.edges[node].sort()
                for edge in self.edges[node]:
                    end, weight = edge
                    print ('{0}->{1}:{2}'.format(node, end, round(weight, 3)))

def NeighborJoining(D,n,node_list=None):
    """ Function implements neighbor joining algorithm
    D: Distance matrix of size NxN
    n: N  
    """
    def remove(i,D):
        """ Remove row and column i from distance matrix D """
        D_new=[]
        for j in range(len(D)):
            if j!=i:
                D_new.append([D[j][k] for k in range(len(D[j])) if k!=i])
        return D_new        
        
    def create_DPrime(Total_distance):
        """ Construct the neighbour joining matrix from distance matrix. Return i, j such that DPrime i,j is a minimum non-diagonal element of DPrime """

        DPrime=[[0]*n for _ in range(n)]

        min_i, min_j, min_D = 0, 0, 10*8    # setting the minimum distance as arbitrary large number

        for i in range(n):
            for j in range(i+1,n):
                DPrime[i][j] = (n-2)*D[i][j] - Total_distance[i] - Total_distance[j]
                DPrime[j][i] = DPrime[i][j]
                if DPrime[i][j] < min_D:
                    min_i, min_j, min_D = i, j, DPrime[i][j]
                    
        return DPrime, min_i, min_j
    
    def create_Delta(Total_distance):
        """ Δ ← (TotalDistanceD(i) - TotalDistanceD(j)) /(n - 2) """
        Delta=[[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1,n):
                Delta[i][j] = (Total_distance[i]-Total_distance[j])/(n-2)
                Delta[j][i] = Delta[i][j]
        return Delta

    if node_list==None:
        node_list=list(range(n))
        
    if n==2:
        """ Base case """
        T=Tree()
        T.link(node_list[0],node_list[1],D[0][1])
        return T

    else:
        Total_distance = [sum(row) for row in D]
        _, i, j = create_DPrime(Total_distance)

        Delta = create_Delta(Total_distance)

        limbLength_i = (D[i][j]+Delta[i][j])/2
        limbLength_j = (D[i][j]-Delta[i][j])/2

        new_row=[0.5*(D[k][i]+D[k][j]-D[i][j]) for k in range(n)]+[0]
 
        D.append(new_row)
        for l in range(n):
            D[l].append(new_row[l])

        m = node_list[-1]+1
        node_list.append(m)

        D = remove(max(i,j),D)
        D = remove(min(i,j),D)
        node_i = node_list[i]
        node_j = node_list[j]

        node_list.remove(node_i)
        node_list.remove(node_j)

        T = NeighborJoining(D, n-1, node_list)

        T.link(node_i,m,limbLength_i)
        T.link(node_j,m,limbLength_j)       
        return T

def main():
    dist_matrix = list()
    flag = 1
    for row in sys.stdin:
        if row == '\n': break
        if flag: 
            flag = 0
            n = int(row)
            continue
        dist_matrix.append([int(i) for i in row.split()])

    NeighborJoining(dist_matrix, n).print_adjacency()

if __name__ == '__main__':
    main()