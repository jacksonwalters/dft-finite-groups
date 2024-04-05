#Theorem 3.3, "Nakayama's Conjecture", pg. 8, Wildon, Vertices of Specht Modules and Blocks of the Symmetric Group:

#Let p be a prime. The p-blocks of the symmetric group Sn are labelled by pairs (\gamma,w), 
#where \gamma is a p-core and w \in \mathbb{N}_0 is the associated weight, such that |\gamma| + wp = n. 
#Thus S^\lambda lies in the block labelled by (\gamma,w) if and only if \lambda has p-core \gamma and weight w.
#
#---
#
#Theorem 2.1, Murphy, The Idempotents of the Symmetric Group and Nakayama's Conjecture
#
#F_j^i = \sum_{\mu,k | t_k^\mu \in T_j^i} E_k^\mu 
#H^i = \sum_j F_j^i 
#
#{F_j^i} is a complete set of orthogonal idempotents of RS_n, similarly for {\overline{F}_j^i}.
#
#where
#
#E_i^\mu = \prod_{c = -n+1}^{n-1} \prod_{u | a_{ui}^\mu \ne c} \frac{c-L_u}{c-\alpha_{ui}^\mu}
#
#L_u = (1,u) + (2,u) + ... + (u-1,u), a sum of transpositions
#
#T_j^i are tableaux belonging to equivalence classes B_i, where B_i are equivalences classes of
#partitions given by \tilde_p
#
#\labmda \tidle_p \mu iff the classes (j-i) of each node have the same residue mod p
#
#The class of the node occupied by u in t_i^\mu is \alpha_{ui}^\mu where t_1^\mu, ..., t_d^\mu are the standard \mu-tableaux 
#in the ordering in pg. 288 of "A new construction of Young's seminormal representation of the symmetric groups".
#
#The "class" of a node (i,j) in a Young diagram \mu is the difference j-i.
#
#d is the number of "standard \mu-tableaux", given by the Hook-length formula.
#
#---
#
#Theorem 2.8 (Nakayama's Conjecture): {\overline{H}^i} is a complete set of primitive orthogonal central idempotents of \overline{K}S_n,
#and S_{\overline{K}}^\mu, S_{\overline{K}}^\lambda belong to the same block of \overline{K}S_n if and only if \mu \tilde_p \lambda.
#
#---

#define the element \alpha_{ui}^\mu
#which is the class of u in the i^th standard Tableaux corresponding to mu
#the class of a node is j-i where i is the row and j is the column
def alpha(mu,u,i):
    tab = StandardTableaux(mu)[i]
    for row in range(len(tab)):
        for col in range(len(tab[row])):
            if u == tab[row][col]:
                return col - row
            
#define the elment L_u = (1,u) + (2,u) + ... + (u-1,u) in K[S_n]
def L(u,n,char=0):
    if char == 0:
        SGA = SymmetricGroupAlgebra(QQ,n)
    else:
        SGA = SymmetricGroupAlgebra(GF(char),n)
    return sum(SGA(PermutationGroupElement(f"({i},{u})")) for i in range(1,u))

#from Murphy
def E(mu,i,n,char=0):
    return prod(prod((c-L(u,n,char))/(c-alpha(mu,u,i)) for u in range(1,n+1) if alpha(mu,u,i) != c) for c in range(-n+1,n))

#helper function to determine if diagram comes from a partition
def diagram_is_from_partition(diag):
    #handle empty diagram
    if len(diag) == 0:
        return True
    from_partition = True
    #check if rows start at 0 and have any gaps
    row_labels = {node[0] for node in diag}
    from_partition &= min(row_labels) == 0
    from_partition &= len(row_labels) == max(row_labels) - min(row_labels) + 1
    #for each row, check if columns have any gaps, and start at 0
    part = []
    for row in row_labels:
        col_labels = {node[1] for node in diag if node[0]==row}
        part.append(len(col_labels))
        from_partition &= min(col_labels) == 0
        from_partition &= len(col_labels) == max(col_labels) - min(col_labels) + 1
    from_partition &= all(part[i] >= part[i+1] for i in range(len(part) - 1))
    return from_partition

#given a diagram which comes from a partition, give the partition
def to_partition(diag):
    assert diagram_is_from_partition(diag)
    return Partition([max({node[1] for node in diag if node[0]==row})+1 for row in {node[0] for node in diag}])

#find the p-core of a partition
#a "rim p-hook" is:
#connected set of p nodes
#in the rim (the node (i+1,j+1) is not in the diagram)
#whose removal leaves a valid young diagram
#this re-implements mu.core(p)
from sage.combinat.diagram import Diagram, NorthwestDiagrams
def core(mu,p):
    #get the diagram associated to the partition
    diag = NorthwestDiagrams().from_partition(mu)
    #find the rim of the diagram/partition
    rim = []
    for node in diag:
        if (node[0]+1,node[1]+1) not in diag:
            rim.append(node)
    #find rim p-hooks
    p_hooks = []
    for start in rim:
        p_hook = [start]
        for step in range(p-1):
            current = p_hook[step]
            down = (current[0]+1,current[1])
            left = (current[0],current[1]-1)
            if (not down in rim) and (not left in rim):
                break
            if down in rim:
                p_hook.append(down)
                continue
            if left in rim:
                p_hook.append(left)
                continue
        if len(p_hook) == p:
            p_hooks.append(p_hook)
    #determine if removal of candidate p_hook results in a diagram
    valid_p_hook_removed_diags = []
    for p_hook in p_hooks:
        diag_minus_p_hook = Diagram([node for node in diag if node not in p_hook])
        if diagram_is_from_partition(diag_minus_p_hook):
            valid_p_hook_removed_diags.append(diag_minus_p_hook)
    #recursive step: if no p-hooks, return original diagram
    if len(valid_p_hook_removed_diags) == 0:
        return to_partition(diag)
    #recursive step: if diagram is unchanged, return, otherwise keep going
    for valid_p_hook_removed_diag in valid_p_hook_removed_diags:
        partition = to_partition(valid_p_hook_removed_diag)
        return core(partition,p)
    
    #determine when two partitions are p-equivalent, i.e. have the same p-core
#look at canonical tableaus (place 1,...,n in top-to-bottom, left-to-right)
#compute p-residues of classes j-i \mod p
#using the method described in Littlewood '51
#just check if the sets have the same residue classes
#BUG: with "residue" method, the same partition with different zero padding will have
#different residue sets
from multiset import *
def p_equiv(mu_1, mu_2,p,method="residue"):
    if method=="core":
        return core(mu_1,p) == core(mu_2,p)
    if method=="residue":
        max_len = max(len(mu_1),len(mu_2))
        #zero pad partitions
        mu_1 += [0]*(max_len - len(mu_1))
        mu_2 += [0]*(max_len - len(mu_2))
        #compute residue sets
        residues_1 = set([(mu_1[i]+max_len-1-i) % p for i in range(max_len)])
        residues_2 = set([(mu_2[i]+max_len-1-i) % p for i in range(max_len)])
        return residues_1 == residues_2
    if method=="multiset":
        diag_1=NorthwestDiagrams().from_partition(mu_1)
        diag_2=NorthwestDiagrams().from_partition(mu_2)
        res_1 = Multiset([(node[1]-node[0]) % p for node in diag_1])
        res_2 = Multiset([(node[1]-node[0]) % p for node in diag_2])
        k=0
        while res_1 >= Multiset({i:k for i in range(5)}):
            k += 1
        reduced_1 = res_1 - Multiset({i:k-1 for i in range(5)})
        k=0
        while res_2 >= Multiset({i:k for i in range(5)}):
            k += 1
        reduced_2 = res_2 - Multiset({i:k-1 for i in range(5)})
        return reduced_1 == reduced_2
    
#equivalence classes of all partitions under \tilde_p
#two partitions are equivalent if they have the same p-core
def B(n,p):
    equiv_classes = []
    for part_1 in Partitions(n):
        equiv_class_found = False
        for equiv_class in equiv_classes:
            for part_2 in equiv_class:
                if p_equiv(part_1,part_2,p,method="core"):
                    equiv_class.append(part_1)
                    equiv_class_found = True
                    break
        if not equiv_class_found:
            equiv_classes.append([part_1])
    return equiv_classes

#in each equivalence class B_i, we have p-equivalent partitions
#each partition in B_i has associated standard tableaux
#we obtain equivalences classes T_j^i by extending ~_p from partitions to tableaux
#define t ~_p t^* if every u <= n occupies a node of the same p-class in both t and t^*
def equiv_tableaux(t1,t2,n,p):
    equiv = True
    #find index of u in t1 and t2
    for u in range(1,n+1):
        for i in range(len(t1)):
            for j in range(len(t1[i])):
                if t1[i][j] == u:
                    res_1 = (j-i) % p
        for i in range(len(t2)):
            for j in range(len(t2[i])):
                if t2[i][j] == u:
                    res_2 = (j-i) % p
        equiv &= res_1 == res_2
    return equiv

#flatten a list
def flatten(l):
    return [item for sublist in l for item in sublist]

#equivalence classes of all partitions under \tilde_p
#two partitions are equivalent if they have the same p-core
#from sage.combinat.tableau import IncreasingTableaux_shape
def T(i,n,p):
    partition_equiv_class = B(n,p)[i]
    equiv_classes_list = []
    #should this be StandardTableaux, or allTableaux?
    all_tableaux = flatten([StandardTableaux(mu) for mu in partition_equiv_class])
    equiv_classes = []
    for tab_1 in all_tableaux:
        equiv_class_found = False
        for equiv_class in equiv_classes:
            for tab_2 in equiv_class:
                if equiv_tableaux(tab_1,tab_2,n,p):
                    equiv_class.append(tab_1)
                    equiv_class_found = True
                    break
        if not equiv_class_found:
            equiv_classes.append([tab_1])
    equiv_classes_list.append(equiv_classes)
    return flatten(equiv_classes_list)

#define idempotents F_i^j
#CLARIFY: flattening and summing over all elements in T(i,n,p)[j] works
#however it's not clear why, i.e. F_ij no longer depends on j
#CLARIFICATION: since we are summing orthogonal idempotents to get another idempotent, they can't be primitive
#therefore we should abandon the j-indices and just sum
def F(i,n,p,j=-1):
    F_ij = 0
    if j >= 0:
        tableaux = T(i,n,p)[j]
    else:
        tableaux = flatten(T(i,n,p))
    for t in tableaux:
        mu = [len(l) for l in t]
        std_tableaux = StandardTableaux(mu)
        k = list(std_tableaux).index(Tableau(t))
        F_ij += E(mu,k,n)
    return F_ij

#compute idempotents for symmetric group over finite field
def central_orthogonal_idempotents(n,p):
    num_p_equiv_partition_classes = len(B(n,p))
    idempotents = []
    for i in range(num_p_equiv_partition_classes):
        idem = sum([GF(p)(item[1])*SymmetricGroupAlgebra(GF(p),n)(item[0]) for item in F(i,n,p,j=-1)])
        idempotents.append(idem)
    return idempotents