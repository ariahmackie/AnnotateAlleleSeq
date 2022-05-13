A = [((148556137, 148567814), 'uc011mxd.2'), ((148556137, 148582727), 'uc011mxg.2'), ((148556137, 148610433), 'uc010nsu.2'), ((148556137, 148611311), 'uc004fcw.4')]
B = [((138804993, 138847077), 'uc004fax.3'), ((138804993, 138910928), 'uc004fba.3')]
C = [((134416301, 134426538), 'uc004eyo.3')]
D = []


#print(len(A)) # 4
#print(len(B)) #2
#print(len(C)) # 1
#print(len(D)) # 0

#print(A[0]) #((148556137, 148567814), 'uc011mxd.2')
#print(A[1]) #((148556137, 148582727), 'uc011mxg.2')
#print(A[0][1]) #'uc011mxd.2'

def get_gene_names(genes_covering_snp):
    if len(genes_covering_snp) == 0:
        return "No Gene Match Found."
    gene_names = "Matches: "
    for i in range(len(genes_covering_snp)):
            if i == len(genes_covering_snp) - 1:
                gene_names = gene_names + genes_covering_snp[i][1]
            else:
                gene_names = gene_names +  genes_covering_snp[i][1] + ", "
    return gene_names

print(get_gene_names(A))
print(get_gene_names(B))
print(get_gene_names(C))
print(get_gene_names(D))
