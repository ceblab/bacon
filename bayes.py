import numpy as np
import itertools
import copy

from regex import W

np.random.seed(314)

class Bayes:
    num_variable = 0
    max_parent_set = 0
    score_disc_list = list()
    


    def __init__(self,num_vari,num_parent):
        self.num_variable = num_vari
        self.max_parent_set = num_parent
        self.score_disc_list = [{(i,) : self.score({i})} for i in range(num_vari)]#key tuple(i) -> parent set ga empty set no  baai 
        self.make_score_disc()


    def output(self):
        print(self.num_variable)
    
    def score(self,subset):
        return np.random.randint(0,1023)

    def make_score_disc(self):
        for num_vari in range(self.num_variable):
            vari_list = [i for i in range(self.num_variable) if i != num_vari]

            for m in range(1, self.max_parent_set + 1):
                for conb in itertools.combinations(vari_list,m ):
                
                    self.score_disc_list[num_vari][conb] = self.score(set(conb))


def make_W_set_of_i(disc_of_i,child_vari,num_variable,max_parent_set):
    W_set_of_i = list()
    W_set_of_i.append((child_vari,))
    score_of_empty_set = disc_of_i[(child_vari,)]
    max_score_list = [score_of_empty_set for i in range(num_variable)]#score of parent set i
    for i in range(num_variable):
        if disc_of_i[(i,)] > score_of_empty_set:
            W_set_of_i.append((i,))
            max_score_list[i] = disc_of_i[(i,)]
    
    vari_list = [i for i in range(num_variable) if i != child_vari]
    
    for conb in itertools.combinations(vari_list,2):
        a,b = conb
        s = disc_of_i[conb]
        if s>max_score_list[a] and s> max_score_list[b]:
            W_set_of_i.append(conb)
    
    for conb in itertools.combinations(vari_list,3):
        a,b,c = conb
        s = disc_of_i[conb]
        if s>max_score_list[a] and s> max_score_list[b] and s> max_score_list[c] and s>disc_of_i[(a,b)] and s>disc_of_i[(a,c)] and s>disc_of_i[(b,c)]:
            W_set_of_i.append(conb)
    return W_set_of_i
            
def decomposition_of_i(w_set_in,child_vari,max_parent_set):
    W_set = {i for i in w_set_in}
    W_set.remove((child_vari,))
    U_set = [i for i in w_set_in if len(i) == 1]
    W_set = W_set - make_union_set(U_set,max_parent_set)
    print("W_set",W_set)#
    V_set = set()
    for w in W_set:
        for i in range (1,len(w)+1):
            for conb in itertools.combinations([j for j in w],i):
                if not (conb in U_set) :
                    V_set.add(conb)
    k = len(V_set)
    print("V_set",V_set)#
    while(len(W_set) > 0):
        W_set_lis = [{i for i in j} for j in W_set]
        V_set_dust = set()
        for i in V_set:
            v_in_count = 0
            v_set_i = {k for k in i}
            for j in W_set_lis:
                if v_set_i <= j:
                    v_in_count=v_in_count+1
            if v_in_count <= 1:
                V_set_dust.add(i)
        V_set = V_set - V_set_dust
        print("V_set",V_set)#
        W_set_dash = set()
        v_dash = ()
        k_dash = 0
        for v in V_set:
            W_set_2dash = W_set & make_union_set_with_v(U_set,max_parent_set,v)
            if len(W_set_2dash) >len(W_set_dash):
                W_set_dash = W_set_2dash
                v_dash = v
                k_dash = len(W_set_2dash)
            if k == k_dash:
                break
        if len(W_set_dash) > 0:
            U_set.append(v_dash)
            V_set.remove(v_dash)
            W_set = W_set - W_set_dash
            W_set_lis = [{i for i in j} for j in W_set]
            v_set_i = {k for k in v_dash}
            v_in_count = 0
            for j in W_set_lis:
                if v_set_i <= j:
                    v_in_count=v_in_count+1
            if v_in_count > 0:
                k = len(W_set_dash) + 1
            else:
                k = len(W_set_dash)
        else:
            U_set.extend([i for i in W_set])
            W_set = set()
    return U_set

def make_union_set(U_set,max_parent_set):
    union_U_set = set()
    l = len(U_set)
    U_set_tuple_to_set = [{i for i in j} for j in U_set]
    for i in range(l):
        for j in range(i,l):
            set1 = U_set_tuple_to_set[i] | U_set_tuple_to_set[j]
            list1 = [i for i in set1]
            if len(list1) < max_parent_set+1:
                set2 = tuple(np.sort(list1))
                union_U_set.add(set2)

            

            
    return union_U_set

def make_union_set_with_v(U_set,max_parent_set,v):
    union_U_set = set()
    l = len(U_set)
    v_set = {i for i in v}
    U_set_tuple_to_set = [{i for i in j} for j in U_set]
    for i in range(l):
        set1 = U_set_tuple_to_set[i] | v_set
        list1 = [i for i in set1]
        if len(list1) < max_parent_set+1:
            set2 = tuple(np.sort(list1))
            union_U_set.add(set2)
            
    return union_U_set

num_variable = 9
max_parent = 3
bayes = Bayes(num_variable,max_parent)
bayes.output()
#s = set((1,2,3))
score = bayes.score_disc_list[1][(1,)]
print(score)
score = bayes.score_disc_list[1][(2,3,8)]
print(score)

w_list = make_W_set_of_i(bayes.score_disc_list[1],1,num_variable,max_parent)
print("w_list",w_list)
#for i in w_list:
#    print(i ,": score" ,bayes.score_disc_list[1][i])


p_list = decomposition_of_i(w_list,1,max_parent)
print("p_list",p_list)