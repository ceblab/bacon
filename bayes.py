from calendar import c
from re import L
import numpy as np
import itertools
import copy
import math
np.set_printoptions(threshold=np.inf)
from regex import W

np.random.seed(314)

class Bayes:
    num_variable = 0
    max_parent_set = 0
    score_disc_list = list()
    u_set_list = list()
    count_len_of_u = list()


    def __init__(self,num_vari,num_parent):
        self.num_variable = num_vari
        self.max_parent_set = num_parent 
        self.u_set_list = [list() for i in range(self.num_variable)]
        self.init_score_disc()

    def init_score_disc(self):
        self.score_disc_list = [dict() for i in range(self.num_variable)]

        
    def output(self):
        print(self.num_variable)
    
    def score(self,subset):
        return np.random.randint(0,1023)

    def make_score_disc_i(self,child_vari):
        self.score_disc_list[child_vari][(child_vari,)] = self.score(set((child_vari,)))
        vari_list = [i for i in range(self.num_variable) if i != child_vari]

        for m in range(1, self.max_parent_set + 1):
            for conb in itertools.combinations(vari_list,m ):
                
                self.score_disc_list[child_vari][conb] = self.score(set(conb))


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
    
    V_set = set()
    for w in W_set:
        for i in range (1,len(w)+1):
            for conb in itertools.combinations([j for j in w],i):
                if not (conb in U_set) :
                    V_set.add(conb)
    k = len(V_set)
    
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

def calc_s_of_i(U_set,W_set,disc_of_i,child_vari):
    s_disc = dict()
    l = len(U_set)
    empty_score = math.log(disc_of_i[(child_vari,)])
    for i in range(l):
        for j in range(i,l):
            if i ==0 & j == 0:
                pair_ij = (child_vari,)
            else:
                set_i= {k for k in U_set[i]}
                set_j = {k for k in U_set[j]}
                set_ij = (set_i | set_j) - {child_vari}            
                pair_ij = tuple(np.sort([k for k in set_ij]))
            if i == 0:
                s_disc[(child_vari,i,j)] = 0


            elif(i == j):
                if(pair_ij in W_set):
                    s_disc[(child_vari,i,j)] = -math.log(disc_of_i[pair_ij]) + empty_score
                else:
                    s_disc[(child_vari,i,j)] = 0
            else:
                if(pair_ij in W_set):
                    if(U_set[i] in W_set):
                        if(U_set[j] in W_set):
                            s_disc[(child_vari,i,j)] = - math.log(disc_of_i[pair_ij]) + math.log(disc_of_i[U_set[i]]) +math.log(disc_of_i[U_set[j]]) -empty_score
                        else:
                            s_disc[(child_vari,i,j)] = - math.log(disc_of_i[pair_ij]) + math.log(disc_of_i[U_set[i]])
                    else:
                        if(U_set[j] in W_set):
                            s_disc[(child_vari,i,j)] = - math.log(disc_of_i[pair_ij]) +math.log(disc_of_i[U_set[j]])
                        else:
                            s_disc[(child_vari,i,j)] = - math.log(disc_of_i[pair_ij]) + empty_score
                else:
                    if(U_set[i] in W_set):
                        if(U_set[j] in W_set):
                            s_disc[(child_vari,i,j)] = math.log(disc_of_i[U_set[i]]) +math.log(disc_of_i[U_set[j]]) - 2 * empty_score
                        else:
                            s_disc[(child_vari,i,j)] = - empty_score + math.log(disc_of_i[U_set[i]])
                    else:
                        if(U_set[j] in W_set):
                            s_disc[(child_vari,i,j)] = - empty_score +math.log(disc_of_i[U_set[j]])
                        else:
                            s_disc[(child_vari,i,j)] = 0


                            
    return s_disc

def making_d_star(U_set,num_vari,child_vari):
    d_ij =[0 for i in range(num_vari)]
    for i in range(num_vari):
        if child_vari == i:
            d_ij[i] = (-1,)
        else:
            p_list = list()
            for j in range(len(U_set)):
                set1 = {k for k in U_set[j]}
                if(i in set1):
                    p_list.append(j)
            d_ij[i]=tuple(p_list)
    return d_ij



def making_QUBO(bayes):
    num_vari = bayes.num_variable
    max_parent_num = bayes.max_parent_set
    len_of_u =[0 for i in range(num_vari)]
    xi = [0 for i in range(num_vari)]
    delta = [0 for i in range(num_vari)]
    d_star =[list() for j in range(num_vari)]
    s_list = [dict() for i in range(num_vari)]
    for num in range(num_vari):
        bayes.make_score_disc_i(num)
        w_i=make_W_set_of_i(bayes.score_disc_list[num],num,num_vari,max_parent_num)
        u_i =decomposition_of_i(w_i,num,max_parent_num)
        bayes.u_set_list[num] = u_i
        #print(num,w_i)
        #print(num,u_i)
        s_list[num] = calc_s_of_i(u_i,w_i,bayes.score_disc_list[num],num)
        min_s = min(s_list[num].values())
        delta[num] = min_s
        xi[num] = -3 * min_s
        d_star[num] = making_d_star(u_i,num_vari,num)
        len_of_u[num] = len(u_i)
    qubo_size = int(sum(len_of_u) + num_vari + num_vari*(num_vari +1)/2)
    qubo = np.zeros((qubo_size,qubo_size))
    count_u = 0
    delta1 =  - min(delta)
    delta2 = (num_vari - 2) * delta1
    count_c=[0 for i in range(num_vari)]
    #print(s_list[0])
    #print(d_star)
    for i in range(num_vari):
        for k , v in s_list[i].items():
            a,b,c = k
            if b == c:
                qubo[b+count_u][c + count_u] += v
            else:
                qubo[b+count_u][c + count_u] += v + xi[i]
        count_c[i] = count_u
        count_u = count_u + len_of_u[i]
    count_b = 0
    for i in range(num_vari):
        qubo[count_u][count_u] += xi[i]
        for j in range(len_of_u[i]):
            qubo[count_u][j + count_b] += - xi[i]
        count_b = count_b + len_of_u[i]
        count_u = count_u + 1
    
    vari_list = [i for i in range(num_vari)]
    for conb in itertools.combinations(vari_list,3):
        i,j,k = conb
        r_ij = i*num_vari - int(i*(i+1)/2) + j - i -1
        r_jk = j*num_vari - int(j*(j+1)/2) + k - j -1
        r_ik = i*num_vari - int(i*(i+1)/2) + k - i -1
        qubo[count_u + r_ij][count_u + r_jk] += delta1
        qubo[count_u + r_ij][count_u + r_ik] += -delta1
        qubo[count_u + r_ik][count_u + r_jk] += -delta1
        qubo[count_u + r_ik][count_u + r_ik] += delta1
    for conb in itertools.combinations(vari_list,2):
        i,j = conb
        r_ij = i*num_vari - int(i*(i+1)/2) + j - i -1
        for d in d_star[i][j]:
            qubo[count_c[i]+d][count_u + r_ij] += delta2
        for d in d_star[j][i]:
            qubo[count_c[j]+d][count_u + r_ij] += -delta2
            qubo[count_c[i]+d][count_c[i]+d] += delta2
    
    bayes.count_len_of_u = len_of_u


    

    return qubo

def recreate(x,bayes):
    parent_set_list = [tuple() for i in bayes.num_variable]
    count_u = 0
    for i in bayes.num_variable:
        tcount = 0
        s1 = set()
    
        for j in range(bayes.count_len_of_u[i]):
            if x[count_u] == 1:
                s1 =s1 | {i for i in bayes.u_set_list[i][j]}
                tcount +=1
        s1 = s1 - {i}
        if tcount >2:
                print('error in x')
        parent_set_list[i] = tuple(np.sort([k for k in s1]))
    
    return parent_set_list


num_variable = 9
max_parent = 3
bayes = Bayes(num_variable,max_parent)

QUBO =making_QUBO(bayes)
np.savetxt('qubo.csv',QUBO)
