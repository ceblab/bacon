from calendar import c
from re import L, M
import numpy as np
import itertools
import copy
import math
import csv
import pprint
import pandas as pd
from amplify import BinaryMatrix
from amplify import BinarySymbolGenerator
from amplify import (
    BinaryPoly,
    BinaryQuadraticModel,
    Solver,
    SymbolGenerator,
)
from amplify.client import FixstarsClient
from pgmpy.readwrite import BIFReader
from pgmpy.sampling import BayesianModelSampling
from scipy.special import gammaln
np.set_printoptions(threshold=np.inf)
#from regex import W
from pgmpy.estimators import BDeuScore
np.random.seed(314)
from pgmpy.models import BayesianNetwork
class Bayes:
    variable_list = list()
    status_list = dict()
    data_set = pd.DataFrame()
    bdeu_a = 1

    num_variable = 0
    max_parent_set = 0
    score_disc_list = list()
    u_set_list = list()
    count_len_of_u = list()
    qubo_size = 0
    max_model = BayesianNetwork()


    def __init__(self,bif_file,num_parent,bdeu_arufa):
        reader = BIFReader(bif_file)
        self.variable_list = reader.get_variables()
        self.num_variable = len(self.variable_list)
        self.status_list = reader.get_states()
        self.data_set =BayesianModelSampling(reader.get_model()).forward_sample(size=int(1e3))
        self.bdeu_a = bdeu_arufa
        self.max_parent_set = num_parent 
        self.u_set_list = [list() for i in range(self.num_variable)]
        self.init_score_disc()
        self.max_model = reader.get_model()
        self.kani_qubo_size = 0

    def init_score_disc(self):
        self.score_disc_list = [dict() for i in range(self.num_variable)]

        
    def output(self):
        print(self.num_variable)

    def score_kara(self,child_vari):
        colum_name = self.variable_list[child_vari]
        child_status_list = self.status_list.get(colum_name)
        ganma = len(child_status_list)
        a_ijk = self.bdeu_a / ganma
        n_list = [0 for i in range(ganma)]
        a_list = [a_ijk for i in range(ganma)]
        s = gammaln(self.bdeu_a) - gammaln(len(self.data_set)+self.bdeu_a)

        for i in range(ganma):
            N_ijk = (self.data_set[colum_name] == child_status_list[i]).sum()
            n_list[i] = N_ijk + a_ijk
        s += (np.sum(gammaln(n_list)) - np.sum(gammaln(a_list)))
        
        return s
    def score_1(self,subset,child_vari):
        a = subset[0]
        parent_colum_name = self.variable_list[a]
        child_colum_name =self.variable_list[child_vari]
        parent_status_list = self.status_list.get(parent_colum_name)
        child_status_list = self.status_list.get(child_colum_name)
        ganma = len(child_status_list)
        beta = len(parent_status_list)
        a_ijk = self.bdeu_a / (ganma*beta)
        a_ij = self.bdeu_a / beta
        n_list = [0 for i in range(ganma)]
        a_list = [a_ijk for i in range(ganma)]
        s = 0
        for j in range(beta):

            d_set_j = (self.data_set[parent_colum_name] == parent_status_list[j])
            N_ij = d_set_j.sum()
            
            s += (gammaln(a_ij) -gammaln(N_ij + a_ij)) 
            d_set_sub = self.data_set[d_set_j]
            for k in range(ganma):

                N_ijk = (d_set_sub[child_colum_name] == child_status_list[k]).sum()
                
                n_list[k] = N_ijk + a_ijk
            s += (np.sum(gammaln(n_list)) - np.sum(gammaln(a_list)))
        return s

    def score_2(self,subset,child_vari):
        a,b= subset
        parent_colum_name_1 = self.variable_list[a]
        parent_colum_name_2 = self.variable_list[b]
        child_colum_name =self.variable_list[child_vari]
        parent_status_list_1 = self.status_list.get(parent_colum_name_1)
        parent_status_list_2 = self.status_list.get(parent_colum_name_2)
        child_status_list = self.status_list.get(child_colum_name)
        ganma = len(child_status_list)
        beta = len(parent_status_list_1) * len(parent_status_list_2)
        a_ijk = self.bdeu_a / (ganma*beta)
        a_ij = self.bdeu_a / beta
        n_list = [0 for i in range(ganma)]
        a_list = [a_ijk for i in range(ganma)]
        s = 0
        for n in range(len(parent_status_list_1)):
            for m in range(len(parent_status_list_2)):



                d_set_j = (self.data_set[parent_colum_name_1] == parent_status_list_1[n]) & (self.data_set[parent_colum_name_2] == parent_status_list_2[m])
                N_ij = d_set_j.sum()
                s += (gammaln(a_ij) -gammaln(N_ij + a_ij) )
                d_set_sub = self.data_set[d_set_j]
                for k in range(ganma):

                    N_ijk = (d_set_sub[child_colum_name] == child_status_list[k]).sum()
                    n_list[k] = N_ijk + a_ijk
                s += (np.sum(gammaln(n_list)) - np.sum(gammaln(a_list)))
        return s

    def score_3(self,subset,child_vari):
        a,b,c= subset
        parent_colum_name_1 = self.variable_list[a]
        parent_colum_name_2 = self.variable_list[b]
        parent_colum_name_3 = self.variable_list[c]
        child_colum_name =self.variable_list[child_vari]
        parent_status_list_1 = self.status_list.get(parent_colum_name_1)
        parent_status_list_2 = self.status_list.get(parent_colum_name_2)
        parent_status_list_3 = self.status_list.get(parent_colum_name_3)
        child_status_list = self.status_list.get(child_colum_name)
        ganma = len(child_status_list)
        beta = len(parent_status_list_1) * len(parent_status_list_2) * len(parent_status_list_3)
        a_ijk = self.bdeu_a / (ganma*beta)
        a_ij = self.bdeu_a / beta
        n_list = [0 for i in range(ganma)]
        a_list = [a_ijk for i in range(ganma)]
        s = 0
        for n in range(len(parent_status_list_1)):
            for m in range(len(parent_status_list_2)):
                for l in range(len(parent_status_list_3)):
            
                    d_set_j = (self.data_set[parent_colum_name_1] == parent_status_list_1[n]) & (self.data_set[parent_colum_name_2] == parent_status_list_2[m]) &(self.data_set[parent_colum_name_3] == parent_status_list_3[l])
                    N_ij = d_set_j.sum()
                    s += (gammaln(a_ij) -gammaln(N_ij + a_ij)) 
                    d_set_sub = self.data_set[d_set_j]
                    for k in range(ganma):

                        N_ijk = (d_set_sub[child_colum_name] == child_status_list[k]).sum()
                        n_list[k] = N_ijk + a_ijk
                    s += (np.sum(gammaln(n_list)) - np.sum(gammaln(a_list)))
        return s


    def score(self,subset,child_vari):
        p_set_len = len(subset)
        s = 0
        if p_set_len == 1:
            s = self.score_1(subset,child_vari)
        elif p_set_len == 2:
            s = self.score_2(subset,child_vari)
        elif p_set_len == 3:
            s = self.score_3(subset,child_vari)

        return s

    def make_score_disc_i(self,child_vari):
        #self.score_disc_list[child_vari][(child_vari,)] = self.score_kara(child_vari)
        #vari_list = [i for i in range(self.num_variable) if i != child_vari]

        #for m in range(1, self.max_parent_set + 1):
        #    for conb in itertools.combinations(vari_list,m ):
        #        
        #        self.score_disc_list[child_vari][conb] = self.score(conb,child_vari)
        ##debug
        bdeu = BDeuScore(self.data_set, equivalent_sample_size=5)
        self.score_disc_list[child_vari][(child_vari,)] = bdeu.local_score(self.variable_list[child_vari], parents=[])
        vari_list = [i for i in range(self.num_variable) if i != child_vari]

        for m in range(1, self.max_parent_set + 1):
            for conb in itertools.combinations(vari_list,m ):
                parent_list = [self.variable_list[i] for i in conb]
                self.score_disc_list[child_vari][conb] = bdeu.local_score(self.variable_list[child_vari], parents=parent_list)


def make_W_set_of_i(disc_of_i,child_vari,num_variable,max_parent_set):
    if max_parent==3:
        wi = make_W_set_of_i_three(disc_of_i,child_vari,num_variable,max_parent_set)
    elif max_parent ==2:
        wi = make_W_set_of_i_two(disc_of_i,child_vari,num_variable,max_parent_set)
    return wi

def make_W_set_of_i_two(disc_of_i,child_vari,num_variable,max_parent_set):
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
    return W_set_of_i

def make_W_set_of_i_three(disc_of_i,child_vari,num_variable,max_parent_set):
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
    empty_score = (disc_of_i[(child_vari,)])
    for i in range(l):
        for j in range(i,l):
            if i ==0 and j == 0:
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
                    s_disc[(child_vari,i,j)] = -(disc_of_i[pair_ij]) + empty_score
                else:
                    s_disc[(child_vari,i,j)] = 0
            else:
                if(pair_ij in W_set):
                    if(U_set[i] in W_set):
                        if(U_set[j] in W_set):
                            s_disc[(child_vari,i,j)] = - (disc_of_i[pair_ij]) + (disc_of_i[U_set[i]]) +(disc_of_i[U_set[j]]) -empty_score
                        else:
                            s_disc[(child_vari,i,j)] = - (disc_of_i[pair_ij]) + (disc_of_i[U_set[i]])
                    else:
                        if(U_set[j] in W_set):
                            s_disc[(child_vari,i,j)] = - (disc_of_i[pair_ij]) +(disc_of_i[U_set[j]])
                        else:
                            s_disc[(child_vari,i,j)] = - (disc_of_i[pair_ij]) + empty_score
                else:
                    if(U_set[i] in W_set):
                        if(U_set[j] in W_set):
                            s_disc[(child_vari,i,j)] = (disc_of_i[U_set[i]]) +(disc_of_i[U_set[j]]) - 2 * empty_score
                        else:
                            s_disc[(child_vari,i,j)] = - empty_score + (disc_of_i[U_set[i]])
                    else:
                        if(U_set[j] in W_set):
                            s_disc[(child_vari,i,j)] = - empty_score +(disc_of_i[U_set[j]])
                        else:
                            s_disc[(child_vari,i,j)] = 0


                            
    return s_disc
'''
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
'''
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
    #debug
    #u_list = [list() for j in range(num_vari)]
    #debug ijou
    for num in range(num_vari):
        print(num,'calc A')
        bayes.make_score_disc_i(num)
        
        w_i=make_W_set_of_i(bayes.score_disc_list[num],num,num_vari,max_parent_num)
        
        #print(w_i) #debug
        u_i =decomposition_of_i(w_i,num,max_parent_num)
        #debug
        #u_list[num] = u_i
        #debug ijou
        bayes.u_set_list[num] = u_i
        
        #print(num,w_i)
        print(num,u_i)
        s_list[num] = calc_s_of_i(u_i,w_i,bayes.score_disc_list[num],num)
        
        min_s = min(s_list[num].values())
        delta[num] = min_s
        xi[num] = -3 * min_s+10 #下限をあげる
        d_star[num] = making_d_star(u_i,num_vari,num)
        len_of_u[num] = len(u_i)
    #print(d_star)
    qubo_size = int(sum(len_of_u)  + num_vari*(num_vari +1)/2)
    bayes.qubo_size = qubo_size
    bayes.kani_qubo_size = int(sum(len_of_u)  + num_vari)
    #print(qubo_size)
    qubo = np.zeros((qubo_size,qubo_size))
    count_u = 0
    delta1 =  - min(delta) + 10
    delta2 = (num_vari - 2) * delta1 +1
    count_c=[0 for i in range(num_vari)]
    #print(s_list[0])
    #print(d_star)
    for i in range(num_vari):
        for k , v in s_list[i].items():
            a,b,c = k
            
            if b == c:
                qubo[b+count_u][c + count_u] += v
                #print((b+count_u,c + count_u),v,u_list[i][b],u_list[i][c])
            else:
                qubo[b+count_u][c + count_u] += v + xi[i]
                #print((b+count_u,c + count_u),v,u_list[i][b],u_list[i][c])
        count_c[i] = count_u
        count_u = count_u + len_of_u[i]
    count_b = 0
    #print(qubo)
    for i in range(num_vari):
        print(i,'calc B')
        qubo[count_u][count_u] += xi[i]
        #print(count_u)
        for j in range(len_of_u[i]):
            #qubo[count_u][j + count_b] += - xi[i]
            qubo[j + count_b][count_u] += - xi[i]
            #print((count_u,j + count_b))
      
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
    #print(22222)
    #print(qubo)
    #print(count_c)
    for conb in itertools.combinations(vari_list,2):
        i,j = conb
        r_ij = i*num_vari - int(i*(i+1)/2) + j - i -1
        for d in d_star[i][j]:
            #print((count_c[i]+d,count_u + r_ij))
            qubo[count_c[i]+d][count_u + r_ij] += delta2
        for d in d_star[j][i]:
            qubo[count_c[j]+d][count_u + r_ij] += -delta2
            qubo[count_c[j]+d][count_c[j]+d] += delta2
    #print(33333)
    #print(qubo)
    bayes.count_len_of_u = len_of_u

    
    

    return qubo

def recreate(x,bayes):
    parent_set_list = [list() for i in range(bayes.num_variable)]
    vari_list = bayes.variable_list
    count_u = 0
    for i in range(bayes.num_variable):
        tcount = 0
        s1 = set()
        
        for j in range(bayes.count_len_of_u[i]):
            if x[count_u] == 1:
                
                s1 =s1 | {i for i in bayes.u_set_list[i][j]}
                tcount +=1
            count_u +=1
        s1 = s1 - {i}
        
        if tcount >2:
                print('error in x',i)
        parent_set_list[i] = list(np.sort([k for k in s1]))
    print(parent_set_list)
    re_list =  [list() for i in range(bayes.num_variable)]
    for i in range(len(parent_set_list)):
        st_list = [vari_list[j] for j in parent_set_list[i]]
        re_list[i] = st_list
    if judge_dag(parent_set_list):
        print('dagである')
    else:
        print('dagではない')

    return re_list

def interpret_parent(parent_set_list,bayes):
    vari_lsit = bayes.variable_list
    #for i in range(len(parent_set_list)):
        


    return
def judge_dag(parent_set_list):
    leng = len(parent_set_list)
    topological = list()
    void_q = list()
    judge = True
    for u in range(leng):
        for i in range(leng):
            if len(parent_set_list[i] )== 0:
                void_q.append(i)
                parent_set_list[i].append(-1)
        if len(void_q) == 0:
            print('dagではない')
            judge = False
            break
        x = void_q.pop(0)
        topological.append(x)
        for v in range(leng):
            if x in parent_set_list[v]:
                parent_set_list[v].remove(x)
        
    
    return judge

def eval_bay(bayes,re_list):
    vari_list = bayes.variable_list
    bdeu = BDeuScore(bayes.data_set, equivalent_sample_size=5)
    #model_max =BayesianNetwork([('asia','tub'),('tub','either'),('smoke','lung'),('smoke','bronc'),('either','xray'),('either','dysp'),('bronc','dysp')])
    model_max = bayes.max_model
    max_bdeu = bdeu.score(model_max)
    calc_model = []
    a = []
    for i in range(len(vari_list)):
        a = [(vari_list[i],re_list[i][j]) for j in range(len(re_list[i]))]
        calc_model.extend(a)
        a = []
    calc_model=BayesianNetwork(calc_model)
    calc_bdeu = bdeu.score(calc_model)
    
    print('max_bdeu',max_bdeu)
    print('calc_bdeu',calc_bdeu)
    print('割合',calc_bdeu/max_bdeu)
    print('ここのスコアは高いほうがよい')
    return
max_parent = 2
bayes = Bayes('alarm.bif',max_parent,1)



QUBO =making_QUBO(bayes)

#print(bayes.score_disc_list)

qubo_amp = BinaryMatrix(QUBO)

# イジングマシンクライアントの設定
client = FixstarsClient()
client.token = "bGezdYx88VisgcbxBgmpP7i0ferFQiZG"
client.parameters.timeout = 1000  # タイムアウト1秒

# ソルバーの実行
solver = Solver(client)
result = solver.solve(qubo_amp)
gen = SymbolGenerator(BinaryPoly) 
q = gen.array(bayes.qubo_size)
# 結果の解析
kai = q.decode(result[0].values)
for solution in result:
    print(f"q = {q.decode(solution.values)}")

kai_int =[int(i) for i in kai]
np.savetxt('qubo.csv',QUBO,delimiter=',' ,fmt="%.i")

with open('qubo_sample.csv','w') as f:
    writer = csv.writer(f)
    for i in range(bayes.qubo_size):
        writer.writerow(QUBO[i])
        #writer.writerow(QUBO[i][i::])


x = [0 for i in range(bayes.qubo_size)]

#with open('minx.txt') as f:
    #for i in range(bayes.qubo_size):
        #x[i] = int(f.read(1))
print(kai_int)
#for i  in range(len(kai_int)):
    #print(kai_int[i],end='')

list =recreate(kai_int,bayes)
#list2=recreate(kai_int,bayes)

print('')
kai_np = np.array(kai_int)
print('score',np.dot(np.dot((kai_np),QUBO ), (kai_np.T)))
for i in range(bayes.num_variable):
    
     print(bayes.variable_list[i],'の親変数は',list[i],'です')
eval_bay(bayes,list)
#print(bayes.variable_list)

#for i in range(bayes.num_variable):
    #print(i,'の親変数は',list2[i],'です_sample')

#print(bayes.score_disc_list)
print(bayes.variable_list)
print('qubo_size',bayes.qubo_size)
print('kani_qubo_size',bayes.kani_qubo_size)