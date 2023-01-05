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
import sys
import threading
import time
#変更前の再帰関数の実行回数の上限を表示
print(sys.getrecursionlimit())

sys.setrecursionlimit(67108864) #64MB
threading.stack_size(2048*2048)  #2の20乗のstackを確保=メモリの確保

#変更後の再帰関数の実行回数の上限を表示
print(sys.getrecursionlimit())
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
    i_jdic = dict()

    def __init__(self,bif_file,num_parent,bdeu_arufa):
        reader = BIFReader(bif_file)
        self.variable_list = reader.get_variables()
        self.num_variable = len(self.variable_list)
        self.status_list = reader.get_states()
        self.data_set =BayesianModelSampling(reader.get_model()).forward_sample(size=int(1e4))
        self.bdeu_a = bdeu_arufa
        self.max_parent_set = num_parent 
        self.u_set_list = [list() for i in range(self.num_variable)]
        self.init_score_disc()
        self.max_model = reader.get_model()
        self.max_w_list = [0 for i in range(self.num_variable)]
        for i in range(self.num_variable):
            for j in range(self.num_variable):
                self.i_jdic[(i,j)] = list()


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
def calc_s_of_i(U_set,W_set,disc_of_i,child_vari,ij_dict):
    weight = 1000
    s_disc = dict()
    l = len(U_set)
    W_set_check_list = dict()
    
    for i in range(len(W_set)):
        W_set_check_list[W_set[i]] = 0
    
    empty_score = (disc_of_i[(child_vari,)])
    for i in range(l):
        for u in U_set[i]:
            ij_dict[(child_vari,u)].append((i,i))
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
                    W_set_check_list[pair_ij] = 1
                    #for u in pair_ij:
                    #    ij_dict[(child_vari,u)].append((0,i))
                else:
                    s_disc[(child_vari,i,j)] = 0
            else:
                if(pair_ij in W_set and W_set_check_list[pair_ij]==0):
                    W_set_check_list[pair_ij] = 1
                    #for u in pair_ij:
                    #    ij_dict[(child_vari,u)].append((i,j))
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
                            s_disc[(child_vari,i,j)] = (disc_of_i[U_set[i]]) +(disc_of_i[U_set[j]]) - 2 * empty_score + weight
                        else:
                            s_disc[(child_vari,i,j)] = - empty_score + (disc_of_i[U_set[i]]) + weight
                    else:
                        if(U_set[j] in W_set):
                            s_disc[(child_vari,i,j)] = - empty_score +(disc_of_i[U_set[j]]) + weight
                        else:
                            s_disc[(child_vari,i,j)] = 0 + weight


                            
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
    i_in_u_list = [list() for i in range(num_vari)]
    max_w_list = [0 for i in range(num_vari)]
    #debug
    #u_list = [list() for j in range(num_vari)]
    #debug ijou
    for num in range(num_vari):
        
        print(num,'calc A')
        print('start make score disc')
        bayes.make_score_disc_i(num)
        max_w_list[num] = max(bayes.score_disc_list[num].values())
        print('end make score disc')
        print('start make_W_set_of_i')
        w_i=make_W_set_of_i(bayes.score_disc_list[num],num,num_vari,max_parent_num)
        print('end W_set')
        #print(w_i) #debug
        print('start decomposition')
        u_i =decomposition_of_i(w_i,num,max_parent_num)
        i_in_u_list[num] = make_i_in_u_set(u_i,num_vari)
        print('end decomposition')
        #debug
        #u_list[num] = u_i
        #debug ijou
        bayes.u_set_list[num] = u_i
        
        #print(num,w_i)
        print(num,u_i)
        s_list[num] = calc_s_of_i(u_i,w_i,bayes.score_disc_list[num],num,bayes.i_jdic)
        
        min_s = min(s_list[num].values())
        delta[num] = min_s
        xi[num] = -3 * min_s+1000 #下限をあげる
        d_star[num] = making_d_star(u_i,num_vari,num)
        len_of_u[num] = len(u_i)
    #print(d_star)
    qubo_size = int(sum(len_of_u)  + num_vari*(num_vari +1)/2)
    bayes.qubo_size = qubo_size
    #print(qubo_size)
    qubo = np.zeros((qubo_size,qubo_size))
    count_u = 0
    delta1 =  - min(delta) +10
    delta2 = (num_vari - 2) * delta1 +1
    count_c=[0 for i in range(num_vari)]
    #print(s_list[0])
    #print(d_star)
    for i in range(num_vari):
        for k , v in s_list[i].items():
            a,b,c = k
            
            if b == 0 and c == 0:
                qubo[b+count_u][c + count_u] += v 
                print(v)
                #print((b+count_u,c + count_u),v,u_list[i][b],u_list[i][c])
            elif b == c:
                qubo[b+count_u][c + count_u] += v #+ xi[i]
            else:
                qubo[b+count_u][c + count_u] += v + xi[i]
                #print((b+count_u,c + count_u),v,u_list[i][b],u_list[i][c])
            #if b == 0 and c>0:
            #    qubo[b+count_u][c + count_u] +=xi[i]

        count_c[i] = count_u
        count_u = count_u + len_of_u[i]
    count_b = 0
    #print(qubo)
    for i in range(num_vari):
        print(i,'calc B')
        qubo[count_u][count_u] += xi[i]#4*xi[i]
        #print(count_u)
        for j in range(len_of_u[i]):
            qubo[j + count_b][count_u] += - xi[i]#- 4*xi[i]
            #print((count_u,j + count_b))
      
        count_b = count_b + len_of_u[i]
        count_u = count_u + 1
    ####################################
    #2cycle wo jogai suru
    for i in range(num_vari):
        u_set_of_i = bayes.u_set_list[i]
        for j in range(1,len(u_set_of_i)):
            for k in u_set_of_i[j]:
                    if i < k:
                        for l in i_in_u_list[k][i]:
                            qubo[count_c[i] + j][count_c[k] +l] = xi[i]
                
    

    ###################################
    '''
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
    '''
    bayes.count_len_of_u = len_of_u
    bayes.max_w_list = max_w_list
    
    

    return qubo

def recreate(x,bayes):
    parent_set_list = [[] for i in range(bayes.num_variable)]
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
        parent_set_list[i] = (np.sort([k for k in s1])).tolist()
        print('親変数',vari_list[i],'の','pの数は',tcount)
    print(parent_set_list)
    re_list =  [[] for i in range(bayes.num_variable)]
    for i in range(len(parent_set_list)):
        st_list = [vari_list[j] for j in parent_set_list[i]]
        re_list[i] = st_list
    if judge_dag(parent_set_list):
        print('dagである')
        eval_bay(bayes,re_list)
    else:
        print('dagではない')
    judge_2cycle(parent_set_list)
    return re_list,parent_set_list

def interpret_parent(parent_set_list,bayes):
    vari_lsit = bayes.variable_list
    #for i in range(len(parent_set_list)):
        


    return
def judge_dag(parent_set_list):
    leng = len(parent_set_list)
    topological = []
    void_q = []
    judge = True
    for u in range(leng):
        for i in range(leng):
            if len(parent_set_list[i] )== 0:
                void_q.append(i)
                parent_set_list[i].append(-1)
        if len(void_q) == 0:
            
            judge = False
            break
        x = void_q.pop(0)
        topological.append(x)
        for v in range(leng):
            if x in parent_set_list[v]:
                parent_set_list[v].remove(x)
        
    
    return judge

def judge_2cycle(parent_set_list):
    judge = True
    leng = len(parent_set_list)
    for i in range(leng):
        parent_set = parent_set_list[i]
        for j in range(len(parent_set)):
            if i in parent_set_list[parent_set[j]]:
                judge = False
    if judge:
        print('2cycle nashi')
    else:
        print('2cycle ari')
    return judge

def make_i_in_u_set(u_set,len_of_vari):
    i_in_u =[list() for i in range(len_of_vari)]
    for i in range(len(u_set)):
        for j in u_set[i]:
            i_in_u[j].append(i)
    return i_in_u
class preli_graph:
    node_status_list = list()
    node_time_stamp= list()
    child_set_list = list()
    parent_set_list = list()
    scc_tree = []
    time=0
    num_node = 0
    def __init__(self,num_nodes,parent_list):
        self.time = 0
        self.scc_tree = []
        self.parent_set_list = parent_list
        self.num_node = num_nodes
        self.node_status_list = [0 for i in range(num_nodes)]##0 -> white, 1 -> gray, 2 -> black
        self.node_time_stamp = [0 for i in range(num_nodes)]
        self.child_set_list = [[] for i in range(num_nodes)]
        for i in range(num_nodes):
            parent_set_of_i = self.parent_set_list[i]
            for j in range(len(parent_set_of_i)):
                self.child_set_list[parent_set_of_i[j]].append(i)
        self.scc()

    def scc(self):
        for i in range(self.num_node):
            if self.node_status_list[i] == 0:
                self.dfs(i)
        ##print(self.node_time_stamp)##
        time_list = np.array(self.node_time_stamp)
        time_list = np.argsort((-time_list))
        self.node_status_list = [0 for i in range(self.num_node)]
        ##print(time_list)
        for i in time_list:
            if self.node_status_list[i] == 0:
                self.dfs_t(i)
                tree_list =[i for i ,x in enumerate(self.node_status_list) if x == 2]
                for j in tree_list:
                    self.node_status_list[j] = 3
                self.scc_tree.append(tree_list)


        return
    def dfs_t(self,node):
        self.node_status_list[node] = 1
        for i in self.parent_set_list[node]:
            if self.node_status_list[i] == 0:
                self.dfs_t(i)
        self.node_status_list[node] = 2
        return
    def dfs(self,node):
        self.node_status_list[node] = 1
        for i in self.child_set_list[node]:
            if self.node_status_list[i] == 0:
                self.dfs(i)
        self.node_status_list[node] = 2
        self.time = self.time + 1
        self.node_time_stamp[node] = self.time
        
        return
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
        #print(bdeu.local_score(vari_list[i], parents=[]))
    #print(bdeu.score(BayesianNetwork()))
    calc_model=BayesianNetwork(calc_model)
    calc_bdeu = bdeu.score(calc_model)
    
    print('max_bdeu',max_bdeu)
    print('calc_bdeu',calc_bdeu)
    print('割合',calc_bdeu/max_bdeu)
    print('ここのスコアは高いほうがよい')
    return
max_parent = 2
bayes = Bayes('asia.bif',max_parent,1)



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

#with open('qubo_rec.csv','w') as f:
#    writer = csv.writer(f)
#    for i in range(bayes.qubo_size):
#         writer.writerow(QUBO[i][i::])
with open('qubo_rec.csv','w') as f:
    writer = csv.writer(f)
    for i in range(bayes.qubo_size):

        writer.writerow(QUBO[i])

x = [0 for i in range(bayes.qubo_size)]

#with open('minx.txt') as f:
    #for i in range(bayes.qubo_size):
        #x[i] = int(f.read(1))
print(kai_int)
#for i  in range(len(kai_int)):
    #print(kai_int[i],end='')

list,parent_set_list =recreate(kai_int,bayes)
#list2=recreate(kai_int,bayes)

print('')
kai_np = np.array(kai_int)
print('score',np.dot(np.dot((kai_np),QUBO ), (kai_np.T)))
for i in range(bayes.num_variable):
    
     print(bayes.variable_list[i],'の親変数は',list[i],'です')
#eval_bay(bayes,list)
#print(bayes.variable_list)

#for i in range(bayes.num_variable):
    #print(i,'の親変数は',list2[i],'です_sample')

#print(bayes.score_disc_list)
#print(bayes.variable_list)
#example_parent_set_list =[[1],[0],[0,7],[1,2],[6],[4],[5],[2,5,6]]
#graph = preli_graph(len(example_parent_set_list),example_parent_set_list)
#print(graph.scc_tree)

class Branch_and_Bound:
    
    uppor_bound =0
    chosen_list = []#eranda mono
    steps = -1#ima nandanme no branch ka
    branch_list = []#erabu list no list
    kai_list = [] #kai no list
    kai_parent_set = []
    weight = 1000000
    preli_kai = []
    kakikae = [] #kaikikaeta list
    seen_set_horizon = set()
    seen_set_horizon_no_steps = set()#miowatta
    run_set = set()#imamiteru
    Judge = True
    qubo_size = 0
    cycle_list = []
    def __init__(self,bayes,Qubo):
        self.QUBO = Qubo
        self.QUBO_origin =np.copy(Qubo)
        self.Bay = bayes
        client = FixstarsClient()
        client.token = "bGezdYx88VisgcbxBgmpP7i0ferFQiZG"
        client.parameters.timeout = 10  # タイムアウト1秒
        self.solver = Solver(client)
        self.count_sum_len_of_u = [0 for i in range(bayes.num_variable)]
        sum = 0
        for i in range(bayes.num_variable):
            self.count_sum_len_of_u[i] = sum
            sum += bayes.count_len_of_u[i]
        self.qubo_size = bayes.qubo_size
        print(self.count_sum_len_of_u)
        print('qubo_size',self.qubo_size)
        self.stlist = [0 for i in range(bayes.num_variable)]
        #debug
        #for i in range(bayes.count_len_of_u[i]):
        #    for j in range(bayes.count_len_of_u[i]):
        #        if i == 0 and j == 0:
        #            self.QUBO[i][j] = 0
        #        elif i <= j:
        #            self.QUBO[i][j] = self.weight + self.Bay.max_w_list[0] - self.Bay.score_disc_list[0][(0,)]

    

    def qubo_result(self,Qubo):
        qubo_amp = BinaryMatrix(Qubo)
        result = self.solver.solve(qubo_amp)
        gen = SymbolGenerator(BinaryPoly) 
        q = gen.array(self.Bay.qubo_size)
        kai = q.decode(result[0].values)
        kai_int =[int(i) for i in kai]
        return kai_int
    
    def make_parent_set_list(self,x):
        parent_set = [[] for i in range(self.Bay.num_variable)]
        count_u = 0
        for i in range(self.Bay.num_variable):
            tcount = 0
            s1 = set()
        
            for j in range(self.Bay.count_len_of_u[i]):
                if x[count_u] == 1:
                    
                    s1 =s1 | {i for i in self.Bay.u_set_list[i][j]}
                    tcount +=1
                count_u +=1
            s1 = s1 - {i}
        
            if tcount >2:
                print('error in x in make parent set_list',i)
                time.sleep(10)
            parent_set[i] = (np.sort([k for k in s1])).tolist()
        return parent_set
    
    def next_qubo(self):
        #print('QUBO',self.QUBO)
        print('nextqubo')
        print('step',self.steps)
        print('seen',self.seen_set_horizon)
        print('seen_no',self.seen_set_horizon_no_steps)
        print('branch',self.chosen_list[self.steps])
        print('root',self.chosen_list[self.steps])
        print('chosen',self.chosen_list)
        cont = True
        if len(self.branch_list[self.steps]) > self.chosen_list[self.steps] +1:
            self.chosen_list[self.steps] += 1
            a,b = self.kakikae.pop()
            for u in self.Bay.i_jdic[(a,b)]:
                    c,d = u
                    offset = self.count_sum_len_of_u[a]
                    self.QUBO[offset + c,offset + d]  = self.QUBO_origin[offset + c,offset + d] 
            #
            print('水平')
            kaki_a = [(g,h) for g,h in self.kakikae if g == a]
            for s in kaki_a:
                for u in self.Bay.i_jdic[s]:
                    a,b = s
                    c,d = u
                    offset = self.count_sum_len_of_u[a]
                    self.QUBO[offset + c,offset + d]  = self.weight
            #
            b_n = self.branch_list[self.steps][self.chosen_list[self.steps]]
            self.check_ban(b_n)#debug
            for u in self.Bay.i_jdic[b_n]:
                j,k = b_n
                c,d = u
                offset = self.count_sum_len_of_u[j]
                self.QUBO[offset + c,offset + d]   += self.weight + self.Bay.max_w_list[j] - self.Bay.score_disc_list[j][(j,)]
            #

            
            self.kakikae.append(b_n)
            if (b_n in self.seen_set_horizon_no_steps):
                print('omit')
                self.next_qubo()
            #
            #
        else:
            print('繰り上がり')
            if self.steps == 0 and len(self.branch_list[self.steps]) == self.chosen_list[self.steps] +1:
                cont = False
            else:
                a,b = self.kakikae.pop()
                for u in self.Bay.i_jdic[(a,b)]:
                    c,d = u
                    offset = self.count_sum_len_of_u[a]
                    self.QUBO[offset + c,offset + d]  = self.QUBO_origin[offset + c,offset + d] 
                
                #
                kaki_a = [(g,h) for g,h in self.kakikae if g == a]
                for s in kaki_a:
                    for u in self.Bay.i_jdic[s]:
                        a,b = s
                        c,d = u
                        offset = self.count_sum_len_of_u[a]
                        self.QUBO[offset + c,offset + d]  = self.weight
                #
                self.chosen_list.pop()
                self.branch_list.pop()
                #self.kai_list.pop()
                #
                e,f= self.kakikae[-1]
                for u in self.Bay.i_jdic[(e,f)]:
                    c,d = u
                    offset = self.count_sum_len_of_u[e]
                    self.QUBO[offset + c,offset + d]  = self.QUBO_origin[offset + c,offset + d] 
                #
                kaki_b = [(g,h) for g,h in self.kakikae if g == e]
                for s in kaki_b:
                    for u in self.Bay.i_jdic[s]:
                        a,b = s
                        c,d = u
                        offset = self.count_sum_len_of_u[a]
                        self.QUBO[offset + c,offset + d]  = self.weight
                #
                self.seen_set_horizon = {(y,z,w) for (y,z,w) in self.seen_set_horizon if y < self.steps}
                self.seen_set_horizon_no_steps = {(z,w) for (y,z,w) in self.seen_set_horizon if y < self.steps}

                self.seen_set_horizon.add((self.steps -1,e,f))
                self.seen_set_horizon_no_steps.add((e,f))
                #
                print('seen',self.seen_set_horizon)
                print('seen_no',self.seen_set_horizon_no_steps)
                
                self.steps = self.steps -1
                self.next_qubo()
                
        print('cont',cont)


        return cont
    def dag_search(self,parents):
        
        leng = len(parents)
        topological = []
        void_q = []
        self.Judge = True
        for u in range(leng):
            for i in range(leng):
                if len(parents[i] )== 0:
                    void_q.append(i)
                    parents[i].append(-1)
            if len(void_q) == 0:
                
                self.Judge = False
                break
            x = void_q.pop(0)
            topological.append(x)
            for v in range(leng):
                if x in parents[v]:
                    parents[v].remove(x)
        
    
        return         
    def i_child(self,num,parent,min_cycle):
        child_set = []
        for i in range(len(parent)):
            if num in parent[i] and i in min_cycle:
                child_set.append(i)
        
        return child_set
            
    def search_cycle(self,min_cycle,parent):
        self.cycle_list.clear()
        #u =np.random.randint(0,len(min_cycle))
        self.stlist = [0 for i in range(bayes.num_variable)]
        for i in range(len(self.stlist)):
            if i in min_cycle:
                self.stlist[i] = 0
            else:
                self.stlist[i] = -1

        start = min_cycle[0]
        self.cycle_list = [start]
        self.df_search(start,parent,start)
        return self.cycle_list
    def df_search(self,vari,parent,start):
        self.stlist[vari] = 1
        con = True
        for i in parent[vari]:
            if i == start:
                con = False
                break

            if self.stlist[i] == 0 and con:
                self.cycle_list.append(i)
                con = con and self.df_search(i,parent,start)
        if con:
            a = self.cycle_list.pop()
            print('pop',vari,a)

        return con
    def depthsf(self):
        print('depthsf')
        print('step',self.steps)
        print('seen',self.seen_set_horizon)
        print('seen_no',self.seen_set_horizon_no_steps)
        if(self.steps >-1):
            print('branch',self.chosen_list[self.steps])
            print('root',self.chosen_list[0])
            print('chosen',self.chosen_list)

        qubo_result =self.qubo_result(self.QUBO)
        parent_set = self.make_parent_set_list(qubo_result)
        #
        print('kai_parent',self.kai_parent_set)
        #
        qubo_result_np = np.array(qubo_result)
        print('qubo_result',qubo_result_np)
        pre_score = np.dot(np.dot((qubo_result_np),self.QUBO), (qubo_result_np.T))
        
        print('pre_score',pre_score)
        min_cycle = []
        min_cycle_num = len(parent_set) +1
        self.dag_search(parent_set)
        parent_set = self.make_parent_set_list(qubo_result)
        print('parent_set',parent_set)
        print('judge_dag',self.Judge)
        if pre_score >= self.uppor_bound:
            if len(self.kakikae) > 0:
                e,f= self.kakikae[-1]
                self.seen_set_horizon.add((self.steps,e,f))
                self.seen_set_horizon_no_steps.add((e,f))
            if self.next_qubo():
                print('prescore>= ')
                self.depthsf()
        else:
            if self.Judge:
                print('judgeee')
                print('parent',parent_set)
                self.uppor_bound = pre_score
                self.preli_kai = qubo_result
                self.kai_parent_set = parent_set
                if len(self.kakikae) > 0:
                    e,f= self.kakikae[-1]
                    self.seen_set_horizon.add((self.steps,e,f))
                    self.seen_set_horizon_no_steps.add((e,f))
                if self.steps >= 0 and self.next_qubo() :
                    print('self judge')
                    self.depthsf()
                
            else:
                print('notttjudge')
                print('parent_set',parent_set)
                self.steps = self.steps + 1
                
                pre_graph = preli_graph(len(parent_set),parent_set)
                woods = pre_graph.scc_tree
                for tree in woods:
                    if min_cycle_num > len(tree) >1:
                        min_cycle = tree
                #
                print('wood',woods)
                #

                #blist = self.i_child(min_cycle[0],parent_set,min_cycle)
                #ban_c = (self.i_child(min_cycle[0],parent_set,min_cycle))[0]

                ban = min_cycle[0]
                #ban_p = 0

                #for i in parent_set[ban]:
                #    if i in min_cycle:
                #        ban_p = i
                cycle_l = self.search_cycle(min_cycle,parent_set)
                blist = []
                for i in range(len(cycle_l)):
                    v = (i +1) % len(cycle_l)
                    blist.append((cycle_l[i],cycle_l[v]))
                #blist = #[(ban_c,ban),(ban,ban_p),(ban_p,ban_c)]#(chi,pa)

                

                m = 0
                
                loop = True
                conti = True
                blist_sin = []
                for b in blist:
                    if (b not in self.seen_set_horizon_no_steps):
                        blist_sin.append(b)
                #
                #if self.steps>18:
                #    blist_sin = []
                #
                if len(blist_sin) == 0:
                    print('omit in depth')
                    self.steps += -1
                    self.next_qubo()
                    conti = False
                    if self.steps == 0 and len(self.branch_list[self.steps]) == self.chosen_list[self.steps] +1:
                        print('終わり')
                    else:
                        self.depthsf()
                blist = blist_sin

#                while loop:
#
#                    b = blist[m]
#                    
#                    if (b in self.seen_set_horizon_no_steps):
#                        m += 1
##                        if m == len(blist):
#
#                            print('omit in depth')
#                            self.steps = self.steps -1
#                            self.next_qubo()
#                            if self.steps == 0 and len(self.branch_list[self.steps]) == self.chosen_list[self.steps] +1:
#                            
#                                print('終わり')
#                                conti = False
#                                m = m -1
#                                loop = False
#                            else:
#                                self.depthsf()
#                                loop = False
#                                conti = False
#
#                    else:
#                        loop = False
                #
                if conti:

                    #ban_domain = self.count_sum_len_of_u[ban]
                    #self.kai_list.append(qubo_result)
                    #kai_ban = qubo_result
                    #kai_1 = [j for j ,x in enumerate(kai_ban) if x == 1 and ban_domain <= j < ban_domain + self.Bay.count_len_of_u[ban]]
                    #kai_1 = np.array(kai_1)
                    #kai_1 = np.sort(kai_1)
                    #u = self.Bay.u_set_list[ban]
                #set1 = {k for k in u[kai_1[0]- self.count_sum_len_of_u[ban]]}
                #set2 = {k for k in u[kai_1[1] - self.count_sum_len_of_u[ban]]}
                #set12 = (set1 | set2) - {ban}
                #ban_par = np.sort([k for k in set12])
                
                #self.branch_list.append(ban_par)
                #debug
                    #ban_par = blist
                    self.branch_list.append(blist)
                    self.chosen_list.append(m)
                #ban = min_cycle[0]
                #ban_domain = self.count_sum_len_of_u[ban]
                    g,h = blist[m]
                    ban_domain = self.count_sum_len_of_u[g]
                
                #kai_ban = self.kai_list[self.steps]
                    kai_1 = [j for j ,x in enumerate(qubo_result) if x == 1 and ban_domain <= j < ban_domain + self.Bay.count_len_of_u[ban]]
                #kai_1 = np.array(kai_1)
                #kai_1 = np.sort(kai_1)
            
                
                    #if len(kai_1) >2:
                    #    print('error occur')
                    print('kai_1',kai_1)
                


                    self.kakikae.append(blist[m])
                    print('ban, banpar',blist,'bm',blist[m])
                    print('seen',self.seen_set_horizon)
                    print('seen_no',self.seen_set_horizon_no_steps)
                #print(self.QUBO[0][0])
                #print(self.QUBO[kai_1[0]][kai_1[0]])
                #print(self.QUBO[kai_1[1]][kai_1[1]])
                    self.check_ban(blist[m])#debug
                    for u in self.Bay.i_jdic[blist[m]]:
                    #print('kakikae')
                        a,b = u
                        parent,child = blist[m]
                        offset = self.count_sum_len_of_u[parent]
                        self.QUBO[offset + a][offset + b]  =self.weight#+= self.weight + self.Bay.max_w_list[parent] - self.Bay.score_disc_list[parent][(parent,)]
                        print((offset+a,offset+b))

                #print('henka',self.Bay.i_jdic[(ban_par[0],ban)])
                #print('qubo,qubo_origin',(self.QUBO == self.QUBO_origin).all())
                    print('branch root',self.branch_list[0])
                    #debug
                    

                    #
                    self.depthsf()


        return
    def check_ban(self,ban):
        l = [self.branch_list[i][self.chosen_list[i]] for i in range(len(self.chosen_list))]
        if len(l) != len(set(l)):
            print('error occur in check_ban')
            print(l)
            print(self.count_sum_len_of_u)
            g,h = ban
            print('ij', self.Bay.i_jdic[ban])
            for u in self.Bay.i_jdic[ban]:
                v,w = u
                print((self.count_sum_len_of_u[g]+v,self.count_sum_len_of_u[g]+w))
            np.savetxt('qubo_error.csv',QUBO,delimiter=',' ,fmt="%.i")
            time.sleep(10)
        return
    def qubo_solve(self):
        for i in range(len(self.count_sum_len_of_u)):
            if i == len(self.count_sum_len_of_u):
                k =self.qubo_size
            else:
                k = self.count_sum_len_of_u[i+1]
            for j in range(self.count_sum_len_of_u[i]+1,k):
                self.QUBO[j,j] = self.weight
            self.depthsf()
            self.QUBO = self.QUBO_origin
#    def qubo_solve(self):
#        for i in range(len(self.count_sum_len_of_u)):
#            for l in range(len(self.count_sum_len_of_u)):
#                print('qubo solve',i,l)
#                offset = self.count_sum_len_of_u[i]
#                for j in range(self.count_sum_len_of_u[i]):
#                    for k in range(j+1,self.count_sum_len_of_u[i]):
#                        self.QUBO[offset+j,offset+k] += self.weight + self.Bay.max_w_list[i] - self.Bay.score_disc_list[i][(i,)]
#                for u in range(len(self.count_sum_len_of_u)):
#                    for q in self.Bay.i_jdic[(u,l)]:
#                        a,b = q
#                        offset = self.count_sum_len_of_u[u]
#                        self.QUBO[offset+a,offset+b] += self.weight + self.Bay.max_w_list[u] - self.Bay.score_disc_list[u][(u,)]
#                self.depthsf()
##               self.QUBO = self.QUBO_origin
#
#
#
#        return

    #def qubo_solve(self):
        #for i in range(len(self.count_sum_len_of_u)):
        # self.QUBO = self.QUBO[self.count_sum_len_of_u[1]:,self.count_sum_len_of_u[1]:]
         #self.count_sum_len_of_u

            
        return
print('ij',bayes.i_jdic)
time_sta = time.time()
bab = Branch_and_Bound(bayes,QUBO)
#bab.qubo_solve()
bab.depthsf()
time_end = time.time()
proc_time = time_end- time_sta
print('解',bab.kai_parent_set)
print('実行時間',proc_time)
list,parent_set_list = recreate(bab.preli_kai,bayes)
print('')

for i in range(bayes.num_variable):
    
     print(bayes.variable_list[i],'の親変数は',list[i],'です')