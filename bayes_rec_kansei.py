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
    weight = 1000
    s_disc = dict()
    l = len(U_set)
    W_set_check_list = dict()
    for i in range(len(W_set)):
        W_set_check_list[W_set[i]] = 0
    
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
                    W_set_check_list[pair_ij] = 1
                else:
                    s_disc[(child_vari,i,j)] = 0
            else:
                if(pair_ij in W_set and W_set_check_list[pair_ij]==0):
                    W_set_check_list[pair_ij] == 1
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
        s_list[num] = calc_s_of_i(u_i,w_i,bayes.score_disc_list[num],num)
        
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
                qubo[b+count_u][c + count_u] += v + xi[i]
            else:
                qubo[b+count_u][c + count_u] += v + 2*xi[i]
                #print((b+count_u,c + count_u),v,u_list[i][b],u_list[i][c])
            if b == 0 and c>0:
                qubo[b+count_u][c + count_u] += xi[i]

        count_c[i] = count_u
        count_u = count_u + len_of_u[i]
    count_b = 0
    #print(qubo)
    for i in range(num_vari):
        print(i,'calc B')
        qubo[count_u][count_u] += 4*xi[i]
        #print(count_u)
        for j in range(len_of_u[i]):
            qubo[j + count_b][count_u] += - 4*xi[i]#qubo[j + count_b][count_u] += - xi[i]
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
    parent_set_list = [list() for i in range(bayes.num_variable)]
    
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
    
    if judge_dag(parent_set_list):
        print('dagである')
        
    else:
        print('dagではない')
    judge_2cycle(parent_set_list)
    return parent_set_list

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
def making_QUBO_C(bayes):
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
max_parent = 3
bayes = Bayes(4,3)

QUBO_C = making_QUBO_C(bayes)

QUBO =making_QUBO(bayes)

#print(bayes.score_disc_list)
qubo_amp = BinaryMatrix(QUBO_C)

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

#print(bayes.variable_list)

#for i in range(bayes.num_variable):
    #print(i,'の親変数は',list2[i],'です_sample')

#print(bayes.score_disc_list)



class Branch_and_Bound:
    
    uppor_bound =0
    chosen_list = []#eranda mono
    steps = -1#ima nandanme no branch ka
    branch_list = []#erabu list no list
    kai_list = [] #kai no list
    kai_parent_set = []
    weight = 1000
    preli_kai = []
    kaikikae = [] #kaikikaeta list
    seen_set_horizon = set()
    seen_set_horizon_no_steps = set()#miowatta
    run_set = set()#imamiteru
    
    qubo_size = 0
    def __init__(self,bayes,Qubo):
        self.QUBO = Qubo
        self.Bay = bayes
        client = FixstarsClient()
        client.token = "bGezdYx88VisgcbxBgmpP7i0ferFQiZG"
        client.parameters.timeout = 1000  # タイムアウト1秒
        self.solver = Solver(client)
        self.count_sum_len_of_u = [0 for i in range(bayes.num_variable)]
        sum = 0
        for i in range(bayes.num_variable):
            self.count_sum_len_of_u[i] = sum
            sum += bayes.count_len_of_u[i]
        self.qubo_size = bayes.qubo_size

        print(self.qubo_size)


        

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
                print('error in x',i)
            parent_set[i] = (np.sort([k for k in s1])).tolist()
        return parent_set
    
    def next_qubo(self):
        print('step',self.steps)
        print('seen',self.seen_set_horizon)
        print('seen_no',self.seen_set_horizon_no_steps)
        cont = True
        if len(self.branch_list[self.steps]) > self.chosen_list[self.steps] +1:
            self.chosen_list[self.steps] += 1
            a,b,value = self.kaikikae.pop()
            self.QUBO[a][b] = value
            #
            
            #
            ban = self.branch_list[self.steps][self.chosen_list[self.steps]]
            ban_domain = self.count_sum_len_of_u[ban]
            kai_ban = self.kai_list[self.steps]
            kai_1 = [j for j ,x in enumerate(kai_ban) if x == 1 and ban_domain <= j < ban_domain + self.Bay.count_len_of_u[ban]]
            kai_1 = np.array(kai_1)
            kai_1 = np.sort(kai_1)
            #

            if len(kai_1) >2:
                print('error occur')
            
            self.kaikikae.append((kai_1[0],kai_1[1],self.QUBO[kai_1[0]][kai_1[1]]))
            self.QUBO[kai_1[0]][kai_1[1]] += self.weight + self.Bay.max_w_list[ban] - self.Bay.score_disc_list[ban][(ban,)]
            #
            if ((kai_1[0],kai_1[1]) in self.seen_set_horizon_no_steps):
                print('omit')
                self.next_qubo()
            #
            #
        else:
            if self.steps == 0 and  len(self.branch_list[self.steps]) == self.chosen_list[self.steps] +1:
                cont = False
            else:
                
                a,b,value = self.kaikikae.pop()
                self.QUBO[a][b] = value
                self.chosen_list.pop()
                self.branch_list.pop()
                self.kai_list.pop()
                #
                c,d,e = self.kaikikae[-1]
                self.QUBO[c][d] = e
                self.seen_set_horizon = {(y,z,w) for (y,z,w) in self.seen_set_horizon if y < self.steps}
                self.seen_set_horizon_no_steps = {(z,w) for (y,z,w) in self.seen_set_horizon if y < self.steps}

                self.seen_set_horizon.add((self.steps -1,c,d))
                self.seen_set_horizon_no_steps.add((c,d))
                #
                print('seen',self.seen_set_horizon)
                print('seen_no',self.seen_set_horizon_no_steps)
                
                self.steps = self.steps -1
                self.next_qubo()
                
        


        return cont
    

    def depthsf(self):
        print('step',self.steps)
        print('seen',self.seen_set_horizon)
        print('seen_no',self.seen_set_horizon_no_steps)
        qubo_result =self.qubo_result(self.QUBO)
        parent_set = self.make_parent_set_list(qubo_result)
        qubo_result_np = np.array(qubo_result)
        pre_score = np.dot(np.dot((qubo_result_np),self.QUBO), (qubo_result_np.T))
        min_cycle = []
        min_cycle_num = len(parent_set) +1
        if pre_score >= self.uppor_bound:
            if self.next_qubo():
                self.depthsf()

        else:
            if judge_dag(parent_set):
                self.uppor_bound = pre_score
                self.preli_kai = qubo_result
                self.kai_parent_set = parent_set
                if self.next_qubo():
                    self.depthsf()
            else:
                self.steps = self.steps + 1
                pre_graph = preli_graph(len(parent_set),parent_set)
                woods = pre_graph.scc_tree
                for tree in woods:
                    if min_cycle_num > len(tree) >1:
                        min_cycle = tree
                self.kai_list.append(qubo_result)
                self.branch_list.append(min_cycle)
                self.chosen_list.append(0)
                ban = min_cycle[0]
                ban_domain = self.count_sum_len_of_u[ban]
                
                print(woods)
                print(min_cycle)
                kai_ban = self.kai_list[self.steps]
                kai_1 = [j for j ,x in enumerate(kai_ban) if x == 1 and ban_domain <= j < ban_domain + self.Bay.count_len_of_u[ban]]
                kai_1 = np.array(kai_1)
                kai_1 = np.sort(kai_1)

                print('lenkai',len(kai_1))
                print(kai_ban)
                print('kai_1',kai_1)
                print('tuika',(kai_1[0],kai_1[1],self.QUBO[kai_1[0]][kai_1[1]]))
                if len(kai_1) >2:
                    print('error occur')
                self.kaikikae.append((kai_1[0],kai_1[1],self.QUBO[kai_1[0]][kai_1[1]]))
                
                
                #print(self.QUBO[0][0])
                #print(self.QUBO[kai_1[0]][kai_1[0]])
                #print(self.QUBO[kai_1[1]][kai_1[1]])
                
                
                self.QUBO[kai_1[0]][kai_1[1]] += self.weight + self.Bay.max_w_list[ban] - self.Bay.score_disc_list[ban][(ban,)]          
                self.depthsf()


        return


bab = Branch_and_Bound(bayes,QUBO)
bab.depthsf()
print(bab.kai_parent_set)