from calendar import c
from re import L, M
import numpy as np
import itertools
import copy
import math
import csv
import pprint
import pandas as pd
from copy import deepcopy

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
    def __init__(self,bayes,Qubo,limit):
        self.QUBO = Qubo
        self.QUBO_origin =np.array(copy.deepcopy(Qubo))
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
        self.step_limit = limit
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
                cont = self.next_qubo()
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
                cont = self.next_qubo()
                
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
            if self.steps > -1 and self.next_qubo():
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
                if self.steps>self.step_limit:
                    blist_sin = []
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
    def qubo_solve_i_sol(self,q):
        if q == len(self.count_sum_len_of_u) -1:
            k = self.qubo_size
        else: 
            k = self.count_sum_len_of_u[(q+1)]
        for u in range(self.count_sum_len_of_u[q]+1,k):
            self.QUBO[u,u] = self.weight
        self.depthsf()
        self.QUBO = np.copy(self.QUBO_origin)
        self.branch_list.clear()
        self.chosen_list.clear()
        self.seen_set_horizon.clear()
        self.seen_set_horizon_no_steps.clear()
        self.kakikae.clear()
        self.steps = -1
    def qubo_solve_i(self):
        k = 0
        for q in range(len(self.count_sum_len_of_u)):
            if k != 0:
                for u in range(self.count_sum_len_of_u[q-1]+1,k):
                    self.QUBO[u,u] = self.QUBO_origin[u,u]
            if q == len(self.count_sum_len_of_u) -1:
                k = self.qubo_size
            else: 
                k = self.count_sum_len_of_u[(q+1)]
            for u in range(self.count_sum_len_of_u[q]+1,k):
                    self.QUBO[u,u] = self.weight
            self.depthsf()
            self.QUBO = np.copy(self.QUBO_origin)
            self.branch_list.clear()
            self.chosen_list.clear()
            self.seen_set_horizon.clear()
            self.seen_set_horizon_no_steps.clear()
            self.kakikae.clear()
            self.steps = -1

    def qubo_solve_ij(self):
        
        
        for q in range(len(self.count_sum_len_of_u)):
            l = [k for k in range(len(self.count_sum_len_of_u)) if k != q]

            for j in l:
                print('l',l)
                print('j',j)
                print('q',q)
                if q == len(self.count_sum_len_of_u) -1:
                    k = self.qubo_size
                else: 
                    k = self.count_sum_len_of_u[(q+1)]
                for u in range(self.count_sum_len_of_u[q]+1,k):
                    self.QUBO[u,u] = self.weight
                for m in range(len(self.count_sum_len_of_u)):
                    if m != j:
                        for s in self.Bay.i_jdic[(m,j)]:
                            a,b = s
                            offset = self.count_sum_len_of_u[m]
                            self.QUBO[offset + a][offset + b]  =self.weight
                self.depthsf()
                self.QUBO = np.copy(self.QUBO_origin)
                self.branch_list.clear()
                self.chosen_list.clear()
                self.seen_set_horizon.clear()
                self.seen_set_horizon_no_steps.clear()
                self.kakikae.clear()
                self.steps = -1

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
bab = Branch_and_Bound(bayes,QUBO,13)
#bab.qubo_solve()
#bab.depthsf()
#bab.qubo_solve_ij()
#bab.qubo_solve_i()
time_end = time.time()
proc_time = time_end- time_sta
print('解',bab.kai_parent_set)
print('実行時間',proc_time)
list,parent_set_list = recreate(bab.preli_kai,bayes)
print('')

for i in range(bayes.num_variable):
    
     print(bayes.variable_list[i],'の親変数は',list[i],'です')
