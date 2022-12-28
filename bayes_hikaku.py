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
    i_to_parent_set = list()
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
        if self.max_parent_set <= 2:
            self.qubo_size = int((3*pow(self.num_variable,2) - self.num_variable)/2)
        elif  self.max_parent_set ==3 and self.num_variable % 2 == 1:
            self.qubo_size = int(self.num_variable*pow(self.num_variable +1,2)/4)
        elif self.max_parent_set == 3 and self.num_variable %2 == 0:
            self.qubo_size = int((self.num_variable*(pow(self.num_variable +1,2) + 1))/4)
        else:
            print('error occur')
        self.QUBO = np.zeros((self.qubo_size,self.qubo_size))
        



    def init_score_disc(self):
        self.score_disc_list = [dict() for i in range(self.num_variable)]

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






    def making_QUBO(self):
        num_vari = self.num_variable
        max_parent_num = self.max_parent_set
        
        xi = [0 for i in range(num_vari)]
        delta = [0 for i in range(num_vari)]
        
        
        #debug
        #u_list = [list() for j in range(num_vari)]
        #debug ijou
        for num in range(num_vari):
            print(num,'calc A')
            bayes.make_score_disc_i(num)
            
            
        
            
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
max_parent = 3
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