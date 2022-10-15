import numpy as np
import itertools

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
            



num_variable = 9
max_parent = 3
bayes = Bayes(num_variable,max_parent)
bayes.output()
#s = set((1,2,3))
score = bayes.score_disc_list[1][(1,)]
print(score)
score = bayes.score_disc_list[1][(2,3,8)]
print(score)

list = make_W_set_of_i(bayes.score_disc_list[1],1,num_variable,max_parent)

for i in list:
    print(i ,": score" ,bayes.score_disc_list[1][i])
