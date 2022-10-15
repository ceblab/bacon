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





bayes = Bayes(9,3)
bayes.output()
#s = set((1,2,3))
score = bayes.score_disc_list[1][(1,)]
print(score)
score = bayes.score_disc_list[1][(2,3,8)]
print(score)
