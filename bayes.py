class Bayes:
    num_variable = 0
    def __init__(self,num_vari):
        self.num_variable = num_vari

    def output(self):
        print(self.num_variable)

bayes = Bayes(10)
bayes.output()