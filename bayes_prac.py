from pgmpy.readwrite import BIFReader
from pgmpy.sampling import BayesianModelSampling
reader = BIFReader("asia.bif")

samples = BayesianModelSampling(reader.get_model()).forward_sample(size=int(1e3))
#print(samples.head())
print((samples['asia'] == 'no'))
print((samples['asia'] == 'yes').sum())
#print(reader.get_variables())#変数名と変数番号の対応表になる
#print(reader .get_states())
