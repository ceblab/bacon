import numpy as np
class preli_graph:
    node_status_list = list()
    node_time_stamp= list()
    child_set_list = list()
    parent_set_list = list()
    scc_tree = []
    cycle_list = []

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
    def search_cycle(self,min_cycle,parent):
        self.cycle_list.clear()
        
        self.stlist = [0 for i in range(self.num_node)]
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
        print(self.cycle_list)
        self.stlist[vari] = 1
        con = True
        for i in parent[vari]:
            print(i)
            if i == start:
                print(self.cycle_list,'false')
                con = False
                break

            if self.stlist[i] == 0 and con:
                self.cycle_list.append(i)
                con = con and self.df_search(i,parent,start)
        if con:
            a = self.cycle_list.pop()
            print('pop',vari,a)

        return con

def dag_search(parents):
        
        leng = len(parents)
        topological = []
        void_q = []
        Judge = True
        for u in range(leng):
            for i in range(leng):
                if len(parents[i] )== 0:
                    void_q.append(i)
                    parents[i].append(-1)
            if len(void_q) == 0:
                
                Judge = False
                break
            x = void_q.pop(0)
            topological.append(x)
            for v in range(leng):
                if x in parents[v]:
                    parents[v].remove(x)
        
    
        return   Judge


list1 =[[], [2], [], [4, 6], [2], [0, 1, 3], [5], [4, 5]] 
list2 = [[], [3, 5], [], [2, 5], [2], [2, 6, 7], [0, 1, 3], [3, 4]]#[[1],[2,3],[4],[0],[1]]

a = preli_graph(len(list2),list2)
print(a.scc_tree)
print(dag_search(list2))
print(a.search_cycle([1,3,5,6,7],list2))