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
        print(void_q)
        print(parent_set_list)
    
    return judge

parent_set_list = [[], [2],[1]]

print(judge_dag(parent_set_list))