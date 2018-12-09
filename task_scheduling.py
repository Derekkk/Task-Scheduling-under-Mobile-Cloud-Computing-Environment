# File:
# -*- coding: utf-8 -*-
# @Time    : 11/29/2018 7:25 PM
# @Author  : Derek Hu
from copy import deepcopy
import time

# define node class
class Node(object):
    def __init__(self, id, parents, children, core_speed, cloud_speed):
        self.id = id
        self.parents = parents # list of Nodes
        self.children = children # list of Nodes
        self.core_speed = core_speed # list: [9, 7, 5] for core1, core2 and core3
        self.cloud_speed = cloud_speed # list [3, 1, 1] cloud speed
        self.ft_l = 0 # local finish time, inf at start
        self.ft_ws = 0
        self.ft_c = 0
        self.ft_wr = 0
        self.ready_time = -1 # local ready time
        self.rt_ws = -1 # cloud ready time
        self.rt_c = -1
        self.rt_wr = -1
        self.is_core = None
        self._local_cloud() # compute self.is_core
        self._computation_cost()
        self.priority_socre = None
        self.assignment = -2 # 0, 1, 2, 3
        self.start_time = [-1, -1, -1, -1] # start time for core1, core2, core3, cloud
        self.is_scheduled = None

    def _local_cloud(self):
        """determine if local or cloud, here assume core3 is faster than others"""
        t_l_min = self.core_speed[2]
        t_c_min = 5 # assume cloud is always 5
        if t_l_min <= t_c_min:
            self.is_core = True
            self.ft_ws = 0
            self.ft_c = 0
            self.ft_wr = 0
        else:
            self.is_core = False
            self.ft_l = 0

    def _computation_cost(self):
        """calculate w_i in section 3"""
        self.w_i = 0
        if self.is_core == True:
            self.w_i = sum(self.core_speed) / len(self.core_speed)
        else:
            self.w_i = 5


def total_T(nodes):
    """compute the total time"""
    total_t = 0
    for node in nodes:
        if len(node.children) == 0:
            total_t = max(node.ft_l, node.ft_wr)
    return total_t

def total_E(nodes, core_cloud_power=[1, 2, 4, 0.5]):
    """compute total energy
        core_cloud_power: [1, 2, 4, 0.5] for core1, core2, core3, cloud sending
    """
    total_energy = 0
    for node in nodes:
        if node.is_core == True:
            current_node_e = node.core_speed[node.assignment] * core_cloud_power[node.assignment]
            total_energy += current_node_e
        if node.is_core == False:
            current_node_e = node.cloud_speed[0] * core_cloud_power[3]
            total_energy += current_node_e
    return total_energy

def primary_assignment(nodes):
    """primary assignment, input is a list of nodes"""
    # [c1, c2, c3, sending, cloud, receiving]
    local_source = [0, 0, 0]
    cloud_source = [0, 0, 0]
    # core and cloud sequence after assignment
    core1_seq = []
    core2_seq = []
    core3_seq = []
    cloud_seq = [] # here we assume cloud cost is same for each node, and the seq of 3 process are always the same

    for i, node in enumerate(nodes):
        if node.is_core == True: # local task
            if len(node.parents) == 0:
                node.ready_time = 0
            else: # equation (3)
                for p in node.parents:
                    p_ft = max(p.ft_l, p.ft_wr)
                    if p_ft > node.ready_time:
                        node.ready_time = p_ft

            core_1_finishtime = max(local_source[0], node.ready_time) + node.core_speed[0]
            core_2_finishtime = max(local_source[1], node.ready_time) + node.core_speed[1]
            core_3_finishtime = max(local_source[2], node.ready_time) + node.core_speed[2]
            # choose the fastest one
            core_assign_id = 0
            core_assign_finishtime = core_1_finishtime
            if core_assign_finishtime > core_2_finishtime:
                core_assign_finishtime = core_2_finishtime
                core_assign_id = 1
            if core_assign_finishtime > core_3_finishtime:
                core_assign_finishtime = core_3_finishtime
                core_assign_id = 2
            node.assignment = core_assign_id
            node.ft_l = core_assign_finishtime
            node.start_time[core_assign_id] = max(local_source[core_assign_id], node.ready_time)
            # update the ready time of the children nodes
            # current_node_ft = max(node.ft_l, node.ft_wr)
            # for child_node in node.children:
            #     if child_node.ready_time < current_node_ft:
            #         child_node.ready_time = current_node_ft
            # local source assignment
            local_source[core_assign_id] = node.ft_l

            # add the node to core seq
            if node.assignment == 0:
                core1_seq.append(node.id)
            if node.assignment == 1:
                core2_seq.append(node.id)
            if node.assignment == 2:
                core3_seq.append(node.id)

            print("node id:{}, assigenment:{}, ready time: {}, local start_time: {}".
                  format(node.id, node.assignment+1, node.ready_time, node.start_time[node.assignment]))
            print(local_source)
            print("-----------")

        if node.is_core == False: # cloud task
            # equation (4)
            for p in node.parents:
                p_ws = max(p.ft_l, p.ft_ws)
                if p_ws > node.rt_ws:
                    node.rt_ws = p_ws
            cloud_ws_finishtime = max(cloud_source[0], node.rt_ws) + node.cloud_speed[0]
            node.ft_ws = cloud_ws_finishtime
            # (5)
            p_max_ft_c = 0
            for p in node.parents:
                if p.ft_c > p_max_ft_c:
                    p_max_ft_c = p.ft_c
            node.rt_c = max(node.ft_ws, p_max_ft_c)
            cloud_c_finishtime = max(cloud_source[1], node.rt_c) + node.cloud_speed[1]
            node.ft_c = cloud_c_finishtime
            #(6)
            node.rt_wr = node.ft_c
            cloud_wr_finishtime = max(cloud_source[2], node.rt_wr) + node.cloud_speed[2]
            node.ft_wr = cloud_wr_finishtime
            node.assignment = 3 # 3 is cloud
            node.start_time[3] = max(cloud_source[0], node.rt_ws) # cloud task start time is sending start time

            cloud_source[0] = cloud_ws_finishtime
            cloud_source[1] = cloud_c_finishtime
            cloud_source[2] = cloud_wr_finishtime

            cloud_seq.append(node.id)
            print("node id:{}, assigenment:{}, ws ready time: {}, c ready time: {}, wr ready time: {}, cloud start time: {}".
                  format(node.id, node.assignment + 1, node.rt_ws, node.rt_c, node.rt_wr, node.start_time[3]))
            print(local_source)
            print("-----------")
    print(cloud_source)
    print("total time: ", total_T(nodes))
    print("total energy: ", total_E(nodes, [1, 2, 4, 0.5]))
    print([i for i in core1_seq])
    print([i for i in core2_seq])
    print([i for i in core3_seq])
    print([i for i in cloud_seq])
    seq = [core1_seq, core2_seq, core3_seq, cloud_seq]
    return seq


def new_squence(nodes, tar_id, k, seq):
    """
    compute new scheduling seq
    :param nodes: node list
    :param tar_id: index of target node
    :param k: migration location: [0, 1, 2, 3] means: core1, core2, core3, cloud
    :param seq: current core sequence: [core1_seq, core2_seq, core3_seq, cloud_seq], each one is a list of node_ids
    :return:
    """
    node_index = {} #{key-node.id: value-index in nodes}
    temp_id = 0
    for _node in nodes:
        node_index[_node.id] = temp_id
        temp_id += 1
        if _node.id == tar_id:
            node_tar = _node
    if node_tar.is_core == True: # calculate tar ready time in (19)
        node_tar_rt = node_tar.ready_time
    if node_tar.is_core == False:
        node_tar_rt = node_tar.rt_ws
    seq[node_tar.assignment].remove(node_tar.id) # original core seq
    s_new = seq[k] # S_new in (19)
    s_new_prim = []
    flag = False
    for _node_id in s_new:
        _node = nodes[node_index[_node_id]]
        if _node.start_time[k] < node_tar_rt:
            s_new_prim.append(_node.id)
        if _node.start_time[k] >= node_tar_rt and flag == False:
            s_new_prim.append(node_tar.id)
            flag = True
        if _node.start_time[k] >= node_tar_rt and flag == True:
            s_new_prim.append(_node.id)
    if flag == False:
        s_new_prim.append(node_tar.id)
    seq[k] = s_new_prim
    node_tar.assignment = k
    if k == 3:
        node_tar.is_core = False
    else:
        node_tar.is_core = True

    return seq


def kernel_algorithm(nodes_new, seq_new):
    """
    kernel algorithm
    :param nodes_new: node list
    :param seq_new: current core sequence: [core1_seq, core2_seq, core3_seq, cloud_seq], each one is a list of nodes
    """

    local_source = [0, 0, 0]
    cloud_source = [0, 0, 0]

    ready1 = [-1]*len(nodes_new) # [-1s] at start, ready1[i] is for node.id==i
    ready2 = [-1]*len(nodes_new)
    ready1[nodes_new[0].id - 1] = 0 # id start from 1
    for each_seq in seq_new:
        if len(each_seq) > 0:
            ready2[each_seq[0] - 1] = 0
    # print(ready1, ready2)

    node_index = {}  # {key-node.id: value-index in nodes}
    temp_id = 0
    for _node in nodes_new:
        node_index[_node.id] = temp_id
        _node.ready_time = -1  # local ready time
        _node.rt_ws = -1  # cloud ready time
        _node.rt_c = -1
        _node.rt_wr = -1
        temp_id += 1

    # start the rescheduling task
    stack = [] # LIFO stack
    stack.append(nodes_new[0])

    while len(stack) != 0: # not empty
        v_i = stack.pop()
        v_i.is_scheduled = "kernel_scheduled" # means is scheduled
        # first, calculate v_i local ready time
        if v_i.is_core == True: # local task
            if len(v_i.parents) == 0:
                v_i.ready_time = 0
            else: # equation (3)
                for p in v_i.parents:
                    p_ft = max(p.ft_l, p.ft_wr)
                    if p_ft > v_i.ready_time:
                        v_i.ready_time = p_ft

        # part 2: schedule on the corresponding core
        if v_i.assignment == 0: # local core1
            v_i.start_time = [-1, -1, -1, -1]
            v_i.start_time[0] = max(local_source[0], v_i.ready_time)
            v_i.ft_l = v_i.start_time[0] + v_i.core_speed[0]
            v_i.ft_ws = -1
            v_i.ft_c = -1
            v_i.ft_wr = -1
            local_source[0] = v_i.ft_l
        if v_i.assignment == 1: # local core2
            v_i.start_time = [-1, -1, -1, -1]
            v_i.start_time[1] = max(local_source[1], v_i.ready_time)
            v_i.ft_l = v_i.start_time[1] + v_i.core_speed[1]
            v_i.ft_ws = -1
            v_i.ft_c = -1
            v_i.ft_wr = -1
            local_source[1] = v_i.ft_l
        if v_i.assignment == 2: # local core3
            v_i.start_time = [-1, -1, -1, -1]
            v_i.start_time[2] = max(local_source[2], v_i.ready_time)
            v_i.ft_l = v_i.start_time[2] + v_i.core_speed[2]
            v_i.ft_ws = -1
            v_i.ft_c = -1
            v_i.ft_wr = -1
            local_source[2] = v_i.ft_l

        if v_i.assignment == 3: # cloud
            if len(v_i.parents) == 0: # 1. sending
                v_i.rt_ws = 0
            else:
                for p in v_i.parents:
                    p_ws = max(p.ft_l, p.ft_ws)
                    if p_ws > v_i.rt_ws:
                        v_i.rt_ws = p_ws
            v_i.ft_ws = max(cloud_source[0], v_i.rt_ws) + v_i.cloud_speed[0]
            v_i.start_time = [-1, -1, -1, -1]
            v_i.start_time[3] = max(cloud_source[0], v_i.rt_ws)
            cloud_source[0] = v_i.ft_ws

            p_max_ft_c = 0 # 2. cloud part
            for p in v_i.parents:
                if p.ft_c > p_max_ft_c:
                    p_max_ft_c = p.ft_c
            v_i.rt_c = max(v_i.ft_ws, p_max_ft_c)
            v_i.ft_c = max(cloud_source[1], v_i.rt_c) + v_i.cloud_speed[1]
            cloud_source[1] = v_i.ft_c

            v_i.rt_wr = v_i.ft_c # 3. receiveing part
            v_i.ft_wr = max(cloud_source[2], v_i.rt_wr) + v_i.cloud_speed[2]
            v_i.ft_l = -1
            cloud_source[2] = v_i.ft_wr


        # if v_i.is_core == True:
        #     print("node id:{}, assigenment:{}, ready time: {}, local start_time: {}, is_scheduled: {}".
        #           format(v_i.id, v_i.assignment + 1, v_i.ready_time, v_i.start_time[v_i.assignment], v_i.is_scheduled))
        #     print(local_source, cloud_source)
        #     print("-----------")
        # else:
        #     print(
        #         "node id:{}, assigenment:{}, ws ready time: {}, c ready time: {}, wr ready time: {}, cloud start time: {}, is_scheduled: {}".
        #         format(v_i.id, v_i.assignment + 1, v_i.rt_ws, v_i.rt_c,v_i.rt_wr, v_i.start_time[3], v_i.is_scheduled))
        #     print(local_source, cloud_source)
        #     print("-----------")


        #update ready1 and ready2
        corresponding_seq = seq_new[v_i.assignment]  # the sequence that current v_i is assigned

        v_i_index = corresponding_seq.index(v_i.id)  # position of v_i in seq list
        if v_i_index != len(corresponding_seq) - 1:
            next_node_id = corresponding_seq[v_i_index + 1]
        else:
            next_node_id = -1 # current node is the last in the seq

        for _node in nodes_new:
            flag = 0
            for p in _node.parents:
                if p.is_scheduled != "kernel_scheduled":
                    flag += 1
                ready1[_node.id - 1] = flag
            if _node.id == next_node_id:
                ready2[_node.id-1] = 0

        for _node in nodes_new:
            # add node into stack if satisfied
            if (ready1[_node.id-1] == 0) and (ready2[_node.id-1] == 0) and (_node.is_scheduled != "kernel_scheduled") and (_node not in stack):
                # print("add stack: ", _node.id)
                stack.append(_node)

    for node in nodes_new:
        node.is_scheduled = None
    return nodes_new

if __name__ == '__main__':

    # node10 = Node(id=10, parents=None, children=[], core_speed=[7, 4, 2], cloud_speed=[3, 1, 1])
    # node9 = Node(id=9, parents=None, children=[node10], core_speed=[5, 3, 2], cloud_speed=[3, 1, 1])
    # node8 = Node(id=8, parents=None, children=[node10], core_speed=[6, 4, 2], cloud_speed=[3, 1, 1])
    # node7 = Node(id=7, parents=None, children=[node10], core_speed=[8, 5, 3], cloud_speed=[3, 1, 1])
    # node6 = Node(id=6, parents=None, children=[node8], core_speed=[7, 6, 4], cloud_speed=[3, 1, 1])
    # node5 = Node(id=5, parents=None, children=[node9], core_speed=[5, 4, 2], cloud_speed=[3, 1, 1])
    # node4 = Node(id=4, parents=None, children=[node8, node9], core_speed=[7, 5, 3], cloud_speed=[3, 1, 1])
    # node3 = Node(id=3, parents=None, children=[node7], core_speed=[6, 5, 4], cloud_speed=[3, 1, 1])
    # node2 = Node(id=2, parents=None, children=[node8, node9], core_speed=[8, 6, 5], cloud_speed=[3, 1, 1])
    # node1 = Node(id=1, parents=None, children=[node2, node3, node4, node5, node6], core_speed=[9, 7, 5], cloud_speed=[3, 1, 1])
    # node1.parents = []
    # node2.parents = [node1]
    # node3.parents = [node1]
    # node4.parents = [node1]
    # node5.parents = [node1]
    # node6.parents = [node1]
    # node7.parents = [node3]
    # node8.parents = [node2, node4, node6]
    # node9.parents = [node2, node4, node5]
    # node10.parents = [node7, node8, node9]
    #
    # node2.is_core = False
    # node2._computation_cost()


    #Test 2
    node10 = Node(id=10, parents=None, children=[], core_speed=[7, 4, 2], cloud_speed=[3, 1, 1])
    node9 = Node(id=9, parents=None, children=[node10], core_speed=[5, 3, 2], cloud_speed=[3, 1, 1])
    node8 = Node(id=8, parents=None, children=[node10], core_speed=[6, 4, 2], cloud_speed=[3, 1, 1])
    node7 = Node(id=7, parents=None, children=[node10], core_speed=[8, 5, 3], cloud_speed=[3, 1, 1])
    node6 = Node(id=6, parents=None, children=[node8, node9], core_speed=[7, 6, 4], cloud_speed=[3, 1, 1])
    node5 = Node(id=5, parents=None, children=[node7, node8], core_speed=[5, 4, 2], cloud_speed=[3, 1, 1])
    node4 = Node(id=4, parents=None, children=[node6], core_speed=[7, 5, 3], cloud_speed=[3, 1, 1])
    node3 = Node(id=3, parents=None, children=[node5, node6], core_speed=[6, 5, 4], cloud_speed=[3, 1, 1])
    node2 = Node(id=2, parents=None, children=[node5], core_speed=[8, 6, 5], cloud_speed=[3, 1, 1])
    node1 = Node(id=1, parents=None, children=[node2, node3, node4], core_speed=[9, 7, 5], cloud_speed=[3, 1, 1])
    node1.parents = []
    node2.parents = [node1]
    node3.parents = [node1]
    node4.parents = [node1]
    node5.parents = [node2, node3]
    node6.parents = [node3, node4]
    node7.parents = [node5]
    node8.parents = [node5, node6]
    node9.parents = [node6]
    node10.parents = [node7, node8, node9]

    start = time.clock()

    node1.ready_time = 0

    # calculate Priority Score
    node_list = [node10, node9, node8, node7, node6, node5, node4, node3, node2, node1]
    for node in node_list:
        priority_socre = node.w_i
        if len(node.children) == 0:
            node.priority_socre = priority_socre
            continue
        child_score = max([i.priority_socre for i in node.children])

        node.priority_socre = priority_socre + child_score

    node_list = sorted(node_list, key=lambda node: node.priority_socre, reverse=True)

    print("compute priority order")
    for node in node_list:
        print(node.id, node.priority_socre)

    sequence = primary_assignment(node_list)
    T_init = total_T(node_list)
    E_init = total_E(node_list, [1, 2, 4, 0.5])
    print("initial time and energy: ", T_init, E_init)

    #############################################
    # start outer loop
    #############################################
    iter_num = 0
    while iter_num < 100:
        # One outer loop
        print("-----" * 20)
        print("iter: ", iter_num)
        print("-----" * 20)
        T_init = total_T(node_list)
        E_init = total_E(node_list, [1, 2, 4, 0.5])
        print("initial time and energy: ", T_init, E_init)
        migration_choice = [[] for i in range(len(node_list))]
        for i in range(len(node_list)):
            if node_list[i].assignment == 3: # cloud node
                current_row_id = node_list[i].id - 1
                current_row_value = [1] * 4  # 4 resources
                migration_choice[current_row_id] = current_row_value
            else:
                current_row_id = node_list[i].id - 1
                current_row_value = [0] * 4 # 4 resources
                current_row_value[node_list[i].assignment] = 1
                migration_choice[current_row_id] = current_row_value
        print("migration")
        # print(migration_choice)

        T_max_constraint = 27
        result_table = [[(-1, -1) for j in range(4)] for i in range(len(node_list))]
        # print(result_table)

        for n in range(len(migration_choice)): # the n-th node
            nth_row = migration_choice[n]
            for k in range(len(nth_row)): # the k-th resource
                if nth_row[k] == 1:
                    continue
                seq_copy = deepcopy(sequence)
                nodes_copy = deepcopy(node_list)
                seq_copy = new_squence(nodes_copy, n+1, k, seq_copy)
                kernel_algorithm(nodes_copy, seq_copy)

                current_T = total_T(nodes_copy)
                current_E = total_E(nodes_copy)
                del nodes_copy
                del seq_copy
                result_table[n][k] = (current_T, current_E)

        # print("after outer loop: ")
        # print(result_table)

        n_best = -1
        k_best = -1
        T_best = T_init
        E_best = E_init
        ration_best = -1
        print(result_table)
        for i in range(len(result_table)):
            for j in range(len(result_table[i])):
                val = result_table[i][j]
                if val == (-1, -1):
                    continue
                if val[0] > 27:
                    continue
                ration = (E_best - val[1]) / abs(val[0] - T_best + 0.00005)
                if ration > ration_best:
                    ration_best = ration
                    n_best = i
                    k_best = j

        if n_best == -1 and k_best == -1:
            break
        n_best += 1
        k_best += 1
        T_best, E_best = result_table[n_best-1][k_best-1]
        print("current migration: task:{}, k: {}, total time: {}, total energy: {}".format(n_best, k_best, T_best, E_best))
        print("update after current outer loop")
        sequence = new_squence(node_list, n_best, k_best-1, sequence)
        kernel_algorithm(node_list, sequence)
        # print("finish time: ", [(node.id, node.start_time, node.ft_wr, node.ft_l, node.assignment) for node in node_list])
        for s in sequence:
            print([i for i in s])
        T_current = total_T(node_list)
        E_current = total_E(node_list, [1, 2, 4, 0.5])
        E_diff = E_init - E_current
        T_diff = abs(T_current - T_init)
        iter_num += 1

        print("time and energy: ", T_current, E_current)
        if E_diff <= 1:
            break

    print("rescheduling finished...")
    for node in node_list:
        if node.is_core == True:
            print("node id:{}, assigenment:{}, ready time: {}, local start_time: {}".
                  format(node.id, node.assignment + 1, node.ready_time, node.start_time[node.assignment]))
            print("-----------")
        else:
            print(
                "node id:{}, assigenment:{}, ws ready time: {}, c ready time: {}, wr ready time: {}, cloud start time: {}".
                format(node.id, node.assignment + 1, node.rt_ws, node.rt_c, node.rt_wr, node.start_time[3]))
            print("-----------")

    elapsed = (time.clock() - start)
    print("Time used:", elapsed)
    print("final sequence: ")
    for s in sequence:
        print([i for i in s])
    # for node in node_list:
    #     print(node.id, node.ft_l, node.ft_wr)
    T_final = total_T(node_list)
    E_final = total_E(node_list, [1, 2, 4, 0.5])
    print("final time: {}, final energy: {}".format(T_current, E_current))



    # print("rescheduling")
    # seq_copy = deepcopy(sequence)
    # nodes_copy = deepcopy(node_list)
    #
    # seq_copy = new_squence(nodes_copy, 1, 3, seq_copy)
    # for s in seq_copy:
    #     print([i for i in s])
    #
    # print("call kernel algorithm")
    # kernel_algorithm(nodes_copy, seq_copy)
    # for s in seq_copy:
    #     print([i for i in s])
    # print("time: ", total_T(nodes_copy), total_E(nodes_copy))
    # print("---"*30)
    # seq_copy = new_squence(nodes_copy, 3, 3, seq_copy)
    # for s in seq_copy:
    #     print([i for i in s])
    #
    # print("call kernel algorithm")
    # kernel_algorithm(nodes_copy, seq_copy)
    # for s in seq_copy:
    #     print([i for i in s])
    # print("time: ", total_T(nodes_copy), total_E(nodes_copy))
    # print("---" * 30)
    # seq_copy = new_squence(nodes_copy, 6, 3, seq_copy)
    # for s in seq_copy:
    #     print([i for i in s])
    #
    # print("call kernel algorithm")
    # kernel_algorithm(nodes_copy, seq_copy)
    # for s in seq_copy:
    #     print([i for i in s])
    # print("time: ", total_T(nodes_copy), total_E(nodes_copy))

