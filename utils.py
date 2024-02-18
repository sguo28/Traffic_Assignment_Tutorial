# Author = "Shuocheng Guo"
# Date = "6/9/2023"

import numpy as np
import pandas as pd
import networkx as nx
from itertools import islice
from collections import defaultdict
from queue import PriorityQueue
from scipy.optimize import fsolve
def find_k_shortest_paths(G, source, target, k, weight=None):
    return list(
        islice(nx.shortest_simple_paths(G, source, target, weight=weight), k)
    )

class Path:
    def __init__(self, path):
        self.path = path
        self.traversed_edges = [(path[i], path[i + 1]) for i in range(len(path) - 1)]
        self.travel_time = None
    def update_travel_time(self, link_time):
        self.travel_time = sum([link_time[edge] for edge in self.traversed_edges])

class TrafficAssignment:
    def __init__(self, file_path, od_path, node_attr, instance_name="sample", algorithm_name ="FW", eps=1e-3):
        G = self.create_network(file_path, node_attr)
        self.graph=G
        od, od_demand = self.get_od(od_path=od_path, od_mask=None)
        self.od=od
        self.od_demand=od_demand
        self.beta = self.get_edge_attributes('beta')
        self.alpha = self.get_edge_attributes('alpha')
        self.fftt = self.get_edge_attributes("travel_time")
        self.capacity = self.get_edge_attributes("capacity")
        self.eps=eps
        self.edges = list(self.graph.edges())
        self.nodes = list(self.graph.nodes())
        self.num_edges = len(self.edges)
        self.num_nodes = len(self.nodes)
        self.link_time = np.zeros(self.num_edges)
        self.link_flow = np.zeros(self.num_edges)
        self.iteration = 0
        self.converged = False
        self.old_gap=1e9
        self.UE_sol = 1e6*np.ones(self.num_edges)
        self.algorithm_name = algorithm_name
        if self.algorithm_name in ["FW","MSA"]:
            self.pqs = None
            self.link_flows = None
        elif self.algorithm_name == "PBA":
            self.pqs = [[] for _ in range(len(self.od))] 
            self.list_flows = [[] for _ in range(len(self.od))]
        else:
            raise ValueError(f"Alogrithm name {self.algorithm_name} not support.")
        if instance_name == 'SiouxFalls':
            self.UE_sol_LB = 4231335.28710744
        elif instance_name == 'sample':
            self.UE_sol_LB = 1719653.819467463
        else:
            raise ValueError("instance name not found")
        
    def get_od(self, od_path, od_mask=None):
        od = pd.read_csv(od_path)
        ods = od[['O', 'D']].values.tolist()
        ods_demand = od['Ton'].values.tolist()
        if od_mask is not None:
            od_mask = od_mask[:od.shape[0]]
            ods = [od for i,od in zip(od_mask,ods) if i==1]
            ods_demand = [od_demand for i,od_demand in zip(od_mask,ods_demand) if i==1]
        return ods, ods_demand

    def get_edge_attributes(self, attributes):
        '''
        get edge attributes: capacity, beta, alpha, travel_time, traffic_flow
        '''
        vals = nx.get_edge_attributes(self.graph,attributes).values()
        return np.array(list(vals))
    
    def get_graph(self):
        return self.graph

    @staticmethod
    def create_network(filepath,node_attr_path):
        '''
        Create a network from a csv file
        '''
        net = pd.read_csv(filepath)
        node_attr = pd.read_csv(node_attr_path)
        G = nx.DiGraph()
        for _, row in net.iterrows():  # iterate over all rows
            edge_id, fn, tn, tt_, capacity_, beta_, alpha_ = row  # unpack your attributes
            G.add_edge(fn, tn, edge_index = edge_id, travel_time=tt_, travel_flow=0, capacity=capacity_, beta=beta_, alpha=alpha_)
        # add node attributes
        for _, row in node_attr.iterrows():
            node, x, y = row
            G.add_node(node, x=x, y=y)

        return G

    def step(self, metric = "relative_gap", search_method = "bisection"):
        if self.algorithm_name == "FW":
            return self.step_fw(metric,search_method)
        elif self.algorithm_name == "MSA":
            return self.step_msa(metric)
        elif self.algorithm_name == "PBA": # Path-Based Algorithm
            return self.step_pba( metric,search_method)
        else:
            raise ValueError("algorithm name not found")
    def step_msa(self, metric = "relative_gap"):
        '''
        one step of traffic assignment
        metric = AEC or relative gap
        '''
        if self.iteration==0:
            init_shortest_path, init_shortest_path_cost=self.get_shortest_path()
            self.link_flow=self.map_path_to_link_flow(init_shortest_path)
        # record the current best link flow
        old_link_flow = self.link_flow
        self.save_UE_cost_function()

        # update travel time
        self.update_travel_time()

        # get the shortest path based on the new travel time
        shortest_path, shortest_path_cost = self.get_shortest_path()
        if metric == "relative_gap":
            # gap = self.get_gap(shortest_path_cost)
            gap = np.sum(self.UE_sol)/self.UE_sol_LB-1
            if gap < self.eps:
            # if self.iteration>100:
                self.converged = True
        else:
            raise ValueError(f"metric {metric} not supported.")
        # map the shortest path to link flow
        new_link_flow = self.map_path_to_link_flow(shortest_path)

        lbd = 1/(self.iteration+1)
        self.link_flow = self.convex_combination(old_link_flow, new_link_flow, lbd)

        self.update_link_flow() # self.link_flow
        self.update_travel_time()

        self.save_shortest_path(shortest_path)

        # update iteration
        self.iteration += 1
        return gap

    def step_fw(self, metric = "relative_gap",search_method = "bisection"):
        '''
        one step of traffic assignment
        metric = AEC or relative gap
        '''
        if self.iteration==0:
            init_shortest_path, init_shortest_path_cost=self.get_shortest_path()
            self.link_flow=self.map_path_to_link_flow(init_shortest_path)
        # record the current best link flow
        old_link_flow = self.link_flow
        self.save_UE_cost_function()

        # update travel time
        self.update_travel_time()
        # print(f"check link_flow: {self.link_flow}\n travel time: {self.link_time}")

        # get the shortest path based on the new travel time
        shortest_path, shortest_path_cost = self.get_shortest_path()
        if metric == "relative_gap":
            # gap = self.get_gap(shortest_path_cost)
            gap = np.sum(self.UE_sol)/self.UE_sol_LB-1
            if gap < self.eps:
            # if self.iteration>100:
                self.converged = True
        else:
            raise ValueError(f"metric {metric} not supported.")
        # map the shortest path to link flow
        new_link_flow = self.map_path_to_link_flow(shortest_path)

        if search_method == "bisection":
            lbd, n_iter = self.bisection_search(old_link_flow, new_link_flow) # optimal lambda
            print('Iteration={},lambda={},bisection iterations={}'.format(self.iteration,lbd,n_iter))
            self.link_flow = self.convex_combination(old_link_flow, new_link_flow, lbd)
        elif search_method == "exact":
            tau = self.fsolve_search(old_link_flow, new_link_flow) # optimal lambda
            self.link_flow=old_link_flow+tau*(new_link_flow-old_link_flow)
            n_iter=-1 # a placeholder, not available for fsolve
        else:
            raise ValueError(f"search method {search_method} not supported.")

        self.update_link_flow() # self.link_flow
        self.update_travel_time()

        self.save_shortest_path(shortest_path)

        # update iteration
        self.iteration += 1
        return gap
    
    def update_travel_time_for_path_set(self, pq, basic_path_traversed_arcs):
        '''
        update travel time for a set of paths
        output: list of travel time for each path, and the basic path id (if has)
        '''
        used_flag = False
        used_path_cost = []
        basic_path_id = None

        i = 0
        # update travel costs for all paths within the same OD.
        for used_traversed_arcs in pq:
            if not used_flag:
                if set(basic_path_traversed_arcs) == set(used_traversed_arcs):
                    used_flag = True
                    basic_path_id = i

            used_path_cost.append(self.get_travel_time_from_path(used_traversed_arcs))
            i += 1
        return used_path_cost, basic_path_id
    
    def step_pba(self, metric = "relative_gap"):
        '''
        path-based algorithm with gradient projection
        '''
        # first, find the shortest path for each OD pair
        if self.iteration==0:
            for od_pair, od_demand, pq, list_flow in zip(self.od, self.od_demand,self.pqs,self.list_flows):
                init_sp, init_sp_cost = self.get_shortest_path_single_od(od_pair)
                init_traversed_arcs = [(i,j) for i,j in zip(init_sp[:-1],init_sp[1:])]
                # put the shortest path to the priority queue with the cost as the priority
                pq.append(init_traversed_arcs) # ([init_sp,init_sp_cost,od_demand]) # append shortest path and demand
                list_flow.append(od_demand)
                self.update_path_to_link_flow_single_od([],init_traversed_arcs, od_demand)
                self.update_travel_time()

        else:
            # next update new shortest path for each OD pair
            for od_pair, od_demand, pq, list_flow in zip(self.od, self.od_demand, self.pqs,self.list_flows):
                basic_path, basic_path_cost = self.get_shortest_path_single_od(od_pair) # based on the updated travel time
                basic_path_traversed_arcs = [(i,j) for i,j in zip(basic_path[:-1],basic_path[1:])]

                used_path_cost, basic_path_id = self.update_travel_time_for_path_set(pq, basic_path_traversed_arcs)

                if basic_path_id is None: # it's a new path:
                    basic_path_flow = 0
                else:
                    basic_path_flow = list_flow[basic_path_id]
                # for each non-basic paths, adjust the flow:
                for i in range(len(pq)):
                    if i != basic_path_id:
                        non_basic_traversed_arcs = pq[i]
                        non_basic_path_cost = used_path_cost[i]
                        # get the shift flow
                        path_cost_diff = non_basic_path_cost-basic_path_cost # non-negative
                        # then get dev_sum
                        dev_sum = self.get_sum_of_gradient(basic_path_traversed_arcs, non_basic_traversed_arcs)
                        flow_to_shift = min(list_flow[i],path_cost_diff/dev_sum)
                        list_flow[i] -= flow_to_shift
                        basic_path_flow += flow_to_shift

                        self.update_path_to_link_flow_single_od(non_basic_traversed_arcs,basic_path_traversed_arcs, flow_to_shift)

                # finally, update the basic path flow
                if basic_path_id is None: # it's a new path:
                    list_flow.append(basic_path_flow)
                    pq.append(basic_path_traversed_arcs)
                else:
                    list_flow[basic_path_id] = basic_path_flow

                # self.update_path_to_link_flow_single_od(used_traversed_arcs, sp_traversed_arcs, shift_flow)
                self.update_travel_time()
        # drop the paths with zero flow
        for od_pair, od_demand, pq, list_flow in zip(self.od, self.od_demand, self.pqs, self.list_flows):
            indices_to_remove = [i for i, flow in enumerate(list_flow) if flow < self.eps]
            for i in sorted(indices_to_remove, reverse=True):
                del pq[i]
                del list_flow[i]
        # check the relative gap
        if metric == "relative_gap":
            gap = np.sum(self.UE_sol)/(self.UE_sol_LB)-1 # initial one is sufficiently large
            if gap < self.eps:
                self.converged = True
        else:
            raise ValueError("metric should be AEC or relative_gap")

        print("iteration: {}, gap: {}".format(self.iteration, gap))

        self.iteration += 1
        # one more thing to do: update the travel time for each path under the same OD.
        self.save_UE_cost_function()
        return gap
    
    def get_travel_time_from_path(self, sp):
        '''
        get the travel time based on the shortest path (list of nodes)
        '''
        t = 0
        for i,j in sp:
            t += self.graph[i][j]['travel_time']
        return t

    def get_sum_of_gradient(self,sp1,sp0):
        '''
        get the sum of derivatives of links.
        '''
        # get the links in sp1 but not in sp0
        link_sp1 = set(sp1)
        link_sp0 = set(sp0)
        # find the links that contribute to either sp1 or sp0
        link_contributed = link_sp1.symmetric_difference(link_sp0)

        # get a masked edge list
        masked_edge_cost = np.zeros(len(self.edges))
        for link in link_contributed:
            masked_edge_cost[self.edges.index(link)] = 1
        dev_sum=np.sum(self.fftt*self.beta*self.alpha/self.capacity*(masked_edge_cost*self.link_flow/self.capacity)**(self.alpha-1))

        return dev_sum

    def update_path_to_link_flow_single_od(self, path_flow_out, path_flow_in, shift_flow):
        '''
        map the shortest path to link flow
        '''
        if len(path_flow_out)>0:
            for i,j in path_flow_out:
                self.link_flow[self.edges.index((i,j))] -= shift_flow


        for i,j in path_flow_in:
            self.link_flow[self.edges.index((i,j))] += shift_flow

    def save_UE_cost_function(self):
        '''
        save the UE cost function to the graph
        '''
        self.UE_sol = self.fftt*self.link_flow + self.capacity/(self.alpha+1)*self.fftt*self.beta*(self.link_flow/self.capacity)**(self.alpha+1)
    
    def save_shortest_path(self,shortest_path):
        '''
        save the shortest path to the graph
        '''
        self.graph._shortest_path = shortest_path
    
    @staticmethod
    def convex_combination(x, y, theta):
        return (1-theta)*x + theta*y


    def update_travel_time(self):
        # update link travel time
        self.link_time = self.fftt*(1+self.beta*(self.link_flow/self.capacity)**self.alpha)
        # then update the link_time to each edge attribute
        for e,tt in zip(self.edges,self.link_time):
            self.graph[e[0]][e[1]]["travel_time"]=tt #update link travel time in the graph

    def get_average_excess_cost(self,shortest_path_cost):
        '''
        get average excess cost
        AEC = \frac{t*x-k*d}{d*1}
        '''
        tx = np.sum(self.link_flow*self.link_time)
        kd = np.sum([self.od_demand[i]*shortest_path_cost[i] for i in range(len(self.od_demand))])
        return (tx-kd)/np.sum(self.od_demand)

    def get_gap(self,shortest_path_cost):
        '''
        get gap
        '''
        # AEC = self.get_average_excess_cost()
        # tx = np.sum(new_travel_time*new_link_flow)
        tx = np.sum(self.link_time*self.link_flow)
        kd = np.sum([self.od_demand[i]*shortest_path_cost[i] for i in range(len(self.od_demand))])
        # TTST = sum(old_travel_time * old_link_flow)
        # SPTT = sum(new_travel_time * new_link_flow)
        # return TTST/SPTT-1

        return tx/kd-1
    
    def update_link_flow(self):
        '''
        updated link flow (according to the shortest path)
        theta: step size (default: 1, bisection method)
        '''
        # update the link_flow to each edge attribute
        for e,f in zip(self.edges,self.link_flow):
            self.graph[e[0]][e[1]]["travel_flow"]=f

    def get_shortest_path_single_od(self,od):
        '''
        get the shortest path for a single od pair
        '''
        sp = nx.dijkstra_path(self.graph, od[0], od[1], weight='travel_time')
        sp_cost = np.sum([self.graph[i][j]["travel_time"] for (i,j) in zip(sp,sp[1:])])
        return sp, sp_cost

    def get_shortest_path(self):
        list_sp = [nx.dijkstra_path(self.graph, self.od[i][0], self.od[i][1], weight='travel_time') for i in range(len(self.od))]
        list_sp_cost = [np.sum([self.graph[i][j]["travel_time"] for (i,j) in zip(sp,sp[1:])]) for sp in list_sp]
        return list_sp, list_sp_cost

    def compute_gradient(self,x1,x2,lbd):
        '''
        x1: the current link flow assignment
        x2: the shortest path flow assignment
        lbd: lambda ratio
        '''
        flow=(1-lbd)*x1+lbd*x2
        link_cost=(self.fftt*(1+self.beta*(flow/self.capacity)**self.alpha))@(x2-x1)
        return link_cost

    def map_path_to_link_flow(self, path):
        '''
        map the shortest path to link flow
        '''
        link_flow = np.zeros(self.num_edges)

        for pid,p in enumerate(path):
            for i in range(len(p) - 1):
                link_flow[self.edges.index((p[i], p[i + 1]))] += self.od_demand[pid]
        return link_flow

    def bisection_search(self,x1,x2,eps=1e-8):
        '''
        bisection search         
        for the optimal step size: lambda
        the goal is to minimize the gap
        '''
        l = 0
        l_ = 1
        n_iter = 0  # iteration
        while abs(l_ - l) > eps:  # stopping criteria
            c = self.compute_gradient(x1, x2, (l + l_) / 2)
            if c > 0:
                l_ = l + (l_ - l) / 2
            else:
                l = l + (l_ - l) / 2
            n_iter += 1
        return (l + l_) / 2, n_iter

    def fsolve_search(self,x1,x2,eps=1e-8):
        '''
        the goal is to minimize the gap
        '''
        def formula():
            '''
            x1: the current link flow assignment
            x2: the shortest path flow assignment
            lbd: lambda ratio
            '''
            # flow = (1 - lbd) * x1 + lbd * x2
            d=x2-x1
            link_cost = (self.fftt * (1 + self.beta * (x1 / self.capacity) ** self.alpha))
            hessian=np.diag(self.fftt*self.beta*self.alpha/self.capacity*(x1/self.capacity)**(self.alpha-1))

            tau= -link_cost.T@d/(d.T@hessian@d)

            return tau

        tau = formula()
        return tau
