# Benders decomposition code for Drone-Truck routing problem.
# Written by Chungmok Lee <chungmok@gmail.com> and Munjeong Kang.
# Last major revision in Nov, 2020
# This code implements the algorithm represented in:
#   An Exact Algorithm for Heterogeneous Drone-Truck Routing Problem, 
#   Munjeong Kang and Chungmok Lee, Transportation Science, 2021
# 
# Works with Python 3.7 and Cplex 12.10



import copy
import cplex
import matplotlib.pyplot as plt
import networkit as nk
import numba
import numpy as np
import random
import re
import sys
import time

from map_util import * 


class MySolve(cplex.callbacks.UserCutCallback):
    def __call__(self):
        self.all_sols = {
            name: self.get_values(name)
            for name in self.all_names
        }

class LazyCutGenCallbackTSP(cplex.callbacks.LazyConstraintCallback):
    def __call__(self):
        self.times_called += 1

        cuts = cut_generation(self.prob, self.get_values(self.var_x_index))

        # print(cuts)

        for (cut, S, k) in cuts:
            self.add(cut, 'G', 0)
            self.generated_cuts.append((cut, S, k))


class CutGenCallback(cplex.callbacks.UserCutCallback):
    def __call__(self):
        self.times_called += 1

        start_time = time.time()

        cuts = cut_generation(self.prob, self.get_values(self.var_x_index))
        for (cut, S, k) in cuts:
            self.add(cut, 'G', 0)
            self.generated_cuts.append((cut, S, k))

        self.times_gcs += len(cuts)
        self.solve_time_gcs += time.time() - start_time


class LazyCutGenCallback(cplex.callbacks.LazyConstraintCallback):
    def __call__(self):
        self.times_called += 1

        start_time = time.time()
        cuts = cut_generation(self.prob, self.get_values(self.var_x_index))

        if len(cuts) > 0:

            self.times_gcs += len(cuts)

            for (cut, S, k) in cuts:
                self.add(cut, 'G', 0)
                self.generated_cuts.append((cut, S, k))

        self.solve_time_gcs += time.time() - start_time

        if len(cuts) == 0:

            start_time = time.time()
            bsp = self.bsp

            bmp_sol = {
                i: self.get_values(self.var_z_index)[i - 1]
                for i in self.prob.N
            }

            for i in self.prob.N:
                bsp.cpx.linear_constraints.set_rhs(bsp.node_visit_const_idx[i], 1 - bmp_sol[i])
                bsp.cpx.linear_constraints.set_rhs(bsp.node_visit_2_const_idx[i], len(self.prob.N) * bmp_sol[i])

            cut, rhs, cut_type = get_cplex_constraint_from_sp(bsp, self.prob, self.get_values(self.var_z_index),
                                                              self.get_values(self.var_W_index), self.w_star_dic,
                                                              cplex.infinity)  # incumbent-cur_tot_vehicle_time)  #cplex.infinity)

            if cut:
                if cut_type == 'inf':
                    self.times_inf += 1
                    self.add(cut, 'G', rhs)
                if cut_type == 'opt':
                    self.times_opt += 1
                    self.add(cut, 'G', rhs)

                self.generated_Benders_cuts.append((cut_type, cut, rhs))

            self.solve_time_sp += time.time() - start_time


def get_gcs(positive_sols, t, epsilon=0.8, num_cuts=-1):
    G = nk.graph.Graph()

    for (i, j), val in positive_sols.items():
        if not G.hasEdge(j, i):
            G.addEdge(i, j, val, addMissing=True)

    cc = nk.components.ConnectedComponents(G)
    cc.run()

    all_S = [set(S) for S in cc.getComponents() if len(S) >= 2 and t not in S]

    all_cuts = []

    for S in all_S:
        for k in S:
            sum_cut_S = sum(val for (i, j), val in positive_sols.items() if i in S and j not in S)
            sum_cut_i = sum(val for (i, j), val in positive_sols.items() if i == k)

            if sum_cut_i - sum_cut_S >= epsilon:
                all_cuts.append((sum_cut_i - sum_cut_S, S, k))

    if num_cuts > 0:
        all_cuts = sorted(all_cuts, reverse=True)
        all_cuts = all_cuts[:num_cuts]
        print(all_cuts)

    return all_cuts

def cut_generation(prob, sol):

    positive_sols = {a: val for (a, val) in zip(prob.A, sol) if val > 0.01}
    all_gcses = get_gcs(positive_sols, prob.t)

    cuts = []

    for violation, S, k in all_gcses:
        cut = get_cplex_constraint_from_gcs(prob, S, k)
        cuts.append((cut, S, k))

    return cuts

def get_cplex_constraint_from_gcs(prob, S, k):

    delta_S = set([(i, j) for (i, j) in prob.A if i in S and j not in S])
    delta_k = set([(i, j) for (i, j) in prob.A if i == k])
    inter_delta = delta_S.intersection(delta_k)
    delta_S -= inter_delta
    delta_k -= inter_delta

    return cplex.SparsePair(
        ind=[f'x_{i}_{j}' for (i, j) in delta_S] + [f'x_{k}_{j}' for (i, j) in delta_k],
        val=[1] * len(delta_S) + [-1] * len(delta_k)
    )

def get_cplex_constraint_from_sp(bsp, prob, z_sol, W_sol, w_star_dic, target_W):

    positive_sols = {a: val for (a, val) in zip(prob.N, z_sol) if val > 0.01}
    zero_sols = {a: val for (a, val) in zip(prob.N, z_sol) if val < 0.01}
    bsp_results = bsp.solve(upper_cutoff=target_W)


    if bsp_results is None or bsp.cpx.solution.get_status() == 103:  # infeasible
        return cplex.SparsePair(
            ind=[f'z_{i}' for i in zero_sols],
            val=[1] * len(zero_sols)
        ), 1, 'inf'

    else:
        bsp_obj = bsp_results['obj'] * prob.scaler
        if W_sol - bsp_obj < -0.1:
            p_star = len(positive_sols)
            omega = 0
            W = int(bsp_obj+ 0.1)
            for add_route_num, visit_node_num in enumerate(range(p_star + 1, int(list(w_star_dic.keys())[-1]) + 1)):
                omega = max(omega, ((W - w_star_dic[visit_node_num]) / (add_route_num + 1)))

            return cplex.SparsePair(
                ind=[f'W'] + [f'z_{i}' for i in zero_sols],
                val=[1] + [omega] * len(zero_sols)
            ), W, 'opt'

        else:
            return None, None, None

class Prob:
    def __init__(self, path, problem_name, alpha, scaler, drone_battery=100, drone_multiplier=0.4, truck_service_time=10, drone_service_time=5):

        self.path = path
        self.problem_name = problem_name
        self.alpha = alpha
        self.B = [round(drone_battery/value,2) for value in alpha]
        self.truck_service_time = truck_service_time
        self.drone_service_time = drone_service_time
        self.drone_multiplier = [round(drone_multiplier/value,2) for value in alpha]
        self.scaler = scaler

        if path.find('solomon') > 0:
            self.filepath = f'{path}/{problem_name}.txt'
            self.nodes_location, self.dist_mat = self.solomon_info()
            self.s, self.t, self.N, self.A1, self.A2, self.A, self.L, self.tv, self.td, self.stv, self.std, self.b, self.M = self.make_parameter()

        if path.find('Augerat') > 0:
            self.filepath = f'{path}/{problem_name}.vrp'
            self.nodes_location, self.dist_mat = self.augerat_info()
            self.s, self.t, self.N, self.A1, self.A2, self.A, self.L, self.tv, self.td, self.stv, self.std, self.b, self.M = self.make_parameter()

        self.info = {
            'num_nodes': len(self.N),
            'problem_name': self.problem_name,
            'alpha': self.alpha,
            'B' : self.B, 
            'truck_service_time' : self.truck_service_time,
            'drone_service_time' : self.drone_service_time,
            'drone_multiplier': self.drone_multiplier,
            'scaler': self.scaler
        }


    def augerat_info(self):

        f = open(self.filepath, 'r')
        
        node_list = []
        demand = {}

        find_node = False

        while True:
            line = f.readline()
            if not line:
                break
            if line.strip() == 'DEMAND_SECTION':
                find_node = False
            if find_node:
                node_list.append(re.split("\W+", line.strip()))
            if line.strip() == 'NODE_COORD_SECTION':
                find_node = True

        f.close()

        nodes_location = {}
        for i, x, y in node_list:
            nodes_location[int(i)-1] = (int(x), int(y))
        nodes_location[len(node_list)] = nodes_location[0]

        dist_mat = np.zeros((len(nodes_location), len(nodes_location)))

        for i in nodes_location:
            for j in nodes_location:
                dist_mat[i][j] = (((nodes_location[i][0] - nodes_location[j][0]) ** 2) + (
                        (nodes_location[i][1] - nodes_location[j][1]) ** 2)) ** 0.5
        

        return nodes_location, dist_mat

    def solomon_info(self):
        filepath = self.filepath

        lineset = []
        nodes_location = {}

        with open(filepath) as fp:
            line = fp.readline()
            while line:
                lineset.append(line.strip())
                line = fp.readline()
                
        n = int(filepath.split('/')[-2])

        for cust_no in lineset[9:9+n+1]:
            element = re.split("\W+", cust_no)
            element = list(map(int, element))
            nodes_location[element[0]] = (element[1], element[2])

        nodes_location[len(nodes_location)] = nodes_location[0]
        dist_mat = np.zeros((len(nodes_location), len(nodes_location)))

        for i in nodes_location:
            for j in nodes_location:
                dist_mat[i][j] = (((nodes_location[i][0] - nodes_location[j][0]) ** 2) + (
                        (nodes_location[i][1] - nodes_location[j][1]) ** 2)) ** 0.5
        

        return nodes_location, dist_mat

    def make_parameter(self):
        nodes_location = self.nodes_location
        dist_mat = self.dist_mat
        scaler = self.scaler
        drone_multiplier = self.drone_multiplier

        s = 0
        
        t = len(nodes_location) - 1

        # Set of customers
        N = list(range(1, t))
        depot_N = [s,t]

        # set of arcs
        A1 = [s] + N
        A2 = N + [t]
        A = [(i, j) for i in A1 for j in A2 if (i != j) and (i != s or j != t)]

        # Number_of_drones
        L = [int(idx) for idx in range(len(self.B))]

        # Travel_time_of_vehicle, scaler를 곱하고 정수로 반올림
        tv = np.around(dist_mat * scaler, decimals=0)

        # Travel_time_of_drone, scaler와 alpha를 곱하고 정수로 반올림
        td = {
            l: np.around(dist_mat * scaler * drone_multiplier[l], decimals=0)
            for l in L
        }

        # Service_time_of_vehicle, scaler를 곱
        stv = {
            i: self.truck_service_time * scaler
            for i in range(len(nodes_location))
        }
        
        for i in depot_N :
            stv[i] = 0

        # Service_time_of_drone
        std = {
            l: {
                i: self.drone_service_time * scaler 
                for i in range(len(nodes_location))
            }
            for l in L
        }
        
        for l in L:
            for i in depot_N:
                std[l][i] = 0

        # Required_battery = 2 * td[l][i][j] + std[l][j]
        b = {
            l: 2 * td[l].copy() + np.array(list(std[l].values()))
            for l in L
        }

        # Big M
        M = 100000

        return s, t, N, A1, A2, A, L, tv, td, stv, std, b, M
        
class Prob2 :
    def __init__(self, path, problem_name, alpha, scaler, drone_battery=1800, truck_speed=8, drone_speed=10, truck_service_time=60, drone_service_time=30, heavy_demand=[]):
        self.path = path
        self.problem_name = problem_name
        self.alpha = alpha
        self.B = [drone_battery/ i for i in alpha]
        self.truck_speed = truck_speed
        self.drone_speed = [drone_speed*i for i in alpha]
        self.truck_service_time = truck_service_time
        self.drone_service_time = drone_service_time
        self.scaler = scaler
        self.heavy_demand = heavy_demand
        
        self.read_vrp()
        self.s, self.t, self.N, self.A1, self.A2, self.A, self.L, self.tv, self.td, self.stv, self.std, self.b, self.M = self.make_parameter()
        
        self.info = {
            'num_nodes': len(self.N),
            'problem_name': self.problem_name,
            'alpha' : self.alpha,
            'B': self.B,
            'truck_speed': self.truck_speed,
            'drone_speed' : self.drone_speed,
            'truck_service_time' : self.truck_service_time,
            'drone_service_time' : self.drone_service_time,
            'scaler': self.scaler,
            'heavy_demand' : self.heavy_demand
        }
        
    def read_vrp(self):
        datafile = f'{self.path}/{self.problem_name}.vrp' 
        
        f = open(datafile, 'r')

        truck_dist = []
        coord = {}

        find_node = False
        find_dist = False

        while True:
            line = f.readline()
            if not line:
                break

            if line.strip() == 'DEPOT_SECTION':
                find_dist = False

            if find_dist:      
                dist = re.split("\W+", line.strip())
                dist = [int(i) for i in dist]
                dist.append(dist[0])
                truck_dist.append(dist)

            if line.strip() == 'EDGE_WEIGHT_SECTION':
                find_node = False
                find_dist = True

            if find_node:
                node_coord = re.split(" ", line.strip())
                coord[int(node_coord[0])] = (float(node_coord[1]), float(node_coord[2]))

            if line.strip() == 'NODE_COORD_SECTION':
                find_node = True

        f.close()

        truck_dist = np.array(truck_dist)
        coord[len(truck_dist)] = coord[0]

        drone_dist = np.zeros((len(coord), len(coord)))

        for i in coord:
            for j in coord:
                drone_dist[i][j] = get_distance_m(coord[i], coord[j])
                
        
        self.truck_dist = truck_dist
        self.coord = coord
        self.drone_dist = drone_dist
        
    def make_parameter(self):
        
        s = 0 

        t = len(self.coord)-1

        depot_N = [s,t]

        N = list(range(1, len(self.coord)-1))

        A1 = [s] + N
        A2 = N + [t]
        A = [(i,j) for i in A1 for j in A2 if (i != j) and (i != s or j != t)]

        L = [int(idx) for idx in range(len(self.B))]

        tv = np.around(self.truck_dist/self.truck_speed * self.scaler, decimals=0)

        td = {
            l: np.around(self.drone_dist/self.drone_speed[l] * self.scaler, decimals=0)
            for l in L
        }


        # Service_time_of_vehicle
        stv = {
            i: self.truck_service_time * self.scaler
            for i in self.coord

        }

        for i in depot_N:
            stv[i] = 0

        # Service_time_of_drone
        std = {
            l: {
                i: self.drone_service_time * self.scaler 
                for i in self.coord
            }
            for l in L
        }

        for l in L:
            for i in depot_N:
                std[l][i] = 0
                
        for l, a in enumerate(self.alpha):
            if a > 1.000001 :
                for i in self.heavy_demand:
                    std[l][i] = 100000 

        # Required_battery = 2 * td[l][i][j] + std[l][j]

        b = {
            l: 2 * td[l].copy() + np.array(list(std[l].values()))
            for l in L
        }

        # Big M
        M = 100000
        
        return s, t, N, A1, A2, A, L, tv, td, stv, std, b, M
    

class MIP:
    def __init__(self, prob):
        self.prob = prob
        self.formulation()

    def formulation(self):
        B = [battery_capacity * self.prob.scaler for battery_capacity in self.prob.B]

        s = self.prob.s
        t = self.prob.t

        N = self.prob.N
        A1 = self.prob.A1
        A2 = self.prob.A2
        A = self.prob.A
        L = self.prob.L

        tv = self.prob.tv
        td = self.prob.td

        stv = self.prob.stv
        std = self.prob.std

        b = self.prob.b
        M = self.prob.M

        # Formulation P
        cpx = cplex.Cplex()

        cpx.variables.add(
            obj=[(1)] * len(A1),
            names=[f'w_{i}'
                   for i in A1]
        )

        cpx.variables.add(
            obj=[(tv[i][j] + stv[j])
                 for (i, j) in A],
            types=["B"] * len(A),
            names=[f'x_{i}_{j}'
                   for (i, j) in A]
        )

        cpx.variables.add(

            types=["B"] * len(L) * len(A),
            names=[f'h_{l}_{i}_{j}'
                   for l in L
                   for (i, j) in A]
        )

        cpx.variables.add(

            names=[f'v_{i}'
                   for i in [s] + N + [t]]
        )

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'x_{s}_{j}' for j in N], ([1] * (len(N)))
                )
            ],
            senses=['E'],
            rhs=[1],
            names=['start'])

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'x_{i}_{t}' for i in N], ([1] * (len(N)))
                )
            ],
            senses=['E'],
            rhs=[1],
            names=['end'])

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'x_{i}_{j}' for j in A2 if j != i] + [f'x_{j}_{i}' for j in A1 if j != i],
                    ([1] * (len(A2) - 1)) + ([-1] * (len(A1) - 1))
                )
                for i in N
            ],
            senses=['E'] * len(N),
            rhs=[0] * len(N),
            names=[f'flow_balance_{i}' for i in N]
        )

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'v_{i}'] + [f'v_{j}'] + [f'x_{i}_{j}'], ([1]) + ([-1]) + ([M])
                )
                for (i, j) in A
            ],
            senses=['L'] * len(A),
            rhs=[M - 1] * len(A),
            names=[f'MTZ_{i}_{j}' for (i, j) in A])

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'x_{i}_{j}' for i in A1 if i != j] + [f'h_{l}_{i}_{j}' for l in L for i in A1 if i != j],
                    ([1 for i in A1 if i != j]) + ([1 for l in L for i in A1 if i != j])
                )
                for j in N
            ],
            senses=['E'] * len(N),
            rhs=[1] * len(N),
            names=[f'node_served_{j}' for j in N])

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'x_{i}_{j}' for j in A2 if j != i and (i != s or j != t)] + [f'h_{l}_{i}_{j}' for l in L for j in N if j != i],
                    ([M for j in A2 if j != i and (i != s or j != t)]) + ([-1 for l in L for j in N if j != i])
                )
                for i in A1
            ],
            senses=['G'] * (len(A1)),
            rhs=[0] * (len(A1)),
            names=[f'vehicle_stop_{i}' for i in A1])

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'h_{l}_{i}_{j}' for j in N for i in A1 if j != i], ([b[l][i][j] for j in N for i in A1 if j != i])
                )
                for l in L

            ],
            senses=['L'] * len(L),
            rhs=[B[l] for l in L],
            names=[f'battery_consumption_{l}' for l in L])

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'h_{l}_{i}_{j}' for j in N if j != i] + [f'w_{i}'],
                    ([-((2 * td[l][i][j]) + std[l][j]) for j in N if j != i]) + ([1])
                )
                for l in L
                for i in A1
            ],
            senses=['G'] * len(L) * len(A1),
            rhs=[0] * len(L) * len(A1),
            names=[f'waiting_time_{i}_{l}' for l in L for i in A1])

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'v_{0}'], ([1])
                )
            ],
            senses=['E'],
            rhs=[0],
            names=[f'MTZ_0'])

        self.cpx = cpx

    def fix_route(self, num):
        cpx = self.cpx
        A1 = self.prob.A1
        A2 = self.prob.A2

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'x_{i}_{j}' for i in A1 for j in A2 if i != j and (i != self.prob.s or j != self.prob.t)],
                    ([1 for i in A1 for j in A2 if i != j and (i != self.prob.s or j != self.prob.t)])
                )
            ],
            senses=['E'],
            rhs=[num],
            names=[f'fix_route'])

        self.cpx = cpx
        self.fix_num = num

    def make_routes(self):
        cpx = self.cpx
        vehicle_nodes = {}

        for (i, j) in self.prob.A:
            val = cpx.solution.get_values(f'x_{i}_{j}')
            if val > 0.001:
                vehicle_nodes[i] = j

        drone_nodes = {}

        for l in self.prob.L:
            drone_nodes[l] = []
            for (i, j) in self.prob.A:
                val = cpx.solution.get_values(f'h_{l}_{i}_{j}')
                if val > 0.001:
                    drone_nodes[l].append([i, j])

        return vehicle_nodes, drone_nodes

    def solve(self, cpu=2, time_limit=3600):
        cpx = self.cpx

        start_time = time.time()
        cpx.parameters.threads.set(cpu)
        cpx.parameters.timelimit.set(time_limit)

        cpx.set_log_stream(None)
        cpx.set_results_stream(None)

        try:
            cpx.solve() # 메모리 부족일 경우가 있으므로 
        except:
            pass

        num_bnb_nodes = cplex._internal._procedural.getnodecnt(cpx._env._e, cpx._lp)
        num_gap = cplex._internal._procedural.getmiprelgap(cpx._env._e, cpx._lp)
        solve_time = time.time() - start_time

        self.num_bnb_nodes = num_bnb_nodes
        self.num_gap = num_gap
        self.solve_time = solve_time

        self.vehicle_nodes, self.drone_nodes = self.make_routes()

        results = {
            'status': cpx.solution.get_status_string(),
            'obj': cpx.solution.get_objective_value() / self.prob.scaler,
            'time': solve_time,
            'num_bnb_nodes': num_bnb_nodes,
            'num_gap': num_gap,
            'cpu': cpu,
            'time_limit': time_limit,
            'vehicle_nodes': self.vehicle_nodes,
            'drone_nodes': self.drone_nodes
        }

        return results

    def solve_relaxation(self, cpu=2, time_limit=3600):
        cpx = self.cpx

        start_time = time.time()
        cpx.parameters.threads.set(cpu)
        cpx.parameters.timelimit.set(time_limit)
        cpx.parameters.mip.limits.nodes.set(0)

        cpx.set_log_stream(None)
        cpx.set_results_stream(None)
        cpx.set_warning_stream(None)

        cb = cpx.register_callback(MySolve)
        cb.all_names = cpx.variables.get_names()
        cpx.solve()

        solve_time = time.time() - start_time
        self.solve_time = solve_time

        results = {
            'status': cpx.solution.get_status_string(),
            'relaxation_obj': cpx.solution.MIP.get_best_objective() / self.prob.scaler,
            'time': solve_time,
        }

        return results

class BendersMasterProb:
    def __init__(self, prob, incumbent, lb, ub, w_star_dic, use_ucb=True):
        self.prob = prob
        self.w_star_dic = {
            i: waiting_time * prob.scaler for i, waiting_time in w_star_dic.items()
        }
        self.opt_cut = 0
        self.inf_cut = 0

        upper_cutoff = incumbent * prob.scaler

        B = [battery_capacity * prob.scaler for battery_capacity in prob.B]

        s = prob.s
        t = prob.t

        N = prob.N
        A = prob.A
        A1 = prob.A1
        A2 = prob.A2
        L = prob.L

        tv = prob.tv
        stv = prob.stv
        td = prob.td
        std = prob.std
        b = prob.b

        cpx = cplex.Cplex()
        cpx.objective.set_sense(cpx.objective.sense.minimize)
        cpx.parameters.mip.tolerances.uppercutoff.set(upper_cutoff)
        # cpx.parameters.mip.tolerances.mipgap.set(0.01)
        # cpx.parameters.mip.tolerances.integrality.set(0.01)

        cpx.set_log_stream(None)
        cpx.set_warning_stream(None)
        cpx.set_results_stream(None)

        cpx.variables.add(
            obj=[tv[i][j] + stv[j] for (i, j) in A],
            types=['B'] * len(A),
            names=[f'x_{i}_{j}' for (i, j) in A]
        )

        cpx.variables.add(
            types=['B'] * len(N),
            names=[f'z_{i}' for i in N]
        )

        self.var_z_index = list(range(len(prob.A), len(prob.A) + len(prob.N)))

        cpx.variables.add(
            names=[f'h_{l}_{i}_{j}' for l in L for (i, j) in A]
        )

        cpx.variables.add(
            obj=[1],
            names=[f'W']
        )

        self.var_W_index = len(prob.A) + len(prob.N) + (len(prob.L) * len(prob.A))

        cpx.variables.add(
            names=[f'w_{i}' for i in A1]
        )

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'x_{s}_{j}' for j in N], [1] * len(N)
                )
            ],
            senses=['E'],
            rhs=[1],
            names=['start']
        )

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'x_{i}_{t}' for i in N], [1] * len(N)
                )
            ],
            senses=['E'],
            rhs=[1],
            names=['end']
        )

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'x_{i}_{j}' for j in A2 if j != i] + [f'x_{j}_{i}' for j in A1 if j != i],
                    ([1] * (len(A2) - 1)) + ([-1] * (len(A1) - 1))
                )
                for i in N
            ],
            senses=['E'] * len(N),
            rhs=[0] * len(N),
            names=[f'flow_balance_{i}' for i in N]
        )

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'x_{i}_{j}' for j in A2 if j != i] + [f'z_{i}'], [1] * (len(A2) - 1) + [-1]
                )
                for i in N
            ],
            senses=['E'] * len(N),
            rhs=[0] * len(N),
            names=[f'vehicle_visit_{i}' for i in N]
        )

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'h_{l}_{i}_{j}' for l in L for i in A1 if j != i] + [f'z_{j}'], [1] * len(L) * (len(A1) - 1) + [1]
                )
                for j in N
            ],
            senses=['E'] * len(N),
            rhs=[1] * len(N),
            names=[f'node_visit_{j}' for j in N]
        )

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'z_{i}'] + [f'h_{l}_{i}_{j}'], [1] + [-1]
                )
                for i in N
                for j in A2
                for l in L
                if i != j
            ],
            senses=['G' for i in N for j in A2 for l in L if i != j],
            rhs=[0 for i in N for j in A2 for l in L if i != j],
            names=[f'vehicle_stop_drone_avail_{l}_{i}_{j}' for i in N for j in A2 for l in L if i != j]
        )

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'h_{l}_{i}_{j}' for i in A1 for j in N if j != i], [b[l][i][j] for i in A1 for j in N if j != i]
                )
                for l in L
            ],
            senses=['L'] * len(L),
            rhs=[B[l] for l in L],
            names=[f'drone_battery{l}' for l in L]
        )

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'W'] + [f'w_{i}' for i in A1], [1] + [-1] * len(A1)
                )
            ],
            senses=['G'],
            rhs=[0],
            names=[f'total_waiting_time']
        )

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'w_{i}'] + [f'h_{l}_{i}_{j}' for j in N if j != i],
                    [1] + [-((2 * td[l][i][j]) + std[l][j]) for j in N if j != i]
                )
                for l in L
                for i in A1
            ],
            senses=['G'] * len(L) * len(A1),
            rhs=[0] * len(L) * len(A1),
            names=[f'waiting_time_{l}_{i}' for l in L for i in A1]

        )

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'w_{i}'] + [f'h_{l}_{i}_{j}' for l in L],
                    ([1]) + ([-min([((2 * td[l][i][j]) + std[l][j]) for l in L]) for l in L])
                )
                for (i, j) in A
            ],
            senses=['G'] * len(A),
            rhs=[0] * len(A),
            names=[f'waiting_time2_{i}_{j}' for i, j in A]
        )

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'x_{i}_{j}' for j in A2 for i in N if j != i], [1] * (len(A2) - 1) * len(N)
                )
            ],
            senses=['G'],
            rhs=[lb],
            names=[f'num_visits_lb']
        )

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'x_{i}_{j}' for j in A2 for i in N if j != i], [1] * (len(A2) - 1) * len(N)
                )
            ],
            senses=['L'],
            rhs=[ub],
            names=[f'num_visits_ub']
        )

        self.cpx = cpx

        bmp_sol = {i: 0 for i in prob.N}

        ucb = cpx.register_callback(CutGenCallback)
        ucb.times_called = 0
        ucb.generated_cuts = []
        ucb.prob = prob
        ucb.var_x_index = list(range(len(prob.A)))
        ucb.solve_time_gcs = 0
        ucb.times_gcs = 0

        if not use_ucb:
            cpx.unregister_callback(CutGenCallback)

        lcb = cpx.register_callback(LazyCutGenCallback)
        lcb.times_called = 0
        lcb.times_inf = 0
        lcb.times_opt = 0
        lcb.times_gcs = 0
        lcb.generated_cuts = []
        lcb.generated_Benders_cuts = []

        lcb.prob = prob
        lcb.var_x_index = list(range(len(prob.A)))
        lcb.var_z_index = list(range(len(prob.A), len(prob.A) + len(prob.N)))
        lcb.var_W_index = len(prob.A) + len(prob.N) + (len(prob.L) * len(prob.A))
        lcb.bsp = BendersSubProb(prob, bmp_sol)
        lcb.solve_time_gcs = 0
        lcb.solve_time_sp = 0
        lcb.gcs_obj = 0
        lcb.w_star_dic = self.w_star_dic

        self.ucb = ucb
        self.lcb = lcb

    def make_routes(self):
        cpx = self.cpx
        vehicle_nodes = {}

        for (i, j) in self.prob.A:
            val = cpx.solution.get_values(f'x_{i}_{j}')
            if val > 0.001:
                vehicle_nodes[i] = j

        bsp = BendersSubProb(self.prob, self.bmp_sol)
        bsp.solve()

        drone_nodes = {}

        for l in self.prob.L:
            drone_nodes[l] = []
            for (i, j) in self.prob.A:
                val = bsp.cpx.solution.get_values(f'h_{l}_{i}_{j}')
                if val > 0.001:
                    drone_nodes[l].append([i, j])

        waiting_time = {}

        for i in self.prob.A1:
            val = bsp.cpx.solution.get_values(f'w_{i}')
            waiting_time[i] = val / self.prob.scaler

        return vehicle_nodes, drone_nodes, waiting_time

    def solve(self, cpu=2, time_limit=3600):
        cpx = self.cpx

        start_time = time.time()
        cpx.parameters.threads.set(cpu)
        cpx.parameters.timelimit.set(time_limit)

        cpx.set_log_stream(None)
        cpx.set_results_stream(None)

        cpx.solve()

        num_bnb_nodes = cplex._internal._procedural.getnodecnt(cpx._env._e, cpx._lp)
        num_gap = cplex._internal._procedural.getmiprelgap(cpx._env._e, cpx._lp)
        solve_time = time.time() - start_time

        self.num_bnb_nodes = num_bnb_nodes
        self.num_gap = num_gap
        self.solve_time = solve_time

        self.bmp_sol = {
            i: self.cpx.solution.get_values(self.var_z_index)[i - 1]
            for i in self.prob.N
        }

        self.vehicle_nodes, self.drone_nodes, self.waiting_time = self.make_routes()

        results = {
            'status': cpx.solution.get_status_string(),
            'obj': cpx.solution.get_objective_value() / self.prob.scaler,
            'tot_solve_time': solve_time,
            'solve_time_gcs': self.lcb.solve_time_gcs + self.ucb.solve_time_gcs,
            'solve_time_sp': self.lcb.solve_time_sp,
            'lcb_times_inf': self.lcb.times_inf,
            'lcb_times_opt': self.lcb.times_opt,
            'lcb_times_gcs': self.lcb.times_gcs,
            'ucb_times_gcs': self.ucb.times_gcs,
            'num_bnb_nodes': num_bnb_nodes,
            'num_gap': num_gap,
            'cpu': cpu,
            'time_limit': time_limit,
            'vehicle_nodes': self.vehicle_nodes,
            'drone_nodes': self.drone_nodes,
            'bsp_waiting_time': self.waiting_time,
            'bmp_waiting_time': cpx.solution.get_values('W') / self.prob.scaler
        }

        return results

class BendersSubProb:
    def __init__(self, prob, bmp_sol):
        self.prob = prob
        self.bmp_sol = bmp_sol

        B = [battery_capacity * prob.scaler for battery_capacity in prob.B]

        N = prob.N
        A = prob.A
        A1 = prob.A1
        L = prob.L

        td = prob.td
        std = prob.std
        b = prob.b

        cpx = cplex.Cplex()
        cpx.objective.set_sense(cpx.objective.sense.minimize)
        cpx.parameters.mip.tolerances.mipgap.set(0.001)
        cpx.parameters.mip.tolerances.absmipgap.set(1)
#         cpx.parameters.preprocessing.symmetry.set(3)
        
        # cpx.parameters.mip.tolerances.integrality.set(0.01)
        # cpx.parameters.mip.tolerances.uppercutoff.set(incumbent)

        cpx.set_log_stream(None)
        cpx.set_warning_stream(None)
        cpx.set_error_stream(None)
        cpx.set_results_stream(None)

        cpx.variables.add(

            types=['B'] * len(L) * len(A),
            names=[f'h_{l}_{i}_{j}' for l in L for (i, j) in A]
        )

        cpx.variables.add(

            obj=[1] * len(A1),
            names=[f'w_{i}' for i in A1]
        )

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'h_{l}_{i}_{j}' for i in A1 for j in N if i != j], [b[l][i][j] for i in A1 for j in N if i != j]
                )
                for l in L
            ],
            senses=['L'] * len(L),
            rhs=[B[l] for l in L],
            names=[f'drone_battery{l}' for l in L]
        )

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'w_{i}'] + [f'h_{l}_{i}_{j}' for j in N if j != i],
                    [1] + [-((2 * td[l][i][j]) + std[l][j]) for j in N if j != i]
                )
                for i in A1
                for l in L
            ],
            senses=['G'] * len(A1) * len(L),
            rhs=[0] * len(A1) * len(L),
            names=[f'waiting_time_{l}_{i}' for i in A1 for l in L]
        )

        const_node_visit_start = cpx.linear_constraints.get_num()

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'h_{l}_{i}_{j}' for l in L for i in A1 if j != i], [1] * len(L) * (len(A1) - 1)
                )
                for j in N
            ],
            senses=['E'] * len(N),
            rhs=[1 - bmp_sol[j] for j in N],
            names=[f'node_visit_{j}' for j in N]
        )

        const_node_visit_end = cpx.linear_constraints.get_num() - 1

        const_node_visit_2_start = cpx.linear_constraints.get_num()

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'h_{l}_{i}_{j}' for l in L for j in N if j != i], [1] * len(L) * (len(N) - 1)
                )
                for i in N
            ],
            senses=['L'] * len(N),
            rhs=[len(N) * bmp_sol[i] for i in N],
            names=[f'node_visit_2_{i}' for i in N]
        )

        const_node_visit_2_end = cpx.linear_constraints.get_num() - 1

        self.cpx = cpx

        self.var_h_index = list(range(0, len(L) * len(A)))
        self.var_w_index = list(range(len(L) * len(A), len(L) * len(A) + len(A1)))

        idx = 0
        h_lij_idx_dic = {}
        for l in L:
            for (i, j) in A:
                h_lij_idx_dic[(l, i, j)] = idx
                idx += 1

        self.h_lij_idx_dic = h_lij_idx_dic

        node_visit_list = list(range(const_node_visit_start, const_node_visit_end + 1))
        node_visit_2_list = list(range(const_node_visit_2_start, const_node_visit_2_end + 1))

        self.node_visit_const_idx = {
            i: node_visit_list[i - 1]
            for i in prob.N
        }

        self.node_visit_2_const_idx = {
            i: node_visit_2_list[i - 1]
            for i in prob.N
        }

    def solve(self, cpu=2, time_limit=3600, upper_cutoff=cplex.infinity):
        cpx = self.cpx

        start_time = time.time()
        cpx.parameters.threads.set(cpu)
        cpx.parameters.timelimit.set(time_limit)

        cpx.parameters.mip.tolerances.uppercutoff.set(upper_cutoff)

        cpx.set_log_stream(None)
        cpx.set_results_stream(None)

        cpx.solve()

        if cpx.solution.get_status() == 103:  # infeasible
            return None

        solve_time = time.time() - start_time
        self.solve_time = solve_time
        self.w_sol = self.cpx.solution.get_values(self.var_w_index)

        results = {
            'status': cpx.solution.get_status_string(),
            'obj': cpx.solution.get_objective_value() / self.prob.scaler,
            'time': solve_time
        }

        return results

class P_median:
    def __init__(self, prob, p_median_num):
        self.prob = prob
        self.p_median_num = p_median_num
        self.formulation()

    def formulation(self):
        B = [battery_capacity * self.prob.scaler for battery_capacity in self.prob.B]

        N = self.prob.N
        A1 = self.prob.A1
        A = self.prob.A
        L = self.prob.L

        td = self.prob.td
        std = self.prob.std
        b = self.prob.b

        p_median_num = self.p_median_num

        # Formulation P_median
        cpx = cplex.Cplex()

        self.var_z_index_list = list(range(len(A1)))

        cpx.variables.add(

            types=["B"] * len(A1),
            names=[f'z_{i}'
                   for i in A1]
        )

        cpx.variables.add(
            obj=[(1)
                 for i in A1],
            names=[f'w_{i}'
                   for i in A1]
        )

        cpx.variables.add(

            types=["B"] * len(A),
            names=[f'y_{i}_{j}'
                   for (i, j) in A]
        )

        cpx.variables.add(

            types=["B"] * len(A) * len(L),
            names=[f'h_{l}_{i}_{j}'
                   for l in L
                   for (i, j) in A
                   ]
        )

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'h_{l}_{i}_{j}' for j in N for i in A1 if j != i],
                    ([b[l][i][j] for j in N for i in A1 if j != i])
                )
                for l in L
            ],
            senses=['L'] * len(L),
            rhs=[B[l] for l in L],
            names=[f'battery_consumption_{l}' for l in L])

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'w_{i}'] + [f'h_{l}_{i}_{j}' for j in N if j != i],
                    ([1]) + ([-((2 * td[l][i][j]) + std[l][j]) for j in N if j != i])
                )
                for i in A1
                for l in L
            ],
            senses=['G'] * len(A1) * len(L),
            rhs=[0] * len(A1) * len(L),
            names=[f'waiting_time_{l}_{i}' for i in A1 for l in L])

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'y_{i}_{j}'] + [f'h_{l}_{i}_{j}' for l in L], ([1]) + ([-1] * (len(L)))
                )
                for (i, j) in A
            ],
            senses=['E'] * len(A),
            rhs=[0] * len(A),
            names=[f'drone_delivery_{i}_{j}' for (i, j) in A])

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'z_{i}'] + [f'y_{i}_{j}'], ([1]) + ([-1])
                )
                for i in A1
                for j in N
                if i != j
            ],
            senses=['G' for i in A1 for j in N if i != j],
            rhs=[0 for i in A1 for j in N if i != j],
            names=[f'drone_{i}_{j}' for i in A1 for j in N if i != j])

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'z_{i}' for i in N], ([1] * len(N))
                )
            ],
            senses=['E'],
            rhs=[p_median_num],
            names=[f'the_number_of_P'])

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'y_{i}_{j}' for i in A1 if i != j] + [f'z_{j}'], ([1 for i in A1 if i != j]) + ([1])
                )
                for j in N
            ],
            senses=['G'] * len(N),
            rhs=[1] * len(N),
            names=[f'visit_nodes_{j}' for j in N])

        self.cpx = cpx

    def make_drone_routes(self):
        prob = self.prob
        cpx = self.cpx
        drone_route = {}

        for l in prob.L:
            drone_route[l] = []
            for (i, j) in prob.A:
                val = cpx.solution.get_values(f'h_{l}_{i}_{j}')
                if val > 0.001:
                    drone_route[l].append((i, j))

        return drone_route

    def solve(self, cpu=2, time_limit=3600):
        cpx = self.cpx

        start_time = time.time()
        cpx.parameters.threads.set(cpu)
        cpx.parameters.timelimit.set(time_limit)

        cpx.set_log_stream(None)
        cpx.set_results_stream(None)
        cpx.set_warning_stream(None)
        cpx.set_error_stream(None)
        cpx.solve()

        num_bnb_nodes = cplex._internal._procedural.getnodecnt(cpx._env._e, cpx._lp)
        num_gap = cplex._internal._procedural.getmiprelgap(cpx._env._e, cpx._lp)
        solve_time = time.time() - start_time

        self.num_bnb_nodes = num_bnb_nodes
        self.num_gap = num_gap
        self.solve_time = solve_time

        val_z_sol = self.cpx.solution.get_values(self.var_z_index_list)
        positive_z_list = [i for i in range(len(val_z_sol)) if val_z_sol[i] > 0.001]
        zero_z_list = [i for i in range(len(val_z_sol)) if val_z_sol[i] < 0.001]

        self.positive_z_list = positive_z_list
        self.zero_z_list = zero_z_list

        drone_route = self.make_drone_routes()

        results = {
            'status': cpx.solution.get_status_string(),
            'obj': cpx.solution.get_objective_value() / self.prob.scaler,
            'time': solve_time,
            'num_bnb_nodes': num_bnb_nodes,
            'num_gap': num_gap,
            'cpu': cpu,
            'time_limit': time_limit,
            'drone_route': drone_route
        }

        return results

    def solve_root(self, cpu=2, time_limit=3600):
        cpx = self.cpx

        start_time = time.time()
        cpx.parameters.threads.set(cpu)
        cpx.parameters.timelimit.set(time_limit)
        cpx.parameters.mip.limits.nodes.set(0)

        cpx.set_log_stream(None)
        cpx.set_results_stream(None)
        cpx.set_warning_stream(None)
        cpx.set_error_stream(None)
        cpx.solve()

        solve_time = time.time() - start_time

        self.solve_time = solve_time

        results = {
            'status': cpx.solution.get_status_string(),
            'obj': cpx.solution.MIP.get_best_objective() / self.prob.scaler,
            'time': solve_time
        }

        return results

class TSP:
    def __init__(self, prob, positive_s_list, zero_s_list):
        s = prob.s
        t = prob.t
        N = prob.N
        A = prob.A
        A1 = prob.A1
        A2 = prob.A2
        tv = prob.tv
        stv = prob.stv

        cpx = cplex.Cplex()
        cpx.objective.set_sense(cpx.objective.sense.minimize)

        cpx.set_log_stream(None)
        cpx.set_warning_stream(None)
        cpx.set_results_stream(None)

        cpx.variables.add(
            obj=[(tv[i][j] + stv[j]) for (i, j) in A],
            types=['B'] * len(A),
            names=[f'x_{i}_{j}' for (i, j) in A]
        )

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'x_{s}_{j}' for j in N], [1] * len(N)
                )
            ],
            senses=['E'],
            rhs=[1],
            names=['start']
        )

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'x_{i}_{t}' for i in N], [1] * len(N)
                )
            ],
            senses=['E'],
            rhs=[1],
            names=['end']
        )

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'x_{i}_{j}' for j in A2 if j != i] + [f'x_{j}_{i}' for j in A1 if j != i],
                    ([1] * (len(A2) - 1)) + ([-1] * (len(A1) - 1))
                )
                for i in N
            ],
            senses=['E'] * len(N),
            rhs=[0] * len(N),
            names=[f'flow_balance_{i}' for i in N]
        )

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'x_{i}_{j}' for j in A2 if j != i], [1] * (len(A2) - 1)
                )
                for i in positive_s_list
            ],
            senses=['E'] * len(positive_s_list),
            rhs=[1] * len(positive_s_list),
            names=[f'single_visit_{i}' for i in positive_s_list]
        )

        cpx.linear_constraints.add(
            lin_expr=[
                cplex.SparsePair(
                    [f'x_{i}_{j}' for j in A2 if j != i], [1] * (len(A2) - 1)
                )
                for i in zero_s_list
            ],
            senses=['E'] * len(zero_s_list),
            rhs=[0] * len(zero_s_list),
            names=[f'single_visit_2_{i}' for i in zero_s_list]
        )

        self.var_x_idx_list = [i for i in range(len(A))]

        ucb = cpx.register_callback(CutGenCallback)
        ucb.times_called = 0
        ucb.generated_cuts = []
        ucb.prob = prob
        ucb.var_x_index = self.var_x_idx_list

        lcb = cpx.register_callback(LazyCutGenCallbackTSP)
        lcb.times_called = 0
        lcb.generated_cuts = []
        lcb.prob = prob
        lcb.var_x_index = self.var_x_idx_list

        self.prob = prob
        self.cpx = cpx
        self.var_names = [f'x_{i}_{j}' for (i, j) in A]

        self.ucb = ucb
        self.lcb = lcb

        self.var_x_idx_dic = {
            idx: (i, j)
            for idx, (i, j) in enumerate(A)
        }

    def sol_dic_to_path(self):
        k = self.prob.s
        path = []
        while k != self.prob.t:
            path.append(k)
            for (i, j) in self.positive_x_sol:
                if i == k:
                    k = j
                    break

        path.append(self.prob.t)
        return path

    def solve(self, cpu=2, time_limit=3600):

        cpx = self.cpx
        cpx.parameters.threads.set(cpu)
        cpx.parameters.timelimit.set(time_limit)

        self.ucb.times_called = 0
        self.ucb.generated_cuts = []
        self.ucb.solve_time_gcs = 0
        self.ucb.times_gcs = 0

        self.lcb.times_called = 0
        self.lcb.generated_cuts = []
        self.lcb.solve_time_gcs = 0
        self.lcb.times_gcs = 0

        cpx.solve()

        var_x_sol_list = cpx.solution.get_values(self.var_x_idx_list)
        self.positive_x_sol = [self.var_x_idx_dic[idx] for idx, val in enumerate(var_x_sol_list) if val > 0.001]

        path = self.sol_dic_to_path()
        vehicle_cost = sum([self.prob.tv[i][j] + self.prob.stv[j] for i, j in zip(path[:-1], path[1:])])

        results = {
            'status': cpx.solution.get_status_string(),
            'obj': cpx.solution.get_objective_value() / self.prob.scaler,
            'vehicle_path': path,
            'vehicle_cost': vehicle_cost
        }

        return results

class Heuristic_Solver:
    def __init__(self, prob):
        self.prob = prob

    def solve(self, cpu_num=2, time_limit_value=3600):
        prob = self.prob
        check_ub_lb = []
        start_time = time.time()

        # 노드의 개수가 증가할수록 오래걸리기 때문에 풀리는 가장 작은 노드 한개만 p-median을 푼다.
        for p_median_num in prob.N:
            p_median_num = int(p_median_num)

            try:
                p_median = P_median(prob, p_median_num)
                p_median_results = p_median.solve(cpu=cpu_num, time_limit=300)
                
                tsp = TSP(prob, p_median.positive_z_list[1:], p_median.zero_z_list)
                tsp_results = tsp.solve(cpu=cpu_num, time_limit=300)
                break

            except:
                pass

        lb_num_r_visits = p_median_num

        fix_num_list = sorted(prob.A2[1:], reverse=True)

        mip_root_node_obj_dic = {}
        new_list = prob.A2[p_median_num:]
        
        iter_num = 0
        fix_num = new_list[round(len(new_list)/3)] # truck length 
        
        while True :
            mip_relaxation = None
            mip_relaxation = MIP(prob)
            mip_relaxation.fix_route(fix_num)
            mip_relaxation_results = mip_relaxation.solve_relaxation(cpu=cpu_num, time_limit=time_limit_value)
            mip_root_node_obj_dic[fix_num - 1] = mip_relaxation_results['relaxation_obj']
            
            if mip_relaxation_results['relaxation_obj'] < p_median_results['obj'] + tsp_results['obj'] and fix_num < prob.t+1 :
                check_ub_lb.append(fix_num - 1)
                fix_num += 1
                iter_num += 1
            elif iter_num == 0 :
                new_list = new_list[:round(len(new_list)/3)]
                if len(new_list) > 0:
                    fix_num = new_list[round(len(new_list)/3)]
                else:
                    check_ub_lb.append(p_median_num)
                    break
            else:
                break
                
#         for fix_num in prob.A2[p_median_num:]:
#             print(f"fix_num : {fix_num}, {time.time() - start_time}")
#             fix_num = int(fix_num)
            
#             mip_relaxation = None
#             mip_relaxation = MIP(prob)
#             mip_relaxation.fix_route(fix_num)
#             mip_relaxation_results = mip_relaxation.solve_relaxation(cpu=cpu_num, time_limit=time_limit_value)
#             mip_root_node_obj_dic[fix_num - 1] = mip_relaxation_results['relaxation_obj']

#             if mip_relaxation_results['relaxation_obj'] < p_median_results['obj'] + tsp_results['obj']:
#                 print("add ub lb")
#                 check_ub_lb.append(fix_num - 1)
#             else:
#                 break

        ub_num_r_visits = max(check_ub_lb)

        w_star_dic = {}
        w_star_dic[p_median_num] = p_median_results['obj']
        
        # lb, ub 사이의 나머지 노드는 p-median을 root 노드만 풀어서 w_star를 구한다.
        for i in range(lb_num_r_visits + 1, ub_num_r_visits + 2):
            p_median_root = P_median(self.prob, i)
            p_median_root_results = p_median_root.solve_root(cpu=cpu_num, time_limit=time_limit_value)
            w_star_dic[i] = p_median_root_results['obj']

        solve_time = time.time() - start_time

        return {
            'w_star_dic': w_star_dic,
            'lb_num_r_visits': lb_num_r_visits,
            'ub_num_r_visits': ub_num_r_visits,
            'incumbent': p_median_results['obj'] + tsp_results['obj'],
            'vehicle_route': tsp_results['vehicle_path'],
            'drone_route': p_median_results['drone_route'],
            'solve_time': solve_time,
            'mip_root_node_obj_dic': mip_root_node_obj_dic,
            'cpu': cpu_num,
            'time_limit': time_limit_value
        }

if __name__ == "__main__":

    import numpy as np

    path = f'./solomon/25'  # solomon 노드개수 지정 (25,50,100)
    problem_name = 'R105'  # Augerat 문제 번호 지정
    service_time = 10 # 모든 문제 서비스 타임은 10
    scaler = 100

    cpu = 2  # cpu 개수 지정

    time_limit = 3600

    B_list = [[30,20,10]]  # 배터리 리스트 지정
    drone_multiplier_list = [[0.6,0.4,0.2]]  # multiplier 지정

    for B in B_list:
        for drone_multiplier in drone_multiplier_list:
            prob = Prob(path,
                        problem_name,
                        B,
                        drone_multiplier,
                        service_time,
                        scaler)

            results = {}

            # 휴리스틱 solve
            heuristic = Heuristic_Solver(prob)
            heuristic_results = heuristic.solve(cpu, time_limit)
            print("finish heuristic")

            # benders에 사용할 preprossing 값
            incumbent = heuristic_results['incumbent']
            lb = heuristic_results['lb_num_r_visits']
            ub = heuristic_results['ub_num_r_visits']
            w_star_dic = heuristic_results['w_star_dic']

            # # MIP solve
            # mip = MIP(prob)
            # mip_results = mip.solve(cpu, time_limit)
            # print("finish mip")
            #
            # # MIP 솔루션 figure 그리기
            # mip_vehicle_nodes = mip_results['vehicle_nodes']
            # mip_drone_nodes = mip_results['drone_nodes']
            # prob.draw_solution(mip_vehicle_nodes, mip_drone_nodes, 'mip', figure_save_path)
            #
            # # MIP Relaxation solve
            # mip_relaxation = MIP(prob)
            # mip_relaxation_results = mip.solve_relaxation(cpu, time_limit)
            # print("finish mip_r")

            # Benders solve
            benders = BendersMasterProb(prob, incumbent, lb, ub, w_star_dic)
            benders_results = benders.solve(cpu, time_limit - heuristic_results['solve_time'])
            print("finish benders")

            # Benders 솔루션 figure 그리기
            # benders_vehicle_nodes = benders_results['vehicle_nodes']
            # benders_drone_nodes = benders_results['drone_nodes']
            # prob.draw_solution(benders_vehicle_nodes, benders_drone_nodes, 'benders', figure_save_path)

            # # 결과값 저장
            # results['prob_info'] = prob.info
            # results['heuristic_results'] = heuristic_results
            # results['mip_results'] = mip_results
            # results['mip_relaxation_results'] = mip_relaxation_results
            # results['benders_results'] = benders_results
            #
            # # json 파일로 저장
            # json.dump(results,
            #           open(f'{json_save_path}/{problem_name}_{len(prob.N)}_{B}_{multiplier}.json', 'w'))


















