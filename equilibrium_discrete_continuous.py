import time
import math
import numpy as np
from scipy.optimize import minimize
import copy
from scipy.integrate import quad


# truncation to ensure boundedness of functions
def trunc(ans, lower_bound, upper_bound):
    if abs(ans) < lower_bound:
        ans = np.sign(ans)*lower_bound
    elif abs(ans) > upper_bound:
        ans = np.sign(ans)*upper_bound
    return ans

# coefficients related to SDE
def bfunc(running_time, running_state, running_control):
    # note that here, we have multiplication of control and state
    # so that when we have no short-selling constraint,
    # u is restricted to be between 0 and 1
    # ans = r*running_state + (alpha - r)*running_control*running_state
    ans = a*running_state + b*running_control
    # ans = a*running_state*running_control
    ans = trunc(ans, 0, Upper_bound_coeff)
    return ans

def sigmafunc(running_time, running_state, running_control):
    # note that here, we have multiplication of control and state
    # so that when we have no short-selling constraint,
    # u is restricted to be between 0 and 1
    # ans = sigma*running_control*running_state
    ans = sigma
    ans = trunc(ans, 0, Upper_bound_coeff)
    return ans

# functions related to objective functional
# psi corresponds to g, terminal cost
def psifunc(running_state):
    ans = trunc(running_state, 0, Upper_bound_func)
    return ans
# phi corresponds to h, running cost
def phifunc(running_state):
    ans = trunc(running_state, 0, Upper_bound_func)
    return ans
# running cost
def hfunc(initial_time, running_time, initial_state, running_state, running_control, mean_field):
    ans = running_control*running_control/2
    ans = trunc(ans, 0, Upper_bound_func)
    return ans
# terminal cost
def gfunc(initial_time, initial_state, running_state, mean_field):
    # ans = running_state - gamma*running_state*running_state/2 + gamma*mean_field*mean_field/2
    if (LQ_or_quartic == 0):
        ans = gamma*(running_state-initial_state)*(running_state-initial_state)/2
    elif (LQ_or_quartic == 1):
        ans = gamma*(running_state-initial_state)*(running_state-initial_state)* \
            (running_state-initial_state)*(running_state-initial_state)/2
    ans = trunc(ans, 0, Upper_bound_func)
    return ans

# probability laws for the markov chain
def prob_up(running_time, running_state, running_control):
    b_coeff = bfunc(running_time, running_state, running_control)
    sigma_coeff = sigmafunc(running_time, running_state, running_control)
    ans = delta_state*b_coeff/(2*Upper_bound_coeff) \
            + sigma_coeff*sigma_coeff/(2*Upper_bound_coeff*Upper_bound_coeff)
    return ans

def prob_down(running_time, running_state, running_control):
    b_coeff = bfunc(running_time, running_state, running_control)
    sigma_coeff = sigmafunc(running_time, running_state, running_control)
    ans = - delta_state*b_coeff/(2*Upper_bound_coeff) \
            + sigma_coeff*sigma_coeff/(2*Upper_bound_coeff*Upper_bound_coeff)
    return ans

def prob_stay(running_time, running_state, running_control):
    sigma_coeff = sigmafunc(running_time, running_state, running_control)
    ans = 1 -sigma_coeff*sigma_coeff/(Upper_bound_coeff*Upper_bound_coeff)
    return ans

# recursion for psi
def psirecursion(to_optim_U, from_time_index, start_state_index):
    temp = np.zeros(2*N*M+1)
    temp_inter = np.zeros(2*M+1)
    from_time_index_N = int(from_time_index / M)
    for tt_N in range(N-1, from_time_index_N-1, -1):
        # when tt_N matches the time we start
        # we should use to_optim_U
        # else, we should use optimal_U
        if tt_N == from_time_index_N:
            to_optim_U_flag = True
        else:
            to_optim_U_flag = False
        temptemp = copy.copy(temp)
        del_N = tt_N - from_time_index_N
        for statediff_index in range(-del_N*M, del_N*M+1):
            # this block is especially for the point at M-1
            # hence, when we are at N-1
            # we initiate temp to be its weighted terminal cost
            if tt_N == N-1:
                for statediff_index_inter in range(-(M-1), (M-1)+1):
                    temptemp_inter_index = M + statediff_index_inter
                    temp_index = start_state_index + statediff_index + statediff_index_inter
                    # temptemp_inter_up = temptemp[temp_index+1]
                    # temptemp_inter_down = temptemp[temp_index-1]
                    # temptemp_inter_stay = temptemp[temp_index]
                    state_up = state_space[temp_index+1]
                    state_down = state_space[temp_index-1]
                    state_stay = state_space[temp_index]
                    tt_now = time_space[tt_N*M+(M-1)]
                    if to_optim_U_flag:
                        uu_now = to_optim_U[(M-1)*(M-1) + M-1 + statediff_index_inter]
                    else:
                        uu_now = optimal_U[start_state_index + statediff_index, \
                                tt_N*M*M+((M-1)*(M-1) + M-1 + statediff_index_inter)]
                    temp_inter[temptemp_inter_index] \
                        = prob_up(tt_now, state_stay, uu_now)*psifunc(state_up) \
                            + prob_down(tt_now, state_stay, uu_now)*psifunc(state_down) \
                            + prob_stay(tt_now, state_stay, uu_now)*psifunc(state_stay)
            # otherwise, we should fill temp_inter using weighted temptemp
            else:
                for statediff_index_inter in range(-(M-1), (M-1)+1):
                    temptemp_inter_index = M + statediff_index_inter
                    temp_index = start_state_index + statediff_index + statediff_index_inter
                    temptemp_inter_up = temptemp[temp_index+1]
                    temptemp_inter_down = temptemp[temp_index-1]
                    temptemp_inter_stay = temptemp[temp_index]
                    state_stay = state_space[temp_index]
                    tt_now = time_space[tt_N*M+(M-1)]
                    if to_optim_U_flag:
                        uu_now = to_optim_U[(M-1)*(M-1) + M-1 + statediff_index_inter]
                    else:
                        uu_now = optimal_U[start_state_index + statediff_index, \
                                tt_N*M*M+((M-1)*(M-1) + M-1 + statediff_index_inter)]
                    temp_inter[temptemp_inter_index] \
                        = prob_up(tt_now, state_stay, uu_now)*temptemp_inter_up \
                            + prob_down(tt_now, state_stay, uu_now)*temptemp_inter_down \
                            + prob_stay(tt_now, state_stay, uu_now)*temptemp_inter_stay
            # else, we update from M-2 till 0
            for tt_N_inter in range(M-2, -1, -1):
                temptemp_inter = copy.copy(temp_inter)
                for statediff_index_inter in range(-tt_N_inter, tt_N_inter+1):
                    # note that here, 0 is left end, M is middle, 2M is right end
                    temptemp_inter_index = M + statediff_index_inter
                    temp_index = start_state_index + statediff_index + statediff_index_inter
                    temptemp_inter_up = temptemp_inter[temptemp_inter_index+1]
                    temptemp_inter_down = temptemp_inter[temptemp_inter_index-1]
                    temptemp_inter_stay = temptemp_inter[temptemp_inter_index]
                    state_stay = state_space[temp_index]
                    tt_now = time_space[tt_N*M+tt_N_inter]
                    if to_optim_U_flag:
                        uu_now = to_optim_U[tt_N_inter*tt_N_inter + tt_N_inter + statediff_index_inter]
                    else:
                        uu_now = optimal_U[start_state_index + statediff_index, \
                                tt_N*M*M+(tt_N_inter*tt_N_inter + tt_N_inter + statediff_index_inter)]
                    temp_inter[temptemp_inter_index] \
                        = prob_up(tt_now, state_stay, uu_now)*temptemp_inter_up \
                            + prob_down(tt_now, state_stay, uu_now)*temptemp_inter_down \
                            + prob_stay(tt_now, state_stay, uu_now)*temptemp_inter_stay
            # after one N-person game is done, we should update the value of temp
            # note that for temp_inter, 0 is left end, M is middle, 2M is right end
            temp[start_state_index + statediff_index] = temp_inter[M]
    return temp[start_state_index]

# recursion for phi
def phirecursion(to_optim_U, from_time_index, to_time_index, start_state_index):
    # when from and to is the same, we do not need to do calculation at all
    if to_time_index == from_time_index:
        ans = phifunc(state_space[start_state_index])
        return(ans)
    temp = np.zeros(2*N*M+1)
    temp_inter = np.zeros(2*M+1)
    from_time_index_N = math.floor(from_time_index / M)
    to_time_index_N = math.floor(to_time_index / M)
    to_time_index_M = to_time_index % M
    tt_N = to_time_index_N
    if tt_N == from_time_index_N:
        to_optim_U_flag = True
    else:
        to_optim_U_flag = False
    del_N = to_time_index_N - from_time_index_N
    # starts from M
    if to_time_index_M == 0:
        for statediff_index in range(-del_N*M, del_N*M+1):
            # at M-1
            for statediff_index_inter in range(-(M-1), (M-1)+1):
                temptemp_inter_index = M + statediff_index_inter
                temp_index = start_state_index + statediff_index + statediff_index_inter
                # temptemp_inter_up = temptemp[temp_index+1]
                # temptemp_inter_down = temptemp[temp_index-1]
                # temptemp_inter_stay = temptemp[temp_index]
                state_up = state_space[temp_index+1]
                state_down = state_space[temp_index-1]
                state_stay = state_space[temp_index]
                tt_now = time_space[tt_N*M+(M-1)]
                if to_optim_U_flag:
                    uu_now = to_optim_U[(M-1)*(M-1) + M-1 + statediff_index_inter]
                else:
                    uu_now = optimal_U[start_state_index + statediff_index, \
                            tt_N*M*M+((M-1)*(M-1) + M-1 + statediff_index_inter)]
                temp_inter[temptemp_inter_index] \
                    = prob_up(tt_now, state_stay, uu_now)\
                            *phifunc(state_up) \
                        + prob_down(tt_now, state_stay, uu_now)\
                            *phifunc(state_down) \
                        + prob_stay(tt_now, state_stay, uu_now)\
                            *phifunc(state_stay)
            # loop from M-2 to the end
            for tt_N_inter in range(M-2, -1, -1):
                temptemp_inter = copy.copy(temp_inter)
                for statediff_index_inter in range(-tt_N_inter, tt_N_inter+1):
                    # note that here, 0 is left end, M is middle, 2M is right end
                    temptemp_inter_index = M + statediff_index_inter
                    temp_index = start_state_index + statediff_index + statediff_index_inter
                    temptemp_inter_up = temptemp_inter[temptemp_inter_index+1]
                    temptemp_inter_down = temptemp_inter[temptemp_inter_index-1]
                    temptemp_inter_stay = temptemp_inter[temptemp_inter_index]
                    state_stay = state_space[temp_index]
                    tt_now = time_space[tt_N*M+tt_N_inter]
                    if to_optim_U_flag:
                        uu_now = to_optim_U[tt_N_inter*tt_N_inter + tt_N_inter + statediff_index_inter]
                    else:
                        uu_now = optimal_U[start_state_index + statediff_index, \
                                tt_N*M*M+(tt_N_inter*tt_N_inter + tt_N_inter + statediff_index_inter)]
                    temp_inter[temptemp_inter_index] \
                        = prob_up(tt_now, state_stay, uu_now)*temptemp_inter_up \
                            + prob_down(tt_now, state_stay, uu_now)*temptemp_inter_down \
                            + prob_stay(tt_now, state_stay, uu_now)*temptemp_inter_stay
            temp[start_state_index + statediff_index] = temp_inter[M]
    else:
        for statediff_index in range(-del_N*to_time_index_M, del_N*to_time_index_M+1):
            # at to_time_index_M-1, where to_time_index_M is not M
            for statediff_index_inter in range(-(to_time_index_M-1), (to_time_index_M-1)+1):
                temptemp_inter_index = M + statediff_index_inter
                temp_index = start_state_index + statediff_index + statediff_index_inter
                # temptemp_inter_up = temptemp[temp_index+1]
                # temptemp_inter_down = temptemp[temp_index-1]
                # temptemp_inter_stay = temptemp[temp_index]
                state_up = state_space[temp_index+1]
                state_down = state_space[temp_index-1]
                state_stay = state_space[temp_index]
                tt_now = time_space[tt_N*M+(to_time_index_M-1)]
                if to_optim_U_flag:
                    uu_now = to_optim_U[(to_time_index_M-1)*(to_time_index_M-1) \
                                + to_time_index_M-1 + statediff_index_inter]
                else:
                    uu_now = optimal_U[start_state_index + statediff_index, \
                            tt_N*M*M+((to_time_index_M-1)*(to_time_index_M-1) \
                                + to_time_index_M-1 + statediff_index_inter)]
                temp_inter[temptemp_inter_index] \
                    = prob_up(tt_now, state_stay, uu_now)\
                            *psifunc(state_up) \
                        + prob_down(tt_now, state_stay, uu_now)\
                            *psifunc(state_down) \
                        + prob_stay(tt_now, state_stay, uu_now)\
                            *psifunc(state_stay)
            # loop from M-2 to the end
            for tt_N_inter in range(to_time_index_M-2, -1, -1):
                temptemp_inter = copy.copy(temp_inter)
                for statediff_index_inter in range(-tt_N_inter, tt_N_inter+1):
                    # note that here, 0 is left end, M is middle, 2M is right end
                    temptemp_inter_index = M + statediff_index_inter
                    temp_index = start_state_index + statediff_index + statediff_index_inter
                    temptemp_inter_up = temptemp_inter[temptemp_inter_index+1]
                    temptemp_inter_down = temptemp_inter[temptemp_inter_index-1]
                    temptemp_inter_stay = temptemp_inter[temptemp_inter_index]
                    state_stay = state_space[temp_index]
                    tt_now = time_space[tt_N*M+tt_N_inter]
                    if to_optim_U_flag:
                        uu_now = to_optim_U[tt_N_inter*tt_N_inter + tt_N_inter + statediff_index_inter]
                    else:
                        uu_now = optimal_U[start_state_index + statediff_index, \
                                tt_N*M*M+(tt_N_inter*tt_N_inter + tt_N_inter + statediff_index_inter)]
                    temp_inter[temptemp_inter_index] \
                        = prob_up(tt_now, state_stay, uu_now)*temptemp_inter_up \
                            + prob_down(tt_now, state_stay, uu_now)*temptemp_inter_down \
                            + prob_stay(tt_now, state_stay, uu_now)*temptemp_inter_stay
            temp[start_state_index + statediff_index] = temp_inter[M]
    for tt_N in range(to_time_index_N-1, from_time_index_N-1, -1):
        # when tt_N matches the time we start
        # we should use to_optim_U
        # else, we should use optimal_U
        if tt_N == from_time_index_N:
            to_optim_U_flag = True
        else:
            to_optim_U_flag = False
        temptemp = copy.copy(temp)
        del_N = tt_N - from_time_index_N
        for statediff_index in range(-del_N*M, del_N*M+1):
            # this block is especially for the point at M-1
            for statediff_index_inter in range(-(M-1), (M-1)+1):
                temptemp_inter_index = M + statediff_index_inter
                temp_index = start_state_index + statediff_index + statediff_index_inter
                temptemp_inter_up = temptemp[temp_index+1]
                temptemp_inter_down = temptemp[temp_index-1]
                temptemp_inter_stay = temptemp[temp_index]
                state_stay = state_space[temp_index]
                tt_now = time_space[tt_N*M+(M-1)]
                if to_optim_U_flag:
                    uu_now = to_optim_U[(M-1)*(M-1) + M-1 + statediff_index_inter]
                else:
                    uu_now = optimal_U[start_state_index + statediff_index, \
                            tt_N*M*M+((M-1)*(M-1) + M-1 + statediff_index_inter)]
                temp_inter[temptemp_inter_index] \
                    = prob_up(tt_now, state_stay, uu_now)*temptemp_inter_up \
                        + prob_down(tt_now, state_stay, uu_now)*temptemp_inter_down \
                        + prob_stay(tt_now, state_stay, uu_now)*temptemp_inter_stay
            # else, we update from M-2 till 0
            for tt_N_inter in range(M-2, -1, -1):
                temptemp_inter = copy.copy(temp_inter)
                for statediff_index_inter in range(-tt_N_inter, tt_N_inter+1):
                    # note that here, 0 is left end, M is middle, 2M is right end
                    temptemp_inter_index = M + statediff_index_inter
                    temp_index = start_state_index + statediff_index + statediff_index_inter
                    temptemp_inter_up = temptemp_inter[temptemp_inter_index+1]
                    temptemp_inter_down = temptemp_inter[temptemp_inter_index-1]
                    temptemp_inter_stay = temptemp_inter[temptemp_inter_index]
                    state_stay = state_space[temp_index]
                    tt_now = time_space[tt_N*M+tt_N_inter]
                    if to_optim_U_flag:
                        uu_now = to_optim_U[tt_N_inter*tt_N_inter + tt_N_inter + statediff_index_inter]
                    else:
                        uu_now = optimal_U[start_state_index + statediff_index, \
                                tt_N*M*M+(tt_N_inter*tt_N_inter + tt_N_inter + statediff_index_inter)]
                    temp_inter[temptemp_inter_index] \
                        = prob_up(tt_now, state_stay, uu_now)*temptemp_inter_up \
                            + prob_down(tt_now, state_stay, uu_now)*temptemp_inter_down \
                            + prob_stay(tt_now, state_stay, uu_now)*temptemp_inter_stay
                # print(temptemp_inter==temp_inter)
            # after one N-person game is done, we should update the value of temp
            # note that for temp_inter, 0 is left end, M is middle, 2M is right end
            temp[start_state_index + statediff_index] = temp_inter[M]
    return temp[start_state_index]

# recursion for g
def grecursion(to_optim_U, from_time_index, start_state_index, initial_time, initial_state):
    # # since no mean field is involved in this problem, we can safely make it 0
    # # to improve the performance
    # mean_field = psirecursion(to_optim_U, from_time_index, start_state_index)
    mean_field = 0
    temp = np.zeros(2*N*M+1)
    temp_inter = np.zeros(2*M+1)
    from_time_index_N = int(from_time_index / M)
    for tt_N in range(N-1, from_time_index_N-1, -1):
        # when tt_N matches the time we start
        # we should use to_optim_U
        # else, we should use optimal_U
        if tt_N == from_time_index_N:
            to_optim_U_flag = True
        else:
            to_optim_U_flag = False
        temptemp = copy.copy(temp)
        del_N = tt_N - from_time_index_N
        for statediff_index in range(-del_N*M, del_N*M+1):
            # this block is especially for the point at M-1
            # hence, when we are at N-1
            # we initiate temp to be its weighted terminal cost
            if tt_N == N-1:
                for statediff_index_inter in range(-(M-1), (M-1)+1):
                    temptemp_inter_index = M + statediff_index_inter
                    temp_index = start_state_index + statediff_index + statediff_index_inter
                    # temptemp_inter_up = temptemp[temp_index+1]
                    # temptemp_inter_down = temptemp[temp_index-1]
                    # temptemp_inter_stay = temptemp[temp_index]
                    state_up = state_space[temp_index+1]
                    state_down = state_space[temp_index-1]
                    state_stay = state_space[temp_index]
                    tt_now = time_space[tt_N*M+(M-1)]
                    if to_optim_U_flag:
                        uu_now = to_optim_U[(M-1)*(M-1) + M-1 + statediff_index_inter]
                    else:
                        uu_now = optimal_U[start_state_index + statediff_index, \
                                tt_N*M*M+((M-1)*(M-1) + M-1 + statediff_index_inter)]
                    temp_inter[temptemp_inter_index] \
                        = prob_up(tt_now, state_stay, uu_now)\
                                *gfunc(initial_time, initial_state, state_up, mean_field) \
                            + prob_down(tt_now, state_stay, uu_now)\
                                *gfunc(initial_time, initial_state, state_down, mean_field) \
                            + prob_stay(tt_now, state_stay, uu_now)\
                                *gfunc(initial_time, initial_state, state_stay, mean_field)
            # otherwise, we should fill temp_inter using weighted temptemp
            else:
                for statediff_index_inter in range(-(M-1), (M-1)+1):
                    temptemp_inter_index = M + statediff_index_inter
                    temp_index = start_state_index + statediff_index + statediff_index_inter
                    temptemp_inter_up = temptemp[temp_index+1]
                    temptemp_inter_down = temptemp[temp_index-1]
                    temptemp_inter_stay = temptemp[temp_index]
                    state_stay = state_space[temp_index]
                    tt_now = time_space[tt_N*M+(M-1)]
                    if to_optim_U_flag:
                        uu_now = to_optim_U[(M-1)*(M-1) + M-1 + statediff_index_inter]
                    else:
                        uu_now = optimal_U[start_state_index + statediff_index, \
                                tt_N*M*M+((M-1)*(M-1) + M-1 + statediff_index_inter)]
                    temp_inter[temptemp_inter_index] \
                        = prob_up(tt_now, state_stay, uu_now)*temptemp_inter_up \
                            + prob_down(tt_now, state_stay, uu_now)*temptemp_inter_down \
                            + prob_stay(tt_now, state_stay, uu_now)*temptemp_inter_stay
            # else, we update from M-2 till 0
            for tt_N_inter in range(M-2, -1, -1):
                temptemp_inter = copy.copy(temp_inter)
                for statediff_index_inter in range(-tt_N_inter, tt_N_inter+1):
                    # note that here, 0 is left end, M is middle, 2M is right end
                    temptemp_inter_index = M + statediff_index_inter
                    temp_index = start_state_index + statediff_index + statediff_index_inter
                    temptemp_inter_up = temptemp_inter[temptemp_inter_index+1]
                    temptemp_inter_down = temptemp_inter[temptemp_inter_index-1]
                    temptemp_inter_stay = temptemp_inter[temptemp_inter_index]
                    state_stay = state_space[temp_index]
                    tt_now = time_space[tt_N*M+tt_N_inter]
                    if to_optim_U_flag:
                        uu_now = to_optim_U[tt_N_inter*tt_N_inter + tt_N_inter + statediff_index_inter]
                    else:
                        uu_now = optimal_U[start_state_index + statediff_index, \
                                tt_N*M*M+(tt_N_inter*tt_N_inter + tt_N_inter + statediff_index_inter)]
                    temp_inter[temptemp_inter_index] \
                        = prob_up(tt_now, state_stay, uu_now)*temptemp_inter_up \
                            + prob_down(tt_now, state_stay, uu_now)*temptemp_inter_down \
                            + prob_stay(tt_now, state_stay, uu_now)*temptemp_inter_stay
                # print(temptemp_inter==temp_inter)
            # after one N-person game is done, we should update the value of temp
            # note that for temp_inter, 0 is left end, M is middle, 2M is right end
            temp[start_state_index + statediff_index] = temp_inter[M]
    return temp[start_state_index]

# recursion for h
def hrecursion(to_optim_U, from_time_index, start_state_index, initial_time, initial_state):
    temp = np.zeros(2*N*M+1)
    temp_inter = np.zeros(2*M+1)
    from_time_index_N = int(from_time_index / M)
    for tt_N in range(N-1, from_time_index_N-1, -1):
        # when tt_N matches the time we start
        # we should use to_optim_U
        # else, we should use optimal_U
        if tt_N == from_time_index_N:
            to_optim_U_flag = True
        else:
            to_optim_U_flag = False
        temptemp = copy.copy(temp)
        del_N = tt_N - from_time_index_N
        for statediff_index in range(-del_N*M, del_N*M+1):
            # this block is especially for the point at M-1
            # hence, when we are at N-1
            # we initiate temp to be its weighted terminal cost
            if tt_N == N-1:
                for statediff_index_inter in range(-(M-1), (M-1)+1):
                    temptemp_inter_index = M + statediff_index_inter
                    temp_index = start_state_index + statediff_index + statediff_index_inter
                    state_stay = state_space[temp_index]
                    tt_now = time_space[tt_N*M+(M-1)]
                    if to_optim_U_flag:
                        uu_now = to_optim_U[(M-1)*(M-1) + M-1 + statediff_index_inter]
                    else:
                        uu_now = optimal_U[start_state_index + statediff_index, \
                                tt_N*M*M+((M-1)*(M-1) + M-1 + statediff_index_inter)]
                    # # since no mean field is involved in this problem, we can safely make it 0
                    # # to improve the performance
                    # mean_field = phirecursion(to_optim_U, from_time_index, tt_N*M+(M-1), start_state_index)
                    mean_field = 0
                    temp_inter[temptemp_inter_index] \
                        = delta_time*hfunc(initial_time, tt_now, initial_state, \
                                   state_stay, uu_now, mean_field)
            # otherwise, we should fill temp_inter using weighted temptemp
            else:
                for statediff_index_inter in range(-(M-1), (M-1)+1):
                    temptemp_inter_index = M + statediff_index_inter
                    temp_index = start_state_index + statediff_index + statediff_index_inter
                    temptemp_inter_up = temptemp[temp_index+1]
                    temptemp_inter_down = temptemp[temp_index-1]
                    temptemp_inter_stay = temptemp[temp_index]
                    state_stay = state_space[temp_index]
                    tt_now = time_space[tt_N*M+(M-1)]
                    if to_optim_U_flag:
                        uu_now = to_optim_U[(M-1)*(M-1) + M-1 + statediff_index_inter]
                    else:
                        uu_now = optimal_U[start_state_index + statediff_index, \
                                tt_N*M*M+((M-1)*(M-1) + M-1 + statediff_index_inter)]
                    # # since no mean field is involved in this problem, we can safely make it 0
                    # # to improve the performance
                    # mean_field = phirecursion(to_optim_U, from_time_index, tt_N*M+(M-1), start_state_index)
                    mean_field = 0
                    temp_inter[temptemp_inter_index] \
                        = prob_up(tt_now, state_stay, uu_now)*temptemp_inter_up \
                            + prob_down(tt_now, state_stay, uu_now)*temptemp_inter_down \
                            + prob_stay(tt_now, state_stay, uu_now)*temptemp_inter_stay
                    temp_inter[temptemp_inter_index] += (delta_time* \
                                hfunc(initial_time, tt_now, initial_state, state_stay, uu_now, mean_field))
            # else, we update from M-2 till 0
            for tt_N_inter in range(M-2, -1, -1):
                temptemp_inter = copy.copy(temp_inter)
                for statediff_index_inter in range(-tt_N_inter, tt_N_inter+1):
                    # note that here, 0 is left end, M is middle, 2M is right end
                    temptemp_inter_index = M + statediff_index_inter
                    temp_index = start_state_index + statediff_index + statediff_index_inter
                    temptemp_inter_up = temptemp_inter[temptemp_inter_index+1]
                    temptemp_inter_down = temptemp_inter[temptemp_inter_index-1]
                    temptemp_inter_stay = temptemp_inter[temptemp_inter_index]
                    state_stay = state_space[temp_index]
                    tt_now = time_space[tt_N*M+tt_N_inter]
                    if to_optim_U_flag:
                        uu_now = to_optim_U[tt_N_inter*tt_N_inter + tt_N_inter + statediff_index_inter]
                    else:
                        uu_now = optimal_U[start_state_index + statediff_index, \
                                tt_N*M*M+(tt_N_inter*tt_N_inter + tt_N_inter + statediff_index_inter)]
                    # # since no mean field is involved in this problem, we can safely make it 0
                    # # to improve the performance
                    # mean_field = phirecursion(to_optim_U, from_time_index, tt_N*M+tt_N_inter, start_state_index)
                    mean_field = 0
                    temp_inter[temptemp_inter_index] \
                        = prob_up(tt_now, state_stay, uu_now)*temptemp_inter_up \
                            + prob_down(tt_now, state_stay, uu_now)*temptemp_inter_down \
                            + prob_stay(tt_now, state_stay, uu_now)*temptemp_inter_stay
                    temp_inter[temptemp_inter_index] += (delta_time* \
                                hfunc(initial_time, tt_now, initial_state, state_stay, uu_now, mean_field))
                # print(temptemp_inter==temp_inter)
            # after one N-person game is done, we should update the value of temp
            # note that for temp_inter, 0 is left end, M is middle, 2M is right end
            temp[start_state_index + statediff_index] = temp_inter[M]
    return temp[start_state_index]

# recursion for objective function
def jrecursion(to_optim_U, from_time_index, start_state_index, initial_time, initial_state):
    jfunc = grecursion(to_optim_U, from_time_index, start_state_index, initial_time, \
        initial_state) + hrecursion(to_optim_U, from_time_index, start_state_index, \
        initial_time, initial_state)
    return jfunc

# recursion for objective function
def jrecursion_quicker(to_optim_U, time_index_N, start_state_index, time_index_M, state_diff):
    tt_now = time_space[time_index_N*M+time_index_M]
    initial_state = state_space[start_state_index]
    state_now = state_space[start_state_index + state_diff]
    state_up = state_space[start_state_index + state_diff + 1]
    state_down = state_space[start_state_index + state_diff - 1]
    uu_now = to_optim_U[0]
    if time_index_M == (M-1):
        if time_index_N == (N-1):
            if (LQ_or_quartic == 0):
                # for LQ
                xsquared_usquared_up = gamma*state_up*state_up/2
                xsquared_usquared_down = gamma*state_down*state_down/2
                xsquared_usquared_stay = gamma*state_now*state_now/2
                x_up = state_up
                x_down = state_down
                x_stay = state_now
            elif (LQ_or_quartic == 1):
                # for quartic
                xfour_usquared_up = gamma*state_up*state_up*state_up*state_up/2
                xfour_usquared_down = gamma*state_down*state_down*state_down*state_down/2
                xfour_usquared_stay = gamma*state_now*state_now*state_now*state_now/2
                xcube_up = state_up*state_up*state_up
                xcube_down = state_down*state_down*state_down
                xcube_stay = state_now*state_now*state_now
                xsquare_up = state_up*state_up
                xsquare_down = state_down*state_down
                xsquare_stay = state_now*state_now
                x_up = state_up
                x_down = state_down
                x_stay = state_now
        else:
            if (LQ_or_quartic == 0):
                # for LQ
                xsquared_usquared_up = xsquared_usquared_iter[start_state_index+state_diff+1, (time_index_N+1)*M*M]
                xsquared_usquared_down = xsquared_usquared_iter[start_state_index+state_diff-1, (time_index_N+1)*M*M]
                xsquared_usquared_stay = xsquared_usquared_iter[start_state_index+state_diff, (time_index_N+1)*M*M]
                x_up = x_iter[start_state_index+state_diff+1, (time_index_N+1)*M*M]
                x_down = x_iter[start_state_index+state_diff-1, (time_index_N+1)*M*M]
                x_stay = x_iter[start_state_index+state_diff, (time_index_N+1)*M*M]
            elif (LQ_or_quartic == 1):
                # for quartic
                xfour_usquared_up = xfour_usquared_iter[start_state_index+state_diff+1, (time_index_N+1)*M*M]
                xfour_usquared_down = xfour_usquared_iter[start_state_index+state_diff-1, (time_index_N+1)*M*M]
                xfour_usquared_stay = xfour_usquared_iter[start_state_index+state_diff, (time_index_N+1)*M*M]
                xcube_up = xcube_iter[start_state_index+state_diff+1, (time_index_N+1)*M*M]
                xcube_down = xcube_iter[start_state_index+state_diff-1, (time_index_N+1)*M*M]
                xcube_stay = xcube_iter[start_state_index+state_diff, (time_index_N+1)*M*M]
                xsquare_up = xsquare_iter[start_state_index+state_diff+1, (time_index_N+1)*M*M]
                xsquare_down = xsquare_iter[start_state_index+state_diff-1, (time_index_N+1)*M*M]
                xsquare_stay = xsquare_iter[start_state_index+state_diff, (time_index_N+1)*M*M]
                x_up = x_iter[start_state_index+state_diff+1, (time_index_N+1)*M*M]
                x_down = x_iter[start_state_index+state_diff-1, (time_index_N+1)*M*M]
                x_stay = x_iter[start_state_index+state_diff, (time_index_N+1)*M*M]
    else:
        if (LQ_or_quartic == 0):
            # for LQ
            xsquared_usquared_up = xsquared_usquared_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter + 1)]
            xsquared_usquared_down = xsquared_usquared_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter - 1)]
            xsquared_usquared_stay = xsquared_usquared_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter)]
            x_up = x_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter + 1)]
            x_down = x_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter - 1)]
            x_stay = x_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter)]
        elif (LQ_or_quartic == 1):
            # for quartic
            xfour_usquared_up = xfour_usquared_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter + 1)]
            xfour_usquared_down = xfour_usquared_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter - 1)]
            xfour_usquared_stay = xfour_usquared_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter)]
            xcube_up = xcube_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter + 1)]
            xcube_down = xcube_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter - 1)]
            xcube_stay = xcube_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter)]
            xsquare_up = xsquare_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter + 1)]
            xsquare_down = xsquare_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter - 1)]
            xsquare_stay = xsquare_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter)]
            x_up = x_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter + 1)]
            x_down = x_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter - 1)]
            x_stay = x_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter)]

    if (LQ_or_quartic == 0):
        # for LQ
        ans = prob_up(tt_now, state_now, uu_now) * (xsquared_usquared_up - gamma*x_up*initial_state) \
            + prob_down(tt_now, state_now, uu_now) * (xsquared_usquared_down - gamma*x_down*initial_state) \
            + prob_stay(tt_now, state_now, uu_now) * (xsquared_usquared_stay - gamma*x_stay*initial_state) \
            + delta_time * uu_now * uu_now / 2 + gamma * initial_state * initial_state / 2
    elif (LQ_or_quartic == 1):
        # for quartic
        ans = prob_up(tt_now, state_now, uu_now) * (xfour_usquared_up \
            - 2*gamma*xcube_up*initial_state + 3*gamma*xsquare_up*initial_state*initial_state \
             - 2*gamma*x_up*initial_state*initial_state*initial_state) \
            + prob_down(tt_now, state_now, uu_now) * (xfour_usquared_down \
            - 2*gamma*xcube_down*initial_state + 3*gamma*xsquare_down*initial_state*initial_state \
             - 2*gamma*x_down*initial_state*initial_state*initial_state) \
            + prob_stay(tt_now, state_now, uu_now) * (xfour_usquared_stay \
            - 2*gamma*xcube_stay*initial_state + 3*gamma*xsquare_stay*initial_state*initial_state \
             - 2*gamma*x_stay*initial_state*initial_state*initial_state) \
            + delta_time * uu_now * uu_now / 2 \
            + gamma * initial_state * initial_state * initial_state * initial_state / 2
    return ans

def jrecursion_fill(time_index_N, start_state_index, time_index_M, state_diff):
    tt_now = time_space[time_index_N*M+time_index_M]
    initial_state = state_space[start_state_index]
    state_now = state_space[start_state_index + state_diff]
    state_up = state_space[start_state_index + state_diff + 1]
    state_down = state_space[start_state_index + state_diff - 1]
    uu_now = optimal_U[start_state_index, time_index_N*M*M+(time_index_M*time_index_M + time_index_M + state_diff)]
    if time_index_M == (M-1):
        if time_index_N == (N-1):
            if (LQ_or_quartic == 0):
                # for LQ
                xsquared_usquared_up = gamma*state_up*state_up/2
                xsquared_usquared_down = gamma*state_down*state_down/2
                xsquared_usquared_stay = gamma*state_now*state_now/2
                x_up = state_up
                x_down = state_down
                x_stay = state_now
            elif (LQ_or_quartic == 1):
                # for quartic
                xfour_usquared_up = gamma*state_up*state_up*state_up*state_up/2
                xfour_usquared_down = gamma*state_down*state_down*state_down*state_down/2
                xfour_usquared_stay = gamma*state_now*state_now*state_now*state_now/2
                xcube_up = state_up*state_up*state_up
                xcube_down = state_down*state_down*state_down
                xcube_stay = state_now*state_now*state_now
                xsquare_up = state_up*state_up
                xsquare_down = state_down*state_down
                xsquare_stay = state_now*state_now
                x_up = state_up
                x_down = state_down
                x_stay = state_now
        else:
            if (LQ_or_quartic == 0):
                # for LQ
                xsquared_usquared_up = xsquared_usquared_iter[start_state_index+state_diff+1, (time_index_N+1)*M*M]
                xsquared_usquared_down = xsquared_usquared_iter[start_state_index+state_diff-1, (time_index_N+1)*M*M]
                xsquared_usquared_stay = xsquared_usquared_iter[start_state_index+state_diff, (time_index_N+1)*M*M]
                x_up = x_iter[start_state_index+state_diff+1, (time_index_N+1)*M*M]
                x_down = x_iter[start_state_index+state_diff-1, (time_index_N+1)*M*M]
                x_stay = x_iter[start_state_index+state_diff, (time_index_N+1)*M*M]
            elif (LQ_or_quartic == 1):
                # for quartic
                xfour_usquared_up = xfour_usquared_iter[start_state_index+state_diff+1, (time_index_N+1)*M*M]
                xfour_usquared_down = xfour_usquared_iter[start_state_index+state_diff-1, (time_index_N+1)*M*M]
                xfour_usquared_stay = xfour_usquared_iter[start_state_index+state_diff, (time_index_N+1)*M*M]
                xcube_up = xcube_iter[start_state_index+state_diff+1, (time_index_N+1)*M*M]
                xcube_down = xcube_iter[start_state_index+state_diff-1, (time_index_N+1)*M*M]
                xcube_stay = xcube_iter[start_state_index+state_diff, (time_index_N+1)*M*M]
                xsquare_up = xsquare_iter[start_state_index+state_diff+1, (time_index_N+1)*M*M]
                xsquare_down = xsquare_iter[start_state_index+state_diff-1, (time_index_N+1)*M*M]
                xsquare_stay = xsquare_iter[start_state_index+state_diff, (time_index_N+1)*M*M]
                x_up = x_iter[start_state_index+state_diff+1, (time_index_N+1)*M*M]
                x_down = x_iter[start_state_index+state_diff-1, (time_index_N+1)*M*M]
                x_stay = x_iter[start_state_index+state_diff, (time_index_N+1)*M*M]
    else:
        if (LQ_or_quartic == 0):
            # for LQ
            xsquared_usquared_up = xsquared_usquared_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter + 1)]
            xsquared_usquared_down = xsquared_usquared_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter - 1)]
            xsquared_usquared_stay = xsquared_usquared_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter)]
            x_up = x_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter + 1)]
            x_down = x_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter - 1)]
            x_stay = x_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter)]
        elif (LQ_or_quartic == 1):
            # for quartic
            xfour_usquared_up = xfour_usquared_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter + 1)]
            xfour_usquared_down = xfour_usquared_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter - 1)]
            xfour_usquared_stay = xfour_usquared_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter)]
            xcube_up = xcube_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter + 1)]
            xcube_down = xcube_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter - 1)]
            xcube_stay = xcube_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter)]
            xsquare_up = xsquare_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter + 1)]
            xsquare_down = xsquare_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter - 1)]
            xsquare_stay = xsquare_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter)]
            x_up = x_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter + 1)]
            x_down = x_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter - 1)]
            x_stay = x_iter[start_state_index, time_index_N*M*M+((time_index_M+1)*(time_index_M+1) + (time_index_M+1) + statediff_index_inter)]

    if (LQ_or_quartic == 0):
        # for LQ
        xsquared_usquared_iter[start_state_index, time_index_N*M*M+(time_index_M*time_index_M + time_index_M + state_diff)] = prob_up(tt_now, state_now, uu_now) * xsquared_usquared_up \
                    + prob_down(tt_now, state_now, uu_now) * xsquared_usquared_down \
                    + prob_stay(tt_now, state_now, uu_now) * xsquared_usquared_stay \
                    + delta_time * uu_now * uu_now / 2

        x_iter[start_state_index, time_index_N*M*M+(time_index_M*time_index_M + time_index_M + state_diff)] = prob_up(tt_now, state_now, uu_now) * x_up \
                    + prob_down(tt_now, state_now, uu_now) * x_down \
                    + prob_stay(tt_now, state_now, uu_now) * x_stay
    elif (LQ_or_quartic == 1):
        # for quartic
        xfour_usquared_iter[start_state_index, time_index_N*M*M+(time_index_M*time_index_M + time_index_M + state_diff)] = prob_up(tt_now, state_now, uu_now) * xfour_usquared_up \
                    + prob_down(tt_now, state_now, uu_now) * xfour_usquared_down \
                    + prob_stay(tt_now, state_now, uu_now) * xfour_usquared_stay \
                    + delta_time * uu_now * uu_now / 2

        xcube_iter[start_state_index, time_index_N*M*M+(time_index_M*time_index_M + time_index_M + state_diff)] = prob_up(tt_now, state_now, uu_now) * xcube_up \
                    + prob_down(tt_now, state_now, uu_now) * xcube_down \
                    + prob_stay(tt_now, state_now, uu_now) * xcube_stay

        xsquare_iter[start_state_index, time_index_N*M*M+(time_index_M*time_index_M + time_index_M + state_diff)] = prob_up(tt_now, state_now, uu_now) * xsquare_up \
                    + prob_down(tt_now, state_now, uu_now) * xsquare_down \
                    + prob_stay(tt_now, state_now, uu_now) * xsquare_stay

        x_iter[start_state_index, time_index_N*M*M+(time_index_M*time_index_M + time_index_M + state_diff)] = prob_up(tt_now, state_now, uu_now) * x_up \
                    + prob_down(tt_now, state_now, uu_now) * x_down \
                    + prob_stay(tt_now, state_now, uu_now) * x_stay


# exact solution...
def beta(s):
    return(gamma*math.exp(a*(T-s)))

def integrand_1(s):
    if a==0:
        ans = math.exp(b*b*gamma*(T-s))*gamma*math.exp(2*a*(T-s))
    else:
        ans = math.exp(b*b*gamma/a*(math.exp(a*(T-s)) - 1))*gamma*math.exp(2*a*(T-s))
    return(ans)

def alpha(s):
    integrated_1 = quad(integrand_1, s, T)[0]
    if a==0:
        w = math.exp(-b*b*gamma*(T-s)) * (1+ b*b*integrated_1)
    else:
        w = math.exp(-b*b*gamma/a*(math.exp(a*(T-s)) - 1)) * (1+ b*b*integrated_1)
    return(gamma*math.exp(2*a*(T-s))/w)

def integrand_2(s):
    return(a+b*b*(beta(s)-alpha(s)))

def gamma_fun(t1, t2):
    integrated_2 = quad(integrand_2, t1, t2)[0]
    return(math.exp(integrated_2))
# here, u is the running time
def integrand_3(u, s):
    return(gamma_fun(u,s)*gamma_fun(u,s))
# here, s is the running time
def integrand_4(s, t, x):
    integrated_3 = quad(integrand_3, t, s, args=(s))[0]
    return((beta(s) - alpha(s))*(beta(s) - alpha(s))*\
        (gamma_fun(t,s)*x*gamma_fun(t,s)*x + sigma*sigma*integrated_3))

def exact_sol_fun(t, x):
    integrated_4 = quad(integrand_4, t, T, args=(t, x))[0]
    integrated_3 = quad(integrand_3, t, T, args=(T))[0]
    return(integrated_4/2 + sigma*sigma*integrated_3*gamma/2)
# TODO: this exact sol is wrong!
# you should have one more term [x*gamma_fun(t,T)]^2-x^2
# now, it won't have any problem because x is 0!

if __name__ == "__main__":
    # 0 for LQ and 1 for quartic
    LQ_or_quartic = 1

    # initial state and initial time
    x0 = 0.0
    t0 = 0.0

    # bounds to ensure boundedness of functions
    Upper_bound_coeff = 2
    Upper_bound_func = Upper_bound_coeff**10
    # Lower_bound = 1 / (x0*100)

    # coefficients: terminal time, risk-free rate,
    #               stock growth rate, stock volatility,
    #               risk aversion rate
    T = 0.1

    # seems like when a not 0, it has problem...
    # for LQ and quartic
    a = 1.0
    b = 1.0
    sigma = 1.0
    gamma = 1.0

    if (LQ_or_quartic == 0):
        filename = "T_" + str(T) + "_quadratic.txt"
    elif (LQ_or_quartic == 1):
        filename = "T_" + str(T) + "_quartic.txt"
    file1 = open(filename,"w")
    file1.write('N,M,tt_now,numerical_sol,runtime\n')
    file1.close()
    if (LQ_or_quartic == 0):
        filename_U = "T_" + str(T) + "quadratic_optimal_U.txt"
    elif (LQ_or_quartic == 1):
        filename_U = "T_" + str(T) + "quartic_optimal_U.txt"
    fileU = open(filename_U,"w")
    fileU.write('N,M,tt_now,state,control,probability\n')


    for iteriter in range(1, 5):
        # iteration from 1 to 5
        # for iter in range(20, 20+1):
        for iter in list(range(1, 11)) + [15, 20]:
            # discretization for N-person differential game
            N = iteriter*5
            # discretization for SDE
            M = iter
            # delta time and delta state
            delta_time = (T-t0)/(N*M)
            delta_state = math.sqrt(delta_time)

            # minimum and maximum of state
            min_space = x0-Upper_bound_coeff*N*M*delta_state
            max_space = x0+Upper_bound_coeff*N*M*delta_state
            # size of state and time
            state_size = 2*N*M+1
            time_size = N*M+1
            # state space and time space
            state_space = np.linspace(min_space, max_space, num=state_size)
            time_space = np.linspace(t0, T, num=time_size)

            # optimal control we are intersted in
            # initialized to be all zeros
            # each row represents different states
            # each column represents different time
            optimal_U = np.zeros((2*N*M+1, N*M*M))
            # initial_guess_U = np.zeros(M*M)
            initial_guess_U = np.zeros(1)
            # NOTE that since optim does not support high dimensional array notation
            # initial_guess_U will be 1-D array
            # the row i represents state_diff from the initial value
            # the column j represents time passed from the initial time
            # hence, this should be (2M-1)*M dimensional
            # here, we adopt the transfomration (i, j)-entry corresponds to
            # i*M+j

            if (LQ_or_quartic == 0):
                # for LQ
                xsquared_usquared_iter = np.zeros((2*N*M+1, N*M*M))
                x_iter = np.zeros((2*N*M+1, N*M*M))
            elif (LQ_or_quartic == 1):
                # for quartic
                xfour_usquared_iter = np.zeros((2*N*M+1, N*M*M))
                xcube_iter = np.zeros((2*N*M+1, N*M*M))
                xsquare_iter = np.zeros((2*N*M+1, N*M*M))
                x_iter = np.zeros((2*N*M+1, N*M*M))

            for tt_N in range(N-1, 0-1, -1):
                start_time = time.time()
                # first update the simulation_information
                # this is to give an idea of where the simulation is
                file2 = open("simulation_information.txt", "w")
                file2.write("N: ")
                file2.write(str(N))
                file2.write("\n")
                file2.write("M: ")
                file2.write(str(M))
                file2.write("\n")
                file2.write("time now: ")
                file2.write(str(time_space[tt_N*M]))
                file2.write("\n")
                file2.close()

                with open(filename, "a") as file1_a:
                    file1_a.write(str(N))
                    file1_a.write(",")
                    file1_a.write(str(M))
                    file1_a.write(",")
                    file1_a.write(str(time_space[tt_N*M]))
                    file1_a.write(",")

                print("N: ", N)
                print("M: ", M)
                print("time now: ", time_space[tt_N*M])
                tt_now = time_space[tt_N*M]
                for statediff_index in range(-tt_N*M, tt_N*M+1):
                    state_index = N*M + statediff_index
                    state_now = state_space[state_index]
                    # -------updated code------------------------------------------------
                    for tt_N_inter in range(M-1, -1, -1):
                        for statediff_index_inter in range(-tt_N_inter, tt_N_inter+1):
                            # note that here, we force the optimal u to be within the bound
                            # (-10, 10), and since BFGS does not support this
                            # we change the method to SLSQP
                            # res = minimize(jrecursion_quicker, initial_guess_U, method='BFGS', \
                            #         args=(tt_N, state_index, tt_N_inter, statediff_index_inter), tol=1e-6, \
                            #         bounds = (-10, 10))
                            bnds = [(-10, 10)]
                            res = minimize(jrecursion_quicker, initial_guess_U, method='SLSQP', \
                                    args=(tt_N, state_index, tt_N_inter, statediff_index_inter), tol=1e-20, \
                                    bounds = bnds)
                            optimal_U[state_index, tt_N*M*M+(tt_N_inter*tt_N_inter + tt_N_inter + statediff_index_inter)] = copy.copy(res.x)
                            initial_guess_U = copy.copy(res.x)
                            jrecursion_fill(tt_N, state_index, tt_N_inter, statediff_index_inter)
                    # res = minimize(jrecursion, initial_guess_U, method='BFGS', \
                    #         args=(tt_N*M, state_index, tt_now, state_now), tol=1e-6)
                    # optimal_U[state_index, (tt_N*M*M):((tt_N+1)*M*M)] = copy.copy(res.x)
                    # initial_guess_U = copy.copy(res.x)
                    if state_index == N*M:
                        numerical_sol = res.fun
                        with open(filename, "a") as file1_a:
                            file1_a.write(str(numerical_sol))
                            file1_a.write(",")
                        print("numerical solution: ", numerical_sol)
                with open(filename, "a") as file1_a:
                    file1_a.write(str(time.time()-start_time))
                    file1_a.write("\n")
                print("time: ", time.time()-start_time)
            if (iteriter == 4) and (iter == 20):
                prob_U = np.zeros((N, 2*N*M+1, M, 2*M+1))
                prob_U[0,N*M,0,M] = 1
                for ii in range(0, N):
                    for kk in range(0, M-1):
                        # the recursion for 0 to M-2
                        for jj in range(-ii*M, ii*M+1):
                            for ll in range(-kk, kk+1):
                                j_index = jj + N*M
                                l_index = ll + M
                                tt_now = time_space[ii*M + kk]
                                x_now = state_space[j_index + ll]
                                u_now = optimal_U[j_index, ii*M*M + kk*kk + kk + ll]
                                prob_now = prob_U[ii, j_index, kk, l_index]
                                fileU.write(str(N))
                                fileU.write(",")
                                fileU.write(str(M))
                                fileU.write(",")
                                fileU.write(str(tt_now))
                                fileU.write(",")
                                fileU.write(str(x_now))
                                fileU.write(",")
                                fileU.write(str(u_now))
                                fileU.write(",")
                                fileU.write(str(prob_now))
                                fileU.write("\n")
                                prob_U[ii, j_index, kk+1, l_index+1] += prob_up(tt_now, x_now, u_now)*prob_now
                                prob_U[ii, j_index, kk+1, l_index-1] += prob_down(tt_now, x_now, u_now)*prob_now
                                prob_U[ii, j_index, kk+1, l_index] += prob_stay(tt_now, x_now, u_now)*prob_now
                    # for the particular case of kk = M-1
                    # when it is M-1, i needs to advance
                    # but when it is at the last point, we do not need to do so!
                    for jj in range(-ii*M, ii*M+1):
                        for ll in range(-(M-1), (M-1)+1):
                            j_index = jj + N*M
                            l_index = ll + M
                            tt_now = time_space[ii*M + M-1]
                            x_now = state_space[j_index + ll]
                            u_now = optimal_U[j_index, ii*M*M + (M-1)*(M-1) + M-1 + ll]
                            prob_now = prob_U[ii, j_index, M-1, l_index]
                            fileU.write(str(N))
                            fileU.write(",")
                            fileU.write(str(M))
                            fileU.write(",")
                            fileU.write(str(tt_now))
                            fileU.write(",")
                            fileU.write(str(x_now))
                            fileU.write(",")
                            fileU.write(str(u_now))
                            fileU.write(",")
                            fileU.write(str(prob_now))
                            fileU.write("\n")
                            if (ii != (N-1)):
                                prob_U[ii+1, j_index+ll+1, 0, M] += prob_up(tt_now, x_now, u_now)*prob_now
                                prob_U[ii+1, j_index+ll-1, 0, M] += prob_down(tt_now, x_now, u_now)*prob_now
                                prob_U[ii+1, j_index+ll, 0, M] += prob_stay(tt_now, x_now, u_now)*prob_now
                # print(optimal_U)
    fileU.close()
