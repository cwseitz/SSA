from ssa._ssa import two_state

end_time = 24.0
k_on = 1.0
k_off = 1.0
time_step = 0.005
x1, x2 = two_state([end_time,time_step,k_on,k_off])
print(x1,x2)
