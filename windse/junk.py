import dolfin

print(dir(dolfin))





# test = "A"

# class A(object):
#     def __init__(self):
#         self.arg = "hey its A"
    

# class B(object):
#     def __init__(self):
#         self.arg = "hey its B"

# class AorB(A if test=="A" else B
#     ):
#     def __init__(self):
#         super(AorB, self).__init__()


# Test_Class = AorB()

# print(Test_Class.arg)




# test = np.float64(1.03)
# test = 4.5
# print(type(test))
# print(test)
# print(float(test))

# N = 100
# x = np.arange(N)
# y = x**2
# from scipy.interpolate import interp1d
# iii = interp1d(x, y, fill_value=(-10, 10), bounds_error=False)
# print(iii(-1))
# print(iii(101))

# exit()

# for t in test:
#     print(t)

# exit()


# test2=test.get("huh",False)

# print(test2)

# print(3*(3,))

# boundary_names = {"inflow":1,"outflow":2,"top":3,"bottom":4}
# boundary_types = {"inflow":    ["inflow","top"],
#                   "no_slip":   ["bottom"],
#                   "no_stress": ["outflow"]}

# bcs_eqns = []
# for bc_type, bs in boundary_types.items():
#     print(bc_type)
#     if bc_type == "inflow":
#         for b in bs:
#             bcs_eqns.append(["INFLOW_EQ",boundary_names[b]])

#     elif bc_type == "no_slip":
#         for b in bs:
#             bcs_eqns.append(["zeros",boundary_names[b]])


#     elif bc_type == "no_stress":
#         for b in bs:
#             bcs_eqns.append([None,boundary_names[b]])

#     else:
#         raise ValueError(bc_type+"is not a recognized boundary type")


# print(bcs_eqns[3][0])