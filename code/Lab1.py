from random import *
from prettytable import PrettyTable

Table_1 = PrettyTable()

a = [randint(0, 20) for i in range(4)]

Table_1.add_column("Number", [i for i in range(4)])
Table_1.add_column("a", a)

def x_generator():
    return [randint(0, 20) for i in range(8)]

X1, X2, X3 = [x_generator() for k in range(3)]

def y_calculator(X1, X2, X3):
    return a[0] + a[1] * X1 + a[2] * X2 + a[3] * X3

Y = [y_calculator(X1[i], X2[i], X3[i]) for i in range(8)]
Y_min = min(Y)

def x0_calculator(X):
    return (max(X) + min(X))/2

X0 = [x0_calculator(i) for i in [X1, X2, X3]]

def dx_calculator(X0, X):
    return X0 - min(X)

x = [X1, X2, X3]
DX = [dx_calculator(X0[i], x[i]) for i in range(3)]

def xn_calculator(X0, DX, x):
    return [round(((i - X0) / DX),3) for i in x]

Xn = [xn_calculator(X0[i], DX[i], x[i]) for i in range(3)]

def Y_average(Y):
    res = 0
    for i in Y:
        res += i
    return res/len(Y)

Y_avg = Y_average(Y)

def Y_optimal(Y_avg, Y):
    opt = []
    for i in range(8):
        opt.append(Y[i] - Y_avg)
    return opt

Y_opt = Y_optimal(Y_avg,Y)

def index(optimal):
    return max((a,i) for i, a in enumerate(Y_opt) if a < 0)[1]

opt_index = index(Y_opt)
opt_result = [X1[opt_index], X2[opt_index], X3[opt_index]]

Table_2 = PrettyTable()
Table_2.add_column("Number", [i for i in range(1, 9)])
Table_2.add_column("X1", X1)
Table_2.add_column("X2", X2)
Table_2.add_column("X3", X3)
Table_2.add_column("Y", Y)
Table_2.add_column("Y_opt", Y_opt)
Table_2.add_column("Xn1", Xn[0])
Table_2.add_column("Xn2", Xn[1])
Table_2.add_column("Xn3", Xn[2])


Table_3 = PrettyTable()
Table_3.field_names = ["Variable", "Value_1", "Value_2", "Value_3"]
Table_3.add_row(["X0", X0[0], X0[1], X0[2]])
Table_3.add_row(["DX", DX[0], DX[1], DX[2]])

print(Table_1)
print(Table_2)
print(Table_3)

print(f"\nЕталонне значення функції: Y = {a[0]} + {a[1]}X0[0] + {a[2]}X0[1] + {a[3]}X0[2]" )
print(f"Функція: Y = {a[0]} + {a[1]}X1 + {a[2]}X2 + {a[3]}X3")
print("Критерій оптимальності: -> Yсереднє")
print("Оптимальна точка плану:  Y({0}, {1}, {2}) = {3}".format(*opt_result, "%.1f" % Y[opt_index]))

print("\nДодаткове завдання(Варіант 217):")
print(f"min(Y):{Y_min}")
