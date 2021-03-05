from sympy import *

N = 2
G = symbols("g")
M,t,r,theta,phi = symbols(" M t r theta phi")
x0,x1 = symbols("x^0 x^1 ")
x0 = theta
x1 = phi 
x = [x0,x1]

def Metric():
    g00 = r*r
    g01 = 0
    g10 = 0
    g11 = r*r*sin(theta)*sin(theta)
    return ([[g00,g01],
             [g10,g11]])


def Metric():
    

G = Matrix(Metric())
Ginv = G.inv(method="LU")

def g(i,j):
    return G[i,j]

def ginv(i,j):
    return Ginv[i,j]

print("Metric has been defined")
print("Metric is:")
for i in range(N):
    for j in range(N):
        print(g(i,j),end='\t')
    print()

print("\n\nInverse Metric is:")
for i in range(N):
    for j in range(N):
        print(ginv(i,j),end='\t')
    print()

print("\n\n\nChristophell Symbols have been calculated")
def T(m,i,j):
    T1 = 0
    T2 = 0
    T3 = 0
    for k in range(N):
        T1 += ginv(m,k)*diff(g(i,k),x[j])  
        T2 += ginv(m,k)*diff(g(k,j),x[i])  
        T3 += ginv(m,k)*diff(g(j,i),x[k]) 
        # print(g(i,k),x[j])
        # print(diff(g(i,k),x[j]))
    return simplify(0.5*(T1 + T2 -T3))

for i in range(N):
    for j in range(N):
        for k in range(N):
            T(i,j,k)
            print(T(i,j,k), end='\t')
        print()
    print("\n\n")


def R(i,j,k,l):
    R = diff(T(i,l,j),x[k]) - diff(T(i,k,j),x[l])
    TT1 = 0
    TT2 = 0
    for n in range(0,N):
        TT1 = TT1 + T(i,k,n)*T(n,l,j)
        TT2 = TT2 + T(i,l,n)*T(n,k,j)
    R = R + TT1 - TT2

    return simplify(R)

print("\n\nReimann Curvature Tensor")
for i in range(N):
    for j in range(N):
        for k in range(N):
            for l in range(N):
                if (R(i,j,k,l)!=0):
                    print(f"Component({i},{j},{k},{l}) = {R(i,j,k,l)}")
print("All other components are zero\n\n")

def covarR(i,j,k,l):
    CovR =0
    for n in range(0,N):
        CovR += g(i,n)*R(n,j,k,l)
    return simplify(CovR)

def RicciTensor(i,j):
    RicT =0
    for m in range(0,N):
        for n in range(0,N):
            RicT += ginv(m,n)*covarR(m,i,n,j)

    return simplify(RicT)

def RicciScalar():
    RicS = 0
    for m in range(0,N):
        for n in range(0,N):
            RicS += ginv(m,n)*RicciTensor(m,n)

    return simplify(RicS)

print("Components of Ricci Tensor")
for i in range(N):
    for j in range(N):
        Rval = RicciTensor(i,j)
        try:
            Rval = float(Rval)
        except:
            pass
        print(Rval,end='\t')
    print()

