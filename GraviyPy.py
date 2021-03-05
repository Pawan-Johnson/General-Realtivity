from sympy import *

class Metric():

    def DefMetric2Sphere(self):
        self.N=2
        self.gij = symbols("g")
        self.M,self.t,self.r,self.theta,self.phi = symbols(" M t r theta phi")
        self.x0,self.x1 = symbols("x^0 x^1")
        self.x0 = self.theta
        self.x1 = self.phi
        self.x = [self.x0, self.x1]
        g00 = self.r**2
        g01 = 0
        g10 = 0
        g11 = self.r*self.r*sin(self.theta)*sin(self.theta)
        self.gij = Matrix([[g00,g01],
                         [g10,g11]])
        
        self.gijinv = self.gij.inv(method="LU")
    
    def DefMetric3D(self):
        self.N = 3
        self.gij = symbols("g")
        self.M,self.t,self.r,self.theta,self.phi = symbols(" M t r theta phi")
        self.x0,self.x1,self.x2,self.x3 = symbols("x^0 x^1 x^2 x^3")
        self.x0 = self.r
        self.x1 = self.theta
        self.x2 = self.phi
        self.x3 = self.phi
        self.x = [self.x0,self.x1,self.x2,self.x3]

    
        g00 = 1
        g01 = 0
        g02 = 0
        g03 = 0
        g10 = 0
        g11 = self.r*self.r
        g12 = 0
        g13 = 0
        g20 = 0
        g21 = 0
        g22 = self.r**2*sin(self.theta)*sin(self.theta)
        g23 = 0
        g30 = 0
        g31 = 0
        g32 = 0
        g33 = self.r**2*sin(self.theta)**2

        self.gij = Matrix([[g00,g01,g02,g03],
                        [g10,g11,g12,g13],
                        [g20,g21,g22,g23],
                        [g30,g31,g32,g33]])
        self.gijinv = self.gij.inv(method="LU")
    
    def DefMetric4D(self):
        pass

    def DefMetricND(self):
        self.N =2
        self.gij = symbols("g")
        self.M,self.t,self.r,self.theta,self.phi = symbols(" M t r theta phi")
        self.x0,self.x1 = symbols("x^0 x^1")
        self.x0 = self.theta
        self.x1 = self.phi
        self.x = [self.x0, self.x1]
        g00 = -1
        g01 = 0
        g10 = 0
        g11 = -1*exp(-self.r*self.x0)
        self.gij = Matrix([[g00,g01],
                         [g10,g11]])
        
        self.gijinv = self.gij.inv(method="LU")

    def DefMetric(self):
        if (N==2):
            DefMetric2D(self)
        elif (N==3) :
            DefMetric3D(self)
        elif (N==4) :
            DefMetric4D(self)
        else:
            DefMetricND(self) 

    def g(self,i,j):
        return self.gij[i,j]  

    def ginv(self,i,j):
        return self.gijinv[i,j]  
    
    def printMetric(self):
        print("\n\nThe metric is defined as:")
        for i in range(self.N):
            for j in range(self.N):
                print(self.g(i,j),end='\t\t\t')
            print()
    
    def printInvMetric(self):
        print("\n\nThe Inverse metric is defined as:")
        for i in range(self.N):
            for j in range(self.N):
                print(self.ginv(i,j),end='\t\t\t')
            print()
    
    def T(self,m,i,j):
        T1 = 0
        T2 = 0
        T3 = 0
        for k in range(self.N):
            T1 += self.ginv(m,k)*diff(self.g(i,k),self.x[j])  
            T2 += self.ginv(m,k)*diff(self.g(k,j),self.x[i])  
            T3 += self.ginv(m,k)*diff(self.g(j,i),self.x[k]) 
            # print(g(i,k),x[j])
            # print(diff(g(i,k),x[j]))
        return simplify(0.5*(T1 + T2 -T3))   

    def printChristophell(self):
        print("\n\nChristophell Symbols")
        for i in range(self.N):
            print(f"{i}th Christophell Matrix:")
            for j in range(self.N):
                for k in range(self.N):
                    self.T(i,j,k)
                    print(self.T(i,j,k), end='\t\t')
                print()
            print()

    def R(self,i,j,k,l):
        R = diff(self.T(i,l,j),self.x[k]) - diff(self.T(i,k,j),self.x[l])
        TT1 = 0
        TT2 = 0
        for n in range(0,self.N):
            TT1 = TT1 + self.T(i,k,n)*self.T(n,l,j)
            TT2 = TT2 + self.T(i,l,n)*self.T(n,k,j)
        R = R + TT1 - TT2
        return simplify(R)

    def printReimannTensor(self):
        print("\n\nReimann Curvature Tensor")
        for i in range(self.N):
            for j in range(self.N):
                for k in range(self.N):
                    for l in range(self.N):
                        if (self.R(i,j,k,l)!=0):
                            print(f"Component({i},{j},{k},{l}) of Reimann Tensor is:  = {self.R(i,j,k,l)}")
        print("All other components are zero\n\n")

    def covarR(self,i,j,k,l):
        CovR =0
        for n in range(0,self.N):
            CovR += self.g(i,n)*self.R(n,j,k,l)
        return simplify(CovR)

    def printCovariantReimannTensor(self):
        print("\n\nCovariant Reimann Curvature Tensor")
        for i in range(self.N):
            for j in range(self.N):
                for k in range(self.N):
                    for l in range(self.N):
                        if (self.R(i,j,k,l)!=0):
                            print(f"Component({i},{j},{k},{l}) of Covariant Reimann Tensor is:  = {self.R(i,j,k,l)}")
        print("All other components are zero\n\n")

    def RicciTensor(self,i,j):
        RicT =0
        for m in range(0,self.N):
            for n in range(0,self.N):
                RicT += self.ginv(m,n)*self.covarR(m,i,n,j)

        return simplify(RicT)

    def printRicciTesnor(self):
        print("Components of Ricci Tensor")
        for i in range(self.N):
            for j in range(self.N):
                Rval = self.RicciTensor(i,j)
                try:
                    Rval = float(Rval)
                except:
                    pass
                print(Rval,end='\t\t\t')
            print()
    
    def RicciScalar(self):
        RicS = 0
        for m in range(0,self.N):
            for n in range(0,self.N):
                RicS += self.ginv(m,n)*self.RicciTensor(m,n)

        return simplify(RicS)
        
    def PrintRicciScalar(self):
        print("\n\nRicci Scalar")
        print(self.RicciScalar())

    def FullyCovariantR(self,a,b,c,d):
        Rabcd = 0
        i = a
        for j in range(self.N):
            for k in range(self.N):
                for l in range(self.N):
                    Rabcd += self.ginv(b,j)*self.ginv(c,k)*self.ginv(d,l)*self.R(i,j,k,l)
        
        return Rabcd

    def KretschmannScalar(self):
        self.K = 0
        for i in range(self.N):
            for j in range(self.N):
                for k in range(self.N):
                    for l in range(self.N):
                        self.K+=self.covarR(i,j,k,l)*self.FullyCovariantR(i,j,k,l)
        self.K = simplify(self.K)

    def printKretschmannScalar(self):
        self.KretschmannScalar()
        print("\n\nKretschmann Scalar:")
        print(self.K)


                        

    
    
        
    


metric = Metric()
metric.DefMetric2Sphere()
metric.printMetric()
metric.printInvMetric()
metric.printChristophell()
metric.printReimannTensor()
metric.printRicciTesnor()
metric.PrintRicciScalar()
metric.printKretschmannScalar()