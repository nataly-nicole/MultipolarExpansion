# coding: utf-8
# Escrito por Nataly Ibarra en noviembre de 2020
# Ante cualquier duda o comentario contactar
# a natalynicole.ibarravera@gmail.com
import itertools
from copy import deepcopy
import sympy as sp
sp.init_printing(use_unicode=True)

# Se representa la permitividad eléctrica del vacío
eps0 = sp.symbols('varepsilon_0', real=True)
eps0

#A continuación se han construido las clases de campos físicos a partir de la clase Tensor disponible en:
#    https://github.com/nataly-nicole/GRpy/blob/master/GRpy/Tensor.py

class formalTensor(object):
    def __init__(self,rank,symbol):
        self.symbol = sp.Symbol(symbol)
        self.rank = rank
    
class Tensor(object):
    '''A class to represent a tensor in a particular basis'''
    def __init__(self,symbol,rank,shape,sym=None,coords=None,**args):
        self.symbol = sp.Symbol(symbol) # The symbol to represent this Tensor
        self.coords = coords # The coordinate system we are using for the representation
        self.shape = shape # The shape of the tensor, for example (-1,1) for k_{a}^{b}
        self.rank = rank # Rank
        self.contr = rank[0]
        self.cov = rank[1]
# We need to know the dimensionality of our spacetime. For the moment
# we will deduce it from the coordinates provide, or if none are given, then
# the assumption will be made that it's stored in the optional arguments
        if coords is not None:
            self.dim = len(self.coords)
        else:
            self.dim = args['dim']
        
        if coords is not None:
            self.allocate(rank)
            self.rep = True
            self.symbolic = formalTensor(rank,symbol)
            
    def __setitem__(self, idx, val):
        self.components[idx] = val
        
    def __getitem__(self, idx):
        return self.components[idx]
    
    def allocate(self, rank):
        '''Allocate the dictionary(hash table) necessary to store the components
        Note that covariant indices are negative! (except for 0 of course)'''
        n = rank[0] + rank[1]
        indc = list(itertools.product(range(self.dim), repeat=n))
        mastr = []
        for i in range(len(indc)):
            temp = []
            for k in range(len(indc[i])):
                if self.shape[k] == -1:
                    temp.append(-indc[i][k])
                else:
                    temp.append(indc[i][k])
        mastr.append(tuple(temp))
        self.components = dict(zip(mastr, [0 for i in range(len(indc))]))
    
    def _dictkeycopy(self, hay):
        keys = hay.keys()
        return dict(zip(keys, [0]*len(keys)))
        
    def getNonZero(self):
        '''Returns only non-zero components of the tensor, if the coordinate system is provided'''
        if self.rep:
            nonzerok = []
            nonzerov = []
            for key in self.components.keys():
                if self.components[key] !=0:
                    nonzerok.append(key)
                    nonzerov.append(self.components[key])
            d = dict(zip(nonzerok,nonzerov))
            keys = d.keys()
            #keys.sort()
            self.nonzero = [(key,d[key]) for key in keys]
            return self.nonzero
        else:
            print("Attempted to get components that have not been initialized!")
    
    def __str__(self):
        '''Print a "nice" human - readable representation of the tensor'''
        self.getNonZero()  
# We will print only non-zero components unless all the components are zero
        ttl=""
        if self.nonzero:
            print(70*'=')
            print('The non-zero components of '+str(self.symbol)+' are:')
            for i in range(len(self.nonzero)):
                ttl = (str(self.nonzero[i][0]) + " : "
                +str(sp.cancel(self.nonzero[i][1])))
                print(ttl)
            print(70*'=')
        else:
            print('All the components of '+str(self.symbol)+' are 0!')
            
class XI(Tensor):
    def __init__(self, coords):
        self.coords = coords
        self.r = sp.sqrt(coords[0]**2 + coords[1]**2 + coords[2]**2)
        super(XI,self).__init__(symbol='x_{i}', rank=(0,1), shape=(-1,), coords=coords)
        for i in range(self.dim):
            self.components[-i] = self.coords[i]
            
class XIn(Tensor):
    def __init__(self, xi, n):
        self.xi = xi
        self.n = n
        super(XIn,self).__init__(symbol='x_{i_1}\cdots x_{i_n} ', rank=(0, self.n), shape=[-i/i for i in range(1,self.n+1)], coords=xi.coords)
        if self.n==0:
            for i1 in range(self.dim):
                self.components[0] = 1
        if self.n==1:
            for i1 in range(self.dim):
                self.components[-i1] = self.xi[-i1]
        if self.n==2:
            for i1 in range(self.dim):
                for i2 in range(self.dim):
                    self.components[-i1,-i2] = xi[-i1]*xi[-i2]
        if self.n==3:
            for i1 in range(self.dim):
                for i2 in range(self.dim):
                    for i3 in range(self.dim):
                        self.components[-i1,-i2,-i3] = xi[-i1]*xi[-i2]*xi[-i3]
        if self.n==4:
            for i1 in range(self.dim):
                for i2 in range(self.dim):
                    for i3 in range(self.dim):
                        for i4 in range(self.dim):
                            self.components[-i1,-i2,-i3,-i4] = xi[-i1]*xi[-i2]*xi[-i3]*xi[-i4]
        if self.n==5:
            for i1 in range(self.dim):
                for i2 in range(self.dim):
                    for i3 in range(self.dim):
                        for i4 in range(self.dim):
                            for i5 in range(self.dim):
                                self.components[-i1,-i2,-i3,-i4,-i5] = xi[-i1]*xi[-i2]*xi[-i3]*xi[-i4]*xi[-i5]
        if self.n==6:
            for i1 in range(self.dim):
                for i2 in range(self.dim):
                    for i3 in range(self.dim):
                        for i4 in range(self.dim):
                            for i5 in range(self.dim):
                                for i6 in range(self.dim):
                                    self.components[-i1,-i2,-i3,-i4,-i5,-i6] = xi[-i1]*xi[-i2]*xi[-i3]*xi[-i4]*xi[-i5]*xi[-i6]
        if self.n==7:
            for i1 in range(self.dim):
                for i2 in range(self.dim):
                    for i3 in range(self.dim):
                        for i4 in range(self.dim):
                            for i5 in range(self.dim):
                                for i6 in range(self.dim):
                                    for i7 in range(self.dim):
                                        self.components[-i1,-i2,-i3,-i4,-i5,-i6,-i7] = xi[-i1]*xi[-i2]*xi[-i3]*xi[-i4]*xi[-i5]*xi[-i6]*xi[-i7]
        if self.n==8:
            for i1 in range(self.dim):
                for i2 in range(self.dim):
                    for i3 in range(self.dim):
                        for i4 in range(self.dim):
                            for i5 in range(self.dim):
                                for i6 in range(self.dim):
                                    for i7 in range(self.dim):
                                        for i8 in range(self.dim):
                                            self.components[-i1,-i2,-i3,-i4,-i5,-i6,-i7,-i8] = xi[-i1]*xi[-i2]*xi[-i3]*xi[-i4]*xi[-i5]*xi[-i6]*xi[-i7]*xi[-i8]
        if self.n==9:
            for i1 in range(self.dim):
                for i2 in range(self.dim):
                    for i3 in range(self.dim):
                        for i4 in range(self.dim):
                            for i5 in range(self.dim):
                                for i6 in range(self.dim):
                                    for i7 in range(self.dim):
                                        for i8 in range(self.dim):
                                            for i9 in range(self.dim):
                                                self.components[-i1,-i2,-i3,-i4,-i5,-i6,-i7,-i8,-i9] = xi[-i1]*xi[-i2]*xi[-i3]*xi[-i4]*xi[-i5]*xi[-i6]*xi[-i7]*xi[-i8]*xi[-i9]                                    
        if self.n==10:
            for i1 in range(self.dim):
                for i2 in range(self.dim):
                    for i3 in range(self.dim):
                        for i4 in range(self.dim):
                            for i5 in range(self.dim):
                                for i6 in range(self.dim):
                                    for i7 in range(self.dim):
                                        for i8 in range(self.dim):
                                            for i9 in range(self.dim):
                                                for i10 in range(self.dim):
                                                    self.components[-i1,-i2,-i3,-i4,-i5,-i6,-i7,-i8,-i9,-i10] = xi[-i1]*xi[-i2]*xi[-i3]*xi[-i4]*xi[-i5]*xi[-i6]*xi[-i7]*xi[-i8]*xi[-i9]*xi[-i10]
        elif self.n>=11:
            print('El orden máximo programado es 10. Favor reduzca el valor de n o programe ordenes superiores.')

class Qin(Tensor):
    def __init__(self, rho, xin, dx, dy, dz):
        self.rho = rho
        self.xin = xin
        self.n = self.xin.n
        self.dx = dx
        self.dy = dy
        self.dz = dz
        super(Qin,self).__init__('Q_{i_{1} \cdots i_{n}}', xin.rank,xin.shape, coords = xin.coords)
#        if self.n==0:
#            self.components[0] = sp.integrate(self.rho,self.dx,self.dy,self.dz)
 #       elif self.n != 0:
        for key in self.xin.components.keys():
            self.components[key] = sp.integrate(self.rho*self.xin.components[key], self.dx, self.dy, self.dz)
            
class Der_partial_n(Tensor):
    def __init__(self, f, xi, n):
        self.xi = xi
        self.n = n
        super(Der_partial_n, self).__init__(symbol='\partial_{i_1}\cdot \partial_{i_n} f ', rank=(0,self.n), shape=[-i/i for i in range(1,self.n+1)], coords=xi.coords)
        if self.n==0:
            self.components[0] = f
        if self.n==1:
            for i1 in range(self.dim):
                self.components[-i1] = sp.diff(f, self.xi[-i1])
#                print i1, self.xi[-i1]
        if self.n==2:
            for i1 in range(self.dim):
                for i2 in range(self.dim):
                    f2 = sp.diff(f, self.xi[-i2])
                    self.components[-i1,-i2] = sp.diff(f2, self.xi[-i1])
        if self.n==3:
            for i1 in range(self.dim):
                for i2 in range(self.dim):
                    for i3 in range(self.dim):
                        f3 = sp.diff(f, self.xi[-i3])
                        f2 = sp.diff(f3, self.xi[-i2])
                        self.components[-i1,-i2,-i3] = sp.diff(f2, self.xi[-i1])
        if self.n==4:
            for i1 in range(self.dim):
                for i2 in range(self.dim):
                    for i3 in range(self.dim):
                        for i4 in range(self.dim):
                            f4 = sp.diff(f, self.xi[-i4])
                            f3 = sp.diff(f4, self.xi[-i3])
                            f2 = sp.diff(f3, self.xi[-i2])
                            self.components[-i1,-i2,-i3,-i4] = sp.diff(f2, self.xi[-i1])
        if self.n==5:
            for i1 in range(self.dim):
                for i2 in range(self.dim):
                    for i3 in range(self.dim):
                        for i4 in range(self.dim):
                            for i5 in range(self.dim):
                                f5 = sp.diff(f, self.xi[-i5])
                                f4 = sp.diff(f5, self.xi[-i4])
                                f3 = sp.diff(f4, self.xi[-i3])
                                f2 = sp.diff(f3, self.xi[-i2])
                                self.components[-i1,-i2,-i3,-i4,-i5] = sp.diff(f2, self.xi[-i1])
        if self.n==6:
            for i1 in range(self.dim):
                for i2 in range(self.dim):
                    for i3 in range(self.dim):
                        for i4 in range(self.dim):
                            for i5 in range(self.dim):
                                for i6 in range(self.dim):
                                    f6 = sp.diff(f, self.xi[-i6])
                                    f5 = sp.diff(f6, self.xi[-i5])
                                    f4 = sp.diff(f5, self.xi[-i4])
                                    f3 = sp.diff(f4, self.xi[-i3])
                                    f2 = sp.diff(f3, self.xi[-i2])
                                    self.components[-i1,-i2,-i3,-i4,-i5,-i6] = sp.diff(f2, self.xi[-i1])
        if self.n==7:
            for i1 in range(self.dim):
                for i2 in range(self.dim):
                    for i3 in range(self.dim):
                        for i4 in range(self.dim):
                            for i5 in range(self.dim):
                                for i6 in range(self.dim):
                                    for i7 in range(self.dim):
                                        f7 = sp.diff(f, self.xi[-i7])
                                        f6 = sp.diff(f7, self.xi[-i6])
                                        f5 = sp.diff(f6, self.xi[-i5])
                                        f4 = sp.diff(f5, self.xi[-i4])
                                        f3 = sp.diff(f4, self.xi[-i3])
                                        f2 = sp.diff(f3, self.xi[-i2])
                                        self.components[-i1,-i2,-i3,-i4,-i5,-i6,-i7] = sp.diff(f2, self.xi[-i1])
        if self.n==8:
            for i1 in range(self.dim):
                for i2 in range(self.dim):
                    for i3 in range(self.dim):
                        for i4 in range(self.dim):
                            for i5 in range(self.dim):
                                for i6 in range(self.dim):
                                    for i7 in range(self.dim):
                                        for i8 in range(self.dim):
                                            f8 = sp.diff(f, self.xi[-i8])
                                            f7 = sp.diff(f8, self.xi[-i7])
                                            f6 = sp.diff(f7, self.xi[-i6])
                                            f5 = sp.diff(f6, self.xi[-i5])
                                            f4 = sp.diff(f5, self.xi[-i4])
                                            f3 = sp.diff(f4, self.xi[-i3])
                                            f2 = sp.diff(f3, self.xi[-i2])
                                            self.components[-i1,-i2,-i3,-i4,-i5,-i6,-i7,-i8] = sp.diff(f2, self.xi[-i1])
        if self.n==9:
            for i1 in range(self.dim):
                for i2 in range(self.dim):
                    for i3 in range(self.dim):
                        for i4 in range(self.dim):
                            for i5 in range(self.dim):
                                for i6 in range(self.dim):
                                    for i7 in range(self.dim):
                                        for i8 in range(self.dim):
                                            for i9 in range(self.dim):
                                                f9 = sp.diff(f, self.xi[-i9])
                                                f8 = sp.diff(f9, self.xi[-i8])
                                                f7 = sp.diff(f8, self.xi[-i7])
                                                f6 = sp.diff(f7, self.xi[-i6])
                                                f5 = sp.diff(f6, self.xi[-i5])
                                                f4 = sp.diff(f5, self.xi[-i4])
                                                f3 = sp.diff(f4, self.xi[-i3])
                                                f2 = sp.diff(f3, self.xi[-i2])
                                                self.components[-i1,-i2,-i3,-i4,-i5,-i6,-i7,-i8,-i9] = sp.diff(f2, self.xi[-i1])
        if self.n==10:
            for i1 in range(self.dim):
                for i2 in range(self.dim):
                    for i3 in range(self.dim):
                        for i4 in range(self.dim):
                            for i5 in range(self.dim):
                                for i6 in range(self.dim):
                                    for i7 in range(self.dim):
                                        for i8 in range(self.dim):
                                            for i9 in range(self.dim):
                                                for i10 in range(self.dim):
                                                    f10 = sp.diff(f, self.xi[-i10])
                                                    f9 = sp.diff(f10, self.xi[-i9])
                                                    f8 = sp.diff(f9, self.xi[-i8])
                                                    f7 = sp.diff(f8, self.xi[-i7])
                                                    f6 = sp.diff(f7, self.xi[-i6])
                                                    f5 = sp.diff(f6, self.xi[-i5])
                                                    f4 = sp.diff(f5, self.xi[-i4])
                                                    f3 = sp.diff(f4, self.xi[-i3])
                                                    f2 = sp.diff(f3, self.xi[-i2])
                                                    self.components[-i1,-i2,-i3,-i4,-i5,-i6,-i7,-i8,-i9,-i10] = sp.diff(f2, self.xi[-i1])
        elif self.n>=11:
            print('El orden máximo programado es 10. Favor reduzca el valor de n o programe ordenes superiores.')

class Potencial_Electrico_n(Tensor):
    def __init__(self, rho, xi, dx, dy, dz, n):
        xin = XIn(xi, n)
        qin = Qin(rho, xin, dx, dy, dz)
        der_n = Der_partial_n(1/xi.r, xi, n)
        super(Potencial_Electrico_n, self).__init__(symbol='\phi_{n} ', rank=(0,0),shape=[], coords=xi.coords)
        if n in range(11):
            suma = 0
#            for  key in list(qin.components)[:-1]:
            for  key in list(qin.components):
                suma += qin.components[key]*der_n.components[key]
            potn = (1/(4*sp.pi*eps0))*((-1)**n/sp.factorial(n))*suma   
        elif n>=11:
            print('El orden máximo programado es 10. Favor reduzca el valor de n o programe ordenes superiores.')
        self.components[0] = potn
            
class Potencial_Electrico(Tensor):
    def __init__(self, rho, xi, dx, dy, dz, nmax):
        self.rho = rho
        self.xi = xi
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.nmax = nmax
        super(Potencial_Electrico, self).__init__(symbol='\phi', rank=(0,0), shape=[], coords=xi.coords)
        if nmax in range(0,11):
            suma = 0
            for n in range(nmax+1):
                potn = Potencial_Electrico_n(rho, xi, dx, dy, dz, n)
                suma += potn.components[0]
            self.components[0] = suma
        elif nmax>=11:
            print('El orden máximo programado es 10. Favor reduzca el valor de n o programe ordenes superiores.')

class Campo_Electrico_n(Tensor):
    def __init__(self, rho, xi, dx, dy, dz, n):
        self.rho = rho
        self.xi = xi
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.n = n
        Potn = Potencial_Electrico_n(rho,xi,dx,dy,dz,n)
        super(Campo_Electrico_n, self).__init__(symbol='E_{i}^{n}', rank=(0,1), shape=(-1,), coords=xi.coords)
        if self.n in range(0,11):
            for i in range(self.dim):
                self.components[-i] = -sp.diff(Potn.components[0], xi[-i])
        elif self.n>=11:
            print('El orden máximo programado es 10. Favor reduzca el valor de n o programe ordenes superiores.')

class Campo_Electrico(Tensor):
    def __init__(self, rho, xi, dx, dy, dz, nmax):
        self.rho = rho
        self.xi = xi
        self.dx = dx
        self.dy = dy
        self.dz = dz
        self.nmax = nmax
        super(Campo_Electrico, self).__init__(symbol='E_{i}', rank=(0,1), shape=(-1,), coords=xi.coords)
        if nmax in range(0,11):
            sumax, sumay, sumaz = 0, 0, 0
            for n in range(nmax+1):
                Ein = Campo_Electrico_n(rho, xi, dx, dy, dz, n)
                sumax += Ein.components[0]
                sumay += Ein.components[-1]
                sumaz += Ein.components[-2]
            self.components[0] = sumax
            self.components[-1] = sumay
            self.components[-2] = sumaz
        elif nmax>=11:
            print('El orden máximo programado es 10. Favor reduzca el valor de n o programe ordenes superiores.')
            
def save_campos_n(rho, xi, dx, dy, dz, n):
    name_file = "./campos_n/campos_n_"+str(n)+".txt"
    Pot_Elec_n = Potencial_Electrico_n(rho, xi, dx, dy, dz, n)
    Campo_Elec_n = Campo_Electrico_n(rho, xi, dx, dy, dz, n)
    campos = [sp.simplify(Pot_Elec_n.components[0]), sp.simplify(Campo_Elec_n.components[0]), 
              sp.simplify(Campo_Elec_n.components[-1]), sp.simplify(Campo_Elec_n.components[-2])]
    file = open(name_file, "w")
    file.write(sp.srepr(campos))
    file.close()
    
def save_campos_nmax(rho, xi, dx, dy, dz, nmax):
    name_file = "./campos_nmax/campos_nmax_"+str(nmax)+".txt"
    Pot_Elec = Potencial_Electrico(rho, xi, dx, dy, dz, nmax)
    Campo_Elec = Campo_Electrico(rho, xi, dx, dy, dz, nmax)
    campos = [sp.simplify(Pot_Elec.components[0]), sp.simplify(Campo_Elec.components[0]),
              sp.simplify(Campo_Elec.components[-1]), sp.simplify(Campo_Elec.components[-2])]
    file = open(name_file, "w")
    file.write(sp.srepr(campos))
    file.close()
    
def read_campos_n(n):
    name_file = "./campos_n/campos_n_"+str(n)+".txt"
    file = open(name_file,"r").read()
    campos_n = sp.sympify(file)
    return campos_n

def read_campos_nmax(nmax):
    name_file = "./campos_nmax/campos_nmax_"+str(nmax)+".txt"
    file = open(name_file,"r").read()
    campos_nmax = sp.sympify(file)
    return campos_nmax
