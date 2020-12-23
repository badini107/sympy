from sympy import *

class Link():
    def __init__(self, x0, y0, x1, y1):
        self.x0 = x0
        self.y0 = y0
        self.x1 = x1
        self.y1 = y1
        self.length =  sqrt((x1 - x0)**2 + (y1 - y0)**2)
        self.mass = 1

class Lagrangian():
    def __init__(self, links):
        self.links = links
        self.g = 10
        self.t = Symbol('t')
    
    def position(self):
        self.angles = []
        for i, _ in enumerate(self.links):
            self.angles.append(Function('theta'+ str(i))(self.t))
        
        self.x = []
        self.y = []      
        for i, link in enumerate(self.links):
            self.x.append(0.5 * sin(self.angles[i]) * link.length)
            self.y.append(-0.5 * cos(self.angles[i]) * link.length)
            for j in range(i):
                self.x[i] += sin(self.angles[j]) * self.links[j].length
                self.y[i] += -cos(self.angles[j]) * self.links[j].length
    
    def derivate(self, expression, derivator):
        d_exp = []
        for exp in expression:
            d_exp.append(diff(exp, derivator))
        return d_exp


links = []
links.append(Link(0, 0, 1, 0))
links.append(Link(1, 0, 2, 0))
lag = Lagrangian(links)
lag.position()
print(lag.x)
vx = lag.derivate(lag.x, 't')
print(vx)
vy = lag.derivate(lag.y, 't')