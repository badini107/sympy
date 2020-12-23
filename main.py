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
        self.angular_velocity = []
        self.angular_acc = []
        for i, _ in enumerate(self.links):
            self.angles.append(Function('theta'+ str(i))(self.t))
            self.angular_velocity.append(Function('theta'+ str(i) + "'")(self.t))
            self.angular_acc.append(Function('theta'+ str(i) + "''")(self.t))
        
        self.x = []
        self.y = []      
        for i, link in enumerate(self.links):
            self.x.append(1 * sin(self.angles[i]) * link.length)
            self.y.append(-1 * cos(self.angles[i]) * link.length)
            for j in range(i):
                self.x[i] += sin(self.angles[j]) * self.links[j].length
                self.y[i] += -cos(self.angles[j]) * self.links[j].length
    
    def derivate(self, expression, derivator):
        d_exp = []
        for i, exp in enumerate(expression):
            exp = diff(exp, derivator)
            for angle, ang_vel, ang_acc in zip(self.angles, self.angular_velocity, self.angular_acc):
                exp = exp.subs(Derivative(angle, self.t), ang_vel)
                exp = exp.subs(Derivative(ang_vel, self.t), ang_acc)
            d_exp.append(exp)
        return d_exp

    def addtion(self, expression_1, expression_2):
        exp_vel  = []
        for vel_1, vel_2 in zip(expression_1, expression_2):
            exp_vel.append(simplify(vel_1 + vel_2))

        return exp_vel
    
    def square(self, expression):
        exp_squared = []
        for exp in expression:
            exp_squared.append(exp**2)
        return  exp_squared
    
    def kinetic(self, velocity):
        self.kinetic_energy = [0]
        for i, vel in enumerate(velocity):
            self.kinetic_energy[0] += vel * self.links[i].mass * 0.5
    
    def potential(self):
        self.potential_energy = [0]
        for i, link in enumerate(self.links):
            self.potential_energy[0] += -1 * link.mass * self.g * link.length * cos(self.angles[i])
            for j in range(i):
                self.potential_energy[0] += -self.links[j].mass * self.g * self.links[j].length * cos(self.angles[j])
    
    def lagrang(self):
        self.L = [self.kinetic_energy[0] - self.potential_energy[0]]

links = []
links.append(Link(0, 0, 1, 0))
links.append(Link(1, 0, 2, 0))
lag = Lagrangian(links)
lag.position()
print('x: ', lag.x)
vx = lag.derivate(lag.x, 't')
print('vx: ',vx)
vy = lag.derivate(lag.y, 't')
v2x = lag.square(vx)
v2y = lag.square(vy)
print('v2x: ',v2x)
vel2 = lag.addtion(v2x, v2y)
print('vel2: ',vel2)
lag.kinetic(vel2)
lag.potential()
print('KE: ',lag.kinetic_energy)
print('PE: ',lag.potential_energy)
lag.lagrang()
print('L: ',lag.L)
dL_dv = lag.derivate(lag.L, lag.angular_velocity[0])
print('dL_dv: ',dL_dv)
dL_dt = [simplify(lag.derivate(dL_dv, 't')[0])]
print('dL_dt: ',dL_dt)
dL_da = lag.derivate(lag.L, lag.angles[0])
print('dL_da: ',dL_da)