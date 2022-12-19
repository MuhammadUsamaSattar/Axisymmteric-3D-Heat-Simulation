import math
import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy

#Sets global parameters
p_wort = 1000
k_wort = 0.6
c_wort = 4184

p_air = 1.2
k_air = 0.025
c_air = 700

R = 0.2
Z = 1

p = [p_wort, p_air]
k = [k_wort, k_air]
c = [c_wort, c_air]

#node class that contains all the values of the node
class node():
    def __init__(self, r, z, T, p, k, c, delta_t, delta_h):
        self.r = r
        self.z = z
        self.T = T
        self.p = p
        self.k = k
        self.c = c
        self.sigma = (self.k/(self.p*self.c))*(delta_t/(delta_h**2))

#grid class that forms a grid of nodes and initiliazes according to set scheme (explicit or implicit)
class grid():
    #Initializes the grid with input parameters
    def __init__(self, boundary_conditions, type, delta_t, delta_h):
        self.nodes = []

        self.delta_t = delta_t
        self.delta_h = delta_h

        self.r_total = int((R/delta_h) + 1)
        self.z_total = int((Z/delta_h) + 1)

        self.boundary_conditions = boundary_conditions
        self.type = type

        self.wort_level = 0.09/(math.pi*R*R)
        self.wort_level = delta_h * round(self.wort_level/delta_h)

        #Generates nodes and sets values for parameters of each node
        for j in range(self.z_total):
            self.nodes.append([])

            choice = 0
            if (j*delta_h/Z) > self.wort_level:
                choice = 1

            for i in range(self.r_total):
                self.nodes[-1].append(node(i*delta_h, j*delta_h, self.boundary_conditions[0][0], p[choice], k[choice], c[choice], delta_t, delta_h))

        #Sets the Dirichlet conditions. Sets the temperature of nodes at certain boundary to some value.
        for i in range(len(self.boundary_conditions[1:])):
            if self.boundary_conditions[i+1][0] == "D":
                if i == 0:
                    for j in range(self.r_total):
                        self.nodes[j][0].T = self.boundary_conditions[i+1][1]

                elif i == 1:
                    for j in range(self.z_total):
                        self.nodes[j][self.r_total-1].T = self.boundary_conditions[i+1][1]

                elif i == 2:
                    for j in range(self.z_total):
                        self.nodes[0][j].T = self.boundary_conditions[i+1][1]

                elif i ==3:
                    for j in range(self.r_total):
                        self.nodes[self.z_total-1][j].T = self.boundary_conditions[i+1][1]

        self.generate_matrix()

    #Generates the factor matrix. The factor matrix does not change at any time so we need to generate it once only
    def generate_matrix(self):
        factor_matrix = []
        constant_matrix = [0] * ((self.r_total*self.z_total))

        for i in range((self.r_total * self.z_total)):
            factor_matrix.append([0]*((self.r_total*self.z_total)))

            node = self.nodes[math.floor(i/self.r_total)][i%self.r_total]

            if math.floor(i/self.r_total) == 0:
                k_z_minus_1 = node.k
                k_z_plus_1 = self.nodes[math.floor(i/self.r_total)+1][i%self.r_total].k

            elif math.floor(i/self.r_total) == (self.z_total-1):
                k_z_plus_1 = node.k
                k_z_minus_1 = self.nodes[math.floor(i/self.r_total)-1][i%self.r_total].k

            else:
                k_z_plus_1 = self.nodes[math.floor(i/self.r_total)+1][i%self.r_total].k
                k_z_minus_1 = self.nodes[math.floor(i/self.r_total)-1][i%self.r_total].k

            if node.r != 0:
                A = node.sigma * ((node.r+delta_h+node.r)/(2*node.r))
                C = node.sigma * ((node.r-delta_h+node.r)/(2*node.r))
            else:
                A = node.sigma
                C = node.sigma

            B = node.sigma * ((k_z_plus_1+node.k)/(2*node.k))
            D = node.sigma * ((k_z_minus_1+node.k)/(2*node.k))
            E = (1 - 2*node.sigma - ((node.sigma*(k_z_plus_1+k_z_minus_1+(2*node.k))/(2*node.k))))

            if self.type == "explicit":
                if (i%self.r_total) == 0 or ((i+1)%self.r_total) == 0 or (math.floor(i/self.r_total)) == 0 or (math.floor(i/self.r_total)+1) == self.z_total:
                    factor_matrix[-1][i] = E

                    if (i%self.r_total) == 0:
                        if self.boundary_conditions[1][0] == "N":
                            if self.boundary_conditions[1][1] == 0:
                                f = 0
                            else:
                                f = 2*self.delta_h*self.boundary_conditions[1][1]*(self.boundary_conditions[1][2]-node.T)/node.k
                            
                            constant_matrix[i] += (-C*f)
                            factor_matrix[-1][i+1] = A+C

                    elif ((i+1)%self.r_total) == 0:
                        if self.boundary_conditions[2][0] == "N":
                            if self.boundary_conditions[2][1] == 0:
                                f = 0
                            else:
                                f = 2*self.delta_h*self.boundary_conditions[2][1]*(self.boundary_conditions[2][2]-node.T)/node.k
                            
                            constant_matrix[i] += (A*f)
                            factor_matrix[-1][i-1] = A+C

                    else:
                        factor_matrix[-1][i+1] = A
                        factor_matrix[-1][i-1] = C

                    if (math.floor(i/self.r_total)) == 0:
                        if self.boundary_conditions[3][0] == "N":
                            if self.boundary_conditions[3][1] == 0:
                                f = 0
                            else:
                                f = 2*self.delta_h*self.boundary_conditions[3][1]*(self.boundary_conditions[3][2]-node.T)/node.k
                            
                            constant_matrix[i] += (-D*f)
                            factor_matrix[-1][i+self.r_total] = B+D

                    elif (math.floor(i/self.r_total)+1) == self.z_total:
                        if self.boundary_conditions[4][0] == "N":
                            if self.boundary_conditions[4][1] == 0:
                                f = 0
                            else:
                                f = 2*self.delta_h*self.boundary_conditions[4][1]*(self.boundary_conditions[4][2]-node.T)/node.k
                            
                            constant_matrix[i] += (B*f)
                            factor_matrix[-1][i-self.r_total] = B+D

                    else:
                        factor_matrix[-1][i+self.r_total] = B
                        factor_matrix[-1][i-self.r_total] = D

                    if (i%self.r_total) == 0:
                        if self.boundary_conditions[1][0] == "D":
                            factor_matrix[-1] = [0]*self.r_total*self.z_total
                            factor_matrix[-1][i] == 1
                            constant_matrix[i] = 0

                    elif ((i+1)%self.r_total) == 0:
                        if self.boundary_conditions[2][0] == "D":
                            factor_matrix[-1] = [0]*(self.r_total*self.z_total)
                            factor_matrix[-1][i] = 1
                            constant_matrix[i] = 0

                    if (math.floor(i/self.r_total)) == 0:
                        if self.boundary_conditions[3][0] == "D":
                            factor_matrix[-1] = [0]*self.r_total*self.z_total
                            factor_matrix[-1][i] = 1
                            constant_matrix[i] = 0

                    elif (math.floor(i/self.r_total)+1) == self.z_total:
                        if self.boundary_conditions[4][0] == "D":
                            factor_matrix[-1] = [0]*self.r_total*self.z_total
                            factor_matrix[-1][i] = 1
                            constant_matrix[i] = 0
                else:
                    factor_matrix[-1][i] = E
                    factor_matrix[-1][i+1] = A
                    factor_matrix[-1][i-1] = C
                    factor_matrix[-1][i+self.r_total] = B
                    factor_matrix[-1][i-self.r_total] = D
                    constant_matrix.append(0)

            elif self.type == "implicit":
                if (i%self.r_total) == 0 or ((i+1)%self.r_total) == 0 or (math.floor(i/self.r_total)) == 0 or (math.floor(i/self.r_total)+1) == self.z_total:
                    factor_matrix[-1][i] = -E+2

                    if (i%self.r_total) == 0:
                        if self.boundary_conditions[1][0] == "N":
                            factor_matrix[-1][i+1] = (A+C)*(-1)

                    elif ((i+1)%self.r_total) == 0:
                        if self.boundary_conditions[2][0] == "N":
                            factor_matrix[-1][i-1] = (A+C)*(-1)

                    else:
                        factor_matrix[-1][i+1] = A*(-1)
                        factor_matrix[-1][i-1] = C*(-1)

                    if (math.floor(i/self.r_total)) == 0:
                        if self.boundary_conditions[3][0] == "N":
                            factor_matrix[-1][i+self.r_total] = (B+D)*(-1)

                    elif (math.floor(i/self.r_total)+1) == self.z_total:
                        if self.boundary_conditions[4][0] == "N":
                            factor_matrix[-1][i-self.r_total] = (B+D)*(-1)

                    else:
                        factor_matrix[-1][i+self.r_total] = B*(-1)
                        factor_matrix[-1][i-self.r_total] = D*(-1)

                    if (i%self.r_total) == 0:
                        if self.boundary_conditions[1][0] == "D":
                            factor_matrix[-1] = [0]*self.r_total*self.z_total
                            factor_matrix[-1][i] == 1

                    elif ((i+1)%self.r_total) == 0:
                        if self.boundary_conditions[2][0] == "D":
                            factor_matrix[-1] = [0]*(self.r_total*self.z_total)
                            factor_matrix[-1][i] = 1

                    if (math.floor(i/self.r_total)) == 0:
                        if self.boundary_conditions[3][0] == "D":
                            factor_matrix[-1] = [0]*self.r_total*self.z_total
                            factor_matrix[-1][i] = 1

                    elif (math.floor(i/self.r_total)+1) == self.z_total:
                        if self.boundary_conditions[4][0] == "D":
                            factor_matrix[-1] = [0]*self.r_total*self.z_total
                            factor_matrix[-1][i] = 1
                else:
                    factor_matrix[-1][i] = -E+2
                    factor_matrix[-1][i+1] = A*(-1)
                    factor_matrix[-1][i-1] = C*(-1)
                    factor_matrix[-1][i+self.r_total] = B*(-1)
                    factor_matrix[-1][i-self.r_total] = D*(-1)

        self.factor_matrix = factor_matrix
        self.constant_matrix = constant_matrix
       
    #Evaluates the flux at a boundary. This is an evaluation of Nuemann boundary condition. Since the flux is dependent on temperature difference of air and outer wall of tank, this changes
    #and has to be evaluated at every time step
    def evaluate_flux(self):
        constant_matrix = [0] * ((self.r_total*self.z_total))

        for i in range((self.r_total * self.z_total)):

            node = self.nodes[math.floor(i/self.r_total)][i%self.r_total]

            if math.floor(i/self.r_total) == 0:
                k_z_minus_1 = node.k
                k_z_plus_1 = self.nodes[math.floor(i/self.r_total)+1][i%self.r_total].k

            elif math.floor(i/self.r_total) == (self.z_total-1):
                k_z_plus_1 = node.k
                k_z_minus_1 = self.nodes[math.floor(i/self.r_total)-1][i%self.r_total].k

            else:
                k_z_plus_1 = self.nodes[math.floor(i/self.r_total)+1][i%self.r_total].k
                k_z_minus_1 = self.nodes[math.floor(i/self.r_total)-1][i%self.r_total].k

            if node.r != 0:
                A = node.sigma * ((node.r+delta_h+node.r)/(2*node.r))
                C = node.sigma * ((node.r-delta_h+node.r)/(2*node.r))
            else:
                A = node.sigma
                C = node.sigma

            B = node.sigma * ((k_z_plus_1+node.k)/(2*node.k))
            D = node.sigma * ((k_z_minus_1+node.k)/(2*node.k))
            E = (1 - 2*node.sigma - ((node.sigma*(k_z_plus_1+k_z_minus_1+(2*node.k))/(2*node.k))))

            if (i%self.r_total) == 0 or ((i+1)%self.r_total) == 0 or (math.floor(i/self.r_total)) == 0 or (math.floor(i/self.r_total)+1) == self.z_total:
                if (i%self.r_total) == 0:
                    if self.boundary_conditions[1][0] == "N":
                        if self.boundary_conditions[1][1] == 0:
                            f = 0
                        else:
                            f = 2*self.delta_h*self.boundary_conditions[1][1]*(self.boundary_conditions[1][2]-node.T)/node.k
                        
                        constant_matrix[i] += (-C*f)

                elif ((i+1)%self.r_total) == 0:
                    if self.boundary_conditions[2][0] == "N":
                        if self.boundary_conditions[2][1] == 0:
                            f = 0
                        else:
                            f = 2*self.delta_h*self.boundary_conditions[2][1]*(self.boundary_conditions[2][2]-node.T)/node.k
                        
                        constant_matrix[i] += (A*f)

                if (math.floor(i/self.r_total)) == 0:
                    if self.boundary_conditions[3][0] == "N":
                        if self.boundary_conditions[3][1] == 0:
                            f = 0
                        else:
                            f = 2*self.delta_h*self.boundary_conditions[3][1]*(self.boundary_conditions[3][2]-node.T)/node.k
                        
                        constant_matrix[i] += (-D*f)

                elif (math.floor(i/self.r_total)+1) == self.z_total:
                    if self.boundary_conditions[4][0] == "N":
                        if self.boundary_conditions[4][1] == 0:
                            f = 0
                        else:
                            f = 2*self.delta_h*self.boundary_conditions[4][1]*(self.boundary_conditions[4][2]-node.T)/node.k
                        
                        constant_matrix[i] += (B*f)
                
                if (i%self.r_total) == 0:
                    if self.boundary_conditions[1][0] == "D":
                        constant_matrix[i] = 0

                elif ((i+1)%self.r_total) == 0:
                    if self.boundary_conditions[2][0] == "D":
                        constant_matrix[i] = 0

                if (math.floor(i/self.r_total)) == 0:
                    if self.boundary_conditions[3][0] == "D":
                        constant_matrix[i] = 0

                elif (math.floor(i/self.r_total)+1) == self.z_total:
                    if self.boundary_conditions[4][0] == "D":
                        constant_matrix[i] = 0

        self.constant_matrix = deepcopy(constant_matrix)

    #Updates the Temperature value of nodes according to the scheme.
    def update(self):
        T_matrix = []

        for i in range(self.r_total * self.z_total):
            node = self.nodes[math.floor(i/self.r_total)][i%self.r_total]
            T_matrix.append(node.T)

        if self.type == "explicit":
            self.evaluate_flux()
            T_matrix = np.matmul(self.factor_matrix, T_matrix)
            T_matrix = np.add(T_matrix, self.constant_matrix)
        elif self.type == "implicit":
            T_matrix = np.linalg.solve(self.factor_matrix, T_matrix)

        for i in range(len(T_matrix)):
            self.nodes[math.floor(i/self.r_total)][i%self.r_total].T = T_matrix[i]

    #Displays the Temperature values of the nodes
    def display(self):
        T_matrix = []

        for row in self.nodes:
            T_matrix.append([])
            for node in row:
                T_matrix[-1].append(round(node.T,3))
        for i in range(len(T_matrix)):
            print(T_matrix[len(T_matrix)-1-i])

        print("--------------")

    #Returns the Temperature values of the nodes
    def get_T_matrix(self):
        T_matrix = []

        for row in self.nodes:
            T_matrix.append([])
            for node in row:
                T_matrix[-1].append(node.T)

        temp = []
        for i in range(len(T_matrix)):
            temp.append(T_matrix[len(T_matrix)-1-i])

        return temp

    #Returns the Temperature value at thermocouple location
    def get_T_thermocouple(self):
        r = 0
        z = math.floor(0.25/delta_h)

        return self.nodes[z][r].T

    #Returns the r, z and T value of a node for countour plot generation
    def get_contour_data(self):
        r_matrix = []
        z_matrix = []
        T_matrix = []

        for row in self.nodes:
            r_matrix.append([])
            z_matrix.append([])
            T_matrix.append([])

            for node in row:
                r_matrix[-1].append(node.r)
                z_matrix[-1].append(node.z)
                T_matrix[-1].append(node.T)

        return [r_matrix, z_matrix, T_matrix]

#########################################################################################################################################################################################
#MAIN BODY OF PROGRAM
#########################################################################################################################################################################################
#Evaluates the problem for explicit scheme
####################################################################
#Set the program parameters from here

delta_h_all = 0.05
delta_t_explicit = 6
delta_t_implicit = 300

delta_h_alternate = 0.05
delta_t_alternate = 6

delta_h_covection = 0.05
delta_t_covection =  1

boundary_conditions_all = [[15], ["N", 0], ["D", 35], ["N", 0], ["N", 0]] #[[t=0 condition], [Type at r=0, value], [Type at r=R, value], [Type at z=0, value], [Type at z=H, value]]
boundary_conditions_convection = [[15], ["N", 0], ["N", 15, 25], ["N", 0], ["N", 0]]

t_total = 24*60*60
t_to_capture = [600, 1800, 3600, 21600, 86400]

####################################################################

boundary_conditions = boundary_conditions_all

t_values = []
t_comp_values = [[],[]]

explicit_data = [[]]
explicit_comp_data = [[],[]]
implicit_data = [[]]

print("--------------")
print("Explicit Scheme")
print("--------------")

t_matrix = deepcopy(t_to_capture)
delta_h = 0.1
delta_t = delta_t_explicit
type = "explicit"
space = grid(boundary_conditions, type, delta_t, delta_h)
space.display()
print("--------------")

steps = math.floor(t_total/delta_t)

for i in range(steps):
    space.update()

    if (((i+1)*100/steps)%5) == 0:
        print("Step |", type, "|", i+1, "/", steps)

    if ((i+1)*delta_t) % 300 == 0:
        t_values.append(((i+1)*delta_t))
        explicit_data[0].append(space.get_T_thermocouple())

    if ((i+1)*delta_t) == t_matrix[0]:
        explicit_data.append(space.get_contour_data())
        t_matrix.pop(0)

    if ((i+1)*delta_t) % 60 == 0:
        if space.get_T_thermocouple() <= 20:
            t_comp_values[0].append(((i+1)*delta_t))
            explicit_comp_data[0].append(space.get_T_thermocouple())

space.display()

#Graphs the solution
for i in range(len(explicit_data[1:])):
    fig, ax = plt.subplots(1,1)
    cp = ax.contourf(explicit_data[i+1][0], explicit_data[i+1][1], explicit_data[i+1][2])
    fig.colorbar(cp)
    ax.set_title("Explicit Method @ " + str(t_to_capture[i]) + " s")
    ax.set_xlabel("r (m)")
    ax.set_ylabel("z (m)")

    plt.show()

plt.plot(t_values, explicit_data[0])
plt.title("Explicit Method at thermocouple location")
plt.xlabel("t (s)")
plt.ylabel("T (C)")
plt.show()

####################################################################
#Evaluates the program from implicit scheme
####################################################################

print("Implicit Scheme")

t_matrix = deepcopy(t_to_capture)
delta_t = delta_t_implicit
delta_h = delta_h_all
type = "implicit"
space = grid(boundary_conditions, type, delta_t, delta_h)

print("--------------")

steps = math.floor(t_total/delta_t)

for i in range(steps):
    space.update()

    if (((i+1)*100/steps)%5) == 0:
        print("Step |", type, "|", i+1, "/", steps)

    if ((i+1)*delta_t) % 300 == 0:
        implicit_data[0].append(space.get_T_thermocouple())

    if ((i+1)*delta_t) == t_matrix[0]:
        implicit_data.append(space.get_contour_data())
        t_matrix.pop(0)

space.display()

#Graphs the solution
for i in range(len(implicit_data[1:])):
    fig, ax = plt.subplots(1,1)
    cp = ax.contourf(implicit_data[i+1][0], implicit_data[i+1][1], implicit_data[i+1][2])
    fig.colorbar(cp)
    ax.set_title("Implicit Method @ " + str(t_to_capture[i]) + " s")
    ax.set_xlabel("r (m)")
    ax.set_ylabel("z (m)")

    plt.show()

plt.plot(t_values, implicit_data[0])
plt.title("Implicit Method at thermocouple location")
plt.xlabel("t (s)")
plt.ylabel("T (C)")
plt.show()

#Graphs the thermocouple temperature variation with time for bot schemes
plt.plot(t_values, explicit_data[0], label = "Explicit Method")
plt.plot(t_values, implicit_data[0], label = "Implicit Method")
plt.title("Explicit vs Implicit Method at thermocouple location")
plt.legend()
plt.xlabel("t (s)")
plt.ylabel("T (C)")
plt.show()

####################################################################
#Evaluates the problem for explicit scheme with different grid and time size from previous explicit scheme
####################################################################

print("Explicit Scheme with different grid size")
print("--------------")

delta_h = delta_h_alternate
delta_t = delta_t_alternate

t_matrix = deepcopy(t_to_capture)
explicit_alt_data = [[]]

type = "explicit"
space = grid(boundary_conditions, type, delta_t, delta_h)
space.display()
print("--------------")

steps = math.floor(t_total/delta_t)

for i in range(steps):
    space.update()
    #print(((i+1)*delta_t))

    if (((i+1)*100/steps)%5) == 0:
        print("Step |", type, "|", i+1, "/", steps)

    if ((i+1)*delta_t) % 60 == 0:
        if space.get_T_thermocouple() <= 20:
            t_comp_values[1].append(((i+1)*delta_t))
            explicit_alt_data[0].append(space.get_T_thermocouple())
            explicit_comp_data[1].append(space.get_T_thermocouple())

    if ((i+1)*delta_t) == t_matrix[0]:
        explicit_alt_data.append(space.get_contour_data())
        t_matrix.pop(0)

space.display()

#Graphs the solution
for i in range(len(explicit_alt_data[1:])):
    fig, ax = plt.subplots(1,1)
    cp = ax.contourf(explicit_alt_data[i+1][0], explicit_alt_data[i+1][1], explicit_alt_data[i+1][2])
    fig.colorbar(cp)
    ax.set_title("Explicit Method with different grid size @ " + str(t_to_capture[i]) + " s")
    ax.set_xlabel("r (m)")
    ax.set_ylabel("z (m)")

    plt.show()

plt.plot(t_comp_values[1], explicit_alt_data[0])
plt.title("Explicit Method with different grid size at thermocouple location")
plt.xlabel("t (s)")
plt.ylabel("T (C)")
plt.show()

#Dispalys the time taken to reach T = 20 at thermocouple location

print("Explicit method with delta_h = ", delta_h_all," and delta_t = ", delta_t_explicit, " takes ", int(t_comp_values[0][-1]), " to reach ", explicit_comp_data[0][-1])
print("Explicit method with delta_h = ", delta_h_alternate," and delta_t = ", delta_t_alternate, " takes ", int(t_comp_values[1][-1]), " to reach ", explicit_comp_data[1][-1])

####################################################################
#Evaluate the problem with convection at outer boundary wall with explicit scheme
####################################################################

print("--------------")
print("Explicit Scheme with convective boundary")
print("--------------")

boundary_conditions = boundary_conditions_convection

delta_h = delta_h_covection
delta_t = delta_t_covection

convection_case_data = [[]]
t_matrix = deepcopy(t_to_capture)
t_values = []

type = "explicit"
space = grid(boundary_conditions, type, delta_t, delta_h)
space.display()
print("--------------")

steps = math.floor(t_total/delta_t)

for i in range(steps):
    space.update()

    if (((i+1)*100/steps)%5) == 0:
        print("Step |", type, "|", i+1, "/", steps)

    if ((i+1)*delta_t) % 300 == 0:
        t_values.append((i+1)*delta_t)
        convection_case_data[0].append(space.get_T_thermocouple())

    if ((i+1)*delta_t) == t_matrix[0]:
        convection_case_data.append(space.get_contour_data())
        t_matrix.pop(0)

space.display()

#Graphs the solution
for i in range(len(convection_case_data[1:])):
    fig, ax = plt.subplots(1,1)
    cp = ax.contourf(convection_case_data[i+1][0], convection_case_data[i+1][1], convection_case_data[i+1][2])
    fig.colorbar(cp)
    ax.set_title("Explicit Method with Convection @ " + str(t_to_capture[i]) + " s")
    ax.set_xlabel("r (m)")
    ax.set_ylabel("z (m)")

    plt.show()

plt.plot(t_values, convection_case_data[0])
plt.title("Explicit Method with Convection at thermocouple location")
plt.xlabel("t (s)")
plt.ylabel("T (C)")
plt.show()