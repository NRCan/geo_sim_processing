# coding: utf-8

# In[21]:

#%matplotlib ipympl
import numpy as np
import random as rd
#rd.seed(5)

import matplotlib.pyplot as plt

temps = 0
L = 50
n_pred0 = 10
n_proie0 = 200
reserves0 = 6.0

# initialisation de la grille
V = np.empty([L, L], dtype=object)
W = np.empty([L, L], dtype=object)


class proie():
    vie = 30
    age_division = 5

    def __init__(self, x, y, age):
        self.x = x
        self.y = y
        self.age = age
        V[x, y] = self

    def __repr__(self):
        #return 'P'+str(self.age)
        return 'P'



    def move(self):
        if self.age >= self.age_division:
            for i in range(5):
                a = rd.randint(1, 9)
                if a == 1:
                    if V[(self.x - 1) % L, (self.y + 1) % L] is None and W[
                        (self.x - 1) % L, (self.y + 1) % L] is None:
                        V[self.x, self.y].age = 0
                        proie((self.x - 1) % L, (self.y + 1) % L, 0)
                        break
                elif a == 2:
                    if V[(self.x) % L, (self.y + 1) % L] is None and W[(self.x) % L, (self.y + 1) % L] is None:
                        V[self.x, self.y].age = 0
                        proie((self.x) % L, (self.y + 1) % L, 0)
                        break
                elif a == 3:
                    if V[(self.x + 1) % L, (self.y + 1) % L] is None and W[
                        (self.x + 1) % L, (self.y + 1) % L] is None:
                        V[self.x, self.y].age = 0
                        proie((self.x + 1) % L, (self.y + 1) % L, 0)
                        break
                elif a == 4:
                    if V[(self.x + 1) % L, (self.y) % L] is None and W[(self.x + 1) % L, (self.y) % L] is None:
                        V[self.x, self.y].age = 0
                        proie((self.x + 1) % L, (self.y) % L, 0)
                        break
                elif a == 5:
                    if V[(self.x + 1) % L, (self.y - 1) % L] is None and W[
                        (self.x + 1) % L, (self.y - 1) % L] is None:
                        V[self.x, self.y].age = 0
                        proie((self.x + 1) % L, (self.y - 1) % L, 0)
                        break
                elif a == 6:
                    if V[(self.x) % L, (self.y - 1) % L] is None and W[(self.x) % L, (self.y - 1) % L] is None:
                        V[self.x, self.y].age = 0
                        proie((self.x) % L, (self.y - 1) % L, 0)
                        break
                elif a == 7:
                    if V[(self.x - 1) % L, (self.y - 1) % L] is None and W[
                        (self.x - 1) % L, (self.y - 1) % L] is None:
                        V[self.x, self.y].age = 0
                        proie((self.x - 1) % L, (self.y - 1) % L, 0)
                        break
                elif a == 8:
                    if V[(self.x - 1) % L, (self.y) % L] is None and W[(self.x - 1) % L, (self.y) % L] is None:
                        V[self.x, self.y].age = 0
                        proie((self.x - 1) % L, (self.y) % L, 0)
                        break

    def temps(self):
        self.age += 1

    def mort(self):
        if self.age >= self.vie:
            V[self.x, self.y] = None
            del self


class pred():
    vie = 100
    age_division = 20
    reserves_division = 10
    reserves_satiete = 20

    def __init__(self, x, y, age, reserves):
        self.x = x
        self.y = y
        self.age = age
        self.reserves = reserves
        self.dep = 0
        W[x, y] = self

    def __repr__(self):
        #return 'V '+str(self.reserves)+" "+" "+str(self.age)
        return 'V'


    def manger(self):
        if self.reserves < self.reserves_satiete:
            if V[self.x, self.y] is not None:
                self.reserves += 1
                V[self.x, self.y] = None

    def deplacement(self):
        for i in range(5):
            if self.dep == 0:
                a = rd.randint(1, 9)
                if a == 1:
                    if V[(self.x - 1) % L, (self.y + 1) % L] is not None and W[ (self.x - 1) % L, (self.y + 1) % L] is None:
                        W[self.x, self.y] = None
                        #---
                        self.x = (self.x - 1) % L
                        self.y = (self.y + 1) % L
                        W[self.x,self.y] = self
                        #---
                        #W[(self.x - 1) % L, (self.y + 1) % L] = self
                        self.dep = 1
                        break
                elif a == 2:
                    if V[(self.x) % L, (self.y + 1) % L] is not None and W[(self.x) % L, (self.y + 1) % L] is None:
                        W[self.x, self.y] = None
                        #W[(self.x) % L, (self.y + 1) % L] = self
                        self.x = (self.x ) % L
                        self.y = (self.y + 1) % L
                        W[self.x, self.y] = self
                        self.dep = 1
                        break
                elif a == 3:
                    if V[(self.x + 1) % L, (self.y + 1) % L] is not None and W[
                        (self.x + 1) % L, (self.y + 1) % L] is None:
                        W[self.x, self.y] = None
                        #W[(self.x + 1) % L, (self.y + 1) % L] = self
                        self.x = (self.x +1) % L
                        self.y = (self.y + 1) % L
                        W[self.x, self.y] = self
                        self.dep = 1
                        break
                elif a == 4:
                    if V[(self.x + 1) % L, (self.y) % L] is not None and W[(self.x + 1) % L, (self.y) % L] is None:
                        W[self.x, self.y] = None
                        #W[(self.x + 1) % L, (self.y) % L] = self
                        self.x = (self.x + 1) % L
                        self.y = (self.y ) % L
                        W[self.x, self.y] = self
                        self.dep = 1
                        break
                elif a == 5:
                    if V[(self.x + 1) % L, (self.y - 1) % L] is not None and W[
                        (self.x + 1) % L, (self.y - 1) % L] is None:
                        W[self.x, self.y] = None
                        #W[(self.x + 1) % L, (self.y - 1) % L] = self
                        self.x = (self.x + 1) % L
                        self.y = (self.y - 1) % L
                        W[self.x,self.y] = self
                        self.dep = 1
                        break
                elif a == 6:
                    if V[(self.x) % L, (self.y - 1) % L] is not None and W[(self.x) % L, (self.y - 1) % L] is None:
                        W[self.x, self.y] = None
                        #W[(self.x) % L, (self.y - 1) % L] = self
                        self.x = (self.x ) % L
                        self.y = (self.y - 1) % L
                        W[self.x,self.y] = self
                        self.dep = 1
                        break
                elif a == 7:
                    if V[(self.x - 1) % L, (self.y - 1) % L] is not None and W[
                        (self.x - 1) % L, (self.y - 1) % L] is None:
                        W[self.x, self.y] = None
                        #W[(self.x - 1) % L, (self.y - 1) % L] = self
                        self.x = (self.x - 1) % L
                        self.y = (self.y - 1) % L
                        W[self.x,self.y] = self
                        self.dep = 1
                        break
                elif a == 8:
                    if V[(self.x - 1) % L, (self.y) % L] is not None and W[(self.x - 1) % L, (self.y) % L] is None:
                        W[self.x, self.y] = None
                        #W[(self.x - 1) % L, (self.y) % L] = self
                        self.x = (self.x - 1) % L
                        self.y = (self.y) % L
                        W[self.x,self.y] = self
                        self.dep = 1
                        break

    def copulation(self):
        if self.age >= self.age_division and self.reserves >= self.reserves_division:
            for i in range(5):
                a = rd.randint(1, 9)
                if a == 1:
                    if W[(self.x - 1) % L, (self.y + 1) % L] is None:
                        W[self.x, self.y].age = 0
                        W[self.x, self.y].reserves = self.reserves / 2
                        pred((self.x - 1) % L, (self.y + 1) % L, 0, self.reserves / 2)
                        break
                elif a == 2:
                    if W[(self.x) % L, (self.y + 1) % L] is None:
                        W[self.x, self.y].age = 0
                        W[self.x, self.y].reserves = self.reserves / 2
                        pred((self.x) % L, (self.y + 1) % L, 0, self.reserves / 2)
                        break
                elif a == 3:
                    if W[(self.x + 1) % L, (self.y + 1) % L] is None:
                        W[self.x, self.y].age = 0
                        W[self.x, self.y].reserves = self.reserves / 2
                        pred((self.x + 1) % L, (self.y + 1) % L, 0, self.reserves / 2)
                        break
                elif a == 4:
                    if W[(self.x + 1) % L, (self.y) % L] is None:
                        W[self.x, self.y].age = 0
                        W[self.x, self.y].reserves = self.reserves / 2
                        pred((self.x + 1) % L, (self.y) % L, 0, self.reserves / 2)
                        break
                elif a == 5:
                    if W[(self.x + 1) % L, (self.y - 1) % L] is None:
                        W[self.x, self.y].age = 0
                        W[self.x, self.y].reserves = self.reserves / 2
                        pred((self.x + 1) % L, (self.y - 1) % L, 0, self.reserves / 2)
                        break
                elif a == 6:
                    if W[(self.x) % L, (self.y - 1) % L] is None:
                        W[self.x, self.y].age = 0
                        W[self.x, self.y].reserves = self.reserves / 2
                        pred((self.x) % L, (self.y - 1) % L, 0, self.reserves / 2)
                        break
                elif a == 7:
                    if W[(self.x - 1) % L, (self.y - 1) % L] is None:
                        W[self.x, self.y].age = 0
                        W[self.x, self.y].reserves = self.reserves / 2
                        pred((self.x - 1) % L, (self.y - 1) % L, 0, self.reserves / 2)
                        break
                elif a == 8:
                    if W[(self.x - 1) % L, (self.y) % L] is None:
                        W[self.x, self.y].age = 0
                        W[self.x, self.y].reserves = self.reserves / 2
                        pred((self.x - 1) % L, (self.y) % L, 0, self.reserves / 2)
                        break

    def effet(self):
        self.age += 1
        self.reserves += -0.5
        self.dep = 0

    def dead(self):
        if self.age >= self.vie or self.reserves <= 0:
            W[self.x, self.y] = None
            del self


# ---------------------------------------------------------------------------
# gille de visualisation
def print_grille():
    for i in range(L):
        ligne = ''
        for j in range(L):
            if V[i, j] is None and W[i, j] is None:
                ligne = ligne + '_'
            elif V[i, j].__repr__() == 'P' and W[i, j].__repr__() == 'V':
                ligne = ligne + 'X'
            elif V[i, j] is not None:
                ligne = ligne + V[i, j].__repr__()
            elif W[i, j] is not None:
                ligne = ligne + W[i, j].__repr__()
            else:
                ligne += '_'
        print(ligne)


# ---------------------------------------------------------------------

# Initialisation du tableau de proies
for i in range(n_proie0):
    x = rd.randint(0, L - 1)
    y = rd.randint(0, L - 1)
    age = rd.randint(0, 31)
    proie(x, y, age)

# Initialisation du tableau de prédateur
for i in range(n_pred0):
    x = rd.randint(0, L - 1)
    y = rd.randint(0, L - 1)
    age = rd.randint(0, 101)
    reserves = reserves0
    pred(x, y, age, reserves)

# Boucle de temps


A = []
B = []

print_grille()
print ("---------------")
for t in range(1000):
    a = 0
    b = 0
    for i in range(L):
        for j in range(L):
            if W[i, j] is not None:
                W[i, j].dead()
            if V[i, j] is not None:
                V[i, j].mort()
    for i in range(L):
        for j in range(L):
            if W[i, j] is not None:
                W[i, j].copulation()
            if V[i, j] is not None:
                V[i, j].move()
    for i in range(L):
        for j in range(L):
            if W[i, j] is not None:
                W[i, j].manger()
    for i in range(L):
        for j in range(L):
            if W[i, j] is not None:
                W[i, j].deplacement()
    for i in range(L):
        for j in range(L):
            if W[i, j] is not None:
                W[i, j].effet()
            if V[i, j] is not None:
                V[i, j].temps()
    for i in range(L):
        for j in range(L):
            if V[i, j] is not None:
                a += 1
            if W[i, j] is not None:
                b += 1
    A.append(a)
    B.append(b)
    print (t)
    #print_grille()
    #print('-------------------------------------------------------------------------')


# for i in range (L):
#    for j in range (L):
#        if W[i,j] is not None:
#            W[i,j].manger()
#

# for i in range (L):
#    for j in range (L):
#        if V[i,j] is not None:
#            V[i,j].move()




#import numpy as np
#import random as rd
#import matplotlib.pyplot as plt

print(A)
print(B)
X = range(0, 1000)

plt.plot(X, A, X, B)

print (10)



