import numpy as np
import time
import matplotlib.pyplot as plt

class task():
    @classmethod
    def CtoK(cls, c): return 273.15 + c
    
class Heater_RTG(task):
    def __init__(self):

        self.heat_emission = 1500  # J/s
        self.init_temp = task.CtoK(-20)  # K
        self.temp_now = self.init_temp
        self.temp_past = -999
        self.density_air = 1.29  # kg/m³
        self.specificHeat_air = 1005  # J/kg/K
        self.outside_temp = task.CtoK(-58)  # K
        self.loss_of_Heat = 18.5  # W/K
        self.time_s = 1

        self.time_data = []
        self.temp_data = []

        plt.ion()
        self.fig, self.ax = plt.subplots()
        self.line, = self.ax.plot([], [], 'r.')
        self.ax.set_xlabel("Time (s)")
        self.ax.set_ylabel("Temperature (°C)")
        self.ax.set_title("RTG Heating Simulation")
        self.ax.grid(True)
        
        self.run_simulation()

    def run_simulation(self):
        while abs(self.temp_now - self.temp_past) > 5e-5:
            self.temp_past = self.temp_now
            exponent = -self.loss_of_Heat * self.time_s / (self.density_air * self.specificHeat_air)
            self.temp_now = self.outside_temp + (self.init_temp - self.outside_temp - (self.heat_emission / self.loss_of_Heat)) * \
                            np.exp(exponent) + (self.heat_emission / self.loss_of_Heat)

            temp_celsius = self.temp_now - 273.15
            print(f"{self.time_s}s: {temp_celsius:.4f} °C")

            self.time_data.append(self.time_s)
            self.temp_data.append(temp_celsius)

            self.line.set_data(self.time_data, self.temp_data)
            self.ax.relim()
            self.ax.autoscale_view()
            plt.draw()
            plt.pause(0.01)

            self.time_s += 1
            time.sleep(0.001)

        plt.ioff()
        plt.show()

class Rendezvous(task):
    def __init__(self):
        pass #

import time
import random
class sandstorm():
    def __init__(self,mode):
        if(mode=="easy"):
            self.stormcome(10)
        if(mode=="intermediate"):
            self.stormcome(30)
        if(mode=="hard"):
            self.stormcome(50)
        
    def stormcome(self, x):
        l=[[0 for i in range(19)] for j in range(19)]
        a,b=map(int,input("Choose your Base Coordinates").split())
        l[a][b]='H'
        for k in range(100):
            num1 = random.randint(0, 18)
            num2 = random.randint(0, 18)
            num3 = random.randint(0,1)
            if(num3):
                for i in range(0,19):
                    l[num1][i]=1
            else:
                for i in range(0,19):
                    l[i][num2]=1
            move=input()        
            time.sleep(0.7)
            for i in range(19):
                for j in range(19):
                    print(l[i][j],end=" ")
                print()
            for i in range(19):
                for j in range(19):
                    if(l[i][j]==1):
                        l[i][j]='S'
                        if(l[a][b]=='S'):
                            print("Oh no The Base Has Been Destroyed")
                            break
                    elif(l[i][j]=='S'):
                        l[i][j]=0

Heater_RTG()
