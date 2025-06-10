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
        pass

Heater_RTG()
