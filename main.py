import numpy as np
import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

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
        self.loss_of_Heat = 18.5  # W/Kdddsfdfn 
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


    mu = np.float64(0.042828e6)  # Mars gravitational parameter (km^3/s^2)

    def __init__(self):
        print("Welcome to the Rendezvous Support System!\nPlease enter data for two spacecraft.")
        c1 = input("Spacecraft 1 position (x y z) [km]: ").replace(',', ' ').split()
        v1 = input("Spacecraft 1 velocity (vx vy vz) [km/s]: ").replace(',', ' ').split()
        c2 = input("Spacecraft 2 position (x y z) [km]: ").replace(',', ' ').split()
        v2 = input("Spacecraft 2 velocity (vx vy vz) [km/s]: ").replace(',', ' ').split()

        self.r1, self.v1 = np.array(c1, dtype=float), np.array(v1, dtype=float)
        self.r2, self.v2 = np.array(c2, dtype=float), np.array(v2, dtype=float)
        print("-" * 40)

        # compute orbital elements and validate ellipticity
        self.elems1 = self._compute_elements(self.r1, self.v1)
        self.elems2 = self._compute_elements(self.r2, self.v2)
        if self.elems1['e'] >= 1 or self.elems2['e'] >= 1:
            raise ValueError("One or both orbits are non-elliptical (e >= 1). Please input e < 1.")

    def _compute_elements(self, r, v):
        tol = 1e-8
        r_norm, v_norm = np.linalg.norm(r), np.linalg.norm(v)
        # semi-major axis
        energy = 0.5 * v_norm**2 - self.mu / r_norm
        a = -self.mu / (2 * energy)

        # eccentricity vector and eccentricity
        h = np.cross(r, v)
        h_norm = np.linalg.norm(h)
        e_vec = np.cross(v, h) / self.mu - r / r_norm
        e = np.linalg.norm(e_vec)

        # inclination
        i = np.arccos(h[2] / h_norm)

        # node vector and RAAN
        k = np.array([0.0, 0.0, 1.0])
        N = np.cross(k, h)
        N_norm = np.linalg.norm(N)
        if N_norm < tol:
            RAAN = 0.0
        else:
            RAAN = np.arccos(np.clip(N[0] / N_norm, -1, 1))
            if N[1] < 0:
                RAAN = 2 * np.pi - RAAN

        # argument of perigee
        if e < tol:
            arg_perigee = 0.0
        else:
            if N_norm < tol:
                # equatorial orbit
                arg_perigee = np.arccos(np.clip(e_vec[0] / e, -1, 1))
                if e_vec[1] < 0:
                    arg_perigee = 2 * np.pi - arg_perigee
            else:
                arg_perigee = np.arccos(np.clip(np.dot(N, e_vec) / (N_norm * e), -1, 1))
                if np.dot(np.cross(N, e_vec), h) < 0:
                    arg_perigee = 2 * np.pi - arg_perigee

        # true anomaly at epoch
        if e < tol:
            # circular orbit
            true_anom0 = np.arccos(np.clip(r[0] / r_norm, -1, 1))
            if r[1] < 0:
                true_anom0 = 2 * np.pi - true_anom0
        else:
            true_anom0 = np.arccos(np.clip(np.dot(e_vec, r) / (e * r_norm), -1, 1))
            if np.dot(r, v) < 0:
                true_anom0 = 2 * np.pi - true_anom0

        # eccentric anomaly and mean anomaly
        E0 = 2 * np.arctan(np.sqrt((1 - e) / (1 + e)) * np.tan(true_anom0 / 2))
        M0 = E0 - e * np.sin(E0)
        n = np.sqrt(self.mu / a**3)  # mean motion

        return {'a': a, 'e': e, 'i': i, 'Ω': RAAN, 'ω': arg_perigee, 'M0': M0, 'n': n}

    def _solve_kepler(self, M, e, tol=1e-8):
        E = M if e < 0.8 else np.pi
        for _ in range(50):
            f = E - e * np.sin(E) - M
            df = 1 - e * np.cos(E)
            dE = -f / df
            E += dE
            if abs(dE) < tol:
                break
        return E

    def _propagate(self, elems, t):
        M = elems['M0'] + elems['n'] * t
        E = self._solve_kepler(M, elems['e'])
        nu = 2 * np.arctan2(np.sqrt(1 + elems['e']) * np.sin(E / 2),
                             np.sqrt(1 - elems['e']) * np.cos(E / 2))

        r_mag = elems['a'] * (1 - elems['e'] * np.cos(E))
        x_orb, y_orb = r_mag * np.cos(nu), r_mag * np.sin(nu)

        p = elems['a'] * (1 - elems['e']**2)
        vx_orb = -np.sqrt(self.mu / p) * np.sin(E)
        vy_orb =  np.sqrt(self.mu / p) * np.sqrt(1 - elems['e']**2) * np.cos(E)

        cos_RAAN, sin_RAAN = np.cos(elems['Ω']), np.sin(elems['Ω'])
        cos_i, sin_i     = np.cos(elems['i']),   np.sin(elems['i'])
        cos_w, sin_w     = np.cos(elems['ω']),   np.sin(elems['ω'])

        R1 = np.array([[ cos_RAAN, -sin_RAAN, 0],
                       [ sin_RAAN,  cos_RAAN, 0],
                       [        0,         0, 1]])
        R2 = np.array([[1,      0,       0],
                       [0,  cos_i, -sin_i],
                       [0,  sin_i,  cos_i]])
        R3 = np.array([[ cos_w, -sin_w, 0],
                       [ sin_w,  cos_w, 0],
                       [     0,      0, 1]])
        R = R1.dot(R2).dot(R3)

        r3d = R.dot([x_orb, y_orb, 0])
        v3d = R.dot([vx_orb, vy_orb, 0])
        return r3d, v3d

    def find_closest_approach(self, t_max=None, steps=5000):
        if t_max is None:
            T1 = 2 * np.pi / self.elems1['n']
            T2 = 2 * np.pi / self.elems2['n']
            t_max = min(T1, T2)

        ts = np.linspace(0, t_max, steps)
        min_d = np.inf
        best = {}
        for t in ts:
            r1, v1 = self._propagate(self.elems1, t)
            r2, v2 = self._propagate(self.elems2, t)
            d = np.linalg.norm(r1 - r2)
            if d < min_d:
                min_d, best = d, {'t': t, 'v1': v1, 'v2': v2}
        rel_speed = np.linalg.norm(best['v1'] - best['v2'])
        return {
            'min_distance_km':    min_d,
            'relative_speed_km_s': rel_speed,
            'time_to_min_s':      best['t']
        }

    def visualize(self, steps=2000):
        # propagate and find best
        T1 = 2 * np.pi / self.elems1['n']
        T2 = 2 * np.pi / self.elems2['n']
        ts = np.linspace(0, min(T1, T2), steps)

        traj1, traj2 = [], []
        min_d, best = np.inf, {}
        for t in ts:
            r1, _ = self._propagate(self.elems1, t)
            r2, _ = self._propagate(self.elems2, t)
            traj1.append(r1); traj2.append(r2)
            d = np.linalg.norm(r1 - r2)
            if d < min_d:
                min_d, best = d, {'t': t, 'r': r1}

        traj1 = np.array(traj1)
        traj2 = np.array(traj2)

        # Mars wireframe
        R_mars = 3390  # km
        u = np.linspace(0, 2*np.pi, 60)
        v = np.linspace(0, np.pi, 30)
        U, V = np.meshgrid(u, v)
        X = R_mars * np.cos(U) * np.sin(V)
        Y = R_mars * np.sin(U) * np.sin(V)
        Z = R_mars * np.cos(V)

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_wireframe(X, Y, Z, rstride=4, cstride=4,
                          alpha=0.2, color='grey')

        ax.plot(traj1[:,0], traj1[:,1], traj1[:,2],
                lw=2, label='Spacecraft 1', color='gold')
        ax.plot(traj2[:,0], traj2[:,1], traj2[:,2],
                lw=2, label='Spacecraft 2', color='orangered')

        rc = best['r']
        ax.scatter(*rc, s=80, marker='X', c='black', label='Closest Approach')
        ax.text(*(rc * 1.05),
                f"d_min={min_d:.1f} km\n t={best['t']:.1f} s",
                fontsize=9, ha='center')

        # equal axis scaling
        all_pts = np.vstack([
            traj1, traj2,
            np.c_[X.flatten(), Y.flatten(), Z.flatten()]
        ])
        center = all_pts.mean(axis=0)
        radius = np.ptp(all_pts, axis=0).max() / 2

        ax.set_xlim(center[0]-radius, center[0]+radius)
        ax.set_ylim(center[1]-radius, center[1]+radius)
        ax.set_zlim(center[2]-radius, center[2]+radius)

        ax.set_xlabel('X (km)')
        ax.set_ylabel('Y (km)')
        ax.set_zlabel('Z (km)')
        ax.set_title('Mars & Spacecraft Orbits with Closest Approach')
        ax.legend()
        plt.tight_layout()
        plt.show()




rv = Rendezvous()
result = rv.find_closest_approach()
print(f"\n최소 접근 거리: {result['min_distance_km']:.3f} km")
print(f"그때 상대 속도: {result['relative_speed_km_s']:.3f} km/s")
print(f"접근까지 걸린 시간: {result['time_to_min_s']:.1f} s")
rv.visualize(steps=3000)


#Heater_RTG()
Rendezvous()
