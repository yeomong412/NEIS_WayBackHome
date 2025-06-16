import numpy as np
import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import tkinter as tk
from tkinter import messagebox, ttk
import random
from datetime import datetime, timedelta
import json


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
        self.elems1 = self.compute_elements(self.r1, self.v1)
        self.elems2 = self.compute_elements(self.r2, self.v2)
        if self.elems1['e'] >= 1 or self.elems2['e'] >= 1:
            raise ValueError("One or both orbits are non-elliptical (e >= 1). Please input e < 1.")

    def compute_elements(self, r, v):
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
        print({'a': a, 'e': e, 'i': i, 'Ω': RAAN, 'ω': arg_perigee, 'M0': M0, 'n': n})

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

         # use a square figure so the subplot itself is square
        fig = plt.figure(figsize=(12,12))
        ax  = fig.add_subplot(111, projection='3d')

        # force equal scaling in x, y, z
        ax.set_box_aspect([1,1,1])
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

class MarsApp:
    def __init__(self, root):
        self.root = root
        self.root.title("화성 생존 시뮬레이터")
        self.root.geometry("740x600")

        self.food_data = {}
        self.consumed = []
        self.produced = []
        self.start_date = None
        self.dday = None
        self.total_calories = 0
        self.daily_calories_needed = 2000
        self.today_calories_consumed = 0
        self.in_farm_screen = False

        self.fire_ready = False
        self.water_ready = False
        self.farming_started = False

        self.load_data()
        self.create_main_menu()

    # ✅ 메뉴 UI
    def create_main_menu(self):
        self.clear_window()
        # create_main_menu() 내 제목 라벨 교체 부분
        title_frame = tk.Frame(self.root, bg="#1a1a1a", padx=20, pady=10)
        title_frame.pack(pady=30)

        tk.Label(
            title_frame,
            text="화성 생존 시뮬레이터",
            font=("Helvetica", 24, "bold"),
            fg="white",
            bg="#1a1a1a"
        ).pack()

        tk.Label(
            title_frame,
            text="Mars Survival Simulator v1.0",
            font=("Helvetica", 12),
            fg="#cccccc",
            bg="#1a1a1a"
        ).pack()

        tk.Button(self.root, text="🍽 식량 관리", width=20, command=self.check_and_open_food).pack(pady=10)
        tk.Button(self.root, text="🌱 농사 관리", width=20, command=self.open_farming_stage_menu).pack(pady=10)

    # ✅ 식량 관리 (원래 코드 그대로 유지)
    def check_and_open_food(self):
        if not self.food_data:
            self.open_food_input()
        else:
            self.open_main_food_screen()

    def open_food_input(self):
        self.clear_window()
        tk.Label(self.root, text="현재 날짜 (YYYY-MM-DD):").pack()
        self.date_entry = tk.Entry(self.root)
        self.date_entry.pack()
        tk.Label(self.root, text="남아있는 식량\n※ (이름:수량:칼로리) 형식으로 입력해주세요.\n예: 감자:20:150, 통조림:5:300").pack()
        self.food_entry = tk.Text(self.root, height=4)
        self.food_entry.pack()
        tk.Label(self.root, text="목표 D-Day (YYYY-MM-DD):").pack()
        self.dday_entry = tk.Entry(self.root)
        self.dday_entry.pack()
        tk.Button(self.root, text="확인", command=self.check_food_inputs).pack(pady=10)
        tk.Button(self.root, text="◀ 돌아가기", command=self.create_main_menu).pack()

    def check_food_inputs(self):
        date = self.date_entry.get().strip()
        food = self.food_entry.get("1.0", tk.END).strip()
        dday = self.dday_entry.get().strip()

        if not date or not food or not dday:
            messagebox.showwarning("입력 오류", "모든 항목을 작성해주세요.")
            return

        try:
            self.start_date = datetime.strptime(date, "%Y-%m-%d").date()
            self.dday = datetime.strptime(dday, "%Y-%m-%d").date()
            self.food_data = {}
            self.total_calories = 0
            for line in food.split(','):
                name, count, cal = line.strip().split(':')
                self.food_data[name.strip()] = [int(count), int(cal)]
                self.total_calories += int(count) * int(cal)
            self.open_main_food_screen()
        except Exception as e:
            messagebox.showerror("형식 오류", f"입력 형식이 잘못되었습니다.\n{e}")

    def open_main_food_screen(self):
        self.in_farm_screen = False
        self.clear_window()
        today = self.start_date
        dday_remaining = (self.dday - self.start_date).days
        tk.Label(self.root, text=f"D-{dday_remaining}", anchor="w").place(x=10, y=10)
        tk.Label(self.root, text=f"{today}", anchor="e").place(x=600, y=10)

        tk.Label(self.root, text="남아있는 식량", font=("Arial", 12)).place(x=20, y=40)
        self.food_tree = ttk.Treeview(self.root, columns=("정보", "총 칼로리"), show="headings", height=10)
        self.food_tree.heading("정보", text="식량 정보")
        self.food_tree.heading("총 칼로리", text="총 칼로리")
        self.food_tree.place(x=20, y=70)

        self.update_food_table()

        right_frame = tk.Frame(self.root, relief=tk.GROOVE, borderwidth=2)
        right_frame.place(x=430, y=50, width=280, height=460)

        tk.Label(right_frame, text="[먹은 식량 기록]", font=("Arial", 10)).pack(pady=2)
        tk.Label(right_frame, text="이름:").pack()
        self.eat_name = tk.Entry(right_frame)
        self.eat_name.pack()
        tk.Label(right_frame, text="수량:").pack()
        self.eat_qty = tk.Entry(right_frame)
        self.eat_qty.pack()
        tk.Button(right_frame, text="섭취 추가", command=self.add_consumed).pack(pady=5)

        tk.Label(right_frame, text="[생산한 식량 기록]", font=("Arial", 10)).pack(pady=2)
        tk.Label(right_frame, text="이름:").pack()
        self.prod_name = tk.Entry(right_frame)
        self.prod_name.pack()
        tk.Label(right_frame, text="수량:").pack()
        self.prod_qty = tk.Entry(right_frame)
        self.prod_qty.pack()
        tk.Label(right_frame, text="칼로리:").pack()
        self.prod_cal = tk.Entry(right_frame)
        self.prod_cal.pack()
        tk.Button(right_frame, text="생산 추가", command=self.add_produced).pack(pady=5)

        tk.Button(self.root, text="💡 생존 팁", command=self.show_tip).place(x=20, y=360)
        tk.Button(self.root, text="🔧 섭취 기준 변경", command=self.change_daily_calories).place(x=120, y=360)
        tk.Button(self.root, text="🏠 홈으로", command=lambda: [self.save_data(), self.create_main_menu()]).place(x=270, y=360)

        self.update_banner()

    def update_food_table(self):
        for row in self.food_tree.get_children():
            self.food_tree.delete(row)
        for name, (qty, cal) in self.food_data.items():
            if qty > 0:
                display = f"{name} ({qty}개, {cal} kcal)"
                self.food_tree.insert("", "end", values=(display, f"{qty * cal} kcal"))

    def update_banner(self):
        survival_days = self.total_calories // self.daily_calories_needed
        remaining_cal = max(0, self.daily_calories_needed - self.today_calories_consumed)
        dday_remain = (self.dday - self.start_date).days
        info_text = f"날짜: {self.start_date} | D-{dday_remain} | 총 보유 칼로리: {self.total_calories} kcal | 예상 생존 가능 일수: {survival_days}일 | 오늘 남은 칼로리: {remaining_cal} kcal | 하루 기준: {self.daily_calories_needed} kcal"

        banner_frame = tk.Frame(self.root, bg="black", height=40)
        banner_frame.place(x=0, y=560, width=740, height=40)
        banner_label = tk.Label(banner_frame, text=info_text, font=("Arial", 12), fg="white", bg="black")
        banner_label.place(x=740, y=5)

        def move():
            x = banner_label.winfo_x()
            if x + banner_label.winfo_width() < 0:
                banner_label.place(x=740, y=5)
            else:
                banner_label.place(x=x - 2, y=5)
            self.root.after(30, move)

        move()

    def add_consumed(self):
        name = self.eat_name.get().strip()
        qty = self.eat_qty.get().strip()
        if name in self.food_data and qty.isdigit():
            qty = int(qty)
            if self.food_data[name][0] >= qty:
                self.food_data[name][0] -= qty
                cal = self.food_data[name][1]
                self.today_calories_consumed += qty * cal
                self.update_food_table()
                self.update_banner()
                messagebox.showinfo("기록됨", f"{name} {qty}개 섭취 기록됨")
            else:
                messagebox.showwarning("재고 부족", "해당 수량만큼 존재하지 않습니다.")
        else:
            messagebox.showerror("입력 오류", "식량이 존재하지 않거나 수량이 잘못되었습니다.")

    def add_produced(self):
        name = self.prod_name.get().strip()
        qty = self.prod_qty.get().strip()
        cal = self.prod_cal.get().strip()
        if name and qty.isdigit() and cal.isdigit():
            qty, cal = int(qty), int(cal)
            if name in self.food_data:
                self.food_data[name][0] += qty
            else:
                self.food_data[name] = [qty, cal]
            self.total_calories += qty * cal
            self.update_food_table()
            self.update_banner()
            messagebox.showinfo("기록됨", f"{name} {qty}개, {cal}kcal 생산 기록됨")

    def change_daily_calories(self):
        def apply():
            new_cal = cal_entry.get()
            if new_cal.isdigit():
                self.daily_calories_needed = int(new_cal)
                top.destroy()
                self.update_banner()
        top = tk.Toplevel(self.root)
        top.title("섭취 기준 변경")
        tk.Label(top, text="하루 필요 칼로리:").pack()
        cal_entry = tk.Entry(top)
        cal_entry.pack()
        tk.Button(top, text="적용", command=apply).pack()

    def show_tip(self):
        tips = [
            "감자는 수확까지 약 70일이 걸립니다.",
            "하루 2000kcal는 성인 기준 생존에 필요한 최소량입니다.",
            "물을 아껴쓰세요. 생존에 필수입니다.",
            "식량은 어두운 곳에 보관하세요.",
            "비상식량은 절대 먼저 먹지 마세요."
        ]
        messagebox.showinfo("생존 팁", random.choice(tips))

    def clear_window(self):
        for widget in self.root.winfo_children():
            widget.destroy()


    def save_data(self):
        data = {
            "food_data": self.food_data,
            "total_calories": self.total_calories,
            "start_date": str(self.start_date),
            "dday": str(self.dday),
            "daily_calories_needed": self.daily_calories_needed,
            "today_calories_consumed": self.today_calories_consumed
        }
        with open("mars_data.json", "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False)

    def load_data(self):
        try:
            with open("mars_data.json", "r", encoding="utf-8") as f:
                data = json.load(f)
                self.food_data = data["food_data"]
                self.total_calories = data["total_calories"]
                self.start_date = datetime.strptime(data["start_date"], "%Y-%m-%d").date()  # 추가
                self.dday = datetime.strptime(data["dday"], "%Y-%m-%d").date()
                self.daily_calories_needed = data["daily_calories_needed"]
                self.today_calories_consumed = data["today_calories_consumed"]
        except Exception:
            pass

    # ✅ 농사 기능 (기존 농사 기능 코드도 그대로 유지)

    def open_farming_stage_menu(self):
        self.clear_window()
        tk.Label(self.root, text="[1단계] 불 피우기").pack(pady=5)
        tk.Button(self.root, text="🔥 불 피우기", command=self.set_fire_ready).pack(pady=5)

        tk.Label(self.root, text="[2단계] 물 만들기").pack(pady=5)
        self.water_btn = tk.Button(self.root, text="💧 물 생성 (불이 켜져야 함)", state=tk.DISABLED, command=self.set_water_ready)
        self.water_btn.pack(pady=5)

        tk.Label(self.root, text="[3단계] 농사 시작").pack(pady=5)
        self.farm_btn = tk.Button(self.root, text="🌱 농사 시작하기", state=tk.DISABLED, command=self.start_or_continue_farming)

        self.farm_btn.pack(pady=5)

        tk.Button(self.root, text="◀ 돌아가기", command=self.create_main_menu).pack(pady=20)
        # ✅ 불과 물이 모두 준비되었고, 농사를 이미 시작했으면 바로 농사 화면 열기
        if self.fire_ready and self.water_ready and self.farming_started:
            self.open_farm_screen()

    def start_or_continue_farming(self):
        if not self.farming_started:
            self.open_farm_input()
        else:
            self.open_farm_screen()


    def open_farm_input(self):
        self.clear_window()
        tk.Label(self.root, text="농사 시작 전 자원을 입력하세요").pack(pady=5)

        self.iridium_entry = tk.Entry(self.root)
        self.fuel_entry = tk.Entry(self.root)
        self.seeds_entry = tk.Entry(self.root)

        tk.Label(self.root, text="이리듐 (일 수):").pack()
        self.iridium_entry.pack()
        tk.Label(self.root, text="연료 (일 수):").pack()
        self.fuel_entry.pack()
        tk.Label(self.root, text="종자 수:").pack()
        self.seeds_entry.pack()

        tk.Label(self.root, text="작물 선택:").pack()
        self.crop_var = tk.StringVar(value="감자")
        tk.OptionMenu(self.root, self.crop_var, "감자", "당근", "배추").pack()

        tk.Button(self.root, text="확인", command=self.start_farming).pack(pady=10)
        tk.Button(self.root, text="◀ 돌아가기", command=self.open_farming_stage_menu).pack()

    def start_farming(self):
        try:
            self.iridium = int(self.iridium_entry.get())
            self.fuel = int(self.fuel_entry.get())
            self.seeds = int(self.seeds_entry.get())
            self.water_days = self.iridium
            self.selected_crop = self.crop_var.get()
            self.crop_info = {
                "감자": {"grow_days": 4},
                "당근": {"grow_days": 5},
                "배추": {"grow_days": 6},
            }
            self.farm_grid = [["o"] * 6 for _ in range(6)]
            self.crop_growth = [[0] * 6 for _ in range(6)]
            self.day_count = 0
            self.open_farm_screen()
        except:
            messagebox.showerror("입력 오류", "숫자를 정확히 입력해주세요.")
        self.farming_started = True  # 플래그 설정

    def open_farm_screen(self):
        self.in_farm_screen = True
        self.clear_window()

        # 상단 상태 표시
        status = f"날짜: {self.start_date} | 작물: {self.selected_crop} | 연료: {self.fuel} | 이리듐: {self.iridium} | 종자: {self.seeds} | 물: {self.water_days}"
        tk.Label(self.root, text=status, font=("Arial", 11)).pack(pady=3)

        # 밭 (6x6)
        frame = tk.Frame(self.root)
        frame.pack()
        self.farm_labels = []
        for i in range(6):
            row = []
            for j in range(6):
                lbl = tk.Label(frame, text=self.farm_grid[i][j], width=2, height=1, borderwidth=1, relief="solid",
                               font=("Courier", 14))
                lbl.grid(row=i, column=j, padx=1, pady=1)
                row.append(lbl)
            self.farm_labels.append(row)

        # 메인 작업 프레임 (좌/우 배치)
        action_frame = tk.Frame(self.root)
        action_frame.pack(pady=10)

        # 왼쪽 작업 (하루, 수확, 다시심기)
        left_frame = tk.Frame(action_frame, relief=tk.GROOVE, borderwidth=2, padx=10, pady=10)
        left_frame.pack(side="left", padx=10)

        tk.Label(left_frame, text="[작물 관리]", font=("Arial", 11)).pack(pady=3)
        tk.Button(left_frame, text="☀ 하루 보내기", command=self.next_day).pack(pady=3)
        tk.Button(left_frame, text="🥕 수확하기", command=self.harvest).pack(pady=3)
        tk.Label(left_frame, text="심을 작물 선택:").pack()
        self.replant_crop_var = tk.StringVar(value="감자")
        tk.OptionMenu(left_frame, self.replant_crop_var, "감자", "당근", "배추").pack()
        tk.Button(left_frame, text="🌱 다시 심기", command=self.replant).pack(pady=3)

        # 오른쪽 작업 (에너지 추가)
        right_frame = tk.Frame(action_frame, relief=tk.GROOVE, borderwidth=2, padx=10, pady=10)
        right_frame.pack(side="right", padx=10)

        tk.Label(right_frame, text="[에너지 추가]", font=("Arial", 11)).pack(pady=3)
        tk.Label(right_frame, text="자원 이름 (연료, 이리듐, 종자 중 택1):").pack()
        self.energy_name = tk.Entry(right_frame)
        self.energy_name.pack()
        tk.Label(right_frame, text="수량:").pack()
        self.energy_amount = tk.Entry(right_frame)
        self.energy_amount.pack()
        tk.Button(right_frame, text="➕ 에너지 추가", command=self.add_energy).pack(pady=3)

        # 하단 홈으로 버튼
        tk.Button(self.root, text="🏠 홈으로", command=self.create_main_menu).pack(pady=10)

    def add_energy(self):
        name = self.energy_name.get().strip()
        amount = self.energy_amount.get().strip()
        if not amount.isdigit():
            messagebox.showwarning("입력 오류", "수량은 숫자로 입력해주세요.")
            return

        amount = int(amount)
        if name == "연료":
            self.fuel += amount
        elif name == "이리듐":
            self.iridium += amount
            self.water_days += amount
        elif name == "종자":
            self.seeds += amount
        else:
            messagebox.showwarning("입력 오류", "연료, 이리듐, 종자 중에서 입력해주세요.")
            return

        self.refresh_farm()
        messagebox.showinfo("추가됨", f"{name} {amount}개 추가됨.")

    def next_day(self):
        self.day_count += 1
        self.fuel -= 1
        self.iridium -= 1
        self.water_days -= 1

        self.start_date += timedelta(days=1)

        if not self.in_farm_screen:  # ← 배너는 식량 화면에서만!
            self.update_banner()

        for i in range(6):
            for j in range(6):
                if self.farm_grid[i][j] in ['o', 'r']:
                    self.crop_growth[i][j] += 1
                    days = self.crop_info[self.selected_crop]["grow_days"]
                    if self.crop_growth[i][j] >= days:
                        self.farm_grid[i][j] = 'Y'
                    else:
                        self.farm_grid[i][j] = 'r'

        if self.water_days <= 0:
            center = 2
            radius = min(3, self.day_count)
            for i in range(6):
                for j in range(6):
                    if abs(i - center) + abs(j - center) <= radius:
                        if self.farm_grid[i][j] in ['o', 'r']:
                            self.farm_grid[i][j] = 'x'

        self.show_warnings()
        self.refresh_farm()

    def harvest(self):
        crop_kcal_map = {"감자": 150, "당근": 120, "배추": 100}
        count = 0
        for i in range(6):
            for j in range(6):
                if self.farm_grid[i][j] == 'Y':
                    self.farm_grid[i][j] = '.'
                    self.crop_growth[i][j] = 0
                    count += 1
                    self.seeds += 1

        # 수확한 작물을 식량 데이터에 반영
        if count > 0:
            cal = crop_kcal_map[self.selected_crop]
            self.total_calories += count * cal
            if self.selected_crop in self.food_data:
                self.food_data[self.selected_crop][0] += count
            else:
                self.food_data[self.selected_crop] = [count, cal]

        messagebox.showinfo("수확 완료", f"{count}개 수확 완료!")
        self.refresh_farm()

    def replant(self):
        crop_symbol = 'o'  # 기본 표시 기호 (여기선 구분용이 아니라 시각용이라 고정)
        selected_crop = self.replant_crop_var.get()
        self.selected_crop = selected_crop  # 선택된 작물 갱신

        for i in range(6):
            for j in range(6):
                if self.farm_grid[i][j] == '.' and self.seeds > 0:
                    self.farm_grid[i][j] = crop_symbol
                    self.crop_growth[i][j] = 0
                    self.seeds -= 1

        self.refresh_farm()

    def refresh_farm(self):
        for i in range(6):
            for j in range(6):
                self.farm_labels[i][j].config(text=self.farm_grid[i][j])

    def show_warnings(self):
        for name, value in [("연료", self.fuel), ("이리듐", self.iridium), ("물", self.water_days), ("종자", self.seeds)]:
            if value in [5, 4, 3, 2, 1]:
                messagebox.showwarning("자원 경고", f"⚠️ {name}이 {value}일치밖에 남지 않았습니다.")

    def set_fire_ready(self):
        self.fire_ready = True
        self.water_btn.config(state=tk.NORMAL)

    def set_water_ready(self):
        if self.fire_ready:
            self.water_ready = True
            self.farm_btn.config(state=tk.NORMAL)
    # ✅ 공통 함수들
    def clear_window(self):
        for widget in self.root.winfo_children():
            widget.destroy()

    def save_data(self):
        data = {
            "food_data": self.food_data,
            "total_calories": self.total_calories,
            "start_date": str(self.start_date),
            "dday": str(self.dday),
            "daily_calories_needed": self.daily_calories_needed,
            "today_calories_consumed": self.today_calories_consumed
        }
        with open("mars_data.json", "w", encoding="utf-8") as f:
            json.dump(data, f, ensure_ascii=False)

    def load_data(self):
        try:
            with open("mars_data.json", "r", encoding="utf-8") as f:
                data = json.load(f)
                self.food_data = data["food_data"]
                self.total_calories = data["total_calories"]
                self.start_date = datetime.strptime(data["start_date"], "%Y-%m-%d").date()
                self.dday = datetime.strptime(data["dday"], "%Y-%m-%d").date()
                self.daily_calories_needed = data["daily_calories_needed"]
                self.today_calories_consumed = data["today_calories_consumed"]
        except Exception:
            pass


#Heater_RTG()


rv = Rendezvous()
result = rv.find_closest_approach()
print(f"\n최소 접근 거리: {result['min_distance_km']:.3f} km")
print(f"상대 속도: {result['relative_speed_km_s']:.3f} km/s")
print(f"접근까지 걸린 시간: {result['time_to_min_s']:.1f} s")
rv.visualize(steps=3000)



#Rendezvous()


# # ✅ 실행부
# if __name__ == "__main__":
#     root = tk.Tk()
#     app = MarsApp(root)
#     root.mainloop()
