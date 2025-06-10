import numpy as np
import time
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import messagebox, ttk
from datetime import datetime
import random

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



class Food_Mars:
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

        self.create_main_menu()

    def create_main_menu(self):
        self.clear_window()
        tk.Label(self.root, text="화성 생존 시뮬레이터", font=("Arial", 16)).pack(pady=20)
        tk.Button(self.root, text="🍽 식량 관리", width=20, command=self.check_and_open_food).pack(pady=10)
        tk.Button(self.root, text="🌱 농사 관리 (준비 중)", width=20, state=tk.DISABLED).pack(pady=10)

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
        self.clear_window()
        today = datetime.now().date()
        dday_remaining = (self.dday - today).days
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
        info_text = f"총 보유 칼로리: {self.total_calories} kcal   |   예상 생존 가능 일수: {survival_days}일   |   오늘 남은 칼로리: {remaining_cal} kcal   |   하루 기준: {self.daily_calories_needed} kcal   "

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

#blabla