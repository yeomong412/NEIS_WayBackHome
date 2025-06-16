# import numpy as np
# import time
# import matplotlib.pyplot as plt

# class task():
#     @classmethod
#     def CtoK(cls, c): return 273.15 + c
    
# class Heater_RTG(task):
#     def __init__(self):

#         self.heat_emission = 1500  # J/s
#         self.init_temp = task.CtoK(-20)  # K
#         self.temp_now = self.init_temp
#         self.temp_past = -999
#         self.density_air = 1.29  # kg/m³
#         self.specificHeat_air = 1005  # J/kg/K
#         self.outside_temp = task.CtoK(-58)  # K
#         self.loss_of_Heat = 18.5  # W/K
#         self.time_s = 1

#         self.time_data = []
#         self.temp_data = []

#         plt.ion()
#         self.fig, self.ax = plt.subplots()
#         self.line, = self.ax.plot([], [], 'r.')
#         self.ax.set_xlabel("Time (s)")
#         self.ax.set_ylabel("Temperature (°C)")
#         self.ax.set_title("RTG Heating Simulation")
#         self.ax.grid(True)
        
#         self.run_simulation()

#     def run_simulation(self):
#         while abs(self.temp_now - self.temp_past) > 5e-5:
#             self.temp_past = self.temp_now
#             exponent = -self.loss_of_Heat * self.time_s / (self.density_air * self.specificHeat_air)
#             self.temp_now = self.outside_temp + (self.init_temp - self.outside_temp - (self.heat_emission / self.loss_of_Heat)) * \
#                             np.exp(exponent) + (self.heat_emission / self.loss_of_Heat)

#             temp_celsius = self.temp_now - 273.15
#             print(f"{self.time_s}s: {temp_celsius:.4f} °C")

#             self.time_data.append(self.time_s)
#             self.temp_data.append(temp_celsius)

#             self.line.set_data(self.time_data, self.temp_data)
#             self.ax.relim()
#             self.ax.autoscale_view()
#             plt.draw()
#             plt.pause(0.01)

#             self.time_s += 1
#             time.sleep(0.001)

#         plt.ioff()
#         plt.show()

# class Rendezvous(task):
#     def __init__(self):
#         pass #

# import time
# import random
# class sandstorm():
#     def __init__(self,mode):
#         if(mode=="easy"):
#             self.stormcome(10)
#         if(mode=="intermediate"):
#             self.stormcome(30)
#         if(mode=="hard"):
#             self.stormcome(50)
        
#     def stormcome(self, x):
#         l=[[0 for i in range(19)] for j in range(19)]
#         a,b=map(int,input("Choose your Base Coordinates").split())
#         l[a-1][b-1]='H'
#         for k in range(100):
#             num1 = random.randint(0, 18)
#             num2 = random.randint(0, 18)
#             num3 = random.randint(0,1)
#             if(num3):
#                 for i in range(0,19):
#                     l[num1][i]=1
#             else:
#                 for i in range(0,19):
#                     l[i][num2]=1
#             if(k>=1):
#                 p,q=map(int,input().split())
#                 if(l[p-1][q-1]==1):
#                     return "Oh no the Base Has Been Destroyed"
#                 l[p-1][q-1]="H"      
#             time.sleep(0.7)
#             for i in range(19):
#                 for j in range(19):#
#                     print(l[i][j],end=" ") 
#                 print()
#             for i in range(19):
#                 for j in range(19):
#                     if(l[i][j]==1):
#                         l[i][j]='S'
#                         if(l[a][b]=='S'):
#                             return "Oh no The Base Has Been Destroyed"
#                     elif(l[i][j]=='S'):
#                         l[i][j]=0


# sandstorm("easy")

# import time
# import random
# import os

# class SandstormGame:
#     def __init__(self, mode="easy"):
#         # 난이도별 스톰 횟수 설정
#         levels = {"easy": 10, "intermediate": 30, "hard": 50}
#         self.turns = levels.get(mode, 10)
#         self.size = 19
#         self.map = [[0]*self.size for _ in range(self.size)]
#         self.place_base()
#         self.play()

#     def place_base(self):
#         while True:
#             try:
#                 a, b = map(int, input(f"기지 좌표를 선택하세요 (1~{self.size}): ").split())
#                 assert 1 <= a <= self.size and 1 <= b <= self.size
#                 self.base = (a-1, b-1)
#                 self.map[a-1][b-1] = 'H'
#                 break
#             except (ValueError, AssertionError):
#                 print("잘못된 입력입니다. 다시 시도해주세요.")

#     def print_map(self):
#         os.system('cls' if os.name=='nt' else 'clear')
#         for row in self.map:
#             print(' '.join(str(cell) for cell in row))
#         print()

#     def play(self):
#         for turn in range(1, self.turns+1):
#             # 새로운 폭풍 행 또는 열 랜덤 선택
#             r = random.randint(0, self.size-1)
#             c = random.randint(0, self.size-1)
#             if random.choice([True, False]):
#                 # 행 폭풍
#                 for j in range(self.size):
#                     self.map[r][j] = 'S'
#             else:
#                 # 열 폭풍
#                 for i in range(self.size):
#                     self.map[i][c] = 'S'

#             # 기지가 폭풍에 휩쓸렸는지 체크
#             if self.map[self.base[0]][self.base[1]] == 'S':
#                 self.print_map()
#                 print(f"턴 {turn}: Oh no, 기지가 파괴되었습니다! Game Over.")
#                 return

#             # 맵 출력
#             self.print_map()
#             print(f"턴 {turn}: 안전합니다. 기지 좌표 {self.base[0]+1},{self.base[1]+1} 보호 중.")

#             # 다음 턴을 위해 기지를 다시 찍어두고, S는 지워줌
#             self.map[self.base[0]][self.base[1]] = 'H'
#             for i in range(self.size):
#                 for j in range(self.size):
#                     if self.map[i][j] == 'S':
#                         self.map[i][j] = 0

#             time.sleep(0.7)

#         self.print_map()
#         print("모든 턴을 생존했습니다! Victory!")

# if __name__ == "__main__":
#     mode = input("난이도(easy, intermediate, hard) 중 선택: ").strip().lower()
#     game = SandstormGame(mode)









# import time, random, os

# class SandstormGame:
#     def __init__(self, mode="easy"):
#         specs = {
#             "easy":        {"turns":10, "storms":(1,1)},
#             "intermediate":{"turns":30, "storms":(5,5)},
#             "hard":        {"turns":50, "storms":(8,8)},
#         }
#         self.spec        = specs.get(mode, specs["easy"])
#         self.size        = 19
#         self.map         = [[0]*self.size for _ in range(self.size)]
#         self.base        = None
#         self.next_storms = []

#         self.place_base()
#         # 첫 턴 전 예고 뿌리기
#         self.plan_next_storm()
#         self.play()

#     def place_base(self):
#         while True:
#             try:
#                 a,b = map(int, input(f"초기 기지 좌표 (1~{self.size}): ").split())
#                 assert 1<=a<=self.size and 1<=b<=self.size
#                 self.base = (a-1, b-1)
#                 self.map[a-1][b-1] = 'H'
#                 break
#             except:
#                 print("잘못된 입력입니다. 다시 시도해주세요.")

#     def print_map(self, title=""):
#         os.system('cls' if os.name=='nt' else 'clear')
#         if title: print("──", title, "──")
#         for row in self.map:
#             print(' '.join(str(cell) for cell in row))
#         print()

#     def apply_storm(self):
#         # 예고된 1 → 실제폭풍 S 로 전환
#         for typ, idx in self.next_storms:
#             if typ=='row':
#                 for j in range(self.size):
#                     self.map[idx][j] = 'S'
#             else:
#                 for i in range(self.size):
#                     self.map[i][idx] = 'S'
#         self.next_storms = []

#     def clear_storm(self):
#         # S를 모두 0으로, H(기지)만 남김
#         bx, by = self.base
#         for i in range(self.size):
#             for j in range(self.size):
#                 if self.map[i][j] == 'S':
#                     self.map[i][j] = 0
#         self.map[bx][by] = 'H'

#     def plan_next_storm(self):
#         # 다음 턴 예고: 랜덤 줄/열 골라서 1로 표시
#         storm_min, storm_max = self.spec["storms"]
#         cnt = random.randint(storm_min, storm_max)
#         plans = []
#         for _ in range(cnt):
#             if random.choice([True, False]):
#                 plans.append(('row', random.randrange(self.size)))
#             else:
#                 plans.append(('col', random.randrange(self.size)))

#         for typ, idx in plans:
#             if typ=='row':
#                 for j in range(self.size):
#                     if self.map[idx][j] != 'H':
#                         self.map[idx][j] = 1
#             else:
#                 for i in range(self.size):
#                     if self.map[i][idx] != 'H':
#                         self.map[i][idx] = 1

#         self.next_storms = plans

#     def play(self):
#         for turn in range(1, self.spec["turns"]+1):
#             # 1) 예고된 '1'을 보여준다
#             self.print_map(f"TURN {turn} – 예고(1): 다음 턴 폭풍 위치")
#             bx, by = self.base

#             # 2) 기지 이동
#             while True:
#                 c = input("기지 이동? (y/n): ").strip().lower()
#                 if c in ('y','n'): break
#                 print("y 또는 n만 입력해주세요.")
#             if c == 'y':
#                 while True:
#                     try:
#                         x, y = map(int, input(f"새 기지 좌표 (1~{self.size}): ").split())
#                         assert 1<=x<=self.size and 1<=y<=self.size
#                         ox, oy = self.base
#                         self.map[ox][oy] = 0
#                         self.base = (x-1, y-1)
#                         self.map[x-1][y-1] = 'H'
#                         break
#                     except:
#                         print("잘못된 입력입니다. 다시 시도해주세요.")

#             # 3) 실제폭풍: 1 → S
#             self.apply_storm()

#             # 4) S 상태 보여주기
#             self.print_map(f"TURN {turn} – 실제폭풍(S): 지금 폭풍이 몰아칩니다!")
#             # 5) 파괴 체크
#             if self.map[self.base[0]][self.base[1]] == 'S':
#                 print("Oh no, 기지가 파괴되었습니다… Game Over")
#                 return

#             # 6) S 제거, 기지 복원
#             self.clear_storm()

#             # 7) 다음 턴 예고 준비 (map에 1 표시, 그러나 화면엔 다음 루프에서만 보임)
#             self.plan_next_storm()

#             time.sleep(0.7)

#         # 모두 생존
#         self.print_map("모든 턴 생존! Victory!")
#         print()

# if __name__ == "__main__":
#     mode = input("난이도(easy, intermediate, hard): ").strip().lower()
#     SandstormGame(mode)









#최종본
# import time, random, os

# class SandstormGame:
#     def __init__(self, mode="easy"):
#         specs = {
#             "easy":        {"turns":10, "storms":(1,1)},
#             "intermediate":{"turns":30, "storms":(5,5)},
#             "hard":        {"turns":50, "storms":(8,8)},
#         }
#         self.spec        = specs.get(mode, specs["easy"])
#         self.size        = 19
#         self.map         = [[0]*self.size for _ in range(self.size)]
#         self.base        = None
#         self.next_storms = []

#         self.place_base()
#         # 첫 턴 전 예고 뿌리기
#         self.plan_next_storm()
#         self.play()

#     def place_base(self):
#         while True:
#             try:
#                 a,b = map(int, input(f"초기 기지 좌표 (1~{self.size}): ").split())
#                 assert 1<=a<=self.size and 1<=b<=self.size
#                 self.base = (a-1, b-1)
#                 self.map[a-1][b-1] = 'H'
#                 break
#             except:
#                 print("잘못된 입력입니다. 다시 시도해주세요.")

#     def print_map(self, title=""):
#         os.system('cls' if os.name=='nt' else 'clear')
#         if title:
#             print(f"── {title} ──")
#         for row in self.map:
#             print(' '.join(str(cell) for cell in row))
#         print()

#     def apply_storm(self):
#         # 예고된 1 → 실제폭풍 S
#         for typ, idx in self.next_storms:
#             if typ == 'row':
#                 for j in range(self.size):
#                     self.map[idx][j] = 'S'
#             else:
#                 for i in range(self.size):
#                     self.map[i][idx] = 'S'
#         self.next_storms = []

#     def clear_storm(self):
#         # S 전부 제거, H만 남기기
#         bx, by = self.base
#         for i in range(self.size):
#             for j in range(self.size):
#                 if self.map[i][j] == 'S':
#                     self.map[i][j] = 0
#         self.map[bx][by] = 'H'

#     def plan_next_storm(self):
#         # 다음 턴 예고: 1로 표시
#         storm_min, storm_max = self.spec["storms"]
#         cnt = random.randint(storm_min, storm_max)
#         plans = []
#         for _ in range(cnt):
#             if random.choice([True, False]):
#                 plans.append(('row', random.randrange(self.size)))
#             else:
#                 plans.append(('col', random.randrange(self.size)))
#         for typ, idx in plans:
#             if typ == 'row':
#                 for j in range(self.size):
#                     if self.map[idx][j] != 'H':
#                         self.map[idx][j] = 1
#             else:
#                 for i in range(self.size):
#                     if self.map[i][idx] != 'H':
#                         self.map[i][idx] = 1
#         self.next_storms = plans

#     def play(self):
#         moves = {'W':(-1,0), 'A':(0,-1), 'S':(1,0), 'D':(0,1)}
#         for turn in range(1, self.spec["turns"]+1):
#             # 1) 예고(1) 표시 상태 보여주기
#             self.print_map(f"TURN {turn} – 예고(1): 다음 턴 폭풍 위치")

#             # 2) 방향 입력
#             bx, by = self.base
#             while True:
#                 cmd = input("이동 방향 (W/A/S/D) 또는 이동 안 함(N): ").strip().upper()
#                 if cmd in ('W','A','S','D','N'):
#                     break
#                 print("W, A, S, D, N 중 하나를 입력하세요.")
#             if cmd != 'N':
#                 dx, dy = moves[cmd]
#                 nx, ny = bx+dx, by+dy
#                 if 0 <= nx < self.size and 0 <= ny < self.size:
#                     # 이동 가능
#                     self.map[bx][by] = 0
#                     self.base = (nx, ny)
#                     self.map[nx][ny] = 'H'
#                 else:
#                     print("그 방향으로는 이동할 수 없습니다.")

#             # 3) 실제폭풍: 1 → S
#             self.apply_storm()

#             # 4) 실제폭풍 상태 보여주기
#             self.print_map(f"TURN {turn} – 실제폭풍(S): 지금 폭풍이 몰아칩니다!")
#             # 파괴 체크
#             bx, by = self.base
#             if self.map[bx][by] == 'S':
#                 print("Oh no, 기지가 파괴되었습니다… Game Over")
#                 return

#             # 5) 폭풍 제거, 기지 복원
#             self.clear_storm()

#             # 6) 다음 턴 예고 준비
#             self.plan_next_storm()

#             time.sleep(0.7)

#         # 모두 생존
#         self.print_map("모든 턴 생존! Victory!")
#         print()

# if __name__ == "__main__":
#     mode = input("난이도(easy, intermediate, hard): ").strip().lower()
#     SandstormGame(mode)



#통신
class MarsMissionApp:
    def __init__(self):
        self.password = "You'll never walk alone"
        self.help_signal = "HELP"
        self.special_text = (
            "Hey Watney, We are doing our best to save you. "
            "Supplies are on your way. Don't worry too much. "
            "The Password is (You'll never walk alone)"
        )

    @staticmethod
    def _encode(msg: str) -> str:
        return ' '.join(str(ord(c)) for c in msg)

    @staticmethod
    def _decode(code_str: str) -> str:
        chars = []
        for token in code_str.split():
            try:
                chars.append(chr(int(token)))
            except ValueError:
                chars.append('?')
        return ''.join(chars)

    def find_pathfinder(self):
        H, W = map(int, input("행과 열 개수 입력 (예: 5 5): ").split())
        matrix = []
        print(f"총 {H}줄을 입력하세요 (각 줄에 {W}개의 0~255 숫자):")
        for _ in range(H):
            row = list(map(int, input().split()))
            if len(row) != W:
                raise ValueError(f"한 줄에 {W}개의 숫자를 입력해야 합니다.")
            matrix.append(row)
        th_in = input("임계값 입력 [기본 200]: ").strip()
        th = int(th_in) if th_in else 200

        xs = []; ys = []
        for y, row in enumerate(matrix):
            for x, val in enumerate(row):
                if val >= th:
                    xs.append(x); ys.append(y)
        if not xs:
            print("밝은 픽셀을 찾지 못했습니다.")
        else:
            print(f"패스파인더 예상 위치: ({sum(xs)//len(xs)}, {sum(ys)//len(ys)})")

    def ascii_communicate(self):
        help_code = self._encode(self.help_signal)
        special_code = self._encode(self.special_text)

        while True:
            cmd = input("\n[H] HELP→ASCII  [SEND] SEND ASCII  [RECV] RECV ASCII  [Q] 뒤로\n>>> ").strip().upper()
            if cmd == 'Q':
                break
            elif cmd == 'H':
                print(f"→ 'HELP' ASCII: {help_code}")
            elif cmd == 'SEND':
                code = input("전송할 ASCII 코드 입력: ").strip()
                if code == help_code:
                    print("\n→ 특수 메시지 ASCII:\n" + special_code)
                else:
                    print("→ 등록된 신호가 아닙니다.")
            elif cmd == 'RECV':
                code = input("받은 ASCII 코드 입력: ").strip()
                text = self._decode(code)
                print(f"→ 복원된 메시지:\n{text}")
            else:
                print("올바른 옵션을 선택하세요.")

    def rover_chat(self):
        pw_input = input("채팅 모드 비밀번호를 입력하세요: ").strip()
        if pw_input != self.password:
            print("❌ 잘못된 비밀번호입니다.")
            return
        print("✅ 채팅 모드 진입 (종료는 'exit')\n")
        while True:
            msg = input("You: ").strip()
            if msg.lower() == 'exit':
                print("채팅 종료.")
                break
            if msg == "Hello":
                print(
                    "Rover: JPL: Mark, this is Vincent Kapoor. "
                    "We've been watching you since SOL54. "
                    "The whole world is rooting for you. "
                    "Amazing job, getting Pathfinder. "
                    "We're working on rescue plans. "
                    "Meantime we're putting together a supply mission "
                    "to keep you fed until Ares 4 arrives.\n"
                )
                continue
            reply = input("Rover 응답을 입력하세요: ").strip()
            print(f"Rover: {reply}\n")

    def run(self):
        while True:
            choice = input(
                "\n=== Mars Communication ===\n"
                "1. 패스파인더 위치 찾기\n"
                "2. ASCII 통신\n"
                "3. Rover 실시간 채팅\n"
                "Q. 종료\n>>> "
            ).strip().upper()
            if choice == '1':
                self.find_pathfinder()
            elif choice == '2':
                self.ascii_communicate()
            elif choice == '3':
                self.rover_chat()
            elif choice == 'Q':
                print("앱을 종료합니다.")
                break
            else:
                print("올바른 번호를 선택하세요.")

if __name__ == "__main__":
    app = MarsMissionApp()
    app.run()

