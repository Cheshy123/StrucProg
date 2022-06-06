#Модуль рассчёта пласта, последовательный рассчёт сразу трёх пластов.
import matplotlib.pyplot as plt #Для построения графиков и их исследований
import math #Востребователен для натуральных логарифмов и квадратных корней
import toml #Нужен для подгрузки данных, таких как  P, T, po и т.д. из готового файла

class Plast:
    def __init__(self):
        pass
    def load(self, filepath):
        #Подключение и интерпретация файла данных Data_HGD_input.toml
        self.DT = toml.load(filepath)
        #Блок ввода исходных параметров системы.
        self.N = self.DT["N1"] #Количество шагов разбиения
        self.h = self.DT["h"] #Мощность пласта, м
        self.k = self.DT["k"] #Проницаемость пласта, м**2
        self.Pc = self.DT["Pc"] #Забойное давление (в начале скважины), Па
        self.Pk = self.DT["Pk"] #Давление на контуре питания, Па
        self.muN0 = self.DT["muN0"] #Динамическая вязкость воды, Па*с
        self.rc = self.DT["rc"] #Радиус скважины, м
        self.rk = self.DT["rk"] #Радиус контура питания, м
        self.alphav = self.DT["alphav"] #Коэффициент обводненности
        self.pi = self.DT["pi"] #Число пи
    def dump(self, filepath):
        self.DT2 = {"Ppl": self.P, "Qpl": self.Q}
        with open(filepath, "w") as io:
            toml.dump(self.DT2, io)
    def solve(self):   
        #Задание списков физических параметров для системы пласта
        self.P = [[0 for _ in range(self.N+1)] for _ in range(3)]
        self.wr = [[0 for _ in range(self.N+1)] for _ in range(3)]
        self.r = [[0 for _ in range(self.N+1)] for _ in range(3)]
        self.Q = [0 for _ in range(3)]
        #Формула Эйнштейна для вязкости эмульсии
        muN = self.muN0 * (1 + 2.5 * self.alphav)

        for i in range (3):
            self.Q[i] = ((2 * self.pi * self.h[i] * self.k[i] * (self.Pk[i] - self.Pc[i])) / (muN * math.log(self.rk[i] / self.rc[i]))) #Формула Дюпюи, расход объемный

            for j in range(0, self.N+1):
                #Увеличение радиуса с каждой итерацией
                if j == 0:
                    self.r[i][j] = self.rc[i]
                else:
                    self.r[i][j] = (self.rk[i] / self.N) * j
                self.P[i][j] = self.Pk[i] - (((self.Pk[i] - self.Pc[i]) / math.log(self.rk[i]/self.rc[i])) * math.log(self.rk[i] / self.r[i][j])) #Распределение давления в пласте
                
                self.wr[i][j] = (self.k[i] * (self.Pc[i] - self.Pk[i])) / (muN * math.log(self.rk[i] / self.rc[i]) * self.r[i][j]) #Распределение скоростей фильтрации

        print("Суммарный дебит, Кг**3/сут:", (self.Q[0] + self.Q[1] + self.Q[2]) * 850 * 60 * 60 *24)
        print("Дебит с 1 скважины, Кг**3/с:", self.Q[0] * 850)
        print("Дебит со 2 скважины, Кг**3/с:", self.Q[1] * 850)
        print("Дебит с 3 скважины, Кг**3/с:", self.Q[2] * 850)



    def plot(self):

        # Построение графиков

        self.P[0] = [self.P[0][_] / 1e+6 for _ in range(0, self.N + 1)]
        self.P[1] = [self.P[1][_] / 1e+6 for _ in range(0, self.N + 1)]
        self.P[2] = [self.P[2][_] / 1e+6 for _ in range(0, self.N + 1)]

        plt.figure(1)
        plt.title("Распределение давления в пласте", fontsize=16) # заголовок
        plt.xlabel("Радиус контура питания, м", fontsize=14) # ось абсцисс
        plt.ylabel("Давление, МПа", fontsize=14) # ось ординат
        plt.grid(which = "major", color = "grey", linewidth = 0.5)      # включение отображение сетки
        plt.plot(self.r[0], self.P[0], c = "orange", linewidth = 2)
        plt.plot(self.r[1], self.P[1], c = "green", linewidth = 2)
        plt.plot(self.r[2], self.P[2], c = "blue", linewidth = 2)
        plt.legend(["Пласт №1", "Пласт №2", "Пласт №3"])
        plt.show()
if __name__ == "__main__":
    PL=Plast()
    PL.load(r"Data_HGD_input.toml")
    PL.solve()
    PL.dump(r"Data_HGD_output.toml")
    PL.plot()
    print(PL)