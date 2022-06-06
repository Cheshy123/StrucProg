# Модуль рассчёта пласта, последовательный рассчёт сразу трёх пластов.

# Для построения графиков и их исследований
import matplotlib.pyplot as plt

# Востребователен для натуральных логарифмов и квадратных корней
import math

# Нужен для подгрузки данных, таких как  P, T, po и т.д. из готового файла
import toml


class Plast:
    def __init__(self):

        pass

    def load(self, filepath):

        # Подключение и интерпретация файла данных Data_HGD_input.toml
        self.DT = toml.load(filepath)

        # Блок ввода исходных параметров системы.

        # Количество шагов разбиения
        self.N = self.DT["N1"]
        # Мощность пласта, м
        self.h = self.DT["h"]
        # Проницаемость пласта, м**2
        self.k = self.DT["k"]
        # Забойное давление (в начале скважины), Па
        self.Pc = self.DT["Pc"]
        # Давление на контуре питания, Па
        self.Pk = self.DT["Pk"]
        # Динамическая вязкость воды, Па*с
        self.muN0 = self.DT["muN0"]
        # Радиус скважины, м
        self.rc = self.DT["rc"]
        # Радиус контура питания, м
        self.rk = self.DT["rk"]
        # Коэффициент обводненности
        self.alphav = self.DT["alphav"]
        # Число пи
        self.pi = self.DT["pi"]

    def dump(self, filepath):

        self.DT2 = {"Ppl": self.P, "Qpl": self.Q}

        with open(filepath, "w") as io:
            toml.dump(self.DT2, io)

    def solve(self):

        # Задание списков физических параметров для системы пласта
        self.P = [[0 for _ in range(self.N + 1)] for _ in range(3)]
        self.wr = [[0 for _ in range(self.N + 1)] for _ in range(3)]
        self.r = [[0 for _ in range(self.N + 1)] for _ in range(3)]
        self.Q = [0 for _ in range(3)]

        # Формула Эйнштейна для вязкости эмульсии
        muN = self.muN0 * (1 + 2.5 * self.alphav)

        # Формула Дюпюи, расход объемный
        for i in range(3):
            self.Q[i] = (
                2 * self.pi * self.h[i] * self.k[i] * (self.Pk[i] - self.Pc[i])
            ) / (muN * math.log(self.rk[i] / self.rc[i]))

            for j in range(0, self.N + 1):
                # Увеличение радиуса с каждой итерацией
                if j == 0:
                    self.r[i][j] = self.rc[i]
                else:
                    self.r[i][j] = (self.rk[i] / self.N) * j

                # Распределение давления в пласте
                self.P[i][j] = self.Pk[i] - (
                    ((self.Pk[i] - self.Pc[i]) / math.log(self.rk[i] / self.rc[i]))
                    * math.log(self.rk[i] / self.r[i][j])
                )

                # Распределение скоростей фильтрации
                self.wr[i][j] = (self.k[i] * (self.Pc[i] - self.Pk[i])) / (
                    muN * math.log(self.rk[i] / self.rc[i]) * self.r[i][j]
                )

        print(
            "Суммарный дебит, Кг**3/сут:",
            (self.Q[0] + self.Q[1] + self.Q[2]) * 850 * 60 * 60 * 24,
        )
        print("Дебит с 1 скважины, Кг**3/с:", self.Q[0] * 850)
        print("Дебит со 2 скважины, Кг**3/с:", self.Q[1] * 850)
        print("Дебит с 3 скважины, Кг**3/с:", self.Q[2] * 850)

    def plot(self):

        # Построение графиков

        self.P[0] = [self.P[0][_] / 1e6 for _ in range(0, self.N + 1)]
        self.P[1] = [self.P[1][_] / 1e6 for _ in range(0, self.N + 1)]
        self.P[2] = [self.P[2][_] / 1e6 for _ in range(0, self.N + 1)]

        plt.figure(1)
        # заголовок
        plt.title("Распределение давления в пласте", fontsize=16)
        # ось абсцисс
        plt.xlabel("Радиус контура питания, м", fontsize=14)
        # ось ординат
        plt.ylabel("Давление, МПа", fontsize=14)
        # включение отображение сетки
        plt.grid(which="major", color="grey", linewidth=0.5)
        plt.plot(self.r[0], self.P[0], c="orange", linewidth=2)
        plt.plot(self.r[1], self.P[1], c="green", linewidth=2)
        plt.plot(self.r[2], self.P[2], c="blue", linewidth=2)
        plt.legend(["Пласт №1", "Пласт №2", "Пласт №3"])
        plt.show()


if __name__ == "__main__":

    PL = Plast()
    PL.load(r"Project_HGD\Data_HGD_input.toml")
    PL.solve()
    PL.dump(r"Project_HGD\Data_HGD_output.toml")
    PL.plot()

    print(PL)
