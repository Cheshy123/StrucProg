# Модуль рассчёта трёх скважин.

# Для построения графиков и их исследований
import matplotlib.pyplot as plt
import math

# Востребователен для натуральных логарифмов и квадратных корней
from numpy import zeros

# Нужен для подгрузки данных, таких как  P, T, po и т.д. из готового файла, чтобы избежать бесконечных импутов
import toml


class SKV:
    def __init__(self):

        pass

    def load(self, filepath1, filepath2):

        # Подключение и интерпретация файла данных Data_HGD_input.toml и Data_HGD_output.toml
        self.DT = toml.load(filepath1)
        self.DT2 = toml.load(filepath2)

        # Блок ввода исходных параметров системы.

        # Количество шагов разбиения
        self.N = self.DT["N3"]
        # Динамическая вязкость воды, Па*с
        self.muN0 = self.DT["muN0"]
        # Диаметр трубопровода, м
        self.Dt = self.DT["Dt"]
        # Коэффициент обводненности
        self.alphav = self.DT["alphav"]
        # Шероховатость
        self.OTsheroh = self.DT["OTsheroh"]
        # Толщина стенки трубопровода
        self.thick = self.DT["thickT"]
        # Эмпирический коэффициент влияния температуры на плотность нефти 1/К
        self.betaN = self.DT["betaN"]
        # Эмпирический коэффициент влияния температуры на плотность воды 1/К
        self.betaV = self.DT["betaV"]
        # Удельная теплоемкость нефти
        self.CN = self.DT["CN"]
        # Удельная теплоемкость воды
        self.CV = self.DT["CV"]
        # Длина трубопроводов (1, 2, 3), м
        self.L = self.DT["L"]
        # Коэффициент теплопроводности грунта Вт/м * К
        self.lyambdaGr = self.DT["lyambdaGr"]
        # Коэффициент теплопрововодности нефти Вт/м * К
        self.lyambdaN0 = self.DT["lyambdaN0"]
        # Коэффициент теплопрововодности воды Вт/м * К
        self.lyambdaV = self.DT["lyambdaV"]
        # Коэффициент теплопроводности материала стенки трубы Вт/м * К
        self.lyambdaSt = self.DT["lyambdaSt"]
        # Удельная теплоёмкость материала стенки трубы
        self.CSt = self.DT["CSt"]
        # Плотность нефти при нормальных условиях
        self.po = self.DT["po"]
        # Коэффициент технического состояний труб
        self.E = self.DT["E"]
        # Коэффициент местных сопротивлений
        self.dzeta = self.DT["dzeta"]
        # Температура грунта на глубине (зависит от глубины пролегания трубопровода H), К
        self.Tgr = self.DT["Tgr"]
        # Глубина пролегания трубопроводов, м
        self.H = self.DT["H"]
        self.Qpl = self.DT2["Qpl"]
        self.Ppl = self.DT2["Ppl"]
        # Ускорение свободного падения
        self.g = self.DT["g"]
        # Критические числа Рейнольдса
        self.Re1 = self.DT["Re1"]
        self.Re2 = self.DT["Re2"]
        # Число пи
        self.pi = self.DT["pi"]
        # Плотность воды
        self.pov = self.DT["pov"]
        # Радиус скважины, м
        self.rc = self.DT["rc"]
        self.Hskv = self.DT["Hskv"]
        self.gradT = self.DT["gradT"]

        # Задание списков физических параметров для системы скважины
        self.v1 = [[0 for _ in range(self.N + 1)] for _ in range(3)]
        self.po1 = [[0 for _ in range(self.N + 1)] for _ in range(3)]
        self.T1 = [[0 for _ in range(self.N + 1)] for _ in range(3)]
        self.P1 = [[0 for _ in range(self.N + 1)] for _ in range(3)]
        self.G1 = [[0 for _ in range(self.N + 1)] for _ in range(3)]
        self.nu = [[0 for _ in range(self.N + 1)] for _ in range(3)]

        # Предварительные расчёты и интерпретация исходных данных
        self.Pc = self.DT["Pc"]

        for i in range(3):
            self.P1[i][0] = self.Pc[i]

        self.Qpl = self.DT2["Qpl"]
        self.Ppl = self.DT2["Ppl"]
        self.h = self.DT["h"]

        self.Dc = [0, 0, 0]
        self.dc = [0, 0, 0]
        self.S1 = [0, 0, 0]
        self.deltaZ = [0, 0, 0]
        self.lyambda2 = [0, 0, 0]
        self.z = []

    def dump(self, filepath):

        self.DT2 = {
            "Pskv": self.P1,
            "Gskv": self.G1,
            "Tskv": self.T1,
            "Ppl": self.Ppl,
            "Qpl": self.Qpl,
        }
        with open(filepath, "w") as io:
            toml.dump(self.DT2, io)

    def solve(self):

        for i in range(3):

            # Расчёт начальных температур
            self.T1[i][0] = self.Tgr + (self.Hskv[i] + self.h[i]) * self.gradT

            # Диаметры скважин

            # Внешний диаметр скважины
            self.Dc[i] = self.rc[i] * 2
            # внутренний диаметр
            self.dc[i] = self.Dc[i] - (2 * self.thick)

            # Площадь сечения скважины
            self.S1[i] = (self.pi * (self.dc[i] ** 2)) / 4

            # Рассчёт начальной плотности смеси
            self.po1[i][0] = (self.po / (1 + (self.betaN * (self.T1[i][0] - 293)))) * (
                1 - self.alphav
            ) + (self.pov / (1 + (self.betaV * (self.T1[i][0] - 293)))) * self.alphav

            # Начальных расход, берётся из рассчёта пласта, пересчёт в кг/с
            self.G1[i][0] = self.Qpl[i] * self.po1[i][0]

            # Рассчёт начальной скорости потока
            self.v1[i][0] = self.G1[i][0] / (self.po1[i][0] * self.S1[i])

            # Рассчёт шага разбиения
            self.deltaZ[i] = self.Hskv[i] / self.N

            # Задание списка глубины скважины
            Z = [self.deltaZ[i] * (self.N - _) for _ in range(self.N + 1)]
            self.z.append(Z)

            # Вычисление lyambda2 для каждой скважины
            self.lyambda2[i] = 0.11 * (((68 / self.Re2) + self.OTsheroh) ** 0.25)

        # Формула вязкости Эйнштейна для смеси воды и нефти
        muN = self.muN0 * (1 + 2.5 * self.alphav)

        lyambda1 = 64 / self.Re1

        # Рассчёт теплоёмкости смеси
        C = self.CN * (1 - self.alphav) + self.CV * self.alphav

        # Рассчёт теплопроводности смеси
        lyambdaN = self.lyambdaN0 * (1 - self.alphav) + self.lyambdaV * self.alphav

        # Рассчёт скважин по сечениям
        for i in range(3):
            for j in range(0, self.N):
                dzeta = 0
                # Формула для вычисления плотности смеси
                self.po1[i][j + 1] = (
                    self.po / (1 + (self.betaN * (self.T1[i][j] - 293)))
                ) * (1 - self.alphav) + (
                    self.pov / (1 + (self.betaV * (self.T1[i][j] - 293)))
                ) * self.alphav

                # Формула для вычисления кинематической вязкости смеси
                self.nu[i][j] = muN / self.po1[i][j]

                # Скорость и расход
                self.v1[i][j + 1] = (self.po1[i][j] * self.v1[i][j]) / (
                    self.po1[i][j + 1]
                )

                # Потери давления, определение числа Рейнольдса
                Re = (self.v1[i][j] * self.dc[i]) / self.nu[i][j]

                # Определение lyambdaT по условиям числа Рейнольдса
                if Re == self.Re2:
                    lyambdaT = self.lyambda2[i]
                # Для ламинарного режима
                elif Re <= self.Re1:
                    lyambdaT = lyambda1
                # Для переходного режима
                elif self.Re1 < Re <= self.Re2:
                    lyambdaT = lyambda1 + (
                        (self.lyambda2[i] - lyambda1) / (self.Re2 - self.Re1)
                    ) * (Re - self.Re1)
                elif self.Re2 < Re <= 500 / self.OTsheroh:
                    lyambdaT = 0.067 * (((158 / Re) + 2 * self.OTsheroh) ** 0.2)
                elif Re > 500 / self.OTsheroh:
                    lyambdaT = 0.067 * ((2.136 * self.OTsheroh) ** 0.2)

                lyambdaTr = (1.05 * lyambdaT) / (self.E**2)

                # Работа сил трения на участке deltaZ, ф. Вейсбаха-Дарси
                deltaPtr = (
                    lyambdaTr
                    * (self.deltaZ[i] / self.dc[i])
                    * self.po1[i][j]
                    * ((self.v1[i][j] ** 2) / 2)
                )

                # потери на колене 90, при выходе из скважины, и на открытой задвижке
                if j == self.N + 1:
                    dzeta = 1.37 + 0.15

                # потери на входе в трубу НКТ
                if j == 0:
                    dzeta = 0.5

                # Местные потери давления
                deltaPmest = dzeta * self.po1[i][j] * ((self.v1[i][j] ** 2) / 2)

                # Суммарные потери давления
                deltaP = deltaPtr + deltaPmest

                # Потери теплоты

                # Расчет температуры грунта по мере уменьшения координаты z
                TgrSkv = self.Tgr + self.z[i][j + 1] * self.gradT

                # Число Прандтля для нефти
                Pr = (muN * C) / lyambdaN

                # Число Прандтля для стенки трубы
                PrSt = 1

                # Число Грасгофа
                Gr = (
                    self.g
                    * ((self.deltaZ[i]) ** 3)
                    * self.betaN
                    * (self.T1[i][j] - TgrSkv)
                ) / self.nu[i][j] ** 2

                # Для ламинарного режима
                alphaL = (
                    0.17
                    * (lyambdaN / self.dc[i])
                    * Re**0.33
                    * Pr**0.43
                    * Gr**0.1
                    * 1
                )

                # Для турбулентного режима
                alphaT = 0.021 * (lyambdaN / self.dc[i]) * Re**0.8 * Pr**0.43 * 1

                if Re <= self.Re1:

                    # Для ламинарного режима
                    alphaSS = alphaL

                    # Коэффициент Кориолиса
                    alphak = 2

                elif Re >= self.Re2:

                    # Для турбулентного режима
                    alphaSS = alphaT

                    # Коэффициент Кориолиса
                    alphak = 1.1

                elif self.Re1 < Re < self.Re2:

                    # Для переходного режима
                    alphaSS = alphaL + ((alphaT - alphaL) / 8000) * (Re - 2000)

                    # Коэффициент Кориолиса
                    alphak = (-1.169e-4 * Re) + 2.269

                # Коэффициент теплопередачи К от теплоносителя в грунт для подземных трубопроводов с учетом стенки трубы, Вт/м**2 * K
                k = 1 / (
                    (1 / alphaSS)
                    + (
                        (self.dc[i] / 2 * self.lyambdaSt)
                        * math.log(self.Dc[i] / self.dc[i])
                    )
                    + ((self.dc[i] / 2 * self.lyambdaGr) * math.log(10))
                )

                # Тепловой поток в окружающую среду, Дж/с
                Qvn = (
                    k * self.pi * self.dc[i] * self.deltaZ[i] * (TgrSkv - self.T1[i][j])
                )

                # Уравнение теплового баланса
                self.T1[i][j + 1] = self.T1[i][j] + (Qvn / (C * self.G1[i][j]))

                # Уравнение Бернулли
                self.P1[i][j + 1] = self.po1[i][j + 1] * (
                    (self.P1[i][j] / self.po1[i][j])
                    + (alphak * ((self.v1[i][j] ** 2) / 2))
                    - (alphak * ((self.v1[i][j + 1] ** 2) / 2))
                    - (self.g * self.deltaZ[i])
                    - (deltaP / self.po1[i][j + 1])
                )

                self.G1[i][j + 1] = self.po1[i][j] * self.v1[i][j] * self.S1[i]

    def plot(self):

        # Построение графика

        self.z[0] = [self.z[0][_] * (-1) for _ in range(0, self.N + 1)]
        self.z[1] = [self.z[1][_] * (-1) for _ in range(0, self.N + 1)]
        self.z[2] = [self.z[2][_] * (-1) for _ in range(0, self.N + 1)]

        self.P1[0] = [self.P1[0][_] / 1e6 for _ in range(0, self.N + 1)]
        self.P1[1] = [self.P1[1][_] / 1e6 for _ in range(0, self.N + 1)]
        self.P1[2] = [self.P1[2][_] / 1e6 for _ in range(0, self.N + 1)]

        fig = plt.figure(1)
        ax = fig.add_subplot(111)
        ax.set_xlim([-(max(self.Hskv) + 100), 100])
        # заголовок
        plt.title("Распределение давления в скважине", fontsize=16)
        # ось абсцисс
        plt.xlabel("Радиус, м", fontsize=14)
        # ось ординат
        plt.ylabel("Давление, МПа", fontsize=14)
        # включение отображение сетки
        plt.grid(which="major", color="grey", linewidth=0.5)
        plt.plot(self.z[0], self.P1[0], c="orange", linewidth=2)
        plt.plot(self.z[1], self.P1[1], c="green", linewidth=2)
        plt.plot(self.z[2], self.P1[2], c="blue", linewidth=2)
        plt.legend(["Скважина №1", "Скважина №2", "Скважина №3"])

        fig = plt.figure(2)
        ax = fig.add_subplot(111)
        ax.set_xlim([-(max(self.Hskv) + 100), 100])
        # заголовок
        plt.title("Распределение температуры в скважине", fontsize=16)
        # ось абсцисс
        plt.xlabel("Радиус, м", fontsize=14)
        # ось ординат
        plt.ylabel("Температура, К", fontsize=14)
        # включение отображение сетки
        plt.grid(which="major", color="grey", linewidth=0.5)
        plt.plot(self.z[0], self.T1[0], c="orange", linewidth=2)
        plt.plot(self.z[1], self.T1[1], c="green", linewidth=2)
        plt.plot(self.z[2], self.T1[2], c="blue", linewidth=2)
        plt.legend(["Скважина №1", "Скважина №2", "Скважина №3"])

        fig = plt.figure(3)
        ax = fig.add_subplot(111)
        ax.set_xlim([-(max(self.Hskv) + 100), 100])
        # заголовок
        plt.title("Распределение плотности в скважине", fontsize=16)
        # ось абсцисс
        plt.xlabel("Радиус, м", fontsize=14)
        # ось ординат
        plt.ylabel("Плотность, кг/м^3", fontsize=14)
        # включение отображение сетки
        plt.grid(which="major", color="grey", linewidth=0.5)
        plt.plot(self.z[0], self.po1[0], c="orange", linewidth=2)
        plt.plot(self.z[1], self.po1[1], c="green", linewidth=2)
        plt.plot(self.z[2], self.po1[2], c="blue", linewidth=2)
        plt.legend(["Скважина №1", "Скважина №2", "Скважина №3"])
        plt.show()


if __name__ == "__main__":

    PL = SKV()
    PL.load(r"Project_HGD\Data_HGD_input.toml", r"Project_HGD\Data_HGD_output.toml")
    PL.solve()
    PL.dump(r"Project_HGD\Data_HGD_output.toml")
    PL.plot()

    print(PL)
