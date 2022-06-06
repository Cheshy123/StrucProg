# Модуль рассчёта трубопроводов от 1, 2, 3 скважин.

# Для построения графиков и их исследований
import matplotlib.pyplot as plt
import math

# Востребователен для натуральных логарифмов и квадратных корней
from numpy import zeros

# Нужен для подгрузки данных, таких как  P, T, po и т.д. из готового файла, чтобы избежать бесконечных импутов
import toml


class TrO:
    def __init__(self):

        pass

    def load(self, filepath1, filepath2):

        # Подключение и интерпретация файла данных Data_HGD_input.toml и Data_HGD_output.toml
        self.DT = toml.load(filepath1)
        self.DT2 = toml.load(filepath2)

        # Блок ввода исходных параметров системы.

        # Количество шагов разбиения
        self.N = self.DT["N4"]
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
        # Длина общего трубопровода, м
        self.LTrO = self.DT["LTrO"]
        # Глубины скважин, м
        self.Hskv = self.DT["Hskv"]
        # Длина трубопровода до узла соединения трубопроводов от скважин, м
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
        # Глубина пролегания трубопровода, м
        self.H = self.DT["H"]
        # Ускорение свободного падения
        self.g = self.DT["g"]
        # Критические числа Рейнольдса
        self.Re1 = self.DT["Re1"]
        self.Re2 = self.DT["Re2"]
        # Число пи
        self.pi = self.DT["pi"]
        # Плотность воды
        self.pov = self.DT["pov"]
        self.Hskv = self.DT["Hskv"]
        self.Pskv = self.DT2["Pskv"]
        self.Gskv = self.DT2["Gskv"]
        self.gradT = self.DT["gradT"]
        self.Tskv = self.DT2["Tskv"]
        self.Ppl = self.DT2["Ppl"]
        self.Qpl = self.DT2["Qpl"]
        self.Ptr = self.DT2["Ptr"]

        # Задание списков физических параметров для системы скважины
        self.v1 = [0 for _ in range(self.N + 1)]
        self.po1 = [0 for _ in range(self.N + 1)]
        self.T1 = [0 for _ in range(self.N + 1)]
        self.P1 = [0 for _ in range(self.N + 1)]
        self.G1 = [0 for _ in range(self.N + 1)]
        self.nu = [0 for _ in range(self.N + 1)]

    def dump(self, filepath):

        self.DT2 = {
            "Pskv": self.Pskv,
            "Gskv": self.Gskv,
            "Tskv": self.Tskv,
            "Ppl": self.Ppl,
            "Qpl": self.Qpl,
            "Ptr": self.P1,
            "Gtr": self.G1,
            "Ttr": self.T1,
        }
        with open(filepath, "w") as io:
            toml.dump(self.DT2, io)

    def solve(self):

        # Внутренний диаметр трубопровода

        # внутренний диаметр
        dt = self.Dt - (2 * self.thick)

        # Площадь сечения трубопроводов
        self.S1 = (self.pi * (dt**2)) / 4

        # Глубина пролегания трубопроводов, до центра сечения трубы, м
        self.H = self.H + (self.Dt / 2)

        # Предварительные расчёты и интерпретация исходных данных
        self.Pskv = self.DT2["Pskv"]
        self.Ptr = self.DT2["Ptr"]
        self.Gtr = self.DT2["Gtr"]

        # Начальных расход, берётся из рассчёта скважин, кг/с
        self.G1[0] = self.Gtr[0][-1] + self.Gtr[1][-1] + self.Gtr[2][-1]
        Ttr = self.DT2["Ttr"]

        # Расчёт начальной температуры, с учётом подвода массы
        self.T1[0] = (
            Ttr[0][-1] * self.Gtr[0][-1]
            + Ttr[1][-1] * self.Gtr[1][-1]
            + Ttr[2][-1] * self.Gtr[2][-1]
        ) / self.G1[0]

        # Расчет итогового давления, получаемого общим трубопроводом
        self.P1[0] = 4.2e6

        # Рассчёт начальной плотности смеси
        self.po1[0] = (self.po / (1 + (self.betaN * (self.T1[0] - 293)))) * (
            1 - self.alphav
        ) + (self.pov / (1 + (self.betaV * (self.T1[0] - 293)))) * self.alphav

        # Рассчёт начальной скорости потока
        self.v1[0] = self.G1[0] / (self.po1[0] * self.S1)

        # Рассчёт шага разбиения
        deltaX = self.LTrO / self.N

        # Задание списка разбиения трубопроводов
        self.x = [deltaX * _ for _ in range(self.N + 1)]

        lyambda2 = 0.11 * (((68 / self.Re2) + self.OTsheroh) ** 0.25)

        # Формула вязкости Эйнштейна для смеси воды и нефти
        muN = self.muN0 * (1 + 2.5 * self.alphav)

        lyambda1 = 64 / self.Re1

        # Рассчёт теплоёмкости смеси
        C = self.CN * (1 - self.alphav) + self.CV * self.alphav

        # Рассчёт теплопроводности смеси
        lyambdaN = self.lyambdaN0 * (1 - self.alphav) + self.lyambdaV * self.alphav

        # Рассчёт трубопроводов по сечениям

        for j in range(0, self.N):

            dzeta = 0

            # Уравнение плотности
            self.po1[j + 1] = (self.po / (1 + (self.betaN * (self.T1[j] - 293)))) * (
                1 - self.alphav
            ) + (self.pov / (1 + (self.betaV * (self.T1[j] - 293)))) * self.alphav

            self.nu[j] = muN / self.po1[j]

            # Скорость и расход
            self.v1[j + 1] = (self.po1[j] * self.v1[j]) / (self.po1[j + 1])

            # Потери давления, определение числа Рейнольдса
            Re = (self.v1[j] * self.Dt) / self.nu[j]

            # Определение lyambdaT по условиям числа Рейнольдса
            if Re == self.Re2:
                lyambdaT = lyambda2

            # Для ламинарного режима
            elif Re <= self.Re1:
                lyambdaT = lyambda1

            # Для переходного режима
            elif self.Re1 < Re <= self.Re2:
                lyambdaT = lyambda1 + (
                    (lyambda2 - lyambda1) / (self.Re2 - self.Re1)
                ) * (Re - self.Re1)
            elif self.Re2 < Re <= 500 / self.OTsheroh:
                lyambdaT = 0.067 * (((158 / Re) + 2 * self.OTsheroh) ** 0.2)
            elif Re > 500 / self.OTsheroh:
                lyambdaT = 0.067 * ((2.136 * self.OTsheroh) ** 0.2)

            lyambdaTr = (1.05 * lyambdaT) / (self.E**2)

            # потери на тройнике
            if j == 0:
                dzeta = 0.23

            # Работа сил трения на участке deltaX, ф. Вейсбаха-Дарси
            deltaPtr = (
                lyambdaTr * (deltaX / self.Dt) * self.po1[j] * ((self.v1[j] ** 2) / 2)
            )

            # Местные потери давления
            deltaPmest = dzeta * self.po1[j] * ((self.v1[j] ** 2) / 2)

            # Суммарные потери давления
            deltaP = deltaPtr + deltaPmest

            # Подвод теплоты

            # Коэффициент теплоотдачи alphaGr от стенки трубопровода к грунту (формула Форхгеймера - Власова), Вт/м**2 * K
            alphaGr = (2 * self.lyambdaGr) / (
                self.Dt
                * (
                    math.log((2 * self.H) / self.Dt)
                    + math.sqrt((((2 * self.H) / self.Dt) ** 2) - 1)
                )
            )

            # Число Прандтля для нефти
            Pr = (muN * C) / lyambdaN

            # Число Грасгофа
            Gr = (
                self.g * ((deltaX) ** 3) * self.betaN * (self.T1[j] - self.Tgr)
            ) / self.nu[j] ** 2

            # Для ламинарного режима
            alphaL = (
                0.17 * (lyambdaN / self.Dt) * Re**0.33 * Pr**0.43 * Gr**0.1 * 1
            )

            # Для турбулентного режима
            alphaT = 0.021 * (lyambdaN / self.Dt) * Re**0.8 * Pr**0.43 * 1

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

            # Коэффициент теплопередачи К от теплоносителя в грунт для подземных трубопроводов, Вт/м**2 * K
            k = 1 / (
                (1 / (alphaSS * self.Dt))
                + ((1 / (2 * lyambdaN)) * math.log(self.Dt / self.Dt))
                + (1 / (alphaGr * self.Dt))
            )

            # Тепловой поток в окружающую среду, Дж/с
            Qvn = k * self.pi * self.Dt * deltaX * (self.Tgr - self.T1[j])

            # Уравнение теплового баланса
            self.T1[j + 1] = self.T1[j] + (Qvn / (C * self.G1[j]))

            # Уравнение Бернулли
            self.P1[j + 1] = self.po1[j + 1] * (
                (self.P1[j] / self.po1[j])
                + (alphak * ((self.v1[j] ** 2) / 2))
                - (alphak * ((self.v1[j + 1] ** 2) / 2))
                - (deltaP / self.po1[j + 1])
            )

            self.G1[j + 1] = self.po1[j] * self.v1[j] * self.S1

    def plot(self):

        # График общего давления системы

        # Задание списка разбиения трубопроводов
        self.Pskv[0] = self.Pskv[0] + self.Ptr[0]
        self.Pskv[0] = self.Pskv[0] + self.P1

        self.Pskv[1] = self.Pskv[1] + self.Ptr[1]
        self.Pskv[1] = self.Pskv[1] + self.P1

        self.Pskv[2] = self.Pskv[2] + self.Ptr[2]
        self.Pskv[2] = self.Pskv[2] + self.P1

        self.Pskv[0] = [self.Pskv[0][_] / 1e6 for _ in range(0, (self.N + 1) * 3)]

        self.Pskv[1] = [self.Pskv[1][_] / 1e6 for _ in range(0, (self.N + 1) * 3)]

        self.Pskv[2] = [self.Pskv[2][_] / 1e6 for _ in range(0, (self.N + 1) * 3)]

        X = [1 * _ for _ in range((self.N + 1) * 3)]
        print(self.Pskv[0])

        plt.figure(0)
        # заголовок
        plt.title("Распределение давления в системе", fontsize=16)
        # ось абсцисс
        plt.xlabel("Номер сечения", fontsize=14)
        # ось ординат
        plt.ylabel("Давление, МПа", fontsize=14)
        # включение отображение сетки
        plt.grid(which="major", color="grey", linewidth=0.5)
        plt.plot(X, self.Pskv[0], c="orange", linewidth=2)
        plt.plot(X, self.Pskv[1], c="green", linewidth=2)
        plt.plot(X, self.Pskv[2], c="blue", linewidth=2)

        # Построение графика
        self.P1 = [self.P1[_] / 1e6 for _ in range(0, self.N + 1)]

        plt.figure(1)
        # заголовок
        plt.title("Распределение давления в трубопроводе", fontsize=16)
        # ось абсцисс
        plt.xlabel("Длина трубопровода, м", fontsize=14)
        # ось ординат
        plt.ylabel("Давление, МПа", fontsize=14)
        # включение отображение сетки
        plt.grid(which="major", color="grey", linewidth=0.5)
        plt.plot(self.x, self.P1, c="green", linewidth=2)
        plt.legend(["Общий трубопровод"])

        plt.figure(2)
        # заголовок
        plt.title("Распределение температуры в трубопроводе", fontsize=16)
        # ось абсцисс
        plt.xlabel("Длина трубопровода, м", fontsize=14)
        # ось ординат
        plt.ylabel("Темпераутра, К", fontsize=14)

        # включение отображение сетки
        plt.grid(which="major", color="grey", linewidth=0.5)
        plt.plot(self.x, self.T1, c="green", linewidth=2)
        plt.legend(["Общий трубопровод"])

        plt.figure(3)
        # заголовок
        plt.title("Распределение расхода в трубопроводе", fontsize=16)
        # ось абсцисс
        plt.xlabel("Длина трубопровода, м", fontsize=14)
        # ось ординат
        plt.ylabel("Расход ,кг/с", fontsize=14)
        # включение отображение сетки
        plt.grid(which="major", color="grey", linewidth=0.5)
        plt.plot(self.x, self.G1, c="green", linewidth=2)
        plt.legend(["Общий трубопровод"])
        plt.show()


if __name__ == "__main__":

    PL = TrO()
    PL.load(r"Project_HGD\Data_HGD_input.toml", r"Project_HGD\Data_HGD_output.toml")
    PL.solve()
    PL.dump(r"Project_HGD\Data_HGD_output.toml")
    PL.plot()

    print(PL)
