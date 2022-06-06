from setuptools import setup, find_packages


setup(
    name = "HGD",
    version='0.15',
    description = "fluid dynamics.",
    long_description='long fluid dynamics.',
    url = "https://github.com/Cheshy123/StrucProg",
    author = "Cheshy123",
    author_email = "kerbal1@mail.ru",
    license = "-",
    packages=find_packages(),
    install_requires = open("requirements.txt", "r").read(),

)