from driverpower import __version__
from setuptools import setup, find_packages

setup(
    name='DriverPower',
    version=__version__,
    author='Shimin Shuai',
    author_email='sshuai@oicr.on.ca',
    url='https://github.com/smshuai/DriverPower',
    license='GPL',
    python_requires='>= 3.5.2',
    packages=find_packages(),
    description='Combined burden and functional test for coding and noncoding cancer drivers',
    py_modules=['driverpower'],
    install_requires=[
        'numpy >= 1.11.1',
        'scipy >= 0.18.1',
        'pandas >= 0.18.1',
        'scikit-learn >= 0.18',
        'statsmodels >= 0.6.1',
        'pytabix >= 0.0.2',
    ],
    entry_points = {
        'console_scripts': ['driverpower=driverpower.cmdline:main'],
    },
)

