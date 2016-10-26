from setuptools import setup, find_packages
setup(
    name='DriverPower',
    version='0.2',
    author='Shimin Shuai',
    author_email='sshuai@oicr.on.ca',
    url='',
    packages=find_packages('src'),
    package_dir = {'driverpower': 'src'},
    py_modules=['driverpower'],
    python_requires='>= 3.5.2',
    install_requires=[
        'numpy >= 1.11.1',
        'scipy >= 0.18.1',
        'pandas >= 0.18.1',
        'scikit-learn >= 0.17',
        'statsmodels >= 0.6.1',

    ],
)

