from setuptools import setup, find_packages

setup(
    name='analysis-pqc',
    version='0.7.0',
    packages=find_packages(exclude=['tests', 'scripts']),
    install_requires=[
        'numpy',
        'scipy'
    ],
    package_data={},
    entry_points={
        'console_scripts': [],
    }
)
