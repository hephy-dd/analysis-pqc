from setuptools import setup, find_packages

setup(
    name='analysis-pqc',
    version='0.1.1',
    packages=find_packages(exclude=['tests']),
    install_requires=[
        'numpy',
        'scipy'
    ],
    package_data={},
    entry_points={
        'console_scripts': [],
    },
    test_suite='tests'
)
