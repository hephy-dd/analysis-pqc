from setuptools import setup, find_packages

setup(
    name='analysis-pqc',
    version='0.1.0',
    packages=find_packages(exclude=['tests']),
    install_requires=[
        'numpy',
        'scipy',
        'matplotlib'
    ],
    package_data={},
    entry_points={
        'console_scripts': [],
    },
    test_suite='tests'
)
