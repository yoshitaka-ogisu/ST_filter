from setuptools import setup, find_packages

setup(
    name='ST_filter',
    version='0.3',
    install_requires=[
        "numpy",
        "scipy",
        "pandas",
    ],
    python_requires='>=3',
    author="Yoshitaka Ogisu",
    author_email="yoshitaka.ogisu@gmail.com",
    description="A tool to find significant ties in a network",
    packages=find_packages()
)

