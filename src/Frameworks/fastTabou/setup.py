from setuptools import setup, find_packages

setup(
    name="fastTabou",
    version="0.1.0",
    author="KEFSI nOURHANE & Mekhazni Ryham",
    author_email="nourhanekefsi@gmail.com",
    description="Framework for Tabu Search with sequential and parallel implementations",
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    url="https://github.com/nourhanekefsi/Master_final_project/src/Framework/fastTabou",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
    install_requires=[
        'numpy>=1.19.0',
    ],
)