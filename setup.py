import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="DynaTMT-py-SB", 
    version="2.9.2",
    author="Kevin Klann",
    author_email="klann@em.uni-frankfurt.de",
    # Updated by = "Süleyman Bozkurt",
    # Updated date = '03/06/2024'
    description="Python package to analyse pSILAC TMT data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/science64/DynaTMT-py-SB",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.8',
)