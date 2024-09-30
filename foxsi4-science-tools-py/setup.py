import setuptools

setuptools.setup(
    name="foxsi4_science_tools_py",
    version="0.0.1",
    description="Software used for FOXSI-4 science.",
    url="https://github.com/foxsi/foxsi4-science-tools",
    install_requires=[
            "setuptools",
            "pytest",
            "numpy", 
            "scipy",
            "matplotlib",
            "astropy",
            "sunpy",
            "pyyaml",
            "bs4",
            "lxml",
            "zeep",
            "drms",
        ],
    packages=setuptools.find_packages(),
    zip_safe=False
)
