from distutils.core import setup


setup(
    name="abedeto",
    version="0.1.0",
    author="Marius Kriegerowski",
    author_email="marius@gfz-potsdam.de",
    description="Array Beam Depth Tool.",
    license="BSD",
    keywords="Array Beam Depth",
    packages=['abedeto'],
    package_dir={'abedeto': 'src'},
    long_description="Array Beam Depth Tool",
    install_requires=[
        "matplotlib",
        "pyrocko",
        "numpy",
    ],
    classifiers=[
        "Development Status :: 1"
        "Topic :: Utilities",
        "License :: OSI Approved :: BSD License",
    ],
    scripts=['src/abedeto'],
)
