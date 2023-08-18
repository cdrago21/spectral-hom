"""Setup file for Hong-Ou-Mandel-Interference module module."""

import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

with open('requirements.txt') as fp:
    install_requires = fp.read()

setuptools.setup(
    name="Hong-Ou-Mandel-Interference",
    description="Hong-Ou-Mandel-Interference",
    long_description=long_description,
    packages=setuptools.find_packages(),
    install_requires=install_requires,
    python_requires='>=3.6'
)
