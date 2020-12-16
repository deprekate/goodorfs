import os
#from distutils.core import setup, Extension
from setuptools import setup, Extension, find_packages


def readme():
	with open("README.md", "r") as fh:
		long_desc = fh.read()
	return long_desc

def get_version():
    with open("VERSION", 'r') as f:
        v = f.readline().strip()
        return v

def main():
	setup (
		name = 'goodorfs',
		version = get_version(),
		author = "Katelyn McNair",
		author_email = "deprekate@gmail.com",
		description = 'A a tool to find coding orfs',
		long_description = readme(),
		long_description_content_type="text/markdown",
		url =  "https://github.com/deprekate/goodorfs",
		scripts=['goodorfs.py'],
		classifiers=[
			"Programming Language :: Python :: 3",
			"License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
			"Operating System :: OS Independent",
		],
		python_requires='>3.5.2',
		packages=find_packages(),
		install_requires=['sklearn','numpy']
	)


if __name__ == "__main__":
	main()
