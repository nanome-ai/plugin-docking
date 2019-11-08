import os
from setuptools import find_packages, setup

with open(os.path.join(os.path.dirname(__file__), "README.md"), 'r') as file:
	README = file.read()

setup(
	name = 'nanome-docking',
	packages=find_packages(),
	version = '0.2.0',
	license='MIT',
	description = 'Nanome Plugin to dock ligands to a receptor',
	long_description = README,
    long_description_content_type = "text/markdown",
	author = 'Nanome',
	author_email = 'hello@nanome.ai',
	url = 'https://github.com/nanome-ai/plugin-docking',
	platforms="any",
	keywords = ['virtual-reality', 'chemistry', 'python', 'api', 'plugin'],
	install_requires=['nanome'],
	entry_points={"console_scripts": ["nanome-docking = nanome_docking.Docking:main"]},
	classifiers=[
		'Development Status :: 3 - Alpha',

		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Chemistry',

		'License :: OSI Approved :: MIT License',

		'Programming Language :: Python :: 2.7',
		'Programming Language :: Python :: 3.5',
		'Programming Language :: Python :: 3.6',
		'Programming Language :: Python :: 3.7',
	],
	package_data={
        "nanome_docking": [
            "*.json",
            "*.png",
			"smina"
        ]
	},
)