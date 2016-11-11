1. How to creat your virtural package so that you can use

`from my_package.my_code import my_class`

	1. Prepare the directory with certain layout:
		*root_package
		(Name doesn't matter)
			*`setup.py`
			*my_package
				(This is also a folder that contains your code)
				** `__init__.py`
					tells Python this is really a package, can be blank
				** `my_code.py`
					your_actually code

	2. In `setup.py`, write:

		`from setuptools import setup`
		`setup(name = 'Mysoftwere', packages=['my_package'])

	3. In the same directory as `setup.py`, run
		
		`python setup.py install`

	4. Now you can import your class anywhere.
