install:
	python3 setup.py install

test:
	python3 -m unittest discover tests -p '*_test.py'

patch: install test
	bumpversion patch
	python3 setup.py sdist bdist_wheel
	python3 -m twine upload dist/* --skip-existing

minor: install test
	bumpversion minor
	python3 setup.py sdist bdist_wheel
	python3 -m twine upload dist/* --skip-existing