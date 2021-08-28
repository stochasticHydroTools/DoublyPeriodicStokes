export cpu = $(CPU)
export clonedir = $(PWD)

all: python

python:
	make -C python_interface
clean:
	make -C python_interface clean





