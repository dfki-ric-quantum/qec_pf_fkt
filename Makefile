all: generator_random_bond Z_to_txt test

generator_random_bond:
	make -C src/generator_random_bond/

Z_to_txt:
	make -C src/Z_to_txt/

test:
	make -C src/test/

clean:
	make -C src/generator_random_bond/ clean
	make -C src/Z_to_txt/ clean
	make -C src/test/ clean
