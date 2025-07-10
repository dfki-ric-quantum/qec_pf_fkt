.PHONY: all clean build generator_random_bond Z_to_txt test

SUBDIRS = generator_random_bond Z_to_txt test

all: generator_random_bond Z_to_txt test

generator_random_bond: | build
	@$(MAKE) --no-print-directory -C src/generator_random_bond

Z_to_txt: | build
	@$(MAKE) --no-print-directory -C src/Z_to_txt

test: | build
	@$(MAKE) --no-print-directory -C src/test

build:
	mkdir -p build/generator_random_bond build/Z_to_txt build/test

clean:
	@for dir in $(SUBDIRS); do \
		$(MAKE) --no-print-directory -C src/$$dir clean; \
	done
	@rm -rf build
	@echo "Cleaned build artifacts."
