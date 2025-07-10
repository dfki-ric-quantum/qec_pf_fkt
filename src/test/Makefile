all:
	$(info Running Tests)
	@../../build/Z_to_txt/isingZToTxt 4096 5 5 42 0.0 0.1 .
	@diff resultsGaussian/0.000000/0.000000/5/5/0.100000/4096/42/Z.txt expectedResults/0.000000/0.000000/5/5/0.100000/4096/Z.txt || \
		(echo "Failed test: 0.1 temperature Z calculation" && exit 1)
	@echo "Passed test: 0.1 temperature Z calculation"

	@../../build/Z_to_txt/isingZToTxt 4096 5 5 42 0.0 99999999999999.000000 .
	@./compare_txt_files.py resultsGaussian/0.000000/0.000000/5/5/99999999999999.000000/4096/42/Z.txt expectedResults/0.000000/0.000000/5/5/99999999999999.000000/4096/Z.txt || \
		(echo "Failed test: High temperature Z calculation" && exit 1)
	@echo "Passed test: High temperature Z calculation"

	@../../build/Z_to_txt/isingZToTxt 2048 4 4 42 0 `cat T_DOS_0.00000001` .
	@./compare_txt_files.py resultsGaussian/0.000000/0.000000/4/4/0.108574/2048/42/Z.txt expectedResults/0.000000/0.000000/4/4/lowT/Z.txt || \
		(echo "Failed test: Low temperature limit Z calculation" && exit 1)
	@echo "Passed test: Low temperature limit Z calculation"

	@../../build/generator_random_bond/isingGeneratorRandomBond 4 4 42 0.000001 . 0.000001
	@diff interactionsGaussian/0.000001/4/4/0.000001/42/interaction_lattice.txt expectedResults/0.000001/4/4/0.000001/42/interaction_lattice.txt || \
		(echo "Failed test: Non uniform noise Z calculation" && exit 1)
	@echo "Passed test: Non uniform noise Z calculation"

	@rm -rf resultsGaussian
	@rm -rf interactionsGaussian/0.000001

clean:
	@rm -f test_*
	@rm -rf tmp_test
