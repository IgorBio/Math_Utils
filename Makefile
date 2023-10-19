BUILD_DIR = ./build
TEST_BUILD_DIR = $(BUILD_DIR)/tests
RUN_TEST=$(TEST_BUILD_DIR)/MathUtilsTests
OS = $(shell uname)

ifeq ($(OS), Linux)
	CHECK_LEAKS=CK_FORK=no valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --log-file=valgrind.log
	OPEN=xdg-open
else
	CHECK_LEAKS=CK_FORK=no leaks --atExit --
	OPEN=open
endif


all: build

build:
	@cmake -S . -B $(BUILD_DIR)
	@cmake --build $(BUILD_DIR)

rebuild: clean build

cppcheck:
	@cmake --build $(BUILD_DIR) --target cppcheck

style:
	@cmake --build $(BUILD_DIR) --target clang-format

tests:
	@cmake -S ./tests -B $(TEST_BUILD_DIR)
	@cmake --build $(TEST_BUILD_DIR)
	@$(RUN_TEST)

gcov_report: tests
	@cmake --build $(TEST_BUILD_DIR) --target coverage

leaks: tests
	@$(CHECK_LEAKS) $(RUN_TEST)

clean:
	@rm -rf $(BUILD_DIR)

.PHONY: all build rebuild clean cppcheck style tests lcov