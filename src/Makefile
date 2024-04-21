CC = gcc
CFLAGS = -Wall -Werror -Wextra -std=c11

OS = $(shell uname)

ifeq ($(OS),Linux)
    det_OS = -lcheck  -lm -lrt -lpthread -lsubunit
else
    det_OS = -lcheck -lm -lpthread
endif

all : clean s21_matrix.a test gcov_report

test: 
	$(CC) $(CFLAGS) --coverage s21_matrix.c test.c -o test $(det_OS)
	./test

s21_matrix.o: s21_matrix.c s21_matrix.h
	$(CC) $(CFLAGS) s21_matrix.c  -c 

s21_matrix.a : s21_matrix.o
	ar rcs s21_matrix.a *.o
	ranlib s21_matrix.a

gcov_report: test
	lcov -t "result" -o res.info -c -d .
	genhtml -o report res.info
	open report/index.html

check_style:
	cp ../materials/linters/.clang-format ./
	clang-format -style=Google -n *.c *.h
	rm .clang-format

fix_style:
	cp ../materials/linters/.clang-format ./
	clang-format -style=Google -i *.c *.h
	rm .clang-format

clean:
	rm -rf report *.gc* *.o *.info *.a test.dSYM test
	
rebuild: clean all

