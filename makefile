release:
	g++ main.cpp -o schedule.exe -O3 -std=c++17 -Wall -Wextra -Wpedantic -Wshadow
debug:
	g++ main.cpp -o schedule.exe -std=c++17 -Wall -Wextra -Wpedantic -Wshadow

clang:
	clang++ main.cpp -o schedule.exe -O3 -std=c++17 -Wall -Wextra -Wpedantic -Wshadow
