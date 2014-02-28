ifndef ANALYZER
$(error $$ANALYZER environment variable not defined)
endif

all:
	analyzer -n -q build.C

clean:
	rm -f *.d *.so
