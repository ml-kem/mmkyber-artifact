#	Makefile

#	these are expected to come from command line
#	make MODE=PKE LEVEL=192 TEST=TESTVEC

#	options: { KEM, PKE }
MODE	?=	KEM

#	options: { 128, 192, 256 }
LEVEL	?=	128

#	options: { BENCH, TESTVEC }
TEST	?=	BENCH

#	configuration
XBIN	?=	xtest
CSRC	=	$(wildcard *.c sym/*.c)
OBJS	= 	$(CSRC:.c=.o)
CC 		?=	gcc
CFLAGS	+=	-Wall -Wextra -Wshadow -march=native -O3
#CFLAGS	+=	-Wall -Wextra -Wshadow -fsanitize=address,undefined -O2 -g
LDLIBS	+=	-lm
CFLAGS	+=	-Isym -DMM_$(MODE) -DMM_$(LEVEL) -D$(TEST)

$(XBIN): $(OBJS)
	$(CC) $(CFLAGS) -o $(XBIN) $(OBJS) $(LDLIBS)

%.o:	%.[csS]
	$(CC) $(CFLAGS) -c $^ -o $@

#	xtest
run:	$(XBIN)
	./$(XBIN)

test:
	./gen_testvec.sh
	cmp testvec.tmp testvec.txt

bench:
	./run_bench.sh

obj-clean:
	$(RM) -f $(XBIN) $(OBJS) gmon.out

clean:	obj-clean
	$(RM) -f testvec.tmp

