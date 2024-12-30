#	Makefile
#	===	Kyber reference benchmarks

#	options: { 128, 192, 256 }
LEVEL	?=	128

#	directory
KYBER	=	kyber/ref

#	configuration
XBIN	?=	xtest

CSRC	=	$(KYBER)/kem.c		$(KYBER)/indcpa.c	$(KYBER)/polyvec.c	\
			$(KYBER)/poly.c		$(KYBER)/ntt.c		$(KYBER)/cbd.c 		\
			$(KYBER)/reduce.c	$(KYBER)/verify.c	$(KYBER)/randombytes.c

CSRC	+=	$(KYBER)/fips202.c	$(KYBER)/symmetric-shake.c

CSRC	+=	$(wildcard *.c)

OBJS	= 	$(CSRC:.c=.o)
CC 		?=	gcc
CFLAGS	+=	-Wall -Wextra -Wshadow -march=native -O3
#CFLAGS	+=	-Wall -Wextra -Wshadow -fsanitize=address,undefined -O2 -g
LDLIBS	+=
CFLAGS	+=	-I$(KYBER) -DKYBER_K="($(LEVEL)/64)"

#	for N=1000 benchmarks
#CFLAGS	+=	-DMM_REP_TOT=100000 -DMM_SET_N=1000

$(XBIN): $(OBJS)
	$(CC) $(CFLAGS) -o $(XBIN) $(OBJS) $(LDLIBS)

%.o:	%.[csS]
	$(CC) $(CFLAGS) -c $^ -o $@

#	xtest
run:	$(XBIN)
	./$(XBIN)

bench:
	./run_bench.sh

obj-clean:
	$(RM) -f $(XBIN) $(OBJS) gmon.out

clean:	obj-clean
	$(RM) -f testvec.tmp
	#cd unit-test && $(MAKE) clean

