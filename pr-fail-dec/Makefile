#	Makefile

XBIN	=	xtest
CC		=	gcc
CFLAGS	=	-Wall -Wextra -O3 -march=native
#CFLAGS	=	-Wall -Wextra -Wshadow -fsanitize=address,undefined -O2 -g 
CSRC	= 	$(wildcard *.c)
OBJS	= 	$(CSRC:.c=.o)
LDLIBS	=	-lmpfr -lgmp -lm

#	Standard Linux C compile

$(XBIN): $(OBJS)
	$(CC) $(CFLAGS) -o $(XBIN) $(OBJS) $(LDLIBS)

test:	$(XBIN)	
	./$(XBIN)

%.o:	%.[cS]
	$(CC) $(CFLAGS) -c $^ -o $@

#	Cleanup
clean:
	$(RM) $(XBIN) $(OBJS)
