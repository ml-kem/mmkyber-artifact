.PHONY: all clean

#       these are expected to come from command line
#       make MODE=PKE LEVEL=192
	
#       options: { KEM, PKE }
MODE    ?=      KEM

#       options: { 128, 192, 256 }
LEVEL   ?=      128


MM_C_DIR	=	../mmKyber-c
LAZER_DIR	=	../../lazer
THIS_DIR	=	$(shell pwd)

CSRC		=	$(wildcard *.c $(MM_C_DIR)/*.c $(MM_C_DIR)/sym/*.c)
OBJS		=	$(CSRC:.c=.o)

CFLAGS		+=	-Wall -Wextra -Wshadow -Wundef -O3 -g
#CFLAGS		+=	-Wall -Wextra -Wshadow -fsanitize=address,undefined -O2 -g
CFLAGS		+=	-I. -I$(LAZER_DIR)
CFLAGS		+=	-I$(MM_C_DIR) -I$(MM_C_DIR)/sym -DMM_$(MODE) -DMM_$(LEVEL)
CFLAGS		+=	-DMM_L2NORM -DMM_BIN_LWE

LDLIBS		=	$(LAZER_DIR)/liblazer.a -lmpfr -lgmp -lm -lstdc++ \
	$(LAZER_DIR)/third_party/hexl-development/build/hexl/lib/libhexl.a

all:		mm_proof

mm_proof:	$(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDLIBS)

%.o:	%.[csS]
	$(CC) $(CFLAGS) -c $^ -o $@

#	running the code generator

clean-all:	clean
	rm	-f $(CODEGEN_H) proof.log

clean:
	rm	-f mm_proof $(OBJS)

