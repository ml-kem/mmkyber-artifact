#	Makefile

clean:
	$(RM) -f $(XBIN) $(OBJS)
	cd mmKyber-c && $(MAKE) clean
	cd mmKyber-py && $(MAKE) clean
	cd mmKyber-pkzk && $(MAKE) clean

