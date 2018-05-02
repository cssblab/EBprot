TOPDIR := $(shell pwd)
TARGET := $(TOPDIR)/bin

all: 
	mkdir -p $(TARGET)
	cd EBprot.MakeGrpData; make
	cd EBprotV2; make

clean:
	cd EBprot.MakeGrpData; rm -rf *.o
	cd EBprotV2; rm -rf *.o
	rm -rf $(TARGET)

