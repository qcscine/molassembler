.PHONY: build push run

build:
	docker build -t gitlab.chab.ethz.ch:4567/sobezj/molassembler/gcc7 .

push:
	docker push gitlab.chab.ethz.ch:4567/sobezj/molassembler/gcc7

run:
	docker run -it gitlab.chab.ethz.ch:4567/sobezj/molassembler/gcc7
